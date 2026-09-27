################################################################################
# AURORA_nerve.R
# 单独脚本，不要和 nerve_TCGA.R / TCGA-BRCA-2026.R 在同一会话 Source
# 数据目录：E:/R/cBioportal breast cancer/AURORA
# 也自动识别子目录 brca_aurora_2023 / brcaaurora_2023（cBioPortal AURORA US, Nat Cancer 2023）
# RStudio 打开后从第一行 Source
#
# 需要的文件（cBioPortal 下载即可；.txt / .txt.gz）：
#   data_clinical_sample.txt
#   data_clinical_patient.txt
#   data_mrna_seq_v2_rsem.txt                  # 优先，原始 RSEM
#   或 data_mrna_seq_v2_rsem_zscores_*.txt     # 仅有 z-score 时也可用
#
# 分析方式固定为 a–d（全部只用原位瘤 RNA）。
# 气泡图纵坐标只写通路英文名，不写 GO 编号。
# a–c 聚焦：
#   GO:2001224  positive regulation of neuron migration
#   GO:2001222  regulation of neuron migration
#   GO:0019227  neuronal action potential propagation
#   GO:0019228  neuronal action potential
#   GO:1902847  regulation of neuronal signal transduction
#   GO:0097492  sympathetic neuron axon guidance
#   GO:0007411  axon guidance
#   GO:0007409  axonogenesis
#   GO:0023041  neuronal signal transduction   # 用户写作 GO00230410
#   GO:1902667  regulation of axon guidance
# d 聚焦：
#   GO:2001222  regulation of neuron migration
#   GO:0031102  neuron projection regeneration
#   GO:0097491  sympathetic neuron projection guidance
#   GO:0097374  sensory neuron axon guidance
#   GO:1902667  regulation of axon guidance
#   GO:0031103  axon regeneration
#   GO:0007411  axon guidance
#   GO:0007409  axonogenesis
#
# 打分：主 z-mean；补充 z-median、ssGSEA
# 全部只用原位瘤 RNA。转移组织只用来判断该患者有没有配对转移灶，不参与打分。
# 转移定义（人数不足的组会跳过并写日志）：
#   a 原位 诊断 M1 vs M0（PRIM_M / YPM）
#   b 原位 Stage IV vs I-III
#   c 原位 N+ vs N0（YPN / PRIM_N）
#   d 主：原位 已转移 vs 未转移
#        已转移 = 诊断 M1 或 Stage IV，或该患者有配对转移组织
#        未转移 = 诊断 M0 且不是 Stage IV，且没有配对转移组织
#
# 负相关基因（原位 RNA；每个神经 GO 单独打分，不合并）：
#   与转移：按 a–d 同一套分组，Spearman r<0 且 p<0.05（主）；另出 r<=-0.15 且 p<0.05（严格）
#   与神经浸润：基因 vs 各神经 GO 分数，同一套阈值
#
# 结果目录：results_AURORA_nerve/
################################################################################

library(data.table)
library(ggplot2)

aurora_out_dir <- "results_AURORA_nerve"
min_set_genes <- 1
min_group_n <- 2
min_expr_frac <- 0.20
neg_pvalue_cutoff <- 0.05
neg_r_cutoff <- 0
strict_r_cutoff <- -0.15

# a–c 与 d 用不同的神经 GO 子集；打分只做并集，每个 GO 仍单独取基因
go_list_abc <- c(
  "GO:2001224",
  "GO:2001222",
  "GO:0019227",
  "GO:0019228",
  "GO:1902847",
  "GO:0097492",
  "GO:0007411",
  "GO:0007409",
  "GO:0023041",
  "GO:1902667"
)
go_list_d <- c(
  "GO:2001222",
  "GO:0031102",
  "GO:0097491",
  "GO:0097374",
  "GO:1902667",
  "GO:0031103",
  "GO:0007411",
  "GO:0007409"
)
go_list <- unique(c(go_list_abc, go_list_d))

gos_for_design <- function(key) {
  if (grepl("^d_", as.character(key))) go_list_d else go_list_abc
}

go_name_map <- c(
  "GO:0023041" = "neuronal signal transduction",
  "GO:1904457" = "positive regulation of neuronal action potential",
  "GO:1904340" = "positive regulation of dopaminergic neuron differentiation",
  "GO:2001224" = "positive regulation of neuron migration",
  "GO:2001222" = "regulation of neuron migration",
  "GO:0019227" = "neuronal action potential propagation",
  "GO:0019228" = "neuronal action potential",
  "GO:1902847" = "regulation of neuronal signal transduction",
  "GO:0031102" = "neuron projection regeneration",
  "GO:0097492" = "sympathetic neuron axon guidance",
  "GO:0097491" = "sympathetic neuron projection guidance",
  "GO:0097374" = "sensory neuron axon guidance",
  "GO:0007158" = "neuron cell-cell adhesion",
  "GO:1902667" = "regulation of axon guidance",
  "GO:0031103" = "axon regeneration",
  "GO:0007411" = "axon guidance",
  "GO:0007409" = "axonogenesis"
)

go_fallback <- list(
  "GO:0023041" = c(
    "MACO1", "NLGN1", "NRXN1", "OLFM1", "RTN4R", "KCNA1", "CLU", "DLG4",
    "GRIN1", "CAMK2A", "SYNGAP1", "MAPK1", "ADCY1", "HOMER1", "GNAO1", "PRKCA",
    "P2RY11"
  ),
  "GO:1904457" = c(
    "GBA1", "SCN1A", "SCN2A", "SCN8A", "ANK3", "FGF12", "KCNA1", "SCN1B",
    "CNTNAP1"
  ),
  "GO:1904340" = c(
    "DKK1", "FGF20", "EN1", "EN2", "LMX1A", "LMX1B", "NR4A2", "PITX3",
    "WNT1", "SHH", "FGF8"
  ),
  "GO:2001224" = c(
    "ADAT2", "ARHGEF2", "DAB2IP", "FLNA", "KIF20B", "MDK", "NIPBL", "NSMF",
    "PLAA", "RAPGEF2", "RELN", "SEMA6A", "SHTN1", "SRGAP2C", "TBC1D24", "WDR62",
    "ZNF609", "DCX", "CDK5", "PAFAH1B1", "DAB1", "NRG1", "CXCL12", "CXCR4"
  ),
  "GO:2001222" = c(
    "ADAT2", "ADGRG1", "ARHGEF2", "CAMK2A", "CAMK2B", "COL3A1", "CTNNA2", "CUL5",
    "CX3CL1", "DAB2IP", "FLNA", "GPR161", "IGSF10", "KIF20B", "KIF26A", "MDK",
    "NEXMIF", "NIPBL", "NSMF", "NTNG1", "NTNG2", "PHACTR1", "PLAA", "PLXNB2",
    "RAPGEF2", "RELN", "RNF7", "SEMA6A", "SHTN1", "SOCS7", "SRGAP2", "TBC1D24",
    "TNN", "ULK4", "VRK1", "WDR62", "ZNF609", "DCX", "CDK5", "PAFAH1B1",
    "CXCL12", "CXCR4", "SEMA3A", "SLIT1", "ROBO1", "GPR56", "NELF", "SRGAP2C"
  ),
  "GO:0019227" = c(
    "CNTNAP1", "SCN1B", "SCN1A", "SCN2A", "SCN8A", "SCN2B", "SCN4B", "ANK3",
    "NFASC", "NRCAM", "CNTN2", "SPTBN4", "KCNA1", "KCNQ2", "KCNQ3"
  ),
  "GO:0019228" = c(
    "ANK3", "ASIC5", "CHRNA1", "FGF12", "FMR1", "GBA1", "GPER1", "GRIA1",
    "KCNA1", "KCNA2", "KCND2", "KCNK2", "KCNK4", "KCNMB1", "KCNMB2", "KCNMB3",
    "KCNMB4", "MTNR1B", "MTOR", "MYH14", "NALCN", "NPR2", "NPY2R", "P2RX1",
    "SCN11A", "SCN1A", "SCN2A", "SCN8A", "SCN4B", "SCN9A", "TRPA1", "UNC79",
    "UNC80", "CNTNAP1", "SCN1B", "NFASC", "NRCAM", "NALF1"
  ),
  "GO:1902847" = c(
    "CLU", "ADCY1", "CAMK2A", "GRIN1", "GRIN2B", "DLG4", "HOMER1", "SYNGAP1",
    "RGS4", "GNAO1", "PRKCA", "MAPK1"
  ),
  "GO:0031102" = c(
    "GAP43", "CNTF", "EPHA4", "MAG", "RTN4", "RTN4R", "L1CAM", "NCAM1",
    "BDNF", "NGF", "NTRK1", "PTEN", "DLG4", "APOA4", "APOD", "CERS2",
    "CSPG5", "CTNNA1", "DAG1", "GRN", "INPP5F", "ISL1", "JAK2", "KIAA0319",
    "KREMEN1", "MAP1B", "MAPK8IP3", "MMP2", "NREP", "OMG", "PTN", "PTPRS",
    "PUM2", "RGMA", "RTCA", "RTN4RL1", "RTN4RL2", "SPP1", "STK24", "THY1",
    "TSPO", "STAT3", "ATF3", "DHFR", "FIGNL2", "FOLR1", "MTR", "SCARF1"
  ),
  "GO:0097492" = c(
    "ECE1", "EDN1", "EDNRA", "PLXNA4", "NRP1", "SEMA3A", "NGF", "NTRK1"
  ),
  "GO:0097491" = c(
    "NRP1", "NRP2", "SEMA3A", "SEMA3F", "ECE1", "EDN1", "EDNRA"
  ),
  "GO:0097374" = c(
    "NRP1", "NRP2", "SEMA3A", "SEMA3F", "PLXNA1", "PLXNA3", "PLXNA4", "NTN1",
    "DCC", "SLIT1", "ROBO1"
  ),
  "GO:0007158" = c(
    "NCAM1", "NCAM2", "L1CAM", "NRCAM", "CNTN1", "CNTN2", "CNTN4", "NRXN1",
    "NRXN2", "NRXN3", "NLGN1", "NLGN2", "NLGN3", "NLGN4X", "CADM1", "ASTN1",
    "ASTN2", "ICAM5", "CDK5R1", "ITGAL", "ITGB2", "NINJ2", "RET", "CDH2",
    "PCDH17", "KIAA0921", "NLGN4Y"
  ),
  "GO:1902667" = c(
    "ATOH7", "KIF21A", "MYCBP2", "NOVA2", "POU4F2", "PTPRO", "ROBO3", "SLIT2",
    "TUBB2B", "YTHDF1"
  ),
  "GO:0031103" = c(
    "GAP43", "PTEN", "MAP3K12", "STAT3", "JAK1", "SOCS3", "SPRR1A", "ATF3",
    "CNTF", "EPHA4", "MAG", "RTN4", "RTN4R", "L1CAM", "NCAM1", "BDNF",
    "NTRK1", "APOA4", "APOD", "CERS2", "CSPG5", "CTNNA1", "DAG1", "GRN",
    "ISL1", "JAK2", "MAP1B", "MMP2", "NREP", "PTN", "PTPRS", "RGMA",
    "RTN4RL1", "SPP1", "STK24", "DHFR", "FIGNL2", "FOLR1", "INPP5F", "KIAA0319",
    "KREMEN1", "MAPK8IP3", "MTR", "PUM2", "RTCA", "RTN4RL2", "SCARF1", "TSPO"
  ),
  "GO:0007411" = c(
    "SEMA3A", "SEMA3C", "SEMA3E", "SEMA3F", "SEMA4D", "SEMA5A", "SEMA6A", "SEMA6D",
    "SLIT1", "SLIT2", "SLIT3", "ROBO1", "ROBO2", "ROBO3", "NTN1", "NTN4",
    "DCC", "UNC5A", "UNC5B", "UNC5C", "UNC5D", "NRP1", "NRP2", "PLXNA1",
    "PLXNA3", "PLXNA4", "PLXNB1", "EFNA1", "EFNA2", "EFNA3", "EFNA4", "EFNA5",
    "EFNB1", "EFNB2", "EFNB3", "EPHA3", "EPHA4", "EPHA5", "EPHA6", "EPHA7",
    "EPHA8", "EPHB1", "EPHB2", "EPHB3", "CXCL12", "DRAXIN", "DSCAM", "ENAH",
    "FYN", "FZD3", "GAP43", "GDNF", "WNT5A", "ADAM17", "ALCAM", "ANOS1",
    "APP", "ARHGAP35", "ARHGEF25", "ARHGEF40", "ARK2C", "ATOH7", "BDNF", "BMPR2",
    "BOC", "BSG", "CDK5R1", "CDK5R2", "CELSR3", "CHN1", "CNTN1", "CNTN2",
    "CNTN4", "CNTN5", "CNTN6", "CSF1R", "CYFIP1", "CYFIP2", "DAG1", "DPYSL5",
    "DSCAML1", "ECE1", "EDN1", "EDN3", "EDNRA", "EMB", "EPHA10", "EPHB6",
    "EVL", "FEZ1", "FEZ2", "FGF8", "FLRT3", "GFRA3", "GLI2", "HMCN2",
    "IGSF9", "KALRN", "KIAA1755", "KIF21A", "KIF5C", "KLF7", "L1CAM", "LGI1",
    "LGR6", "LHX1", "LHX9", "LMO4", "LYPLA2", "MEGF8", "MYCBP2", "MYOT",
    "MYPN", "NCAM1", "NECTIN1", "NELL2", "NEO1", "NEXN", "NFASC", "NFIB",
    "NOTCH1", "NOTCH2", "NOTCH3", "NOVA2", "NPTN", "NRCAM", "NRXN1", "NRXN3",
    "NTN3", "NTRK1", "OPHN1", "OTX2", "PALLD", "PRKCQ", "PTCH1", "PTK2",
    "PTK7", "PTPRH", "PTPRJ", "PTPRM", "PTPRO", "RAC1", "RAC3", "RELN",
    "RET", "RIC1", "ROBO4", "RPS6KA5", "RYK", "SCN1B", "SEMA3B", "SEMA3D",
    "SEMA3G", "SEMA4A", "SEMA4B", "SEMA4C", "SEMA4F", "SEMA4G", "SEMA5B", "SEMA6B",
    "SEMA6C", "SEMA7A", "SHH", "SIAH1", "SMO", "SOS1", "TENM1", "TENM2",
    "TENM3", "TENM4", "TNR", "TRIO", "TUBB2B", "TUBB3", "USP33", "VANGL2",
    "VASP", "VEGFA", "WNT7B", "YTHDF1"
  ),
  "GO:0007409" = c(
    "ABL1", "ADGRB1", "ADNP", "ALCAM", "AMIGO1", "ANK3", "ANOS1", "APBB1",
    "APLP1", "APOE", "APP", "ATL1", "ATOH7", "AUTS2", "BAIAP2", "BDNF",
    "BMPR2", "BRSK1", "BRSK2", "CDH2", "CDK5", "CDK5R1", "CELSR1", "CELSR2",
    "CELSR3", "CHN1", "CNTN1", "CNTN2", "CNTNAP1", "CTNNA2", "CXCL12", "DCC",
    "DCLK1", "DISC1", "DOCK7", "DPYSL5", "DRAXIN", "DSCAM", "EFNA1", "EFNA2",
    "EFNA3", "EFNA4", "EFNA5", "EFNB1", "EFNB2", "EFNB3", "ENAH", "EPHA3",
    "EPHA4", "EPHA5", "EPHA6", "EPHA7", "EPHA8", "EPHB1", "EPHB2", "EPHB3",
    "FEZ1", "FEZ2", "FGF13", "FYN", "FZD3", "GAP43", "GDNF", "GSK3B",
    "ISL1", "ISL2", "KALRN", "KIF21A", "KIF5C", "KLF7", "L1CAM", "LHX1",
    "LIMK1", "LRRC4C", "MACF1", "MAP1B", "MAP2", "MAPT", "MYCBP2", "NCAM1",
    "NEFH", "NEO1", "NFASC", "NRCAM", "NRP1", "NRP2", "NRXN1", "NTN1",
    "NTN4", "NTNG1", "NTRK1", "NTRK2", "PAFAH1B1", "PAK1", "PLXNA3", "PLXNA4",
    "PLXNB1", "POU4F1", "POU4F2", "PTEN", "PTK2", "PTPRO", "RELN", "RET",
    "ROBO1", "ROBO2", "ROBO3", "RTN4", "RTN4R", "SEMA3A", "SEMA3C", "SEMA3E",
    "SEMA3F", "SEMA4D", "SEMA5A", "SEMA6A", "SEMA6D", "SHH", "SLIT1", "SLIT2",
    "SLIT3", "SLITRK1", "SPAST", "SPTBN4", "TENM1", "TENM2", "TENM3", "TENM4",
    "TUBB3", "ULK1", "UNC5A", "UNC5B", "UNC5C", "UNC5D", "WNT5A", "WNT7A",
    "ACTB", "ACTBL2", "ACTG1", "ACTL8", "ADAM17", "ADCY10", "AFG3L2", "ANAPC2",
    "APLP2", "ARHGAP35", "ARHGAP4", "ARHGEF25", "ARHGEF40", "ARK2C", "B4GALT5", "B4GALT6",
    "BOC", "BSG", "CCK", "CDH1", "CDHR2", "CDK5R2", "CDKL3", "CDKL5",
    "CHODL", "CHRNB2", "CNTN4", "CNTN5", "CNTN6", "COBL", "CSF1R", "CTTN",
    "CYFIP1", "CYFIP2", "DAG1", "DCHS2", "DIP2B", "DNM2", "DRD2", "DSCAML1",
    "ECE1", "EDN1", "EDN2", "EDN3", "EDNRA", "EMB", "EPHA10", "EPHB6",
    "EVL", "FAT4", "FGF8", "FGFR2", "FLRT3", "FOXB1", "GDI1", "GFRA3",
    "GLI2", "GOLGA4", "HEL113", "HMCN2", "IGSF9", "ISLR", "ISLR2", "ITGA4",
    "JUP", "KIAA0319", "KIAA1755", "KIF13B", "KIFBP", "LGI1", "LGR6", "LHX9",
    "LLGL1", "LMO4", "LRP4", "LYPLA2", "MAG", "MAP1A", "MAP1S", "MAP3K13",
    "MAP6", "MARK2", "MCF2", "MEGF8", "METRN", "MT3", "MYOT", "MYPN",
    "NECTIN1", "NELL2", "NEXN", "NFIB", "NLGN3", "NOLC1", "NOTCH1", "NOTCH2",
    "NOTCH3", "NOVA2", "NPR2", "NPTN", "NPTX1", "NRDC", "NRXN3", "NTN3",
    "NTNG2", "OLFM1", "OPHN1", "OTX2", "PAK2", "PAK3", "PALLD", "PARD3",
    "PAX2", "POTEE", "POTEF", "POTEI", "POTEJ", "PRKCA", "PRKCQ", "PRKG1",
    "PTCH1", "PTK7", "PTPRH", "PTPRJ", "PTPRM", "PTPRS", "RAB10", "RAB21",
    "RAB3A", "RAB8A", "RAC1", "RAC3", "RGMA", "RIC1", "RNF6", "ROBO4",
    "RPS6KA5", "RUFY3", "RYK", "S100B", "SCN11A", "SCN1B", "SEMA3B", "SEMA3D",
    "SEMA3G", "SEMA4A", "SEMA4B", "SEMA4C", "SEMA4F", "SEMA4G", "SEMA5B", "SEMA6B",
    "SEMA6C", "SEMA7A", "SHOX2", "SHTN1", "SIAH1", "SIN3A", "SIPA1L1", "SKIL",
    "SLC9A6", "SLITRK2", "SLITRK3", "SLITRK4", "SLITRK5", "SLITRK6", "SMN1", "SMO",
    "SMURF1", "SOS1", "SPG11", "SPP1", "SSNA1", "STK11", "STXBP1", "SZT2",
    "TAOK2", "THY1", "TIAM1", "TIAM2", "TNR", "TRAK1", "TRAK2", "TRIM46",
    "TRIO", "TRPC5", "TRPV2", "TSKU", "TUBB2B", "TWF2", "ULK2", "USP33",
    "USP9X", "VANGL2", "VASP", "VEGFA", "VIM", "WNT7B", "YTHDF1", "ZDHHC17",
    "ZFYVE27"
  )
)
go_fallback <- lapply(go_fallback, function(x) unique(trimws(x)))

go_title <- function(go_id) {
  go_id <- as.character(go_id)
  out <- unname(go_name_map[go_id])
  miss <- is.na(out) | !nzchar(out)
  out[miss] <- go_id[miss]
  out
}
# 图上只标英文通路名，不写 GO 编号
go_lab <- function(go_id) go_title(go_id)
safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

# ==============================================================================
# 工具函数（全部先定义，读完表达矩阵再分析）
# ==============================================================================
first_present <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) NA_character_ else hit[1]
}
norm_id <- function(x) {
  x <- toupper(gsub("[._ ]", "-", as.character(x)))
  gsub("-+", "-", x)
}
looks_like_cbioportal <- function(dir) {
  if (!dir.exists(dir)) return(FALSE)
  hits <- list.files(dir, pattern = "data_clinical_sample|data_mrna", ignore.case = TRUE)
  length(hits) > 0
}
resolve_aurora_dir <- function() {
  env <- Sys.getenv("AURORA_NERVE_DIR", unset = "")
  cands <- unique(c(
    env,
    "E:/R/cBioportal breast cancer/AURORA",
    "E:/R/cBioportal breast cancer/AURORA/brca_aurora_2023",
    "E:/R/cBioportal breast cancer/AURORA/brcaaurora_2023",
    file.path(getwd(), "brca_aurora_2023"),
    file.path(getwd(), "brcaaurora_2023"),
    getwd()
  ))
  cands <- cands[nzchar(cands)]
  extra <- unlist(lapply(cands, function(d) {
    if (!dir.exists(d)) return(character())
    subs <- list.dirs(d, recursive = FALSE, full.names = TRUE)
    named <- file.path(d, list.files(d, pattern = "aurora", ignore.case = TRUE, include.dirs = TRUE))
    c(subs, named)
  }))
  cands <- unique(c(cands, extra))
  hit <- cands[vapply(cands, looks_like_cbioportal, logical(1))]
  if (length(hit) == 0) {
    stop(
      "找不到 AURORA / cBioPortal 数据。请把脚本放到数据目录，或设置 AURORA_NERVE_DIR。\n",
      "已试：", paste(cands, collapse = " | ")
    )
  }
  normalizePath(hit[1], winslash = "/", mustWork = FALSE)
}

fread_cbioportal <- function(path, ...) {
  if (!file.exists(path)) stop("找不到文件：", path)
  npeek <- 30L
  peek <- tryCatch(readLines(path, n = npeek, warn = FALSE), error = function(e) character())
  n_hash <- sum(grepl("^#", peek))
  out <- tryCatch(
    fread(path, sep = "\t", header = TRUE, skip = n_hash, fill = TRUE, showProgress = TRUE, ...),
    error = function(e) e
  )
  if (!inherits(out, "error")) return(out)
  alt <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    sub("\\.gz$", "", path, ignore.case = TRUE)
  } else {
    paste0(path, ".gz")
  }
  if (file.exists(alt)) {
    message("读 ", path, " 失败，改读 ", alt)
    peek2 <- tryCatch(readLines(alt, n = npeek, warn = FALSE), error = function(e) character())
    n_hash2 <- sum(grepl("^#", peek2))
    return(fread(alt, sep = "\t", header = TRUE, skip = n_hash2, fill = TRUE, showProgress = TRUE, ...))
  }
  stop("读不了 ", path, "：", conditionMessage(out))
}

list_data_files <- function(dir, pattern) {
  hits <- list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  hits <- hits[!is.na(file.info(hits)$isdir) & !file.info(hits)$isdir]
  hits
}

pick_expression_file <- function(dir) {
  all_m <- list_data_files(dir, "^data_mrna.*\\.(txt|tsv)(\\.gz)?$")
  if (length(all_m) == 0) {
    stop("找不到表达矩阵：", dir, " 下应有 data_mrna_seq_*.txt")
  }
  rawish <- all_m[!grepl("zscore", basename(all_m), ignore.case = TRUE)]
  pick <- if (length(rawish) > 0) rawish else all_m
  info <- file.info(pick)
  pick <- pick[order(-info$size)]
  message("表达候选：\n", paste(sprintf("  %s  %.1fMB", pick, file.info(pick)$size / 1024^2), collapse = "\n"))
  pick[1]
}

pick_clin_file <- function(dir, stem) {
  hits <- list_data_files(dir, paste0("^", stem, "\\.(txt|tsv)(\\.gz)?$"))
  if (length(hits) == 0) return(NA_character_)
  hits[which.max(file.info(hits)$size)]
}

as_symbol_matrix_cbioportal <- function(expr_dt) {
  expr_dt <- as.data.table(expr_dt)
  nms <- names(expr_dt)
  gene_col <- first_present(nms, c(
    "Hugo_Symbol", "HUGO_SYMBOL", "hugo_symbol", "Gene_Symbol",
    "gene_symbol", "SYMBOL", "Gene", "gene"
  ))
  drop <- intersect(nms, c(
    gene_col, "Entrez_Gene_Id", "ENTREZ_GENE_ID", "entrez_gene_id",
    "Entrez", "ENSEMBL", "Ensembl", "ensembl"
  ))
  if (is.na(gene_col)) {
    gene_col <- nms[1]
    drop <- nms[1]
  }
  samp_cols <- setdiff(nms, drop)
  if (length(samp_cols) < 5) stop("表达表样本列太少：", length(samp_cols))
  sym <- as.character(expr_dt[[gene_col]])
  mat <- as.matrix(expr_dt[, samp_cols, with = FALSE])
  storage.mode(mat) <- "double"
  colnames(mat) <- norm_id(colnames(mat))
  keep <- !is.na(sym) & nzchar(sym) & !grepl("^NA$|^-$|^\\.$", sym)
  mat <- mat[keep, , drop = FALSE]
  sym <- sym[keep]
  if (anyDuplicated(sym)) {
    dt <- data.table(symbol = sym, as.data.table(mat))
    dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = symbol]
    mat <- as.matrix(dt[, -1, with = FALSE])
    rownames(mat) <- dt$symbol
  } else {
    rownames(mat) <- sym
  }
  if (anyDuplicated(colnames(mat))) {
    mat <- mat[, !duplicated(colnames(mat)), drop = FALSE]
  }
  mat
}

detect_already_z <- function(mat, path) {
  if (grepl("zscore", basename(path), ignore.case = TRUE)) return(TRUE)
  finite <- mat[is.finite(mat)]
  if (length(finite) < 100) return(FALSE)
  isTRUE(max(abs(finite), na.rm = TRUE) < 40 && abs(stats::median(finite, na.rm = TRUE)) < 1)
}

classify_sample_class <- function(sample_type, sample_id) {
  out <- rep(NA_character_, length(sample_type))
  xl <- toupper(as.character(sample_type))
  sid <- toupper(as.character(sample_id))
  out[grepl("METASTA", xl)] <- "Metastatic"
  out[grepl("PRIMARY|PRIMARY SOLID", xl)] <- "Primary"
  out[is.na(out) & grepl("-TTM|\\bTTM", sid)] <- "Metastatic"
  out[is.na(out) & grepl("-TTP|\\bTTP", sid)] <- "Primary"
  out
}
classify_m <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("M1", x)] <- "M1"
  out[is.na(out) & grepl("M0", x)] <- "M0"
  out[grepl("MX|UNKNOWN|NOT AVAILABLE|NOT REPORTED", x)] <- NA_character_
  out
}
classify_n <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("N[1-3]", x)] <- "Nplus"
  out[is.na(out) & grepl("N0", x)] <- "N0"
  out[grepl("NX|UNKNOWN|NOT AVAILABLE|NOT REPORTED", x)] <- NA_character_
  out
}
classify_stage <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("\\bIV\\b|STAGE\\s*4", x)] <- "Stage IV"
  out[is.na(out) & grepl("III|II|I[ABC]?\\b|STAGE\\s*[123]", x)] <- "Stage I-III"
  out[grepl("X\\b|UNKNOWN|NOT AVAILABLE|NOT REPORTED", x)] <- NA_character_
  out
}

# 通路活性。already_z=TRUE 时不再做 z（cBioPortal 已是全样本 z-score）。变量不要叫 score
pathway_zstat <- function(expr_mat, genes, how = "mean", already_z = FALSE) {
  genes <- unique(intersect(as.character(genes), rownames(expr_mat)))
  if (length(genes) < min_set_genes) return(NULL)
  sub <- as.matrix(expr_mat[genes, , drop = FALSE])
  storage.mode(sub) <- "double"
  if (isTRUE(already_z)) {
    z <- sub
    z[!is.finite(z)] <- 0
  } else {
    gene_mean <- rowMeans(sub, na.rm = TRUE)
    gene_sd <- sqrt(rowMeans((sub - gene_mean)^2, na.rm = TRUE))
    gene_sd[!is.finite(gene_sd) | gene_sd < 1e-12] <- 1
    z <- (sub - gene_mean) / gene_sd
    z[!is.finite(z)] <- 0
  }
  if (identical(how, "median")) {
    set_score <- apply(z, 2, stats::median, na.rm = TRUE)
  } else {
    set_score <- colMeans(z, na.rm = TRUE)
  }
  names(set_score) <- colnames(sub)
  attr(set_score, "n_genes") <- length(genes)
  attr(set_score, "genes") <- genes
  set_score
}
pathway_zmean <- function(expr_mat, genes, already_z = FALSE) {
  pathway_zstat(expr_mat, genes, "mean", already_z)
}
pathway_zmedian <- function(expr_mat, genes, already_z = FALSE) {
  pathway_zstat(expr_mat, genes, "median", already_z)
}

ssgsea_via_gsva <- function(expr_mat, gene_sets) {
  if (!requireNamespace("GSVA", quietly = TRUE)) return(NULL)
  gsets <- lapply(gene_sets, function(g) unique(intersect(as.character(g), rownames(expr_mat))))
  gsets <- Filter(function(g) length(g) >= min_set_genes, gsets)
  if (length(gsets) == 0) return(NULL)
  mat <- as.matrix(expr_mat)
  storage.mode(mat) <- "double"
  scored <- tryCatch({
    if (exists("ssgseaParam", envir = asNamespace("GSVA"), inherits = FALSE)) {
      param <- GSVA::ssgseaParam(mat, gsets, normalize = TRUE)
      GSVA::gsva(param, verbose = FALSE)
    } else {
      GSVA::gsva(mat, gsets, method = "ssgsea", ssgsea.norm = TRUE, verbose = FALSE)
    }
  }, error = function(e) {
    message("GSVA ssGSEA 失败，改用脚本内实现：", conditionMessage(e))
    NULL
  })
  if (is.null(scored)) return(NULL)
  scored <- as.matrix(scored)
  # GSVA: rows = gene sets, columns = samples
  out <- lapply(rownames(scored), function(nm) {
    v <- as.numeric(scored[nm, ])
    names(v) <- colnames(scored)
    attr(v, "n_genes") <- length(gsets[[nm]])
    attr(v, "genes") <- gsets[[nm]]
    v
  })
  names(out) <- rownames(scored)
  out
}

ssgsea_builtin <- function(expr_mat, gene_sets, tau = 0.25) {
  mat <- as.matrix(expr_mat)
  storage.mode(mat) <- "double"
  genes_all <- rownames(mat)
  ng <- nrow(mat)
  ns <- ncol(mat)
  set_idx <- lapply(gene_sets, function(g) unique(which(genes_all %in% unique(as.character(g)))))
  keep <- vapply(set_idx, length, integer(1)) >= min_set_genes
  set_idx <- set_idx[keep]
  if (length(set_idx) == 0) return(list())
  out_mat <- matrix(NA_real_, ns, length(set_idx),
                    dimnames = list(colnames(mat), names(set_idx)))
  pos_tau <- seq_len(ng)^tau
  for (j in seq_len(ns)) {
    o <- order(mat[, j], decreasing = TRUE, na.last = TRUE)
    for (k in seq_along(set_idx)) {
      hit <- o %in% set_idx[[k]]
      n_hit <- sum(hit)
      n_miss <- ng - n_hit
      if (n_hit < min_set_genes || n_miss < 1) next
      hit_w <- numeric(ng)
      hit_w[hit] <- pos_tau[hit]
      s <- sum(hit_w)
      if (!is.finite(s) || s <= 0) next
      walk <- cumsum(hit_w / s - ifelse(hit, 0, 1 / n_miss))
      out_mat[j, k] <- sum(walk)
    }
  }
  for (k in seq_len(ncol(out_mat))) {
    mx <- max(abs(out_mat[, k]), na.rm = TRUE)
    if (is.finite(mx) && mx > 0) out_mat[, k] <- out_mat[, k] / mx
  }
  out <- lapply(colnames(out_mat), function(nm) {
    v <- as.numeric(out_mat[, nm])
    names(v) <- rownames(out_mat)
    attr(v, "n_genes") <- length(set_idx[[nm]])
    attr(v, "genes") <- genes_all[set_idx[[nm]]]
    v
  })
  names(out) <- colnames(out_mat)
  out
}

score_ssgsea_sets <- function(expr_mat, gene_sets) {
  via <- ssgsea_via_gsva(expr_mat, gene_sets)
  if (!is.null(via) && length(via) > 0) {
    message("  ssGSEA 后端：GSVA")
    return(via)
  }
  message("  ssGSEA 后端：脚本内置（未安装 GSVA）")
  ssgsea_builtin(expr_mat, gene_sets)
}

get_go_genes <- function(go_id) {
  mapped <- NULL
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    pick <- function(keytype) {
      tryCatch(
        AnnotationDbi::select(
          org.Hs.eg.db::org.Hs.eg.db, keys = go_id, keytype = keytype,
          columns = c("SYMBOL", "ENSEMBL", "ENTREZID")
        ),
        error = function(e) NULL
      )
    }
    mapped <- suppressMessages(pick("GOALL"))
    if (is.null(mapped) || nrow(mapped) == 0) mapped <- suppressMessages(pick("GO"))
    extra <- tryCatch({
      go2eg <- as.list(org.Hs.eg.db::org.Hs.egGO2ALLEGS)
      eg <- unique(as.character(go2eg[[go_id]]))
      if (length(eg) == 0) NULL else {
        AnnotationDbi::select(
          org.Hs.eg.db::org.Hs.eg.db, keys = eg, keytype = "ENTREZID",
          columns = c("SYMBOL", "ENSEMBL")
        )
      }
    }, error = function(e) NULL)
    if (!is.null(extra) && nrow(extra) > 0) {
      extra <- as.data.table(extra)
      extra[, GO := go_id]
      mapped <- if (is.null(mapped) || nrow(mapped) == 0) extra else {
        rbind(as.data.table(mapped), extra, fill = TRUE)
      }
    }
  }
  if (!is.null(mapped) && nrow(mapped) > 0) {
    mapped <- as.data.table(mapped)
    if ("GOALL" %in% names(mapped)) setnames(mapped, "GOALL", "GO", skip_absent = TRUE)
    if (!"GO" %in% names(mapped)) mapped[, GO := go_id]
    syms <- unique(mapped[!is.na(SYMBOL) & SYMBOL != "", SYMBOL])
    if (length(syms) > 0) return(list(genes = syms, source = "org.Hs.eg.db"))
  }
  fb <- go_fallback[[go_id]]
  if (is.null(fb) || length(fb) == 0) return(list(genes = character(), source = "none"))
  list(genes = unique(fb), source = "fallback")
}

compare_groups <- function(value_vec, group, pos, neg, grouping) {
  df <- data.frame(
    nerve_value = as.numeric(value_vec),
    group = as.character(group),
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$nerve_value) & df$group %in% c(pos, neg), ]
  n_pos <- sum(df$group == pos)
  n_neg <- sum(df$group == neg)
  if (n_pos < 1 && n_neg < 1) return(NULL)
  pval <- NA_real_
  if (n_pos >= min_group_n && n_neg >= min_group_n) {
    wt <- suppressWarnings(stats::wilcox.test(nerve_value ~ group, data = df))
    pval <- wt$p.value
  }
  med_pos <- if (n_pos > 0) stats::median(df$nerve_value[df$group == pos], na.rm = TRUE) else NA_real_
  med_neg <- if (n_neg > 0) stats::median(df$nerve_value[df$group == neg], na.rm = TRUE) else NA_real_
  data.table(
    grouping = grouping, pos_level = pos, neg_level = neg,
    n_pos = n_pos, n_neg = n_neg,
    median_pos = med_pos, median_neg = med_neg,
    delta_median = med_pos - med_neg,
    pvalue = pval, test = "wilcoxon_unpaired"
  )
}

spearman_vs_go_score <- function(mat, go_score_vec) {
  common <- intersect(colnames(mat), names(go_score_vec))
  if (length(common) < 5) return(data.table())
  mat <- mat[, common, drop = FALSE]
  go_score_vec <- go_score_vec[common]
  keep <- apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
  mat <- mat[keep, , drop = FALSE]
  n <- ncol(mat)
  r <- as.numeric(cor(
    base::t(as.matrix(mat)), go_score_vec,
    method = "spearman", use = "pairwise.complete.obs"
  ))
  names(r) <- rownames(mat)
  r <- pmin(pmax(r, -0.999999), 0.999999)
  tstat <- r * sqrt((n - 2) / pmax(1e-12, 1 - r^2))
  p <- 2 * stats::pt(-abs(tstat), df = n - 2)
  data.table(feature = names(r), spearman_r = r, pvalue = p,
             fdr = p.adjust(p, method = "BH"), n = n)
}

spearman_vs_binary <- function(mat, group, pos, neg) {
  g <- as.character(group)
  names(g) <- names(group)
  common <- intersect(colnames(mat), names(g))
  g <- g[common]
  keep_s <- g %in% c(pos, neg)
  if (sum(g[keep_s] == pos) < min_group_n || sum(g[keep_s] == neg) < min_group_n) {
    return(data.table())
  }
  y <- ifelse(g[keep_s] == pos, 1, 0)
  mat <- mat[, names(y), drop = FALSE]
  keep <- apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
  mat <- mat[keep, , drop = FALSE]
  n <- ncol(mat)
  r <- as.numeric(cor(
    base::t(as.matrix(mat)), y,
    method = "spearman", use = "pairwise.complete.obs"
  ))
  names(r) <- rownames(mat)
  r <- pmin(pmax(r, -0.999999), 0.999999)
  tstat <- r * sqrt((n - 2) / pmax(1e-12, 1 - r^2))
  p <- 2 * stats::pt(-abs(tstat), df = n - 2)
  data.table(
    feature = names(r), spearman_r = r, pvalue = p,
    fdr = p.adjust(p, method = "BH"), n = n,
    n_pos = sum(y == 1), n_neg = sum(y == 0),
    pos_level = pos, neg_level = neg
  )
}

head_dt <- function(dt, n) {
  if (is.null(dt) || nrow(dt) == 0L || n <= 0L) return(dt[0])
  dt[seq_len(min(as.integer(n), nrow(dt)))]
}

pick_volcano_labels <- function(plot_dt, n_neg = 12L, n_pos = 6L) {
  neg <- head_dt(plot_dt[significant_neg == TRUE], n_neg)
  pos <- head_dt(plot_dt[spearman_r > 0][order(-spearman_r)], n_pos)
  out <- unique(rbindlist(list(neg, pos), fill = TRUE))
  if (nrow(out) == 0L) return(out)
  out[is.finite(spearman_r) & is.finite(neglogp) & !is.na(feature) & nzchar(as.character(feature))]
}

# 不用 ggrepel：ggplot2 4.x 下会报 depth(NULL)
add_volcano_labels <- function(p, lab_dt) {
  if (is.null(lab_dt) || nrow(lab_dt) == 0L) return(p)
  p + ggplot2::geom_text(
    data = as.data.frame(lab_dt),
    aes(x = spearman_r, y = neglogp, label = feature),
    size = 2.3, vjust = -0.55, inherit.aes = FALSE, check_overlap = TRUE
  )
}

plot_gene_volcano <- function(tab, title, subtitle, path_stub, n_neg = 12L, n_pos = 6L) {
  if (is.null(tab) || nrow(tab) == 0) return(invisible(NULL))
  plot_dt <- copy(as.data.table(tab))
  plot_dt[, neglogp := pmin(12, -log10(pmax(pvalue, 1e-12)))]
  plot_dt[, col := ifelse(significant_neg, "Negative",
                          ifelse(spearman_r > 0 & pvalue < neg_pvalue_cutoff, "Positive", "NS"))]
  top_lab <- pick_volcano_labels(plot_dt, n_neg, n_pos)
  p_vol <- tryCatch({
    p <- ggplot(plot_dt, aes(x = spearman_r, y = neglogp, color = col)) +
      geom_point(alpha = 0.45, size = 0.7) +
      geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
      geom_hline(yintercept = -log10(neg_pvalue_cutoff), linetype = 2, color = "grey50") +
      scale_color_manual(values = c("Negative" = "#3C5488", "Positive" = "#E64B35", "NS" = "grey75")) +
      labs(title = title, subtitle = subtitle,
           x = "Spearman r", y = expression(-log[10](p)), color = NULL) +
      theme_bw()
    add_volcano_labels(p, top_lab)
  }, error = function(e) {
    message("火山图失败（已跳过）：", conditionMessage(e))
    NULL
  })
  if (!is.null(p_vol)) save_plot(p_vol, path_stub, 10, 7)
  invisible(p_vol)
}

compare_paired_patients <- function(value_vec, ann, grouping) {
  dt <- data.table(
    sample = names(value_vec),
    nerve_value = as.numeric(value_vec)
  )
  dt <- merge(dt, ann[, .(sample, patient, sample_class)], by = "sample", all.x = TRUE)
  dt <- dt[is.finite(nerve_value) & sample_class %in% c("Primary", "Metastatic")]
  if (nrow(dt) == 0) return(NULL)
  pat <- dt[, .(
    primary_mean = mean(nerve_value[sample_class == "Primary"], na.rm = TRUE),
    met_mean = mean(nerve_value[sample_class == "Metastatic"], na.rm = TRUE)
  ), by = patient]
  pat <- pat[is.finite(primary_mean) & is.finite(met_mean)]
  n <- nrow(pat)
  if (n < min_group_n) return(NULL)
  pval <- suppressWarnings(stats::wilcox.test(pat$met_mean, pat$primary_mean, paired = TRUE)$p.value)
  data.table(
    grouping = grouping, pos_level = "Metastatic", neg_level = "Primary",
    n_pos = n, n_neg = n,
    median_pos = stats::median(pat$met_mean, na.rm = TRUE),
    median_neg = stats::median(pat$primary_mean, na.rm = TRUE),
    delta_median = stats::median(pat$met_mean - pat$primary_mean, na.rm = TRUE),
    pvalue = pval, test = "wilcoxon_paired"
  )
}

save_plot <- function(p, path_stub, width = 11, height = 8) {
  tryCatch(
    ggsave(paste0(path_stub, ".pdf"), p, width = width, height = height),
    error = function(e) message("保存 PDF 失败：", conditionMessage(e))
  )
  tryCatch(
    ggsave(paste0(path_stub, ".png"), p, width = width, height = height, dpi = 150),
    error = function(e) message("保存 PNG 失败：", conditionMessage(e))
  )
  if (interactive()) {
    tryCatch(print(p), error = function(e) {
      message("预览图失败（已跳过）：", conditionMessage(e))
    })
  }
}

expand_two_cols <- function(stat_dt) {
  d <- copy(as.data.table(stat_dt))
  if (nrow(d) == 0) return(d)
  rbindlist(list(
    d[, .(
      GO, GO_name, grouping, panel,
      side = "Non-metastatic", x_lab = neg_lab,
      n = n_neg, median_value = median_neg, pvalue
    )],
    d[, .(
      GO, GO_name, grouping, panel,
      side = "Metastatic", x_lab = pos_lab,
      n = n_pos, median_value = median_pos, pvalue
    )]
  ), fill = TRUE)
}

plot_bubble_two_cols <- function(stat_dt, title, subtitle, path_stub, facet = FALSE,
                                 go_order = NULL) {
  long <- expand_two_cols(stat_dt)
  long <- long[is.finite(median_value)]
  if (nrow(long) == 0) return(invisible(NULL))
  present <- unique(as.character(stat_dt$GO))
  if (!is.null(go_order) && length(go_order) > 0) {
    go_ids <- c(go_order[go_order %in% present], present[!present %in% go_order])
  } else {
    go_ids <- present
  }
  go_lv <- unique(go_lab(go_ids))
  long[, y_lab := factor(go_lab(GO), levels = rev(go_lv))]
  long[, x_lab := factor(x_lab, levels = unique(c(stat_dt$neg_lab, stat_dt$pos_lab)))]
  long[, neglogp := ifelse(is.finite(pvalue), pmin(10, -log10(pmax(pvalue, 1e-12))), 0.5)]
  fill_lim <- max(abs(long$median_value), na.rm = TRUE)
  if (!is.finite(fill_lim) || fill_lim < 1e-8) fill_lim <- 0.1
  p <- ggplot(long, aes(x = x_lab, y = y_lab)) +
    geom_point(
      aes(size = neglogp, fill = median_value),
      shape = 21, color = "black", stroke = 0.5 / ggplot2::.pt
    ) +
    scale_fill_gradientn(
      colours = c("#3C5488", "#5B7FA6", "#FFFFFF", "#EE8A7A", "#E64B35"),
      values = c(0, 0.47, 0.50, 0.53, 1),
      limits = c(-fill_lim, fill_lim),
      oob = scales::squish,
      name = "Pathway score\nmedian"
    ) +
    scale_size_continuous(range = c(3, 11), name = expression(-log[10](p))) +
    scale_x_discrete(expand = expansion(add = 0.85)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1, size = 11),
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      strip.text = element_text(size = 10),
      panel.spacing.x = grid::unit(0.55, "lines")
    )
  if (isTRUE(facet) && "panel" %in% names(long) && uniqueN(long$panel) > 1) {
    # a–c 与 d 的 GO 子集不同，分面允许各自 y 轴
    p <- p + facet_wrap(~ panel, nrow = 1, scales = "free")
  }
  n_panel <- if (isTRUE(facet)) max(1, uniqueN(long$panel)) else 1
  n_y <- uniqueN(long$y_lab)
  fig_w <- if (isTRUE(facet)) max(10.2, 2.9 * n_panel + 3.6) else 8.2
  fig_h <- max(6.4, 0.48 * n_y + 2.8)
  save_plot(p, path_stub, fig_w, fig_h)
}

build_aurora_annotation <- function(sample_ids, clin_sample, clin_patient) {
  ann <- data.table(sample = norm_id(sample_ids))
  ann <- ann[!is.na(sample) & sample != ""]
  ann <- ann[!duplicated(sample)]

  cs <- as.data.table(clin_sample)
  names(cs) <- toupper(names(cs))
  sid <- first_present(names(cs), c("SAMPLE_ID", "SAMPLE", "SAMPLEID"))
  pid <- first_present(names(cs), c("PATIENT_ID", "PATIENT", "PATIENTID"))
  if (is.na(sid)) stop("临床样本表没有 SAMPLE_ID")
  cs[, sample := norm_id(cs[[sid]])]
  cs[, patient := if (!is.na(pid)) as.character(cs[[pid]]) else sub("-(TTP|TTM).*$", "", sample, ignore.case = TRUE)]
  st <- first_present(names(cs), c("SAMPLE_TYPE", "TUMOR_TYPE", "TISSUE_TYPE"))
  ms <- first_present(names(cs), c("METASTATIC_SITE", "METASTATIC.SITE"))
  keep_s <- unique(c("sample", "patient", na.omit(c(st, ms))))
  extra_s <- cs[, keep_s, with = FALSE]
  extra_s <- extra_s[!duplicated(sample)]
  ann <- merge(ann, extra_s, by = "sample", all.x = TRUE)

  if (!is.null(clin_patient) && nrow(clin_patient) > 0) {
    cp <- as.data.table(clin_patient)
    names(cp) <- toupper(names(cp))
    pp <- first_present(names(cp), c("PATIENT_ID", "PATIENT", "PATIENTID"))
    if (!is.na(pp)) {
      cp[, patient := as.character(cp[[pp]])]
      want <- c(
        "PRIM_M", "YPM", "PRIM_N", "YPN", "PRIM_STAGE_DX",
        "YSTAGE_AT_PSTAGING", "DFS_STATUS", "OS_STATUS"
      )
      keep_p <- unique(c("patient", intersect(want, names(cp))))
      extra_p <- cp[, keep_p, with = FALSE]
      extra_p <- extra_p[!duplicated(patient)]
      ann <- merge(ann, extra_p, by = "patient", all.x = TRUE)
    }
  }

  raw_type <- if (!is.na(st) && st %in% names(ann)) ann[[st]] else rep(NA_character_, nrow(ann))
  type_vec <- classify_sample_class(raw_type, ann$sample)
  m_raw <- if ("PRIM_M" %in% names(ann)) ann$PRIM_M else rep(NA_character_, nrow(ann))
  if ("YPM" %in% names(ann)) {
    miss <- is.na(classify_m(m_raw))
    m_raw[miss] <- ann$YPM[miss]
  }
  n_raw <- if ("YPN" %in% names(ann)) ann$YPN else rep(NA_character_, nrow(ann))
  if ("PRIM_N" %in% names(ann)) {
    miss <- is.na(classify_n(n_raw))
    n_raw[miss] <- ann$PRIM_N[miss]
  }
  st_raw <- if ("PRIM_STAGE_DX" %in% names(ann)) ann$PRIM_STAGE_DX else rep(NA_character_, nrow(ann))
  if ("YSTAGE_AT_PSTAGING" %in% names(ann)) {
    miss <- is.na(classify_stage(st_raw))
    st_raw[miss] <- ann$YSTAGE_AT_PSTAGING[miss]
  }

  ann[, sample_class := factor(type_vec, levels = c("Primary", "Metastatic"))]
  ann[, distant_M := factor(classify_m(m_raw), levels = c("M0", "M1"))]
  ann[, node_N := factor(classify_n(n_raw), levels = c("N0", "Nplus"))]
  ann[, stage_IV := factor(classify_stage(st_raw), levels = c("Stage I-III", "Stage IV"))]
  if (is.na(ann$patient[1]) || !nzchar(ann$patient[1])) {
    ann[, patient := sub("-(TTP|TTM).*$", "", sample, ignore.case = TRUE)]
  }
  message(
    "选用临床：SAMPLE_TYPE=", st,
    "  M=PRIM_M/YPM  N=YPN/PRIM_N  stage=PRIM_STAGE_DX/YSTAGE"
  )
  ann
}

# 原位瘤转移标签：临床 M1/Stage IV，或该患者有配对转移组织
label_primary_met_status <- function(ann, clin_sample = NULL) {
  met_pats <- unique(na.omit(as.character(ann$patient[ann$sample_class == "Metastatic"])))
  if (!is.null(clin_sample) && nrow(clin_sample) > 0) {
    cs <- as.data.table(clin_sample)
    names(cs) <- toupper(names(cs))
    pid <- first_present(names(cs), c("PATIENT_ID", "PATIENT", "PATIENTID"))
    st <- first_present(names(cs), c("SAMPLE_TYPE", "TUMOR_TYPE", "TISSUE_TYPE"))
    sid <- first_present(names(cs), c("SAMPLE_ID", "SAMPLE", "SAMPLEID"))
    if (!is.na(pid) && !is.na(st)) {
      raw_id <- if (!is.na(sid)) cs[[sid]] else cs[[pid]]
      cls <- classify_sample_class(cs[[st]], raw_id)
      met_pats <- unique(c(met_pats, as.character(cs[[pid]])[cls == "Metastatic"]))
    }
  }
  met_pats <- unique(met_pats[nzchar(met_pats) & !is.na(met_pats)])
  ann[, has_paired_met := patient %in% met_pats]
  m1 <- as.character(ann$distant_M) == "M1"
  m0 <- as.character(ann$distant_M) == "M0"
  st4 <- as.character(ann$stage_IV) == "Stage IV"
  status <- rep(NA_character_, nrow(ann))
  status[m1 | st4 | ann$has_paired_met] <- "Metastasized"
  status[is.na(status) & m0 & !st4 & !ann$has_paired_met] <- "Non-metastasized"
  ann[, primary_met_status := factor(status, levels = c("Non-metastasized", "Metastasized"))]
  fwrite(data.table(
    field = c("primary_met_status", "has_paired_met"),
    meaning = c(
      "Primary RNA only. Metastasized = diagnosis M1 or Stage IV or patient has a paired metastatic tissue; Non-metastasized = M0, not Stage IV, and no paired metastatic tissue",
      "Patient has at least one metastatic sample in clinical/RNA (used only as a label, met RNA is not scored)"
    )
  ), file.path(aurora_out_dir, "00_primary_met_definition.csv"))
  ann
}

# ==============================================================================
# 读数据（分析函数在 expr 建好之后才调用）
# ==============================================================================
aurora_work_dir <- resolve_aurora_dir()
setwd(aurora_work_dir)
message("数据目录：", aurora_work_dir)

expr_file <- pick_expression_file(aurora_work_dir)
samp_file <- pick_clin_file(aurora_work_dir, "data_clinical_sample")
pat_file <- pick_clin_file(aurora_work_dir, "data_clinical_patient")
if (is.na(samp_file)) stop("找不到 data_clinical_sample.txt")

message("表达：", expr_file)
message("样本临床：", samp_file)
if (!is.na(pat_file)) message("患者临床：", pat_file)

expr_dt <- fread_cbioportal(expr_file)
clin_sample <- fread_cbioportal(samp_file)
clin_patient <- if (!is.na(pat_file)) fread_cbioportal(pat_file) else NULL

dir.create(aurora_out_dir, showWarnings = FALSE, recursive = TRUE)

message("正在把表达表映射为基因符号矩阵…")
aurora_expr_all <- as_symbol_matrix_cbioportal(expr_dt)
aurora_already_z <- detect_already_z(aurora_expr_all, expr_file)
mx <- suppressWarnings(max(aurora_expr_all, na.rm = TRUE))
if (!isTRUE(aurora_already_z) && is.finite(mx) && mx > 50) {
  message("检测到原始表达（max=", round(mx, 2), "），做 log2(x+1)")
  aurora_expr_all <- log2(pmax(aurora_expr_all, 0) + 1)
} else if (isTRUE(aurora_already_z)) {
  message("表达已是 z-score（", basename(expr_file), "），z-mean/z-median 不再重复标准化")
} else {
  message("表达值范围较小（max=", round(mx, 2), "），视为已 log 转换")
}

aurora_ann <- build_aurora_annotation(colnames(aurora_expr_all), clin_sample, clin_patient)
aurora_ann <- aurora_ann[sample %in% colnames(aurora_expr_all)]
aurora_ann <- label_primary_met_status(aurora_ann, clin_sample)
fwrite(aurora_ann, file.path(aurora_out_dir, "00_sample_annotation.csv"))

is_prim <- aurora_ann$sample_class == "Primary"
message(
  "样本：原发=", sum(is_prim, na.rm = TRUE),
  "  转移组织（只用于配对标签）=", sum(aurora_ann$sample_class == "Metastatic", na.rm = TRUE),
  "  原位已转移=", sum(is_prim & aurora_ann$primary_met_status == "Metastasized", na.rm = TRUE),
  "  原位未转移=", sum(is_prim & aurora_ann$primary_met_status == "Non-metastasized", na.rm = TRUE),
  "  原位中有配对转移组织=", sum(is_prim & aurora_ann$has_paired_met, na.rm = TRUE),
  "  原位中 M1=", sum(is_prim & aurora_ann$distant_M == "M1", na.rm = TRUE),
  "  原位中 Stage IV=", sum(is_prim & aurora_ann$stage_IV == "Stage IV", na.rm = TRUE),
  "  原位中 N+=", sum(is_prim & aurora_ann$node_N == "Nplus", na.rm = TRUE)
)
if (sum(is_prim & aurora_ann$primary_met_status == "Non-metastasized", na.rm = TRUE) < min_group_n) {
  message("注意：AURORA 几乎都是转移性乳腺癌，未转移原位瘤可能很少；主比较 d 人数不足会跳过，仍会尝试 a/b/c")
}

keep_g <- rowMeans(is.finite(aurora_expr_all), na.rm = TRUE) >= min_expr_frac
if (!isTRUE(aurora_already_z)) {
  keep_g <- keep_g & rowMeans(is.finite(aurora_expr_all) & aurora_expr_all > 0, na.rm = TRUE) >= min_expr_frac
}
aurora_expr <- aurora_expr_all[keep_g, , drop = FALSE]
primary_ids <- aurora_ann$sample[aurora_ann$sample_class == "Primary"]
met_ids <- aurora_ann$sample[aurora_ann$sample_class == "Metastatic"]
aurora_expr_primary <- aurora_expr[, intersect(primary_ids, colnames(aurora_expr)), drop = FALSE]
message(
  "表达矩阵：全部 ", ncol(aurora_expr), " 样本 x ", nrow(aurora_expr),
  " 基因；原发 ", ncol(aurora_expr_primary)
)
if (ncol(aurora_expr_primary) < 5) stop("有表达的原位样本太少：", ncol(aurora_expr_primary))

# ==============================================================================
# 主分析（必须在 expr 建好之后）
# ==============================================================================
run_aurora_nerve <- function() {
  if (!exists("aurora_expr_primary", inherits = TRUE)) stop("还没有 aurora_expr_primary，请从脚本开头 Source")

  mat_from_list <- function(lst) {
    if (length(lst) == 0) return(NULL)
    common <- Reduce(intersect, lapply(lst, names))
    if (length(common) == 0) return(NULL)
    mat <- do.call(cbind, lapply(lst, function(x) x[common]))
    colnames(mat) <- names(lst)
    rownames(mat) <- common
    mat
  }

  message("收集聚焦神经 GO 基因（a–c 与 d 子集不同；每个 GO 单独，不合并）")
  go_map <- lapply(go_list, get_go_genes)
  names(go_map) <- go_list
  go_map <- Filter(function(x) length(x$genes) >= min_set_genes, go_map)
  if (length(go_map) == 0) stop("没有任何神经 GO 能打分")
  gene_sets <- lapply(go_map, function(x) x$genes)
  fwrite(rbindlist(lapply(names(go_map), function(g) {
    data.table(
      GO = g, GO_name = go_title(g), gene_source = go_map[[g]]$source,
      n_genes = length(go_map[[g]]$genes),
      used_in_abc = g %in% go_list_abc,
      used_in_d = g %in% go_list_d,
      genes = paste(go_map[[g]]$genes, collapse = ";")
    )
  }), fill = TRUE), file.path(aurora_out_dir, "00_GO_genes_used.csv"))
  fwrite(data.table(
    design = c(rep("a-c", length(go_list_abc)), rep("d", length(go_list_d))),
    GO = c(go_list_abc, go_list_d),
    GO_name = go_title(c(go_list_abc, go_list_d))
  ), file.path(aurora_out_dir, "00_GO_focus_by_design.csv"))
  fwrite(data.table(
    method = c("zmean", "zmedian", "ssgsea"),
    role = c("primary", "supplement", "supplement"),
    description = c(
      "Gene-wise z-score then mean (skipped if matrix is already z-scored)",
      "Gene-wise z-score then median (skipped if matrix is already z-scored)",
      "ssGSEA (Barbie 2009); GSVA if installed, otherwise built-in"
    )
  ), file.path(aurora_out_dir, "00_scoring_methods.csv"))

  score_z_method <- function(expr_mat, how) {
    lst <- lapply(names(gene_sets), function(g) {
      fn <- if (identical(how, "median")) pathway_zmedian else pathway_zmean
      set_score <- fn(expr_mat, gene_sets[[g]], already_z = aurora_already_z)
      if (is.null(set_score)) {
        message("  ", g, " 映射基因不足，跳过")
        return(NULL)
      }
      message("  ", g, "  ", go_title(g), "  基因数=", attr(set_score, "n_genes"))
      set_score
    })
    names(lst) <- names(gene_sets)
    Filter(Negate(is.null), lst)
  }

  plot_one_method <- function(method, sm_p) {
    if (is.null(sm_p) || nrow(sm_p) == 0) {
      message("  ", method$title, " 没有可画的原位分数，跳过")
      return(invisible(NULL))
    }
    mdir <- file.path(aurora_out_dir, method$id)
    dir.create(mdir, showWarnings = FALSE, recursive = TRUE)
    fwrite(data.table(sample = rownames(sm_p), as.data.table(sm_p)),
           file.path(mdir, "01_pathway_scores_primary.csv"))

    ann_p <- aurora_ann[sample %in% rownames(sm_p) & sample_class == "Primary"]

    subset_score_mat <- function(sm, gos) {
      if (is.null(sm) || ncol(sm) == 0) return(sm)
      keep <- intersect(gos, colnames(sm))
      if (length(keep) == 0) return(sm[, 0, drop = FALSE])
      sm[, keep, drop = FALSE]
    }
    designs <- list(
      list(key = "a_distant_M", title = "Distant M at first diagnosis", panel = "1a Distant M",
           group = setNames(as.character(ann_p$distant_M), ann_p$sample),
           pos = "M1", neg = "M0", pos_lab = "M1", neg_lab = "M0",
           gos = go_list_abc,
           score_mat = subset_score_mat(sm_p, go_list_abc)),
      list(key = "b_AJCC_stageIV", title = "AJCC stage at diagnosis", panel = "1b AJCC stage",
           group = setNames(as.character(ann_p$stage_IV), ann_p$sample),
           pos = "Stage IV", neg = "Stage I-III",
           pos_lab = "Stage IV", neg_lab = "Stage I-III",
           gos = go_list_abc,
           score_mat = subset_score_mat(sm_p, go_list_abc)),
      list(key = "c_node_N", title = "Lymph node", panel = "1c Lymph node",
           group = setNames(as.character(ann_p$node_N), ann_p$sample),
           pos = "Nplus", neg = "N0", pos_lab = "N+", neg_lab = "N0",
           gos = go_list_abc,
           score_mat = subset_score_mat(sm_p, go_list_abc)),
      list(key = "d_primary_met_status", title = "Primary later / known metastasis",
           panel = "1d Primary met status",
           group = setNames(as.character(ann_p$primary_met_status), ann_p$sample),
           pos = "Metastasized", neg = "Non-metastasized",
           pos_lab = "Metastasized primary", neg_lab = "Non-metastasized primary",
           gos = go_list_d,
           score_mat = subset_score_mat(sm_p, go_list_d))
    )

    all_stat <- list()
    for (ds in designs) {
      message("气泡图：", method$title, " / ", ds$title)
      stat_rows <- list()
      sm <- ds$score_mat
      if (is.null(sm) || nrow(sm) == 0) {
        message("  无分数矩阵，跳过 ", ds$key)
        next
      }
      for (g in colnames(sm)) {
        value_vec <- as.numeric(sm[, g])
        names(value_vec) <- rownames(sm)
        one <- compare_groups(value_vec, ds$group[names(value_vec)], ds$pos, ds$neg, ds$key)
        if (is.null(one)) next
        if (one$n_pos < min_group_n || one$n_neg < min_group_n) {
          next
        }
        one[, `:=`(
          scoring = method$id, GO = g, GO_name = go_title(g),
          panel = ds$panel, pos_lab = ds$pos_lab, neg_lab = ds$neg_lab
        )]
        stat_rows[[g]] <- one
      }
      stat_dt <- rbindlist(stat_rows, fill = TRUE)
      if (nrow(stat_dt) == 0) {
        message("  分组人数不足，跳过 ", method$id, " / ", ds$key)
        next
      }
      stat_dt[, fdr := p.adjust(pvalue, method = "BH")]
      fwrite(stat_dt, file.path(mdir, paste0("02_", ds$key, "_GO_vs_metastasis.csv")))
      fwrite(expand_two_cols(stat_dt),
             file.path(mdir, paste0("02_", ds$key, "_GO_vs_metastasis_two_cols.csv")))
      all_stat[[ds$key]] <- stat_dt
      plot_bubble_two_cols(
        stat_dt,
        title = paste0(ds$title, ": neural GO (", method$title, ")"),
        subtitle = paste0(
          "AURORA US; scoring = ", method$title,
          "; Y = pathway English name (scored separately); X = two groups; fill = median; size = -log10(p)"
        ),
        path_stub = file.path(mdir, paste0("02_", ds$key, "_bubble")),
        go_order = ds$gos
      )
    }

    if (length(all_stat) == 0) return(invisible(NULL))
    bubble <- rbindlist(all_stat, fill = TRUE)
    fwrite(bubble, file.path(mdir, "02_summary_GO_vs_metastasis.csv"))
    fwrite(expand_two_cols(bubble),
           file.path(mdir, "02_summary_GO_vs_metastasis_two_cols.csv"))
    abc_keys <- c("a_distant_M", "b_AJCC_stageIV", "c_node_N")
    bubble_abc <- bubble[grouping %in% abc_keys]
    bubble_d <- bubble[grouping == "d_primary_met_status"]
    if (nrow(bubble_abc) > 0) {
      plot_bubble_two_cols(
        bubble_abc,
        title = paste0("AURORA a-c neural GO versus metastasis (", method$title, ")"),
        subtitle = paste0("Scoring = ", method$title, "; focused a-c GO subset; Y = English name"),
        path_stub = file.path(mdir, "02_summary_abc_bubble_GO_vs_metastasis"),
        facet = TRUE,
        go_order = go_list_abc
      )
    }
    if (nrow(bubble_d) > 0) {
      plot_bubble_two_cols(
        bubble_d,
        title = paste0("AURORA d neural GO versus metastasis (", method$title, ")"),
        subtitle = paste0("Scoring = ", method$title, "; focused d GO subset; Y = English name"),
        path_stub = file.path(mdir, "02_summary_d_bubble_GO_vs_metastasis"),
        go_order = go_list_d
      )
    }
    plot_bubble_two_cols(
      bubble,
      title = paste0("AURORA neural GO versus metastasis (", method$title, ")"),
      subtitle = paste0(
        "Scoring = ", method$title,
        "; a-c and d use different focused GO subsets; not pooled"
      ),
      path_stub = file.path(mdir, "02_summary_bubble_GO_vs_metastasis"),
      facet = TRUE
    )
    if (isTRUE(method$primary)) {
      fwrite(bubble, file.path(aurora_out_dir, "02_summary_GO_vs_metastasis.csv"))
      if (nrow(bubble_abc) > 0) {
        plot_bubble_two_cols(
          bubble_abc,
          title = "AURORA primary tumors a-c: neural GO versus metastasis (z-mean)",
          subtitle = "Primary RNA only; Y = English pathway name",
          path_stub = file.path(aurora_out_dir, "02_summary_abc_bubble_GO_vs_metastasis"),
          facet = TRUE,
          go_order = go_list_abc
        )
      }
      if (nrow(bubble_d) > 0) {
        plot_bubble_two_cols(
          bubble_d,
          title = "AURORA primary tumors d: neural GO versus metastasis (z-mean)",
          subtitle = "Primary RNA only; Metastasized = M1 or Stage IV or paired metastatic tissue",
          path_stub = file.path(aurora_out_dir, "02_summary_d_bubble_GO_vs_metastasis"),
          go_order = go_list_d
        )
      }
      plot_bubble_two_cols(
        bubble,
        title = "AURORA primary tumors: neural GO versus metastasis (z-mean)",
        subtitle = "Primary RNA only; a-c and d use different focused GO subsets",
        path_stub = file.path(aurora_out_dir, "02_summary_bubble_GO_vs_metastasis"),
        facet = TRUE
      )
    }
    bubble
  }

  score_methods <- list(
    list(id = "zmean", title = "z-mean", primary = TRUE),
    list(id = "zmedian", title = "z-median", primary = FALSE),
    list(id = "ssgsea", title = "ssGSEA", primary = FALSE)
  )
  score_primary_list <- NULL
  all_method_stats <- list()
  for (method in score_methods) {
    message("打分：", method$title, if (isTRUE(method$primary)) "（主）" else "（补充）")
    if (identical(method$id, "ssgsea")) {
      lst <- score_ssgsea_sets(aurora_expr_primary, gene_sets)
    } else {
      how <- if (identical(method$id, "zmedian")) "median" else "mean"
      lst <- score_z_method(aurora_expr_primary, how)
    }
    sm_p <- mat_from_list(lst)
    if (isTRUE(method$primary)) score_primary_list <- lst
    all_method_stats[[method$id]] <- plot_one_method(method, sm_p)
  }

  keep <- Filter(function(x) is.data.table(x) && nrow(x) > 0, all_method_stats)
  if (length(keep) > 0) {
    fwrite(rbindlist(keep, fill = TRUE),
           file.path(aurora_out_dir, "02_all_scoring_methods_vs_metastasis.csv"))
  }

  # ---- 2) 原位肿瘤内，按 a–d 分组，与转移负相关的基因 ----
  ann_p <- aurora_ann[sample %in% colnames(aurora_expr_primary) & sample_class == "Primary"]
  met_defs <- list(
    list(key = "a_distant_M", title = "2a Genes negatively correlated with distant M1 (primary tumors)",
         group = setNames(as.character(ann_p$distant_M), ann_p$sample),
         pos = "M1", neg = "M0"),
    list(key = "b_AJCC_stageIV", title = "2b Genes negatively correlated with Stage IV (primary tumors)",
         group = setNames(as.character(ann_p$stage_IV), ann_p$sample),
         pos = "Stage IV", neg = "Stage I-III"),
    list(key = "c_node_N", title = "2c Genes negatively correlated with N+ (primary tumors)",
         group = setNames(as.character(ann_p$node_N), ann_p$sample),
         pos = "Nplus", neg = "N0"),
    list(key = "d_primary_met_status",
         title = "2d Genes negatively correlated with later / known metastasis (primary tumors)",
         group = setNames(as.character(ann_p$primary_met_status), ann_p$sample),
         pos = "Metastasized", neg = "Non-metastasized")
  )
  met_neg_summary <- list()
  for (md in met_defs) {
    message("全基因组相关：", md$title)
    tab <- spearman_vs_binary(aurora_expr_primary, md$group, md$pos, md$neg)
    if (nrow(tab) == 0) {
      message("  分组人数不足，跳过 ", md$key)
      next
    }
    tab[, `:=`(
      design = md$key,
      significant_neg = spearman_r < neg_r_cutoff & pvalue < neg_pvalue_cutoff,
      strict_neg = spearman_r <= strict_r_cutoff & pvalue < neg_pvalue_cutoff
    )]
    setorder(tab, spearman_r)
    fwrite(tab, file.path(aurora_out_dir, paste0("03_", md$key, "_genes_vs_metastasis_all.csv")))
    fwrite(tab[significant_neg == TRUE],
           file.path(aurora_out_dir, paste0("03_", md$key, "_genes_NEG_vs_metastasis.csv")))
    fwrite(tab[strict_neg == TRUE],
           file.path(aurora_out_dir, paste0("03_", md$key, "_genes_NEG_strict_vs_metastasis.csv")))
    message("  负相关基因 ", sum(tab$significant_neg),
            "（严格 r<=", strict_r_cutoff, "：", sum(tab$strict_neg), "）")
    met_neg_summary[[md$key]] <- data.table(
      design = md$key, title = md$title,
      n_tested = nrow(tab), n_pos = tab$n_pos[1], n_neg_group = tab$n_neg[1],
      n_neg = sum(tab$significant_neg), n_neg_strict = sum(tab$strict_neg)
    )
    plot_gene_volcano(
      tab,
      title = md$title,
      subtitle = paste0(
        "Spearman: gene vs ", md$pos, " (1) / ", md$neg, " (0); blue = negative vs metastasis"
      ),
      path_stub = file.path(aurora_out_dir, paste0("03_", md$key, "_volcano_genes_vs_metastasis"))
    )
  }
  if (length(met_neg_summary) > 0) {
    fwrite(rbindlist(met_neg_summary, fill = TRUE),
           file.path(aurora_out_dir, "03_summary_neg_genes_vs_metastasis.csv"))
  }

  # ---- 3) 原位肿瘤内，与每个神经 GO 分数负相关的基因 ----
  dir.create(file.path(aurora_out_dir, "04_neg_vs_neural_per_GO"), showWarnings = FALSE)
  summary_neg <- list()
  if (is.null(score_primary_list) || length(score_primary_list) == 0) {
    message("没有可用的原位神经 GO 分数，跳过与神经浸润负相关的基因分析")
  } else {
    for (g in names(score_primary_list)) {
      message("与神经浸润负相关：", go_title(g))
      go_score_vec <- score_primary_list[[g]]
      tab <- spearman_vs_go_score(aurora_expr_primary, go_score_vec)
      if (nrow(tab) == 0) next
      tab[, `:=`(
        GO = g, GO_name = go_title(g),
        used_in_abc = g %in% go_list_abc,
        used_in_d = g %in% go_list_d,
        significant_neg = spearman_r < neg_r_cutoff & pvalue < neg_pvalue_cutoff,
        strict_neg = spearman_r <= strict_r_cutoff & pvalue < neg_pvalue_cutoff
      )]
      setorder(tab, spearman_r)
      gdir <- file.path(aurora_out_dir, "04_neg_vs_neural_per_GO",
                        paste0("GO_", safe_name(sub("GO:", "", g))))
      dir.create(gdir, showWarnings = FALSE)
      fwrite(tab, file.path(gdir, "genes_vs_neural_GO_all.csv"))
      fwrite(tab[significant_neg == TRUE], file.path(gdir, "genes_NEG_vs_neural_GO.csv"))
      fwrite(tab[strict_neg == TRUE], file.path(gdir, "genes_NEG_strict_vs_neural_GO.csv"))
      summary_neg[[g]] <- data.table(
        GO = g, GO_name = go_title(g),
        used_in_abc = g %in% go_list_abc,
        used_in_d = g %in% go_list_d,
        n_pathway_genes = attr(go_score_vec, "n_genes"),
        n_tested = nrow(tab),
        n_neg = sum(tab$significant_neg),
        n_neg_strict = sum(tab$strict_neg)
      )
      plot_gene_volcano(
        tab,
        title = "Genes negatively correlated with neural invasion",
        subtitle = paste0(go_title(g), "; blue = negative vs this GO score"),
        path_stub = file.path(gdir, "volcano_neg_vs_neural_GO"),
        n_neg = 10L, n_pos = 5L
      )
    }
  }
  if (length(summary_neg) > 0) {
    sum_dt <- rbindlist(summary_neg, fill = TRUE)
    fwrite(sum_dt, file.path(aurora_out_dir, "04_summary_neg_genes_vs_each_neural_GO.csv"))
    p_n <- ggplot(sum_dt, aes(x = n_neg, y = reorder(GO_name, n_neg))) +
      geom_col(fill = "#3C5488", width = 0.7) +
      labs(title = "Number of genes negatively correlated with each neural GO",
           subtitle = paste0("Spearman r < 0 and p < ", neg_pvalue_cutoff, "; GO sets not pooled"),
           x = "Number of negative genes", y = NULL) +
      theme_bw()
    save_plot(p_n, file.path(aurora_out_dir, "04_summary_neg_gene_counts"), 10, 6.5)
  }

  message("完成。结果目录：", normalizePath(aurora_out_dir, winslash = "/", mustWork = FALSE))
  message("主图（原位已转移 vs 未转移）：zmean/02_d_primary_met_status_bubble.png")
  message("转移负相关基因：03_*_genes_NEG_vs_metastasis.csv")
  message("神经 GO 负相关基因：04_neg_vs_neural_per_GO/")
  invisible(TRUE)
}

if (exists("aurora_expr_primary") && ncol(aurora_expr_primary) > 5) {
  run_aurora_nerve()
} else {
  stop("原位表达矩阵未建好，请从第一行完整 Source")
}
