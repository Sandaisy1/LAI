################################################################################
# nerve_TCGA.R
# 单独脚本，不要和 TCGA-BRCA-2026.R 在同一会话 Source
# 放到 E:/R/TCGA-BRCA-2026，RStudio 打开后从第一行 Source
#
# 数据（.tsv / .tsv.gz 均可；空的 .tsv 会跳过，改读 .gz）：
#   TCGA-BRCA.star_fpkm.tsv(.gz)
#   TCGA-BRCA.clinical.tsv(.gz)
#   TCGA-BRCA.survival.tsv(.gz)          # 可选
#   gencode.v36.annotation.gtf.gene.probemap
#
# 神经浸润按 6 大类分别打分（类与类不合并，类内各条目也不合并）：
#   1 病理 PNI（Xena 临床若无此列则跳过并写说明）
#   2 神经信号 GO（5 个 GO 各自打分）
#   3 施旺细胞 marker（核心 / 扩展 两套）
#   4 神经营养因子（配体 / 受体 两套）
#   5 轴突与神经元结构（神经丝 / 结构蛋白 两套）
#   6 交感 vs 感觉神经
#
# 每类对 4 种转移定义做 Wilcoxon，并画两列气泡图（Non-metastatic | Metastatic）
#   a 原位肿瘤 M1 vs M0
#   b 原位肿瘤 Stage IV vs I-III
#   c 原位肿瘤 N+ vs N0
#   d 转移组织 vs 原位肿瘤
#
# 结果目录：results_nerve_TCGA/
################################################################################

library(data.table)
library(ggplot2)

nerve_work_dir <- Sys.getenv(
  "NERVE_TCGA_DIR",
  unset = Sys.getenv("TCGA_BRCA_2026_DIR", unset = "E:/R/TCGA-BRCA-2026")
)
nerve_out_dir <- "results_nerve_TCGA"
min_set_genes <- 1
min_group_n <- 2
min_expr_frac <- 0.20

# ==============================================================================
# 6 类定义（基因集合写死；GO 有 org.Hs.eg.db 时用数据库，否则用后备列表）
# ==============================================================================
go_name_map <- c(
  "GO:0019227" = "neuronal action potential propagation",
  "GO:1902847" = "regulation of neuronal signal transduction",
  "GO:0097374" = "sensory neuron axon guidance",
  "GO:1902667" = "regulation of axon guidance",
  "GO:0007409" = "axonogenesis"
)

go_fallback <- list(
  "GO:0019227" = c(
    "CNTNAP1", "SCN1B", "SCN1A", "SCN2A", "SCN8A", "SCN2B", "SCN4B",
    "ANK3", "NFASC", "NRCAM", "CNTN2", "SPTBN4", "KCNA1", "KCNQ2", "KCNQ3"
  ),
  "GO:1902847" = c(
    "CLU", "ADCY1", "CAMK2A", "GRIN1", "GRIN2B", "DLG4",
    "HOMER1", "SYNGAP1", "RGS4", "GNAO1", "PRKCA", "MAPK1"
  ),
  "GO:0097374" = c(
    "NRP1", "NRP2", "SEMA3A", "SEMA3F", "PLXNA1", "PLXNA3", "PLXNA4",
    "NTN1", "DCC", "SLIT1", "ROBO1"
  ),
  "GO:1902667" = c(
    "ATOH7", "KIF21A", "MYCBP2", "NOVA2", "POU4F2", "PTPRO",
    "ROBO3", "SLIT2", "TUBB2B", "YTHDF1"
  ),
  "GO:0007409" = strsplit(paste(
    "ABL1,ADGRB1,ADNP,ALCAM,AMIGO1,ANK3,ANOS1,APBB1,APLP1,APOE,APP,ATL1,ATOH7,AUTS2,",
    "BAIAP2,BDNF,BMPR2,BRSK1,BRSK2,CDH2,CDK5,CDK5R1,CELSR1,CELSR2,CELSR3,CHN1,",
    "CNTN1,CNTN2,CNTNAP1,CTNNA2,CXCL12,DCC,DCLK1,DISC1,DOCK7,DPYSL5,DRAXIN,DSCAM,",
    "EFNA1,EFNA2,EFNA3,EFNA4,EFNA5,EFNB1,EFNB2,EFNB3,ENAH,EPHA3,EPHA4,EPHA5,EPHA6,",
    "EPHA7,EPHA8,EPHB1,EPHB2,EPHB3,FEZ1,FEZ2,FGF13,FYN,FZD3,GAP43,GDNF,GSK3B,",
    "ISL1,ISL2,KALRN,KIF21A,KIF5C,KLF7,L1CAM,LHX1,LIMK1,LRRC4C,MACF1,MAP1B,MAP2,",
    "MAPT,MYCBP2,NCAM1,NEFH,NEO1,NFASC,NRCAM,NRP1,NRP2,NRXN1,NTN1,NTN4,NTNG1,",
    "NTRK1,NTRK2,PAFAH1B1,PAK1,PLXNA3,PLXNA4,PLXNB1,POU4F1,POU4F2,PTEN,PTK2,",
    "PTPRO,RELN,RET,ROBO1,ROBO2,ROBO3,RTN4,RTN4R,SEMA3A,SEMA3C,SEMA3E,SEMA3F,",
    "SEMA4D,SEMA5A,SEMA6A,SEMA6D,SHH,SLIT1,SLIT2,SLIT3,SLITRK1,SPAST,SPTBN4,",
    "TENM1,TENM2,TENM3,TENM4,TUBB3,ULK1,UNC5A,UNC5B,UNC5C,UNC5D,WNT5A,WNT7A"
  ), ",", fixed = TRUE)[[1]]
)
go_fallback <- lapply(go_fallback, function(x) unique(trimws(x)))

marker_sets <- list(
  Schwann_core = c("S100B", "SOX10", "MPZ", "MBP", "PMP22", "PLP1", "MAG", "PRX", "NGFR"),
  Schwann_extended = c(
    "GFAP", "NCAM1", "L1CAM", "CDH19", "ERBB3", "NRG1", "GAP43",
    "SCN7A", "EGR2", "POU3F1", "MAL", "GJB1", "NES"
  ),
  Neurotrophin_ligands = c(
    "NGF", "BDNF", "NTF3", "NTF4", "GDNF", "NRTN", "ARTN", "PSPN", "CNTF"
  ),
  Neurotrophin_receptors = c(
    "NTRK1", "NTRK2", "NTRK3", "NGFR", "RET",
    "GFRA1", "GFRA2", "GFRA3", "GFRA4", "CNTFR"
  ),
  Neurofilament = c("NEFL", "NEFM", "NEFH", "INA", "PRPH"),
  Neuron_structural = c("TUBB3", "MAP2", "MAPT", "GAP43", "UCHL1", "SYP"),
  Sympathetic = c("TH", "DBH", "SLC6A2", "PNMT"),
  Sensory = c("TAC1", "CALCA", "CALCB", "SCN9A", "SCN10A", "TRPV1", "NTRK1")
)

scheme_catalog <- list(
  list(
    id = "01_pathology",
    title = "1 Pathology PNI",
    type = "pathology",
    fill_name = "PNI-positive\nrate",
    features = list(list(id = "Pathology_PNI", name = "Pathology perineural invasion"))
  ),
  list(
    id = "02_neural_GO",
    title = "2 Neural-signal GO",
    type = "go",
    fill_name = "Pathway score\nmedian",
    features = lapply(names(go_name_map), function(g) {
      list(id = g, name = unname(go_name_map[[g]]))
    })
  ),
  list(
    id = "03_schwann",
    title = "3 Schwann-cell markers",
    type = "marker",
    fill_name = "Marker score\nmedian",
    features = list(
      list(id = "Schwann_core", name = "Schwann core", genes = marker_sets$Schwann_core),
      list(id = "Schwann_extended", name = "Schwann extended", genes = marker_sets$Schwann_extended)
    )
  ),
  list(
    id = "04_neurotrophin",
    title = "4 Neurotrophin ligands and receptors",
    type = "marker",
    fill_name = "Marker score\nmedian",
    features = list(
      list(id = "Neurotrophin_ligands", name = "Neurotrophin ligands", genes = marker_sets$Neurotrophin_ligands),
      list(id = "Neurotrophin_receptors", name = "Neurotrophin receptors", genes = marker_sets$Neurotrophin_receptors)
    )
  ),
  list(
    id = "05_axon_neuron",
    title = "5 Axon and neuron structural genes",
    type = "marker",
    fill_name = "Marker score\nmedian",
    features = list(
      list(id = "Neurofilament", name = "Neurofilament", genes = marker_sets$Neurofilament),
      list(id = "Neuron_structural", name = "Neuron structural", genes = marker_sets$Neuron_structural)
    )
  ),
  list(
    id = "06_autonomic_sensory",
    title = "6 Sympathetic vs sensory",
    type = "marker",
    fill_name = "Marker score\nmedian",
    features = list(
      list(id = "Sympathetic", name = "Sympathetic", genes = marker_sets$Sympathetic),
      list(id = "Sensory", name = "Sensory", genes = marker_sets$Sensory)
    )
  )
)

# ==============================================================================
# 工具函数（全部先定义，读完表达矩阵再分析）
# ==============================================================================
normalize_barcode <- function(x) {
  x <- toupper(gsub("\\.", "-", as.character(x)))
  x <- sub("A$", "", x)
  ifelse(nchar(x) >= 15, substr(x, 1, 15), x)
}
patient_id <- function(x) {
  x <- normalize_barcode(x)
  ifelse(nchar(x) >= 12, substr(x, 1, 12), x)
}
sample_type_code <- function(x) {
  x <- normalize_barcode(x)
  ifelse(nchar(x) >= 15, substr(x, 14, 15), x)
}
first_present <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) NA_character_ else hit[1]
}
safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

existing_files <- function(stems) {
  cands <- unique(unlist(lapply(stems, function(s) c(s, paste0(s, ".gz")))))
  cands[file.exists(cands) & !is.na(file.info(cands)$isdir) & !file.info(cands)$isdir]
}
probe_table_file <- function(f) {
  info <- file.info(f)
  if (is.null(info) || is.na(info$size) || info$size < 200) {
    return(data.table(file = f, size = 0, n_col = 0L, readable = FALSE))
  }
  hdr <- tryCatch(names(fread(f, nrows = 0, fill = TRUE, showProgress = FALSE)),
                  error = function(e) character())
  data.table(file = f, size = as.numeric(info$size), n_col = length(hdr),
             readable = length(hdr) >= 1L)
}
pick_best_table <- function(stems, must = TRUE, min_cols = 1L, label = "table") {
  hits <- existing_files(stems)
  if (length(hits) == 0) {
    if (must) stop("找不到", label, "：", paste(stems, collapse = " / "))
    return(NA_character_)
  }
  tab <- rbindlist(lapply(hits, probe_table_file), fill = TRUE)
  message(label, "候选：\n",
          paste(sprintf("  %s  大小=%.1fMB  列=%d  可读=%s",
                        tab$file, tab$size / 1024^2, tab$n_col, tab$readable),
                collapse = "\n"))
  ok <- tab[readable == TRUE & n_col >= min_cols]
  if (nrow(ok) == 0) {
    if (must) stop(label, "存在但读不了。请删掉空的 .tsv，保留 .tsv.gz 后再 Source。")
    return(NA_character_)
  }
  setorder(ok, -n_col, -size)
  ok$file[1]
}
safe_fread <- function(path, ...) {
  out <- tryCatch(fread(path, showProgress = TRUE, ...), error = function(e) e)
  if (!inherits(out, "error")) return(out)
  alt <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    sub("\\.gz$", "", path, ignore.case = TRUE)
  } else {
    paste0(path, ".gz")
  }
  if (file.exists(alt) && isTRUE(file.info(alt)$size > 200)) {
    message("读 ", path, " 失败（", conditionMessage(out), "），改读 ", alt)
    return(fread(alt, showProgress = TRUE, ...))
  }
  stop("读不了 ", path, "：", conditionMessage(out))
}
pick_best_clinical <- function() {
  f <- pick_best_table(
    c("TCGA-BRCA.clinical.tsv", "TCGA-BRCA.GDC_phenotype.tsv"),
    must = TRUE, min_cols = 5L, label = "临床表"
  )
  hdr <- names(fread(f, nrows = 0, fill = TRUE, showProgress = FALSE))
  message("选用临床：", f, "  列=", length(hdr))
  f
}
pick_clin_col <- function(dt, patterns) {
  nms <- names(dt)
  skip <- grepl(
    "^(sample|patient|barcode|submitter|case_id|id|project|age_|days_|year_|uuid)",
    nms, ignore.case = TRUE
  )
  nms2 <- nms[!skip]
  for (p in patterns) {
    hit <- grep(p, nms2, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

simplify_stage <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("IV|STAGE.?4", x)] <- "Stage IV"
  out[is.na(out) & grepl("III|STAGE.?3", x)] <- "Stage III"
  out[is.na(out) & grepl("II|STAGE.?2", x)] <- "Stage II"
  out[is.na(out) & grepl("I|STAGE.?1", x)] <- "Stage I"
  out[grepl("\\bX\\b|NOT AVAILABLE|UNKNOWN|NOT REPORTED|NOT APPLICABLE", x)] <- NA_character_
  out
}
classify_m <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("M1[ABC]?|\\bM1\\b", x)] <- "M1"
  out[is.na(out) & grepl("\\bM0\\b|CM0", x)] <- "M0"
  out[grepl("MX|NOT AVAILABLE|UNKNOWN|NOT REPORTED", x)] <- NA_character_
  out
}
classify_n <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("N[1-3]", x)] <- "Nplus"
  out[is.na(out) & grepl("\\bN0\\b", x)] <- "N0"
  out[grepl("NX|NOT AVAILABLE|UNKNOWN|NOT REPORTED", x)] <- NA_character_
  out
}
classify_sample_type <- function(x, barcode) {
  out <- rep(NA_character_, length(x))
  xl <- toupper(as.character(x))
  out[grepl("PRIMARY", xl)] <- "PrimaryTumor"
  out[grepl("METASTATIC", xl)] <- "MetastaticTissue"
  out[grepl("NORMAL", xl)] <- "Normal"
  code <- sample_type_code(barcode)
  out[is.na(out) & code == "01"] <- "PrimaryTumor"
  out[is.na(out) & code %in% c("06", "07")] <- "MetastaticTissue"
  out[is.na(out) & code == "11"] <- "Normal"
  out
}
classify_pni <- function(x) {
  x <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[grepl("^(YES|Y|TRUE|PRESENT|POSITIVE|POS|1)\\b", x)] <- "PNI_pos"
  out[grepl("^(NO|N|FALSE|ABSENT|NEGATIVE|NEG|0)\\b", x)] <- "PNI_neg"
  out[grepl("NOT AVAILABLE|UNKNOWN|NOT REPORTED|NOT APPLICABLE|INDETERMINATE", x)] <- NA_character_
  out
}

find_pni_column <- function(clin) {
  nms <- names(clin)
  pat <- "perineural|peri[_ .]?neural|neural[_ .]?invasion|nerve[_ .]?invasion|(^|_)pni(_|$)"
  hit <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
  if (length(hit) == 0) {
    return(list(available = FALSE, column = NA_character_, n_pos = 0L, n_neg = 0L))
  }
  col <- hit[1]
  lab <- classify_pni(clin[[col]])
  list(
    available = TRUE, column = col,
    n_pos = sum(lab == "PNI_pos", na.rm = TRUE),
    n_neg = sum(lab == "PNI_neg", na.rm = TRUE),
    labels = lab
  )
}

# 通路/基因集活性：基因 z 后对样本取均值。不要用 scale()/t()，变量不要叫 score
pathway_zmean <- function(expr_mat, genes) {
  genes <- unique(intersect(as.character(genes), rownames(expr_mat)))
  if (length(genes) < min_set_genes) return(NULL)
  sub <- as.matrix(expr_mat[genes, , drop = FALSE])
  storage.mode(sub) <- "double"
  gene_mean <- rowMeans(sub, na.rm = TRUE)
  gene_sd <- sqrt(rowMeans((sub - gene_mean)^2, na.rm = TRUE))
  gene_sd[!is.finite(gene_sd) | gene_sd < 1e-12] <- 1
  z <- (sub - gene_mean) / gene_sd
  z[!is.finite(z)] <- 0
  set_score <- colMeans(z, na.rm = TRUE)
  names(set_score) <- colnames(sub)
  attr(set_score, "n_genes") <- length(genes)
  attr(set_score, "genes") <- genes
  set_score
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
    if (length(syms) > 0) {
      return(list(genes = syms, source = "org.Hs.eg.db"))
    }
  }
  fb <- go_fallback[[go_id]]
  if (is.null(fb) || length(fb) == 0) {
    return(list(genes = character(), source = "none"))
  }
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
    pvalue = pval
  )
}

as_symbol_matrix <- function(fpkm_data, probe_annot) {
  fpkm_data <- as.data.table(fpkm_data)
  id_col <- names(fpkm_data)[1]
  rid <- as.character(fpkm_data[[id_col]])
  ens <- sub("\\..*$", "", rid)
  mat <- as.matrix(fpkm_data[, -1, with = FALSE])
  storage.mode(mat) <- "double"
  colnames(mat) <- normalize_barcode(colnames(mat))
  if (anyDuplicated(colnames(mat))) {
    mat <- mat[, !duplicated(colnames(mat)), drop = FALSE]
  }
  probe_annot <- as.data.table(probe_annot)
  pid <- first_present(names(probe_annot), c("id", "ensembl", "Ensembl", "ENSEMBL", names(probe_annot)[1]))
  psym <- first_present(names(probe_annot), c("gene", "symbol", "Gene", "SYMBOL", names(probe_annot)[2]))
  map_ens <- sub("\\..*$", "", as.character(probe_annot[[pid]]))
  map_sym <- as.character(probe_annot[[psym]])
  if (mean(grepl("^ENSG", ens, ignore.case = TRUE), na.rm = TRUE) >= 0.5) {
    sym <- map_sym[match(ens, map_ens)]
    keep <- !is.na(sym) & sym != "" & !grepl("^XLOC_", sym)
    mat <- mat[keep, , drop = FALSE]
    sym <- sym[keep]
  } else {
    sym <- rid
  }
  if (anyDuplicated(sym)) {
    dt <- data.table(symbol = sym, as.data.table(mat))
    dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = symbol]
    mat <- as.matrix(dt[, -1, with = FALSE])
    rownames(mat) <- dt$symbol
  } else {
    rownames(mat) <- sym
  }
  mat
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
      scheme_id, scheme_title, feature_id, feature_name, grouping, panel,
      side = "Non-metastatic", x_lab = neg_lab,
      n = n_neg, median_value = median_neg, pvalue
    )],
    d[, .(
      scheme_id, scheme_title, feature_id, feature_name, grouping, panel,
      side = "Metastatic", x_lab = pos_lab,
      n = n_pos, median_value = median_pos, pvalue
    )]
  ), fill = TRUE)
}

plot_bubble_two_cols <- function(stat_dt, title, subtitle, path_stub,
                                 fill_name = "Score\nmedian", facet = FALSE) {
  long <- expand_two_cols(stat_dt)
  long <- long[is.finite(median_value)]
  if (nrow(long) == 0) return(invisible(NULL))
  feat_lv <- unique(stat_dt$feature_name)
  long[, y_lab := factor(feature_name, levels = rev(feat_lv))]
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
      name = fill_name
    ) +
    scale_size_continuous(range = c(3, 11), name = expression(-log[10](p))) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1, size = 11),
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      strip.text = element_text(size = 10)
    )
  if (isTRUE(facet) && "panel" %in% names(long) && uniqueN(long$panel) > 1) {
    p <- p + facet_wrap(~ panel, nrow = 1, scales = "free_x")
  }
  n_panel <- if (isTRUE(facet)) max(1, uniqueN(long$panel)) else 1
  save_plot(p, path_stub, max(8, 4.2 * n_panel + 4), max(6, 0.42 * uniqueN(long$y_lab) + 2.6))
}

build_annotation <- function(sample_ids, clin) {
  ann <- data.table(sample = normalize_barcode(sample_ids))
  ann <- ann[!is.na(sample) & sample != ""]
  ann <- ann[!duplicated(sample)]
  ann[, patient := patient_id(sample)]
  ann[, barcode_type := sample_type_code(sample)]

  clin <- as.data.table(clin)
  idc <- first_present(names(clin), c(
    "sample", "sampleID", "sample_id", "submitter_id.samples",
    "bcr_sample_barcode", "bcr_patient_barcode", "submitter_id", names(clin)[1]
  ))
  clin[, sample_std := normalize_barcode(clin[[idc]])]
  clin[, patient_std := patient_id(clin[[idc]])]

  st <- pick_clin_col(clin, c(
    "ajcc_pathologic_stage", "ajcc_pathologic_tumor_stage",
    "pathologic_stage", "tumor_stage", "clinical_stage"
  ))
  mc <- pick_clin_col(clin, c(
    "ajcc_pathologic_m", "ajcc_metastasis_pathologic_pm", "pathologic_m", "clinical_m"
  ))
  nc <- pick_clin_col(clin, c(
    "ajcc_pathologic_n", "ajcc_nodes_pathologic_pn", "pathologic_n", "clinical_n"
  ))
  tc <- pick_clin_col(clin, c("sample_type\\.samples", "sample_type", "tumor_descriptor"))
  pni_info <- find_pni_column(clin)
  message(
    "选用临床列：stage=", st, "  M=", mc, "  N=", nc,
    "  sample_type=", tc, "  PNI=", pni_info$column
  )

  keep_cols <- unique(c("sample_std", "patient_std", na.omit(c(st, mc, nc, tc, pni_info$column))))
  extra <- clin[, keep_cols, with = FALSE]
  extra_s <- extra[!duplicated(sample_std)]
  ann <- merge(ann, extra_s, by.x = "sample", by.y = "sample_std", all.x = TRUE)
  extra_p <- extra[!duplicated(patient_std)]
  for (nm in setdiff(names(extra_p), c("sample_std", "patient_std", names(ann)))) {
    src <- extra_p[[nm]][match(ann$patient, extra_p$patient_std)]
    cur <- if (nm %in% names(ann)) ann[[nm]] else rep(NA, nrow(ann))
    miss <- is.na(cur) | as.character(cur) %in% c("", "NA")
    cur[miss] <- src[miss]
    ann[, (nm) := cur]
  }

  n <- nrow(ann)
  st_vec <- if (!is.na(st) && st %in% names(ann)) simplify_stage(ann[[st]]) else rep(NA_character_, n)
  m_vec  <- if (!is.na(mc) && mc %in% names(ann)) classify_m(ann[[mc]]) else rep(NA_character_, n)
  n_vec  <- if (!is.na(nc) && nc %in% names(ann)) classify_n(ann[[nc]]) else rep(NA_character_, n)
  raw_type <- if (!is.na(tc) && tc %in% names(ann)) ann[[tc]] else rep(NA_character_, n)
  type_vec <- classify_sample_type(raw_type, ann$sample)
  stage_iv <- rep(NA_character_, n)
  stage_iv[!is.na(st_vec) & st_vec == "Stage IV"] <- "Stage IV"
  stage_iv[!is.na(st_vec) & st_vec %in% c("Stage I", "Stage II", "Stage III")] <- "Stage I-III"
  pni_vec <- if (isTRUE(pni_info$available) && !is.na(pni_info$column) && pni_info$column %in% names(ann)) {
    classify_pni(ann[[pni_info$column]])
  } else {
    rep(NA_character_, n)
  }

  ann[, stage_simplified := st_vec]
  ann[, distant_M := factor(m_vec, levels = c("M0", "M1"))]
  ann[, node_N := factor(n_vec, levels = c("N0", "Nplus"))]
  ann[, stage_IV := factor(stage_iv, levels = c("Stage I-III", "Stage IV"))]
  ann[, sample_class := factor(type_vec, levels = c("PrimaryTumor", "MetastaticTissue", "Normal"))]
  ann[, pni_status := factor(pni_vec, levels = c("PNI_neg", "PNI_pos"))]
  attr(ann, "pni_info") <- pni_info
  ann
}

# ==============================================================================
# 读数据（分析函数在 expr 建好之后才调用）
# ==============================================================================
if (dir.exists(nerve_work_dir)) {
  setwd(nerve_work_dir)
} else {
  message("未找到 ", nerve_work_dir, " ，改用当前目录：", getwd())
}

fpkm_file <- pick_best_table(c("TCGA-BRCA.star_fpkm.tsv"), must = TRUE, min_cols = 20L, label = "FPKM")
clin_file <- pick_best_clinical()
probe_file <- pick_best_table(c(
  "gencode.v36.annotation.gtf.gene.probemap",
  "gencode.v36.annotation.gtf.gene.probemap.tsv"
), must = TRUE, min_cols = 2L, label = "基因注释")
surv_file <- pick_best_table(c("TCGA-BRCA.survival.tsv"), must = FALSE, min_cols = 2L, label = "生存表")

message("FPKM：", fpkm_file)
message("临床：", clin_file)
message("注释：", probe_file)
if (!is.na(surv_file)) message("生存：", surv_file)

fpkm_data <- safe_fread(fpkm_file)
probe_annot <- safe_fread(probe_file)
clinical_data <- safe_fread(clin_file)

dir.create(nerve_out_dir, showWarnings = FALSE, recursive = TRUE)

message("正在把 FPKM 映射为基因符号矩阵…")
nerve_expr_all <- as_symbol_matrix(fpkm_data, probe_annot)
mx <- suppressWarnings(max(nerve_expr_all, na.rm = TRUE))
if (is.finite(mx) && mx > 50) {
  message("检测到原始 FPKM（max=", round(mx, 2), "），做 log2(x+1)")
  nerve_expr_all <- log2(nerve_expr_all + 1)
} else {
  message("表达值范围较小（max=", round(mx, 2), "），视为已 log 转换")
}

nerve_ann <- build_annotation(colnames(nerve_expr_all), clinical_data)
nerve_ann <- nerve_ann[sample %in% colnames(nerve_expr_all)]
fwrite(nerve_ann, file.path(nerve_out_dir, "00_sample_annotation.csv"))
pni_info <- attr(nerve_ann, "pni_info")
if (is.null(pni_info)) pni_info <- list(available = FALSE, column = NA_character_, n_pos = 0L, n_neg = 0L)

message(
  "样本：原位=", sum(nerve_ann$sample_class == "PrimaryTumor", na.rm = TRUE),
  "  转移组织=", sum(nerve_ann$sample_class == "MetastaticTissue", na.rm = TRUE),
  "  正常=", sum(nerve_ann$sample_class == "Normal", na.rm = TRUE),
  "  M1=", sum(nerve_ann$distant_M == "M1", na.rm = TRUE),
  "  Stage IV=", sum(nerve_ann$stage_simplified == "Stage IV", na.rm = TRUE),
  "  N+=", sum(nerve_ann$node_N == "Nplus", na.rm = TRUE),
  "  PNI+=", sum(nerve_ann$pni_status == "PNI_pos", na.rm = TRUE)
)

primary_ids <- nerve_ann$sample[nerve_ann$sample_class == "PrimaryTumor"]
met_ids <- nerve_ann$sample[nerve_ann$sample_class == "MetastaticTissue"]
tumor_ids <- unique(c(primary_ids, met_ids))
if (length(primary_ids) < 20) stop("原位肿瘤样本太少：", length(primary_ids))

nerve_expr_tumor <- nerve_expr_all[, intersect(tumor_ids, colnames(nerve_expr_all)), drop = FALSE]
keep_g <- rowMeans(is.finite(nerve_expr_tumor) & nerve_expr_tumor > 0, na.rm = TRUE) >= min_expr_frac
nerve_expr_tumor <- nerve_expr_tumor[keep_g, , drop = FALSE]
nerve_expr_primary <- nerve_expr_tumor[, intersect(primary_ids, colnames(nerve_expr_tumor)), drop = FALSE]
message(
  "表达矩阵：原位 ", ncol(nerve_expr_primary), " 样本 x ", nrow(nerve_expr_primary),
  " 基因；肿瘤(含转移组织) ", ncol(nerve_expr_tumor)
)

# ==============================================================================
# 主分析（必须在 expr 建好之后）
# ==============================================================================
run_nerve_tcga <- function() {
  if (!exists("nerve_expr_primary", inherits = TRUE)) {
    stop("还没有 nerve_expr_primary，请从脚本开头 Source")
  }

  score_one_set <- function(expr_mat, genes, label) {
    set_score <- pathway_zmean(expr_mat, genes)
    if (is.null(set_score)) {
      message("  ", label, " 映射基因不足，跳过")
      return(NULL)
    }
    message("  ", label, "  基因数=", attr(set_score, "n_genes"))
    set_score
  }

  mat_from_list <- function(lst) {
    if (length(lst) == 0) return(NULL)
    common <- Reduce(intersect, lapply(lst, names))
    if (length(common) == 0) return(NULL)
    mat <- do.call(cbind, lapply(lst, function(x) x[common]))
    colnames(mat) <- names(lst)
    rownames(mat) <- common
    mat
  }

  metastasis_designs <- function(ann_p, ann_t) {
    list(
      list(key = "a_distant_M", title = "Distant metastasis", panel = "1a Distant M",
           group = setNames(as.character(ann_p$distant_M), ann_p$sample),
           pos = "M1", neg = "M0", pos_lab = "M1", neg_lab = "M0",
           use = "primary"),
      list(key = "b_AJCC_stageIV", title = "AJCC stage", panel = "1b AJCC stage",
           group = setNames(as.character(ann_p$stage_IV), ann_p$sample),
           pos = "Stage IV", neg = "Stage I-III",
           pos_lab = "Stage IV", neg_lab = "Stage I-III",
           use = "primary"),
      list(key = "c_node_N", title = "Lymph node", panel = "1c Lymph node",
           group = setNames(as.character(ann_p$node_N), ann_p$sample),
           pos = "Nplus", neg = "N0", pos_lab = "N+", neg_lab = "N0",
           use = "primary"),
      list(key = "d_sample_type", title = "Sample type", panel = "1d Sample type",
           group = setNames(as.character(ann_t$sample_class), ann_t$sample),
           pos = "MetastaticTissue", neg = "PrimaryTumor",
           pos_lab = "Metastatic tissue", neg_lab = "Primary tumor",
           use = "tumor")
    )
  }

  run_one_scheme <- function(scheme, primary_list, tumor_list, ann_p, ann_t) {
    sdir <- file.path(nerve_out_dir, scheme$id)
    dir.create(sdir, showWarnings = FALSE, recursive = TRUE)
    sm_p <- mat_from_list(primary_list)
    sm_t <- mat_from_list(tumor_list)
    if (is.null(sm_p) && is.null(sm_t)) {
      writeLines("No scorable features in this scheme.", file.path(sdir, "SKIPPED.txt"))
      return(data.table())
    }
    if (!is.null(sm_p)) {
      fwrite(
        data.table(sample = rownames(sm_p), as.data.table(sm_p)),
        file.path(sdir, "scores_primary.csv")
      )
    }
    if (!is.null(sm_t)) {
      fwrite(
        data.table(sample = rownames(sm_t), as.data.table(sm_t)),
        file.path(sdir, "scores_tumor.csv")
      )
    }

    feat_name <- setNames(
      vapply(scheme$features, function(f) f$name, character(1)),
      vapply(scheme$features, function(f) f$id, character(1))
    )
    designs <- metastasis_designs(ann_p, ann_t)
    all_stat <- list()
    for (ds in designs) {
      sm <- if (identical(ds$use, "tumor")) sm_t else sm_p
      if (is.null(sm)) next
      stat_rows <- list()
      for (fid in colnames(sm)) {
        value_vec <- as.numeric(sm[, fid])
        names(value_vec) <- rownames(sm)
        one <- compare_groups(value_vec, ds$group[names(value_vec)], ds$pos, ds$neg, ds$key)
        if (is.null(one)) next
        one[, `:=`(
          scheme_id = scheme$id, scheme_title = scheme$title,
          feature_id = fid,
          feature_name = unname(if (fid %in% names(feat_name)) feat_name[[fid]] else fid),
          panel = ds$panel, pos_lab = ds$pos_lab, neg_lab = ds$neg_lab
        )]
        stat_rows[[fid]] <- one
      }
      stat_dt <- rbindlist(stat_rows, fill = TRUE)
      if (nrow(stat_dt) == 0) {
        message("  分组人数不足，跳过 ", scheme$id, " / ", ds$key)
        next
      }
      stat_dt[, fdr := p.adjust(pvalue, method = "BH")]
      fwrite(stat_dt, file.path(sdir, paste0(ds$key, "_vs_metastasis.csv")))
      fwrite(expand_two_cols(stat_dt), file.path(sdir, paste0(ds$key, "_vs_metastasis_two_cols.csv")))
      all_stat[[ds$key]] <- stat_dt
      plot_bubble_two_cols(
        stat_dt,
        title = paste0(scheme$title, " — ", ds$title),
        subtitle = "Y = neural-invasion feature (scored separately); X = non-metastatic | metastatic; fill = median; size = -log10(Wilcoxon p)",
        path_stub = file.path(sdir, paste0(ds$key, "_bubble")),
        fill_name = scheme$fill_name
      )
    }
    if (length(all_stat) == 0) return(data.table())
    bubble <- rbindlist(all_stat, fill = TRUE)
    fwrite(bubble, file.path(sdir, "summary_vs_metastasis.csv"))
    fwrite(expand_two_cols(bubble), file.path(sdir, "summary_vs_metastasis_two_cols.csv"))
    plot_bubble_two_cols(
      bubble,
      title = paste0(scheme$title, " versus metastasis"),
      subtitle = "Each panel: non-metastatic | metastatic; features within this scheme are not pooled",
      path_stub = file.path(sdir, "summary_bubble"),
      fill_name = scheme$fill_name,
      facet = TRUE
    )
    bubble
  }

  gene_log <- list()
  scheme_rows <- list()
  all_scheme_stats <- list()

  write_readme <- function() {
    txt <- c(
      "nerve_TCGA results",
      "Six neural-invasion schemes are scored separately and never pooled.",
      "Bubble plots: two x-columns (non-metastatic | metastatic); fill = median; size = -log10(p).",
      "1 Pathology PNI: gold-standard label from clinical table if a PNI column exists.",
      "2 Neural-signal GO: GO:0019227, GO:1902847, GO:0097374, GO:1902667, GO:0007409.",
      "3 Schwann-cell markers: core and extended sets (SOX10/S100B also mark some basal/myoepithelial cells).",
      "4 Neurotrophin ligands vs receptors (not mixed).",
      "5 Neurofilament vs neuron-structural genes.",
      "6 Sympathetic (TH/DBH/SLC6A2/PNMT) vs sensory (TAC1/CALCA/SCN9A/NTRK1).",
      paste0("PNI column available: ", isTRUE(pni_info$available),
             if (!is.na(pni_info$column)) paste0(" (", pni_info$column, ")") else "")
    )
    writeLines(txt, file.path(nerve_out_dir, "00_README.txt"))
  }
  write_readme()

  # ---- 1 病理 PNI ----
  s1 <- scheme_catalog[[1]]
  s1dir <- file.path(nerve_out_dir, s1$id)
  dir.create(s1dir, showWarnings = FALSE, recursive = TRUE)
  if (!isTRUE(pni_info$available) || sum(!is.na(nerve_ann$pni_status)) < 4) {
    writeLines(c(
      "Pathology PNI is not available in this Xena/GDC clinical table.",
      "TCGA-BRCA.clinical.tsv has no perineural-invasion column.",
      "This scheme is skipped. Schemes 2-6 still run.",
      paste0("Scanned columns: ", paste(names(clinical_data), collapse = ", "))
    ), file.path(s1dir, "NOT_AVAILABLE.txt"))
    message("方案 1 病理 PNI：临床表无此列，已跳过")
    scheme_rows[[s1$id]] <- data.table(
      scheme_id = s1$id, scheme_title = s1$title, status = "not_available",
      n_features = 0L
    )
  } else {
    message("方案 1 病理 PNI：列=", pni_info$column)
    pni_num <- ifelse(as.character(nerve_ann$pni_status) == "PNI_pos", 1,
                      ifelse(as.character(nerve_ann$pni_status) == "PNI_neg", 0, NA_real_))
    names(pni_num) <- nerve_ann$sample
    pni_p <- pni_num[intersect(names(pni_num), colnames(nerve_expr_primary))]
    pni_t <- pni_num[intersect(names(pni_num), colnames(nerve_expr_tumor))]
    attr(pni_p, "n_genes") <- 0L
    attr(pni_p, "genes") <- "clinical_PNI"
    attr(pni_t, "n_genes") <- 0L
    attr(pni_t, "genes") <- "clinical_PNI"
    gene_log[[s1$id]] <- data.table(
      scheme_id = s1$id, feature_id = "Pathology_PNI",
      feature_name = "Pathology perineural invasion",
      gene_source = paste0("clinical:", pni_info$column),
      n_genes = 0L, genes = pni_info$column
    )
    all_scheme_stats[[s1$id]] <- run_one_scheme(
      s1, list(Pathology_PNI = pni_p), list(Pathology_PNI = pni_t),
      nerve_ann[sample %in% names(pni_p)],
      nerve_ann[sample %in% names(pni_t)]
    )
    scheme_rows[[s1$id]] <- data.table(
      scheme_id = s1$id, scheme_title = s1$title, status = "ok", n_features = 1L
    )
  }

  # ---- 2 神经 GO ----
  s2 <- scheme_catalog[[2]]
  message("方案 2 神经信号 GO（每个 GO 单独，不合并）")
  go_primary <- list()
  go_tumor <- list()
  for (ft in s2$features) {
    gm <- get_go_genes(ft$id)
    gene_log[[ft$id]] <- data.table(
      scheme_id = s2$id, feature_id = ft$id, feature_name = ft$name,
      gene_source = gm$source, n_genes = length(gm$genes),
      genes = paste(gm$genes, collapse = ";")
    )
    go_primary[[ft$id]] <- score_one_set(nerve_expr_primary, gm$genes, paste(ft$id, ft$name))
    go_tumor[[ft$id]] <- score_one_set(nerve_expr_tumor, gm$genes, paste(ft$id, "tumor"))
  }
  go_primary <- Filter(Negate(is.null), go_primary)
  go_tumor <- Filter(Negate(is.null), go_tumor)
  all_scheme_stats[[s2$id]] <- run_one_scheme(
    s2, go_primary, go_tumor,
    nerve_ann[sample %in% colnames(nerve_expr_primary)],
    nerve_ann[sample %in% colnames(nerve_expr_tumor)]
  )
  scheme_rows[[s2$id]] <- data.table(
    scheme_id = s2$id, scheme_title = s2$title, status = "ok",
    n_features = length(go_primary)
  )

  # ---- 3–6 marker 方案 ----
  for (scheme in scheme_catalog[3:6]) {
    message("方案 ", scheme$title, "（各条目单独打分，不合并）")
    primary_list <- list()
    tumor_list <- list()
    for (ft in scheme$features) {
      gene_log[[ft$id]] <- data.table(
        scheme_id = scheme$id, feature_id = ft$id, feature_name = ft$name,
        gene_source = "curated_markers", n_genes = length(ft$genes),
        genes = paste(ft$genes, collapse = ";")
      )
      primary_list[[ft$id]] <- score_one_set(nerve_expr_primary, ft$genes, ft$name)
      tumor_list[[ft$id]] <- score_one_set(nerve_expr_tumor, ft$genes, paste(ft$name, "tumor"))
    }
    primary_list <- Filter(Negate(is.null), primary_list)
    tumor_list <- Filter(Negate(is.null), tumor_list)
    all_scheme_stats[[scheme$id]] <- run_one_scheme(
      scheme, primary_list, tumor_list,
      nerve_ann[sample %in% colnames(nerve_expr_primary)],
      nerve_ann[sample %in% colnames(nerve_expr_tumor)]
    )
    scheme_rows[[scheme$id]] <- data.table(
      scheme_id = scheme$id, scheme_title = scheme$title, status = "ok",
      n_features = length(primary_list)
    )
  }

  gene_dt <- rbindlist(gene_log, fill = TRUE)
  fwrite(gene_dt, file.path(nerve_out_dir, "00_genes_used.csv"))
  fwrite(rbindlist(scheme_rows, fill = TRUE), file.path(nerve_out_dir, "00_scheme_index.csv"))

  keep_stats <- Filter(
    function(x) is.data.table(x) && nrow(x) > 0,
    all_scheme_stats[names(all_scheme_stats) != "01_pathology"]
  )
  if (length(keep_stats) > 0) {
    overview <- rbindlist(keep_stats, fill = TRUE)
    fwrite(overview, file.path(nerve_out_dir, "07_all_schemes_vs_metastasis.csv"))
    fwrite(expand_two_cols(overview), file.path(nerve_out_dir, "07_all_schemes_vs_metastasis_two_cols.csv"))
    plot_bubble_two_cols(
      overview,
      title = "Six neural-invasion schemes versus metastasis",
      subtitle = "Schemes are not pooled; each row is one feature inside one scheme",
      path_stub = file.path(nerve_out_dir, "07_all_schemes_summary_bubble"),
      fill_name = "Score\nmedian",
      facet = TRUE
    )
  }

  message("完成。结果目录：", normalizePath(nerve_out_dir, winslash = "/", mustWork = FALSE))
  message("每类气泡图：01_pathology/ … 06_autonomic_sensory/ 下的 summary_bubble.png")
  message("总览：07_all_schemes_summary_bubble.png")
  invisible(TRUE)
}

if (exists("nerve_expr_primary") && ncol(nerve_expr_primary) > 10) {
  run_nerve_tcga()
} else {
  stop("表达矩阵未建好，请从第一行完整 Source")
}
