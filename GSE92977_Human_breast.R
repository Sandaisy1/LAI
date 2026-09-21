#!/usr/bin/env Rscript
# =============================================================================
# GSE92977_Human_breast.R
# Cejalvo 2017 ConvertHER：配对原发 vs 肺/骨转移 NanoString（GSE92977）
#
# 数据目录（脚本与结果默认同目录）：
#   E:/R/Human breast cancer/GSE92977
#
# 四个问题（只分析「原发相对配对转移更低」= 转移/原发 FC > 1）：
#   1) 配对原发 vs 肺转移、配对原发 vs 骨转移
#   2) 器官特异：只肺或只骨（本队列每人只有一个转移部位）
#   3) 轴突导向 / 施旺细胞 / 神经营养 三套签名（受 ~105 基因面板限制）
#   4) 三种神经浸润分别与肺转移的关系
#
# 阈值：p < 0.05，FC > 1 与 FC > 1.25。不做 top50–300。
# 配对：Patient N 原发 对 Patient N 该器官转移（热图列顺序 Primary_i, Met_i）。
#       以 title「Patient N primary/metastatic tumor」为准（个别 description/individual 有笔误）。
#
# 本数据集：123 对 NanoString（管家基因标准化后约 105 基因），不是全转录组。
#   肺转移 7 对；骨转移 16 对；无人同时有肺和骨。
# 不要改 TG_RNAseq_*.R。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")
options(clusterProfiler.download.method = "auto")

# -----------------------------------------------------------------------------
# 0. 包
# -----------------------------------------------------------------------------
cran_required <- c(
  "dplyr", "tidyr", "tibble", "stringr", "ggplot2", "ggrepel",
  "pheatmap", "RColorBrewer", "matrixStats", "writexl"
)
cran_optional <- c("ggvenn", "GSVA")
bioc_required <- c(
  "limma", "clusterProfiler", "org.Hs.eg.db",
  "enrichplot", "AnnotationDbi", "fgsea"
)
bioc_optional <- c("ReactomePA", "msigdbr", "pathview")

install_if_missing <- function(pkgs, bioc = FALSE, required = TRUE) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) == 0) return(invisible(TRUE))
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    tryCatch(
      BiocManager::install(miss, update = FALSE, ask = FALSE),
      error = function(e) message("Bioconductor install failed: ", e$message)
    )
  } else {
    tryCatch(
      install.packages(miss, repos = "https://cloud.r-project.org"),
      error = function(e) message("CRAN install failed: ", e$message)
    )
  }
  still <- miss[!vapply(miss, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still) > 0 && required) {
    stop("缺少必需 R 包: ", paste(still, collapse = ", "))
  }
  if (length(still) > 0) {
    message("可选包未安装，相关分析将跳过: ", paste(still, collapse = ", "))
  }
  invisible(TRUE)
}

install_if_missing(cran_required, bioc = FALSE, required = TRUE)
install_if_missing(cran_optional, bioc = FALSE, required = FALSE)
install_if_missing(bioc_required, bioc = TRUE, required = TRUE)
install_if_missing(bioc_optional, bioc = TRUE, required = FALSE)

safe_library <- function(pkgs) {
  for (p in pkgs) {
    if (requireNamespace(p, quietly = TRUE)) {
      suppressPackageStartupMessages(library(p, character.only = TRUE))
    }
  }
}
safe_library(c(cran_required, cran_optional, bioc_required, bioc_optional))
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

# -----------------------------------------------------------------------------
# 1. 路径
# -----------------------------------------------------------------------------
resolve_project_dir <- function() {
  env_dir <- Sys.getenv("GSE92977_DIR", unset = "")
  script_dir <- tryCatch({
    of <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(of)) dirname(normalizePath(sub("^--file=", "", of[1]))) else NA_character_
  }, error = function(e) NA_character_)
  candidates <- c(
    env_dir,
    "E:/R/Human breast cancer/GSE92977",
    "E:\\R\\Human breast cancer\\GSE92977",
    script_dir,
    file.path(getwd(), "GSE92977"),
    getwd()
  )
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  hit_file <- function(d) {
    if (!dir.exists(d)) return(FALSE)
    ff <- list.files(d, full.names = FALSE, ignore.case = TRUE)
    any(grepl("GSE92977|raw_data|series_matrix", ff, ignore.case = TRUE))
  }
  for (d in candidates) {
    if (hit_file(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

project_dir <- resolve_project_dir()
result_dir  <- file.path(project_dir, "results_GSE92977_Human_breast")
log_dir     <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("GSE92977_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_msg("Working directory: ", project_dir)
log_msg("Results: ", result_dir)

p_cutoff  <- 0.05
fc_cutoffs <- c("FC_1" = 1, "FC_1.25" = 1.25)

# -----------------------------------------------------------------------------
# 2. 读 GEO 文件
# -----------------------------------------------------------------------------
find_first <- function(patterns) {
  ff <- list.files(project_dir, full.names = TRUE, recursive = TRUE)
  ff <- ff[!grepl("results_GSE92977", ff, ignore.case = TRUE)]
  base <- basename(ff)
  for (pat in patterns) {
    hit <- ff[grepl(pat, base, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

strip_quotes <- function(x) gsub('^\"|\"$', "", x)

# 探针 ID -> 符号（series matrix 行名是数字；raw_data 第二列是 Gene.Name）
builtin_probe_map <- c(
  "1" = "ACTB", "2" = "ACTR3B", "3" = "ADM", "4" = "ANGPTL4", "5" = "ANLN",
  "6" = "AXL", "7" = "BAG1", "8" = "BCL2", "9" = "BIRC5", "10" = "BLVRA",
  "11" = "CAV1", "12" = "CCNB1", "13" = "CCNE1", "14" = "CD24", "15" = "CDC20",
  "16" = "CDC6", "17" = "CDH3", "18" = "CENPF", "19" = "CEP55", "20" = "CLDN3",
  "21" = "CLDN4", "22" = "CLDN7", "23" = "CRYAB", "24" = "CXXC5", "25" = "DDIT4",
  "26" = "DDR1", "27" = "DSP", "28" = "EGFR", "29" = "EMP3", "30" = "EPCAM",
  "31" = "ERBB2", "32" = "ERBB3", "33" = "ERBB4", "34" = "ESR1", "35" = "ESRP1",
  "36" = "EVI2A", "37" = "EXO1", "38" = "F11R", "39" = "FABP5", "40" = "FBN1",
  "41" = "FGFR4", "42" = "FLVCR2", "43" = "FOXA1", "44" = "FOXC1", "45" = "GAL",
  "46" = "GATA3", "47" = "GNG11", "48" = "GPR160", "49" = "GRB7", "50" = "GRHL2",
  "51" = "JUP", "52" = "KIF2C", "53" = "KRT14", "54" = "KRT17", "55" = "KRT19",
  "56" = "KRT5", "57" = "KRT8", "58" = "LEPRE1", "59" = "LHFP", "60" = "MAPT",
  "61" = "MDM2", "62" = "MELK", "63" = "MET", "64" = "MIA", "65" = "MKI67",
  "66" = "MLPH", "67" = "MMP11", "68" = "MPP1", "69" = "MRPL19", "70" = "MYBL2",
  "71" = "MYC", "72" = "NAT1", "73" = "NDC80", "74" = "NDRG1", "75" = "NT5E",
  "76" = "NUF2", "77" = "ORC6L", "78" = "PGR", "79" = "PHGDH", "80" = "PIK3CA",
  "81" = "PLOD1", "82" = "PNP", "83" = "PSMC4", "84" = "PTEN", "85" = "PTTG1",
  "86" = "PVRL3", "87" = "RAB25", "88" = "RAD17", "89" = "RAD50", "90" = "RB1",
  "91" = "RPLP0", "92" = "RRAGD", "93" = "RRM2", "94" = "SF3A1", "95" = "SFRP1",
  "96" = "SH2B3", "97" = "SLC16A3", "98" = "SLC39A6", "99" = "SPINT1", "100" = "SPINT2",
  "101" = "TMEM158", "102" = "TMEM45B", "103" = "TYMS", "104" = "UBE2C", "105" = "UBE2T",
  "106" = "UCHL1", "107" = "VAMP8", "108" = "VEGFA", "109" = "VIM", "110" = "ZEB1"
)
hk_genes <- c("ACTB", "MRPL19", "PSMC4", "RPLP0", "SF3A1")

read_probe_map_from_raw <- function(path) {
  df <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    utils::read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  id_col <- names(df)[1]
  gene_col <- if ("Gene.Name" %in% names(df)) "Gene.Name" else names(df)[2]
  ids <- as.character(df[[id_col]])
  sym <- as.character(df[[gene_col]])
  keep <- nzchar(ids) & nzchar(sym) & !is.na(ids) & !is.na(sym)
  stats::setNames(sym[keep], ids[keep])
}

canonical_site <- function(x) {
  u <- toupper(trimws(as.character(x)))
  u[is.na(u)] <- ""
  out <- rep("Other", length(u))
  out[u == ""] <- NA_character_
  out[grepl("LUNG", u)] <- "Lung"
  out[grepl("BONE", u)] <- "Bone"
  out[grepl("PLEURA", u)] <- "Pleura"
  out[grepl("LIVER", u)] <- "Liver"
  out[grepl("BRAIN", u)] <- "Brain"
  out[grepl("LYMPH", u)] <- "LN"
  out[grepl("SKIN", u)] <- "Skin"
  out[grepl("OVAR", u)] <- "Ovary"
  out[grepl("^BREAST$", u)] <- "BreastMet"
  out[grepl("PERICARD", u)] <- "Pericardium"
  out[grepl("ENDOMET", u)] <- "Endometrium"
  out[grepl("BOWEL", u)] <- "Bowel"
  out[grepl("PERITON", u)] <- "Peritoneum"
  out
}

# 读 series matrix：管家基因标准化表达 + 表型
read_gse92977 <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  titles <- gsms <- desc <- NULL
  chars <- list()
  in_table <- FALSE
  header <- NULL
  gene <- character(0)
  rows <- list()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "!series_matrix_table_begin")) {
      in_table <- TRUE
      next
    }
    if (startsWith(line, "!series_matrix_table_end")) break
    if (in_table) {
      parts <- strip_quotes(strsplit(line, "\t", fixed = TRUE)[[1]])
      if (is.null(header)) {
        header <- parts
      } else {
        gene <- c(gene, parts[1])
        num <- suppressWarnings(as.numeric(parts[-1]))
        rows[[length(rows) + 1]] <- num
      }
      next
    }
    if (startsWith(line, "!Sample_title\t")) {
      titles <- strip_quotes(strsplit(line, "\t", fixed = TRUE)[[1]][-1])
    } else if (startsWith(line, "!Sample_geo_accession\t")) {
      gsms <- strip_quotes(strsplit(line, "\t", fixed = TRUE)[[1]][-1])
    } else if (startsWith(line, "!Sample_description\t")) {
      desc <- strip_quotes(strsplit(line, "\t", fixed = TRUE)[[1]][-1])
    } else if (startsWith(line, "!Sample_characteristics_ch1\t")) {
      vals <- strip_quotes(strsplit(line, "\t", fixed = TRUE)[[1]][-1])
      field <- sub(":.*$", "", vals[1])
      chars[[field]] <- sub("^[^:]+:\\s*", "", vals)
    }
  }
  if (is.null(titles) || is.null(gsms)) stop("series matrix 缺 sample 注释: ", path)
  if (length(rows) == 0) stop("series matrix 没有表达表: ", path)
  mat <- do.call(rbind, rows)
  rownames(mat) <- gene
  colnames(mat) <- header[-1]
  n <- length(titles)
  pick <- function(nm, default = NA_character_) {
    if (nm %in% names(chars)) chars[[nm]] else rep(default, n)
  }
  patient <- sub("(?i)^.*patient\\s+([0-9]+).*$", "\\1", titles, perl = TRUE)
  unmatched <- !grepl("(?i)patient\\s+[0-9]+", titles, perl = TRUE)
  if (any(unmatched)) patient[unmatched] <- paste0("P", which(unmatched))
  tissue_raw <- pick("tissue")
  kind <- ifelse(grepl("primary", tissue_raw, ignore.case = TRUE), "primary",
                 ifelse(grepl("metastat", tissue_raw, ignore.case = TRUE), "met", "other"))
  site <- canonical_site(pick("metastatic site"))
  site[kind == "primary"] <- "Breast"
  pheno <- data.frame(
    sample = gsms,
    gsm = gsms,
    title = titles,
    barcode = if (is.null(desc)) rep(NA_character_, n) else desc,
    patient = patient,
    kind = kind,
    tissue = site,
    met_site_raw = pick("metastatic site"),
    subtype = pick("intrinsic subtype breast cancer"),
    age = pick("age_dx metastatic tumor"),
    stringsAsFactors = FALSE
  )
  pheno <- pheno[match(colnames(mat), pheno$sample), ]
  rownames(pheno) <- pheno$sample
  list(mat = mat, pheno = pheno)
}

pheno_path <- find_first(c("GSE92977_series_matrix\\.txt(\\.gz)?$", "series_matrix"))
raw_path <- find_first(c("GSE92977_raw_data\\.txt(\\.gz)?$", "raw_data"))
if (is.na(pheno_path) || !file.exists(pheno_path)) {
  stop("找不到 GSE92977_series_matrix.txt.gz，请放在: ", project_dir,
       "\n同时建议保留 GSE92977_raw_data.txt.gz（探针 ID 到基因名）。")
}

obj <- read_gse92977(pheno_path)
raw <- obj$mat
pheno <- obj$pheno
probe_map <- builtin_probe_map
if (!is.na(raw_path) && file.exists(raw_path)) {
  pm <- tryCatch(read_probe_map_from_raw(raw_path), error = function(e) NULL)
  if (!is.null(pm) && length(pm) > 0) {
    probe_map <- pm
    log_msg("Probe map from ", basename(raw_path), ": ", length(probe_map), " probes")
  }
} else {
  log_msg("WARNING: 未找到 raw_data，使用内置探针-基因名对照")
}
sym <- unname(probe_map[rownames(raw)])
sym[is.na(sym) | !nzchar(sym)] <- paste0("PROBE_", rownames(raw)[is.na(sym) | !nzchar(sym)])
if (any(duplicated(sym))) {
  means <- rowMeans(raw, na.rm = TRUE)
  ord <- order(means, decreasing = TRUE)
  raw <- raw[ord, , drop = FALSE]
  sym <- sym[ord]
  keep <- !duplicated(sym)
  raw <- raw[keep, , drop = FALSE]
  sym <- sym[keep]
}
rownames(raw) <- sym
drop_hk <- rownames(raw) %in% hk_genes
if (any(drop_hk)) {
  log_msg("Drop housekeeping probes still in matrix: ", paste(rownames(raw)[drop_hk], collapse = ", "))
  raw <- raw[!drop_hk, , drop = FALSE]
}

storage.mode(raw) <- "double"
raw[!is.finite(raw)] <- NA
log_msg("NanoString matrix: ", nrow(raw), " genes x ", ncol(raw), " samples from ", basename(pheno_path))
log_msg("Kind: ", paste(names(table(pheno$kind)), table(pheno$kind), sep = "=", collapse = "; "))
log_msg("Met sites: ", paste(names(table(pheno$tissue[pheno$kind == "met"])),
                             table(pheno$tissue[pheno$kind == "met"]), sep = "=", collapse = "; "))

mx <- stats::quantile(as.numeric(raw), 0.99, na.rm = TRUE)
if (is.finite(mx) && mx > 50) {
  log_msg("Values look like raw counts (q99=", signif(mx, 3), "); log2(x+1)")
  logmat <- log2(pmax(raw, 0) + 1)
} else {
  log_msg("Values look log-scale / HK-normalized (q99=", signif(mx, 3), "); use as-is")
  logmat <- raw
}
keep_gene <- rowSums(is.finite(logmat)) >= 3 &
  matrixStats::rowSds(logmat, na.rm = TRUE) > 1e-6
log_msg("Low-variance filter: keep ", sum(keep_gene), " / ", nrow(logmat))
logmat <- logmat[keep_gene, , drop = FALSE]
pheno_use <- pheno
heat_mat <- t(scale(t(logmat)))
heat_mat[!is.finite(heat_mat)] <- 0
utils::write.csv(pheno_use, file.path(result_dir, "00_logs", "sample_annotation.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. 一一对应：Patient N 原发 vs Patient N 该器官转移
# -----------------------------------------------------------------------------
pick_primary <- function(rows) {
  if (nrow(rows) == 1) return(rows$sample[1])
  rows$sample[1]
}

build_pairs <- function(organ, organ_labels) {
  mets <- pheno_use[pheno_use$kind == "met" & pheno_use$tissue %in% organ_labels, , drop = FALSE]
  prim <- pheno_use[pheno_use$kind == "primary", , drop = FALSE]
  out <- list()
  for (pid in intersect(unique(mets$patient), unique(prim$patient))) {
    p_rows <- prim[prim$patient == pid, , drop = FALSE]
    m_rows <- mets[mets$patient == pid, , drop = FALSE]
    p_id <- pick_primary(p_rows)
    m_ids <- m_rows$sample
    out[[length(out) + 1]] <- data.frame(
      pair_id = pid,
      patient = pid,
      organ = organ,
      primary = p_id,
      metastasis = if (length(m_ids) == 1) m_ids else paste(m_ids, collapse = ";"),
      n_met_samples = length(m_ids),
      primary_title = p_rows$title[match(p_id, p_rows$sample)],
      met_title = paste(m_rows$title, collapse = ";"),
      met_tissue = paste(unique(m_rows$tissue), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) {
    return(data.frame(
      pair_id = character(), patient = character(), organ = character(),
      primary = character(), metastasis = character(), n_met_samples = integer(),
      primary_title = character(), met_title = character(),
      met_tissue = character(), stringsAsFactors = FALSE
    ))
  }
  dplyr::bind_rows(out)
}

pairs_lung <- build_pairs("Lung", "Lung")
pairs_bone <- build_pairs("Bone", c("Bone", "Rib", "Spine", "Spinal", "Skull"))
pair_dir <- file.path(result_dir, "00_sample_pairing")
dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(pairs_lung, file.path(pair_dir, "pairs_primary_vs_lung_1to1.csv"), row.names = FALSE)
utils::write.csv(pairs_bone, file.path(pair_dir, "pairs_primary_vs_bone_1to1.csv"), row.names = FALSE)
log_msg("Paired primary-lung: ", nrow(pairs_lung), " patients")
log_msg("Paired primary-bone: ", nrow(pairs_bone), " patients")
if (nrow(pairs_lung) > 0) log_msg("Lung pairs: ", paste(pairs_lung$patient, collapse = ", "))
if (nrow(pairs_bone) > 0) log_msg("Bone pairs: ", paste(pairs_bone$patient, collapse = ", "))
writeLines(
  c("Pairing uses title 'Patient N primary/metastatic tumor', not the description barcode.",
    "GEO quirks: Patient 81 descriptions are 8IP/8IM; Patient 119 both labeled 119M;",
    "Patient 96M individual field wrongly says Patient 95. Title+tissue are the source of truth.",
    "Each patient has exactly one metastatic site in this series."),
  file.path(pair_dir, "NOTE_pairing_rules.txt")
)

# 每个配对患者：原发向量 vs 该器官转移（多灶则均值）
pair_vectors <- function(pairs) {
  prim_mat <- matrix(NA_real_, nrow = nrow(logmat), ncol = nrow(pairs),
                     dimnames = list(rownames(logmat), pairs$patient))
  met_mat <- prim_mat
  for (i in seq_len(nrow(pairs))) {
    p <- pairs$primary[i]
    m <- strsplit(pairs$metastasis[i], ";", fixed = TRUE)[[1]]
    m <- intersect(m, colnames(logmat))
    prim_mat[, i] <- logmat[, p]
    met_mat[, i] <- if (length(m) == 1) logmat[, m] else rowMeans(logmat[, m, drop = FALSE])
  }
  list(primary = prim_mat, met = met_mat, patients = pairs$patient)
}

# -----------------------------------------------------------------------------
# 5. 配对 limma：转移 vs 原发；FC = 转移/原发（>1 表示原发更低）
# -----------------------------------------------------------------------------
paired_limma <- function(pairs, label) {
  if (nrow(pairs) < 2) {
    log_msg(label, ": fewer than 2 pairs, skip limma")
    return(NULL)
  }
  if (nrow(pairs) < 3) {
    log_msg("WARNING ", label, ": only ", nrow(pairs),
            " pairs. p-values have very few residual df; interpret cautiously.")
  }
  vec <- pair_vectors(pairs)
  expr <- cbind(vec$primary, vec$met)
  colnames(expr) <- c(paste0(vec$patients, "_Primary"), paste0(vec$patients, "_Met"))
  patient <- factor(c(vec$patients, vec$patients))
  tissue  <- factor(c(rep("Primary", ncol(vec$primary)), rep("Met", ncol(vec$met))),
                    levels = c("Primary", "Met"))
  design <- stats::model.matrix(~ patient + tissue)
  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit, coef = "tissueMet", number = Inf, sort.by = "none")
  fc_pairs <- 2^(vec$met - vec$primary)
  colnames(fc_pairs) <- paste0("FC_", vec$patients)
  de <- data.frame(
    gene = rownames(tt),
    log2FC = tt$logFC,
    FC = 2^tt$logFC,
    AveExpr = tt$AveExpr,
    t = tt$t,
    pvalue = tt$P.Value,
    padj = tt$adj.P.Val,
    n_pairs = nrow(pairs),
    n_pairs_FC_gt_1 = rowSums(fc_pairs > 1, na.rm = TRUE),
    n_pairs_FC_gt_1.25 = rowSums(fc_pairs > 1.25, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  de <- cbind(de, as.data.frame(fc_pairs))
  de <- de[order(de$pvalue, -de$log2FC), ]
  rownames(de) <- NULL
  attr(de, "pair_expr") <- vec
  attr(de, "pairs") <- pairs
  attr(de, "label") <- label
  de
}

select_up_in_met <- function(de, fc_min, p_min = p_cutoff) {
  if (is.null(de) || !is.data.frame(de) || nrow(de) == 0) {
    return(data.frame())
  }
  de[!is.na(de$pvalue) & de$pvalue < p_min & de$FC > fc_min, , drop = FALSE]
}

# -----------------------------------------------------------------------------
# 6. 绘图 / ORA / GSEA
# -----------------------------------------------------------------------------
save_gg <- function(plot, path_stub, width = 8, height = 6) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".pdf"), plot, width = width, height = height),
           error = function(e) log_msg("pdf ggsave failed: ", e$message))
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".png"), plot, width = width, height = height, dpi = 300),
           error = function(e) log_msg("png ggsave failed: ", e$message))
}

note_empty <- function(stub, msg) writeLines(msg, paste0(stub, "_EMPTY.txt"))

try_save_plot <- function(fun, stub, width = 9, height = 7) {
  p <- tryCatch(fun(), error = function(e) {
    log_msg("Plot failed (", basename(stub), "): ", e$message)
    NULL
  })
  if (is.null(p)) return(invisible(FALSE))
  save_gg(p, stub, width = width, height = height)
}

plot_volcano <- function(de, highlight, title, outfile, fc_line = 1) {
  df <- de
  df$y <- -log10(pmax(df$pvalue, 1e-300))
  df$set <- ifelse(df$gene %in% highlight, "selected", "other")
  df$label <- ifelse(df$gene %in% utils::head(highlight, 15), df$gene, NA)
  lfc_line <- log2(fc_line)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = y, color = set)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.4) +
    ggplot2::scale_color_manual(values = c(other = "grey70", selected = "#D62828")) +
    ggplot2::geom_vline(xintercept = c(-lfc_line, lfc_line), linetype = 2, color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = 2, color = "grey40") +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = title,
      x = "log2FC (matched metastasis / primary); >0 = lower in primary",
      y = "-log10(p value)", color = NULL
    )
  save_gg(p, outfile)
}

plot_de_bar <- function(sub, title, outfile) {
  if (nrow(sub) == 0) return(invisible(NULL))
  df <- sub[order(sub$log2FC, decreasing = TRUE), , drop = FALSE]
  if (nrow(df) > 60) df <- rbind(utils::head(df, 30), utils::tail(df, 30))
  df$gene <- factor(df$gene, levels = rev(unique(df$gene)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = gene, y = log2FC)) +
    ggplot2::geom_col(fill = "#D62828", width = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title, x = NULL, y = "log2FC (met / primary)")
  save_gg(p, outfile, width = 8, height = max(5, min(16, 0.22 * nrow(df) + 2)))
}

group_heatmap <- function(mat, samples, groups, genes, title, outfile) {
  genes <- intersect(genes, rownames(mat))
  samples <- intersect(samples, colnames(mat))
  if (length(genes) > 200) genes <- genes[seq_len(200)]
  if (length(genes) < 2 || length(samples) < 2) {
    log_msg("Group heatmap skipped: ", title)
    return(invisible(NULL))
  }
  sub <- mat[genes, samples, drop = FALSE]
  ann <- data.frame(Group = groups[match(samples, names(groups))], row.names = samples)
  default_pal <- c(
    low = "#4C78A8", high = "#D62828",
    Primary = "#4C78A8", Metastasis = "#D62828",
    lung_only = "#D62828", pleura_only = "#4C78A8", bone_only = "#2A9D8F"
  )
  ug <- unique(as.character(ann$Group))
  pal_g <- default_pal[names(default_pal) %in% ug]
  missing <- setdiff(ug, names(pal_g))
  if (length(missing) > 0) {
    extra_cols <- RColorBrewer::brewer.pal(max(3, length(missing)), "Set2")[seq_along(missing)]
    pal_g <- c(pal_g, stats::setNames(extra_cols, missing))
  }
  pal <- list(Group = pal_g)
  draw <- function() {
    pheatmap::pheatmap(
      sub, scale = "row", annotation_col = ann, annotation_colors = pal,
      cluster_cols = TRUE, cluster_rows = TRUE,
      show_rownames = nrow(sub) <= 80, fontsize_row = 6, main = title,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100)
    )
  }
  grDevices::pdf(paste0(outfile, ".pdf"), width = 8, height = max(6, min(18, 0.18 * nrow(sub) + 3)))
  on.exit({
    while (grDevices::dev.cur() > 1) grDevices::dev.off()
  }, add = TRUE)
  draw()
  grDevices::dev.off()
  grDevices::png(paste0(outfile, ".png"), width = 2400, height = max(1800, 40 * nrow(sub) + 400), res = 300)
  draw()
  grDevices::dev.off()
}

paired_heatmap <- function(de_obj, genes, title, outfile) {
  vec <- attr(de_obj, "pair_expr")
  genes <- intersect(genes, rownames(vec$primary))
  if (length(genes) > 200) genes <- genes[seq_len(200)]
  if (length(genes) < 2) {
    log_msg("Heatmap skipped (<2 genes): ", title)
    return(invisible(NULL))
  }
  mat <- matrix(NA_real_, nrow = length(genes), ncol = 2 * length(vec$patients),
                dimnames = list(genes, NULL))
  cn <- character(0)
  ann <- data.frame(Patient = character(0), Tissue = character(0), stringsAsFactors = FALSE)
  col_i <- 1
  for (j in seq_along(vec$patients)) {
    mat[, col_i] <- vec$primary[genes, j]
    mat[, col_i + 1] <- vec$met[genes, j]
    cn <- c(cn, paste0(vec$patients[j], "_P"), paste0(vec$patients[j], "_M"))
    ann <- rbind(ann, data.frame(
      Patient = vec$patients[j], Tissue = "Primary", stringsAsFactors = FALSE
    ))
    ann <- rbind(ann, data.frame(
      Patient = vec$patients[j], Tissue = "Metastasis", stringsAsFactors = FALSE
    ))
    col_i <- col_i + 2
  }
  colnames(mat) <- cn
  rownames(ann) <- cn
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  pal <- list(
    Tissue = c(Primary = "#4C78A8", Metastasis = "#D62828")
  )
  draw <- function() {
    pheatmap::pheatmap(
      z, scale = "none", annotation_col = ann, annotation_colors = pal,
      cluster_cols = FALSE, cluster_rows = TRUE,
      show_rownames = nrow(z) <= 80, fontsize_row = 6, main = title,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100)
    )
  }
  grDevices::pdf(paste0(outfile, ".pdf"), width = max(8, 0.45 * ncol(z) + 4),
                 height = max(6, min(18, 0.18 * nrow(z) + 3)))
  on.exit({
    while (grDevices::dev.cur() > 1) grDevices::dev.off()
  }, add = TRUE)
  draw()
  grDevices::dev.off()
  grDevices::png(paste0(outfile, ".png"), width = max(1800, 80 * ncol(z)),
                 height = max(1800, 40 * nrow(z) + 400), res = 300)
  draw()
  grDevices::dev.off()
}

pick_official_symbol <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.[0-9]+$", "", x)
  if (length(x) != 1 || is.na(x) || x %in% c("", "-", ".", "NA")) return(NA_character_)
  parts <- unlist(strsplit(x, "[,;|/]+"))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) return(NA_character_)
  parts[1]
}

map_to_entrez <- function(symbols) {
  symbols <- unique(vapply(symbols, pick_official_symbol, character(1), USE.NAMES = FALSE))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  if (length(symbols) == 0) return(data.frame(gene = character(), entrez = character()))
  m <- tryCatch(
    clusterProfiler::bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
  )
  if (nrow(m) == 0) return(data.frame(gene = character(), entrez = character()))
  m <- m[!duplicated(m[[1]]), ]
  data.frame(gene = m[[1]], entrez = m[[2]], stringsAsFactors = FALSE)
}

ranked_entrez <- function(de) {
  mp <- map_to_entrez(de$gene)
  de2 <- merge(de, mp, by = "gene")
  de2 <- de2[!is.na(de2$entrez) & !is.na(de2$log2FC), ]
  de2 <- de2[order(abs(de2$log2FC), decreasing = TRUE), ]
  de2 <- de2[!duplicated(de2$entrez), ]
  stats <- de2$log2FC
  names(stats) <- de2$entrez
  sort(stats, decreasing = TRUE)
}

enrich_or_relax <- function(strict_fun, relax_fun, label) {
  obj <- tryCatch(strict_fun(), error = function(e) {
    log_msg(label, " strict failed: ", e$message)
    NULL
  })
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    attr(obj, "relaxed") <- FALSE
    return(obj)
  }
  obj2 <- tryCatch(relax_fun(), error = function(e) {
    log_msg(label, " relaxed failed: ", e$message)
    NULL
  })
  if (!is.null(obj2)) attr(obj2, "relaxed") <- TRUE
  obj2
}

plot_ora_object <- function(x, stub, title, fold_change = NULL) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) {
    note_empty(stub, "no enrichment terms")
    return(invisible(NULL))
  }
  df <- as.data.frame(x)
  utils::write.csv(df, paste0(stub, ".csv"), row.names = FALSE)
  nshow <- min(15, nrow(df))
  try_save_plot(function() enrichplot::dotplot(x, showCategory = nshow) + ggplot2::ggtitle(title),
                paste0(stub, "_dotplot"), 9, 7)
  try_save_plot(function() enrichplot::barplot(x, showCategory = nshow) + ggplot2::ggtitle(title),
                paste0(stub, "_barplot"), 9, 7)
}

run_ora_plots <- function(genes, de_sub, outdir, label, tag) {
  go_dir <- file.path(outdir, "GO")
  pw_dir <- file.path(outdir, "Pathway")
  kg_dir <- file.path(outdir, "KEGG")
  dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(kg_dir, recursive = TRUE, showWarnings = FALSE)
  pref <- paste0(tag, "_")
  mp <- map_to_entrez(genes)
  entrez <- unique(mp$entrez)
  fc_sym <- setNames(de_sub$log2FC, de_sub$gene)
  if (length(entrez) < 3) {
    writeLines(paste("mapped_entrez", length(entrez)), file.path(outdir, paste0(pref, "ORA_skipped.txt")))
    note_empty(file.path(go_dir, paste0(pref, "ORA_GO")), "too few mapped genes")
    note_empty(file.path(pw_dir, paste0(pref, "ORA_Pathway")), "too few mapped genes")
    note_empty(file.path(kg_dir, paste0(pref, "ORA_KEGG")), "too few mapped genes")
    return(invisible(NULL))
  }
  for (ont in c("BP", "MF", "CC")) {
    ego <- enrich_or_relax(
      function() clusterProfiler::enrichGO(
        gene = entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
        pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
      ),
      function() clusterProfiler::enrichGO(
        gene = entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
        pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
      ),
      paste("enrichGO", ont)
    )
    plot_ora_object(ego, file.path(go_dir, paste0(pref, "ORA_GO_", ont)),
                    paste(label, "| ORA GO", ont), fold_change = fc_sym)
  }
  ek <- enrich_or_relax(
    function() clusterProfiler::enrichKEGG(
      gene = entrez, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.2
    ),
    function() clusterProfiler::enrichKEGG(
      gene = entrez, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1
    ),
    "enrichKEGG"
  )
  if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) {
    ek <- tryCatch(clusterProfiler::setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
                   error = function(e) ek)
  }
  plot_ora_object(ek, file.path(kg_dir, paste0(pref, "ORA_KEGG")),
                  paste(label, "| ORA KEGG"), fold_change = fc_sym)
  if (has_pkg("ReactomePA")) {
    er <- enrich_or_relax(
      function() ReactomePA::enrichPathway(
        gene = entrez, organism = "human", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
      ),
      function() ReactomePA::enrichPathway(
        gene = entrez, organism = "human", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
      ),
      "enrichPathway"
    )
    plot_ora_object(er, file.path(pw_dir, paste0(pref, "ORA_Reactome_pathway")),
                    paste(label, "| ORA Reactome"), fold_change = fc_sym)
  } else {
    note_empty(file.path(pw_dir, paste0(pref, "ORA_Reactome_pathway")), "ReactomePA not installed")
  }
  writeLines(
    c("GO/Pathway/KEGG 是 ORA，不是 GSEA。",
      "GSEA 在同级 GSEA/ 目录，文件名以 GSEA_ 开头。"),
    file.path(outdir, paste0(pref, "00_ORA_is_not_GSEA.txt"))
  )
}

run_gsea_plots <- function(full_de, sub, outdir, tag, label) {
  gsea_dir <- file.path(outdir, "GSEA")
  dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
  pref <- paste0(tag, "_")
  stats <- ranked_entrez(full_de)
  if (length(stats) < 8) {
    note_empty(file.path(gsea_dir, paste0(pref, "GSEA")), "too few ranked genes")
    return(invisible(NULL))
  }
  kegg <- enrich_or_relax(
    function() clusterProfiler::gseKEGG(
      geneList = stats, organism = "hsa", minGSSize = 5, maxGSSize = 500,
      pvalueCutoff = 0.05, verbose = FALSE, eps = 0
    ),
    function() clusterProfiler::gseKEGG(
      geneList = stats, organism = "hsa", minGSSize = 3, maxGSSize = 500,
      pvalueCutoff = 1, verbose = FALSE, eps = 0
    ),
    paste("gseKEGG", tag)
  )
  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
    kegg <- tryCatch(clusterProfiler::setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
                     error = function(e) kegg)
  }
  plot_ora_object(kegg, file.path(gsea_dir, paste0(pref, "GSEA_KEGG")),
                  paste(label, "| GSEA KEGG"))
  if (has_pkg("msigdbr") && has_pkg("fgsea")) {
    hm <- tryCatch({
      m <- tryCatch(
        msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
        error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H")
      )
      gs <- split(as.character(m$entrez_gene), m$gs_name)
      fgsea::fgsea(pathways = gs, stats = stats, minSize = 5, maxSize = 500)
    }, error = function(e) {
      log_msg("fgsea hallmark failed: ", e$message)
      NULL
    })
    if (!is.null(hm) && nrow(hm) > 0) {
      hm <- hm[order(hm$pval), ]
      utils::write.csv(as.data.frame(hm), file.path(gsea_dir, paste0(pref, "GSEA_Hallmark.csv")),
                       row.names = FALSE)
    }
  }
}

emit_subset <- function(comp_name, de_full, sub, tag, title, outdir, fc_line) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(paste("comparison:", comp_name),
      paste("subset:", tag),
      paste("title:", title),
      paste("n_genes:", if (is.null(sub)) 0 else nrow(sub)),
      "FC = matched metastasis / matched primary (same patient).",
      "Selected genes are LOW in the primary relative to the paired metastasis."),
    file.path(outdir, paste0("00_", tag, "_THIS_FOLDER.txt"))
  )
  if (is.null(sub) || nrow(sub) == 0) {
    writeLines("no genes", file.path(outdir, paste0(tag, "_EMPTY.txt")))
    log_msg(comp_name, " ", tag, ": EMPTY")
    return(invisible(NULL))
  }
  utils::write.csv(sub, file.path(outdir, paste0(tag, "_DE_selected_genes.csv")), row.names = FALSE)
  tryCatch(writexl::write_xlsx(sub, file.path(outdir, paste0(tag, "_DE_selected_genes.xlsx"))),
           error = function(e) log_msg("xlsx write failed: ", e$message))
  log_msg(comp_name, " ", tag, ": n = ", nrow(sub))
  tryCatch(plot_de_bar(sub, paste(title, "| genes low in primary"), file.path(outdir, paste0(tag, "_DE_log2FC_barplot"))),
           error = function(e) log_msg("barplot failed: ", e$message))
  tryCatch(plot_volcano(de_full, sub$gene, title, file.path(outdir, paste0(tag, "_volcano")), fc_line = fc_line),
           error = function(e) log_msg("volcano failed: ", e$message))
  hm_mode <- attr(de_full, "heatmap_mode")
  if (identical(hm_mode, "group")) {
    tryCatch(group_heatmap(
      attr(de_full, "heat_mat"), attr(de_full, "heat_samples"),
      attr(de_full, "heat_groups"), sub$gene, title,
      file.path(outdir, paste0(tag, "_heatmap_high_vs_low"))
    ), error = function(e) {
      while (grDevices::dev.cur() > 1) grDevices::dev.off()
      log_msg("group heatmap failed: ", e$message)
    })
  } else {
    tryCatch(paired_heatmap(de_full, sub$gene, title, file.path(outdir, paste0(tag, "_heatmap_paired_1to1"))),
             error = function(e) {
               while (grDevices::dev.cur() > 1) grDevices::dev.off()
               log_msg("heatmap failed: ", e$message)
             })
  }
  tryCatch(run_ora_plots(sub$gene, sub, outdir, title, tag),
           error = function(e) log_msg("ORA failed: ", e$message))
  tryCatch(run_gsea_plots(de_full, sub, outdir, tag, title),
           error = function(e) log_msg("GSEA failed: ", e$message))
}

analyze_paired <- function(comp_name, pairs, folder) {
  base <- file.path(result_dir, folder)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  de <- paired_limma(pairs, comp_name)
  if (is.null(de)) {
    writeLines("not enough pairs", file.path(base, "SKIPPED_not_enough_pairs.txt"))
    return(NULL)
  }
  utils::write.csv(de, file.path(base, paste0(comp_name, "_full_genome_paired_DE.csv")), row.names = FALSE)
  tryCatch(writexl::write_xlsx(de, file.path(base, paste0(comp_name, "_full_genome_paired_DE.xlsx"))),
           error = function(e) log_msg(e$message))
  gsea_all <- file.path(base, "00_GSEA_all_genes_NOT_FC_or_topN")
  dir.create(gsea_all, recursive = TRUE, showWarnings = FALSE)
  tryCatch(run_gsea_plots(de, de, gsea_all, "ALLgenes", paste(comp_name, "all genes")),
           error = function(e) log_msg("all-gene GSEA failed: ", e$message))
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    sub <- select_up_in_met(de, fc_min = fc)
    outdir <- file.path(base, "FoldChange", nm)
    title <- paste0(comp_name, " | p<", p_cutoff, " & FC>", fc, " (low in primary)")
    emit_subset(comp_name, de, sub, nm, title, outdir, fc_line = fc)
  }
  de
}

# -----------------------------------------------------------------------------
# Q1 配对肺 / 骨
# -----------------------------------------------------------------------------
de_lung <- analyze_paired("paired_primary_vs_lung", pairs_lung, "01_lung_paired_met_vs_primary")
de_bone <- analyze_paired("paired_primary_vs_bone", pairs_bone, "02_bone_paired_met_vs_primary")

# -----------------------------------------------------------------------------
# Q2 器官特异：肺对过阈值但骨对不过，或反过来
# -----------------------------------------------------------------------------
organ_dir <- file.path(result_dir, "03_organ_specific")
dir.create(organ_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  c("Each ConvertHER patient has exactly one metastatic site.",
    "Lung-specific = paired lung DE genes not in paired bone DE.",
    "Bone-specific = the reverse. Unpaired tropism compares lung-only vs bone-only primaries."),
  file.path(organ_dir, "NOTE_one_met_site_per_patient.txt")
)

write_setdiff <- function(a, b, name_a, name_b, out_stub, de_ref, folder_tag) {
  genes_a <- if (is.null(a) || nrow(a) == 0) character(0) else a$gene
  genes_b <- if (is.null(b) || nrow(b) == 0) character(0) else b$gene
  only_a <- setdiff(genes_a, genes_b)
  empty <- if (is.null(a)) data.frame() else a[0, ]
  sub <- if (is.null(a) || length(only_a) == 0) empty else a[a$gene %in% only_a, , drop = FALSE]
  title <- paste0(name_a, " specific vs ", name_b)
  emit_subset(title, de_ref, sub, folder_tag, title, out_stub, fc_line = 1)
  invisible(sub)
}

if (!is.null(de_lung) && !is.null(de_bone)) {
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    lung_sub <- select_up_in_met(de_lung, fc)
    bone_sub <- select_up_in_met(de_bone, fc)
    write_setdiff(
      lung_sub, bone_sub, "Lung", "Bone",
      file.path(organ_dir, nm, "lung_specific_not_bone"),
      de_lung, paste0(nm, "_lung_specific")
    )
    write_setdiff(
      bone_sub, lung_sub, "Bone", "Lung",
      file.path(organ_dir, nm, "bone_specific_not_lung"),
      de_bone, paste0(nm, "_bone_specific")
    )
    both <- intersect(lung_sub$gene, bone_sub$gene)
    overlap <- if (length(both) == 0) lung_sub[0, ] else lung_sub[lung_sub$gene %in% both, ]
    emit_subset(
      "lung_and_bone_shared", de_lung, overlap, paste0(nm, "_shared"),
      paste0("Shared lung & bone | p<", p_cutoff, " FC>", fc),
      file.path(organ_dir, nm, "shared_lung_and_bone"), fc
    )
    if (has_pkg("ggvenn")) {
      vdf <- list(Lung = lung_sub$gene, Bone = bone_sub$gene)
      p <- ggvenn::ggvenn(vdf, fill_color = c("#4C78A8", "#D62828")) +
        ggplot2::ggtitle(paste("Organ-specific DE", nm, "p<", p_cutoff))
      save_gg(p, file.path(organ_dir, nm, paste0(nm, "_venn_lung_vs_bone")))
    }
  }
} else {
  writeLines("Need both lung and bone paired DE to call organ-specific genes.",
             file.path(organ_dir, "SKIPPED.txt"))
}

lung_only_pid <- setdiff(pairs_lung$patient, pairs_bone$patient)
bone_only_pid <- setdiff(pairs_bone$patient, pairs_lung$patient)
tropism <- data.frame(
  patient = c(lung_only_pid, bone_only_pid),
  tropism = c(rep("lung_only", length(lung_only_pid)),
              rep("bone_only", length(bone_only_pid))),
  stringsAsFactors = FALSE
)
utils::write.csv(tropism, file.path(organ_dir, "patients_lung_only_vs_bone_only.csv"), row.names = FALSE)
log_msg("Lung-only patients: ", paste(lung_only_pid, collapse = ", "))
log_msg("Bone-only patients: ", paste(bone_only_pid, collapse = ", "))

if (length(lung_only_pid) >= 2 && length(bone_only_pid) >= 2) {
  prim_ids <- c(
    pairs_lung$primary[match(lung_only_pid, pairs_lung$patient)],
    pairs_bone$primary[match(bone_only_pid, pairs_bone$patient)]
  )
  grp <- factor(c(rep("lung_only", length(lung_only_pid)),
                  rep("bone_only", length(bone_only_pid))),
                levels = c("bone_only", "lung_only"))
  design <- stats::model.matrix(~ grp)
  fit <- limma::lmFit(logmat[, prim_ids, drop = FALSE], design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit, coef = "grplung_only", number = Inf)
  trop_de <- data.frame(
    gene = rownames(tt), log2FC = tt$logFC, FC = 2^tt$logFC,
    pvalue = tt$P.Value, padj = tt$adj.P.Val, AveExpr = tt$AveExpr,
    stringsAsFactors = FALSE
  )
  trop_de <- trop_de[order(trop_de$pvalue), ]
  utils::write.csv(trop_de, file.path(organ_dir, "primary_lung_only_vs_bone_only_full_DE.csv"),
                   row.names = FALSE)
  attr(trop_de, "heatmap_mode") <- "group"
  attr(trop_de, "heat_mat") <- logmat
  attr(trop_de, "heat_samples") <- prim_ids
  attr(trop_de, "heat_groups") <- stats::setNames(as.character(grp), prim_ids)
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    sub <- trop_de[!is.na(trop_de$pvalue) & trop_de$pvalue < p_cutoff & trop_de$FC < (1 / fc), ]
    if (nrow(sub) > 0) sub$FC_low_in_lung_primary <- 1 / sub$FC
    emit_subset(
      "primary_low_in_lung_only", trop_de, sub, nm,
      paste0("Primary lower in lung-only vs bone-only | p<", p_cutoff, " FC>", fc),
      file.path(organ_dir, "primary_unpaired_lung_only_low", nm), fc
    )
  }
}

# -----------------------------------------------------------------------------
# Q3 神经浸润签名（原发）
# GEO 无 PNI 病理；用轴突导向 / 施旺 / 神经营养 三套基因在原发打分
# -----------------------------------------------------------------------------
pni_curated <- list(
  # ConvertHER ~105-gene panel: classic SEMA/SLIT/NTRK/SOX10 mostly absent
  axon_guidance = c("DDR1", "EGFR", "MET", "SFRP1", "VEGFA", "ADM", "AXL"),
  schwann = c("ERBB2", "ERBB3", "ERBB4", "CRYAB", "NDRG1", "MIA", "VIM", "EGFR"),
  neurotrophic = c("VEGFA", "GAL", "MET", "EGFR", "FGFR4", "ADM")
)

expand_msig <- function(key_words) {
  if (!has_pkg("msigdbr")) return(character(0))
  m <- tryCatch({
    tryCatch(
      msigdbr::msigdbr(species = "Homo sapiens", collection = "C5"),
      error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "C5")
    )
  }, error = function(e) NULL)
  if (is.null(m) || !("gs_name" %in% names(m))) return(character(0))
  hit <- m[grepl(key_words, m$gs_name, ignore.case = TRUE), ]
  unique(as.character(hit$gene_symbol))
}

pni_sets <- list(
  axon_guidance = unique(intersect(c(pni_curated$axon_guidance, expand_msig("AXON_GUIDANCE")), rownames(logmat))),
  schwann = unique(intersect(c(pni_curated$schwann, expand_msig("SCHWANN")), rownames(logmat))),
  neurotrophic = unique(intersect(c(pni_curated$neurotrophic, expand_msig("NEUROTROPH")), rownames(logmat)))
)

zscore_sig <- function(mat, genes) {
  genes <- intersect(unique(genes), rownames(mat))
  if (length(genes) < 3) {
    return(list(score = setNames(rep(NA_real_, ncol(mat)), colnames(mat)), genes = genes))
  }
  z <- t(scale(t(mat[genes, , drop = FALSE])))
  z[!is.finite(z)] <- NA
  list(score = colMeans(z, na.rm = TRUE), genes = genes)
}

prim_samples <- pheno_use$sample[pheno_use$kind == "primary"]
pni_dir <- file.path(result_dir, "04_PNI_primary")
dir.create(pni_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  c("GSE92977 is a ~105-gene NanoString panel (housekeeping-normalized), not RNA-seq.",
    "Classic SEMA/SLIT/NTRK/SOX10 genes are mostly absent.",
    "PNI scores use panel genes that map to axon/Schwann/neurotrophic pathways.",
    "If a signature has <3 genes in the matrix, that arm is skipped."),
  file.path(pni_dir, "NOTE_panel_limited_PNI_signatures.txt")
)

score_df <- data.frame(
  sample = prim_samples,
  patient = pheno_use$patient[match(prim_samples, pheno_use$sample)],
  has_paired_lung = pheno_use$patient[match(prim_samples, pheno_use$sample)] %in% pairs_lung$patient,
  has_paired_bone = pheno_use$patient[match(prim_samples, pheno_use$sample)] %in% pairs_bone$patient,
  stringsAsFactors = FALSE
)
used_genes <- list()
for (nm in names(pni_sets)) {
  sc <- zscore_sig(logmat[, prim_samples, drop = FALSE], pni_sets[[nm]])
  score_df[[paste0("score_", nm)]] <- sc$score[prim_samples]
  used_genes[[nm]] <- sc$genes
  writeLines(sc$genes, file.path(pni_dir, paste0("signature_genes_", nm, ".txt")))
  log_msg("PNI signature ", nm, ": ", length(sc$genes), " genes in matrix")
}
utils::write.csv(score_df, file.path(pni_dir, "primary_PNI_signature_scores.csv"), row.names = FALSE)

unpaired_limma <- function(expr, group, coef_name) {
  group <- factor(group)
  design <- stats::model.matrix(~ group)
  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  coef <- grep(coef_name, colnames(design), value = TRUE)
  if (length(coef) == 0) coef <- colnames(design)[ncol(design)]
  tt <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
  data.frame(
    gene = rownames(tt), log2FC = tt$logFC, FC = 2^tt$logFC,
    AveExpr = tt$AveExpr, t = tt$t, pvalue = tt$P.Value, padj = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
}

# 高神经浸润 vs 低：基因在高分组更低 => log2FC(high/low) < 0，FC_low/high > 阈值
pni_de_list <- list()
for (nm in names(pni_sets)) {
  base <- file.path(pni_dir, nm)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  if (length(pni_sets[[nm]]) < 3) {
    writeLines(paste("too few genes on NanoString panel:", paste(pni_sets[[nm]], collapse = ", ")),
               file.path(base, "SKIPPED_too_few_signature_genes.txt"))
    log_msg("PNI ", nm, " skipped: only ", length(pni_sets[[nm]]), " panel genes")
    next
  }
  sc <- score_df[[paste0("score_", nm)]]
  ok <- prim_samples[!is.na(sc)]
  sc2 <- sc[match(ok, prim_samples)]
  hi <- ok[sc2 >= stats::median(sc2, na.rm = TRUE)]
  lo <- ok[sc2 < stats::median(sc2, na.rm = TRUE)]
  log_msg("PNI ", nm, " high n=", length(hi), " low n=", length(lo))
  expr <- logmat[, c(lo, hi), drop = FALSE]
  grp <- factor(c(rep("low", length(lo)), rep("high", length(hi))), levels = c("low", "high"))
  de <- unpaired_limma(expr, grp, "high")
  de <- de[order(de$pvalue), ]
  # 高 PNI 原发里更低：FC_high/low < 1
  de$FC_low_over_high <- 1 / de$FC
  pni_de_list[[nm]] <- de
  utils::write.csv(de, file.path(base, paste0(nm, "_high_vs_low_full_DE.csv")), row.names = FALSE)
  grp_ids <- c(lo, hi)
  grp_lab <- setNames(c(rep("low", length(lo)), rep("high", length(hi))), grp_ids)
  for (fcnm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[fcnm]])
    sub <- de[!is.na(de$pvalue) & de$pvalue < p_cutoff & de$FC_low_over_high > fc, ]
    if (nrow(sub) > 0) {
      sub$log2FC <- log2(sub$FC_low_over_high)
      sub$FC <- sub$FC_low_over_high
    }
    full_for_volcano <- de
    full_for_volcano$log2FC <- log2(pmax(full_for_volcano$FC_low_over_high, 1e-8))
    full_for_volcano$FC <- full_for_volcano$FC_low_over_high
    attr(full_for_volcano, "heatmap_mode") <- "group"
    attr(full_for_volcano, "heat_mat") <- logmat
    attr(full_for_volcano, "heat_samples") <- grp_ids
    attr(full_for_volcano, "heat_groups") <- grp_lab
    emit_subset(
      paste0("PNI_", nm), full_for_volcano, sub, fcnm,
      paste0("PNI ", nm, " | genes lower in high-score primaries | p<", p_cutoff, " FC>", fc),
      file.path(base, "FoldChange", fcnm), fc
    )
  }
}

# -----------------------------------------------------------------------------
# Q4 三种神经浸润 vs 肺转移（分别做）
# -----------------------------------------------------------------------------
q4_dir <- file.path(result_dir, "05_PNI_vs_lung_met")
dir.create(q4_dir, recursive = TRUE, showWarnings = FALSE)

score_long <- tidyr::pivot_longer(
  score_df,
  cols = dplyr::starts_with("score_"),
  names_to = "signature",
  values_to = "score"
)
score_long$signature <- sub("^score_", "", score_long$signature)
score_long$lung_status <- ifelse(score_long$has_paired_lung,
                                 "primary_with_paired_lung",
                                 "primary_without_paired_lung")

p <- ggplot2::ggplot(score_long, ggplot2::aes(x = lung_status, y = score, fill = signature)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  ggplot2::geom_jitter(width = 0.12, size = 1.4, alpha = 0.8) +
  ggplot2::facet_wrap(~ signature, scales = "free_y") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)) +
  ggplot2::labs(
    title = "Primary PNI signature scores vs paired lung metastasis",
    x = NULL, y = "z-score (primary tumor)"
  )
save_gg(p, file.path(q4_dir, "PNI_scores_by_lung_pair_status"), width = 11, height = 5)

wilcox_rows <- list()
for (nm in names(pni_sets)) {
  sub <- score_df[!is.na(score_df[[paste0("score_", nm)]]), ]
  a <- sub[[paste0("score_", nm)]][sub$has_paired_lung]
  b <- sub[[paste0("score_", nm)]][!sub$has_paired_lung]
  wt <- tryCatch(stats::wilcox.test(a, b, exact = FALSE), error = function(e) NULL)
  wilcox_rows[[nm]] <- data.frame(
    signature = nm,
    n_primary_with_paired_lung = length(a),
    n_primary_without_paired_lung = length(b),
    median_with_lung = stats::median(a, na.rm = TRUE),
    median_without_lung = stats::median(b, na.rm = TRUE),
    wilcox_p = if (is.null(wt)) NA_real_ else unname(wt$p.value),
    stringsAsFactors = FALSE
  )
}
wilcox_tab <- dplyr::bind_rows(wilcox_rows)
utils::write.csv(wilcox_tab, file.path(q4_dir, "PNI_score_wilcox_paired_lung_vs_not.csv"), row.names = FALSE)

# 配对肺患者：原发 vs 肺转移 的三种签名分数（同一患者一一对应）
if (nrow(pairs_lung) >= 2) {
  vec <- pair_vectors(pairs_lung)
  paired_score <- list()
  for (nm in names(pni_sets)) {
    sc_p <- zscore_sig(vec$primary, pni_sets[[nm]])$score
    sc_m <- zscore_sig(vec$met, pni_sets[[nm]])$score
    df <- data.frame(
      patient = names(sc_p),
      signature = nm,
      score_primary = unname(sc_p),
      score_lung_met = unname(sc_m),
      delta_met_minus_primary = unname(sc_m - sc_p),
      stringsAsFactors = FALSE
    )
    paired_score[[nm]] <- df
    wt <- tryCatch(stats::wilcox.test(df$score_lung_met, df$score_primary, paired = TRUE, exact = FALSE),
                   error = function(e) NULL)
    log_msg("Paired PNI ", nm, " lung vs primary Wilcoxon p=",
            if (is.null(wt)) "NA" else signif(wt$p.value, 3))
    pdf <- tidyr::pivot_longer(df, cols = c("score_primary", "score_lung_met"),
                               names_to = "tissue", values_to = "score")
    pdf$tissue <- ifelse(pdf$tissue == "score_primary", "Primary", "Lung met")
    pp <- ggplot2::ggplot(pdf, ggplot2::aes(x = tissue, y = score, group = patient)) +
      ggplot2::geom_line(color = "grey60") +
      ggplot2::geom_point(ggplot2::aes(color = tissue), size = 2.5) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(
        title = paste0("1-to-1 PNI ", nm, " : primary vs paired lung"),
        x = NULL, y = "signature z-score"
      )
    save_gg(pp, file.path(q4_dir, nm, paste0(nm, "_paired_primary_vs_lung_score")))
  }
  utils::write.csv(dplyr::bind_rows(paired_score),
                   file.path(q4_dir, "paired_PNI_scores_primary_vs_lung.csv"),
                   row.names = FALSE)
}

# 基因交集：原发相对肺转移更低 ∩ 高神经浸润原发里更低
if (!is.null(de_lung)) {
  for (nm in names(pni_de_list)) {
    for (fcnm in names(fc_cutoffs)) {
      fc <- unname(fc_cutoffs[[fcnm]])
      lung_sub <- select_up_in_met(de_lung, fc)
      pni_de <- pni_de_list[[nm]]
      pni_sub <- pni_de[!is.na(pni_de$pvalue) & pni_de$pvalue < p_cutoff & pni_de$FC_low_over_high > fc, ]
      both <- intersect(lung_sub$gene, pni_sub$gene)
      ov <- if (length(both) == 0) lung_sub[0, ] else lung_sub[lung_sub$gene %in% both, ]
      emit_subset(
        paste0("PNI_", nm, "_AND_lung"),
        de_lung, ov, fcnm,
        paste0(nm, " PNI-low genes also low in primary vs paired lung | ", fcnm),
        file.path(q4_dir, nm, "overlap_with_lung_paired", fcnm), fc
      )
    }
  }
}

# -----------------------------------------------------------------------------
# 总表
# -----------------------------------------------------------------------------
summarize_n <- function(de, fc) {
  if (is.null(de) || nrow(de) == 0) return(0L)
  nrow(select_up_in_met(de, fc))
}
summary_tab <- data.frame(
  question = c(
    "Q1_lung_paired", "Q1_bone_paired",
    "Q2_lung_specific_FC1.25", "Q2_bone_specific_FC1.25",
    "Q3_axon_guidance", "Q3_schwann", "Q3_neurotrophic"
  ),
  n_pairs_or_primaries = c(
    nrow(pairs_lung), nrow(pairs_bone),
    nrow(pairs_lung), nrow(pairs_bone),
    sum(!is.na(score_df$score_axon_guidance)),
    sum(!is.na(score_df$score_schwann)),
    sum(!is.na(score_df$score_neurotrophic))
  ),
  n_genes_p0.05_FC_gt_1 = c(
    summarize_n(de_lung, 1), summarize_n(de_bone, 1),
    NA, NA, NA, NA, NA
  ),
  n_genes_p0.05_FC_gt_1.25 = c(
    summarize_n(de_lung, 1.25), summarize_n(de_bone, 1.25),
    NA, NA, NA, NA, NA
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(summary_tab, file.path(result_dir, "00_run_summary.csv"), row.names = FALSE)
writeLines(
  c(
    "GSE92977_Human_breast.R finished.",
    paste("Pairs lung:", nrow(pairs_lung), paste(pairs_lung$patient, collapse = ", ")),
    paste("Pairs bone:", nrow(pairs_bone), paste(pairs_bone$patient, collapse = ", ")),
    "Heatmaps keep patient order: Primary_i next to Met_i (cluster_cols = FALSE).",
    "No Top50-300. Thresholds: p < 0.05 and FC > 1 / 1.25 only.",
    "NanoString ~105-gene panel; PNI signatures are panel-limited.",
    "Each of 123 patients has exactly one matched met site; pairing uses Patient N in the title."
  ),
  file.path(result_dir, "00_README_results.txt")
)
log_msg("Done. Results in ", result_dir)
