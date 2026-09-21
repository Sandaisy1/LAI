#!/usr/bin/env Rscript
# =============================================================================
# MSK 2018 Cancer Cell cBioPortal 乳腺癌：原发灶低剂量 vs 配对肺/骨转移，以及神经浸润
# 文件名：MSK_2018_cancer_cell_cbioportal_breast.R
# （即“MSK_2018 cancer cell_cbioportal breast”）
#
# 数据目录（脚本会自动寻找）：
#   E:/R/cBioportal breast cancer/MSK_2018 cancer cell
#   或其下的 breast_msk_2018/ 解压子目录
#
# 运行：
#   setwd("E:/R/cBioportal breast cancer/MSK_2018 cancer cell")
#   source("MSK_2018_cancer_cell_cbioportal_breast.R")
# 也可设置环境变量 MSK2018_DIR 指向解压后的 breast_msk_2018 目录。
#
# 四个问题：
#   1) 同一患者原发灶 vs 配对肺转移 / 骨转移：原发灶相对转移灶“低表达/低剂量”的基因
#      （1 号原发只配 1 号该器官转移，绝不把不同患者混成一组）
#   2) 在 1) 基础上取器官特异：只配肺、不配骨；或只配骨、不配肺
#   3) 用轴突导向 / 施旺细胞 / 神经营养 三套基因集给原发灶打神经浸润分，
#      再找原发灶里低剂量、与高神经浸润表型相关的基因（三套分开）
#   4) 三种神经浸润评分分别与是否发生肺转移的关系
#
# 阈值：p < 0.05，且 FC >= 1 与 FC >= 1.25 两档。不做 top50–300。
#
# 重要：MSK 2018（Razavi et al. Cancer Cell 2018, breast_msk_2018）公开包是
# MSK-IMPACT 靶向测序，没有全转录组 RNA。本脚本：
#   - 若目录里出现 mRNA 矩阵则优先用 RNA（与 AURORA 脚本同一套 FC 定义）
#   - 否则用离散 GISTIC 拷贝数作为剂量代理（-2/-1/0/+1/+2）
#     FC = 2^(CNA_met - CNA_primary)
#   配对人数通常很小（肺约 3，骨约 10）；n < 3 不伪造 p 值。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")
options(clusterProfiler.download.method = "auto")

# -----------------------------------------------------------------------------
# 0. 依赖
# -----------------------------------------------------------------------------
cran_required <- c(
  "dplyr", "tidyr", "tibble", "stringr", "ggplot2", "ggrepel",
  "pheatmap", "RColorBrewer", "matrixStats"
)
cran_optional <- c("writexl")
bioc_required <- character(0)
bioc_optional <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "AnnotationDbi", "msigdbr")

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
  if (length(still) > 0) message("可选包未安装，相关步骤将跳过: ", paste(still, collapse = ", "))
  invisible(TRUE)
}

install_if_missing(cran_required, bioc = FALSE, required = TRUE)
install_if_missing(cran_optional, bioc = FALSE, required = FALSE)
still_bioc <- bioc_optional[!vapply(bioc_optional, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_bioc) > 0) {
  message("可选 Bioconductor 包未安装，相关步骤将跳过: ", paste(still_bioc, collapse = ", "))
}

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
script_dir <- tryCatch({
  ofile <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(ofile) == 1) {
    dirname(normalizePath(sub("^--file=", "", ofile), winslash = "/", mustWork = FALSE))
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    dirname(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE))
  } else {
    NA_character_
  }
}, error = function(e) NA_character_)

dir_has_msk2018 <- function(d) {
  if (!nzchar(d) || !dir.exists(d)) return(FALSE)
  hits <- list.files(d, pattern = "data_clinical_sample\\.txt$", recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) return(FALSE)
  # 优先确认是 MSK 2018；若只有临床表也接受（用户指定了该目录）
  meta <- list.files(d, pattern = "meta_study\\.txt$", recursive = TRUE, full.names = TRUE)
  if (length(meta) > 0) {
    txt <- paste(readLines(meta[1], warn = FALSE), collapse = "\n")
    if (grepl("breast_msk_2018|MSK, Cancer Cell 2018", txt, ignore.case = TRUE)) return(TRUE)
  }
  TRUE
}

resolve_msk2018_dir <- function() {
  env_dir <- Sys.getenv("MSK2018_DIR", unset = "")
  candidates <- c(
    env_dir,
    "E:/R/cBioportal breast cancer/MSK_2018 cancer cell",
    "E:\\R\\cBioportal breast cancer\\MSK_2018 cancer cell",
    file.path("E:/R/cBioportal breast cancer/MSK_2018 cancer cell", "breast_msk_2018"),
    getwd(),
    file.path(getwd(), "breast_msk_2018"),
    file.path(getwd(), "MSK_2018 cancer cell"),
    script_dir,
    file.path(script_dir, "MSK_2018 cancer cell"),
    file.path(script_dir, "breast_msk_2018")
  )
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  for (d in candidates) {
    if (dir_has_msk2018(d)) {
      clin <- list.files(d, pattern = "data_clinical_sample\\.txt$", recursive = TRUE, full.names = TRUE)
      return(normalizePath(dirname(clin[1]), winslash = "/", mustWork = FALSE))
    }
  }
  stop(
    "未找到 MSK 2018 临床表 data_clinical_sample.txt。请把解压后的文件放在：",
    "E:/R/cBioportal breast cancer/MSK_2018 cancer cell （或设置环境变量 MSK2018_DIR）"
  )
}

project_dir <- resolve_msk2018_dir()
result_dir  <- file.path(project_dir, "results_MSK_2018_cbioportal_breast")
log_dir     <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("MSK2018_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
log_msg("MSK 2018 data directory: ", project_dir)
log_msg("Results: ", result_dir)

p_cutoff   <- 0.05
fc_cutoffs <- c("FC_1" = 1, "FC_1.25" = 1.25)

save_tbl <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, path, row.names = FALSE, na = "")
  if (has_pkg("writexl") && grepl("\\.csv$", path)) {
    tryCatch(
      writexl::write_xlsx(df, sub("\\.csv$", ".xlsx", path)),
      error = function(e) NULL
    )
  }
}

save_gg <- function(plot, path_stub, width = 8, height = 6) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".pdf"), plot, width = width, height = height),
           error = function(e) log_msg("pdf ggsave failed: ", e$message))
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".png"), plot, width = width, height = height, dpi = 300),
           error = function(e) log_msg("png ggsave failed: ", e$message))
}

# -----------------------------------------------------------------------------
# 2. 读临床表 + 表达/拷贝数矩阵
# -----------------------------------------------------------------------------
read_cbioportal_clinical <- function(path) {
  raw <- readLines(path, warn = FALSE, encoding = "UTF-8")
  hit <- which(grepl("^PATIENT_ID\\t", raw))
  if (length(hit) == 0) stop("临床表没有 PATIENT_ID 表头: ", path)
  utils::read.delim(path, skip = hit[1] - 1, check.names = FALSE, stringsAsFactors = FALSE)
}

find_one <- function(patterns) {
  allf <- list.files(project_dir, recursive = TRUE, full.names = TRUE)
  base <- basename(allf)
  for (p in patterns) {
    hit <- allf[grepl(p, base, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NULL
}

matrix_from_portal <- function(path) {
  raw <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  gene_col <- if ("Hugo_Symbol" %in% names(raw)) "Hugo_Symbol" else names(raw)[1]
  entrez_col <- if ("Entrez_Gene_Id" %in% names(raw)) "Entrez_Gene_Id" else NULL
  sample_cols <- setdiff(names(raw), c(gene_col, entrez_col))
  mat <- as.matrix(raw[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "double"
  genes <- trimws(as.character(raw[[gene_col]]))
  keep <- !is.na(genes) & nzchar(genes) & genes != "NA"
  mat <- mat[keep, , drop = FALSE]
  genes <- genes[keep]
  if (any(duplicated(genes))) {
    means <- matrixStats::rowMeans2(mat, na.rm = TRUE)
    ord <- order(means, decreasing = TRUE, na.last = TRUE)
    mat <- mat[ord, , drop = FALSE]
    genes <- genes[ord]
    uniq <- !duplicated(genes)
    mat <- mat[uniq, , drop = FALSE]
    genes <- genes[uniq]
  }
  rownames(mat) <- genes
  mat
}

clin_sample_file <- find_one("^data_clinical_sample\\.txt$")
clin_patient_file <- find_one("^data_clinical_patient\\.txt$")
if (is.null(clin_sample_file)) stop("缺少 data_clinical_sample.txt")
clin <- read_cbioportal_clinical(clin_sample_file)
patient <- if (!is.null(clin_patient_file)) read_cbioportal_clinical(clin_patient_file) else NULL
log_msg("Clinical samples: ", nrow(clin), " from ", basename(clin_sample_file))

# 优先 RNA；MSK 2018 公开包没有 RNA，回退到离散 CNA
expr_file <- find_one(c(
  "^data_mrna_seq_v2_rsem\\.txt$",
  "^data_mrna_seq_rsem\\.txt$",
  "^data_mrna.*rsem\\.txt$",
  "^data_mrna.*counts\\.txt$",
  "^data_mrna.*fpkm\\.txt$"
))
expr_is_zscore <- FALSE
assay_mode <- "RNA"
if (is.null(expr_file)) {
  zf <- find_one(c(
    "^data_mrna_seq_v2_rsem_zscores.*\\.txt$",
    "^data_mrna.*zscores.*\\.txt$"
  ))
  if (!is.null(zf)) {
    expr_file <- zf
    expr_is_zscore <- TRUE
    assay_mode <- "RNA_zscore"
  }
}
if (is.null(expr_file)) {
  expr_file <- find_one(c("^data_cna\\.txt$", "^data_linear_cna\\.txt$"))
  assay_mode <- if (!is.null(expr_file) && grepl("linear", basename(expr_file), ignore.case = TRUE)) {
    "CNA_linear"
  } else {
    "CNA_discrete"
  }
}
if (is.null(expr_file)) stop("未找到 mRNA 或拷贝数矩阵（data_cna.txt）")
log_msg("Assay file: ", basename(expr_file), " | mode=", assay_mode, " | z-score=", expr_is_zscore)
if (grepl("^CNA", assay_mode)) {
  log_msg("NOTE: MSK 2018 公开包无 RNA。用 IMPACT 离散拷贝数当剂量代理。",
          "低表达 ≈ 原发灶拷贝数低于配对转移灶；FC = 2^(CNA_met - CNA_primary)。",
          "这不是转录组倍数，解读时必须写明。")
}

mat <- matrix_from_portal(expr_file)
log_msg("Matrix: ", nrow(mat), " genes x ", ncol(mat), " samples")

# -----------------------------------------------------------------------------
# 3. 器官归并 + 同一患者 1 对 1 配对
# -----------------------------------------------------------------------------
# MSK 2018 用 SAMPLE_SITE（原发为 Treatment Naive Primary 等；转移为 Lung/Bone/...）
site_col <- if ("SAMPLE_SITE" %in% names(clin)) {
  "SAMPLE_SITE"
} else if ("METASTATIC_SITE" %in% names(clin)) {
  "METASTATIC_SITE"
} else {
  stop("临床表没有 SAMPLE_SITE / METASTATIC_SITE")
}

classify_organ <- function(site, sample_type) {
  s <- tolower(trimws(as.character(site)))
  s[is.na(s)] <- ""
  st <- tolower(trimws(as.character(sample_type)))
  out <- rep(NA_character_, length(s))
  # 原发灶部位不是转移器官
  is_met <- grepl("metasta", st)
  lung <- is_met & grepl("lung|lobe", s) & !grepl("liver", s)
  bone <- is_met & grepl("bone|vertebra|rib|iliac|spine|femur|sternum|calvarium|perioste", s)
  brain <- is_met & grepl("brain|cereb|dura|occipital", s)
  liver <- is_met & grepl("liver", s)
  out[liver] <- "Liver"
  out[brain] <- "Brain"
  out[bone] <- "Bone"
  out[lung] <- "Lung"
  out[is_met & s != "" & is.na(out)] <- "Other"
  out
}

clin$SAMPLE_ID <- as.character(clin$SAMPLE_ID)
clin$PATIENT_ID <- as.character(clin$PATIENT_ID)
clin$SAMPLE_TYPE <- as.character(clin$SAMPLE_TYPE)
clin$SAMPLE_TYPE[clin$SAMPLE_TYPE %in% c("Metastasis", "metastasis")] <- "Metastatic"
clin$SITE_RAW <- as.character(clin[[site_col]])
clin$organ <- classify_organ(clin$SITE_RAW, clin$SAMPLE_TYPE)
clin$has_assay <- clin$SAMPLE_ID %in% colnames(mat)

log_msg("Assay by sample type: ",
        paste(names(table(clin$SAMPLE_TYPE[clin$has_assay])),
              table(clin$SAMPLE_TYPE[clin$has_assay]), sep = "=", collapse = ", "))
log_msg("Metastatic organs with assay: ",
        paste(names(table(clin$organ[clin$has_assay & clin$SAMPLE_TYPE == "Metastatic"], useNA = "ifany")),
              table(clin$organ[clin$has_assay & clin$SAMPLE_TYPE == "Metastatic"], useNA = "ifany"),
              sep = "=", collapse = ", "))

site_tab <- as.data.frame(table(SAMPLE_TYPE = clin$SAMPLE_TYPE, SITE = clin$SITE_RAW), stringsAsFactors = FALSE)
site_tab <- site_tab[site_tab$Freq > 0, ]
site_tab$organ <- classify_organ(site_tab$SITE, site_tab$SAMPLE_TYPE)
save_tbl(site_tab[order(-site_tab$Freq), ], file.path(result_dir, "01_pairing", "sample_site_to_organ.csv"))

pick_one_primary <- function(rows) {
  if (nrow(rows) == 1) return(rows$SAMPLE_ID[1])
  naive <- grepl("treatment naive|treatment-naive|naive primary", rows$SITE_RAW, ignore.case = TRUE)
  if (any(naive)) rows <- rows[naive, , drop = FALSE]
  if (nrow(rows) == 1) return(rows$SAMPLE_ID[1])
  nuc <- if ("PERCENT_TUMOR_NUCLEI" %in% names(rows)) {
    suppressWarnings(as.numeric(rows$PERCENT_TUMOR_NUCLEI))
  } else {
    rep(NA_real_, nrow(rows))
  }
  if (all(is.na(nuc))) return(rows$SAMPLE_ID[1])
  rows$SAMPLE_ID[which.max(replace(nuc, is.na(nuc), -Inf))]
}

# 同一患者、同一器官：多枚转移灶先平均，保证 1 个原发 : 1 个该器官转移轮廓
build_pairs <- function(organ_name) {
  prim <- clin[which(clin$SAMPLE_TYPE == "Primary" & clin$has_assay), , drop = FALSE]
  mets <- clin[which(clin$SAMPLE_TYPE == "Metastatic" & clin$has_assay &
                       !is.na(clin$organ) & clin$organ == organ_name), , drop = FALSE]
  both <- intersect(unique(prim$PATIENT_ID), unique(mets$PATIENT_ID))
  if (length(both) == 0) {
    return(data.frame(
      patient_id = character(), primary_sample = character(), met_sample = character(),
      primary_site = character(), met_site = character(), n_met = integer(),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(both, function(pid) {
    pr <- prim[prim$PATIENT_ID == pid, , drop = FALSE]
    mt <- mets[mets$PATIENT_ID == pid, , drop = FALSE]
    data.frame(
      patient_id = pid,
      primary_sample = pick_one_primary(pr),
      met_sample = paste(mt$SAMPLE_ID, collapse = ";"),
      primary_site = paste(pr$SITE_RAW[match(pick_one_primary(pr), pr$SAMPLE_ID)], collapse = ";"),
      met_site = paste(mt$SITE_RAW, collapse = ";"),
      n_met = nrow(mt),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

pair_expr <- function(pairs) {
  if (nrow(pairs) == 0) {
    return(list(
      primary = matrix(numeric(0), nrow = nrow(mat), ncol = 0, dimnames = list(rownames(mat), NULL)),
      met = matrix(numeric(0), nrow = nrow(mat), ncol = 0, dimnames = list(rownames(mat), NULL))
    ))
  }
  prim_mat <- mat[, pairs$primary_sample, drop = FALSE]
  met_mat <- sapply(seq_len(nrow(pairs)), function(i) {
    ids <- strsplit(pairs$met_sample[i], ";", fixed = TRUE)[[1]]
    ids <- intersect(ids, colnames(mat))
    if (length(ids) == 1) return(mat[, ids])
    matrixStats::rowMeans2(mat[, ids, drop = FALSE], na.rm = TRUE)
  })
  if (is.null(dim(met_mat))) met_mat <- matrix(met_mat, ncol = 1)
  colnames(prim_mat) <- pairs$patient_id
  colnames(met_mat) <- pairs$patient_id
  rownames(met_mat) <- rownames(mat)
  list(primary = prim_mat, met = met_mat)
}

pairs_lung <- build_pairs("Lung")
pairs_bone <- build_pairs("Bone")
save_tbl(pairs_lung, file.path(result_dir, "01_pairing", "pairs_primary_vs_lung.csv"))
save_tbl(pairs_bone, file.path(result_dir, "01_pairing", "pairs_primary_vs_bone.csv"))
log_msg("Paired primary–lung patients: ", nrow(pairs_lung),
        " | ", paste(pairs_lung$patient_id, collapse = ", "))
log_msg("Paired primary–bone patients: ", nrow(pairs_bone),
        " | ", paste(pairs_bone$patient_id, collapse = ", "))
if (nrow(pairs_lung) < 3) {
  log_msg("WARNING: 肺转移配对 < 3，不能稳定估计 p 值，不伪造 p。")
} else if (nrow(pairs_lung) == 3) {
  log_msg("WARNING: 肺转移配对仅 3 人，p 值可算但不稳定，结果只作探索。")
}
if (nrow(pairs_bone) < 3) {
  log_msg("WARNING: 骨转移配对 < 3，不能稳定估计 p 值，不伪造 p。")
}

writeLines(
  c(
    "配对规则：同一 PATIENT_ID；1 个原发灶对应该患者该器官的转移灶。",
    "同一器官多枚转移灶：先对基因剂量取平均，再与原发比较。",
    "原发若有多枚：优先 Treatment Naive Primary。",
    paste0("当前矩阵模式: ", assay_mode),
    paste0("肺配对 n=", nrow(pairs_lung), "  骨配对 n=", nrow(pairs_bone)),
    "不要把不同患者的原发/转移混成两组做非配对检验来回答问题 1。"
  ),
  file.path(result_dir, "01_pairing", "PAIRING_README.txt")
)

# -----------------------------------------------------------------------------
# 4. 配对差异：原发低剂量 = 转移灶相对原发灶升高
# -----------------------------------------------------------------------------
paired_de <- function(pairs, label) {
  pe <- pair_expr(pairs)
  n <- ncol(pe$primary)
  genes_all <- rownames(pe$primary)
  if (n == 0) {
    log_msg(label, ": 0 pairs, skip")
    return(data.frame(gene = genes_all, n_pairs = 0, mean_primary = NA_real_,
                      mean_met = NA_real_, log2FC = NA_real_, FC = NA_real_,
                      pvalue = NA_real_, padj = NA_real_, direction = NA_character_,
                      stringsAsFactors = FALSE))
  }
  delta <- pe$met - pe$primary
  mean_p <- matrixStats::rowMeans2(pe$primary, na.rm = TRUE)
  mean_m <- matrixStats::rowMeans2(pe$met, na.rm = TRUE)
  log2fc <- matrixStats::rowMeans2(delta, na.rm = TRUE)
  n_ok <- rowSums(is.finite(delta))
  pval <- rep(NA_real_, length(genes_all))
  if (n >= 3) {
    pval <- apply(delta, 1, function(x) {
      x <- x[is.finite(x)]
      if (length(x) < 3) return(NA_real_)
      if (stats::sd(x) == 0) return(if (abs(mean(x)) < 1e-12) 1 else NA_real_)
      tryCatch(stats::t.test(x)$p.value, error = function(e) NA_real_)
    })
  } else {
    log_msg(label, ": n_pairs=", n, " < 3，不估计 p 值（不伪造）")
  }
  padj <- if (all(is.na(pval))) rep(NA_real_, length(pval)) else p.adjust(pval, method = "BH")
  data.frame(
    gene = genes_all,
    n_pairs = n_ok,
    mean_primary = mean_p,
    mean_met = mean_m,
    log2FC = log2fc,
    FC = 2^log2fc,
    pvalue = pval,
    padj = padj,
    direction = ifelse(log2fc > 0, "low_in_primary_vs_matched_met",
                       ifelse(log2fc < 0, "high_in_primary_vs_matched_met", "unchanged")),
    stringsAsFactors = FALSE
  )
}

select_low_in_primary <- function(de, fc_min, p_min = p_cutoff) {
  ok_fc <- is.finite(de$FC) & de$FC >= fc_min & de$log2FC > 0
  has_p <- any(is.finite(de$pvalue))
  if (has_p) {
    de[ok_fc & is.finite(de$pvalue) & de$pvalue < p_min, , drop = FALSE]
  } else {
    log_msg("无 p 值：不把 FC 列表当成显著基因，避免 1 对样本刷出上千“阳性”")
    de[0, , drop = FALSE]
  }
}

xlab_fc <- if (grepl("^CNA", assay_mode)) {
  "paired log2FC (matched met CNA - primary CNA)"
} else {
  "paired log2FC (matched met / primary)"
}

plot_volcano <- function(de, highlight, title, outfile, fc_line = 1) {
  df <- de
  has_p <- "pvalue" %in% names(df) && any(is.finite(df$pvalue))
  if (has_p) {
    df$y <- -log10(pmax(df$pvalue, 1e-300))
    ylab <- "-log10(p value)"
    hline <- -log10(p_cutoff)
  } else {
    df$y <- abs(df$log2FC)
    ylab <- "|paired log2FC| (p not estimated)"
    hline <- NULL
  }
  df$set <- ifelse(df$gene %in% highlight, "selected", "other")
  df$label <- ifelse(df$gene %in% utils::head(highlight, 20), df$gene, NA)
  lfc_line <- log2(fc_line)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = y, color = set)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.4) +
    ggplot2::scale_color_manual(values = c(other = "grey70", selected = "#D62828")) +
    ggplot2::geom_vline(xintercept = c(-lfc_line, lfc_line), linetype = 2, color = "grey40") +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, x = xlab_fc, y = ylab, color = NULL)
  if (!is.null(hline)) p <- p + ggplot2::geom_hline(yintercept = hline, linetype = 2, color = "grey40")
  save_gg(p, outfile)
}

plot_pair_heatmap <- function(pairs, genes, title, outfile) {
  genes <- intersect(genes, rownames(mat))
  if (length(genes) > 200) genes <- genes[seq_len(200)]
  if (length(genes) < 2 || nrow(pairs) == 0) {
    log_msg("Heatmap skipped: ", title)
    return(invisible(NULL))
  }
  pe <- pair_expr(pairs)
  heat <- cbind(pe$primary[genes, , drop = FALSE], pe$met[genes, , drop = FALSE])
  colnames(heat) <- c(paste0(pairs$patient_id, "_Primary"), paste0(pairs$patient_id, "_Met"))
  interleaved <- as.vector(rbind(
    paste0(pairs$patient_id, "_Primary"),
    paste0(pairs$patient_id, "_Met")
  ))
  heat <- heat[, interleaved, drop = FALSE]
  ann <- data.frame(
    Tissue = ifelse(grepl("_Primary$", colnames(heat)), "Primary", "Matched_met"),
    Patient = gsub("_(Primary|Met)$", "", colnames(heat)),
    row.names = colnames(heat)
  )
  pal <- list(Tissue = c(Primary = "#4C78A8", Matched_met = "#D62828"))
  use_scale <- if (grepl("^CNA_discrete", assay_mode)) "none" else "row"
  grDevices::pdf(paste0(outfile, ".pdf"), width = max(8, 0.45 * ncol(heat) + 4),
                 height = max(6, min(18, 0.18 * nrow(heat) + 3)))
  pheatmap::pheatmap(
    heat, scale = use_scale, annotation_col = ann, annotation_colors = pal,
    show_rownames = nrow(heat) <= 80, fontsize_row = 6, main = title,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100),
    cluster_cols = FALSE
  )
  grDevices::dev.off()
}

map_to_entrez <- function(symbols) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  if (length(symbols) == 0 || !has_pkg("org.Hs.eg.db") || !has_pkg("AnnotationDbi")) {
    return(data.frame(gene = character(), entrez = character()))
  }
  mp <- tryCatch(
    AnnotationDbi::select(org.Hs.eg.db, keys = symbols, keytype = "SYMBOL",
                          columns = c("ENTREZID", "SYMBOL")),
    error = function(e) NULL
  )
  if (is.null(mp) || nrow(mp) == 0) return(data.frame(gene = character(), entrez = character()))
  mp <- mp[!is.na(mp$ENTREZID), , drop = FALSE]
  data.frame(gene = mp$SYMBOL, entrez = as.character(mp$ENTREZID), stringsAsFactors = FALSE)
}

run_ora <- function(genes, de_sub, outdir, label, tag) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(as.character(genes), file.path(outdir, paste0(tag, "_gene_list.txt")))
  if (!has_pkg("clusterProfiler") || !has_pkg("org.Hs.eg.db")) {
    writeLines("clusterProfiler / org.Hs.eg.db 未安装，已跳过 ORA。",
               file.path(outdir, paste0(tag, "_ORA_skipped.txt")))
    return(invisible(NULL))
  }
  mp <- map_to_entrez(genes)
  entrez <- unique(mp$entrez)
  if (length(entrez) < 5) {
    writeLines(paste("mapped_entrez", length(entrez)), file.path(outdir, paste0(tag, "_ORA_skipped.txt")))
    return(invisible(NULL))
  }
  ego <- tryCatch(
    clusterProfiler::enrichGO(
      gene = entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
    ),
    error = function(e) NULL
  )
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    save_tbl(as.data.frame(ego), file.path(outdir, paste0(tag, "_ORA_GO_BP.csv")))
    if (has_pkg("enrichplot")) {
      tryCatch({
        p <- enrichplot::dotplot(ego, showCategory = 15) +
          ggplot2::ggtitle(paste(label, "| ORA GO BP"))
        save_gg(p, file.path(outdir, paste0(tag, "_ORA_GO_BP_dotplot")), width = 9, height = 7)
      }, error = function(e) NULL)
    }
  }
  ek <- tryCatch(
    clusterProfiler::enrichKEGG(gene = entrez, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
    error = function(e) NULL
  )
  if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) {
    if (has_pkg("AnnotationDbi")) {
      ek <- tryCatch(clusterProfiler::setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"), error = function(e) ek)
    }
    save_tbl(as.data.frame(ek), file.path(outdir, paste0(tag, "_ORA_KEGG.csv")))
    if (has_pkg("enrichplot")) {
      tryCatch({
        p <- enrichplot::dotplot(ek, showCategory = 15) +
          ggplot2::ggtitle(paste(label, "| ORA KEGG"))
        save_gg(p, file.path(outdir, paste0(tag, "_ORA_KEGG_dotplot")), width = 9, height = 7)
      }, error = function(e) NULL)
    }
  }
}

export_de_set <- function(de, pairs, organ, fc_min, tag, out_root) {
  hit <- select_low_in_primary(de, fc_min)
  hit <- hit[order(hit$pvalue, -hit$FC, na.last = TRUE), , drop = FALSE]
  outdir <- file.path(out_root, tag)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ranked <- de[order(de$pvalue, -de$FC, na.last = TRUE), ]
  save_tbl(ranked, file.path(outdir, paste0(tag, "_all_genes.csv")))
  save_tbl(hit, file.path(outdir, paste0(tag, "_low_in_primary.csv")))
  if (!any(is.finite(de$pvalue))) {
    fc_only <- ranked[is.finite(ranked$FC) & ranked$FC >= fc_min & ranked$log2FC > 0, , drop = FALSE]
    save_tbl(fc_only, file.path(outdir, paste0(tag, "_FC_ranked_NO_PVALUE.csv")))
    writeLines(
      "配对 < 3，p 值未估计。本表仅按配对 FC 排序，不能当作显著差异基因。",
      file.path(outdir, paste0(tag, "_NO_PVALUE_README.txt"))
    )
  }
  log_msg(organ, " ", tag, ": ", nrow(hit), " genes (low in paired primary, FC>=", fc_min,
          if (any(is.finite(de$pvalue))) paste0(", p<", p_cutoff) else ", no p", ")")
  plot_volcano(
    de, hit$gene,
    paste0("MSK2018 | paired primary vs ", organ, " met | ", tag, " | ", assay_mode),
    file.path(outdir, paste0(tag, "_volcano")),
    fc_line = fc_min
  )
  plot_pair_heatmap(
    pairs, hit$gene,
    paste0("MSK2018 paired ", organ, " | ", tag, " low in primary"),
    file.path(outdir, paste0(tag, "_heatmap_paired"))
  )
  run_ora(hit$gene, de, file.path(outdir, "ORA"), paste("MSK2018", organ, tag), tag)
  hit
}

log_msg("==== Q1: paired primary vs lung / bone, genes low in primary ====")
de_lung <- paired_de(pairs_lung, "Lung")
de_bone <- paired_de(pairs_bone, "Bone")
save_tbl(de_lung, file.path(result_dir, "02_lung", "lung_paired_all_genes.csv"))
save_tbl(de_bone, file.path(result_dir, "03_bone", "bone_paired_all_genes.csv"))

lung_sets <- list()
bone_sets <- list()
for (nm in names(fc_cutoffs)) {
  lung_sets[[nm]] <- export_de_set(de_lung, pairs_lung, "Lung", fc_cutoffs[[nm]], nm,
                                   file.path(result_dir, "02_lung"))
  bone_sets[[nm]] <- export_de_set(de_bone, pairs_bone, "Bone", fc_cutoffs[[nm]], nm,
                                   file.path(result_dir, "03_bone"))
}

# -----------------------------------------------------------------------------
# 5. Q2 器官特异：只肺或只骨
# -----------------------------------------------------------------------------
log_msg("==== Q2: organ-specific low dosage in paired primary ====")
org_dir <- file.path(result_dir, "04_organ_specific")
dir.create(org_dir, recursive = TRUE, showWarnings = FALSE)

patient_organs <- clin %>%
  dplyr::filter(.data$SAMPLE_TYPE == "Metastatic") %>%
  dplyr::group_by(.data$PATIENT_ID) %>%
  dplyr::summarise(
    ever_lung = any(.data$organ == "Lung", na.rm = TRUE),
    ever_bone = any(.data$organ == "Bone", na.rm = TRUE),
    .groups = "drop"
  )

for (nm in names(fc_cutoffs)) {
  g_lung <- lung_sets[[nm]]$gene
  g_bone <- bone_sets[[nm]]$gene
  only_lung <- setdiff(g_lung, g_bone)
  only_bone <- setdiff(g_bone, g_lung)
  both <- intersect(g_lung, g_bone)
  tab <- data.frame(
    gene = c(only_lung, only_bone, both),
    organ_specificity = c(
      rep("lung_only", length(only_lung)),
      rep("bone_only", length(only_bone)),
      rep("lung_and_bone", length(both))
    ),
    stringsAsFactors = FALSE
  )
  if (nrow(tab) > 0) {
    tab$FC_lung <- de_lung$FC[match(tab$gene, de_lung$gene)]
    tab$p_lung <- de_lung$pvalue[match(tab$gene, de_lung$gene)]
    tab$FC_bone <- de_bone$FC[match(tab$gene, de_bone$gene)]
    tab$p_bone <- de_bone$pvalue[match(tab$gene, de_bone$gene)]
  }
  save_tbl(tab, file.path(org_dir, paste0(nm, "_organ_specific_genes.csv")))
  save_tbl(data.frame(gene = only_lung), file.path(org_dir, paste0(nm, "_lung_ONLY.csv")))
  save_tbl(data.frame(gene = only_bone), file.path(org_dir, paste0(nm, "_bone_ONLY.csv")))
  log_msg(nm, " lung-only=", length(only_lung), " bone-only=", length(only_bone),
          " both=", length(both))
  run_ora(only_lung, de_lung, file.path(org_dir, paste0(nm, "_lung_ONLY_ORA")),
          paste("MSK2018 lung-only", nm), paste0(nm, "_lung_ONLY"))
  run_ora(only_bone, de_bone, file.path(org_dir, paste0(nm, "_bone_ONLY_ORA")),
          paste("MSK2018 bone-only", nm), paste0(nm, "_bone_ONLY"))
}

# -----------------------------------------------------------------------------
# 6. 神经浸润基因集（轴突导向 / 施旺 / 神经营养）
# -----------------------------------------------------------------------------
curated_sets <- list(
  axon_guidance = c(
    "ROBO1", "ROBO2", "ROBO3", "SLIT1", "SLIT2", "SLIT3",
    "SEMA3A", "SEMA3B", "SEMA3C", "SEMA3D", "SEMA3E", "SEMA3F", "SEMA3G",
    "SEMA4A", "SEMA4B", "SEMA4C", "SEMA4D", "SEMA4F", "SEMA4G",
    "SEMA5A", "SEMA5B", "SEMA6A", "SEMA6B", "SEMA6C", "SEMA6D", "SEMA7A",
    "NRP1", "NRP2", "PLXNA1", "PLXNA2", "PLXNA3", "PLXNA4",
    "PLXNB1", "PLXNB2", "PLXNB3", "PLXNC1", "PLXND1",
    "EPHA1", "EPHA2", "EPHA3", "EPHA4", "EPHA5", "EPHA7", "EPHA8",
    "EPHB1", "EPHB2", "EPHB3", "EPHB4", "EPHB6",
    "EFNA1", "EFNA2", "EFNA3", "EFNA4", "EFNA5", "EFNB1", "EFNB2", "EFNB3",
    "NTN1", "NTN4", "DCC", "UNC5A", "UNC5B", "UNC5C", "UNC5D", "NEO1",
    "CXCL12", "CXCR4", "MET", "PTK2", "FYN", "RAC1", "CDC42", "RHOA",
    "PAK1", "MAPK1", "MAPK3", "GSK3B", "CDK5", "ABL1", "NCK1", "NCK2",
    "SRGAP1", "SRGAP2", "SRGAP3", "LIMK1", "LIMK2", "CFL1", "ROCK1", "ROCK2",
    "ITGB1", "L1CAM", "NCAM1", "CNTN2"
  ),
  schwann = c(
    "SOX10", "S100B", "MPZ", "MBP", "PMP22", "MAG", "PRX", "EGR2",
    "POU3F1", "POU3F2", "ERBB2", "ERBB3", "NRG1", "NGFR", "GFAP", "GAP43",
    "NCAM1", "L1CAM", "DHH", "CDH19", "MAL", "PLP1", "GJB1", "MPZL1",
    "SCN7A", "NES", "FOXD3", "PAX3", "TFAP2A", "SOX2", "CDH2", "ITGA6",
    "LAMA2", "LAMA4", "LAMB2", "DRP2"
  ),
  neurotrophic = c(
    "NGF", "BDNF", "NTF3", "NTF4", "GDNF", "NRTN", "ARTN", "PSPN", "CNTF",
    "NTRK1", "NTRK2", "NTRK3", "NGFR", "RET", "GFRA1", "GFRA2", "GFRA3", "GFRA4",
    "SORT1", "IRS1", "PIK3CA", "AKT1", "PLCG1", "MTOR", "MAPK1", "MAPK3",
    "CREB1", "BCL2", "BAX", "TP53", "FRS2", "SHC1", "GRB2", "SOS1",
    "VGF", "GAL", "NPY"
  )
)

expand_with_msig <- function(sets) {
  if (!has_pkg("msigdbr")) return(sets)
  ms <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens"), error = function(e) NULL)
  if (is.null(ms) || nrow(ms) == 0) return(sets)
  name_col <- if ("gs_name" %in% names(ms)) "gs_name" else names(ms)[grepl("name", names(ms), ignore.case = TRUE)][1]
  gene_col <- if ("gene_symbol" %in% names(ms)) "gene_symbol" else if ("symbol" %in% names(ms)) "symbol" else "gene_symbol"
  if (is.na(name_col) || !gene_col %in% names(ms)) return(sets)
  pick <- function(keys) {
    hit <- ms[grepl(keys, ms[[name_col]], ignore.case = TRUE), , drop = FALSE]
    unique(as.character(hit[[gene_col]]))
  }
  sets$axon_guidance <- unique(c(sets$axon_guidance, pick("AXON_GUIDANCE|AXON_GUIDE")))
  sets$schwann <- unique(c(sets$schwann, pick("SCHWANN")))
  sets$neurotrophic <- unique(c(
    sets$neurotrophic,
    pick("NEUROTROPHIN|NEUROTROPHIC|NGF_SIGNALLING|BDNF")
  ))
  sets
}

gene_sets_full <- expand_with_msig(curated_sets)
gene_sets <- lapply(gene_sets_full, function(g) intersect(unique(g), rownames(mat)))
log_msg("Gene set sizes on panel: axon=", length(gene_sets$axon_guidance),
        "/", length(unique(gene_sets_full$axon_guidance)),
        " schwann=", length(gene_sets$schwann),
        "/", length(unique(gene_sets_full$schwann)),
        " neurotrophic=", length(gene_sets$neurotrophic),
        "/", length(unique(gene_sets_full$neurotrophic)))
if (grepl("^CNA", assay_mode)) {
  log_msg("IMPACT 面板只有约 474 基因，神经浸润基因集大部分不在芯片上；评分只用面板交集。")
}

score_signature <- function(expr_mat, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) < 3) return(rep(NA_real_, ncol(expr_mat)))
  sub <- expr_mat[genes, , drop = FALSE]
  as.numeric(matrixStats::colMeans2(sub, na.rm = TRUE))
}

prim_ids <- clin$SAMPLE_ID[clin$SAMPLE_TYPE == "Primary" & clin$has_assay]
prim_mat <- mat[, prim_ids, drop = FALSE]
prim_ann <- clin[match(prim_ids, clin$SAMPLE_ID), , drop = FALSE]
prim_ann$score_axon <- score_signature(prim_mat, gene_sets$axon_guidance)
prim_ann$score_schwann <- score_signature(prim_mat, gene_sets$schwann)
prim_ann$score_neurotrophic <- score_signature(prim_mat, gene_sets$neurotrophic)
prim_ann$has_paired_lung <- prim_ann$PATIENT_ID %in% pairs_lung$patient_id
prim_ann$has_paired_bone <- prim_ann$PATIENT_ID %in% pairs_bone$patient_id
if (nrow(patient_organs) > 0) {
  prim_ann$ever_lung <- patient_organs$ever_lung[match(prim_ann$PATIENT_ID, patient_organs$PATIENT_ID)]
  prim_ann$ever_bone <- patient_organs$ever_bone[match(prim_ann$PATIENT_ID, patient_organs$PATIENT_ID)]
} else {
  prim_ann$ever_lung <- NA
  prim_ann$ever_bone <- NA
}
prim_ann$ever_lung[is.na(prim_ann$ever_lung)] <- FALSE
prim_ann$ever_bone[is.na(prim_ann$ever_bone)] <- FALSE
prim_ann$mark_axon <- ifelse(prim_ann$score_axon >= stats::median(prim_ann$score_axon, na.rm = TRUE),
                             "High_axon_guidance", "Low_axon_guidance")
prim_ann$mark_schwann <- ifelse(prim_ann$score_schwann >= stats::median(prim_ann$score_schwann, na.rm = TRUE),
                                "High_Schwann", "Low_Schwann")
prim_ann$mark_neurotrophic <- ifelse(prim_ann$score_neurotrophic >= stats::median(prim_ann$score_neurotrophic, na.rm = TRUE),
                                     "High_neurotrophic", "Low_neurotrophic")
save_tbl(prim_ann, file.path(result_dir, "05_neural_invasion", "primary_neural_invasion_marks.csv"))
log_msg("Marked ", nrow(prim_ann), " assay primaries with 3 neural-invasion scores (median split)")

# 高神经浸润原发灶 vs 低：找低剂量基因（高浸润组更低）
unpaired_low_in_high <- function(high_ids, low_ids, label) {
  high_ids <- intersect(high_ids, colnames(mat))
  low_ids <- intersect(low_ids, colnames(mat))
  if (length(high_ids) < 2 || length(low_ids) < 2) {
    log_msg(label, ": too few samples high=", length(high_ids), " low=", length(low_ids))
    return(data.frame(gene = rownames(mat), n_high = length(high_ids), n_low = length(low_ids),
                      mean_high = NA_real_, mean_low = NA_real_, log2FC = NA_real_,
                      FC = NA_real_, pvalue = NA_real_, padj = NA_real_,
                      stringsAsFactors = FALSE))
  }
  h <- mat[, high_ids, drop = FALSE]
  l <- mat[, low_ids, drop = FALSE]
  mean_h <- matrixStats::rowMeans2(h, na.rm = TRUE)
  mean_l <- matrixStats::rowMeans2(l, na.rm = TRUE)
  log2fc <- mean_l - mean_h  # 正值 = 高浸润组更低 = 低剂量伴随神经浸润
  pval <- apply(mat, 1, function(x) {
    a <- x[high_ids]; b <- x[low_ids]
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a) < 2 || length(b) < 2) return(NA_real_)
    if (stats::sd(a) == 0 && stats::sd(b) == 0) {
      return(if (abs(mean(a) - mean(b)) < 1e-12) 1 else NA_real_)
    }
    tryCatch(stats::t.test(a, b)$p.value, error = function(e) NA_real_)
  })
  data.frame(
    gene = rownames(mat),
    n_high = length(high_ids),
    n_low = length(low_ids),
    mean_high = mean_h,
    mean_low = mean_l,
    log2FC = log2fc,
    FC = 2^log2fc,
    pvalue = pval,
    padj = p.adjust(pval, method = "BH"),
    stringsAsFactors = FALSE
  )
}

neural_specs <- list(
  axon_guidance = list(mark = "mark_axon", high = "High_axon_guidance", score = "score_axon", genes = gene_sets$axon_guidance),
  schwann = list(mark = "mark_schwann", high = "High_Schwann", score = "score_schwann", genes = gene_sets$schwann),
  neurotrophic = list(mark = "mark_neurotrophic", high = "High_neurotrophic", score = "score_neurotrophic", genes = gene_sets$neurotrophic)
)

log_msg("==== Q3: genes low in high neural-invasion primaries (3 signatures) ====")
neural_de <- list()
for (nm in names(neural_specs)) {
  sp <- neural_specs[[nm]]
  ok <- !is.na(prim_ann[[sp$mark]])
  high_ids <- prim_ann$SAMPLE_ID[ok & prim_ann[[sp$mark]] == sp$high]
  low_ids <- prim_ann$SAMPLE_ID[ok & prim_ann[[sp$mark]] != sp$high]
  de_n <- unpaired_low_in_high(high_ids, low_ids, nm)
  neural_de[[nm]] <- de_n
  ndir <- file.path(result_dir, "05_neural_invasion", nm)
  save_tbl(de_n, file.path(ndir, paste0(nm, "_all_genes_high_vs_low_PNI.csv")))
  for (fcnm in names(fc_cutoffs)) {
    hit <- de_n[is.finite(de_n$FC) & de_n$FC >= fc_cutoffs[[fcnm]] & de_n$log2FC > 0 &
                  is.finite(de_n$pvalue) & de_n$pvalue < p_cutoff, , drop = FALSE]
    hit <- hit[order(hit$pvalue, -hit$FC), ]
    save_tbl(hit, file.path(ndir, paste0(fcnm, "_low_in_high_", nm, ".csv")))
    in_set <- hit[hit$gene %in% sp$genes, , drop = FALSE]
    save_tbl(in_set, file.path(ndir, paste0(fcnm, "_", nm, "_geneset_only.csv")))
    log_msg("Q3 ", nm, " ", fcnm, ": genome-wide ", nrow(hit),
            " | in gene set ", nrow(in_set))
    plot_volcano(
      de_n, hit$gene,
      paste0("Primary high vs low ", nm, " | ", fcnm),
      file.path(ndir, paste0(fcnm, "_volcano")),
      fc_line = fc_cutoffs[[fcnm]]
    )
    run_ora(hit$gene, de_n, file.path(ndir, paste0(fcnm, "_ORA")),
            paste("PNI", nm, fcnm), paste0(fcnm, "_", nm))
  }
  gs_tbl <- data.frame(
    gene = unique(gene_sets_full[[nm]]),
    in_matrix = unique(gene_sets_full[[nm]]) %in% rownames(mat),
    stringsAsFactors = FALSE
  )
  if (nrow(gs_tbl) > 0) {
    gs_tbl$mean_high_PNI <- de_n$mean_high[match(gs_tbl$gene, de_n$gene)]
    gs_tbl$mean_low_PNI <- de_n$mean_low[match(gs_tbl$gene, de_n$gene)]
    gs_tbl$FC_low_over_high <- de_n$FC[match(gs_tbl$gene, de_n$gene)]
    gs_tbl$pvalue <- de_n$pvalue[match(gs_tbl$gene, de_n$gene)]
    gs_tbl$log2FC <- de_n$log2FC[match(gs_tbl$gene, de_n$gene)]
  }
  save_tbl(gs_tbl, file.path(ndir, paste0(nm, "_geneset_membership.csv")))
}

# -----------------------------------------------------------------------------
# 7. Q4：三种神经浸润 vs 肺转移
# -----------------------------------------------------------------------------
log_msg("==== Q4: three neural-invasion scores vs lung metastasis ====")
q4_dir <- file.path(result_dir, "06_neural_vs_lung")
dir.create(q4_dir, recursive = TRUE, showWarnings = FALSE)

score_long <- prim_ann %>%
  dplyr::select("PATIENT_ID", "SAMPLE_ID", "has_paired_lung", "ever_lung",
                "score_axon", "score_schwann", "score_neurotrophic") %>%
  tidyr::pivot_longer(
    cols = c("score_axon", "score_schwann", "score_neurotrophic"),
    names_to = "signature", values_to = "score"
  )
score_long$signature <- dplyr::recode(
  score_long$signature,
  score_axon = "axon_guidance",
  score_schwann = "schwann",
  score_neurotrophic = "neurotrophic"
)
save_tbl(score_long, file.path(q4_dir, "primary_scores_long.csv"))

assoc_one <- function(nm, grp) {
  sc <- switch(nm,
               axon_guidance = prim_ann$score_axon,
               schwann = prim_ann$score_schwann,
               neurotrophic = prim_ann$score_neurotrophic)
  ok <- is.finite(sc) & !is.na(grp)
  wt <- tryCatch(stats::wilcox.test(sc[ok & grp], sc[ok & !grp])$p.value, error = function(e) NA_real_)
  tt <- tryCatch(stats::t.test(sc[ok & grp], sc[ok & !grp])$p.value, error = function(e) NA_real_)
  data.frame(
    signature = nm,
    n_yes = sum(ok & grp),
    n_no = sum(ok & !grp),
    mean_score_yes = mean(sc[ok & grp], na.rm = TRUE),
    mean_score_no = mean(sc[ok & !grp], na.rm = TRUE),
    wilcoxon_p = wt,
    ttest_p = tt,
    stringsAsFactors = FALSE
  )
}

assoc_paired <- dplyr::bind_rows(lapply(
  c("axon_guidance", "schwann", "neurotrophic"),
  function(nm) {
    out <- assoc_one(nm, prim_ann$has_paired_lung)
    names(out)[names(out) == "n_yes"] <- "n_paired_lung_primary"
    names(out)[names(out) == "n_no"] <- "n_other_primary"
    names(out)[names(out) == "mean_score_yes"] <- "mean_score_paired_lung"
    names(out)[names(out) == "mean_score_no"] <- "mean_score_other"
    out
  }
))
save_tbl(assoc_paired, file.path(q4_dir, "neural_score_vs_paired_lung_wilcoxon.csv"))
log_msg("Neural score vs paired-lung primary: ",
        paste(assoc_paired$signature, "p=", signif(assoc_paired$wilcoxon_p, 3), collapse = "; "))

# 临床是否出现过肺转移（不要求该患者同时有原发测序+肺转移测序）
assoc_ever <- dplyr::bind_rows(lapply(
  c("axon_guidance", "schwann", "neurotrophic"),
  function(nm) {
    out <- assoc_one(nm, prim_ann$ever_lung)
    names(out)[names(out) == "n_yes"] <- "n_ever_lung"
    names(out)[names(out) == "n_no"] <- "n_never_lung"
    names(out)[names(out) == "mean_score_yes"] <- "mean_score_ever_lung"
    names(out)[names(out) == "mean_score_no"] <- "mean_score_never_lung"
    out
  }
))
save_tbl(assoc_ever, file.path(q4_dir, "neural_score_vs_ever_lung_wilcoxon.csv"))

pbox <- ggplot2::ggplot(
  score_long[!is.na(score_long$has_paired_lung), ],
  ggplot2::aes(x = ifelse(.data$has_paired_lung, "Paired lung met", "No paired lung met"),
               y = .data$score, fill = ifelse(.data$has_paired_lung, "Paired lung met", "No paired lung met"))
) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  ggplot2::geom_jitter(width = 0.15, size = 1.4, alpha = 0.8) +
  ggplot2::facet_wrap(~signature, scales = "free_y") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)) +
  ggplot2::labs(
    title = "MSK 2018 primary neural-invasion scores vs paired lung metastasis",
    x = NULL,
    y = if (grepl("^CNA", assay_mode)) "Signature score (mean CNA of gene set)" else "Signature score (mean of gene set)"
  )
save_gg(pbox, file.path(q4_dir, "boxplot_neural_score_vs_paired_lung"), width = 10, height = 5)

pbox2 <- ggplot2::ggplot(
  score_long[!is.na(score_long$ever_lung), ],
  ggplot2::aes(x = ifelse(.data$ever_lung, "Ever lung met (clinical)", "No lung met recorded"),
               y = .data$score, fill = ifelse(.data$ever_lung, "Ever lung met (clinical)", "No lung met recorded"))
) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  ggplot2::geom_jitter(width = 0.15, size = 0.8, alpha = 0.35) +
  ggplot2::facet_wrap(~signature, scales = "free_y") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)) +
  ggplot2::labs(
    title = "MSK 2018 primary neural-invasion scores vs any recorded lung metastasis",
    x = NULL,
    y = if (grepl("^CNA", assay_mode)) "Signature score (mean CNA of gene set)" else "Signature score (mean of gene set)"
  )
save_gg(pbox2, file.path(q4_dir, "boxplot_neural_score_vs_ever_lung"), width = 10, height = 5)

for (nm in names(neural_specs)) {
  sp <- neural_specs[[nm]]
  keep <- !is.na(prim_ann[[sp$mark]])
  tab <- table(
    High = prim_ann[[sp$mark]][keep] == sp$high,
    PairedLung = prim_ann$has_paired_lung[keep]
  )
  ft <- tryCatch(stats::fisher.test(tab), error = function(e) NULL)
  cap <- data.frame(
    signature = nm,
    comparison = "paired_lung",
    fisher_p = if (is.null(ft)) NA_real_ else ft$p.value,
    odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate),
    stringsAsFactors = FALSE
  )
  save_tbl(as.data.frame.matrix(tab), file.path(q4_dir, paste0(nm, "_highLow_vs_paired_lung_counts.csv")))
  save_tbl(cap, file.path(q4_dir, paste0(nm, "_fisher_highPNI_vs_paired_lung.csv")))
  log_msg("Fisher high ", nm, " vs paired lung: p=",
          if (is.null(ft)) "NA" else signif(ft$p.value, 3))

  tab2 <- table(
    High = prim_ann[[sp$mark]][keep] == sp$high,
    EverLung = prim_ann$ever_lung[keep]
  )
  ft2 <- tryCatch(stats::fisher.test(tab2), error = function(e) NULL)
  cap2 <- data.frame(
    signature = nm,
    comparison = "ever_lung",
    fisher_p = if (is.null(ft2)) NA_real_ else ft2$p.value,
    odds_ratio = if (is.null(ft2)) NA_real_ else unname(ft2$estimate),
    stringsAsFactors = FALSE
  )
  save_tbl(as.data.frame.matrix(tab2), file.path(q4_dir, paste0(nm, "_highLow_vs_ever_lung_counts.csv")))
  save_tbl(cap2, file.path(q4_dir, paste0(nm, "_fisher_highPNI_vs_ever_lung.csv")))
}

for (fcnm in names(fc_cutoffs)) {
  g_lung <- lung_sets[[fcnm]]$gene
  ov_rows <- lapply(names(neural_specs), function(nm) {
    de_n <- neural_de[[nm]]
    g_pni <- de_n$gene[is.finite(de_n$FC) & de_n$FC >= fc_cutoffs[[fcnm]] &
                         de_n$log2FC > 0 & is.finite(de_n$pvalue) & de_n$pvalue < p_cutoff]
    inter <- intersect(g_lung, g_pni)
    in_set <- intersect(inter, neural_specs[[nm]]$genes)
    data.frame(
      signature = nm,
      n_lung_low = length(g_lung),
      n_pni_low = length(g_pni),
      n_overlap = length(inter),
      n_overlap_in_signature_geneset = length(in_set),
      overlap_genes = paste(inter, collapse = ";"),
      overlap_in_geneset = paste(in_set, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  ov <- dplyr::bind_rows(ov_rows)
  save_tbl(ov, file.path(q4_dir, paste0(fcnm, "_overlap_lung_low_AND_PNI_low.csv")))
  log_msg(fcnm, " overlap lung-low & PNI-low: ",
          paste(ov$signature, ov$n_overlap, sep = "=", collapse = ", "))
}

if (nrow(pairs_lung) >= 3) {
  pe <- pair_expr(pairs_lung)
  sig_fc <- lapply(names(gene_sets), function(nm) {
    g <- intersect(gene_sets[[nm]], rownames(pe$primary))
    if (length(g) < 3) return(rep(NA_real_, nrow(pairs_lung)))
    d <- pe$met[g, , drop = FALSE] - pe$primary[g, , drop = FALSE]
    matrixStats::colMeans2(d, na.rm = TRUE)
  })
  names(sig_fc) <- paste0("paired_log2FC_", names(gene_sets))
  cor_df <- cbind(pairs_lung, as.data.frame(sig_fc))
  cor_df$score_axon <- prim_ann$score_axon[match(cor_df$patient_id, prim_ann$PATIENT_ID)]
  cor_df$score_schwann <- prim_ann$score_schwann[match(cor_df$patient_id, prim_ann$PATIENT_ID)]
  cor_df$score_neurotrophic <- prim_ann$score_neurotrophic[match(cor_df$patient_id, prim_ann$PATIENT_ID)]
  save_tbl(cor_df, file.path(q4_dir, "paired_lung_patient_neural_score_vs_geneset_FC.csv"))
  cor_tests <- data.frame(
    signature = c("axon_guidance", "schwann", "neurotrophic"),
    pearson_r = c(
      suppressWarnings(stats::cor(cor_df$score_axon, cor_df$paired_log2FC_axon_guidance, use = "complete.obs")),
      suppressWarnings(stats::cor(cor_df$score_schwann, cor_df$paired_log2FC_schwann, use = "complete.obs")),
      suppressWarnings(stats::cor(cor_df$score_neurotrophic, cor_df$paired_log2FC_neurotrophic, use = "complete.obs"))
    ),
    stringsAsFactors = FALSE
  )
  cor_tests$pearson_p <- c(
    tryCatch(stats::cor.test(cor_df$score_axon, cor_df$paired_log2FC_axon_guidance)$p.value, error = function(e) NA_real_),
    tryCatch(stats::cor.test(cor_df$score_schwann, cor_df$paired_log2FC_schwann)$p.value, error = function(e) NA_real_),
    tryCatch(stats::cor.test(cor_df$score_neurotrophic, cor_df$paired_log2FC_neurotrophic)$p.value, error = function(e) NA_real_)
  )
  save_tbl(cor_tests, file.path(q4_dir, "correlation_primary_neural_score_vs_paired_lung_FC.csv"))
}

# -----------------------------------------------------------------------------
# 8. 总览
# -----------------------------------------------------------------------------
summary_n <- data.frame(
  item = c(
    "assay_mode", "n_genes", "n_primary", "n_metastasis",
    "paired_primary_lung", "paired_primary_bone",
    paste0("lung_", names(fc_cutoffs), "_n_genes"),
    paste0("bone_", names(fc_cutoffs), "_n_genes")
  ),
  value = c(
    assay_mode, as.character(nrow(mat)),
    as.character(sum(clin$SAMPLE_TYPE == "Primary" & clin$has_assay)),
    as.character(sum(clin$SAMPLE_TYPE == "Metastatic" & clin$has_assay)),
    as.character(nrow(pairs_lung)), as.character(nrow(pairs_bone)),
    as.character(vapply(lung_sets, nrow, integer(1))),
    as.character(vapply(bone_sets, nrow, integer(1)))
  ),
  stringsAsFactors = FALSE
)
save_tbl(summary_n, file.path(result_dir, "00_logs", "analysis_summary_counts.csv"))
log_msg("Done. Results in: ", result_dir)
log_msg("Read 01_pairing first to confirm 1 primary : 1 organ-met per patient.")
if (grepl("^CNA", assay_mode)) {
  log_msg("NOTE: 当前是拷贝数剂量，不是 RNA。FC = 2^(CNA_met - CNA_primary)。")
}
if (expr_is_zscore) {
  log_msg("NOTE: 当前矩阵是 z-score。FC = 2^(z_met - z_primary)。")
}
