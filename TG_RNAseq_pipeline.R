#!/usr/bin/env Rscript
# =============================================================================
# TG BRCA 细胞 RNA-seq 分析流程（下调方向，六种设计）
# 组别：NTC_rep0、NTC_rep1、TG_sh1、TG_sh5
# 六种比较（各自单独作图；仅设计 2 合并两个 NTC）：
#   1) 单独 1-vs-1：TG_sh1 vs NTC_rep0，TG_sh5 vs NTC_rep0，
#                  TG_sh1 vs NTC_rep1，TG_sh5 vs NTC_rep1
#   2) (TG_sh1 + TG_sh5)/2 vs NTC组均值(NTC_rep0, NTC_rep1)
#   3) 共同下调：TG_sh1 vs NTC_rep0 与 TG_sh5 vs NTC_rep0 的交集
#   4) 共同下调：TG_sh1 vs NTC_rep1 与 TG_sh5 vs NTC_rep1 的交集
#   5) (TG_sh1 + TG_sh5)/2 vs NTC_rep0
#   6) (TG_sh1 + TG_sh5)/2 vs NTC_rep1
# 预处理：过滤低表达 + 标准化消除技术偏差
# 子集策略：
#   A) 下调 FC <= 1 / 0.8 / 0.6 / 0.5
#   B) 下调排名 top 50 / 75 / 100 / 150 / 200 / 250 / 300（log2FC 最负）
# 每个比较 × 每个 FC 阈值 × 每个 topN 都必须出图：
#   差异基因表/柱状图、火山图、热图、GO图、通路富集图、KEGG图、GSEA图
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")
options(clusterProfiler.download.method = "auto")

# -----------------------------------------------------------------------------
# 0. 依赖包（必需包失败才中止；pathview / ReactomePA 缺失时跳过对应图）
# -----------------------------------------------------------------------------
cran_required <- c(
  "readxl", "writexl", "dplyr", "tidyr", "tibble", "stringr", "ggplot2",
  "ggrepel", "pheatmap", "RColorBrewer", "matrixStats", "cowplot",
  "ggridges", "ggnewscale", "igraph"
)
cran_optional <- c("ggvenn")
bioc_required <- c(
  "DESeq2", "edgeR", "limma", "clusterProfiler", "org.Hs.eg.db",
  "enrichplot", "DOSE", "AnnotationDbi", "fgsea", "msigdbr",
  "SummarizedExperiment"
)
bioc_optional <- c("ReactomePA", "pathview")

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
  if (length(still) > 0) message("可选包未安装，相关分析将跳过: ", paste(still, collapse = ", "))
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
# 1. 路径与分析参数
# -----------------------------------------------------------------------------
resolve_project_dir <- function() {
  env_dir <- Sys.getenv("TG_RNASEQ_DIR", unset = "")
  candidates <- c(
    env_dir,
    "E:/R/TG_BRCA/TG",
    "E:\\R\\TG_BRCA\\TG",
    getwd()
  )
  candidates <- unique(candidates[nzchar(candidates)])
  for (d in candidates) {
    if (dir.exists(d) && (
      file.exists(file.path(d, "genes.read_group_tracking")) ||
      file.exists(file.path(d, "genes.fpkm_tracking")) ||
      file.exists(file.path(d, "genes.count_tracking")) ||
      file.exists(file.path(d, "read_groups.info"))
    )) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

project_dir <- resolve_project_dir()
result_dir  <- file.path(project_dir, "results")
log_dir     <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

padj_cutoff <- 0.05
# 下调倍数：treat/control <= 1 / 0.8 / 0.6 / 0.5
fc_cutoffs <- c("FC_1" = 1, "FC_0.8" = 0.8, "FC_0.6" = 0.6, "FC_0.5" = 0.5)
top_ns     <- c(50, 75, 100, 150, 200, 250, 300)
down_fill  <- "#1D4E89"

# -----------------------------------------------------------------------------
# 2. 样本名识别
# -----------------------------------------------------------------------------
classify_sample <- function(name) {
  n <- toupper(gsub("[^A-Za-z0-9]", "", name))
  if (grepl("SH5|SHRNA5|TGSH5|SHTG5", n)) return("TG_sh5")
  if (grepl("SH1|SHRNA1|TGSH1|SHTG1", n)) return("TG_sh1")
  if (grepl("SHTG|TGS", n) && grepl("5", n)) return("TG_sh5")
  if (grepl("SHTG|TGS", n)) return("TG_sh1")
  if (grepl("NTC|SHNC|NEGCTRL|CTRL|CONTROL", n)) return("NTC")
  if (grepl("^NC[0-9]*$", n)) return("NTC")
  NA_character_
}

# 两个 NTC 样品单独标记为 NTC_rep0 / NTC_rep1，不在 1-vs-1 比较里合并
add_ntc_ids <- function(sample_info) {
  sample_info$ntc_id <- NA_character_
  ntc <- which(sample_info$group == "NTC")
  if (length(ntc) == 0) return(sample_info)
  labs <- as.character(sample_info$sample[ntc])
  ids <- rep(NA_character_, length(labs))
  ids[grepl("rep1|[_-]1$", labs, ignore.case = TRUE)] <- "NTC_rep1"
  ids[grepl("rep0|[_-]0$", labs, ignore.case = TRUE)] <- "NTC_rep0"
  if (any(is.na(ids)) && length(ntc) == 2) {
    ids[order(labs)] <- c("NTC_rep0", "NTC_rep1")
  }
  if (length(ntc) == 1 && is.na(ids[1])) ids[1] <- "NTC_rep0"
  sample_info$ntc_id[ntc] <- ids
  sample_info
}

find_sample <- function(sample_info, group, ntc_id = NULL) {
  if (!is.null(ntc_id)) {
    hit <- sample_info$sample[sample_info$group == "NTC" & sample_info$ntc_id == ntc_id]
  } else {
    hit <- sample_info$sample[sample_info$group == group]
  }
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

guess_gene_column <- function(df) {
  nms <- names(df)
  low <- tolower(nms)
  keys <- c("gene_short_name", "gene_name", "gene_symbol", "symbol",
            "geneid", "gene_id", "tracking_id", "id", "genes")
  for (k in keys) {
    hit <- which(low == k | grepl(paste0("^", k, "$"), low))
    if (length(hit) > 0) return(nms[hit[1]])
  }
  non_num <- nms[vapply(df, function(x) !is.numeric(x) && !is.integer(x), logical(1))]
  if (length(non_num) > 0) return(non_num[1])
  nms[1]
}

pick_official_symbol <- function(x) {
  x <- trimws(as.character(x))
  if (length(x) != 1 || is.na(x) || x %in% c("", "-", ".", "NA")) return(NA_character_)
  parts <- unlist(strsplit(x, "[,;|/]+"))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts) & !parts %in% c("-", ".", "NA")]
  if (length(parts) == 0) return(NA_character_)
  is_fusion <- vapply(parts, function(t) {
    bits <- strsplit(t, "-", fixed = TRUE)[[1]]
    length(bits) == 2 && (bits[1] %in% parts || bits[2] %in% parts)
  }, logical(1))
  if (any(!is_fusion)) parts <- parts[!is_fusion]
  score <- vapply(parts, function(s) {
    if (grepl("^(XLOC|TCONS|CUFF)_", s, ignore.case = TRUE)) return(0)
    if (grepl("^MIR[0-9]", s, ignore.case = TRUE)) return(1)
    if (grepl("^LOC[0-9]+$", s, ignore.case = TRUE)) return(2)
    3
  }, numeric(1))
  parts[which.max(score)]
}

clean_gene_names <- function(symbols, tracking_ids = NULL, nearest_ref = NULL) {
  out <- vapply(symbols, pick_official_symbol, character(1), USE.NAMES = FALSE)
  if (!is.null(nearest_ref)) {
    need <- is.na(out) | grepl("^(XLOC|TCONS|CUFF)_", out, ignore.case = TRUE)
    ref <- vapply(nearest_ref, pick_official_symbol, character(1), USE.NAMES = FALSE)
    use_ref <- need & !is.na(ref) & !grepl("^(XLOC|TCONS|CUFF|NM_|NR_|ENST)", ref, ignore.case = TRUE)
    out[use_ref] <- ref[use_ref]
  }
  if (!is.null(tracking_ids)) {
    still <- is.na(out) | !nzchar(out)
    out[still] <- as.character(tracking_ids[still])
  }
  out
}

apply_gene_labels <- function(mat, symbols, tracking_ids = NULL, nearest_ref = NULL) {
  raw <- as.character(symbols)
  genes <- clean_gene_names(raw, tracking_ids = tracking_ids, nearest_ref = nearest_ref)
  n_fused <- sum(grepl("[,;|/]", raw), na.rm = TRUE)
  mat <- collapse_by_gene(mat, genes)
  n_xloc <- sum(grepl("^(XLOC|TCONS|CUFF)_", rownames(mat), ignore.case = TRUE))
  n_comma <- sum(grepl(",", rownames(mat), fixed = TRUE))
  log_msg(
    "Gene labels: ", nrow(mat), " unique after collapsing; ",
    "Cufflinks fused/comma names cleaned: ", n_fused, "; ",
    "remaining XLOC (novel/unannotated): ", n_xloc, "; ",
    "remaining comma names: ", n_comma
  )
  mat
}

collapse_by_gene <- function(mat, genes) {
  genes[is.na(genes) | genes == "" | genes == "-"] <- NA
  keep <- !is.na(genes)
  mat <- mat[keep, , drop = FALSE]
  genes <- genes[keep]
  if (any(duplicated(genes))) {
    means <- rowMeans(mat, na.rm = TRUE)
    ord <- order(means, decreasing = TRUE)
    mat <- mat[ord, , drop = FALSE]
    genes <- genes[ord]
    keep_u <- !duplicated(genes)
    mat <- mat[keep_u, , drop = FALSE]
    genes <- genes[keep_u]
  }
  rownames(mat) <- genes
  mat
}

# -----------------------------------------------------------------------------
# 3. 读入表达矩阵（Cuffdiff tracking；不再使用已删除的 Excel）
# -----------------------------------------------------------------------------
read_read_group_tracking <- function(path) {
  rg <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  need <- c("tracking_id", "condition", "replicate")
  if (!all(need %in% names(rg))) {
    log_msg("genes.read_group_tracking 缺少必要列: ", paste(setdiff(need, names(rg)), collapse = ", "))
    return(NULL)
  }
  value_col <- if ("raw_frags" %in% names(rg)) {
    "raw_frags"
  } else if ("external_scaled_frags" %in% names(rg)) {
    "external_scaled_frags"
  } else if ("FPKM" %in% names(rg)) {
    "FPKM"
  } else {
    return(NULL)
  }
  log_msg("Cuffdiff conditions: ", paste(unique(as.character(rg$condition)), collapse = ", "))
  log_msg("Cuffdiff value column: ", value_col)
  rg$sample <- paste(rg$condition, rg$replicate, sep = "_rep")
  rg$group <- vapply(as.character(rg$condition), classify_sample, character(1))
  if (all(is.na(rg$group))) {
    rg$group <- vapply(rg$sample, classify_sample, character(1))
  }
  if (all(is.na(rg$group))) {
    log_msg("无法把 Cuffdiff condition 映射到 NTC/TG_sh1/TG_sh5: ",
            paste(unique(as.character(rg$condition)), collapse = ", "))
    return(NULL)
  }
  dropped <- unique(as.character(rg$condition[is.na(rg$group)]))
  if (length(dropped) > 0) log_msg("Unmapped conditions dropped: ", paste(dropped, collapse = ", "))
  rg <- rg[!is.na(rg$group), , drop = FALSE]
  if (nrow(rg) == 0) return(NULL)
  wide <- tidyr::pivot_wider(
    rg[, c("tracking_id", "sample", value_col)],
    names_from = "sample",
    values_from = value_col,
    values_fn = mean
  )
  gene_map <- NULL
  fpkm_file <- file.path(dirname(path), "genes.fpkm_tracking")
  if (file.exists(fpkm_file)) {
    fp <- utils::read.delim(fpkm_file, check.names = FALSE, stringsAsFactors = FALSE)
    if (all(c("tracking_id", "gene_short_name") %in% names(fp))) {
      gene_map <- fp
    }
  }
  genes <- wide$tracking_id
  nearest <- NULL
  if (!is.null(gene_map)) {
    hit <- match(wide$tracking_id, gene_map$tracking_id)
    genes <- gene_map$gene_short_name[hit]
    if ("nearest_ref_id" %in% names(gene_map)) {
      nearest <- gene_map$nearest_ref_id[hit]
    }
  }
  mat <- as.matrix(wide[, setdiff(names(wide), "tracking_id"), drop = FALSE])
  storage.mode(mat) <- "double"
  mat[is.na(mat)] <- 0
  mat <- apply_gene_labels(mat, genes, tracking_ids = wide$tracking_id, nearest_ref = nearest)
  sample_info <- unique(rg[, c("sample", "group")])
  sample_info <- sample_info[match(colnames(mat), sample_info$sample), ]
  log_msg("Sample mapping: ", paste(paste0(sample_info$sample, "=", sample_info$group), collapse = "; "))
  list(mat = mat, sample_info = sample_info, source = basename(path), value_col = value_col)
}

read_tracking_matrix <- function(path, value_pattern) {
  tr <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  gene_col <- if ("gene_short_name" %in% names(tr)) "gene_short_name" else if ("gene_id" %in% names(tr)) "gene_id" else names(tr)[1]
  val_cols <- grep(value_pattern, names(tr), value = TRUE, ignore.case = TRUE)
  val_cols <- val_cols[!grepl("variance|conf|status|dispersion|uncertainty", val_cols, ignore.case = TRUE)]
  if (length(val_cols) == 0) return(NULL)
  groups <- vapply(val_cols, classify_sample, character(1))
  if (all(is.na(groups))) {
    rg_info <- file.path(dirname(path), "read_groups.info")
    if (file.exists(rg_info)) {
      info <- utils::read.delim(rg_info, check.names = FALSE, stringsAsFactors = FALSE)
      cond <- if ("condition" %in% names(info)) info$condition else info[[2]]
      if (length(cond) == length(val_cols)) {
        groups <- vapply(as.character(cond), classify_sample, character(1))
        names(val_cols) <- paste(cond, seq_along(cond), sep = "_rep")
        colnames_keep <- names(val_cols)
      }
    }
  }
  keep <- !is.na(groups)
  if (sum(keep) < 2) return(NULL)
  mat <- as.matrix(tr[, val_cols[keep], drop = FALSE])
  storage.mode(mat) <- "double"
  mat[is.na(mat)] <- 0
  if (exists("colnames_keep")) colnames(mat) <- colnames_keep[keep]
  tid <- if ("tracking_id" %in% names(tr)) tr$tracking_id else NULL
  nearest <- if ("nearest_ref_id" %in% names(tr)) tr$nearest_ref_id else NULL
  mat <- apply_gene_labels(mat, tr[[gene_col]], tracking_ids = tid, nearest_ref = nearest)
  list(
    mat = mat,
    sample_info = data.frame(sample = colnames(mat), group = unname(groups[keep]), stringsAsFactors = FALSE),
    source = basename(path)
  )
}

load_expression <- function(project_dir) {
  rg <- file.path(project_dir, "genes.read_group_tracking")
  if (file.exists(rg)) {
    log_msg("Reading Cuffdiff replicate file: genes.read_group_tracking")
    obj <- tryCatch(read_read_group_tracking(rg), error = function(e) {
      log_msg("read_group_tracking import failed: ", e$message)
      NULL
    })
    if (!is.null(obj)) return(obj)
  }
  ct <- file.path(project_dir, "genes.count_tracking")
  if (file.exists(ct)) {
    log_msg("Reading Cuffdiff count file: genes.count_tracking")
    obj <- tryCatch(read_tracking_matrix(ct, "_count$|^q[0-9]+_count$"), error = function(e) {
      log_msg("count_tracking import failed: ", e$message)
      NULL
    })
    if (!is.null(obj)) return(obj)
  }
  fp <- file.path(project_dir, "genes.fpkm_tracking")
  if (file.exists(fp)) {
    log_msg("Reading Cuffdiff FPKM file: genes.fpkm_tracking")
    obj <- tryCatch(read_tracking_matrix(fp, "_FPKM$|^q[0-9]+_FPKM$"), error = function(e) {
      log_msg("fpkm_tracking import failed: ", e$message)
      NULL
    })
    if (!is.null(obj)) return(obj)
  }
  stop("未找到 Cuffdiff 表达文件。请确认目录中有 genes.read_group_tracking / genes.count_tracking / genes.fpkm_tracking: ", project_dir)
}

detect_value_type <- function(mat) {
  x <- as.numeric(mat)
  x <- x[is.finite(x)]
  frac_int <- mean(abs(x - round(x)) < 1e-6)
  if (frac_int > 0.85 && stats::quantile(x, 0.95, na.rm = TRUE) > 50) {
    return("counts")
  }
  "fpkm"
}

# -----------------------------------------------------------------------------
# 4. 过滤低表达 + 标准化
# -----------------------------------------------------------------------------
filter_low_expression <- function(mat, sample_info, value_type) {
  min_n <- max(2, min(table(sample_info$group)))
  if (value_type == "counts") {
    keep <- rowSums(mat >= 10, na.rm = TRUE) >= min_n
  } else {
    keep <- rowSums(mat > 1, na.rm = TRUE) >= min_n
  }
  if (sum(keep) < 200) {
    keep <- rowSums(mat > 0, na.rm = TRUE) >= min_n
    log_msg("Strict filter left too few genes; fallback to expressed-in-", min_n, "-samples")
  }
  log_msg("Low-expression filter: keep ", sum(keep), " / ", nrow(mat), " genes")
  mat[keep, , drop = FALSE]
}

normalize_expression <- function(mat, sample_info, value_type) {
  group <- factor(sample_info$group, levels = c("NTC", "TG_sh1", "TG_sh5"))
  group <- droplevels(group)
  if (value_type == "counts") {
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(mat),
      colData = data.frame(row.names = colnames(mat), group = group),
      design = ~ group
    )
    dds <- DESeq2::estimateSizeFactors(dds)
    log_msg("DESeq2 size-factor normalization")
    norm_counts <- DESeq2::counts(dds, normalized = TRUE)
    log_mat <- log2(norm_counts + 1)
    vsd <- tryCatch({
      SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))
    }, error = function(e) {
      log_msg("vst failed, use log2(norm+1): ", e$message)
      log_mat
    })
    list(dds = dds, log_mat = log_mat, heat_mat = vsd, group = group, method = "DESeq2")
  } else {
    dge <- edgeR::DGEList(counts = pmax(mat, 0), group = group)
    keep <- edgeR::filterByExpr(dge, group = group)
    if (sum(keep) >= 200) {
      mat <- mat[keep, , drop = FALSE]
      dge <- dge[keep, , keep.lib.sizes = FALSE]
    }
    log_msg("Quantile normalize log2(FPKM+1)")
    log_mat <- limma::normalizeBetweenArrays(log2(pmax(mat, 0) + 1), method = "quantile")
    list(dds = NULL, log_mat = log_mat, heat_mat = log_mat, group = group, method = "quantile_logFPKM")
  }
}

# -----------------------------------------------------------------------------
# 5. 差异分析
#   1-vs-1 比较：标准化后的 log 值直接相减，无 P 值
#   合并比较：两个 knockdown 等权平均 vs 两个 NTC 的组均值
# -----------------------------------------------------------------------------
pairwise_de <- function(log_mat, treat_sample, ntc_sample, comp_name) {
  if (is.na(treat_sample) || is.na(ntc_sample)) return(NULL)
  if (!all(c(treat_sample, ntc_sample) %in% colnames(log_mat))) return(NULL)
  log_msg(comp_name, " : ", treat_sample, " vs ", ntc_sample, " (1-vs-1, FC only)")
  log2FC <- log_mat[, treat_sample] - log_mat[, ntc_sample]
  ave <- (log_mat[, treat_sample] + log_mat[, ntc_sample]) / 2
  data.frame(
    gene = rownames(log_mat),
    log2FC = as.numeric(log2FC),
    AveExpr = as.numeric(ave),
    pvalue = NA_real_,
    padj = NA_real_,
    treat_sample = treat_sample,
    ntc_sample = ntc_sample,
    stringsAsFactors = FALSE
  )
}

mean_kd_vs_ntc_de <- function(log_mat, sample_info) {
  ntc <- sample_info$sample[sample_info$group == "NTC"]
  sh1 <- find_sample(sample_info, "TG_sh1")
  sh5 <- find_sample(sample_info, "TG_sh5")
  if (length(ntc) < 1 || is.na(sh1) || is.na(sh5)) return(NULL)
  log_msg("TGsh_mean_vs_NTC : mean(", sh1, ", ", sh5, ") vs mean(", paste(ntc, collapse = ", "), ")")
  ntc_mean <- rowMeans(log_mat[, ntc, drop = FALSE])
  sh_mean <- (log_mat[, sh1] + log_mat[, sh5]) / 2
  data.frame(
    gene = rownames(log_mat),
    log2FC = as.numeric(sh_mean - ntc_mean),
    AveExpr = as.numeric((sh_mean + ntc_mean) / 2),
    pvalue = NA_real_,
    padj = NA_real_,
    stringsAsFactors = FALSE
  )
}

mean_kd_vs_one_ntc_de <- function(log_mat, sample_info, ntc_id, comp_name) {
  sh1 <- find_sample(sample_info, "TG_sh1")
  sh5 <- find_sample(sample_info, "TG_sh5")
  ntc <- find_sample(sample_info, "NTC", ntc_id)
  if (is.na(sh1) || is.na(sh5) || is.na(ntc)) return(NULL)
  if (!all(c(sh1, sh5, ntc) %in% colnames(log_mat))) return(NULL)
  log_msg(comp_name, " : mean(", sh1, ", ", sh5, ") vs ", ntc, " (KD mean vs one NTC, FC only)")
  sh_mean <- (log_mat[, sh1] + log_mat[, sh5]) / 2
  ntc_val <- log_mat[, ntc]
  data.frame(
    gene = rownames(log_mat),
    log2FC = as.numeric(sh_mean - ntc_val),
    AveExpr = as.numeric((sh_mean + ntc_val) / 2),
    pvalue = NA_real_,
    padj = NA_real_,
    treat_sample = paste0("mean(", sh1, ",", sh5, ")"),
    ntc_sample = ntc,
    stringsAsFactors = FALSE
  )
}

empty_common_de <- function() {
  data.frame(
    gene = character(), log2FC = numeric(), log2FC_sh1 = numeric(), log2FC_sh5 = numeric(),
    AveExpr = numeric(), pvalue = numeric(), padj = numeric()
  )
}

build_common_down <- function(a, b) {
  if (is.null(a) || is.null(b)) return(NULL)
  a <- a[!is.na(a$log2FC), ]
  b <- b[!is.na(b$log2FC), ]
  common <- intersect(a$gene[a$log2FC < 0], b$gene[b$log2FC < 0])
  if (length(common) == 0) return(empty_common_de())
  aa <- a[match(common, a$gene), ]
  bb <- b[match(common, b$gene), ]
  data.frame(
    gene = common,
    log2FC = (aa$log2FC + bb$log2FC) / 2,
    log2FC_sh1 = aa$log2FC,
    log2FC_sh5 = bb$log2FC,
    AveExpr = (aa$AveExpr + bb$AveExpr) / 2,
    pvalue = NA_real_,
    padj = NA_real_,
    stringsAsFactors = FALSE
  )
}

# 保留旧名，避免外部脚本调用时报错
build_common_up <- function(a, b) build_common_down(a, b)

full_rank_two <- function(a, b) {
  if (is.null(a) || is.null(b)) return(NULL)
  both <- merge(
    a[, c("gene", "log2FC", "AveExpr", "pvalue", "padj")],
    b[, c("gene", "log2FC", "AveExpr", "pvalue", "padj")],
    by = "gene", suffixes = c("_sh1", "_sh5")
  )
  both$log2FC <- (both$log2FC_sh1 + both$log2FC_sh5) / 2
  both$AveExpr <- (both$AveExpr_sh1 + both$AveExpr_sh5) / 2
  both$pvalue <- NA_real_
  both$padj <- NA_real_
  both
}

passes_padj <- function(padj, have_pvalue) {
  if (!have_pvalue) return(rep(TRUE, length(padj)))
  !is.na(padj) & padj < padj_cutoff
}

order_down <- function(de) {
  if (is.null(de) || nrow(de) == 0) return(de)
  de[order(de$log2FC, decreasing = FALSE, na.last = TRUE), , drop = FALSE]
}

select_by_fc <- function(de, fc, have_pvalue) {
  keep <- !is.na(de$log2FC) & (2^de$log2FC <= fc) & passes_padj(de$padj, have_pvalue)
  if ("log2FC_sh1" %in% names(de)) {
    keep <- keep & (2^de$log2FC_sh1 <= fc) & (2^de$log2FC_sh5 <= fc)
  }
  order_down(de[keep, , drop = FALSE])
}

select_by_topn <- function(de, n, have_pvalue) {
  x <- de[!is.na(de$log2FC) & de$log2FC < 0, , drop = FALSE]
  if ("log2FC_sh1" %in% names(x)) {
    x <- x[x$log2FC_sh1 < 0 & x$log2FC_sh5 < 0, , drop = FALSE]
  }
  sig <- x[passes_padj(x$padj, have_pvalue), , drop = FALSE]
  if (nrow(sig) == 0) sig <- x
  utils::head(order_down(sig), n)
}

# -----------------------------------------------------------------------------
# 6. 基因 ID 转换
# -----------------------------------------------------------------------------
map_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (length(symbols) == 0) return(data.frame(gene = character(), entrez = character()))
  m <- tryCatch(
    clusterProfiler::bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
  )
  if (nrow(m) == 0) {
    m <- tryCatch(
      clusterProfiler::bitr(symbols, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
      error = function(e) data.frame(ENSEMBL = character(), ENTREZID = character())
    )
    if (nrow(m) > 0) names(m)[1] <- "SYMBOL"
  }
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

# -----------------------------------------------------------------------------
# 7. 绘图
# -----------------------------------------------------------------------------
save_gg <- function(plot, path_stub, width = 8, height = 6) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".pdf"), plot, width = width, height = height),
           error = function(e) log_msg("pdf ggsave failed: ", e$message))
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".png"), plot, width = width, height = height, dpi = 300),
           error = function(e) log_msg("png ggsave failed: ", e$message))
}

plot_volcano <- function(de, highlight, title, outfile, fc_line = 1) {
  df <- de
  has_p <- "pvalue" %in% names(df) && any(!is.na(df$pvalue))
  if (has_p) {
    df$y <- -log10(pmax(df$pvalue, 1e-300))
    ylab <- "-log10(p value)"
    hline <- -log10(0.05)
  } else {
    df$y <- df$AveExpr
    ylab <- "Average expression (1-vs-1, no p-value)"
    hline <- NULL
  }
  df$set <- ifelse(df$gene %in% highlight, "selected", "other")
  df$label <- ifelse(df$gene %in% utils::head(highlight, 15), df$gene, NA)
  lfc_abs <- abs(log2(fc_line))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = y, color = set)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.4) +
    ggplot2::scale_color_manual(values = c(other = "grey70", selected = down_fill)) +
    ggplot2::geom_vline(xintercept = c(-lfc_abs, lfc_abs), linetype = 2, color = "grey40") +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, x = "log2 Fold Change", y = ylab, color = NULL)
  if (!is.null(hline)) p <- p + ggplot2::geom_hline(yintercept = hline, linetype = 2, color = "grey40")
  save_gg(p, outfile)
}

plot_scatter_common <- function(de, highlight, title, outfile) {
  if (!all(c("log2FC_sh1", "log2FC_sh5") %in% names(de))) return(invisible(NULL))
  df <- de
  df$set <- ifelse(df$gene %in% highlight, "common_down", "other")
  df$label <- ifelse(df$gene %in% utils::head(highlight, 15), df$gene, NA)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC_sh1, y = log2FC_sh5, color = set)) +
    ggplot2::geom_point(alpha = 0.75, size = 1.6) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
    ggplot2::scale_color_manual(values = c(other = "grey70", common_down = down_fill)) +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, x = "log2FC TG_sh1 vs this NTC", y = "log2FC TG_sh5 vs this NTC", color = NULL)
  save_gg(p, outfile)
}

plot_heatmap <- function(heat_mat, sample_info, genes, title, outfile) {
  genes <- intersect(genes, rownames(heat_mat))
  if (length(genes) > 200) {
    log_msg("Heatmap truncated to 200 genes from ", length(genes), ": ", title)
    genes <- genes[seq_len(200)]
  }
  if (length(genes) < 2) {
    log_msg("Heatmap skipped (<2 genes): ", title)
    return(invisible(NULL))
  }
  sub <- heat_mat[genes, , drop = FALSE]
  ann <- data.frame(Group = sample_info$group, row.names = sample_info$sample)
  ann <- ann[colnames(sub), , drop = FALSE]
  pal <- c(NTC = "#4C78A8", TG_sh1 = "#F58518", TG_sh5 = "#54A24B")
  draw_hm <- function() {
    args <- list(
      mat = sub, scale = "row", annotation_col = ann,
      annotation_colors = list(Group = pal[names(pal) %in% unique(ann$Group)]),
      show_rownames = nrow(sub) <= 80, fontsize_row = 6, main = title,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100),
      clustering_distance_cols = "euclidean"
    )
    tryCatch(
      do.call(pheatmap::pheatmap, c(args, list(clustering_distance_rows = "correlation"))),
      error = function(e) do.call(pheatmap::pheatmap, c(args, list(clustering_distance_rows = "euclidean")))
    )
  }
  grDevices::pdf(paste0(outfile, ".pdf"), width = 8, height = max(6, min(18, 0.18 * nrow(sub) + 3)))
  on.exit({
    while (grDevices::dev.cur() > 1) grDevices::dev.off()
  }, add = TRUE)
  draw_hm()
  grDevices::dev.off()
  grDevices::png(paste0(outfile, ".png"), width = 2400, height = max(1800, 40 * nrow(sub) + 400), res = 300)
  draw_hm()
  grDevices::dev.off()
}

plot_pca <- function(heat_mat, sample_info, outfile) {
  if (ncol(heat_mat) < 2) {
    log_msg("PCA skipped: fewer than 2 samples")
    return(invisible(NULL))
  }
  pca <- tryCatch(stats::prcomp(t(heat_mat), scale. = TRUE), error = function(e) {
    log_msg("PCA failed: ", e$message)
    NULL
  })
  if (is.null(pca)) return(invisible(NULL))
  df <- data.frame(pca$x[, 1:2], group = sample_info$group[match(rownames(pca$x), sample_info$sample)],
                   sample = rownames(pca$x))
  varp <- summary(pca)$importance[2, 1:2] * 100
  p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, color = group, label = sample)) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = "PCA after normalization",
      x = sprintf("PC1 (%.1f%%)", varp[1]),
      y = sprintf("PC2 (%.1f%%)", varp[2])
    )
  save_gg(p, outfile)
}

plot_de_bar <- function(sub, title, outfile) {
  if (nrow(sub) == 0) return(invisible(NULL))
  df <- order_down(sub)
  if (nrow(df) > 60) df <- utils::head(df, 60)
  df$gene <- factor(df$gene, levels = rev(unique(df$gene)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = gene, y = log2FC)) +
    ggplot2::geom_col(fill = down_fill, width = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title, x = NULL, y = "log2 Fold Change")
  save_gg(p, outfile, width = 8, height = max(5, min(16, 0.22 * nrow(df) + 2)))
}

# -----------------------------------------------------------------------------
# 8. 富集分析与绘图（每个子集都必须出图）
# -----------------------------------------------------------------------------
try_save_plot <- function(fun, stub, width = 9, height = 7) {
  p <- tryCatch(fun(), error = function(e) {
    log_msg("Plot failed (", basename(stub), "): ", e$message)
    NULL
  })
  if (is.null(p)) return(invisible(FALSE))
  ok <- tryCatch({
    save_gg(p, stub, width = width, height = height)
    TRUE
  }, error = function(e) {
    log_msg("ggsave failed (", basename(stub), "): ", e$message)
    FALSE
  })
  ok
}

note_empty <- function(stub, msg) {
  writeLines(msg, paste0(stub, "_EMPTY.txt"))
}

plot_ora_object <- function(x, stub, title, fold_change = NULL, also_export_focus = TRUE) {
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
  x2 <- tryCatch(enrichplot::pairwise_termsim(x), error = function(e) NULL)
  if (!is.null(x2) && nrow(df) >= 2) {
    try_save_plot(function() enrichplot::emapplot(x2, showCategory = min(20, nrow(df))) + ggplot2::ggtitle(title),
                  paste0(stub, "_emapplot"), 10, 8)
    try_save_plot(function() enrichplot::treeplot(x2, showCategory = min(20, nrow(df))) + ggplot2::ggtitle(title),
                  paste0(stub, "_treeplot"), 11, 8)
  }
  try_save_plot(function() enrichplot::cnetplot(
    x, showCategory = min(8, nshow), foldChange = fold_change, circular = FALSE
  ) + ggplot2::ggtitle(title), paste0(stub, "_cnetplot"), 10, 8)
  try_save_plot(function() enrichplot::heatplot(x, showCategory = nshow, foldChange = fold_change) +
                  ggplot2::ggtitle(title), paste0(stub, "_heatplot"), 11, 6)
  if (isTRUE(also_export_focus)) export_focus_terms(x, stub, title)
}

plot_gsea_object <- function(x, stub, title, also_export_focus = TRUE) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) {
    note_empty(stub, "no GSEA terms")
    return(invisible(NULL))
  }
  df <- as.data.frame(x)
  utils::write.csv(df, paste0(stub, ".csv"), row.names = FALSE)
  nshow <- min(15, nrow(df))
  try_save_plot(function() {
    p <- enrichplot::dotplot(x, showCategory = nshow, split = ".sign")
    p <- tryCatch(p + ggplot2::facet_grid(. ~ .sign) + ggplot2::ggtitle(title), error = function(e) p + ggplot2::ggtitle(title))
    p
  }, paste0(stub, "_dotplot"), 10, 7)
  try_save_plot(function() enrichplot::ridgeplot(x, showCategory = nshow) + ggplot2::ggtitle(title),
                paste0(stub, "_ridgeplot"), 10, 8)
  ncurve <- min(5, nrow(df))
  try_save_plot(function() enrichplot::gseaplot2(x, geneSetID = seq_len(ncurve), pvalue_table = TRUE, title = title),
                paste0(stub, "_gseaplot"), 10, 8)
  for (i in seq_len(min(3, nrow(df)))) {
    desc <- as.character(df$Description[i])
    try_save_plot(function() enrichplot::gseaplot2(x, geneSetID = i, title = desc),
                  paste0(stub, "_gseaplot_top", i), 8, 6)
  }
  x2 <- tryCatch(enrichplot::pairwise_termsim(x), error = function(e) NULL)
  if (!is.null(x2) && nrow(df) >= 2) {
    try_save_plot(function() enrichplot::emapplot(x2, showCategory = min(20, nrow(df))) + ggplot2::ggtitle(title),
                  paste0(stub, "_emapplot"), 10, 8)
    try_save_plot(function() enrichplot::cnetplot(x, showCategory = min(8, nshow)) + ggplot2::ggtitle(title),
                  paste0(stub, "_cnetplot"), 10, 8)
  }
  if (isTRUE(also_export_focus)) export_focus_terms(x, stub, title)
}

plot_gsea_selected_ids <- function(x, ids, stub, title) {
  if (is.null(x) || length(ids) == 0) {
    note_empty(stub, "no overlapping GSEA terms")
    return(invisible(NULL))
  }
  ids <- ids[ids %in% as.data.frame(x)$ID]
  if (length(ids) == 0) {
    note_empty(stub, "no overlapping GSEA terms")
    return(invisible(NULL))
  }
  ids <- utils::head(ids, 5)
  try_save_plot(function() enrichplot::gseaplot2(x, geneSetID = ids, pvalue_table = TRUE, title = title),
                stub, 10, 8)
}

enrich_or_relax <- function(fun_strict, fun_relax, label) {
  obj <- tryCatch(fun_strict(), error = function(e) {
    log_msg(label, " strict failed: ", e$message)
    NULL
  })
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    attr(obj, "relaxed") <- FALSE
    return(obj)
  }
  obj <- tryCatch(fun_relax(), error = function(e) {
    log_msg(label, " relaxed failed: ", e$message)
    NULL
  })
  if (!is.null(obj)) attr(obj, "relaxed") <- TRUE
  obj
}

title_maybe_relaxed <- function(obj, base) {
  if (isTRUE(attr(obj, "relaxed"))) paste0(base, " (relaxed cutoff)") else base
}

msig_hallmark_map <- function() {
  msig <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  )
  gene_col <- intersect(c("ncbi_gene", "entrez_gene"), names(msig))[1]
  msig[, c("gs_name", gene_col)]
}

# -----------------------------------------------------------------------------
# 8b. 细胞骨架运动 / 线粒体专项富集（不改全库 GO/KEGG 的 p 值）
# -----------------------------------------------------------------------------
.focus_env <- new.env(parent = emptyenv())

focus_keyword_patterns <- function() {
  list(
    cytoskeleton_motility = paste(
      "cytoskelet", "\\bactin\\b", "microtubul", "\\bmyosin\\b", "kinesin", "dynein",
      "lamellipod", "filopod", "stress fiber", "focal adhesion",
      "cell migration", "cell motility", "cell locomotion", "chemotaxis",
      "actin filament", "actin cytoskeleton", "microfilament",
      "rho gtpase", "\\bcdc42\\b", "pseudopod", "podosome", "invadopod",
      "adherens junction", "tight junction", "gap junction",
      "ameboid", "amoeboid", "cell leading edge", "cortical actin",
      "regulation of actin", "myofibril", "sarcomere", "ruffle",
      sep = "|"
    ),
    mitochondria = paste(
      "mitochondr", "oxidative phosphorylat", "respiratory chain",
      "electron transport", "citric acid", "tca cycle", "krebs",
      "mitophag", "oxphos", "respiratory electron", "atp synthase",
      "mitochondrial translation", "cristae", "complex i",
      "inner mitochondrial", "mitochondrial gene", "mitochondrial respir",
      "fatty acid beta-oxidation", "proton-transporting atp",
      "thermogenesis",
      sep = "|"
    )
  )
}

focus_kegg_ids <- function() {
  list(
    cytoskeleton_motility = c(
      "hsa04810", "hsa04510", "hsa04520", "hsa04530", "hsa04540",
      "hsa04512", "hsa04670", "04810", "04510", "04520", "04530", "04540", "04512", "04670"
    ),
    mitochondria = c(
      "hsa00190", "hsa00020", "hsa04137", "hsa04714", "hsa00071", "hsa01212",
      "00190", "00020", "04137", "04714", "00071", "01212"
    )
  )
}

classify_focus_term <- function(id, desc) {
  id <- as.character(id)[1]
  desc <- as.character(desc)[1]
  if (is.na(id)) id <- ""
  if (is.na(desc)) desc <- ""
  kid <- focus_kegg_ids()
  if (id %in% kid$mitochondria || grepl("^MITO_", id)) return("mitochondria")
  if (id %in% kid$cytoskeleton_motility || grepl("^CYTO_", id)) return("cytoskeleton_motility")
  txt <- tolower(paste(id, desc))
  pats <- focus_keyword_patterns()
  is_mito <- grepl(pats$mitochondria, txt, perl = TRUE, ignore.case = TRUE)
  is_cyto <- grepl(pats$cytoskeleton_motility, txt, perl = TRUE, ignore.case = TRUE)
  if (is_cyto && is_mito) return("both")
  if (is_mito) return("mitochondria")
  if (is_cyto) return("cytoskeleton_motility")
  NA_character_
}

plot_focus_term_bar <- function(df, stub, title) {
  if (nrow(df) == 0) return(invisible(NULL))
  lab <- if ("Description" %in% names(df)) as.character(df$Description) else as.character(df$ID)
  df$lab <- paste0(df$focus_class, " | ", lab)
  if ("NES" %in% names(df) && any(is.finite(df$NES))) {
    df <- df[order(abs(df$NES), decreasing = TRUE), , drop = FALSE]
    df <- utils::head(df, 20)
    df$lab <- factor(df$lab, levels = rev(unique(df$lab)))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = NES, y = lab, fill = focus_class)) +
      ggplot2::geom_col() +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::labs(title = title, y = NULL, fill = NULL)
  } else {
    yv <- if ("p.adjust" %in% names(df)) df$p.adjust else df$pvalue
    df$neglog <- -log10(pmax(as.numeric(yv), 1e-300))
    df <- df[order(yv), , drop = FALSE]
    df <- utils::head(df, 20)
    df$lab <- factor(df$lab, levels = rev(unique(df$lab)))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = neglog, y = lab, fill = focus_class)) +
      ggplot2::geom_col() +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::labs(title = title, x = "-log10(p.adjust)", y = NULL, fill = NULL)
  }
  save_gg(p, stub, width = 12, height = max(5, min(12, 0.38 * nrow(df) + 2)))
}

export_focus_terms <- function(x, stub, title) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) return(invisible(NULL))
  df <- as.data.frame(x)
  desc <- if ("Description" %in% names(df)) df$Description else df$ID
  df$genome_wide_rank <- seq_len(nrow(df))
  df$focus_class <- vapply(seq_len(nrow(df)), function(i) {
    classify_focus_term(df$ID[i], desc[i])
  }, character(1))
  hit <- df[!is.na(df$focus_class), , drop = FALSE]
  if (nrow(hit) == 0) {
    note_empty(paste0(stub, "_FOCUS_cytoskeleton_mito"), "no cytoskeleton/mitochondria terms in this result")
    return(invisible(NULL))
  }
  utils::write.csv(hit, paste0(stub, "_FOCUS_cytoskeleton_mito.csv"), row.names = FALSE)
  plot_focus_term_bar(
    hit, paste0(stub, "_FOCUS_cytoskeleton_mito_barplot"),
    paste(title, "| cytoskeleton / mitochondria terms (original p-values)")
  )
}

entrez_for_go <- function(go_id) {
  ids <- tryCatch(
    AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = go_id, column = "ENTREZID",
      keytype = "GOALL", multiVals = "list"
    )[[go_id]],
    error = function(e) {
      tryCatch(
        AnnotationDbi::mapIds(
          org.Hs.eg.db, keys = go_id, column = "ENTREZID",
          keytype = "GO", multiVals = "list"
        )[[go_id]],
        error = function(e2) character()
      )
    }
  )
  unique(as.character(ids[!is.na(ids)]))
}

entrez_for_kegg_path <- function(path_id) {
  pid <- sub("^hsa", "", as.character(path_id))
  hits <- tryCatch(
    AnnotationDbi::select(org.Hs.eg.db, keys = pid, columns = "ENTREZID", keytype = "PATH"),
    error = function(e) NULL
  )
  if (is.null(hits) || nrow(hits) == 0) return(character())
  unique(as.character(hits$ENTREZID[!is.na(hits$ENTREZID)]))
}

append_term2gene <- function(lst, name, entrez) {
  entrez <- unique(as.character(entrez))
  entrez <- entrez[nzchar(entrez) & !is.na(entrez)]
  if (length(entrez) < 8) return(lst)
  lst[[length(lst) + 1]] <- data.frame(
    gs_name = name, entrez = entrez, stringsAsFactors = FALSE
  )
  lst
}

get_focus_term2gene <- function() {
  if (exists("term2gene", envir = .focus_env, inherits = FALSE)) {
    return(.focus_env$term2gene)
  }
  rows <- list()
  go_sets <- c(
    CYTO_GO_cytoskeleton = "GO:0005856",
    CYTO_GO_actin_cytoskeleton = "GO:0015629",
    CYTO_GO_cytoskeleton_organization = "GO:0007010",
    CYTO_GO_actin_cytoskeleton_organization = "GO:0030036",
    CYTO_GO_actin_filament_based_process = "GO:0030029",
    CYTO_GO_actin_filament_organization = "GO:0007015",
    CYTO_GO_regulation_of_actin_cytoskeleton = "GO:0032956",
    CYTO_GO_microtubule_based_process = "GO:0007017",
    CYTO_GO_microtubule = "GO:0005874",
    CYTO_GO_cell_motility = "GO:0048870",
    CYTO_GO_cell_migration = "GO:0016477",
    CYTO_GO_regulation_of_cell_migration = "GO:0030334",
    CYTO_GO_cell_leading_edge = "GO:0031252",
    CYTO_GO_lamellipodium = "GO:0030027",
    CYTO_GO_focal_adhesion = "GO:0005925",
    CYTO_GO_stress_fiber = "GO:0001725",
    CYTO_GO_actin_binding = "GO:0003779",
    CYTO_GO_cytoskeletal_protein_binding = "GO:0008092",
    CYTO_GO_adherens_junction = "GO:0005912",
    MITO_GO_mitochondrion = "GO:0005739",
    MITO_GO_mitochondrial_envelope = "GO:0005740",
    MITO_GO_mitochondrial_inner_membrane = "GO:0005743",
    MITO_GO_mitochondrial_matrix = "GO:0005759",
    MITO_GO_mitochondrion_organization = "GO:0007005",
    MITO_GO_oxidative_phosphorylation = "GO:0006119",
    MITO_GO_electron_transport_chain = "GO:0022900",
    MITO_GO_mito_ATP_synthesis_ETC = "GO:0042775",
    MITO_GO_TCA_cycle = "GO:0006099",
    MITO_GO_mitochondrial_translation = "GO:0032543",
    MITO_GO_mitophagy = "GO:0000422",
    MITO_GO_mitochondrial_transport = "GO:0006839",
    MITO_GO_mitochondrial_respiratory_chain = "GO:0005746",
    MITO_GO_respiratory_chain_complex_assembly = "GO:0033108",
    MITO_GO_fatty_acid_beta_oxidation = "GO:0006635"
  )
  log_msg("Building focused cytoskeleton / mitochondria gene sets from GO")
  for (nm in names(go_sets)) {
    rows <- append_term2gene(rows, nm, entrez_for_go(go_sets[[nm]]))
  }
  kegg_sets <- c(
    CYTO_KEGG_regulation_of_actin_cytoskeleton = "04810",
    CYTO_KEGG_focal_adhesion = "04510",
    CYTO_KEGG_adherens_junction = "04520",
    CYTO_KEGG_tight_junction = "04530",
    CYTO_KEGG_gap_junction = "04540",
    CYTO_KEGG_ECM_receptor_interaction = "04512",
    CYTO_KEGG_leukocyte_transendothelial_migration = "04670",
    MITO_KEGG_oxidative_phosphorylation = "00190",
    MITO_KEGG_citrate_cycle_TCA = "00020",
    MITO_KEGG_mitophagy = "04137",
    MITO_KEGG_thermogenesis = "04714",
    MITO_KEGG_fatty_acid_degradation = "00071"
  )
  for (nm in names(kegg_sets)) {
    rows <- append_term2gene(rows, nm, entrez_for_kegg_path(kegg_sets[[nm]]))
  }
  if (!any(grepl("^CYTO_KEGG_|^MITO_KEGG_", vapply(rows, function(x) x$gs_name[1], character(1))))) {
    msig_kegg <- tryCatch({
      tryCatch(
        msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG"),
        error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
      )
    }, error = function(e) {
      tryCatch(msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG"), error = function(e2) NULL)
    })
    if (!is.null(msig_kegg) && nrow(msig_kegg) > 0) {
      gcol <- intersect(c("ncbi_gene", "entrez_gene"), names(msig_kegg))[1]
      kegg_name_map <- c(
        KEGG_REGULATION_OF_ACTIN_CYTOSKELETON = "CYTO_KEGG_regulation_of_actin_cytoskeleton",
        KEGG_FOCAL_ADHESION = "CYTO_KEGG_focal_adhesion",
        KEGG_ADHERENS_JUNCTION = "CYTO_KEGG_adherens_junction",
        KEGG_TIGHT_JUNCTION = "CYTO_KEGG_tight_junction",
        KEGG_GAP_JUNCTION = "CYTO_KEGG_gap_junction",
        KEGG_ECM_RECEPTOR_INTERACTION = "CYTO_KEGG_ECM_receptor_interaction",
        KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION = "CYTO_KEGG_leukocyte_transendothelial_migration",
        KEGG_OXIDATIVE_PHOSPHORYLATION = "MITO_KEGG_oxidative_phosphorylation",
        KEGG_CITRATE_CYCLE_TCA_CYCLE = "MITO_KEGG_citrate_cycle_TCA",
        KEGG_FATTY_ACID_DEGRADATION = "MITO_KEGG_fatty_acid_degradation"
      )
      for (old in names(kegg_name_map)) {
        hit <- grepl(paste0("^", old, "$"), msig_kegg$gs_name)
        if (!any(hit)) hit <- grepl(old, msig_kegg$gs_name, ignore.case = TRUE)
        rows <- append_term2gene(rows, kegg_name_map[[old]], msig_kegg[[gcol]][hit])
      }
    }
  }
  hm <- tryCatch(msig_hallmark_map(), error = function(e) NULL)
  if (!is.null(hm) && nrow(hm) > 0) {
    hall_map <- c(
      HALLMARK_OXIDATIVE_PHOSPHORYLATION = "MITO_HALLMARK_OXIDATIVE_PHOSPHORYLATION",
      HALLMARK_FATTY_ACID_METABOLISM = "MITO_HALLMARK_FATTY_ACID_METABOLISM",
      HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = "MITO_HALLMARK_REACTIVE_OXYGEN_SPECIES",
      HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "CYTO_HALLMARK_EMT",
      HALLMARK_APICAL_JUNCTION = "CYTO_HALLMARK_APICAL_JUNCTION",
      HALLMARK_MYOGENESIS = "CYTO_HALLMARK_MYOGENESIS"
    )
    gcol <- names(hm)[2]
    for (old in names(hall_map)) {
      rows <- append_term2gene(rows, hall_map[[old]], hm[[gcol]][hm[[1]] == old])
    }
  }
  if (length(rows) == 0) {
    log_msg("WARNING: focused cytoskeleton/mito gene sets are empty")
    .focus_env$term2gene <- data.frame(gs_name = character(), entrez = character())
    return(.focus_env$term2gene)
  }
  t2g <- do.call(rbind, rows)
  t2g <- unique(t2g)
  log_msg(
    "Focused gene sets: ", length(unique(t2g$gs_name)),
    " terms, ", length(unique(t2g$entrez)), " unique Entrez genes"
  )
  .focus_env$term2gene <- t2g
  t2g
}

plot_focus_gene_heatmap <- function(de, t2g, prefix, heat_mat, sample_info, outfile, title) {
  sets <- unique(t2g$gs_name[startsWith(t2g$gs_name, prefix)])
  entrez <- unique(t2g$entrez[t2g$gs_name %in% sets])
  mp <- map_to_entrez(de$gene)
  genes <- unique(mp$gene[mp$entrez %in% entrez])
  genes <- intersect(genes, rownames(heat_mat))
  if (length(genes) < 2) {
    note_empty(outfile, "fewer than 2 mapped genes in this focused set")
    return(invisible(NULL))
  }
  lfc <- de$log2FC[match(genes, de$gene)]
  genes <- genes[order(abs(lfc), decreasing = TRUE, na.last = TRUE)]
  utils::write.csv(
    data.frame(gene = genes, log2FC = lfc[order(abs(lfc), decreasing = TRUE, na.last = TRUE)], stringsAsFactors = FALSE),
    paste0(outfile, "_genes.csv"),
    row.names = FALSE
  )
  plot_heatmap(heat_mat, sample_info, utils::head(genes, 80), title, outfile)
}

run_focused_ora <- function(entrez, de_sub, outdir, label, tag, fc_sym) {
  t2g <- get_focus_term2gene()
  if (is.null(t2g) || nrow(t2g) < 8) return(invisible(NULL))
  fdir <- file.path(outdir, "Focused_cytoskeleton_mito")
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  pref <- paste0(tag, "_")
  obj <- enrich_or_relax(
    function() clusterProfiler::enricher(
      entrez, TERM2GENE = t2g, minGSSize = 5, maxGSSize = 2500,
      pvalueCutoff = 0.05, qvalueCutoff = 0.2
    ),
    function() clusterProfiler::enricher(
      entrez, TERM2GENE = t2g, minGSSize = 5, maxGSSize = 2500,
      pvalueCutoff = 1, qvalueCutoff = 1
    ),
    "focused ORA cytoskeleton/mito"
  )
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    obj <- tryCatch(
      clusterProfiler::setReadable(obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
      error = function(e) obj
    )
  }
  plot_ora_object(
    obj, file.path(fdir, paste0(pref, "ORA_focused_cytoskeleton_mito")),
    title_maybe_relaxed(obj, paste(label, "| ORA focused cytoskeleton / mitochondria")),
    fold_change = fc_sym, also_export_focus = FALSE
  )
}

run_focused_gsea <- function(stats, de, heat_mat, sample_info, outdir, label) {
  t2g <- get_focus_term2gene()
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c("这不是改全库 GO/KEGG 的 p 值或排名。",
      "全基因组 GO/KEGG/GSEA 仍按原统计量排序。",
      "本文件夹只检验细胞骨架运动、细胞迁移、线粒体相关基因集，",
      "所以这些条目会排在这里的结果前面。",
      "全库结果旁边的 *_FOCUS_cytoskeleton_mito.csv 保留原始 p 值和 genome_wide_rank。",
      "若专项 GSEA 仍不显著，说明这些通路在本数据里没有协同变化，不能人为抬到全库第一。"),
    file.path(outdir, "00_README.txt")
  )
  if (is.null(t2g) || nrow(t2g) < 8 || length(stats) < 10) {
    note_empty(file.path(outdir, "GSEA_focused_cytoskeleton_mito"), "too few genes or empty gene sets")
    return(invisible(NULL))
  }
  obj <- enrich_or_relax(
    function() clusterProfiler::GSEA(
      geneList = stats, TERM2GENE = t2g, minGSSize = 8, maxGSSize = 2500,
      pvalueCutoff = 0.05, eps = 0, verbose = FALSE
    ),
    function() clusterProfiler::GSEA(
      geneList = stats, TERM2GENE = t2g, minGSSize = 5, maxGSSize = 2500,
      pvalueCutoff = 1, eps = 0, verbose = FALSE
    ),
    "focused GSEA cytoskeleton/mito"
  )
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    obj <- tryCatch(
      clusterProfiler::setReadable(obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
      error = function(e) obj
    )
  }
  plot_gsea_object(
    obj, file.path(outdir, "GSEA_focused_cytoskeleton_mito"),
    paste(label, "| GSEA focused cytoskeleton / mitochondria"),
    also_export_focus = FALSE
  )
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    df <- as.data.frame(obj)
    df$focus_class <- vapply(seq_len(nrow(df)), function(i) {
      classify_focus_term(df$ID[i], df$Description[i])
    }, character(1))
    cyto <- df[df$focus_class %in% c("cytoskeleton_motility", "both"), , drop = FALSE]
    mito <- df[df$focus_class %in% c("mitochondria", "both"), , drop = FALSE]
    if (nrow(cyto) > 0) {
      plot_focus_term_bar(cyto, file.path(outdir, "GSEA_cytoskeleton_motility_barplot"),
                          paste(label, "| cytoskeleton / motility (focused GSEA)"))
    }
    if (nrow(mito) > 0) {
      plot_focus_term_bar(mito, file.path(outdir, "GSEA_mitochondria_barplot"),
                          paste(label, "| mitochondria (focused GSEA)"))
    }
  }
  plot_focus_gene_heatmap(
    de, t2g, "CYTO_", heat_mat, sample_info,
    file.path(outdir, "heatmap_cytoskeleton_motility_genes"),
    paste(label, "| cytoskeleton / motility genes")
  )
  plot_focus_gene_heatmap(
    de, t2g, "MITO_", heat_mat, sample_info,
    file.path(outdir, "heatmap_mitochondria_genes"),
    paste(label, "| mitochondria genes")
  )
}

build_gsea_cache <- function(de) {
  stats <- ranked_entrez(de)
  out <- list(stats = stats)
  if (length(stats) < 10) return(out)
  gsea_one <- function(fun, label) {
    enrich_or_relax(
      function() fun(pvalueCutoff = 0.05, minGSSize = 10),
      function() fun(pvalueCutoff = 1, minGSSize = 5),
      label
    )
  }
  out$GO_BP <- gsea_one(function(pvalueCutoff, minGSSize) {
    clusterProfiler::gseGO(
      geneList = stats, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID",
      minGSSize = minGSSize, maxGSSize = 500, pvalueCutoff = pvalueCutoff,
      verbose = FALSE, eps = 0
    )
  }, "gseGO_BP")
  out$GO_MF <- gsea_one(function(pvalueCutoff, minGSSize) {
    clusterProfiler::gseGO(
      geneList = stats, OrgDb = org.Hs.eg.db, ont = "MF", keyType = "ENTREZID",
      minGSSize = minGSSize, maxGSSize = 500, pvalueCutoff = pvalueCutoff,
      verbose = FALSE, eps = 0
    )
  }, "gseGO_MF")
  out$GO_CC <- gsea_one(function(pvalueCutoff, minGSSize) {
    clusterProfiler::gseGO(
      geneList = stats, OrgDb = org.Hs.eg.db, ont = "CC", keyType = "ENTREZID",
      minGSSize = minGSSize, maxGSSize = 500, pvalueCutoff = pvalueCutoff,
      verbose = FALSE, eps = 0
    )
  }, "gseGO_CC")
  out$KEGG <- gsea_one(function(pvalueCutoff, minGSSize) {
    clusterProfiler::gseKEGG(
      geneList = stats, organism = "hsa", minGSSize = minGSSize, maxGSSize = 500,
      pvalueCutoff = pvalueCutoff, verbose = FALSE, eps = 0
    )
  }, "gseKEGG")
  if (has_pkg("ReactomePA")) {
    out$Reactome <- gsea_one(function(pvalueCutoff, minGSSize) {
      ReactomePA::gsePathway(
        stats, organism = "human", minGSSize = minGSSize, maxGSSize = 500,
        pvalueCutoff = pvalueCutoff, verbose = FALSE, eps = 0
      )
    }, "gsePathway")
  } else {
    log_msg("ReactomePA not installed; skip Reactome GSEA")
  }
  term2gene <- tryCatch(msig_hallmark_map(), error = function(e) NULL)
  if (!is.null(term2gene)) {
    out$Hallmark <- gsea_one(function(pvalueCutoff, minGSSize) {
      clusterProfiler::GSEA(
        geneList = stats, TERM2GENE = term2gene, minGSSize = minGSSize,
        maxGSSize = 500, pvalueCutoff = pvalueCutoff, eps = 0, verbose = FALSE
      )
    }, "GSEA_Hallmark")
  }
  for (nm in c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Hallmark")) {
    if (!is.null(out[[nm]]) && nrow(as.data.frame(out[[nm]])) > 0) {
      out[[nm]] <- tryCatch(
        clusterProfiler::setReadable(out[[nm]], OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
        error = function(e) out[[nm]]
      )
    }
  }
  out
}

plot_fgsea_hallmark <- function(stats, outdir, title, prefix = "") {
  stub <- function(x) file.path(outdir, paste0(prefix, x))
  if (length(stats) < 5) {
    note_empty(stub("GSEA_Hallmark_fgsea"), "too few ranked genes")
    return(invisible(NULL))
  }
  term2gene <- tryCatch(msig_hallmark_map(), error = function(e) NULL)
  if (is.null(term2gene)) return(invisible(NULL))
  pathways <- split(as.character(term2gene[[2]]), term2gene[[1]])
  fg <- tryCatch(fgsea::fgsea(pathways = pathways, stats = stats, minSize = 5, maxSize = 500), error = function(e) {
    log_msg("fgsea Hallmark failed: ", e$message)
    NULL
  })
  if (is.null(fg) || nrow(fg) == 0) {
    note_empty(stub("GSEA_Hallmark_fgsea"), "no fgsea terms")
    return(invisible(NULL))
  }
  fg <- as.data.frame(fg)
  fg <- fg[order(fg$pval), ]
  utils::write.csv(fg, stub("GSEA_Hallmark_fgsea.csv"), row.names = FALSE)
  plot_df <- utils::head(fg, 15)
  plot_df$pathway <- factor(plot_df$pathway, levels = rev(plot_df$pathway))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = NES, y = pathway, fill = padj < 0.05)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#D62828", "FALSE" = "grey70")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title, y = NULL, fill = "padj < 0.05")
  save_gg(p, stub("GSEA_Hallmark_fgsea_barplot"), 10, 7)
  top_ids <- utils::head(fg$pathway[is.finite(fg$NES)], 3)
  for (i in seq_along(top_ids)) {
    pid <- top_ids[i]
    pe <- tryCatch(fgsea::plotEnrichment(pathways[[pid]], stats) + ggplot2::labs(title = paste(title, "|", pid)), error = function(e) NULL)
    if (!is.null(pe)) save_gg(pe, stub(paste0("GSEA_Hallmark_enrichment_top", i)), 8, 5)
  }
}

gsea_ids_overlapping_genes <- function(gsea_obj, symbols, entrez) {
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) return(character())
  df <- as.data.frame(gsea_obj)
  if (!"core_enrichment" %in% names(df)) return(character())
  keep <- vapply(df$core_enrichment, function(s) {
    gs <- unlist(strsplit(as.character(s), "/"))
    any(gs %in% symbols) || any(gs %in% entrez)
  }, logical(1))
  df$ID[keep]
}

plot_kegg_pathview <- function(kegg_obj, stats, outdir) {
  if (is.null(kegg_obj) || nrow(as.data.frame(kegg_obj)) == 0) return(invisible(NULL))
  if (!has_pkg("pathview")) {
    log_msg("pathview not installed; skip KEGG pathway maps")
    return(invisible(NULL))
  }
  ids <- utils::head(as.character(as.data.frame(kegg_obj)$ID), 3)
  old <- getwd()
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  setwd(outdir)
  on.exit(setwd(old), add = TRUE)
  for (id in ids) {
    pid <- sub("^hsa", "", id)
    tryCatch(
      pathview::pathview(gene.data = stats, pathway.id = pid, species = "hsa", kegg.native = TRUE),
      error = function(e) log_msg("pathview failed for ", id, ": ", e$message)
    )
  }
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
  fc_entrez <- setNames(de_sub$log2FC[match(mp$gene, de_sub$gene)], mp$entrez)
  if (length(entrez) < 3) {
    log_msg("ORA skipped, mapped genes < 3: ", outdir)
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
                    title_maybe_relaxed(ego, paste(label, "| ORA GO", ont, "(not GSEA)")), fold_change = fc_sym)
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
    ek <- tryCatch(clusterProfiler::setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"), error = function(e) ek)
  }
  plot_ora_object(ek, file.path(kg_dir, paste0(pref, "ORA_KEGG")),
                  title_maybe_relaxed(ek, paste(label, "| ORA KEGG (not GSEA)")), fold_change = fc_sym)
  plot_kegg_pathview(ek, fc_entrez, kg_dir)

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
                    title_maybe_relaxed(er, paste(label, "| ORA Reactome pathway (not GSEA)")), fold_change = fc_sym)
  } else {
    note_empty(file.path(pw_dir, paste0(pref, "ORA_Reactome_pathway")), "ReactomePA not installed")
  }

  hm <- enrich_or_relax(
    function() {
      term2gene <- msig_hallmark_map()
      clusterProfiler::enricher(entrez, TERM2GENE = term2gene, pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    },
    function() {
      term2gene <- msig_hallmark_map()
      clusterProfiler::enricher(entrez, TERM2GENE = term2gene, pvalueCutoff = 1, qvalueCutoff = 1)
    },
    "Hallmark"
  )
  if (!is.null(hm) && nrow(as.data.frame(hm)) > 0) {
    hm <- tryCatch(clusterProfiler::setReadable(hm, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"), error = function(e) hm)
  }
  plot_ora_object(hm, file.path(pw_dir, paste0(pref, "ORA_MSigDB_Hallmark_pathway")),
                  title_maybe_relaxed(hm, paste(label, "| ORA Hallmark pathway (not GSEA)")), fold_change = fc_sym)
  tryCatch(
    run_focused_ora(entrez, de_sub, outdir, label, tag, fc_sym),
    error = function(e) log_msg("focused ORA failed: ", e$message)
  )
  writeLines(
    c("This GO/Pathway/KEGG folder is ORA (over-representation), NOT GSEA.",
      "GSEA files are in the sibling folder named GSEA/ and start with GSEA_."),
    file.path(outdir, paste0(pref, "00_ORA_is_not_GSEA.txt"))
  )
}
run_gsea_plots <- function(sub, gsea_cache, outdir, tag, label) {
  gsea_dir <- file.path(outdir, "GSEA")
  dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
  pref <- paste0(tag, "_")
  sub_stats <- ranked_entrez(sub)
  plot_fgsea_hallmark(sub_stats, gsea_dir, paste(label, "| GSEA Hallmark"), prefix = pref)

  if (length(sub_stats) >= 8) {
    term2gene <- tryCatch(msig_hallmark_map(), error = function(e) NULL)
    if (!is.null(term2gene)) {
      hm <- enrich_or_relax(
        function() clusterProfiler::GSEA(
          geneList = sub_stats, TERM2GENE = term2gene, minGSSize = 5,
          maxGSSize = 500, pvalueCutoff = 0.05, eps = 0, verbose = FALSE
        ),
        function() clusterProfiler::GSEA(
          geneList = sub_stats, TERM2GENE = term2gene, minGSSize = 3,
          maxGSSize = 500, pvalueCutoff = 1, eps = 0, verbose = FALSE
        ),
        paste("subset Hallmark GSEA", tag)
      )
      plot_gsea_object(hm, file.path(gsea_dir, paste0(pref, "GSEA_Hallmark")),
                       paste(label, "| GSEA Hallmark (subset ranked)"))
    }
    kegg <- enrich_or_relax(
      function() clusterProfiler::gseKEGG(
        geneList = sub_stats, organism = "hsa", minGSSize = 5, maxGSSize = 500,
        pvalueCutoff = 0.05, verbose = FALSE, eps = 0
      ),
      function() clusterProfiler::gseKEGG(
        geneList = sub_stats, organism = "hsa", minGSSize = 3, maxGSSize = 500,
        pvalueCutoff = 1, verbose = FALSE, eps = 0
      ),
      paste("subset KEGG GSEA", tag)
    )
    if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
      kegg <- tryCatch(clusterProfiler::setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"), error = function(e) kegg)
    }
    plot_gsea_object(kegg, file.path(gsea_dir, paste0(pref, "GSEA_KEGG")),
                     paste(label, "| GSEA KEGG (subset ranked)"))
  }

  mp <- map_to_entrez(sub$gene)
  for (nm in c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Hallmark")) {
    ids <- gsea_ids_overlapping_genes(gsea_cache[[nm]], sub$gene, mp$entrez)
    plot_gsea_selected_ids(
      gsea_cache[[nm]], ids,
      file.path(gsea_dir, paste0(pref, "GSEA_fullrank_overlap_", nm)),
      paste(label, "| GSEA", nm, "(full-rank overlap)")
    )
  }
}

emit_subset_analysis <- function(comp_name, sub, tag, title, outdir, full_de_for_volcano,
                                 heat_mat, sample_info, gsea_cache, fc_line = 1) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(paste("comparison:", comp_name),
      paste("subset:", tag),
      paste("title:", title),
      paste("n_genes:", nrow(sub))),
    file.path(outdir, paste0("00_", tag, "_THIS_FOLDER.txt"))
  )
  utils::write.csv(sub, file.path(outdir, paste0(tag, "_DE_selected_genes.csv")), row.names = FALSE)
  tryCatch(writexl::write_xlsx(sub, file.path(outdir, paste0(tag, "_DE_selected_genes.xlsx"))),
           error = function(e) log_msg("xlsx write failed: ", e$message))
  log_msg(comp_name, " ", tag, ": n = ", nrow(sub))
  if (nrow(sub) == 0) {
    writeLines("no genes", file.path(outdir, paste0(tag, "_EMPTY.txt")))
    return(invisible(NULL))
  }
  tryCatch(plot_de_bar(sub, paste0(title, " | DE genes"), file.path(outdir, paste0(tag, "_DE_log2FC_barplot"))),
           error = function(e) log_msg("DE barplot failed: ", e$message))
  tryCatch(plot_volcano(full_de_for_volcano, sub$gene, title, file.path(outdir, paste0(tag, "_volcano")), fc_line = fc_line),
           error = function(e) log_msg("volcano failed: ", e$message))
  if ("log2FC_sh1" %in% names(full_de_for_volcano)) {
    tryCatch(plot_scatter_common(full_de_for_volcano, sub$gene, title, file.path(outdir, paste0(tag, "_scatter_sh1_sh5"))),
             error = function(e) log_msg("scatter failed: ", e$message))
  }
  tryCatch(plot_heatmap(heat_mat, sample_info, sub$gene, title, file.path(outdir, paste0(tag, "_heatmap"))),
           error = function(e) {
             while (grDevices::dev.cur() > 1) grDevices::dev.off()
             log_msg("heatmap failed: ", e$message)
           })
  tryCatch(run_ora_plots(sub$gene, sub, outdir, title, tag),
           error = function(e) log_msg("ORA/plots failed: ", e$message))
  tryCatch(run_gsea_plots(sub, gsea_cache, outdir, tag, title),
           error = function(e) log_msg("GSEA/plots failed: ", e$message))
}

# -----------------------------------------------------------------------------
# 9. 对单个比较执行全部子集分析（4 个下调 FC x 7 个 topN，全部出图）
# -----------------------------------------------------------------------------
analyze_one_comparison <- function(comp_name, de, full_de_for_volcano, heat_mat, sample_info, have_pvalue, gsea_de = NULL) {
  if (is.null(gsea_de)) gsea_de <- de
  base <- file.path(result_dir, comp_name)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(de, file.path(base, "DE_full.csv"), row.names = FALSE)
  tryCatch(writexl::write_xlsx(de, file.path(base, "DE_full.xlsx")),
           error = function(e) log_msg("DE_full xlsx failed: ", e$message))

  # 先建好 FC / topN 目录，避免只看到 GSEA 文件夹
  fc_dirs <- file.path(base, "FoldChange", names(fc_cutoffs))
  top_dirs <- file.path(base, "TopRank", paste0("top", top_ns))
  invisible(lapply(c(fc_dirs, top_dirs), dir.create, recursive = TRUE, showWarnings = FALSE))
  writeLines(
    c("请打开下面两个文件夹，不要只看 GSEA：",
      "  1_FoldChange_and_TopRank_are_here",
      "  FoldChange/FC_1  FC_0.8  FC_0.6  FC_0.5  （下调 FC <= 阈值）",
      "  TopRank/top50 ... top300  （下调排名，log2FC 最负）",
      "每个子文件夹里：火山图、热图、ORA_GO、ORA通路、ORA_KEGG、以及 GSEA。",
      "细胞骨架运动 / 线粒体专项结果在 Focused_cytoskeleton_mito/（不是改全库排名）。",
      "00_GSEA_all_genes_NOT_FC_or_topN 只是全基因 GSEA，不是分层图。"),
    file.path(base, "00_READ_ME_先看这里.txt")
  )

  gsea_cache <- list()
  log_msg("Subset plots first (volcano/heatmap/ORA), GSEA-all-genes later: ", comp_name)

  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    sub <- select_by_fc(de, fc, have_pvalue)
    tryCatch(
      emit_subset_analysis(
        comp_name, sub, nm, paste0(comp_name, " | down FC <= ", fc),
        file.path(base, "FoldChange", nm),
        full_de_for_volcano, heat_mat, sample_info, gsea_cache, fc_line = fc
      ),
      error = function(e) log_msg("ERROR subset ", comp_name, " ", nm, ": ", e$message)
    )
  }

  for (n in top_ns) {
    tag <- paste0("top", n)
    sub <- select_by_topn(de, n, have_pvalue)
    tryCatch(
      emit_subset_analysis(
        comp_name, sub, tag, paste0(comp_name, " | downregulated top ", n),
        file.path(base, "TopRank", tag),
        full_de_for_volcano, heat_mat, sample_info, gsea_cache, fc_line = 1
      ),
      error = function(e) log_msg("ERROR subset ", comp_name, " ", tag, ": ", e$message)
    )
  }

  tryCatch({
    log_msg("Building full-list GSEA after subset plots: ", comp_name)
    gsea_cache <- build_gsea_cache(gsea_de)
    full_gsea_dir <- file.path(base, "00_GSEA_all_genes_NOT_FC_or_topN")
    dir.create(full_gsea_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines("全基因 GSEA，不是 FC/topN 分层结果。分层图在 FoldChange/ 和 TopRank/。",
               file.path(full_gsea_dir, "00_README.txt"))
    for (nm in c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Hallmark")) {
      plot_gsea_object(gsea_cache[[nm]], file.path(full_gsea_dir, paste0("allGenes_GSEA_", nm)),
                       paste("GSEA", nm, "|", comp_name, "| ALL genes, NOT FC/topN"))
    }
    plot_fgsea_hallmark(gsea_cache$stats, full_gsea_dir,
                        paste("GSEA Hallmark |", comp_name, "| ALL genes"), prefix = "allGenes_")
  }, error = function(e) log_msg("full-list GSEA failed for ", comp_name, ": ", e$message))

  tryCatch({
    log_msg("Focused cytoskeleton / mitochondria GSEA: ", comp_name)
    focus_stats <- if (exists("gsea_cache", inherits = FALSE) && !is.null(gsea_cache$stats)) {
      gsea_cache$stats
    } else {
      ranked_entrez(gsea_de)
    }
    run_focused_gsea(
      focus_stats, gsea_de, heat_mat, sample_info,
      file.path(base, "Focused_cytoskeleton_mito"),
      paste(comp_name, "| all genes")
    )
  }, error = function(e) log_msg("focused GSEA failed for ", comp_name, ": ", e$message))
}

plot_venn_down <- function(de_a, de_b, outdir, label_a, label_b, title_prefix) {
  if (is.null(de_a) || is.null(de_b)) return(invisible(NULL))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  have_pvalue <- FALSE
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    lst <- list(
      x = select_by_fc(de_a, fc, have_pvalue)$gene,
      y = select_by_fc(de_b, fc, have_pvalue)$gene
    )
    names(lst) <- c(label_a, label_b)
    p <- tryCatch({
      if (!has_pkg("ggvenn")) stop("ggvenn not installed")
      ggvenn::ggvenn(lst, fill_color = c("#F58518", "#54A24B")) +
        ggplot2::labs(title = paste0(title_prefix, " | FC <= ", fc))
    }, error = function(e) {
      log_msg("venn failed: ", e$message)
      NULL
    })
    vdir <- file.path(outdir, "FoldChange", nm)
    dir.create(vdir, recursive = TRUE, showWarnings = FALSE)
    if (!is.null(p)) save_gg(p, file.path(vdir, paste0("venn_", nm)), width = 7, height = 6)
  }
  for (n in top_ns) {
    tag <- paste0("top", n)
    lst <- list(
      x = select_by_topn(de_a, n, have_pvalue)$gene,
      y = select_by_topn(de_b, n, have_pvalue)$gene
    )
    names(lst) <- c(label_a, label_b)
    p <- tryCatch({
      if (!has_pkg("ggvenn")) stop("ggvenn not installed")
      ggvenn::ggvenn(lst, fill_color = c("#F58518", "#54A24B")) +
        ggplot2::labs(title = paste0(title_prefix, " | ", tag))
    }, error = function(e) {
      log_msg("venn failed: ", e$message)
      NULL
    })
    vdir <- file.path(outdir, "TopRank", tag)
    dir.create(vdir, recursive = TRUE, showWarnings = FALSE)
    if (!is.null(p)) save_gg(p, file.path(vdir, paste0("venn_", tag)), width = 7, height = 6)
  }
}

plot_venn_up <- function(...) plot_venn_down(...)

# -----------------------------------------------------------------------------
# 10. 主流程
# -----------------------------------------------------------------------------
log_msg("Project dir: ", project_dir)
expr <- load_expression(project_dir)
expr$sample_info <- add_ntc_ids(expr$sample_info)
log_msg("Loaded from ", expr$source, " | genes=", nrow(expr$mat), " samples=", ncol(expr$mat))
print(expr$sample_info)
utils::write.csv(expr$sample_info, file.path(log_dir, "sample_info.csv"), row.names = FALSE)

present <- unique(expr$sample_info$group)
log_msg("Detected groups: ", paste(present, collapse = ", "))
ntc_ids <- unique(stats::na.omit(expr$sample_info$ntc_id))
log_msg("NTC samples kept separate: ", paste(ntc_ids, collapse = ", "))
if (!"NTC" %in% present) stop("未检测到 NTC 对照样本")
if (!"TG_sh1" %in% present) log_msg("WARNING: 未检测到 TG_sh1")
if (!"TG_sh5" %in% present) log_msg("WARNING: 未检测到 TG_sh5")

value_type <- detect_value_type(expr$mat)
log_msg("Value type inferred as: ", value_type)

filt <- filter_low_expression(expr$mat, expr$sample_info, value_type)
norm <- normalize_expression(filt, expr$sample_info, value_type)
expr$sample_info <- expr$sample_info[match(colnames(norm$log_mat), expr$sample_info$sample), ]
expr$sample_info <- add_ntc_ids(expr$sample_info)
utils::write.csv(
  cbind(gene = rownames(norm$log_mat), as.data.frame(norm$log_mat)),
  file.path(result_dir, "normalized_log_matrix.csv"),
  row.names = FALSE
)
plot_pca(norm$heat_mat, expr$sample_info, file.path(result_dir, "00_QC_PCA"))

log_mat <- norm$log_mat
si <- expr$sample_info
sh1 <- find_sample(si, "TG_sh1")
sh5 <- find_sample(si, "TG_sh5")
ntc0 <- find_sample(si, "NTC", "NTC_rep0")
ntc1 <- find_sample(si, "NTC", "NTC_rep1")

de_list <- list()
# 设计1：四个 1-vs-1，各自单独作图，不用 NTC 均值
de_list$TG_sh1_vs_NTC_rep0 <- pairwise_de(log_mat, sh1, ntc0, "TG_sh1_vs_NTC_rep0")
de_list$TG_sh5_vs_NTC_rep0 <- pairwise_de(log_mat, sh5, ntc0, "TG_sh5_vs_NTC_rep0")
de_list$TG_sh1_vs_NTC_rep1 <- pairwise_de(log_mat, sh1, ntc1, "TG_sh1_vs_NTC_rep1")
de_list$TG_sh5_vs_NTC_rep1 <- pairwise_de(log_mat, sh5, ntc1, "TG_sh5_vs_NTC_rep1")
# 设计2：两个 knockdown 等权平均 vs 两个 NTC 的组均值
de_list$TGsh_mean_vs_NTC <- mean_kd_vs_ntc_de(log_mat, si)
# 设计3 / 4：分别相对同一个 NTC 样品的共同下调
de_list$common_down_vs_NTC_rep0 <- build_common_down(de_list$TG_sh1_vs_NTC_rep0, de_list$TG_sh5_vs_NTC_rep0)
de_list$common_down_vs_NTC_rep1 <- build_common_down(de_list$TG_sh1_vs_NTC_rep1, de_list$TG_sh5_vs_NTC_rep1)
# 设计5 / 6：两个 knockdown 等权平均 vs 单个 NTC
de_list$TGsh_mean_vs_NTC_rep0 <- mean_kd_vs_one_ntc_de(log_mat, si, "NTC_rep0", "TGsh_mean_vs_NTC_rep0")
de_list$TGsh_mean_vs_NTC_rep1 <- mean_kd_vs_one_ntc_de(log_mat, si, "NTC_rep1", "TGsh_mean_vs_NTC_rep1")
de_list <- de_list[!vapply(de_list, is.null, logical(1))]

gsea_rank <- list(
  common_down_vs_NTC_rep0 = full_rank_two(de_list$TG_sh1_vs_NTC_rep0, de_list$TG_sh5_vs_NTC_rep0),
  common_down_vs_NTC_rep1 = full_rank_two(de_list$TG_sh1_vs_NTC_rep1, de_list$TG_sh5_vs_NTC_rep1)
)

for (nm in names(de_list)) {
  volcano_df <- de_list[[nm]]
  gsea_de <- de_list[[nm]]
  if (nm %in% names(gsea_rank) && !is.null(gsea_rank[[nm]])) {
    volcano_df <- gsea_rank[[nm]]
    gsea_de <- gsea_rank[[nm]]
  }
  have_p <- any(!is.na(de_list[[nm]]$padj))
  tryCatch(
    analyze_one_comparison(
      nm, de_list[[nm]], volcano_df, norm$heat_mat, si, have_p, gsea_de
    ),
    error = function(e) log_msg("ERROR in comparison ", nm, ": ", e$message)
  )
}

tryCatch(plot_venn_down(
  de_list$TG_sh1_vs_NTC_rep0, de_list$TG_sh5_vs_NTC_rep0,
  file.path(result_dir, "common_down_vs_NTC_rep0"),
  "TG_sh1_vs_NTC_rep0", "TG_sh5_vs_NTC_rep0", "Common down vs NTC_rep0"
), error = function(e) log_msg("venn NTC_rep0 error: ", e$message))
tryCatch(plot_venn_down(
  de_list$TG_sh1_vs_NTC_rep1, de_list$TG_sh5_vs_NTC_rep1,
  file.path(result_dir, "common_down_vs_NTC_rep1"),
  "TG_sh1_vs_NTC_rep1", "TG_sh5_vs_NTC_rep1", "Common down vs NTC_rep1"
), error = function(e) log_msg("venn NTC_rep1 error: ", e$message))

base::writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
log_msg("All done. Results in: ", result_dir)
