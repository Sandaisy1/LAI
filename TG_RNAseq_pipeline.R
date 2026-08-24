#!/usr/bin/env Rscript
# =============================================================================
# TG BRCA 细胞 RNA-seq：按 metastasis_custom_genes.txt 中列出的 GO 通路分析表达
#
# 比较 1–4（本文件）：
#   1) 四个 1-vs-1：TG_sh1 vs NTC_rep0，TG_sh5 vs NTC_rep0，
#                  TG_sh1 vs NTC_rep1，TG_sh5 vs NTC_rep1（各自单独作图）
#   2) mean(TG_sh1, TG_sh5) vs mean(NTC_rep0, NTC_rep1)
#   3) 相对 NTC_rep0 的共同上调（sh1 与 sh5 交集）
#   4) 相对 NTC_rep1 的共同上调（sh1 与 sh5 交集）
# 比较 5–6 在 TG_RNAseq_TGsh_mean_vs_NTC_reps.R，不要为加 5–6 而改本文件主流程。
#
# 分层只两种（均 p<0.05）：上调 FC>=1/1.25/1.5/2，以及上调 top 50–300。
# 气泡图：先全基因组 enrichGO，再抽出列出 GO 的 GeneRatio / p.adjust / Count。
# 不在自选通路上重新校正 p。最终气泡图：p.adjust < 0.2，再按 GeneRatio 取前 15 与前 20。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")

# -----------------------------------------------------------------------------
# 0. 依赖包
# -----------------------------------------------------------------------------
cran_required <- c(
  "writexl", "dplyr", "tidyr", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer"
)
cran_optional <- c("ggvenn")
bioc_required <- c(
  "DESeq2", "edgeR", "limma", "clusterProfiler", "org.Hs.eg.db",
  "AnnotationDbi", "SummarizedExperiment"
)
bioc_optional <- c("GO.db", "GSVA")

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
  candidates <- c(env_dir, "E:/R/TG_BRCA/TG", "E:\\R\\TG_BRCA\\TG", getwd())
  candidates <- unique(candidates[nzchar(candidates)])
  for (d in candidates) {
    if (dir.exists(d) && (
      file.exists(file.path(d, "genes.read_group_tracking")) ||
      file.exists(file.path(d, "genes.fpkm_tracking")) ||
      file.exists(file.path(d, "genes.count_tracking")) ||
      file.exists(file.path(d, "read_groups.info")) ||
      file.exists(file.path(d, "metastasis_custom_genes.txt"))
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

p_cutoffs <- c("p0.05" = 0.05)
fc_cutoffs <- c("FC_1" = 1, "FC_1.25" = 1.25, "FC_1.5" = 1.5, "FC_2" = 2)
top_ns     <- c(50, 75, 100, 150, 200, 250, 300)
bubble_top_ns <- c(15, 20)
padj_plot_cutoff <- 0.2

# -----------------------------------------------------------------------------
# 2. 样本名与基因名
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

add_ntc_ids <- function(sample_info) {
  sample_info$ntc_id <- NA_character_
  ntc <- which(sample_info$group == "NTC")
  if (length(ntc) == 0) return(sample_info)
  labs <- as.character(sample_info$sample[ntc])
  ids <- rep(NA_character_, length(labs))
  ids[grepl("rep1|[_-]1$", labs, ignore.case = TRUE)] <- "NTC_rep1"
  ids[grepl("rep0|[_-]0$", labs, ignore.case = TRUE)] <- "NTC_rep0"
  if (any(is.na(ids)) && length(ntc) == 2) ids[order(labs)] <- c("NTC_rep0", "NTC_rep1")
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

apply_gene_labels <- function(mat, symbols, tracking_ids = NULL, nearest_ref = NULL) {
  raw <- as.character(symbols)
  genes <- clean_gene_names(raw, tracking_ids = tracking_ids, nearest_ref = nearest_ref)
  mat <- collapse_by_gene(mat, genes)
  log_msg(
    "Gene labels: ", nrow(mat), " unique; remaining XLOC: ",
    sum(grepl("^(XLOC|TCONS|CUFF)_", rownames(mat), ignore.case = TRUE))
  )
  mat
}

# -----------------------------------------------------------------------------
# 3. 读入 Cuffdiff
# -----------------------------------------------------------------------------
read_read_group_tracking <- function(path) {
  rg <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  need <- c("tracking_id", "condition", "replicate")
  if (!all(need %in% names(rg))) return(NULL)
  value_col <- if ("raw_frags" %in% names(rg)) {
    "raw_frags"
  } else if ("external_scaled_frags" %in% names(rg)) {
    "external_scaled_frags"
  } else if ("FPKM" %in% names(rg)) {
    "FPKM"
  } else {
    return(NULL)
  }
  log_msg("Cuffdiff value column: ", value_col)
  rg$sample <- paste(rg$condition, rg$replicate, sep = "_rep")
  rg$group <- vapply(as.character(rg$condition), classify_sample, character(1))
  if (all(is.na(rg$group))) rg$group <- vapply(rg$sample, classify_sample, character(1))
  if (all(is.na(rg$group))) return(NULL)
  rg <- rg[!is.na(rg$group), , drop = FALSE]
  if (nrow(rg) == 0) return(NULL)
  wide <- tidyr::pivot_wider(
    rg[, c("tracking_id", "sample", value_col)],
    names_from = "sample", values_from = value_col, values_fn = mean
  )
  gene_map <- NULL
  fpkm_file <- file.path(dirname(path), "genes.fpkm_tracking")
  if (file.exists(fpkm_file)) {
    fp <- utils::read.delim(fpkm_file, check.names = FALSE, stringsAsFactors = FALSE)
    if (all(c("tracking_id", "gene_short_name") %in% names(fp))) gene_map <- fp
  }
  genes <- wide$tracking_id
  nearest <- NULL
  if (!is.null(gene_map)) {
    hit <- match(wide$tracking_id, gene_map$tracking_id)
    genes <- gene_map$gene_short_name[hit]
    if ("nearest_ref_id" %in% names(gene_map)) nearest <- gene_map$nearest_ref_id[hit]
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
  colnames_keep <- NULL
  if (all(is.na(groups))) {
    rg_info <- file.path(dirname(path), "read_groups.info")
    if (file.exists(rg_info)) {
      info <- utils::read.delim(rg_info, check.names = FALSE, stringsAsFactors = FALSE)
      cond <- if ("condition" %in% names(info)) info$condition else info[[2]]
      if (length(cond) == length(val_cols)) {
        groups <- vapply(as.character(cond), classify_sample, character(1))
        colnames_keep <- paste(cond, seq_along(cond), sep = "_rep")
      }
    }
  }
  keep <- !is.na(groups)
  if (sum(keep) < 2) return(NULL)
  mat <- as.matrix(tr[, val_cols[keep], drop = FALSE])
  storage.mode(mat) <- "double"
  mat[is.na(mat)] <- 0
  if (!is.null(colnames_keep)) colnames(mat) <- colnames_keep[keep]
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
    obj <- tryCatch(read_tracking_matrix(ct, "_count$|^q[0-9]+_count$"), error = function(e) NULL)
    if (!is.null(obj)) return(obj)
  }
  fp <- file.path(project_dir, "genes.fpkm_tracking")
  if (file.exists(fp)) {
    log_msg("Reading Cuffdiff FPKM file: genes.fpkm_tracking")
    obj <- tryCatch(read_tracking_matrix(fp, "_FPKM$|^q[0-9]+_FPKM$"), error = function(e) NULL)
    if (!is.null(obj)) return(obj)
  }
  stop("未找到 Cuffdiff 表达文件: ", project_dir)
}

detect_value_type <- function(mat) {
  x <- as.numeric(mat)
  x <- x[is.finite(x)]
  frac_int <- mean(abs(x - round(x)) < 1e-6)
  if (frac_int > 0.85 && stats::quantile(x, 0.95, na.rm = TRUE) > 50) return("counts")
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
    log_mat <- log2(DESeq2::counts(dds, normalized = TRUE) + 1)
    vsd <- tryCatch(
      SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE)),
      error = function(e) log_mat
    )
    list(log_mat = log_mat, heat_mat = vsd, method = "DESeq2")
  } else {
    log_msg("Quantile normalize log2(FPKM+1)")
    log_mat <- limma::normalizeBetweenArrays(log2(pmax(mat, 0) + 1), method = "quantile")
    list(log_mat = log_mat, heat_mat = log_mat, method = "quantile_logFPKM")
  }
}

# -----------------------------------------------------------------------------
# 5. 差异分析（pvalue，不用 padj 顶替）
# -----------------------------------------------------------------------------
empty_de <- function() {
  data.frame(
    gene = character(), log2FC = numeric(), AveExpr = numeric(),
    pvalue = numeric(), padj = numeric(), stringsAsFactors = FALSE
  )
}

pairwise_de <- function(log_mat, treat_sample, ntc_sample, comp_name) {
  if (is.na(treat_sample) || is.na(ntc_sample)) return(NULL)
  if (!all(c(treat_sample, ntc_sample) %in% colnames(log_mat))) return(NULL)
  log_msg(comp_name, " : ", treat_sample, " vs ", ntc_sample, " (1-vs-1, FC only, no p)")
  log2FC <- log_mat[, treat_sample] - log_mat[, ntc_sample]
  data.frame(
    gene = rownames(log_mat),
    log2FC = as.numeric(log2FC),
    AveExpr = as.numeric((log_mat[, treat_sample] + log_mat[, ntc_sample]) / 2),
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
  de <- data.frame(
    gene = rownames(log_mat),
    log2FC = as.numeric(sh_mean - ntc_mean),
    AveExpr = as.numeric((sh_mean + ntc_mean) / 2),
    pvalue = NA_real_,
    padj = NA_real_,
    stringsAsFactors = FALSE
  )
  if (length(ntc) >= 2) {
    kd <- c(sh1, sh5)
    grp <- factor(ifelse(colnames(log_mat) %in% ntc, "NTC",
                         ifelse(colnames(log_mat) %in% kd, "KD", NA_character_)),
                  levels = c("NTC", "KD"))
    keep <- !is.na(grp)
    if (sum(keep) >= 4 && nlevels(droplevels(grp[keep])) == 2) {
      design <- stats::model.matrix(~ grp[keep])
      fit <- limma::eBayes(limma::lmFit(log_mat[, keep, drop = FALSE], design))
      tt <- limma::topTable(fit, coef = ncol(design), number = Inf, sort.by = "none")
      de$pvalue <- tt$P.Value[match(de$gene, rownames(tt))]
      de$padj <- tt$adj.P.Val[match(de$gene, rownames(tt))]
      log_msg("TGsh_mean_vs_NTC : limma p-values from 2-vs-2")
    }
  }
  de
}

mean_kd_vs_one_ntc_de <- function(log_mat, sample_info, ntc_id, comp_name) {
  sh1 <- find_sample(sample_info, "TG_sh1")
  sh5 <- find_sample(sample_info, "TG_sh5")
  ntc <- find_sample(sample_info, "NTC", ntc_id)
  if (is.na(sh1) || is.na(sh5) || is.na(ntc)) return(NULL)
  if (!all(c(sh1, sh5, ntc) %in% colnames(log_mat))) return(NULL)
  log_msg(comp_name, " : mean(", sh1, ", ", sh5, ") vs ", ntc, " (2-vs-1, FC only, no p)")
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

build_common <- function(a, b, direction = c("up", "down")) {
  direction <- match.arg(direction)
  if (is.null(a) || is.null(b)) return(NULL)
  a <- a[!is.na(a$log2FC), ]
  b <- b[!is.na(b$log2FC), ]
  if (direction == "up") {
    ga <- a$gene[a$log2FC > 0]
    gb <- b$gene[b$log2FC > 0]
  } else {
    ga <- a$gene[a$log2FC < 0]
    gb <- b$gene[b$log2FC < 0]
  }
  common <- intersect(ga, gb)
  if (length(common) == 0) {
    return(data.frame(
      gene = character(), log2FC = numeric(), log2FC_sh1 = numeric(), log2FC_sh5 = numeric(),
      AveExpr = numeric(), pvalue = numeric(), padj = numeric(), stringsAsFactors = FALSE
    ))
  }
  aa <- a[match(common, a$gene), ]
  bb <- b[match(common, b$gene), ]
  pv <- if (all(c("pvalue") %in% names(aa)) && any(!is.na(aa$pvalue)) && any(!is.na(bb$pvalue))) {
    pmax(aa$pvalue, bb$pvalue, na.rm = TRUE)
  } else {
    NA_real_
  }
  data.frame(
    gene = common,
    log2FC = (aa$log2FC + bb$log2FC) / 2,
    log2FC_sh1 = aa$log2FC,
    log2FC_sh5 = bb$log2FC,
    AveExpr = (aa$AveExpr + bb$AveExpr) / 2,
    pvalue = pv,
    padj = NA_real_,
    stringsAsFactors = FALSE
  )
}

full_rank_two <- function(a, b) {
  if (is.null(a) || is.null(b)) return(NULL)
  both <- merge(
    a[, c("gene", "log2FC", "AveExpr", "pvalue")],
    b[, c("gene", "log2FC", "AveExpr", "pvalue")],
    by = "gene", suffixes = c("_sh1", "_sh5")
  )
  both$log2FC <- (both$log2FC_sh1 + both$log2FC_sh5) / 2
  both$AveExpr <- (both$AveExpr_sh1 + both$AveExpr_sh5) / 2
  both$pvalue <- NA_real_
  both$padj <- NA_real_
  both
}

has_real_pvalue <- function(de) {
  !is.null(de) && "pvalue" %in% names(de) && any(!is.na(de$pvalue))
}

passes_p <- function(pvalue, p_cutoff, have_pvalue) {
  if (!have_pvalue) return(rep(TRUE, length(pvalue)))
  !is.na(pvalue) & pvalue < p_cutoff
}

select_by_fc <- function(de, fc, p_cutoff, have_pvalue) {
  keep <- !is.na(de$log2FC) & (2^de$log2FC >= fc) & passes_p(de$pvalue, p_cutoff, have_pvalue)
  if ("log2FC_sh1" %in% names(de)) {
    keep <- keep & (2^de$log2FC_sh1 >= fc) & (2^de$log2FC_sh5 >= fc)
  }
  de[keep, , drop = FALSE]
}

select_by_topn <- function(de, n, p_cutoff, have_pvalue) {
  x <- de[!is.na(de$log2FC) & de$log2FC > 0, , drop = FALSE]
  if ("log2FC_sh1" %in% names(x)) {
    x <- x[x$log2FC_sh1 > 0 & x$log2FC_sh5 > 0, , drop = FALSE]
  }
  sig <- x[passes_p(x$pvalue, p_cutoff, have_pvalue), , drop = FALSE]
  if (nrow(sig) == 0) sig <- x
  sig <- sig[order(sig$log2FC, decreasing = TRUE), , drop = FALSE]
  utils::head(sig, n)
}

select_by_direction <- function(de, direction = c("up", "down", "all"), p_cutoff, have_pvalue) {
  direction <- match.arg(direction)
  keep <- !is.na(de$log2FC) & passes_p(de$pvalue, p_cutoff, have_pvalue)
  has_arms <- all(c("log2FC_sh1", "log2FC_sh5") %in% names(de))
  if (direction == "up") {
    keep <- keep & de$log2FC > 0
    if (has_arms) keep <- keep & de$log2FC_sh1 > 0 & de$log2FC_sh5 > 0
  } else if (direction == "down") {
    keep <- keep & de$log2FC < 0
    if (has_arms) keep <- keep & de$log2FC_sh1 < 0 & de$log2FC_sh5 < 0
  } else {
    keep <- keep & de$log2FC != 0
    if (has_arms) {
      keep <- keep & (
        (de$log2FC_sh1 > 0 & de$log2FC_sh5 > 0) |
          (de$log2FC_sh1 < 0 & de$log2FC_sh5 < 0)
      )
    }
  }
  out <- de[keep, , drop = FALSE]
  if (direction == "down") {
    out <- out[order(out$log2FC, decreasing = FALSE), , drop = FALSE]
  } else {
    out <- out[order(out$log2FC, decreasing = TRUE), , drop = FALSE]
  }
  out
}

# -----------------------------------------------------------------------------
# 6. 列出的 GO 通路：读取、映射基因、通路表达分数
# -----------------------------------------------------------------------------
find_custom_go_file <- function(project_dir) {
  candidates <- c(
    file.path(project_dir, "metastasis_custom_genes.txt"),
    file.path(getwd(), "metastasis_custom_genes.txt")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop(
      "找不到 metastasis_custom_genes.txt。请在 ", project_dir,
      " 放置该文件，每行一条通路：GO 名称 + GO 号（例如 cell migration\\tGO:0016477）。",
      "不要只写基因符号。"
    )
  }
  hit[1]
}

go_official_name <- function(go_id) {
  if (!has_pkg("GO.db")) return(NA_character_)
  tryCatch({
    trm <- AnnotationDbi::Term(GO.db::GOTERM[[go_id]])
    if (is.null(trm) || length(trm) == 0) NA_character_ else as.character(trm)
  }, error = function(e) NA_character_)
}

parse_custom_go_file <- function(path) {
  raw <- readLines(path, warn = FALSE, encoding = "UTF-8")
  raw <- sub("\ufeff", "", raw, fixed = TRUE)
  keep <- trimws(raw)
  keep <- keep[nzchar(keep) & !startsWith(keep, "#")]
  if (length(keep) == 0) {
    stop("metastasis_custom_genes.txt 没有有效行（全是空行或注释）。需要 GO 通路名称和 GO 号。")
  }
  rows <- lapply(keep, function(line) {
    ids <- unique(regmatches(line, gregexpr("GO[:_][0-9]{5,7}", line, ignore.case = TRUE))[[1]])
    ids <- toupper(gsub("_", ":", ids, fixed = TRUE))
    ids <- vapply(ids, function(id) {
      num <- sub("^GO:", "", id)
      sprintf("GO:%07d", as.integer(num))
    }, character(1), USE.NAMES = FALSE)
    name <- trimws(gsub("GO[:_][0-9]{5,7}", "", line, ignore.case = TRUE))
    name <- gsub("[,;|/]+$", "", name)
    name <- gsub("^[,;|/]+", "", name)
    name <- gsub("\\s+", " ", name)
    if (length(ids) == 0) return(NULL)
    data.frame(go_id = ids, name_in_file = name, stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(tab) || nrow(tab) == 0) {
    stop(
      "metastasis_custom_genes.txt 里没有解析到 GO 号。",
      "请写成「通路名称 + GO:0007411」，不要只放基因符号。"
    )
  }
  tab <- tab[!duplicated(tab$go_id), , drop = FALSE]
  tab$name <- vapply(seq_len(nrow(tab)), function(i) {
    official <- go_official_name(tab$go_id[i])
    file_nm <- tab$name_in_file[i]
    if (nzchar(file_nm)) file_nm else if (!is.na(official)) official else tab$go_id[i]
  }, character(1))
  tab$safe_id <- gsub(":", "_", tab$go_id, fixed = TRUE)
  log_msg("Custom GO pathways: ", nrow(tab), " from ", path)
  tab
}

map_go_to_symbols <- function(go_ids, expressed) {
  sets <- tryCatch(
    AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db, keys = go_ids, column = "SYMBOL",
      keytype = "GOALL", multiVals = "list"
    ),
    error = function(e) {
      log_msg("GOALL mapping failed: ", e$message)
      setNames(vector("list", length(go_ids)), go_ids)
    }
  )
  out <- lapply(go_ids, function(id) {
    g <- unique(as.character(sets[[id]]))
    g <- g[!is.na(g) & nzchar(g)]
    intersect(g, expressed)
  })
  names(out) <- go_ids
  n <- vapply(out, length, integer(1))
  log_msg("GO gene mapping (expressed genes): ", paste(paste0(go_ids, "=", n), collapse = "; "))
  out
}

all_pathway_genes <- function(go_sets) unique(unlist(go_sets, use.names = FALSE))

pathway_score_matrix <- function(log_mat, go_sets) {
  z <- t(scale(t(log_mat)))
  z[!is.finite(z)] <- 0
  score <- vapply(go_sets, function(genes) {
    genes <- intersect(genes, rownames(z))
    if (length(genes) == 0) return(rep(NA_real_, ncol(z)))
    colMeans(z[genes, , drop = FALSE], na.rm = TRUE)
  }, numeric(ncol(z)))
  if (is.null(dim(score))) {
    score <- matrix(score, nrow = ncol(z), dimnames = list(colnames(z), names(go_sets)))
  } else {
    score <- t(score)
    colnames(score) <- colnames(z)
  }
  if (has_pkg("GSVA") && length(go_sets) > 0) {
    gs <- lapply(go_sets, function(g) intersect(g, rownames(log_mat)))
    gs <- gs[vapply(gs, length, integer(1)) >= 2]
    if (length(gs) > 0) {
      gsva_mat <- tryCatch({
        expr <- as.matrix(log_mat)
        if (utils::packageVersion("GSVA") >= "1.50.0" && exists("ssgseaParam", envir = asNamespace("GSVA"), inherits = FALSE)) {
          GSVA::gsva(GSVA::ssgseaParam(expr, gs), verbose = FALSE)
        } else {
          GSVA::gsva(expr, gs, method = "ssgsea", kcdf = "Gaussian", verbose = FALSE)
        }
      }, error = function(e) {
        log_msg("ssGSEA failed, use mean z-score: ", e$message)
        NULL
      })
      if (!is.null(gsva_mat)) {
        log_msg("Pathway scores: ssGSEA")
        return(list(score = gsva_mat, method = "ssgsea"))
      }
    }
  }
  log_msg("Pathway scores: mean z-score of genes in each GO")
  list(score = score, method = "mean_z")
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

write_table <- function(df, stub) {
  dir.create(dirname(stub), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, paste0(stub, ".csv"), row.names = FALSE)
  tryCatch(writexl::write_xlsx(df, paste0(stub, ".xlsx")),
           error = function(e) log_msg("xlsx write failed: ", e$message))
}

plot_volcano <- function(de, highlight, title, outfile, fc_line = 1, p_line = 0.05) {
  df <- de
  has_p <- has_real_pvalue(df)
  if (has_p) {
    df$y <- -log10(pmax(df$pvalue, 1e-300))
    ylab <- "-log10(p value)"
    hline <- -log10(p_line)
  } else {
    df$y <- df$AveExpr
    ylab <- "Average expression (no p-value)"
    hline <- NULL
  }
  df$set <- ifelse(df$gene %in% highlight, "selected", "other")
  df$label <- ifelse(df$gene %in% utils::head(highlight, 15), df$gene, NA)
  lfc_line <- log2(fc_line)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = y, color = set)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.4) +
    ggplot2::scale_color_manual(values = c(other = "grey70", selected = "#D62828")) +
    ggplot2::geom_vline(xintercept = c(-lfc_line, lfc_line), linetype = 2, color = "grey40") +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, x = "log2 Fold Change", y = ylab, color = NULL)
  if (!is.null(hline)) p <- p + ggplot2::geom_hline(yintercept = hline, linetype = 2, color = "grey40")
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
    writeLines("heatmap skipped: fewer than 2 genes", paste0(outfile, "_EMPTY.txt"))
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
  if (ncol(heat_mat) < 2) return(invisible(NULL))
  pca <- tryCatch(stats::prcomp(t(heat_mat), scale. = TRUE), error = function(e) NULL)
  if (is.null(pca)) return(invisible(NULL))
  df <- data.frame(
    pca$x[, 1:2],
    group = sample_info$group[match(rownames(pca$x), sample_info$sample)],
    sample = rownames(pca$x)
  )
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

map_to_entrez <- function(symbols) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  empty <- data.frame(gene = character(), entrez = character(), stringsAsFactors = FALSE)
  if (length(symbols) == 0) return(empty)
  m <- tryCatch(
    clusterProfiler::bitr(
      symbols, fromType = "SYMBOL", toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    ),
    error = function(e) NULL
  )
  if (is.null(m) || nrow(m) == 0) return(empty)
  m <- m[!duplicated(m$SYMBOL), , drop = FALSE]
  data.frame(gene = m$SYMBOL, entrez = as.character(m$ENTREZID), stringsAsFactors = FALSE)
}

parse_gene_ratio <- function(x) {
  vapply(as.character(x), function(s) {
    p <- strsplit(s, "/", fixed = TRUE)[[1]]
    if (length(p) < 2) return(suppressWarnings(as.numeric(s)))
    as.numeric(p[1]) / as.numeric(p[2])
  }, numeric(1), USE.NAMES = FALSE)
}

ora_df_from_enrich <- function(ego, ont) {
  if (is.null(ego)) return(NULL)
  df <- as.data.frame(ego)
  if (nrow(df) == 0) return(NULL)
  df$ONTOLOGY <- ont
  df$genome_wide_rank <- seq_len(nrow(df))
  df$GeneRatio_num <- parse_gene_ratio(df$GeneRatio)
  df
}

# 与原先正常 GO 分析相同：clusterProfiler::enrichGO，不限制 universe，
# BH 在全库 GO 条目上校正。pvalueCutoff=1 只为留下自选通路，不改 p.adjust。
run_genome_enrichGO <- function(genes, go_dir, tag) {
  dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(
      "全库 GO 气泡图在本文件夹，不要只看 xlsx：",
      "  *_ORA_GO_BP_dotplot_top15.pdf / top20.pdf",
      "  *_ORA_GO_CC_dotplot_top15.pdf / top20.pdf",
      "  *_ORA_GO_MF_dotplot_top15.pdf / top20.pdf"
    ),
    file.path(go_dir, "00_READ_ME_气泡图在这里.txt")
  )
  entrez <- unique(map_to_entrez(genes)$entrez)
  if (length(entrez) < 3) {
    writeLines(
      paste("mapped_entrez", length(entrez)),
      file.path(go_dir, paste0(tag, "_ORA_skipped.txt"))
    )
    return(NULL)
  }
  pieces <- list()
  for (ont in c("BP", "MF", "CC")) {
    ego <- tryCatch(
      clusterProfiler::enrichGO(
        gene = entrez,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 10,
        maxGSSize = 500,
        readable = TRUE
      ),
      error = function(e) {
        log_msg("enrichGO ", ont, " failed: ", e$message)
        NULL
      }
    )
    df <- ora_df_from_enrich(ego, ont)
    if (!is.null(df)) {
      write_table(df, file.path(go_dir, paste0(tag, "_ORA_GO_", ont)))
      for (n_show in bubble_top_ns) {
        top_tag <- paste0("top", n_show)
        tryCatch(
          plot_ora_bubble(
            df,
            paste0(tag, " | ORA GO ", ont, " ", top_tag,
                   " (p.adjust<", padj_plot_cutoff, ", by GeneRatio)"),
            file.path(go_dir, paste0(tag, "_ORA_GO_", ont, "_dotplot_", top_tag)),
            n_show = n_show
          ),
          error = function(e) log_msg("GO ", ont, " bubble failed (", top_tag, "): ", e$message)
        )
      }
      pieces[[ont]] <- df
    }
  }
  if (length(pieces) == 0) return(NULL)
  dplyr::bind_rows(pieces)
}

extract_listed_go_ora <- function(genome_ora, go_tab) {
  missing <- go_tab$go_id
  empty <- list(all_listed = NULL, missing = missing)
  if (is.null(genome_ora) || nrow(genome_ora) == 0) return(empty)
  id_col <- if ("ID" %in% names(genome_ora)) "ID" else "go_id"
  hit <- genome_ora[genome_ora[[id_col]] %in% go_tab$go_id, , drop = FALSE]
  if (nrow(hit) == 0) return(empty)
  hit <- hit[order(hit$p.adjust, hit$pvalue, -hit$Count), , drop = FALSE]
  hit <- hit[!duplicated(hit[[id_col]]), , drop = FALSE]
  hit$go_id <- hit[[id_col]]
  if (!"GeneRatio_num" %in% names(hit)) hit$GeneRatio_num <- parse_gene_ratio(hit$GeneRatio)
  hit$listed_rank <- seq_len(nrow(hit))
  list(
    all_listed = hit,
    missing = setdiff(go_tab$go_id, hit$go_id)
  )
}

plot_ora_bubble <- function(ora, title, outfile, n_show) {
  df <- ora
  if (is.null(df) || nrow(df) == 0) {
    writeLines("no GO terms to plot", paste0(outfile, "_EMPTY.txt"))
    return(invisible(NULL))
  }
  if (!"GeneRatio_num" %in% names(df)) {
    df$GeneRatio_num <- if (is.numeric(df$GeneRatio)) df$GeneRatio else parse_gene_ratio(df$GeneRatio)
  }
  df <- df[is.finite(df$GeneRatio_num), , drop = FALSE]
  df <- df[is.finite(df$p.adjust) & df$p.adjust < padj_plot_cutoff, , drop = FALSE]
  if (nrow(df) == 0) {
    writeLines(
      paste0("no GO terms with p.adjust < ", padj_plot_cutoff),
      paste0(outfile, "_EMPTY.txt")
    )
    return(invisible(NULL))
  }
  df <- df[order(-df$GeneRatio_num, df$p.adjust, df$pvalue), , drop = FALSE]
  df <- utils::head(df, n_show)
  gid <- if ("go_id" %in% names(df)) df$go_id else df$ID
  df$label <- paste0(df$Description, " (", gid, ")")
  df$label <- factor(df$label, levels = rev(unique(as.character(df$label))))
  df$p_adjust <- df$p.adjust
  df$p_adjust[!is.finite(df$p_adjust)] <- 1
  df$p_adjust <- pmax(df$p_adjust, 1e-300)
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = GeneRatio_num, y = label, size = Count, fill = p_adjust)
  ) +
    ggplot2::geom_point(shape = 21, color = "grey30", stroke = 0.5) +
    ggplot2::scale_size_continuous(range = c(6, 18)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.16))) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.55)) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(6, 10, 6, 6),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = title,
      x = "GeneRatio",
      y = NULL,
      fill = "p.adjust",
      size = "Count"
    )
  rng <- range(df$p_adjust[df$p_adjust > 0], na.rm = TRUE)
  diverging <- c("blue", "white", "red")
  if (is.finite(rng[1]) && rng[1] > 0 && rng[2] / rng[1] >= 10) {
    p <- p + ggplot2::scale_fill_gradientn(colours = diverging, trans = "log10")
  } else {
    p <- p + ggplot2::scale_fill_gradientn(colours = diverging)
  }
  h <- max(5.2, min(10.5, 0.32 * nrow(df) + 2.4))
  save_gg(p, outfile, width = 9, height = h)
}

plot_pathway_mean_fc <- function(de, go_tab, go_sets, title, outfile) {
  rows <- lapply(seq_len(nrow(go_tab)), function(i) {
    id <- go_tab$go_id[i]
    genes <- intersect(go_sets[[id]], de$gene)
    fc <- de$log2FC[match(genes, de$gene)]
    data.frame(
      go_id = id,
      Description = go_tab$name[i],
      n_genes = length(genes),
      mean_log2FC = if (length(fc) == 0) NA_real_ else mean(fc, na.rm = TRUE),
      median_log2FC = if (length(fc) == 0) NA_real_ else stats::median(fc, na.rm = TRUE),
      n_up = sum(fc > 0, na.rm = TRUE),
      n_down = sum(fc < 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  write_table(df, paste0(outfile, "_table"))
  df <- df[is.finite(df$mean_log2FC), , drop = FALSE]
  if (nrow(df) == 0) return(invisible(df))
  df$label <- paste0(df$Description, " (", df$go_id, ")")
  df$label <- factor(df$label, levels = df$label[order(df$mean_log2FC)])
  p <- ggplot2::ggplot(df, ggplot2::aes(x = mean_log2FC, y = label, fill = mean_log2FC > 0)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c("#4C78A8", "#D62828"), guide = "none") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, x = "Mean log2FC of genes in GO", y = NULL)
  save_gg(p, outfile, width = 10, height = max(5, min(16, 0.35 * nrow(df) + 3)))
  df
}

plot_score_heatmap <- function(score, sample_info, go_tab, title, outfile) {
  if (is.null(score) || nrow(score) < 1 || ncol(score) < 2) return(invisible(NULL))
  lab <- go_tab$name[match(rownames(score), go_tab$go_id)]
  lab[is.na(lab) | !nzchar(lab)] <- rownames(score)
  rownames(score) <- make.unique(paste0(lab, " (", rownames(score), ")"))
  ann <- data.frame(Group = sample_info$group, row.names = sample_info$sample)
  ann <- ann[colnames(score), , drop = FALSE]
  pal <- c(NTC = "#4C78A8", TG_sh1 = "#F58518", TG_sh5 = "#54A24B")
  draw <- function() {
    pheatmap::pheatmap(
      score, scale = "none", annotation_col = ann,
      annotation_colors = list(Group = pal[names(pal) %in% unique(ann$Group)]),
      main = title, fontsize_row = 8,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100)
    )
  }
  grDevices::pdf(paste0(outfile, ".pdf"), width = 9, height = max(5, min(16, 0.28 * nrow(score) + 3)))
  on.exit({
    while (grDevices::dev.cur() > 1) grDevices::dev.off()
  }, add = TRUE)
  draw()
  grDevices::dev.off()
  grDevices::png(paste0(outfile, ".png"), width = 2400, height = max(1600, 50 * nrow(score) + 400), res = 300)
  draw()
  grDevices::dev.off()
}

# -----------------------------------------------------------------------------
# 8. 每个子集：差异表、火山图、热图、列出 GO 的气泡图
# -----------------------------------------------------------------------------
emit_subset_analysis <- function(comp_name, sub, tag, title, outdir, full_de,
                                 heat_mat, sample_info, go_tab, go_sets,
                                 pathway_genes, fc_line = 1, p_line = 0.05) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(paste("comparison:", comp_name), paste("subset:", tag),
      paste("title:", title), paste("n_genes:", nrow(sub))),
    file.path(outdir, paste0("00_", tag, "_THIS_FOLDER.txt"))
  )
  write_table(sub, file.path(outdir, paste0(tag, "_DE_selected_genes")))
  log_msg(comp_name, " ", tag, ": n = ", nrow(sub))
  if (nrow(sub) == 0) {
    writeLines("no genes", file.path(outdir, paste0(tag, "_EMPTY.txt")))
    return(invisible(NULL))
  }
  tryCatch(
    plot_volcano(full_de, sub$gene, title, file.path(outdir, paste0(tag, "_volcano")),
                 fc_line = fc_line, p_line = p_line),
    error = function(e) log_msg("volcano failed: ", e$message)
  )
  tryCatch(
    plot_heatmap(heat_mat, sample_info, sub$gene, paste(title, "| DE genes"),
                 file.path(outdir, paste0(tag, "_heatmap"))),
    error = function(e) {
      while (grDevices::dev.cur() > 1) grDevices::dev.off()
      log_msg("heatmap failed: ", e$message)
    }
  )
  in_path <- intersect(sub$gene, pathway_genes)
  tryCatch(
    plot_heatmap(
      heat_mat, sample_info, in_path,
      paste(title, "| listed GO genes"),
      file.path(outdir, paste0(tag, "_pathway_gene_heatmap"))
    ),
    error = function(e) log_msg("pathway heatmap failed: ", e$message)
  )
  cg_dir <- file.path(outdir, "CustomGO")
  go_dir <- file.path(outdir, "GO")
  dir.create(cg_dir, recursive = TRUE, showWarnings = FALSE)
  genome_ora <- tryCatch(
    run_genome_enrichGO(sub$gene, go_dir, tag),
    error = function(e) {
      log_msg("genome enrichGO failed: ", e$message)
      NULL
    }
  )
  extracted <- extract_listed_go_ora(genome_ora, go_tab)
  if (!is.null(extracted$all_listed)) {
    write_table(extracted$all_listed, file.path(cg_dir, paste0(tag, "_ORA_CustomGO")))
  }
  if (length(extracted$missing) > 0) {
    writeLines(
      c(
        "下列自选 GO 未出现在全基因组 enrichGO 结果中（可能超出 minGSSize/maxGSSize，或该子集没有注释到）。",
        "气泡图使用全库 GeneRatio / p.adjust / Count，未重新计算这些缺失项。",
        extracted$missing
      ),
      file.path(cg_dir, paste0(tag, "_listed_GO_not_in_genome_ORA.txt"))
    )
  }
  for (n_show in bubble_top_ns) {
    top_tag <- paste0("top", n_show)
    if (!is.null(extracted$all_listed)) {
      ranked <- extracted$all_listed
      ranked <- ranked[is.finite(ranked$p.adjust) & ranked$p.adjust < padj_plot_cutoff, , drop = FALSE]
      ranked <- ranked[order(-ranked$GeneRatio_num, ranked$p.adjust, ranked$pvalue), , drop = FALSE]
      write_table(
        utils::head(ranked, n_show),
        file.path(cg_dir, paste0(tag, "_ORA_CustomGO_", top_tag))
      )
    }
    tryCatch(
      plot_ora_bubble(
        extracted$all_listed,
        paste0(title, " | ORA listed GO ", top_tag,
               " (p.adjust<", padj_plot_cutoff, ", by GeneRatio)"),
        file.path(cg_dir, paste0(tag, "_ORA_CustomGO_dotplot_", top_tag)),
        n_show = n_show
      ),
      error = function(e) log_msg("bubble failed (", top_tag, "): ", e$message)
    )
  }
}

run_two_tracks <- function(comp_name, de, full_de, heat_mat, sample_info,
                            go_tab, go_sets, pathway_genes) {
  base <- file.path(result_dir, comp_name)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  write_table(de, file.path(base, "DE_full"))
  have_p <- has_real_pvalue(de)
  writeLines(
    c(
      "本比较只跑两种分层（均要求 p<0.05，无法估 p 时不伪造 p，仍按 FC/排名）：",
      "  1) FoldChange 上调 FC >= 1 / 1.25 / 1.5 / 2",
      "  2) 上调排名 top 50 / 75 / 100 / 150 / 200 / 250 / 300",
      "每个非空子集：差异表、火山图、热图；GO/ 全库 BP/CC/MF 气泡图（top15 与 top20）。",
      "最终气泡图只保留 p.adjust<0.2，再按 GeneRatio 取前15和前20。",
      paste("p-value estimated:", have_p)
    ),
    file.path(base, "00_READ_ME.txt")
  )
  if (!have_p) {
    writeLines(
      "1-vs-1 或 2-vs-1 无法估计 p，未伪造 p 值。p0.05 目录里是 FC/topN 分层（未用 p 过滤）。",
      file.path(base, "NO_PVALUE.txt")
    )
  }

  p_tag <- "p0.05"
  pcut <- unname(p_cutoffs[[p_tag]])
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    sub <- select_by_fc(de, fc, pcut, have_p)
    if (nrow(sub) > 0) sub <- sub[order(sub$log2FC, decreasing = TRUE), , drop = FALSE]
    emit_subset_analysis(
      comp_name, sub, paste0(p_tag, "_", nm),
      paste0(comp_name, " | ", p_tag, " | up FC >= ", fc),
      file.path(base, p_tag, "FoldChange", nm),
      full_de, heat_mat, sample_info, go_tab, go_sets, pathway_genes,
      fc_line = fc, p_line = pcut
    )
  }
  for (n in top_ns) {
    tag <- paste0("top", n)
    sub <- select_by_topn(de, n, pcut, have_p)
    emit_subset_analysis(
      comp_name, sub, paste0(p_tag, "_", tag),
      paste0(comp_name, " | ", p_tag, " | upregulated top ", n),
      file.path(base, p_tag, "TopRank", tag),
      full_de, heat_mat, sample_info, go_tab, go_sets, pathway_genes,
      fc_line = 1, p_line = pcut
    )
  }

  tryCatch(
    plot_pathway_mean_fc(
      de, go_tab, go_sets,
      paste(comp_name, "| mean log2FC of genes in each listed GO"),
      file.path(base, "CustomGO", paste0(comp_name, "_pathway_mean_log2FC"))
    ),
    error = function(e) log_msg("pathway mean FC failed: ", e$message)
  )
}

write_pathway_expression <- function(log_mat, heat_mat, sample_info, go_tab, go_sets) {
  out <- file.path(result_dir, "00_PathwayExpression")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  set_rows <- do.call(rbind, lapply(seq_len(nrow(go_tab)), function(i) {
    data.frame(
      go_id = go_tab$go_id[i],
      name = go_tab$name[i],
      n_expressed_genes = length(go_sets[[go_tab$go_id[i]]]),
      genes = paste(go_sets[[go_tab$go_id[i]]], collapse = "/"),
      stringsAsFactors = FALSE
    )
  }))
  write_table(set_rows, file.path(out, "listed_GO_gene_sets"))
  scored <- pathway_score_matrix(log_mat, setNames(go_sets, names(go_sets)))
  score <- scored$score
  write.csv(
    cbind(go_id = rownames(score), as.data.frame(score)),
    file.path(out, paste0("pathway_scores_", scored$method, ".csv")),
    row.names = FALSE
  )
  plot_score_heatmap(
    score, sample_info, go_tab,
    paste0("Listed GO pathway expression (", scored$method, ")"),
    file.path(out, "pathway_score_heatmap")
  )
  per <- file.path(out, "per_GO")
  dir.create(per, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(go_tab))) {
    id <- go_tab$go_id[i]
    genes <- go_sets[[id]]
    plot_heatmap(
      heat_mat, sample_info, genes,
      paste0(go_tab$name[i], " (", id, ")"),
      file.path(per, paste0(go_tab$safe_id[i], "_gene_heatmap"))
    )
  }
}

plot_venn_up <- function(de_a, de_b, outdir, label_a, label_b, title_prefix,
                         p_cutoff, p_tag, have_pvalue) {
  if (is.null(de_a) || is.null(de_b) || !has_pkg("ggvenn")) return(invisible(NULL))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    lst <- list(
      x = select_by_fc(de_a, fc, p_cutoff, have_pvalue)$gene,
      y = select_by_fc(de_b, fc, p_cutoff, have_pvalue)$gene
    )
    names(lst) <- c(label_a, label_b)
    p <- tryCatch(
      ggvenn::ggvenn(lst, fill_color = c("#F58518", "#54A24B")) +
        ggplot2::labs(title = paste0(title_prefix, " | ", p_tag, " | FC >= ", fc)),
      error = function(e) NULL
    )
    vdir <- file.path(outdir, p_tag, "FoldChange", nm)
    dir.create(vdir, recursive = TRUE, showWarnings = FALSE)
    if (!is.null(p)) save_gg(p, file.path(vdir, paste0("venn_", p_tag, "_", nm)), width = 7, height = 6)
  }
}

analyze_one_comparison <- function(comp_name, de, full_de, heat_mat, sample_info,
                                   go_tab, go_sets, pathway_genes) {
  run_two_tracks(comp_name, de, full_de, heat_mat, sample_info, go_tab, go_sets, pathway_genes)
}

# -----------------------------------------------------------------------------
# 9. 主流程：比较 1–4
# -----------------------------------------------------------------------------
run_comparisons_1_to_4 <- function() {
  log_msg("Project dir: ", project_dir)
  go_path <- find_custom_go_file(project_dir)
  go_tab <- parse_custom_go_file(go_path)
  expr <- load_expression(project_dir)
  expr$sample_info <- add_ntc_ids(expr$sample_info)
  log_msg("Loaded from ", expr$source, " | genes=", nrow(expr$mat), " samples=", ncol(expr$mat))
  utils::write.csv(expr$sample_info, file.path(log_dir, "sample_info.csv"), row.names = FALSE)

  if (!"NTC" %in% expr$sample_info$group) stop("未检测到 NTC 对照样本")
  value_type <- detect_value_type(expr$mat)
  log_msg("Value type inferred as: ", value_type)
  filt <- filter_low_expression(expr$mat, expr$sample_info, value_type)
  norm <- normalize_expression(filt, expr$sample_info, value_type)
  si <- expr$sample_info[match(colnames(norm$log_mat), expr$sample_info$sample), ]
  si <- add_ntc_ids(si)
  utils::write.csv(
    cbind(gene = rownames(norm$log_mat), as.data.frame(norm$log_mat)),
    file.path(result_dir, "normalized_log_matrix.csv"),
    row.names = FALSE
  )
  plot_pca(norm$heat_mat, si, file.path(result_dir, "00_QC_PCA"))

  go_sets <- map_go_to_symbols(go_tab$go_id, rownames(norm$log_mat))
  names(go_sets) <- go_tab$go_id
  pathway_genes <- all_pathway_genes(go_sets)
  write_pathway_expression(norm$log_mat, norm$heat_mat, si, go_tab, go_sets)

  log_mat <- norm$log_mat
  sh1 <- find_sample(si, "TG_sh1")
  sh5 <- find_sample(si, "TG_sh5")
  ntc0 <- find_sample(si, "NTC", "NTC_rep0")
  ntc1 <- find_sample(si, "NTC", "NTC_rep1")

  de_list <- list()
  de_list$TG_sh1_vs_NTC_rep0 <- pairwise_de(log_mat, sh1, ntc0, "TG_sh1_vs_NTC_rep0")
  de_list$TG_sh5_vs_NTC_rep0 <- pairwise_de(log_mat, sh5, ntc0, "TG_sh5_vs_NTC_rep0")
  de_list$TG_sh1_vs_NTC_rep1 <- pairwise_de(log_mat, sh1, ntc1, "TG_sh1_vs_NTC_rep1")
  de_list$TG_sh5_vs_NTC_rep1 <- pairwise_de(log_mat, sh5, ntc1, "TG_sh5_vs_NTC_rep1")
  de_list$TGsh_mean_vs_NTC <- mean_kd_vs_ntc_de(log_mat, si)
  de_list$common_up_vs_NTC_rep0 <- build_common(
    de_list$TG_sh1_vs_NTC_rep0, de_list$TG_sh5_vs_NTC_rep0, "up"
  )
  de_list$common_up_vs_NTC_rep1 <- build_common(
    de_list$TG_sh1_vs_NTC_rep1, de_list$TG_sh5_vs_NTC_rep1, "up"
  )
  de_list <- de_list[!vapply(de_list, is.null, logical(1))]

  gsea_rank <- list(
    common_up_vs_NTC_rep0 = full_rank_two(de_list$TG_sh1_vs_NTC_rep0, de_list$TG_sh5_vs_NTC_rep0),
    common_up_vs_NTC_rep1 = full_rank_two(de_list$TG_sh1_vs_NTC_rep1, de_list$TG_sh5_vs_NTC_rep1)
  )

  assign("log_mat", log_mat, envir = .GlobalEnv)
  assign("si", si, envir = .GlobalEnv)
  assign("norm", norm, envir = .GlobalEnv)
  assign("go_tab", go_tab, envir = .GlobalEnv)
  assign("go_sets", go_sets, envir = .GlobalEnv)
  assign("pathway_genes", pathway_genes, envir = .GlobalEnv)

  for (nm in names(de_list)) {
    volcano_df <- de_list[[nm]]
    if (nm %in% names(gsea_rank) && !is.null(gsea_rank[[nm]])) volcano_df <- gsea_rank[[nm]]
    tryCatch(
      analyze_one_comparison(
        nm, de_list[[nm]], volcano_df, norm$heat_mat, si, go_tab, go_sets, pathway_genes
      ),
      error = function(e) log_msg("ERROR in comparison ", nm, ": ", e$message)
    )
  }

  for (p_tag in names(p_cutoffs)) {
    tryCatch(plot_venn_up(
      de_list$TG_sh1_vs_NTC_rep0, de_list$TG_sh5_vs_NTC_rep0,
      file.path(result_dir, "common_up_vs_NTC_rep0"),
      "TG_sh1_vs_NTC_rep0", "TG_sh5_vs_NTC_rep0", "Common up vs NTC_rep0",
      unname(p_cutoffs[[p_tag]]), p_tag, has_real_pvalue(de_list$TG_sh1_vs_NTC_rep0)
    ), error = function(e) log_msg("venn NTC_rep0 error: ", e$message))
    tryCatch(plot_venn_up(
      de_list$TG_sh1_vs_NTC_rep1, de_list$TG_sh5_vs_NTC_rep1,
      file.path(result_dir, "common_up_vs_NTC_rep1"),
      "TG_sh1_vs_NTC_rep1", "TG_sh5_vs_NTC_rep1", "Common up vs NTC_rep1",
      unname(p_cutoffs[[p_tag]]), p_tag, has_real_pvalue(de_list$TG_sh1_vs_NTC_rep1)
    ), error = function(e) log_msg("venn NTC_rep1 error: ", e$message))
  }

  base::writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
  log_msg("Comparisons 1-4 done. Results in: ", result_dir)
}

if (!isTRUE(getOption("tg.rnaseq.functions_only", FALSE))) {
  run_comparisons_1_to_4()
}
