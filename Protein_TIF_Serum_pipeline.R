#!/usr/bin/env Rscript
# =============================================================================
# TIF / 血清蛋白质组：T vs N
# 输入：DIA-NN TIF_report.pg_matrix、Serum_report.pg_matrix
# 分析：
#   1) 组织间质液 T vs N：差异蛋白、火山图、上调 GO、上调 KEGG
#   2) 血清 T vs N：同上
#   3) TIF T vs N 有、血清 T vs N 无的蛋白，并绘制排名图
# T6 不并入 T vs N；TIF 与血清分开标准化，不混样本。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")
options(clusterProfiler.download.method = "auto")

# -----------------------------------------------------------------------------
# 0. 依赖包
# -----------------------------------------------------------------------------
cran_required <- c(
  "dplyr", "tidyr", "tibble", "stringr", "ggplot2", "ggrepel",
  "pheatmap", "RColorBrewer", "writexl"
)
bioc_required <- c(
  "limma", "clusterProfiler", "org.Hs.eg.db", "enrichplot",
  "AnnotationDbi"
)
bioc_optional <- c("pathview")

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
install_if_missing(bioc_required, bioc = TRUE, required = TRUE)
install_if_missing(bioc_optional, bioc = TRUE, required = FALSE)

safe_library <- function(pkgs) {
  for (p in pkgs) {
    if (requireNamespace(p, quietly = TRUE)) {
      suppressPackageStartupMessages(library(p, character.only = TRUE))
    }
  }
}
safe_library(c(cran_required, bioc_required, bioc_optional))
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

# -----------------------------------------------------------------------------
# 1. 路径与参数
# -----------------------------------------------------------------------------
resolve_project_dir <- function() {
  env_dir <- Sys.getenv("PROTEIN_TIF_SERUM_DIR", unset = "")
  candidates <- c(
    env_dir,
    "E:/R/Protein TIF serum",
    "E:\\R\\Protein TIF serum",
    getwd()
  )
  candidates <- unique(candidates[nzchar(candidates)])
  needed <- c("TIF_report.pg_matrix", "Serum_report.pg_matrix")
  for (d in candidates) {
    if (dir.exists(d) && all(file.exists(file.path(d, needed)))) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  for (d in candidates) {
    if (dir.exists(d) && any(file.exists(file.path(d, needed)))) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

project_dir <- resolve_project_dir()
result_dir  <- file.path(project_dir, "results_protein")
log_dir     <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

p_cutoff  <- 0.01
fc_ora    <- 1.5
rank_plot_n <- 50

anno_cols <- c(
  "Protein.Group", "Protein.Groups", "Protein.Ids", "Protein.ID",
  "Protein.Names", "Protein.Name", "Genes", "Gene", "Gene.Name",
  "Gene.Names", "First.Protein.Description", "Description",
  "N.Sequences", "N.Proteotypic.Sequences", "N.Proteotypic.Peptides",
  "Protein.Q.Value", "Global.Q.Value", "Global.PG.Q.Value",
  "PG.Q.Value", "Protein.Group.Q.Value"
)

# -----------------------------------------------------------------------------
# 2. 样本名识别：先 T6，再 T；去掉 TIF/SERUM 前缀
# -----------------------------------------------------------------------------
tokenize_sample <- function(colname) {
  b <- gsub("\\\\", "/", colname)
  b <- basename(b)
  b <- sub("\\.(raw|mzml|wiff|dia|parquet|tsv|txt)$", "", b, ignore.case = TRUE)
  s <- toupper(gsub("[^A-Za-z0-9]+", "_", b))
  s <- gsub("^_|_$", "", s)
  tokens <- unlist(strsplit(s, "_", fixed = TRUE))
  tokens <- tokens[nzchar(tokens)]
  drop <- c(
    "TIF", "SERUM", "REPORT", "PG", "PR", "MATRIX", "DIANN", "DIA", "NN",
    "SAMPLE", "REP", "REPLICATE", "INTENSITY", "LFQ", "MAXLFQ", "PROTEIN"
  )
  tokens[!tokens %in% drop]
}

classify_sample <- function(colname) {
  tokens <- tokenize_sample(colname)
  if (length(tokens) == 0) {
    return(list(group = NA_character_, replicate = NA_character_))
  }

  t6 <- grep("^T6", tokens, value = TRUE)
  if (length(t6) > 0) {
    lab <- t6[[1]]
    replicate <- "T6"
    if (lab %in% c("T61", "T6_1") || grepl("T6.*1$", lab)) replicate <- "T6-1"
    if (lab %in% c("T62", "T6_2") || grepl("T6.*2$", lab)) replicate <- "T6-2"
    if (lab %in% c("T63", "T6_3") || grepl("T6.*3$", lab)) replicate <- "T6-3"
    idx <- match(lab, tokens)
    if (lab == "T6" && idx < length(tokens) && tokens[idx + 1] %in% c("1", "2", "3")) {
      replicate <- paste0("T6-", tokens[idx + 1])
    }
    return(list(group = "T6", replicate = replicate))
  }

  ntok <- grep("^N([137])?$", tokens, value = TRUE)
  if (length(ntok) > 0) {
    lab <- ntok[[1]]
    replicate <- if (identical(lab, "N")) "N" else lab
    return(list(group = "N", replicate = replicate))
  }

  ttok <- grep("^T([135])?$", tokens, value = TRUE)
  if (length(ttok) > 0) {
    lab <- ttok[[1]]
    replicate <- if (identical(lab, "T")) "T" else lab
    return(list(group = "T", replicate = replicate))
  }

  list(group = NA_character_, replicate = NA_character_)
}

apply_optional_sample_map <- function(sample_info, assay) {
  map_path <- file.path(project_dir, "sample_map.csv")
  if (!file.exists(map_path)) return(sample_info)
  mp <- tryCatch(
    utils::read.csv(map_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(mp) || !all(c("file", "assay", "group") %in% names(mp))) {
    log_msg("sample_map.csv 存在但缺少 file/assay/group 列，忽略")
    return(sample_info)
  }
  mp <- mp[toupper(mp$assay) == toupper(assay), , drop = FALSE]
  if (nrow(mp) == 0) return(sample_info)
  for (i in seq_len(nrow(sample_info))) {
    hit <- mp$file[vapply(mp$file, function(x) grepl(x, sample_info$column[i], ignore.case = TRUE), logical(1))]
    if (length(hit) == 0) next
    row <- mp[match(hit[[1]], mp$file), ]
    sample_info$group[i] <- as.character(row$group)
    if ("replicate" %in% names(row) && nzchar(as.character(row$replicate))) {
      sample_info$replicate[i] <- as.character(row$replicate)
    }
  }
  sample_info
}

# -----------------------------------------------------------------------------
# 3. 读入 DIA-NN pg_matrix
# -----------------------------------------------------------------------------
first_symbol <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  parts <- unlist(strsplit(as.character(x), "[;|,/]+"))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  parts <- parts[!grepl("^XLOC_", parts, ignore.case = TRUE)]
  if (length(parts) == 0) return(NA_character_)
  parts[[1]]
}

first_uniprot <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  parts <- unlist(strsplit(as.character(x), "[;|,]+"))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) return(NA_character_)
  parts[[1]]
}

pick_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(NULL)
  df[[hit[[1]]]]
}

read_table_auto <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE)
  sep <- if (grepl("\t", first, fixed = TRUE)) "\t" else ","
  utils::read.delim(
    path, sep = sep, check.names = FALSE, stringsAsFactors = FALSE,
    comment.char = "", na.strings = c("", "NA", "NaN", "Filtered")
  )
}

is_mostly_numeric <- function(x) {
  xn <- suppressWarnings(as.numeric(as.character(x)))
  mean(!is.na(xn) | is.na(x) | as.character(x) %in% c("", "0")) > 0.7
}

read_pg_matrix <- function(path, assay) {
  if (!file.exists(path)) stop("找不到文件: ", path)
  raw <- read_table_auto(path)
  log_msg(assay, " raw dim: ", nrow(raw), " x ", ncol(raw), " from ", basename(path))

  sample_idx <- vapply(seq_along(raw), function(i) {
    nm <- names(raw)[i]
    if (nm %in% anno_cols) return(FALSE)
    if (grepl("q\\.value|qvalue|unique\\.peptides|sequences", nm, ignore.case = TRUE)) {
      return(FALSE)
    }
    is_mostly_numeric(raw[[i]])
  }, logical(1))

  if (!any(sample_idx)) stop(assay, ": 未识别到样品强度列")

  intens <- as.data.frame(lapply(raw[sample_idx], function(x) {
    xn <- suppressWarnings(as.numeric(as.character(x)))
    xn[is.finite(xn) & xn <= 0] <- NA_real_
    xn
  }), check.names = FALSE)

  protein_group <- pick_col(raw, c("Protein.Group", "Protein.Groups", "Protein.Ids", "Protein.ID"))
  protein_ids   <- pick_col(raw, c("Protein.Ids", "Protein.ID", "Protein.Group"))
  genes_raw     <- pick_col(raw, c("Genes", "Gene", "Gene.Names", "Gene.Name"))
  desc          <- pick_col(raw, c("First.Protein.Description", "Description", "Protein.Names", "Protein.Name"))

  if (is.null(protein_group)) protein_group <- paste0("PG_", seq_len(nrow(raw)))
  if (is.null(protein_ids)) protein_ids <- protein_group
  if (is.null(genes_raw)) genes_raw <- rep(NA_character_, nrow(raw))
  if (is.null(desc)) desc <- rep(NA_character_, nrow(raw))

  gene <- vapply(genes_raw, first_symbol, character(1), USE.NAMES = FALSE)
  uniprot <- vapply(protein_ids, first_uniprot, character(1), USE.NAMES = FALSE)
  display <- ifelse(!is.na(gene) & nzchar(gene), gene, uniprot)
  display <- ifelse(!is.na(display) & nzchar(display), display, as.character(protein_group))

  dup <- duplicated(display) | duplicated(display, fromLast = TRUE)
  display[dup] <- paste0(display[dup], "_", uniprot[dup])
  display <- make.unique(display)

  rownames(intens) <- display
  sample_info <- data.frame(
    column = names(intens),
    assay = assay,
    stringsAsFactors = FALSE
  )
  cls <- lapply(sample_info$column, classify_sample)
  sample_info$group <- vapply(cls, function(x) x$group, character(1))
  sample_info$replicate <- vapply(cls, function(x) x$replicate, character(1))
  sample_info$sample <- ifelse(
    !is.na(sample_info$replicate) & nzchar(sample_info$replicate),
    paste(assay, sample_info$replicate, sep = "_"),
    paste(assay, seq_len(nrow(sample_info)), sep = "_")
  )
  sample_info <- apply_optional_sample_map(sample_info, assay)
  sample_info$sample <- ifelse(
    !is.na(sample_info$replicate) & nzchar(sample_info$replicate),
    paste(assay, sample_info$replicate, sep = "_"),
    sample_info$sample
  )
  sample_info$sample <- make.unique(sample_info$sample)
  colnames(intens) <- sample_info$sample

  anno <- data.frame(
    protein = display,
    gene = gene,
    protein_group = as.character(protein_group),
    uniprot = uniprot,
    description = as.character(desc),
    stringsAsFactors = FALSE
  )
  rownames(anno) <- display

  list(mat = as.matrix(intens), sample_info = sample_info, anno = anno, assay = assay)
}

# -----------------------------------------------------------------------------
# 4. 过滤、log2、中位数标准化（每个矩阵单独做）
# -----------------------------------------------------------------------------
filter_and_normalize <- function(obj) {
  mat <- obj$mat
  si <- obj$sample_info
  keep_si <- si$group %in% c("N", "T", "T6")
  if (!any(keep_si)) {
    log_msg(obj$assay, ": 未能从列名识别 N/T/T6，保留全部样品但无法做 T vs N")
  } else {
    mat <- mat[, keep_si, drop = FALSE]
    si <- si[keep_si, , drop = FALSE]
  }

  n_valid <- rowSums(is.finite(mat))
  keep <- n_valid >= max(1, ceiling(0.3 * ncol(mat)))
  log_msg(obj$assay, ": filter proteins ", sum(keep), " / ", length(keep),
          " (>=30% valid values)")
  mat <- mat[keep, , drop = FALSE]
  anno <- obj$anno[rownames(mat), , drop = FALSE]

  log_mat <- log2(mat)
  col_med <- apply(log_mat, 2, stats::median, na.rm = TRUE)
  global_med <- stats::median(col_med, na.rm = TRUE)
  log_mat <- sweep(log_mat, 2, col_med - global_med, "-")

  list(log_mat = log_mat, sample_info = si, anno = anno, assay = obj$assay)
}

impute_left_na <- function(log_mat) {
  out <- log_mat
  for (i in seq_len(nrow(out))) {
    miss <- !is.finite(out[i, ])
    if (!any(miss) || all(miss)) next
    obs <- out[i, !miss]
    fill <- min(obs, na.rm = TRUE) - 1
    out[i, miss] <- fill
  }
  out
}

# -----------------------------------------------------------------------------
# 5. T vs N 差异分析
# -----------------------------------------------------------------------------
tn_samples <- function(sample_info, group) {
  sample_info$sample[sample_info$group == group]
}

de_t_vs_n <- function(norm) {
  si <- norm$sample_info
  t_samp <- tn_samples(si, "T")
  n_samp <- tn_samples(si, "N")
  if (length(t_samp) == 0 || length(n_samp) == 0) {
    log_msg(norm$assay, ": 缺少 T 或 N 样品，跳过差异分析。当前组别: ",
            paste(unique(si$group), collapse = ", "))
    return(NULL)
  }
  if (any(si$group == "T6")) {
    log_msg(norm$assay, ": T6 样品已识别但不进入 T vs N: ",
            paste(tn_samples(si, "T6"), collapse = ", "))
  }

  use <- c(t_samp, n_samp)
  mat <- norm$log_mat[, use, drop = FALSE]
  min_t <- if (length(t_samp) >= 3) 2 else 1
  min_n <- if (length(n_samp) >= 3) 2 else 1
  keep <- rowSums(is.finite(mat[, t_samp, drop = FALSE])) >= min_t &
    rowSums(is.finite(mat[, n_samp, drop = FALSE])) >= min_n
  mat <- mat[keep, , drop = FALSE]
  anno <- norm$anno[rownames(mat), , drop = FALSE]
  log_msg(norm$assay, " T vs N proteins after valid-value filter: ", nrow(mat),
          " | T=", paste(t_samp, collapse = ","),
          " | N=", paste(n_samp, collapse = ","))

  mean_t <- rowMeans(mat[, t_samp, drop = FALSE], na.rm = TRUE)
  mean_n <- rowMeans(mat[, n_samp, drop = FALSE], na.rm = TRUE)
  log2FC <- mean_t - mean_n
  ave <- (mean_t + mean_n) / 2
  have_rep <- length(t_samp) >= 2 && length(n_samp) >= 2

  pvalue <- rep(NA_real_, nrow(mat))
  padj <- rep(NA_real_, nrow(mat))
  tstat <- rep(NA_real_, nrow(mat))

  if (have_rep) {
    imp <- impute_left_na(mat)
    group <- factor(si$group[match(colnames(imp), si$sample)], levels = c("N", "T"))
    design <- stats::model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    fit <- limma::lmFit(imp, design)
    cont <- limma::makeContrasts(T_vs_N = T - N, levels = design)
    fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont), trend = TRUE, robust = TRUE)
    tt <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "none")
    tt <- tt[rownames(imp), , drop = FALSE]
    pvalue <- tt$P.Value
    padj <- tt$adj.P.Val
    tstat <- tt$t
    log_msg(norm$assay, " T vs N: limma with ", length(t_samp), " T vs ", length(n_samp), " N")
  } else {
    log_msg(norm$assay, " T vs N: 无生物学重复，只算 FC，不伪造 p 值")
  }

  de <- data.frame(
    protein = rownames(mat),
    gene = anno$gene,
    uniprot = anno$uniprot,
    protein_group = anno$protein_group,
    description = anno$description,
    log2FC = as.numeric(log2FC),
    FoldChange = 2^as.numeric(log2FC),
    AveExpr = as.numeric(ave),
    mean_T = as.numeric(mean_t),
    mean_N = as.numeric(mean_n),
    pvalue = pvalue,
    padj = padj,
    t = tstat,
    stringsAsFactors = FALSE
  )
  de$gene_key <- ifelse(!is.na(de$gene) & nzchar(de$gene), de$gene, de$protein)
  de <- de[order(de$log2FC, decreasing = TRUE), ]
  attr(de, "have_pvalue") <- have_rep
  attr(de, "t_samples") <- t_samp
  attr(de, "n_samples") <- n_samp
  attr(de, "log_mat") <- mat
  attr(de, "sample_info") <- si[si$sample %in% use, , drop = FALSE]
  de
}

select_up <- function(de) {
  have_p <- isTRUE(attr(de, "have_pvalue"))
  ok <- !is.na(de$log2FC) & de$log2FC > 0
  if (have_p) {
    ok <- ok & !is.na(de$pvalue) & de$pvalue < p_cutoff
  } else {
    ok <- ok & (de$FoldChange >= fc_ora)
    log_msg("无 p 值，上调蛋白改用 FC >= ", fc_ora)
  }
  de[ok, , drop = FALSE]
}

# -----------------------------------------------------------------------------
# 6. 绘图
# -----------------------------------------------------------------------------
save_gg <- function(plot, path_stub, width = 8, height = 6) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".pdf"), plot, width = width, height = height),
           error = function(e) log_msg("pdf ggsave failed: ", e$message))
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".png"), plot, width = width, height = height, dpi = 300),
           error = function(e) log_msg("png ggsave failed: ", e$message))
}

note_empty <- function(stub, msg) {
  dir.create(dirname(stub), recursive = TRUE, showWarnings = FALSE)
  writeLines(msg, paste0(stub, "_EMPTY.txt"))
}

plot_volcano <- function(de, highlight, title, outfile) {
  df <- de
  have_p <- "pvalue" %in% names(df) && any(!is.na(df$pvalue))
  if (have_p) {
    df$y <- -log10(pmax(df$pvalue, 1e-300))
    ylab <- "-log10(p value)"
    hline <- -log10(p_cutoff)
  } else {
    df$y <- df$AveExpr
    ylab <- "Average log2 intensity (no p-value)"
    hline <- NULL
  }
  df$set <- ifelse(df$protein %in% highlight, "up", "other")
  lab_src <- utils::head(highlight, 15)
  df$label <- ifelse(df$protein %in% lab_src, df$gene_key, NA)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = y, color = set)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.4) +
    ggplot2::scale_color_manual(values = c(other = "grey70", up = "#D62828")) +
    ggplot2::geom_vline(xintercept = c(-log2(fc_ora), log2(fc_ora)), linetype = 2, color = "grey40") +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, x = "log2 Fold Change (T / N)", y = ylab, color = NULL)
  if (!is.null(hline)) p <- p + ggplot2::geom_hline(yintercept = hline, linetype = 2, color = "grey40")
  save_gg(p, outfile)
}

plot_heatmap <- function(de, up, title, outfile) {
  mat <- attr(de, "log_mat")
  si <- attr(de, "sample_info")
  genes <- intersect(up$protein, rownames(mat))
  if (length(genes) > 200) genes <- genes[seq_len(200)]
  if (length(genes) < 2) {
    log_msg("Heatmap skipped (<2 proteins): ", title)
    return(invisible(NULL))
  }
  sub <- mat[genes, , drop = FALSE]
  rownames(sub) <- de$gene_key[match(rownames(sub), de$protein)]
  ann <- data.frame(Group = si$group, row.names = si$sample)
  pal <- c(N = "#4C78A8", T = "#E45756", T6 = "#54A24B")
  draw_hm <- function() {
    pheatmap::pheatmap(
      sub, scale = "row", annotation_col = ann,
      annotation_colors = list(Group = pal[names(pal) %in% unique(ann$Group)]),
      show_rownames = nrow(sub) <= 80, fontsize_row = 6, main = title,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100)
    )
  }
  grDevices::pdf(paste0(outfile, ".pdf"), width = 8, height = max(6, min(18, 0.18 * nrow(sub) + 3)))
  on.exit({
    while (grDevices::dev.cur() > 1) grDevices::dev.off()
  }, add = TRUE)
  tryCatch(draw_hm(), error = function(e) log_msg("heatmap pdf failed: ", e$message))
  grDevices::dev.off()
  grDevices::png(paste0(outfile, ".png"), width = 2400, height = max(1800, 40 * nrow(sub) + 400), res = 300)
  tryCatch(draw_hm(), error = function(e) log_msg("heatmap png failed: ", e$message))
  grDevices::dev.off()
}

plot_rank <- function(df, title, outfile, n = rank_plot_n) {
  if (nrow(df) == 0) {
    note_empty(outfile, "no TIF-specific proteins")
    return(invisible(NULL))
  }
  x <- df[order(df$log2FC, decreasing = TRUE), , drop = FALSE]
  x$rank <- seq_len(nrow(x))
  plot_df <- utils::head(x, n)
  plot_df$label <- factor(plot_df$gene_key, levels = rev(plot_df$gene_key))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = log2FC, y = label)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = log2FC, y = label, yend = label), color = "grey70") +
    ggplot2::geom_point(size = 2.6, color = "#D62828") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Showing top ", nrow(plot_df), " / ", nrow(x), " by TIF log2FC"),
      x = "TIF log2 Fold Change (T / N)",
      y = "Ranked protein"
    )
  save_gg(p, outfile, width = 9, height = max(6, min(18, 0.22 * nrow(plot_df) + 2)))
  p2 <- ggplot2::ggplot(x, ggplot2::aes(x = rank, y = log2FC)) +
    ggplot2::geom_line(color = "#D62828") +
    ggplot2::geom_point(size = 1.1, color = "#D62828") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = paste(title, "| rank vs log2FC"), x = "Rank (1 = highest TIF FC)", y = "TIF log2FC")
  save_gg(p2, paste0(outfile, "_rank_vs_log2FC"), width = 8, height = 5)
}

plot_pca <- function(log_mat, sample_info, outfile) {
  if (ncol(log_mat) < 2) return(invisible(NULL))
  keep <- rowSums(is.finite(log_mat)) == ncol(log_mat)
  if (sum(keep) < 3) {
    tmp <- impute_left_na(log_mat)
  } else {
    tmp <- log_mat[keep, , drop = FALSE]
  }
  pca <- tryCatch(stats::prcomp(t(tmp), scale. = TRUE), error = function(e) NULL)
  if (is.null(pca)) return(invisible(NULL))
  df <- data.frame(
    pca$x[, 1:2, drop = FALSE],
    group = sample_info$group[match(rownames(pca$x), sample_info$sample)],
    sample = rownames(pca$x)
  )
  varp <- summary(pca)$importance[2, 1:min(2, ncol(summary(pca)$importance))] * 100
  p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, color = group, label = sample)) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = "PCA after log2 median-centering",
      x = sprintf("PC1 (%.1f%%)", varp[1]),
      y = sprintf("PC2 (%.1f%%)", varp[2])
    )
  save_gg(p, outfile)
}

# -----------------------------------------------------------------------------
# 7. 上调蛋白 GO / KEGG（ORA）
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
      clusterProfiler::bitr(symbols, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
      error = function(e) data.frame(UNIPROT = character(), ENTREZID = character())
    )
    if (nrow(m) > 0) names(m)[1] <- "SYMBOL"
  }
  if (nrow(m) == 0) return(data.frame(gene = character(), entrez = character()))
  m <- m[!duplicated(m[[1]]), ]
  data.frame(gene = m[[1]], entrez = m[[2]], stringsAsFactors = FALSE)
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

try_save_plot <- function(fun, stub, width = 9, height = 7) {
  p <- tryCatch(fun(), error = function(e) {
    log_msg("Plot failed (", basename(stub), "): ", e$message)
    NULL
  })
  if (is.null(p)) return(invisible(FALSE))
  tryCatch({
    save_gg(p, stub, width = width, height = height)
    TRUE
  }, error = function(e) {
    log_msg("ggsave failed (", basename(stub), "): ", e$message)
    FALSE
  })
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

run_up_ora <- function(up, outdir, label) {
  go_dir <- file.path(outdir, "GO")
  kg_dir <- file.path(outdir, "KEGG")
  dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(kg_dir, recursive = TRUE, showWarnings = FALSE)

  genes <- unique(c(up$gene[!is.na(up$gene) & nzchar(up$gene)], up$uniprot))
  mp <- map_to_entrez(genes)
  entrez <- unique(mp$entrez)
  fc_sym <- setNames(up$log2FC, up$gene_key)
  if (length(entrez) < 3) {
    log_msg("ORA skipped, mapped genes < 3: ", label)
    note_empty(file.path(go_dir, "ORA_GO"), "too few mapped upregulated proteins")
    note_empty(file.path(kg_dir, "ORA_KEGG"), "too few mapped upregulated proteins")
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
    plot_ora_object(
      ego, file.path(go_dir, paste0("ORA_GO_", ont)),
      title_maybe_relaxed(ego, paste(label, "| upregulated ORA GO", ont)),
      fold_change = fc_sym
    )
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
    ek <- tryCatch(
      clusterProfiler::setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
      error = function(e) ek
    )
  }
  plot_ora_object(
    ek, file.path(kg_dir, "ORA_KEGG"),
    title_maybe_relaxed(ek, paste(label, "| upregulated ORA KEGG")),
    fold_change = fc_sym
  )
}

write_de_tables <- function(de, up, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(de, file.path(outdir, "DE_full.csv"), row.names = FALSE)
  utils::write.csv(up, file.path(outdir, "DE_upregulated.csv"), row.names = FALSE)
  tryCatch({
    writexl::write_xlsx(
      list(DE_full = de, upregulated = up),
      file.path(outdir, "DE_T_vs_N.xlsx")
    )
  }, error = function(e) log_msg("xlsx write failed: ", e$message))
}

analyze_tn <- function(norm, comp_name) {
  base <- file.path(result_dir, comp_name)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(paste("comparison:", comp_name),
      "T6 samples are excluded from T vs N.",
      "GO/KEGG are ORA on upregulated proteins only."),
    file.path(base, "00_README.txt")
  )
  plot_pca(norm$log_mat, norm$sample_info, file.path(base, "00_QC_PCA"))
  utils::write.csv(norm$sample_info, file.path(base, "sample_info.csv"), row.names = FALSE)

  de <- de_t_vs_n(norm)
  if (is.null(de)) {
    writeLines("missing T or N samples", file.path(base, "SKIPPED.txt"))
    return(NULL)
  }
  up <- select_up(de)
  log_msg(comp_name, " upregulated n = ", nrow(up), " / ", nrow(de))
  write_de_tables(de, up, base)
  plot_volcano(de, up$protein, paste(comp_name, "| volcano (up highlighted)"), file.path(base, "volcano_T_vs_N"))
  tryCatch(plot_heatmap(de, up, paste(comp_name, "| upregulated heatmap"), file.path(base, "heatmap_up_T_vs_N")),
           error = function(e) {
             while (grDevices::dev.cur() > 1) grDevices::dev.off()
             log_msg("heatmap failed: ", e$message)
           })
  run_up_ora(up, base, paste(comp_name, "| upregulated"))
  list(de = de, up = up, norm = norm)
}

# -----------------------------------------------------------------------------
# 8. TIF 有、血清无
# -----------------------------------------------------------------------------
tif_specific_proteins <- function(tif_res, serum_res) {
  outdir <- file.path(result_dir, "TIF_specific_vs_Serum")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(tif_res) || is.null(serum_res)) {
    writeLines("need both TIF and Serum T vs N results", file.path(outdir, "SKIPPED.txt"))
    return(invisible(NULL))
  }

  tif_de <- tif_res$de
  serum_de <- serum_res$de
  tif_up <- tif_res$up
  serum_up <- serum_res$up

  serum_keys <- unique(c(serum_de$gene_key, serum_de$gene, serum_de$uniprot))
  serum_keys <- serum_keys[!is.na(serum_keys) & nzchar(serum_keys)]
  serum_up_keys <- unique(c(serum_up$gene_key, serum_up$gene, serum_up$uniprot))
  serum_up_keys <- serum_up_keys[!is.na(serum_up_keys) & nzchar(serum_up_keys)]

  key_match <- function(df, pool) {
    if (nrow(df) == 0) return(logical(0))
    apply(df[, c("gene_key", "gene", "uniprot"), drop = FALSE], 1, function(v) {
      any(v[!is.na(v) & nzchar(v)] %in% pool)
    })
  }

  tif_up_not_serum_up <- tif_up[!key_match(tif_up, serum_up_keys), , drop = FALSE]
  tif_up_not_serum_detected <- tif_up[!key_match(tif_up, serum_keys), , drop = FALSE]
  tif_detected_not_serum <- tif_de[!key_match(tif_de, serum_keys), , drop = FALSE]

  if (nrow(tif_up_not_serum_up) > 0) {
    tif_up_not_serum_up <- tif_up_not_serum_up[order(tif_up_not_serum_up$log2FC, decreasing = TRUE), ]
    tif_up_not_serum_up$rank <- seq_len(nrow(tif_up_not_serum_up))
  }

  utils::write.csv(tif_up_not_serum_up, file.path(outdir, "TIF_up_not_in_Serum_up.csv"), row.names = FALSE)
  utils::write.csv(tif_up_not_serum_detected, file.path(outdir, "TIF_up_not_detected_in_Serum.csv"), row.names = FALSE)
  utils::write.csv(tif_detected_not_serum, file.path(outdir, "TIF_detected_not_in_Serum.csv"), row.names = FALSE)
  tryCatch(writexl::write_xlsx(
    list(
      TIF_up_not_Serum_up = tif_up_not_serum_up,
      TIF_up_not_detected_in_Serum = tif_up_not_serum_detected,
      TIF_detected_not_in_Serum = tif_detected_not_serum
    ),
    file.path(outdir, "TIF_specific_vs_Serum.xlsx")
  ), error = function(e) log_msg("TIF-specific xlsx failed: ", e$message))

  writeLines(
    c("TIF-specific definition:",
      "1) TIF_up_not_in_Serum_up: TIF T vs N upregulated, gene not upregulated in Serum T vs N (main ranking plot).",
      "2) TIF_up_not_detected_in_Serum: TIF upregulated and gene not quantified in Serum T vs N matrix.",
      "3) TIF_detected_not_in_Serum: quantified in TIF T/N but not in Serum T/N."),
    file.path(outdir, "00_README.txt")
  )

  plot_rank(
    tif_up_not_serum_up,
    "TIF-specific upregulated proteins (not upregulated in Serum T vs N)",
    file.path(outdir, "rankplot_TIF_specific_up")
  )
  if (nrow(tif_up_not_serum_up) == 0 && nrow(tif_detected_not_serum) > 0) {
    plot_rank(
      tif_detected_not_serum,
      "TIF-detected proteins absent from Serum T vs N",
      file.path(outdir, "rankplot_TIF_detected_not_in_Serum")
    )
  }
  log_msg("TIF-specific upregulated (not Serum up): ", nrow(tif_up_not_serum_up))
  log_msg("TIF upregulated not detected in Serum: ", nrow(tif_up_not_serum_detected))
}

# -----------------------------------------------------------------------------
# 9. 主流程
# -----------------------------------------------------------------------------
log_msg("Project dir: ", project_dir)
tif_path <- file.path(project_dir, "TIF_report.pg_matrix")
serum_path <- file.path(project_dir, "Serum_report.pg_matrix")
if (!file.exists(tif_path)) stop("缺少 TIF_report.pg_matrix: ", tif_path)
if (!file.exists(serum_path)) stop("缺少 Serum_report.pg_matrix: ", serum_path)

tif_raw <- read_pg_matrix(tif_path, "TIF")
serum_raw <- read_pg_matrix(serum_path, "Serum")
utils::write.csv(tif_raw$sample_info, file.path(log_dir, "TIF_sample_info_raw.csv"), row.names = FALSE)
utils::write.csv(serum_raw$sample_info, file.path(log_dir, "Serum_sample_info_raw.csv"), row.names = FALSE)

if (any(is.na(tif_raw$sample_info$group))) {
  log_msg("WARNING: TIF 有未识别样品，请检查列名或提供 sample_map.csv")
}
if (any(is.na(serum_raw$sample_info$group))) {
  log_msg("WARNING: 血清有未识别样品，请检查列名或提供 sample_map.csv")
}

tif_norm <- filter_and_normalize(tif_raw)
serum_norm <- filter_and_normalize(serum_raw)

tif_res <- tryCatch(analyze_tn(tif_norm, "TIF_T_vs_N"), error = function(e) {
  log_msg("TIF T vs N failed: ", e$message)
  NULL
})
serum_res <- tryCatch(analyze_tn(serum_norm, "Serum_T_vs_N"), error = function(e) {
  log_msg("Serum T vs N failed: ", e$message)
  NULL
})
tryCatch(tif_specific_proteins(tif_res, serum_res), error = function(e) {
  log_msg("TIF-specific comparison failed: ", e$message)
})

base::writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
log_msg("All done. Results in: ", result_dir)
