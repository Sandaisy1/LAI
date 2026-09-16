#!/usr/bin/env Rscript
# =============================================================================
# TIF / 血清蛋白质组：T vs N（独立脚本，不依赖 Protein_TIF_Serum_pipeline.R）
#
# 输入：E:/R/Protein TIF serum 下的 DIA-NN 蛋白矩阵
#   TIF_report.pg_matrix / TIF_report.pg_matrix.tsv
#   Serum_report.pg_matrix / Serum_report.pg_matrix.tsv
# 不要用 *.pr_matrix。不要把 T6 并进 T vs N。TIF 与血清分开标准化。
#
# 三件事：
#   1) 组织间质液 T vs N：差异蛋白、火山图、上调 GO、上调 KEGG
#   2) 血清 T vs N：同上
#   3) TIF T vs N 有、血清 T vs N 全表（上调+下调）都没有的蛋白 → 排名图
#
# 在 R / RStudio 控制台运行（不要输入 Rscript）：
#   setwd("E:/R/Protein TIF serum")
#   source("Protein_TIF_Serum_TVsN_standalone.R")
#
# 结果写到 results_protein_standalone/ ，不会覆盖 results_protein/ 。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")
options(clusterProfiler.download.method = "auto")

# -----------------------------------------------------------------------------
# 0. 依赖包
# -----------------------------------------------------------------------------
cran_required <- c("ggplot2", "ggrepel", "writexl")
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
install_if_missing("limma", bioc = TRUE, required = TRUE)
if (Sys.getenv("PROTEIN_SKIP_OPTIONAL_BIOC", unset = "0") != "1") {
  install_if_missing(
    c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "AnnotationDbi", "GO.db"),
    bioc = TRUE, required = FALSE
  )
}

safe_library <- function(pkgs) {
  for (p in pkgs) {
    if (requireNamespace(p, quietly = TRUE)) {
      suppressPackageStartupMessages(library(p, character.only = TRUE))
    }
  }
}
safe_library(c(cran_required, "limma", "clusterProfiler", "org.Hs.eg.db",
               "enrichplot", "AnnotationDbi", "GO.db"))
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

# -----------------------------------------------------------------------------
# 1. 路径与参数
# -----------------------------------------------------------------------------
THIS_SCRIPT <- "Protein_TIF_Serum_TVsN_standalone.R"

script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]),
                                 winslash = "/", mustWork = FALSE)))
  }
  n <- tryCatch(sys.nframe(), error = function(e) 0L)
  if (is.finite(n) && n >= 1) {
    for (i in n:1) {
      ofile <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
      if (!is.null(ofile) && nzchar(ofile) && file.exists(ofile)) {
        return(dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)))
      }
    }
  }
  if (file.exists(THIS_SCRIPT)) {
    return(normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
})

find_pg_matrix <- function(dir, assay = c("TIF", "Serum")) {
  assay <- toupper(match.arg(assay))
  if (!dir.exists(dir)) return(NA_character_)
  stems <- c(
    paste0(assay, "_report.pg_matrix"),
    paste0(tolower(assay), "_report.pg_matrix"),
    paste0(assay, ".pg_matrix")
  )
  exts <- c("", ".tsv", ".txt", ".csv")
  for (stem in stems) {
    for (ext in exts) {
      p <- file.path(dir, paste0(stem, ext))
      if (file.exists(p) && !grepl("pr[_\\.]?matrix", basename(p), ignore.case = TRUE)) {
        return(normalizePath(p, winslash = "/", mustWork = FALSE))
      }
    }
  }
  files <- list.files(dir, full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) return(NA_character_)
  bn <- basename(files)
  hit <- grepl("pg([._])?matrix", bn, ignore.case = TRUE) &
    !grepl("pr([._])?matrix", bn, ignore.case = TRUE) &
    grepl(assay, bn, ignore.case = TRUE)
  if (!any(hit)) return(NA_character_)
  prefer <- hit & grepl("report", bn, ignore.case = TRUE)
  pick <- if (any(prefer)) which(prefer)[1] else which(hit)[1]
  normalizePath(files[pick], winslash = "/", mustWork = FALSE)
}

resolve_project_dir <- function() {
  env_dir <- Sys.getenv("PROTEIN_TIF_SERUM_DIR", unset = "")
  candidates <- c(
    env_dir,
    "E:/R/Protein TIF serum",
    "E:\\R\\Protein TIF serum",
    getwd(),
    script_dir,
    file.path(script_dir, "demo_protein_tif_serum"),
    file.path(getwd(), "demo_protein_tif_serum")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  score_dir <- function(d) {
    if (!dir.exists(d)) return(0L)
    as.integer(!is.na(find_pg_matrix(d, "TIF"))) + as.integer(!is.na(find_pg_matrix(d, "Serum")))
  }
  scores <- vapply(candidates, score_dir, integer(1))
  if (any(scores == 2L)) {
    return(normalizePath(candidates[which(scores == 2L)[1]], winslash = "/", mustWork = FALSE))
  }
  if (any(scores == 1L)) {
    return(normalizePath(candidates[which(scores == 1L)[1]], winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

project_dir <- resolve_project_dir()
result_dir  <- file.path(project_dir, "results_protein_standalone")
log_dir     <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("standalone_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(...))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

p_cutoff    <- 0.01
fc_ora      <- 1.5
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
# 6. 火山图
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

plot_rank <- function(df, title, outfile, n = rank_plot_n) {
  if (is.null(df) || nrow(df) == 0) {
    note_empty(outfile, "no TIF-specific proteins for ranking plot")
    return(invisible(NULL))
  }
  x <- df[order(df$log2FC, decreasing = TRUE), , drop = FALSE]
  x$rank <- seq_len(nrow(x))
  x$direction <- ifelse(!is.na(x$log2FC) & x$log2FC >= 0, "TIF T > N", "TIF T < N")
  plot_df <- utils::head(x, n)
  plot_df$label <- factor(plot_df$gene_key, levels = rev(unique(plot_df$gene_key)))
  col_vals <- c("TIF T > N" = "#D62828", "TIF T < N" = "#1D4E89")
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = log2FC, y = label, color = direction)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = log2FC, y = label, yend = label),
      linewidth = 0.6
    ) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_color_manual(values = col_vals) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Showing top ", nrow(plot_df), " / ", nrow(x),
                        " by TIF log2FC; serum T vs N (up or down) absent"),
      x = "TIF log2 Fold Change (T / N)",
      y = "Ranked protein (1 = highest TIF FC at top)",
      color = NULL
    )
  save_gg(p, outfile, width = 9, height = max(6, min(18, 0.22 * nrow(plot_df) + 2)))

  p2 <- ggplot2::ggplot(x, ggplot2::aes(x = rank, y = log2FC, color = direction)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    ggplot2::geom_line(color = "grey55") +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_color_manual(values = col_vals) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = paste(title, "| rank vs log2FC"),
      x = "Rank (1 = highest TIF FC)",
      y = "TIF log2FC (T / N)",
      color = NULL
    )
  save_gg(p2, paste0(outfile, "_rank_vs_log2FC"), width = 8, height = 5)
  invisible(x)
}

# -----------------------------------------------------------------------------
# 7. 上调蛋白 GO / KEGG（ORA）
# -----------------------------------------------------------------------------
map_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  empty <- data.frame(gene = character(), entrez = character())
  if (length(symbols) == 0 || !has_pkg("org.Hs.eg.db") || !has_pkg("AnnotationDbi")) return(empty)
  m <- tryCatch(
    AnnotationDbi::select(org.Hs.eg.db, keys = symbols, keytype = "SYMBOL", columns = "ENTREZID"),
    error = function(e) NULL
  )
  if (is.null(m) || !nrow(m) || all(is.na(m$ENTREZID))) {
    m <- tryCatch(
      AnnotationDbi::select(org.Hs.eg.db, keys = symbols, keytype = "UNIPROT", columns = "ENTREZID"),
      error = function(e) NULL
    )
    if (!is.null(m) && nrow(m) > 0) names(m)[names(m) %in% c("UNIPROT", "SYMBOL")][1] <- "SYMBOL"
  }
  if (is.null(m) || nrow(m) == 0) return(empty)
  m <- m[!is.na(m$ENTREZID) & nzchar(m$ENTREZID), , drop = FALSE]
  if (nrow(m) == 0) return(empty)
  m <- m[!duplicated(m[[1]]), ]
  data.frame(gene = as.character(m[[1]]), entrez = as.character(m$ENTREZID), stringsAsFactors = FALSE)
}

ora_hyper <- function(query, term2gene, min_size = 5, max_size = 500) {
  universe <- unique(as.character(term2gene$gene))
  query <- unique(intersect(as.character(query), universe))
  k <- length(query)
  n_uni <- length(universe)
  if (k < 3 || n_uni < 20) return(NULL)
  split_terms <- split(as.character(term2gene$gene), as.character(term2gene$term))
  rows <- lapply(names(split_terms), function(term) {
    genes <- unique(intersect(split_terms[[term]], universe))
    m <- length(genes)
    if (m < min_size || m > max_size) return(NULL)
    overlap <- intersect(query, genes)
    x <- length(overlap)
    if (x < 2) return(NULL)
    p <- stats::phyper(x - 1, m, n_uni - m, k, lower.tail = FALSE)
    data.frame(
      ID = term, Count = x, Size = m,
      GeneRatio = x / k, BgRatio = m / n_uni, pvalue = p,
      geneID = paste(overlap, collapse = "/"),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$p.adjust <- stats::p.adjust(df$pvalue, method = "BH")
  df[order(df$pvalue), , drop = FALSE]
}

go_term2gene <- function(ont) {
  xx <- as.list(org.Hs.egGO2ALLEGS)
  ont_map <- AnnotationDbi::Ontology(GO.db::GOTERM)
  keep <- names(xx)[!is.na(ont_map[names(xx)]) & ont_map[names(xx)] == ont]
  term <- rep(keep, lengths(xx[keep]))
  gene <- unlist(xx[keep], use.names = FALSE)
  data.frame(term = term, gene = as.character(gene), stringsAsFactors = FALSE)
}

kegg_term2gene <- function() {
  cache <- file.path(log_dir, "kegg_hsa_pathway_link.tsv")
  if (!file.exists(cache)) {
    ok <- tryCatch({
      utils::download.file("https://rest.kegg.jp/link/hsa/pathway", cache, quiet = TRUE, mode = "wb")
      TRUE
    }, error = function(e) {
      log_msg("KEGG download failed: ", e$message)
      FALSE
    })
    if (!ok) return(NULL)
  }
  raw <- tryCatch(utils::read.delim(cache, header = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(raw) || ncol(raw) < 2) return(NULL)
  names(raw)[1:2] <- c("term", "gene")
  raw$term <- sub("^path:", "", raw$term)
  raw$gene <- sub("^hsa:", "", raw$gene)
  raw
}

kegg_term_names <- function(ids) {
  cache <- file.path(log_dir, "kegg_hsa_pathway_list.tsv")
  if (!file.exists(cache)) {
    tryCatch(
      utils::download.file("https://rest.kegg.jp/list/pathway/hsa", cache, quiet = TRUE, mode = "wb"),
      error = function(e) log_msg("KEGG list download failed: ", e$message)
    )
  }
  nms <- setNames(ids, ids)
  if (!file.exists(cache)) return(nms)
  raw <- tryCatch(utils::read.delim(cache, header = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(raw) || ncol(raw) < 2) return(nms)
  raw[[1]] <- sub("^path:", "", raw[[1]])
  hit <- match(ids, raw[[1]])
  nms[!is.na(hit)] <- raw[[2]][hit[!is.na(hit)]]
  nms
}

go_term_names <- function(ids) {
  nms <- tryCatch(AnnotationDbi::Term(GO.db::GOTERM[ids]), error = function(e) NULL)
  if (is.null(nms)) return(setNames(ids, ids))
  out <- as.character(nms)
  names(out) <- names(nms)
  out[ids]
}

plot_simple_ora <- function(df, stub, title) {
  if (is.null(df) || nrow(df) == 0) {
    note_empty(stub, "no enrichment terms")
    return(invisible(NULL))
  }
  utils::write.csv(df, paste0(stub, ".csv"), row.names = FALSE)
  top <- utils::head(df, 15)
  top$Description <- factor(top$Description, levels = rev(unique(top$Description)))
  p_bar <- ggplot2::ggplot(top, ggplot2::aes(x = Count, y = Description, fill = -log10(pmax(p.adjust, 1e-300)))) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_gradient(low = "#FEE0D2", high = "#CB181D") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, x = "Gene count", y = NULL, fill = "-log10(padj)")
  save_gg(p_bar, paste0(stub, "_barplot"), 9, 7)
  p_dot <- ggplot2::ggplot(top, ggplot2::aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "#CB181D", high = "#2171B5") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, x = "GeneRatio", y = NULL)
  save_gg(p_dot, paste0(stub, "_dotplot"), 9, 7)
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
  if (!has_pkg("org.Hs.eg.db") || !has_pkg("AnnotationDbi")) {
    log_msg("ORA skipped (need org.Hs.eg.db): ", label)
    note_empty(file.path(go_dir, "ORA_GO"), "org.Hs.eg.db not installed")
    note_empty(file.path(kg_dir, "ORA_KEGG"), "org.Hs.eg.db not installed")
    return(invisible(NULL))
  }

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

  if (has_pkg("clusterProfiler") && has_pkg("enrichplot")) {
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
    return(invisible(TRUE))
  }

  log_msg("clusterProfiler 不可用，改用 org.Hs.eg.db 超几何检验做上调 GO/KEGG: ", label)
  if (!has_pkg("GO.db")) {
    note_empty(file.path(go_dir, "ORA_GO"), "GO.db not installed")
  } else {
    for (ont in c("BP", "MF", "CC")) {
      t2g <- tryCatch(go_term2gene(ont), error = function(e) {
        log_msg("GO term map failed (", ont, "): ", e$message)
        NULL
      })
      df <- if (is.null(t2g)) NULL else ora_hyper(entrez, t2g)
      if (!is.null(df) && nrow(df) > 0) {
        nms <- go_term_names(df$ID)
        df$Description <- ifelse(is.na(nms), df$ID, unname(nms))
      }
      plot_simple_ora(df, file.path(go_dir, paste0("ORA_GO_", ont)),
                      paste(label, "| upregulated ORA GO", ont))
    }
  }
  t2g_k <- tryCatch(kegg_term2gene(), error = function(e) {
    log_msg("KEGG term map failed: ", e$message)
    NULL
  })
  dfk <- if (is.null(t2g_k)) NULL else ora_hyper(entrez, t2g_k)
  if (!is.null(dfk) && nrow(dfk) > 0) {
    nms <- kegg_term_names(dfk$ID)
    dfk$Description <- ifelse(is.na(nms), dfk$ID, unname(nms))
  }
  plot_simple_ora(dfk, file.path(kg_dir, "ORA_KEGG"),
                  paste(label, "| upregulated ORA KEGG"))
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
      "Standalone script: Protein_TIF_Serum_TVsN_standalone.R",
      "T6 samples are excluded from T vs N.",
      "GO/KEGG are ORA on upregulated proteins only."),
    file.path(base, "00_README.txt")
  )
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
  run_up_ora(up, base, paste(comp_name, "| upregulated"))
  list(de = de, up = up, norm = norm)
}

# -----------------------------------------------------------------------------
# 8. 第3条：TIF T vs N 有、血清 T vs N 全表（上调+下调）无 → 排名图
# -----------------------------------------------------------------------------
norm_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("[-][0-9]+$", "", x)
  x[is.na(x) | !nzchar(x) | x %in% c("NA", "NAN")] <- NA_character_
  x
}

protein_id_pool <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(character(0))
  ids <- c(
    norm_id(df$gene),
    norm_id(df$gene_key),
    norm_id(df$uniprot),
    norm_id(df$protein_group)
  )
  unique(ids[!is.na(ids)])
}

row_in_pool <- function(df, pool) {
  if (nrow(df) == 0) return(logical(0))
  vapply(seq_len(nrow(df)), function(i) {
    any(protein_id_pool(df[i, , drop = FALSE]) %in% pool)
  }, logical(1))
}

tif_specific_rank <- function(tif_res, serum_res) {
  outdir <- file.path(result_dir, "TIF_specific_vs_Serum")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  unlink(list.files(outdir, full.names = TRUE))
  if (is.null(tif_res) || is.null(serum_res)) {
    writeLines("need both TIF and Serum T vs N tables", file.path(outdir, "SKIPPED.txt"))
    return(invisible(NULL))
  }

  tif_de <- tif_res$de
  serum_de <- serum_res$de
  serum_de_pool <- protein_id_pool(serum_de)
  in_serum <- row_in_pool(tif_de, serum_de_pool)
  tif_only <- tif_de[!in_serum, , drop = FALSE]
  tif_and_serum <- tif_de[in_serum, , drop = FALSE]
  if (nrow(tif_only) > 0) {
    tif_only <- tif_only[order(tif_only$log2FC, decreasing = TRUE), ]
    tif_only$rank <- seq_len(nrow(tif_only))
    tif_only$in_Serum_T_vs_N <- FALSE
    tif_only$TIF_direction <- ifelse(tif_only$log2FC >= 0, "TIF_T_gt_N", "TIF_T_lt_N")
  }
  if (nrow(tif_and_serum) > 0) {
    tif_and_serum$in_Serum_T_vs_N <- TRUE
  }

  utils::write.csv(tif_de, file.path(outdir, "00_TIF_T_vs_N_DE.csv"), row.names = FALSE)
  utils::write.csv(serum_de, file.path(outdir, "00_Serum_T_vs_N_DE.csv"), row.names = FALSE)
  utils::write.csv(tif_only, file.path(outdir, "TIF_present_not_in_Serum_TN.csv"), row.names = FALSE)
  utils::write.csv(tif_and_serum, file.path(outdir, "TIF_also_in_Serum_TN_excluded.csv"), row.names = FALSE)
  tryCatch(writexl::write_xlsx(
    list(
      TIF_T_vs_N_DE = tif_de,
      Serum_T_vs_N_DE = serum_de,
      TIF_present_not_in_Serum_TN = tif_only,
      overlap_excluded = tif_and_serum
    ),
    file.path(outdir, "TIF_present_not_in_Serum_TN.xlsx")
  ), error = function(e) log_msg("TIF-specific xlsx failed: ", e$message))

  writeLines(
    c("Item 3 (standalone ranking script):",
      "Set = proteins present in TIF T vs N DE table",
      "      minus proteins present in Serum T vs N DE table",
      "      (Serum exclusion uses the FULL table: up AND down).",
      "Match by gene / UniProt / protein group.",
      "Ranking plots: rank_TIF_present_absent_from_Serum_TN*.",
      "T6 is never used. TIF and Serum are never mixed for DE.",
      paste("TIF T vs N proteins:", nrow(tif_de)),
      paste("Serum T vs N proteins (up+down):", nrow(serum_de)),
      paste("Present in both (excluded):", nrow(tif_and_serum)),
      paste("TIF present, absent from Serum T vs N:", nrow(tif_only))),
    file.path(outdir, "00_README.txt")
  )

  tryCatch(
    plot_rank(
      tif_only,
      "TIF T vs N present, absent from Serum T vs N (up or down)",
      file.path(outdir, "rank_TIF_present_absent_from_Serum_TN")
    ),
    error = function(e) log_msg("ranking plot failed: ", e$message)
  )

  log_msg("Item 3: TIF T vs N n=", nrow(tif_de),
          " Serum T vs N n=", nrow(serum_de),
          " overlap excluded=", nrow(tif_and_serum),
          " TIF present not in Serum T vs N=", nrow(tif_only))
}

# -----------------------------------------------------------------------------
# 9. 主流程
# -----------------------------------------------------------------------------
log_msg("Standalone script (does not source Protein_TIF_Serum_pipeline.R)")
log_msg("Project dir: ", project_dir)
log_msg("Results dir: ", result_dir)
tif_path <- find_pg_matrix(project_dir, "TIF")
serum_path <- find_pg_matrix(project_dir, "Serum")
listed <- paste(list.files(project_dir), collapse = ", ")
if (is.na(tif_path) || !file.exists(tif_path)) {
  stop(
    "找不到 TIF 蛋白矩阵 (pg_matrix)。已查目录: ", project_dir, "\n",
    "当前文件: ", listed, "\n",
    "资源管理器若显示“TSV 文件”，实际文件名往往是 TIF_report.pg_matrix.tsv ，脚本现已支持该后缀。",
    call. = FALSE
  )
}
if (is.na(serum_path) || !file.exists(serum_path)) {
  stop(
    "找不到 Serum 蛋白矩阵 (pg_matrix)。已查目录: ", project_dir, "\n",
    "当前文件: ", listed, "\n",
    "资源管理器若显示“TSV 文件”，实际文件名往往是 Serum_report.pg_matrix.tsv 。",
    call. = FALSE
  )
}
log_msg("TIF matrix: ", tif_path)
log_msg("Serum matrix: ", serum_path)

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
tryCatch(tif_specific_rank(tif_res, serum_res), error = function(e) {
  log_msg("TIF-specific ranking failed: ", e$message)
})

base::writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
log_msg("All done. Results in: ", result_dir)
log_msg("R console: setwd(\"E:/R/Protein TIF serum\"); source(\"", THIS_SCRIPT, "\")")
