#!/usr/bin/env Rscript
# =============================================================================
# PCY_RNA：shP1 / shP2 相对 NTC 的基因差异表达
#
# 输入目录：E:/R/PCY_RNA
#   genes.read_group_tracking  （首选，每个重复的 raw_frags）
#   genes.count_tracking
#   genes.fpkm_tracking        （基因符号；没有 count 时也作为表达量）
#
# 比较（同一组别里的重复样品作为生物学重复，不拆开）：
#   1. shP1 vs NTC
#   2. shP2 vs NTC
# 火山图用全部基因。下调且 P < 0.05 的基因再做 GO、KEGG、Reactome 通路，
# 并画气泡图、富集热图，以及这组基因的表达热图。
#
# 运行：
#   source("E:/R/PCY_RNA/PCY_RNA.R", encoding = "UTF-8")
#   或
#   Rscript PCY_RNA.R "E:/R/PCY_RNA"
# 物种：人（org.Hs.eg.db）。log2FC = log2(sh / NTC)，负值表示 sh 低于 NTC。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")
options(clusterProfiler.download.method = "auto")

p_cutoff <- 0.05
lfc_cutoff <- 0

# -----------------------------------------------------------------------------
# 日志与小工具
# -----------------------------------------------------------------------------
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  message(msg)
  if (exists("log_file", envir = .GlobalEnv)) {
    cat(msg, "\n", file = get("log_file", envir = .GlobalEnv), append = TRUE)
  }
}

has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

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
    stop("缺少必需 R 包: ", paste(still, collapse = ", "), call. = FALSE)
  }
  if (length(still) > 0) {
    log_msg("可选包未安装，对应分析会跳过: ", paste(still, collapse = ", "))
  }
  invisible(TRUE)
}

this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(f) >= 1) {
    return(dirname(normalizePath(f[1], winslash = "/", mustWork = FALSE)))
  }
  for (i in rev(seq_len(sys.nframe()))) {
    env <- sys.frame(i)
    if (exists("ofile", envir = env, inherits = FALSE)) {
      ofile <- get("ofile", envir = env, inherits = FALSE)
      if (!is.null(ofile) && nzchar(ofile)) {
        return(dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)))
      }
    }
  }
  NA_character_
}

has_expression_file <- function(d) {
  if (!nzchar(d) || is.na(d) || !dir.exists(d)) return(FALSE)
  any(file.exists(file.path(d, c(
    "genes.read_group_tracking",
    "genes.count_tracking",
    "genes.fpkm_tracking"
  ))))
}

resolve_project_dir <- function() {
  cmd <- commandArgs(trailingOnly = TRUE)
  candidates <- c(
    Sys.getenv("PCY_RNA_DIR", unset = ""),
    if (length(cmd) >= 1) cmd[[1]] else "",
    "E:/R/PCY_RNA",
    "E:\\R\\PCY_RNA",
    this_script_dir(),
    getwd()
  )
  candidates <- unique(candidates[nzchar(candidates) & !is.na(candidates)])
  for (d in candidates) {
    if (has_expression_file(d)) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  stop(
    "未找到 Cuffdiff 表达文件。请把 genes.read_group_tracking、genes.count_tracking、",
    "genes.fpkm_tracking 放在 E:/R/PCY_RNA，或运行 Rscript PCY_RNA.R \"你的目录\"。已查找: ",
    paste(candidates, collapse = " | "),
    call. = FALSE
  )
}

save_gg <- function(plot, path_stub, width = 8, height = 6) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(path_stub, ".pdf"), plot, width = width, height = height)
  ggplot2::ggsave(paste0(path_stub, ".png"), plot, width = width, height = height, dpi = 300)
}

write_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, path, row.names = FALSE, fileEncoding = "UTF-8")
}

note_empty <- function(path, msg) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(msg, path)
}

# -----------------------------------------------------------------------------
# 样品名 → shP1 / shP2 / NTC
# -----------------------------------------------------------------------------
classify_pcy <- function(name) {
  raw <- toupper(trimws(as.character(name)))
  n <- gsub("[^A-Z0-9]", "", raw)
  if (!nzchar(n)) return(NA_character_)
  if (grepl("SHP2|SHRNAP2", n) || n %in% c("P2", "SH2")) return("shP2")
  if (grepl("SHP1|SHRNAP1", n) || n %in% c("P1", "SH1")) return("shP1")
  if (grepl("NTC|SHNC|NEGCTRL|CTRL|CONTROL", n) || grepl("^NC[0-9]*$", n)) return("NTC")
  NA_character_
}

# -----------------------------------------------------------------------------
# 基因符号：复合名取一个官方符号，没有符号的 XLOC 保留
# -----------------------------------------------------------------------------
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
  genes[is.na(genes) | genes == "" | genes == "-"] <- NA_character_
  keep <- !is.na(genes)
  mat <- mat[keep, , drop = FALSE]
  genes <- genes[keep]
  if (nrow(mat) == 0) return(mat)
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
  n_fused <- sum(grepl("[,;|/]", raw), na.rm = TRUE)
  mat <- collapse_by_gene(mat, genes)
  n_xloc <- sum(grepl("^(XLOC|TCONS|CUFF)_", rownames(mat), ignore.case = TRUE))
  log_msg(
    "基因名整理后 ", nrow(mat), " 个；复合名 ", n_fused,
    " 个；仍无官方符号的 XLOC ", n_xloc, " 个（留在差异表里，GO/KEGG 通常对不上）"
  )
  mat
}

# -----------------------------------------------------------------------------
# 读入
# -----------------------------------------------------------------------------
read_read_groups_info <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- lines[!grepl("^\\s*$", lines)]
  if (length(lines) == 0) return(NULL)
  lines <- sub("^#+\\s*", "", lines)
  tab <- strsplit(lines, "\t", fixed = TRUE)
  ncol_max <- max(lengths(tab))
  first <- tolower(tab[[1]])
  has_header <- any(first %in% c("file", "group", "condition", "replicate", "sample", "rep_name"))
  if (has_header) {
    hdr <- make.unique(tolower(gsub("[^a-z0-9]+", "_", tab[[1]])))
    rows <- tab[-1]
  } else {
    hdr <- paste0("V", seq_len(ncol_max))
    rows <- tab
  }
  if (length(rows) == 0) return(NULL)
  mat <- do.call(rbind, lapply(rows, function(x) {
    length(x) <- ncol_max
    x
  }))
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  names(df) <- hdr[seq_len(ncol(df))]
  df
}

condition_from_info <- function(info) {
  if (is.null(info) || nrow(info) == 0) return(NULL)
  nms <- names(info)
  if ("condition" %in% nms) return(as.character(info$condition))
  scored <- vapply(nms, function(nm) {
    sum(!is.na(vapply(info[[nm]], classify_pcy, character(1))))
  }, numeric(1))
  if (length(scored) == 0 || max(scored) == 0) return(NULL)
  as.character(info[[names(which.max(scored))]])
}

pivot_tracking <- function(id, sample, value) {
  value <- suppressWarnings(as.numeric(value))
  value[!is.finite(value)] <- 0
  key <- paste(id, sample, sep = "\r")
  if (any(duplicated(key))) {
    log_msg("tracking_id 与样品有重复行，同一格取平均")
    ord <- order(key)
    id <- id[ord]
    sample <- sample[ord]
    value <- value[ord]
    key <- key[ord]
    keep <- !duplicated(key)
    value <- ave(value, key, FUN = mean)[keep]
    id <- id[keep]
    sample <- sample[keep]
  }
  genes <- unique(id)
  samples <- unique(sample)
  mat <- matrix(0, nrow = length(genes), ncol = length(samples), dimnames = list(genes, samples))
  mat[cbind(match(id, genes), match(sample, samples))] <- value
  storage.mode(mat) <- "double"
  mat
}

read_read_group_tracking <- function(path) {
  rg <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  need <- c("tracking_id", "condition", "replicate")
  if (!all(need %in% names(rg))) {
    log_msg("genes.read_group_tracking 缺少列: ", paste(setdiff(need, names(rg)), collapse = ", "))
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
  log_msg("表达量列: ", value_col)
  rg$group <- vapply(as.character(rg$condition), classify_pcy, character(1))
  dropped <- unique(as.character(rg$condition[is.na(rg$group)]))
  if (length(dropped) > 0) {
    log_msg("未能识别为 shP1/shP2/NTC 的 condition: ", paste(dropped, collapse = ", "))
  }
  rg <- rg[!is.na(rg$group), , drop = FALSE]
  if (nrow(rg) == 0) return(NULL)
  rep_lab <- as.character(rg$replicate)
  rep_lab[!nzchar(rep_lab) | is.na(rep_lab)] <- "0"
  rg$sample <- make.unique(paste(rg$condition, rep_lab, sep = "_rep"), sep = "_")
  mat <- pivot_tracking(rg$tracking_id, rg$sample, rg[[value_col]])
  gene_map <- read_gene_map(file.path(dirname(path), "genes.fpkm_tracking"))
  mat <- label_matrix(mat, gene_map)
  sample_info <- unique(rg[, c("sample", "group")])
  sample_info <- sample_info[match(colnames(mat), sample_info$sample), , drop = FALSE]
  list(mat = mat, sample_info = sample_info, source = basename(path), value_col = value_col)
}

read_gene_map <- function(path) {
  if (!file.exists(path)) return(NULL)
  fp <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (!all(c("tracking_id", "gene_short_name") %in% names(fp))) return(NULL)
  fp
}

label_matrix <- function(mat, gene_map) {
  genes <- rownames(mat)
  nearest <- NULL
  if (!is.null(gene_map)) {
    hit <- match(rownames(mat), gene_map$tracking_id)
    genes <- gene_map$gene_short_name[hit]
    genes[is.na(hit)] <- rownames(mat)[is.na(hit)]
    if ("nearest_ref_id" %in% names(gene_map)) nearest <- gene_map$nearest_ref_id[hit]
  }
  apply_gene_labels(mat, genes, tracking_ids = rownames(mat), nearest_ref = nearest)
}

read_tracking_matrix <- function(path, value_pattern, value_col_name) {
  tr <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  val_cols <- grep(value_pattern, names(tr), value = TRUE, ignore.case = TRUE)
  val_cols <- val_cols[!grepl("variance|conf|status|dispersion|uncertainty", val_cols, ignore.case = TRUE)]
  if (length(val_cols) == 0) return(NULL)
  groups <- vapply(val_cols, classify_pcy, character(1))
  sample_names <- val_cols
  if (all(is.na(groups)) || sum(!is.na(groups)) < 2) {
    info <- read_read_groups_info(file.path(dirname(path), "read_groups.info"))
    cond <- condition_from_info(info)
    q_cols <- grep("^q[0-9]+_", val_cols, ignore.case = TRUE)
    if (!is.null(cond) && length(cond) == length(val_cols)) {
      groups <- vapply(cond, classify_pcy, character(1))
      sample_names <- make.unique(paste(cond, seq_along(cond) - 1L, sep = "_rep"), sep = "_")
    } else if (!is.null(cond) && length(q_cols) == length(cond) && length(q_cols) == length(val_cols)) {
      groups <- vapply(cond, classify_pcy, character(1))
      sample_names <- make.unique(paste(cond, seq_along(cond) - 1L, sep = "_rep"), sep = "_")
    }
  } else {
    sample_names <- make.unique(gsub(value_pattern, "", val_cols, ignore.case = TRUE), sep = "_")
  }
  keep <- !is.na(groups)
  if (sum(keep) < 2) return(NULL)
  mat <- as.matrix(tr[, val_cols[keep], drop = FALSE])
  storage.mode(mat) <- "double"
  mat[!is.finite(mat)] <- 0
  colnames(mat) <- sample_names[keep]
  sample_info <- data.frame(sample = colnames(mat), group = unname(groups[keep]), stringsAsFactors = FALSE)
  tid <- if ("tracking_id" %in% names(tr)) tr$tracking_id else rownames(tr)
  rownames(mat) <- tid
  gene_col <- if ("gene_short_name" %in% names(tr)) tr$gene_short_name else tid
  nearest <- if ("nearest_ref_id" %in% names(tr)) tr$nearest_ref_id else NULL
  mat <- apply_gene_labels(mat, gene_col, tracking_ids = tid, nearest_ref = nearest)
  list(mat = mat, sample_info = sample_info, source = basename(path), value_col = value_col_name)
}

load_expression <- function(project_dir) {
  rg <- file.path(project_dir, "genes.read_group_tracking")
  if (file.exists(rg)) {
    log_msg("读取 ", rg)
    obj <- tryCatch(read_read_group_tracking(rg), error = function(e) {
      log_msg("genes.read_group_tracking 读取失败: ", e$message)
      NULL
    })
    if (!is.null(obj)) return(obj)
  }
  ct <- file.path(project_dir, "genes.count_tracking")
  if (file.exists(ct)) {
    log_msg("读取 ", ct)
    obj <- tryCatch(
      read_tracking_matrix(ct, "_count$|^q[0-9]+_count$", "count"),
      error = function(e) {
        log_msg("genes.count_tracking 读取失败: ", e$message)
        NULL
      }
    )
    if (!is.null(obj)) return(obj)
  }
  fp <- file.path(project_dir, "genes.fpkm_tracking")
  if (file.exists(fp)) {
    log_msg("读取 ", fp)
    obj <- tryCatch(
      read_tracking_matrix(fp, "_FPKM$|^q[0-9]+_FPKM$", "FPKM"),
      error = function(e) {
        log_msg("genes.fpkm_tracking 读取失败: ", e$message)
        NULL
      }
    )
    if (!is.null(obj)) return(obj)
  }
  stop("没有可用的基因表达文件。", call. = FALSE)
}

detect_value_type <- function(obj) {
  if (obj$value_col %in% c("raw_frags", "count")) return("counts")
  x <- as.numeric(obj$mat)
  x <- x[is.finite(x)]
  frac_int <- mean(abs(x - round(x)) < 1e-6)
  if (frac_int > 0.85 && stats::quantile(x, 0.95, na.rm = TRUE) > 50) return("counts")
  "fpkm"
}

prepare_matrix <- function(mat) {
  mat[!is.finite(mat)] <- 0
  n_neg <- sum(mat < 0, na.rm = TRUE)
  if (n_neg > 0) {
    log_msg("负值记为 0: ", n_neg, " 个")
    mat[mat < 0] <- 0
  }
  rs <- rowSums(mat)
  cs <- colSums(mat)
  if (any(cs <= 0)) log_msg("去掉总表达为 0 的样品: ", paste(colnames(mat)[cs <= 0], collapse = ", "))
  mat <- mat[rs > 0, cs > 0, drop = FALSE]
  mat
}

filter_low_expression <- function(mat, sample_info, value_type) {
  min_n <- max(2L, as.integer(min(table(sample_info$group))))
  if (value_type == "counts") {
    keep <- rowSums(mat >= 10, na.rm = TRUE) >= min_n
  } else {
    keep <- rowSums(mat > 1, na.rm = TRUE) >= min_n
  }
  if (sum(keep) < 200) {
    keep <- rowSums(mat > 0, na.rm = TRUE) >= min_n
    log_msg("严格过滤后基因过少，改为至少在 ", min_n, " 个样品中有表达")
  }
  log_msg("低表达过滤: 保留 ", sum(keep), " / ", nrow(mat), " 个基因")
  mat[keep, , drop = FALSE]
}

align_samples <- function(mat, sample_info) {
  sample_info <- sample_info[match(colnames(mat), sample_info$sample), , drop = FALSE]
  if (any(is.na(sample_info$sample))) stop("样品列和分组对不上。", call. = FALSE)
  rownames(sample_info) <- sample_info$sample
  sample_info$group <- factor(sample_info$group, levels = c("NTC", "shP1", "shP2"))
  sample_info
}

# -----------------------------------------------------------------------------
# 差异分析
# -----------------------------------------------------------------------------
residual_df <- function(sample_info) {
  g <- droplevels(sample_info$group)
  ncol_s <- nrow(sample_info)
  nlev <- nlevels(g)
  list(n = ncol_s, n_group = nlev, df = ncol_s - nlev)
}

standardize_de <- function(df) {
  if (!"mean_expr" %in% names(df)) {
    if ("baseMean" %in% names(df)) df$mean_expr <- df$baseMean
    else if ("AveExpr" %in% names(df)) df$mean_expr <- df$AveExpr
    else df$mean_expr <- NA_real_
  }
  df$regulation <- "NS"
  ok <- !is.na(df$pvalue) & !is.na(df$log2FC)
  df$regulation[ok & df$pvalue < p_cutoff & df$log2FC < -lfc_cutoff] <- "Down"
  df$regulation[ok & df$pvalue < p_cutoff & df$log2FC > lfc_cutoff] <- "Up"
  df$gene <- as.character(df$gene)
  df
}

de_export <- function(df, comparison, method) {
  df <- standardize_de(df)
  df$comparison <- comparison
  df$method <- method
  df <- df[order(df$pvalue, df$log2FC), , drop = FALSE]
  keep <- intersect(
    c("gene", "comparison", "log2FC", "pvalue", "padj", "mean_expr", "regulation", "method", "stat", "lfcSE"),
    names(df)
  )
  df[, keep, drop = FALSE]
}

run_deseq <- function(mat, sample_info) {
  counts <- round(mat)
  counts[counts > .Machine$integer.max] <- .Machine$integer.max
  storage.mode(counts) <- "integer"
  group <- droplevels(sample_info$group)
  coldata <- data.frame(group = group, row.names = colnames(counts))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ group)
  dds <- if ("quiet" %in% names(formals(DESeq2::DESeq))) {
    DESeq2::DESeq(dds, quiet = TRUE)
  } else {
    suppressMessages(DESeq2::DESeq(dds))
  }
  sf <- DESeq2::sizeFactors(dds)
  log_msg("DESeq2 size factor: ", paste(paste0(names(sf), "=", signif(sf, 3)), collapse = ", "))
  heat <- tryCatch(
    SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE)),
    error = function(e) {
      log_msg("vst 失败，热图改用 log2(标准化 count + 1): ", e$message)
      log2(DESeq2::counts(dds, normalized = TRUE) + 1)
    }
  )
  list(dds = dds, heat = heat, method = "DESeq2")
}

contrast_from_deseq <- function(dds, treat, control) {
  res <- DESeq2::results(dds, contrast = c("group", treat, control), alpha = p_cutoff)
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df$log2FC <- df$log2FoldChange
  df
}

run_limma <- function(mat, sample_info, treat, control) {
  log_mat <- limma::normalizeBetweenArrays(log2(pmax(mat, 0) + 1), method = "quantile")
  group <- droplevels(sample_info$group)
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  fit <- limma::lmFit(log_mat, design)
  cont <- limma::makeContrasts(contrasts = paste0(treat, "-", control), levels = design)
  cfit <- limma::contrasts.fit(fit, cont)
  fit2 <- tryCatch(
    limma::eBayes(cfit, trend = TRUE, robust = TRUE),
    error = function(e) limma::eBayes(cfit, trend = TRUE)
  )
  tt <- limma::topTable(fit2, number = Inf, sort.by = "none")
  df <- data.frame(
    gene = rownames(tt),
    log2FC = tt$logFC,
    AveExpr = tt$AveExpr,
    pvalue = tt$P.Value,
    padj = tt$adj.P.Val,
    stat = tt$t,
    stringsAsFactors = FALSE
  )
  list(de = df, heat = log_mat, method = "limma_trend_quantile")
}

find_gene_exp_diff <- function(project_dir) {
  candidates <- c(
    "gene_exp.diff", "gene_exp.diff.txt", "gene_exp",
    file.path("cuffdiff", "gene_exp.diff")
  )
  for (f in candidates) {
    p <- file.path(project_dir, f)
    if (file.exists(p) && !dir.exists(p)) return(p)
  }
  hits <- list.files(project_dir, pattern = "^gene_exp\\.diff", full.names = TRUE)
  if (length(hits) >= 1) return(hits[[1]])
  NA_character_
}

read_cuffdiff_contrast <- function(path, treat, control) {
  gd <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  need <- c("sample_1", "sample_2", "p_value")
  if (!all(need %in% names(gd))) {
    stop("gene_exp.diff 缺少 sample_1/sample_2/p_value 列。", call. = FALSE)
  }
  fc_col <- grep("^log2", names(gd), value = TRUE)[1]
  if (is.na(fc_col)) stop("gene_exp.diff 缺少 log2 fold change 列。", call. = FALSE)
  g1 <- vapply(as.character(gd$sample_1), classify_pcy, character(1))
  g2 <- vapply(as.character(gd$sample_2), classify_pcy, character(1))
  forward <- !is.na(g1) & !is.na(g2) & g1 == control & g2 == treat
  reverse <- !is.na(g1) & !is.na(g2) & g1 == treat & g2 == control
  use <- forward | reverse
  if (!any(use)) {
    stop("gene_exp.diff 里没有 ", treat, " 与 ", control, " 的比较。样品名: ",
         paste(unique(c(gd$sample_1, gd$sample_2)), collapse = ", "), call. = FALSE)
  }
  sub <- gd[use, , drop = FALSE]
  forward <- forward[use]
  gene <- if ("gene" %in% names(sub)) as.character(sub$gene) else as.character(sub$test_id)
  gene <- clean_gene_names(gene, tracking_ids = if ("test_id" %in% names(sub)) sub$test_id else gene)
  lfc <- suppressWarnings(as.numeric(sub[[fc_col]]))
  lfc[!forward] <- -lfc[!forward]
  p <- suppressWarnings(as.numeric(sub$p_value))
  padj <- if ("q_value" %in% names(sub)) suppressWarnings(as.numeric(sub$q_value)) else NA_real_
  if ("status" %in% names(sub)) {
    bad <- toupper(as.character(sub$status)) != "OK"
    p[bad] <- NA_real_
    padj[bad] <- NA_real_
  }
  df <- data.frame(gene = gene, log2FC = lfc, pvalue = p, padj = padj, stringsAsFactors = FALSE)
  df <- df[!is.na(df$gene) & nzchar(df$gene), , drop = FALSE]
  if (any(duplicated(df$gene))) {
    log_msg("gene_exp.diff 里同一基因有多行，log2FC 取平均，P 取较大值")
    sp <- split(df, df$gene)
    df <- do.call(rbind, lapply(sp, function(x) {
      data.frame(
        gene = x$gene[1],
        log2FC = mean(x$log2FC, na.rm = TRUE),
        pvalue = if (all(is.na(x$pvalue))) NA_real_ else max(x$pvalue, na.rm = TRUE),
        padj = if (all(is.na(x$padj))) NA_real_ else max(x$padj, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
    rownames(df) <- NULL
  }
  df
}

# -----------------------------------------------------------------------------
# 火山图、表达热图、气泡图、富集热图
# -----------------------------------------------------------------------------
plot_volcano <- function(df, title, outfile) {
  plot_df <- df[!is.na(df$log2FC) & !is.na(df$pvalue) & df$pvalue > 0, , drop = FALSE]
  if (nrow(plot_df) < 2) {
    log_msg("火山图跳过，可用基因少于 2 个: ", title)
    return(invisible(NULL))
  }
  plot_df$y <- -log10(pmax(plot_df$pvalue, 1e-300))
  plot_df$regulation <- factor(plot_df$regulation, levels = c("Down", "Up", "NS"))
  label_of <- function(flag, n) {
    sub <- plot_df[plot_df$regulation == flag, , drop = FALSE]
    if (nrow(sub) == 0) return(character())
    sub <- sub[order(sub$pvalue, -abs(sub$log2FC)), , drop = FALSE]
    utils::head(sub$gene, n)
  }
  labs <- c(label_of("Down", 10), label_of("Up", 5))
  plot_df$label <- ifelse(plot_df$gene %in% labs, plot_df$gene, NA_character_)
  n_down <- sum(plot_df$regulation == "Down")
  n_up <- sum(plot_df$regulation == "Up")
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = log2FC, y = y, color = regulation)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.3) +
    ggplot2::scale_color_manual(
      values = c(Down = "#2166AC", Up = "#B2182B", NS = "grey75"),
      labels = c(Down = "下调", Up = "上调", NS = "不显著"),
      drop = FALSE
    ) +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = 2, color = "grey40") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("下调 P<%s: %d；上调 P<%s: %d", p_cutoff, n_down, p_cutoff, n_up),
      x = "log2 Fold Change (sh / NTC)",
      y = "-log10(P)",
      color = NULL
    )
  if (any(!is.na(plot_df$label))) {
    p <- p + ggrepel::geom_text_repel(
      data = plot_df[!is.na(plot_df$label), , drop = FALSE],
      ggplot2::aes(label = label),
      size = 3, max.overlaps = 40, show.legend = FALSE
    )
  }
  save_gg(p, outfile, width = 8, height = 6.5)
}

plot_expression_heatmap <- function(heat, sample_info, genes, title, outfile) {
  genes <- unique(genes[genes %in% rownames(heat)])
  if (length(genes) > 80) {
    log_msg("表达热图只画 P 值最小的 80 个下调基因（共 ", length(genes), " 个）")
    genes <- genes[seq_len(80)]
  }
  if (length(genes) < 2) {
    log_msg("表达热图跳过，下调基因少于 2 个")
    note_empty(paste0(outfile, "_EMPTY.txt"), "downregulated genes < 2")
    return(invisible(NULL))
  }
  sub <- heat[genes, , drop = FALSE]
  sds <- apply(sub, 1, stats::sd, na.rm = TRUE)
  sub <- sub[is.finite(sds) & sds > 0, , drop = FALSE]
  if (nrow(sub) < 2) {
    log_msg("表达热图跳过，基因方差为 0")
    return(invisible(NULL))
  }
  ann <- data.frame(Group = as.character(sample_info$group), row.names = sample_info$sample)
  ann <- ann[colnames(sub), , drop = FALSE]
  pal <- c(NTC = "#4C78A8", shP1 = "#E45756", shP2 = "#54A24B")
  ann_col <- list(Group = pal[names(pal) %in% unique(ann$Group)])
  draw <- function(dist_rows) {
    pheatmap::pheatmap(
      sub, scale = "row", annotation_col = ann, annotation_colors = ann_col,
      show_rownames = nrow(sub) <= 80, fontsize_row = 6, main = title,
      color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      clustering_distance_rows = dist_rows,
      clustering_distance_cols = "euclidean",
      border_color = NA
    )
  }
  height <- max(6, min(16, 0.16 * nrow(sub) + 3))
  draw_safe <- function() {
    tryCatch(draw("correlation"), error = function(e) draw("euclidean"))
  }
  grDevices::pdf(paste0(outfile, ".pdf"), width = 8, height = height)
  on.exit(while (grDevices::dev.cur() > 1) grDevices::dev.off(), add = TRUE)
  draw_safe()
  grDevices::dev.off()
  grDevices::png(
    paste0(outfile, ".png"),
    width = 8, height = height, units = "in", res = 300
  )
  draw_safe()
  grDevices::dev.off()
}

parse_ratio <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) {
    if (length(z) != 2) return(NA_real_)
    a <- suppressWarnings(as.numeric(z[1]))
    b <- suppressWarnings(as.numeric(z[2]))
    if (!is.finite(a) || !is.finite(b) || b == 0) NA_real_ else a / b
  }, numeric(1))
}

plot_bubble <- function(df, title, outfile, top_n = 15) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  df <- df[order(df$pvalue, df$p.adjust), , drop = FALSE]
  df <- utils::head(df, top_n)
  df$GeneRatioNum <- parse_ratio(df$GeneRatio)
  desc <- as.character(df$Description)
  desc[nchar(desc) > 55] <- paste0(substr(desc[nchar(desc) > 55], 1, 52), "...")
  desc <- make.unique(desc, sep = " ")
  df$Description <- factor(desc, levels = rev(desc))
  color_val <- df$p.adjust
  color_val[!is.finite(color_val)] <- df$pvalue[!is.finite(color_val)]
  df$color_val <- color_val
  p <- ggplot2::ggplot(df, ggplot2::aes(x = GeneRatioNum, y = Description, size = Count)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title, x = "Gene Ratio", y = NULL, size = "Count", color = "p.adjust")
  if (length(unique(df$color_val[is.finite(df$color_val)])) >= 2) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = color_val), alpha = 0.9) +
      ggplot2::scale_color_gradient(low = "#B2182B", high = "#2166AC")
  } else {
    p <- p + ggplot2::geom_point(color = "#B2182B", alpha = 0.9)
  }
  save_gg(p, outfile, width = 9, height = max(4.5, 0.38 * nrow(df) + 1.8))
}

plot_enrich_heatmap <- function(df, fc_symbol, fc_entrez, title, outfile) {
  if (is.null(df) || nrow(df) == 0 || !"geneID" %in% names(df)) return(invisible(NULL))
  df <- df[order(df$pvalue, df$p.adjust), , drop = FALSE]
  df <- utils::head(df, 12)
  tokens <- unlist(strsplit(as.character(df$geneID), "/"))
  use_symbol <- mean(tokens %in% names(fc_symbol), na.rm = TRUE) >= mean(tokens %in% names(fc_entrez), na.rm = TRUE)
  fc <- if (isTRUE(use_symbol)) fc_symbol else fc_entrez
  term_genes <- lapply(strsplit(as.character(df$geneID), "/"), function(gs) intersect(gs, names(fc)))
  genes <- unique(unlist(term_genes))
  if (length(genes) < 2 || nrow(df) < 1) {
    log_msg("富集热图跳过，可画的基因或条目不足: ", title)
    return(invisible(NULL))
  }
  if (length(genes) > 40) {
    score <- table(unlist(term_genes))
    genes <- names(sort(score[genes], decreasing = TRUE))[seq_len(40)]
    term_genes <- lapply(term_genes, function(gs) intersect(gs, genes))
  }
  desc <- as.character(df$Description)
  desc[nchar(desc) > 40] <- paste0(substr(desc[nchar(desc) > 40], 1, 37), "...")
  desc <- make.unique(desc, sep = " ")
  mat <- matrix(NA_real_, nrow = length(genes), ncol = length(desc), dimnames = list(genes, desc))
  for (j in seq_along(desc)) {
    gs <- term_genes[[j]]
    mat[gs, j] <- fc[gs]
  }
  if (sum(is.finite(mat)) < 2) return(invisible(NULL))
  gene_order <- names(sort(rowMeans(mat, na.rm = TRUE)))
  mat <- mat[gene_order, , drop = FALSE]
  height <- max(6, min(16, 0.22 * nrow(mat) + 3))
  width <- max(7, min(14, 0.45 * ncol(mat) + 4))
  draw <- function() {
    pheatmap::pheatmap(
      mat,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      na_col = "grey92",
      color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      main = title,
      fontsize_row = 7,
      fontsize_col = 8,
      border_color = NA
    )
  }
  grDevices::pdf(paste0(outfile, ".pdf"), width = width, height = height)
  on.exit(while (grDevices::dev.cur() > 1) grDevices::dev.off(), add = TRUE)
  tryCatch(draw(), error = function(e) log_msg("富集热图失败: ", e$message))
  grDevices::dev.off()
  grDevices::png(paste0(outfile, ".png"), width = width, height = height, units = "in", res = 300)
  tryCatch(draw(), error = function(e) log_msg("富集热图失败: ", e$message))
  grDevices::dev.off()
}

# -----------------------------------------------------------------------------
# GO / KEGG / Reactome
# -----------------------------------------------------------------------------
map_to_entrez <- function(symbols) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  symbols <- symbols[!grepl("^(XLOC|TCONS|CUFF)_", symbols, ignore.case = TRUE)]
  empty <- data.frame(gene = character(), entrez = character(), stringsAsFactors = FALSE)
  if (length(symbols) == 0) return(empty)
  m <- tryCatch(
    clusterProfiler::bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) {
      log_msg("SYMBOL 转 ENTREZ 失败: ", e$message)
      empty
    }
  )
  if (nrow(m) == 0) return(empty)
  m <- m[!duplicated(m$SYMBOL), , drop = FALSE]
  data.frame(gene = m$SYMBOL, entrez = as.character(m$ENTREZID), stringsAsFactors = FALSE)
}

enrich_pair <- function(strict_fun, relax_fun, label) {
  got <- tryCatch(strict_fun(), error = function(e) {
    log_msg(label, " 失败: ", e$message)
    NULL
  })
  relaxed <- FALSE
  n0 <- if (is.null(got)) 0L else nrow(as.data.frame(got))
  if (n0 == 0) {
    got <- tryCatch(relax_fun(), error = function(e) {
      log_msg(label, " 放宽阈值后仍失败: ", e$message)
      NULL
    })
    relaxed <- TRUE
  }
  list(obj = got, relaxed = relaxed)
}

as_enrich_df <- function(obj, relaxed) {
  if (is.null(obj)) return(NULL)
  df <- as.data.frame(obj)
  if (nrow(df) == 0) return(NULL)
  df <- df[order(df$pvalue, df$p.adjust), , drop = FALSE]
  if (relaxed) df <- utils::head(df, 50)
  df
}

save_enrichment <- function(obj, relaxed, out_stub, title, fc_symbol, fc_entrez) {
  df <- as_enrich_df(obj, relaxed)
  if (is.null(df)) {
    note_empty(paste0(out_stub, "_EMPTY.txt"), "no enrichment terms")
    log_msg("没有富集条目: ", title)
    return(invisible(NULL))
  }
  write_csv(df, paste0(out_stub, ".csv"))
  plot_title <- title
  if (relaxed) {
    plot_title <- paste0(title, "\n（没有 P<0.05 的条目，图中是按 P 值靠前的条目）")
    log_msg("富集没有 P<0.05 的条目，气泡图改为展示排序靠前的条目: ", title)
  }
  tryCatch(plot_bubble(df, plot_title, paste0(out_stub, "_bubble")), error = function(e) {
    log_msg("气泡图失败: ", e$message)
  })
  tryCatch(
    plot_enrich_heatmap(df, fc_symbol, fc_entrez, plot_title, paste0(out_stub, "_heatmap")),
    error = function(e) log_msg("富集热图失败: ", e$message)
  )
}

call_with_universe <- function(fun, args, universe) {
  if ("universe" %in% names(formals(fun))) args$universe <- universe
  do.call(fun, args)
}

run_downstream <- function(de, heat, sample_info, comparison, outdir, id_map) {
  down <- de[de$regulation == "Down", , drop = FALSE]
  write_csv(down, file.path(outdir, "down_P0.05.csv"))
  log_msg(comparison, " 下调且 P<", p_cutoff, ": ", nrow(down), " 个基因")
  tryCatch(
    plot_expression_heatmap(
      heat, sample_info, down$gene,
      paste0(comparison, " 下调基因 (P<", p_cutoff, ")"),
      file.path(outdir, "down_P0.05_expression_heatmap")
    ),
    error = function(e) log_msg("表达热图失败: ", e$message)
  )
  if (nrow(down) < 3) {
    note_empty(file.path(outdir, "enrichment_EMPTY.txt"), "downregulated genes < 3; skip GO/KEGG/pathway")
    log_msg(comparison, " 下调基因少于 3 个，跳过 GO/KEGG/通路")
    return(invisible(NULL))
  }

  all_genes <- de$gene[!grepl("^(XLOC|TCONS|CUFF)_", de$gene)]
  mp_all <- id_map[id_map$gene %in% all_genes, , drop = FALSE]
  mp_down <- id_map[id_map$gene %in% down$gene, , drop = FALSE]
  universe <- unique(mp_all$entrez)
  entrez <- unique(mp_down$entrez)
  entrez <- intersect(entrez, universe)
  log_msg(comparison, " 下调基因映射到 Entrez: ", length(entrez), " / ", nrow(down))
  if (length(entrez) < 3) {
    note_empty(file.path(outdir, "enrichment_EMPTY.txt"), "fewer than 3 genes mapped to Entrez")
    return(invisible(NULL))
  }

  fc_symbol <- setNames(down$log2FC, down$gene)
  fc_entrez <- setNames(down$log2FC[match(mp_down$gene, down$gene)], mp_down$entrez)
  fc_entrez <- fc_entrez[!duplicated(names(fc_entrez)) & nzchar(names(fc_entrez))]

  go_dir <- file.path(outdir, "GO")
  kegg_dir <- file.path(outdir, "KEGG")
  pw_dir <- file.path(outdir, "Pathway")
  dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(kegg_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pw_dir, recursive = TRUE, showWarnings = FALSE)

  for (ont in c("BP", "CC", "MF")) {
    ont_name <- c(BP = "Biological Process", CC = "Cellular Component", MF = "Molecular Function")[[ont]]
    res <- enrich_pair(
      function() clusterProfiler::enrichGO(
        gene = entrez, universe = universe, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
        ont = ont, pAdjustMethod = "BH", pvalueCutoff = p_cutoff, qvalueCutoff = 1,
        minGSSize = 5, maxGSSize = 500, readable = TRUE
      ),
      function() clusterProfiler::enrichGO(
        gene = entrez, universe = universe, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
        ont = ont, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
        minGSSize = 5, maxGSSize = 500, readable = TRUE
      ),
      paste("GO", ont)
    )
    save_enrichment(
      res$obj, res$relaxed,
      file.path(go_dir, paste0("GO_", ont)),
      paste0(comparison, " 下调基因 GO ", ont_name),
      fc_symbol, fc_entrez
    )
  }

  res_k <- enrich_pair(
    function() call_with_universe(
      clusterProfiler::enrichKEGG,
      list(
        gene = entrez, organism = "hsa", pvalueCutoff = p_cutoff, qvalueCutoff = 1,
        minGSSize = 5, maxGSSize = 500
      ),
      universe
    ),
    function() call_with_universe(
      clusterProfiler::enrichKEGG,
      list(
        gene = entrez, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1,
        minGSSize = 5, maxGSSize = 500
      ),
      universe
    ),
    "KEGG"
  )
  if (!is.null(res_k$obj) && nrow(as.data.frame(res_k$obj)) > 0) {
    res_k$obj <- tryCatch(
      clusterProfiler::setReadable(res_k$obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
      error = function(e) res_k$obj
    )
  }
  save_enrichment(
    res_k$obj, res_k$relaxed,
    file.path(kegg_dir, "KEGG"),
    paste0(comparison, " 下调基因 KEGG"),
    fc_symbol, fc_entrez
  )

  if (!has_pkg("ReactomePA")) {
    note_empty(file.path(pw_dir, "Reactome_EMPTY.txt"), "ReactomePA is not installed")
    log_msg("未安装 ReactomePA，跳过通路富集")
    return(invisible(NULL))
  }
  res_p <- enrich_pair(
    function() call_with_universe(
      ReactomePA::enrichPathway,
      list(
        gene = entrez, organism = "human", pvalueCutoff = p_cutoff, qvalueCutoff = 1,
        minGSSize = 5, maxGSSize = 500, readable = TRUE
      ),
      universe
    ),
    function() call_with_universe(
      ReactomePA::enrichPathway,
      list(
        gene = entrez, organism = "human", pvalueCutoff = 1, qvalueCutoff = 1,
        minGSSize = 5, maxGSSize = 500, readable = TRUE
      ),
      universe
    ),
    "Reactome"
  )
  save_enrichment(
    res_p$obj, res_p$relaxed,
    file.path(pw_dir, "Reactome"),
    paste0(comparison, " 下调基因 Reactome 通路"),
    fc_symbol, fc_entrez
  )
}

write_readme <- function(result_dir) {
  txt <- c(
    "PCY_RNA 结果说明",
    "log2FC = log2(sh / NTC)。负值表示该 sh 低于 NTC（下调）。",
    "下调基因：P < 0.05 且 log2FC < 0。火山图里的横线是 P = 0.05。",
    "同一组的重复样品作为生物学重复，shP1、shP2、NTC 各算一组。",
    "",
    "shP1_vs_NTC/ 与 shP2_vs_NTC/ 结构相同：",
    "  DEG_all.csv                         全部基因",
    "  volcano.pdf / .png                  火山图",
    "  down_P0.05.csv                      下调且 P<0.05",
    "  down_P0.05_expression_heatmap       这些基因的表达热图",
    "  GO/GO_BP_bubble、GO_CC_bubble、GO_MF_bubble    GO 气泡图",
    "  GO/GO_*_heatmap                     GO 富集热图（颜色是 log2FC）",
    "  KEGG/KEGG_bubble、KEGG_heatmap      KEGG",
    "  Pathway/Reactome_bubble、Reactome_heatmap   通路（Reactome）",
    "",
    "物种注释是人 org.Hs.eg.db。没有官方基因符号的 XLOC 会留在差异表，但通常进不了 GO/KEGG。"
  )
  writeLines(txt, file.path(result_dir, "00_README.txt"))
}

# -----------------------------------------------------------------------------
# 主程序
# -----------------------------------------------------------------------------
main <- function() {
  project_dir <- resolve_project_dir()
  result_dir <- file.path(project_dir, "PCY_RNA_results")
  log_dir <- file.path(result_dir, "00_logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <<- file.path(log_dir, paste0("PCY_RNA_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
  log_msg("项目目录: ", project_dir)
  log_msg("P 阈值: ", p_cutoff, "；下调还要求 log2FC < ", -lfc_cutoff)

  expr <- load_expression(project_dir)
  mat <- prepare_matrix(expr$mat)
  sample_info <- expr$sample_info
  sample_info <- sample_info[sample_info$sample %in% colnames(mat), , drop = FALSE]
  sample_info <- align_samples(mat, sample_info)
  log_msg("样品: ", paste(paste0(sample_info$sample, "=", sample_info$group), collapse = "; "))
  counts <- table(sample_info$group)
  log_msg("每组样品数: ", paste(paste0(names(counts), "=", as.integer(counts)), collapse = ", "))
  miss <- setdiff(c("NTC", "shP1", "shP2"), as.character(sample_info$group))
  if (length(miss) > 0) {
    stop("缺少组别 ", paste(miss, collapse = ", "),
         "。当前样品: ", paste(paste0(sample_info$sample, "=", sample_info$group), collapse = "; "),
         call. = FALSE)
  }
  write_csv(sample_info, file.path(result_dir, "00_sample_info.csv"))
  if (identical(Sys.getenv("PCY_RNA_STOP_AFTER", unset = ""), "load")) {
    log_msg("PCY_RNA_STOP_AFTER=load，读入检查结束")
    return(invisible(list(mat = mat, sample_info = sample_info, expr = expr)))
  }

  install_if_missing(c("ggplot2", "ggrepel", "pheatmap"), bioc = FALSE, required = TRUE)
  install_if_missing(c("DESeq2", "limma", "clusterProfiler", "org.Hs.eg.db"), bioc = TRUE, required = TRUE)
  install_if_missing(c("ReactomePA"), bioc = TRUE, required = FALSE)
  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
    library(pheatmap)
    library(DESeq2)
    library(limma)
    library(clusterProfiler)
    library(org.Hs.eg.db)
  })
  if (has_pkg("ReactomePA")) suppressPackageStartupMessages(library(ReactomePA))

  value_type <- detect_value_type(expr)
  log_msg("数值类型: ", value_type, " （来自 ", expr$source, " / ", expr$value_col, "）")
  mat <- filter_low_expression(mat, sample_info, value_type)
  sample_info <- align_samples(mat, sample_info)
  df_info <- residual_df(sample_info)
  log_msg("差异模型自由度: 样品 ", df_info$n, "，组别 ", df_info$n_group, "，残差 df ", df_info$df)

  comparisons <- list(
    list(name = "shP1_vs_NTC", treat = "shP1", control = "NTC", title = "shP1 vs NTC"),
    list(name = "shP2_vs_NTC", treat = "shP2", control = "NTC", title = "shP2 vs NTC")
  )
  method <- NULL
  heat <- NULL
  de_sets <- list()

  if (df_info$df > 0 && value_type == "counts") {
    fit <- run_deseq(mat, sample_info)
    heat <- fit$heat
    method <- fit$method
    for (comp in comparisons) {
      de_sets[[comp$name]] <- contrast_from_deseq(fit$dds, comp$treat, comp$control)
    }
  } else if (df_info$df > 0) {
    heat <- NULL
    method <- "limma_trend_quantile"
    for (comp in comparisons) {
      fit <- run_limma(mat, sample_info, comp$treat, comp$control)
      de_sets[[comp$name]] <- fit$de
      heat <- fit$heat
    }
  } else {
    diff_path <- find_gene_exp_diff(project_dir)
    if (is.na(diff_path)) {
      stop(
        "每个组都只有 1 个样品，无法从重复里估计 P 值。",
        "请在同一目录提供 Cuffdiff 的 gene_exp.diff，或保证 shP1、shP2、NTC 至少有一组包含生物学重复。",
        call. = FALSE
      )
    }
    log_msg("没有残差自由度，改用 Cuffdiff 的 P 值: ", diff_path)
    method <- "Cuffdiff_gene_exp.diff"
    for (comp in comparisons) {
      de_sets[[comp$name]] <- read_cuffdiff_contrast(diff_path, comp$treat, comp$control)
    }
    if (value_type == "counts") {
      coldata <- data.frame(group = droplevels(sample_info$group), row.names = colnames(mat))
      dds <- DESeq2::DESeqDataSetFromMatrix(round(mat), coldata, design = ~ 1)
      dds <- DESeq2::estimateSizeFactors(dds)
      heat <- log2(DESeq2::counts(dds, normalized = TRUE) + 1)
    } else {
      heat <- limma::normalizeBetweenArrays(log2(pmax(mat, 0) + 1), method = "quantile")
    }
  }

  write_readme(result_dir)
  id_map <- map_to_entrez(rownames(mat))
  log_msg("全基因映射到 Entrez: ", nrow(id_map), " / ", nrow(mat))

  summary_rows <- list()
  for (comp in comparisons) {
    log_msg("开始 ", comp$title)
    outdir <- file.path(result_dir, comp$name)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    de <- tryCatch(
      de_export(de_sets[[comp$name]], comp$name, method),
      error = function(e) {
        log_msg(comp$name, " 差异表失败: ", e$message)
        NULL
      }
    )
    if (is.null(de)) next
    write_csv(de, file.path(outdir, "DEG_all.csv"))
    tryCatch(
      plot_volcano(de, comp$title, file.path(outdir, "volcano")),
      error = function(e) log_msg("火山图失败: ", e$message)
    )
    tryCatch(
      run_downstream(de, heat, sample_info, comp$title, outdir, id_map),
      error = function(e) log_msg(comp$name, " 下游分析失败: ", e$message)
    )
    summary_rows[[comp$name]] <- data.frame(
      comparison = comp$name,
      method = method,
      genes = nrow(de),
      up_p = sum(de$regulation == "Up"),
      down_p = sum(de$regulation == "Down"),
      stringsAsFactors = FALSE
    )
  }
  if (length(summary_rows) > 0) {
    write_csv(do.call(rbind, summary_rows), file.path(result_dir, "00_summary.csv"))
  }
  log_msg("完成。结果在: ", result_dir)
  invisible(result_dir)
}

main()
