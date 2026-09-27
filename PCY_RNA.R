#!/usr/bin/env Rscript
# =============================================================================
# PCY_RNA
# 四个样品单独保留：NTC_rep0、NTC_rep1、PCY_sh1、PCY_sh4
#
# 六类比较（各自出表、出图；两个 NTC 不在 1-vs-1 里合并）：
#   1. PCY_sh1 vs NTC_rep0、PCY_sh4 vs NTC_rep0、
#      PCY_sh1 vs NTC_rep1、PCY_sh4 vs NTC_rep1
#   2. mean(PCY_sh1, PCY_sh4) vs mean(NTC_rep0, NTC_rep1)
#   3. 相对 NTC_rep0，sh1 与 sh4 的共同上调
#   4. 相对 NTC_rep1，sh1 与 sh4 的共同上调
#   5. mean(PCY_sh1, PCY_sh4) vs NTC_rep0
#   6. mean(PCY_sh1, PCY_sh4) vs NTC_rep1
#
# 只分析上调。有 P 值时先保留 P < 0.05，再按 FC > 1 和 FC > 1.25 分层。
# 每一层做 GO、KEGG、Reactome 通路，并画气泡图、富集热图和表达热图。
# 1-vs-1 没有重复时不编造 P 值；若同目录有 gene_exp.diff，则用其中的 P 值。
#
# 囊泡运输专项单独出结果，不改全基因组 GO 的 p 值：
#   内体运输与内体运动、囊泡运输、外泌体分泌，以及名称上属于这三类的 GO。
#   每个比较的 FoldChange 下有 Focused_vesicle_transport/。
#   六类比较放在一起的热图和相关性在 00_vesicle_transport_across_6/。
#
# 运行：
#   source("E:/R/PCY_RNA/PCY_RNA.R", encoding = "UTF-8")
# 物种：人（org.Hs.eg.db）。log2FC = log2(PCY / NTC)，正值表示 PCY 高于 NTC。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")
options(clusterProfiler.download.method = "auto")

p_cutoff <- 0.05
fc_levels <- c(FC_1 = 1, FC_1.25 = 1.25)
sample_levels <- c("NTC_rep0", "NTC_rep1", "PCY_sh1", "PCY_sh4")

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
# 样品名 → NTC_rep0 / NTC_rep1 / PCY_sh1 / PCY_sh4
# -----------------------------------------------------------------------------
canonical_one <- function(condition, replicate = NA_character_) {
  cond <- toupper(trimws(as.character(condition)))
  rep <- toupper(trimws(as.character(replicate)))
  if (is.na(cond)) cond <- ""
  if (is.na(rep) || !nzchar(rep) || rep == "NA") rep <- ""
  cond_n <- gsub("[^A-Z0-9]", "", cond)
  blob <- gsub("[^A-Z0-9]", "", paste0(cond, rep))
  if (!nzchar(cond_n)) return(NA_character_)
  # 只认 PCY_sh1 / PCY_sh4，不把 TG_sh1、TG_sh5 改名进来
  if (grepl("PCY", cond_n) && grepl("SH4", cond_n)) return("PCY_sh4")
  if (grepl("PCY", cond_n) && grepl("SH1", cond_n)) return("PCY_sh1")
  is_ntc <- grepl("NTC", cond_n) || grepl("CTRL|CONTROL|SHNC", cond_n) || grepl("^NC[0-9]*$", cond_n)
  if (!is_ntc) return(NA_character_)
  if (grepl("REP0", blob) || rep %in% c("0", "REP0")) return("NTC_rep0")
  if (grepl("REP1", blob) || rep %in% c("1", "REP1")) return("NTC_rep1")
  if (grepl("REP0", cond_n)) return("NTC_rep0")
  if (grepl("REP1", cond_n)) return("NTC_rep1")
  NA_character_
}

canonical_sample <- function(condition, replicate = NA_character_) {
  cond <- as.character(condition)
  rep <- as.character(replicate)
  if (length(rep) == 1L && length(cond) > 1L) rep <- rep(rep, length(cond))
  if (length(rep) != length(cond)) rep <- rep(NA_character_, length(cond))
  vapply(seq_along(cond), function(i) canonical_one(cond[[i]], rep[[i]]), character(1))
}

classify_pcy <- function(name) {
  canonical_one(name, NA_character_)
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
  rg$sample <- canonical_sample(rg$condition, rg$replicate)
  dropped <- unique(as.character(rg$condition[is.na(rg$sample)]))
  if (length(dropped) > 0) {
    log_msg("未识别的 condition: ", paste(dropped, collapse = ", "))
  }
  mapped <- unique(data.frame(
    condition = as.character(rg$condition),
    replicate = as.character(rg$replicate),
    sample = rg$sample,
    stringsAsFactors = FALSE
  ))
  mapped <- mapped[!is.na(mapped$sample), , drop = FALSE]
  if (nrow(mapped) > 0) {
    log_msg("组别对应: ", paste(paste0(mapped$condition, " / rep ", mapped$replicate, " -> ", mapped$sample), collapse = "; "))
  }
  rg <- rg[!is.na(rg$sample), , drop = FALSE]
  if (nrow(rg) == 0) return(NULL)
  rg$group <- rg$sample
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
  bare <- gsub(value_pattern, "", val_cols, ignore.case = TRUE)
  groups <- canonical_sample(bare)
  if (sum(!is.na(groups)) < 2) groups <- canonical_sample(val_cols)
  sample_names <- groups
  if (all(is.na(groups)) || sum(!is.na(groups)) < 2) {
    info <- read_read_groups_info(file.path(dirname(path), "read_groups.info"))
    cond <- condition_from_info(info)
    if (!is.null(cond) && length(cond) == length(val_cols)) {
      groups <- canonical_sample(cond, seq_along(cond) - 1L)
      sample_names <- groups
    }
  }
  keep <- !is.na(groups)
  if (sum(keep) < 2) return(NULL)
  mat <- as.matrix(tr[, val_cols[keep], drop = FALSE])
  storage.mode(mat) <- "double"
  mat[!is.finite(mat)] <- 0
  colnames(mat) <- sample_names[keep]
  mat <- collapse_duplicate_samples(mat)
  sample_info <- data.frame(sample = colnames(mat), group = colnames(mat), stringsAsFactors = FALSE)
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

collapse_duplicate_samples <- function(mat) {
  ids <- colnames(mat)
  if (!any(duplicated(ids))) return(mat)
  log_msg("同一组别有多列，取平均: ", paste(unique(ids[duplicated(ids)]), collapse = ", "))
  out_ids <- unique(ids)
  out <- matrix(NA_real_, nrow(mat), length(out_ids), dimnames = list(rownames(mat), out_ids))
  for (id in out_ids) out[, id] <- rowMeans(mat[, ids == id, drop = FALSE], na.rm = TRUE)
  out
}

align_samples <- function(mat, sample_info) {
  sample_info <- sample_info[match(colnames(mat), sample_info$sample), , drop = FALSE]
  if (any(is.na(sample_info$sample))) stop("样品列和分组对不上。", call. = FALSE)
  rownames(sample_info) <- sample_info$sample
  sample_info$group <- factor(as.character(sample_info$group), levels = sample_levels)
  sample_info
}

sample_column <- function(sample_info, sample_id) {
  hit <- sample_info$sample[as.character(sample_info$group) == sample_id]
  if (length(hit) != 1) return(NA_character_)
  hit
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

has_real_p <- function(p) {
  any(!is.na(p))
}

passes_p <- function(p) {
  if (!has_real_p(p)) return(rep(TRUE, length(p)))
  !is.na(p) & p < p_cutoff
}

passes_fc <- function(log2fc, fc) {
  !is.na(log2fc) & log2fc > log2(fc)
}

standardize_de <- function(df) {
  if (!"mean_expr" %in% names(df)) {
    if ("baseMean" %in% names(df)) df$mean_expr <- df$baseMean
    else if ("AveExpr" %in% names(df)) df$mean_expr <- df$AveExpr
    else df$mean_expr <- NA_real_
  }
  df$gene <- as.character(df$gene)
  up <- !is.na(df$log2FC) & df$log2FC > 0 & passes_p(df$pvalue)
  if (all(c("log2FC_sh1", "log2FC_sh4") %in% names(df))) {
    up <- up & df$log2FC_sh1 > 0 & df$log2FC_sh4 > 0
    if (has_real_p(df$pvalue_sh1)) up <- up & passes_p(df$pvalue_sh1)
    if (has_real_p(df$pvalue_sh4)) up <- up & passes_p(df$pvalue_sh4)
  }
  df$regulation <- ifelse(up, "Up", "NS")
  df$regulation[!is.na(df$log2FC) & df$log2FC < 0] <- "Down"
  df
}

de_export <- function(df, comparison, method) {
  df <- standardize_de(df)
  df$comparison <- comparison
  df$method <- method
  df <- df[order(df$pvalue, -df$log2FC), , drop = FALSE]
  keep <- intersect(
    c(
      "gene", "comparison", "log2FC", "log2FC_sh1", "log2FC_sh4",
      "pvalue", "pvalue_sh1", "pvalue_sh4", "padj", "mean_expr",
      "regulation", "method", "stat", "lfcSE"
    ),
    names(df)
  )
  df[, keep, drop = FALSE]
}

normalize_expression <- function(mat, value_type) {
  if (value_type == "counts") {
    counts <- round(pmax(mat, 0))
    counts[counts > .Machine$integer.max] <- .Machine$integer.max
    storage.mode(counts) <- "integer"
    coldata <- data.frame(sample = colnames(counts), row.names = colnames(counts))
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ 1)
    dds <- DESeq2::estimateSizeFactors(dds)
    sf <- DESeq2::sizeFactors(dds)
    log_msg("DESeq2 size factor: ", paste(paste0(names(sf), "=", signif(sf, 3)), collapse = ", "))
    log_mat <- log2(DESeq2::counts(dds, normalized = TRUE) + 1)
    heat <- tryCatch(
      SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE)),
      error = function(e) {
        log_msg("vst 失败，热图改用 log2(标准化 count + 1): ", e$message)
        log_mat
      }
    )
    return(list(log_mat = log_mat, heat = heat, method = "DESeq2_size_factor"))
  }
  log_msg("分位数标准化 log2(FPKM+1)")
  log_mat <- limma::normalizeBetweenArrays(log2(pmax(mat, 0) + 1), method = "quantile")
  list(log_mat = log_mat, heat = log_mat, method = "quantile_logFPKM")
}

limma_pvalues <- function(log_mat, treat_cols, control_cols) {
  cols <- unique(c(control_cols, treat_cols))
  cols <- cols[!is.na(cols) & cols %in% colnames(log_mat)]
  if (length(cols) <= 2) return(NULL)
  sub <- log_mat[, cols, drop = FALSE]
  grp <- ifelse(colnames(sub) %in% treat_cols, "T", "C")
  if (length(unique(grp)) < 2) return(NULL)
  design <- stats::model.matrix(~ 0 + factor(grp, levels = c("C", "T")))
  colnames(design) <- c("C", "T")
  fit <- limma::lmFit(sub, design)
  cont <- limma::makeContrasts(contrasts = "T-C", levels = design)
  cfit <- limma::contrasts.fit(fit, cont)
  fit2 <- tryCatch(
    limma::eBayes(cfit, trend = TRUE, robust = TRUE),
    error = function(e) limma::eBayes(cfit, trend = TRUE)
  )
  tt <- limma::topTable(fit2, number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), pvalue = tt$P.Value, padj = tt$adj.P.Val, stringsAsFactors = FALSE)
}

fc_from_columns <- function(log_mat, treat_cols, control_cols) {
  treat_mean <- rowMeans(log_mat[, treat_cols, drop = FALSE])
  control_mean <- rowMeans(log_mat[, control_cols, drop = FALSE])
  data.frame(
    gene = rownames(log_mat),
    log2FC = as.numeric(treat_mean - control_mean),
    mean_expr = as.numeric((treat_mean + control_mean) / 2),
    pvalue = NA_real_,
    padj = NA_real_,
    stringsAsFactors = FALSE
  )
}

attach_limma_p <- function(de, log_mat, treat_cols, control_cols) {
  got <- tryCatch(
    limma_pvalues(log_mat, treat_cols, control_cols),
    error = function(e) {
      log_msg("limma P 值失败: ", e$message)
      NULL
    }
  )
  if (is.null(got)) return(de)
  hit <- match(de$gene, got$gene)
  de$pvalue <- got$pvalue[hit]
  de$padj <- got$padj[hit]
  de
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
    log_msg("gene_exp.diff 没有 ", treat, " vs ", control, "。样品名: ",
            paste(unique(c(gd$sample_1, gd$sample_2)), collapse = ", "))
    return(NULL)
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

attach_cuffdiff_p <- function(de, project_dir, treat, control) {
  path <- find_gene_exp_diff(project_dir)
  if (is.na(path)) return(de)
  got <- tryCatch(read_cuffdiff_contrast(path, treat, control), error = function(e) {
    log_msg("读取 gene_exp.diff 失败: ", e$message)
    NULL
  })
  if (is.null(got) || !has_real_p(got$pvalue)) return(de)
  hit <- match(de$gene, got$gene)
  de$pvalue <- got$pvalue[hit]
  de$padj <- got$padj[hit]
  log_msg(treat, " vs ", control, " 使用 gene_exp.diff 的 P 值，匹配到 ", sum(!is.na(hit)), " 个基因")
  de
}

select_up <- function(de, fc) {
  if (all(c("log2FC_sh1", "log2FC_sh4") %in% names(de))) {
    keep <- passes_fc(de$log2FC_sh1, fc) & passes_fc(de$log2FC_sh4, fc)
    if (has_real_p(de$pvalue_sh1)) keep <- keep & passes_p(de$pvalue_sh1)
    if (has_real_p(de$pvalue_sh4)) keep <- keep & passes_p(de$pvalue_sh4)
  } else {
    keep <- passes_fc(de$log2FC, fc) & passes_p(de$pvalue)
  }
  out <- de[keep, , drop = FALSE]
  out[order(-out$log2FC), , drop = FALSE]
}

build_common_up <- function(a, b) {
  cols <- c("gene", "log2FC", "pvalue", "padj", "mean_expr")
  m <- merge(a[, cols], b[, cols], by = "gene", suffixes = c("_sh1", "_sh4"))
  m$log2FC <- (m$log2FC_sh1 + m$log2FC_sh4) / 2
  m$mean_expr <- (m$mean_expr_sh1 + m$mean_expr_sh4) / 2
  both_p <- has_real_p(m$pvalue_sh1) && has_real_p(m$pvalue_sh4)
  m$pvalue <- if (both_p) pmax(m$pvalue_sh1, m$pvalue_sh4) else NA_real_
  m$padj <- if (both_p) pmax(m$padj_sh1, m$padj_sh4) else NA_real_
  m
}

# -----------------------------------------------------------------------------
# 火山图、表达热图、气泡图、富集热图
# -----------------------------------------------------------------------------
plot_volcano <- function(df, title, outfile) {
  plot_df <- df[!is.na(df$log2FC), , drop = FALSE]
  if (nrow(plot_df) < 2) {
    log_msg("火山图跳过，可用基因少于 2 个: ", title)
    return(invisible(NULL))
  }
  use_p <- has_real_p(plot_df$pvalue)
  if (use_p) {
    plot_df <- plot_df[!is.na(plot_df$pvalue) & plot_df$pvalue > 0, , drop = FALSE]
    plot_df$y <- -log10(pmax(plot_df$pvalue, 1e-300))
    ylab <- "-log10(P)"
  } else {
    plot_df$y <- plot_df$mean_expr
    ylab <- "Average log2 expression"
  }
  plot_df$show <- ifelse(plot_df$regulation == "Up", "Up", "Other")
  up_genes <- plot_df$gene[plot_df$show == "Up"]
  up_genes <- up_genes[order(plot_df$log2FC[match(up_genes, plot_df$gene)], decreasing = TRUE)]
  labs <- utils::head(up_genes, 10)
  plot_df$label <- ifelse(plot_df$gene %in% labs, plot_df$gene, NA_character_)
  n_up <- sum(plot_df$show == "Up")
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = log2FC, y = y, color = show)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.3) +
    ggplot2::scale_color_manual(values = c(Up = "#B2182B", Other = "grey75"), drop = FALSE) +
    ggplot2::geom_vline(xintercept = c(0, log2(1.25)), linetype = 2, color = "grey40") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = title,
      subtitle = if (use_p) {
        sprintf("上调且 P<%s: %d。竖线为 FC = 1 和 FC = 1.25", p_cutoff, n_up)
      } else {
        sprintf("本比较没有 P 值，红点为上调基因: %d。竖线为 FC = 1 和 FC = 1.25", n_up)
      },
      x = "log2 Fold Change (PCY / NTC)",
      y = ylab,
      color = NULL
    )
  if (use_p) p <- p + ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = 2, color = "grey40")
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
    log_msg("表达热图只画前 80 个上调基因（共 ", length(genes), " 个）")
    genes <- genes[seq_len(80)]
  }
  if (length(genes) < 2) {
    log_msg("表达热图跳过，上调基因少于 2 个")
    note_empty(paste0(outfile, "_EMPTY.txt"), "upregulated genes < 2")
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
  pal <- c(NTC_rep0 = "#4C78A8", NTC_rep1 = "#72B7B2", PCY_sh1 = "#E45756", PCY_sh4 = "#F58518")
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

spread_constant_matrix <- function(mat) {
  vals <- mat[is.finite(mat)]
  if (length(vals) < 2 || length(unique(vals)) >= 2) return(mat)
  idx <- which(is.finite(mat))[1]
  mat[idx] <- vals[1] + max(1e-3, abs(vals[1]) * 1e-3)
  mat
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
  mat <- spread_constant_matrix(mat)
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

# -----------------------------------------------------------------------------
# 胞内囊泡运输专项（不改全基因组 GO/KEGG 的 p 值）
# -----------------------------------------------------------------------------
.vesicle_env <- new.env(parent = emptyenv())

vesicle_class_labels <- c(
  endosome_transport = "内体运输与内体运动",
  vesicle_transport = "囊泡运输",
  exosome_secretion = "外泌体分泌"
)

vesicle_catalog <- function() {
  if (exists("catalog", envir = .vesicle_env, inherits = FALSE)) {
    return(get("catalog", envir = .vesicle_env, inherits = FALSE))
  }
  txt <- "
go_id	term	ontology	vesicle_class
GO:0006897	endocytosis	BP	endosome_transport
GO:0006898	receptor-mediated endocytosis	BP	endosome_transport
GO:0072583	clathrin-dependent endocytosis	BP	endosome_transport
GO:0016197	endosomal transport	BP	endosome_transport
GO:0007032	endosome organization	BP	endosome_transport
GO:0032456	endocytic recycling	BP	endosome_transport
GO:0045022	early endosome to late endosome transport	BP	endosome_transport
GO:0061502	early endosome to recycling endosome transport	BP	endosome_transport
GO:0034498	early endosome to Golgi transport	BP	endosome_transport
GO:0006895	Golgi to endosome transport	BP	endosome_transport
GO:0034499	late endosome to Golgi transport	BP	endosome_transport
GO:0008333	endosome to lysosome transport	BP	endosome_transport
GO:0099638	endosome to plasma membrane protein transport	BP	endosome_transport
GO:0042147	retrograde transport, endosome to Golgi	BP	endosome_transport
GO:0032509	endosome transport via multivesicular body sorting pathway	BP	endosome_transport
GO:0034058	endosomal vesicle fusion	BP	endosome_transport
GO:0005768	endosome	CC	endosome_transport
GO:0005769	early endosome	CC	endosome_transport
GO:0005770	late endosome	CC	endosome_transport
GO:0055037	recycling endosome	CC	endosome_transport
GO:0010008	endosome membrane	CC	endosome_transport
GO:0030139	endocytic vesicle	CC	endosome_transport
GO:0030136	clathrin-coated vesicle	CC	endosome_transport
GO:0045334	clathrin-coated endocytic vesicle	CC	endosome_transport
GO:0016192	vesicle-mediated transport	BP	vesicle_transport
GO:0098927	vesicle-mediated transport between endosomal compartments	BP	vesicle_transport
GO:0098876	vesicle-mediated transport to the plasma membrane	BP	vesicle_transport
GO:0006888	endoplasmic reticulum to Golgi vesicle-mediated transport	BP	vesicle_transport
GO:0006901	vesicle coating	BP	vesicle_transport
GO:0048278	vesicle docking	BP	vesicle_transport
GO:0006906	vesicle fusion	BP	vesicle_transport
GO:0047496	vesicle transport along microtubule	BP	vesicle_transport
GO:0030050	vesicle transport along actin filament	BP	vesicle_transport
GO:0030133	transport vesicle	CC	vesicle_transport
GO:0031410	cytoplasmic vesicle	CC	vesicle_transport
GO:0030137	COPI-coated vesicle	CC	vesicle_transport
GO:0030134	COPII-coated ER to Golgi transport vesicle	CC	vesicle_transport
GO:1990182	exosomal secretion	BP	exosome_secretion
GO:1903541	regulation of exosomal secretion	BP	exosome_secretion
GO:1903551	regulation of extracellular exosome assembly	BP	exosome_secretion
GO:1903552	negative regulation of extracellular exosome assembly	BP	exosome_secretion
GO:1903553	positive regulation of extracellular exosome assembly	BP	exosome_secretion
GO:0071971	extracellular exosome assembly	BP	exosome_secretion
GO:0097734	extracellular exosome biogenesis	BP	exosome_secretion
GO:0140112	extracellular vesicle biogenesis	BP	exosome_secretion
GO:0071985	multivesicular body sorting pathway	BP	exosome_secretion
GO:0036258	multivesicular body assembly	BP	exosome_secretion
GO:0036257	multivesicular body organization	BP	exosome_secretion
GO:0070062	extracellular exosome	CC	exosome_secretion
GO:1903561	extracellular vesicle	CC	exosome_secretion
GO:0005771	multivesicular body	CC	exosome_secretion
GO:0032585	multivesicular body membrane	CC	exosome_secretion
"
  catalog <- utils::read.delim(text = txt, stringsAsFactors = FALSE, strip.white = TRUE)
  catalog <- catalog[nzchar(catalog$go_id), , drop = FALSE]
  catalog$class_label <- vesicle_class_labels[catalog$vesicle_class]
  .vesicle_env$catalog <- catalog
  catalog
}

classify_vesicle_term <- function(id, term) {
  catalog <- vesicle_catalog()
  id <- as.character(id)
  term_l <- tolower(as.character(term))
  out <- rep(NA_character_, length(id))
  known <- match(id, catalog$go_id)
  out[!is.na(known)] <- catalog$vesicle_class[known[!is.na(known)]]
  need <- is.na(out) & !is.na(term_l)
  bad <- grepl(
    paste(
      "viral", "virus", "vitellogenesis", "melanosome", "pigment granule",
      "neurotransmitter", "rnase", "rna polymerase", "nuclear exosome",
      "nucleolar exosome", "bacterial extracellular", "synaptic", "acrosome",
      "acrosomal", "pattern recognition",
      sep = "|"
    ),
    term_l, perl = TRUE
  )
  exo <- grepl("exosomal secretion|extracellular exosome|extracellular vesicle|multivesicular body|via exosome", term_l, perl = TRUE)
  endo <- grepl("endocyt|endosome|endosomal", term_l, perl = TRUE)
  ves <- grepl(
    "vesicle-mediated transport|vesicle transport along|vesicle docking|vesicle fusion|vesicle coating|transport vesicle|cytoplasmic vesicle|copi-coated vesicle|copii-coated",
    term_l, perl = TRUE
  )
  bad[is.na(bad)] <- FALSE
  exo[is.na(exo)] <- FALSE
  endo[is.na(endo)] <- FALSE
  ves[is.na(ves)] <- FALSE
  out[need & exo & !bad] <- "exosome_secretion"
  need <- is.na(out)
  out[need & endo & !bad] <- "endosome_transport"
  need <- is.na(out)
  out[need & ves & !bad] <- "vesicle_transport"
  out
}

enrich_df_ranked <- function(obj) {
  if (is.null(obj)) return(NULL)
  df <- as.data.frame(obj)
  if (nrow(df) == 0) return(NULL)
  df <- df[order(df$pvalue, df$p.adjust), , drop = FALSE]
  df$genome_wide_rank <- seq_len(nrow(df))
  df
}

export_go_vesicle_focus <- function(full_df, out_stub) {
  if (is.null(full_df) || nrow(full_df) == 0) {
    note_empty(paste0(out_stub, "_EMPTY.txt"), "no genome-wide GO terms")
    return(NULL)
  }
  desc <- if ("Description" %in% names(full_df)) full_df$Description else full_df$ID
  full_df$vesicle_class <- classify_vesicle_term(full_df$ID, desc)
  hit <- full_df[!is.na(full_df$vesicle_class), , drop = FALSE]
  if (nrow(hit) == 0) {
    note_empty(paste0(out_stub, "_EMPTY.txt"), "no vesicle-transport GO terms in the genome-wide result")
    return(NULL)
  }
  hit$class_label <- vesicle_class_labels[hit$vesicle_class]
  write_csv(hit, paste0(out_stub, ".csv"))
  hit
}

vesicle_term2gene <- function() {
  if (exists("t2g", envir = .vesicle_env, inherits = FALSE)) {
    got <- get("t2g", envir = .vesicle_env, inherits = FALSE)
    if (!is.null(got)) return(got)
  }
  empty <- data.frame(go_id = character(), entrez = character(), stringsAsFactors = FALSE)
  catalog <- vesicle_catalog()
  if (!requireNamespace("GO.db", quietly = TRUE)) {
    log_msg("未安装 GO.db，跳过囊泡专项基因集")
    .vesicle_env$t2g <- empty
    return(empty)
  }
  mapped <- tryCatch(
    AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = catalog$go_id, column = "ENTREZID",
      keytype = "GOALL", multiVals = "list"
    ),
    error = function(e) {
      log_msg("GOALL 映射失败，改用 GO: ", e$message)
      tryCatch(
        AnnotationDbi::mapIds(
          org.Hs.eg.db, keys = catalog$go_id, column = "ENTREZID",
          keytype = "GO", multiVals = "list"
        ),
        error = function(e2) NULL
      )
    }
  )
  if (is.null(mapped)) {
    .vesicle_env$t2g <- empty
    return(empty)
  }
  rows <- lapply(names(mapped), function(gid) {
    eg <- unique(as.character(mapped[[gid]]))
    eg <- eg[!is.na(eg) & nzchar(eg)]
    if (length(eg) == 0) return(NULL)
    data.frame(go_id = gid, entrez = eg, stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) {
    .vesicle_env$t2g <- empty
    return(empty)
  }
  t2g <- do.call(rbind, rows)
  log_msg(
    "囊泡 GO 基因集: ", length(unique(t2g$go_id)), " 个条目, ",
    length(unique(t2g$entrez)), " 个 Entrez"
  )
  .vesicle_env$t2g <- t2g
  t2g
}

run_vesicle_focus <- function(entrez, universe, up, id_map, heat, sample_info, comparison, comparison_id, outdir, fc_label, fc_symbol, fc_entrez, genome_hits) {
  fdir <- file.path(outdir, "Focused_vesicle_transport")
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(
      "这是胞内囊泡运输专项，不改全基因组 GO 的 p 值或排名。",
      "三类：内体运输与内体运动、囊泡运输、外泌体分泌。",
      "ORA_vesicle_transport.csv 只检验这批 GO。p 值是这批条目内部的超几何检验。",
      "GO/ORA_GO_*_FOCUS_vesicle_transport.csv 来自全基因组 GO，保留原始 p 和 genome_wide_rank。",
      "基因数超过全基因组 maxGSSize=500 的大条目，只会出现在本文件夹。"
    ),
    file.path(fdir, "00_README.txt")
  )
  t2g <- vesicle_term2gene()
  catalog <- vesicle_catalog()
  if (nrow(t2g) == 0 || length(entrez) < 3) {
    note_empty(file.path(fdir, "ORA_vesicle_transport_EMPTY.txt"), "too few genes or empty vesicle GO sets")
    return(invisible(NULL))
  }
  obj <- tryCatch(
    clusterProfiler::enricher(
      gene = entrez, universe = universe, TERM2GENE = t2g[, c("go_id", "entrez")],
      TERM2NAME = catalog[, c("go_id", "term")],
      pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
      minGSSize = 5, maxGSSize = 8000
    ),
    error = function(e) {
      log_msg("囊泡专项 ORA 失败: ", e$message)
      NULL
    }
  )
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    obj <- tryCatch(
      clusterProfiler::setReadable(obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
      error = function(e) obj
    )
  }
  df <- enrich_df_ranked(obj)
  if (is.null(df)) {
    note_empty(file.path(fdir, "ORA_vesicle_transport_EMPTY.txt"), "no vesicle GO term passed the gene-set size filter")
    log_msg(comparison, " ", fc_label, " 囊泡专项没有可检验的 GO")
    return(invisible(NULL))
  }
  df$genome_wide_rank <- NULL
  meta <- catalog[match(df$ID, catalog$go_id), c("vesicle_class", "class_label", "ontology", "term")]
  df$vesicle_class <- meta$vesicle_class
  df$class_label <- meta$class_label
  df$ontology <- meta$ontology
  df$Description <- meta$term
  df <- df[!is.na(df$vesicle_class), , drop = FALSE]
  genome_hits <- Filter(Negate(is.null), genome_hits)
  gw <- NULL
  if (length(genome_hits) > 0) {
    cols <- Reduce(intersect, lapply(genome_hits, names))
    gw <- do.call(rbind, lapply(genome_hits, function(x) x[, cols, drop = FALSE]))
  }
  if (!is.null(gw) && nrow(gw) > 0) {
    gw <- gw[!duplicated(gw$ID), , drop = FALSE]
    hit <- match(df$ID, gw$ID)
    df$genome_wide_pvalue <- gw$pvalue[hit]
    df$genome_wide_padj <- gw$p.adjust[hit]
    df$genome_wide_rank <- gw$genome_wide_rank[hit]
  } else {
    df$genome_wide_pvalue <- NA_real_
    df$genome_wide_padj <- NA_real_
    df$genome_wide_rank <- NA_integer_
  }
  df$comparison <- comparison_id
  df$fc_label <- fc_label
  write_csv(df, file.path(fdir, "ORA_vesicle_transport.csv"))
  log_msg(
    comparison, " ", fc_label, " 囊泡 GO: ", nrow(df),
    " 个条目，其中 P<", p_cutoff, " 的有 ", sum(df$pvalue < p_cutoff, na.rm = TRUE), " 个"
  )
  plot_df <- df
  plot_df$Description <- paste(plot_df$class_label, plot_df$Description, sep = " | ")
  sig <- plot_df[!is.na(plot_df$pvalue) & plot_df$pvalue < p_cutoff, , drop = FALSE]
  relaxed <- nrow(sig) == 0
  show <- if (relaxed) plot_df else sig
  save_enrichment(
    show, relaxed,
    file.path(fdir, "ORA_vesicle_transport"),
    paste0(comparison, " ", fc_label, " 囊泡运输相关 GO"),
    fc_symbol, fc_entrez
  )
  # save_enrichment already wrote the csv of the plotted subset over a different stub.
  # The full table remains ORA_vesicle_transport.csv; plots use the same stub and would
  # overwrite that csv with only the plotted rows. Write the full table again.
  write_csv(df, file.path(fdir, "ORA_vesicle_transport.csv"))

  core_ids <- c("GO:0006897", "GO:0016197", "GO:0016192", "GO:1990182", "GO:0071985", "GO:0047496", "GO:0030050")
  sig_bp <- df$ID[!is.na(df$pvalue) & df$pvalue < p_cutoff & df$ontology == "BP"]
  sig_ids <- df$ID[!is.na(df$pvalue) & df$pvalue < p_cutoff]
  use_ids <- if (length(sig_bp) > 0) sig_bp else if (length(sig_ids) > 0) sig_ids else core_ids
  t2g_use <- t2g[t2g$go_id %in% use_ids, , drop = FALSE]
  up_map <- id_map[id_map$gene %in% up$gene, , drop = FALSE]
  up_hit <- up_map[up_map$entrez %in% t2g_use$entrez, , drop = FALSE]
  if (nrow(up_hit) > 0) {
    class_of <- tapply(t2g_use$go_id, t2g_use$entrez, function(gids) {
      paste(unique(catalog$class_label[match(gids, catalog$go_id)]), collapse = ";")
    })
    gene_tbl <- data.frame(
      gene = up_hit$gene,
      entrez = up_hit$entrez,
      log2FC = up$log2FC[match(up_hit$gene, up$gene)],
      vesicle_class_label = unname(class_of[up_hit$entrez]),
      stringsAsFactors = FALSE
    )
    gene_tbl <- gene_tbl[order(-gene_tbl$log2FC), , drop = FALSE]
    gene_tbl <- gene_tbl[!duplicated(gene_tbl$gene), , drop = FALSE]
    write_csv(gene_tbl, file.path(fdir, "vesicle_up_genes.csv"))
    tryCatch(
      plot_expression_heatmap(
        heat, sample_info, gene_tbl$gene,
        paste0(comparison, " ", fc_label, " 囊泡运输相关上调基因"),
        file.path(fdir, "vesicle_up_expression_heatmap")
      ),
      error = function(e) log_msg("囊泡表达热图失败: ", e$message)
    )
  } else {
    note_empty(file.path(fdir, "vesicle_up_genes_EMPTY.txt"), "no upregulated gene mapped to the vesicle GO sets")
  }

  keep <- c("comparison", "fc_label", "ID", "Description", "vesicle_class", "class_label", "ontology", "pvalue", "p.adjust", "Count")
  # Description in df is still the GO term name; plot_df changed a copy.
  slim <- df[, intersect(keep, names(df)), drop = FALSE]
  .vesicle_env$focused[[length(.vesicle_env$focused) + 1L]] <- slim
  invisible(df)
}

write_vesicle_across <- function(result_dir) {
  outdir <- file.path(result_dir, "00_vesicle_transport_across_6")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(
      "六类比较里，胞内囊泡运输相关 GO 是否一起变化。",
      "列的顺序是四组 1 对 1、两个 knockdown 的均值对两个 NTC、",
      "相对每个 NTC 的共同上调、均值分别对 NTC_rep0 和 NTC_rep1。",
      "数值是专项检验的 -log10(p.adjust)，不是改过的全基因组 p 值。",
      "比较之间的相关性是这些 -log10(p.adjust) 的 Pearson 相关。",
      "FC_1 和 FC_1.25 分开，因为用的上调基因集合不同。"
    ),
    file.path(outdir, "00_README.txt")
  )
  hits <- .vesicle_env$focused
  if (length(hits) == 0) {
    note_empty(file.path(outdir, "vesicle_across_EMPTY.txt"), "no focused vesicle GO results")
    log_msg("没有囊泡专项结果，跳过六类比较汇总")
    return(invisible(NULL))
  }
  all <- do.call(rbind, hits)
  write_csv(all, file.path(outdir, "vesicle_focused_GO_all_comparisons.csv"))
  comps <- .vesicle_env$comparison_order
  if (is.null(comps) || length(comps) == 0) comps <- unique(all$comparison)
  draw_heat <- function(mat, title, stub) {
    if (nrow(mat) < 1 || ncol(mat) < 1 || sum(is.finite(mat)) < 1) {
      note_empty(paste0(stub, "_EMPTY.txt"), "nothing to plot")
      return(invisible(NULL))
    }
    plot_mat <- spread_constant_matrix(mat)
    plot_mat[is.finite(plot_mat) & plot_mat > 20] <- 20
    diverging <- isTRUE(min(plot_mat, na.rm = TRUE) < 0 && max(plot_mat, na.rm = TRUE) <= 1)
    if (diverging) plot_mat[is.finite(plot_mat) & plot_mat < -1] <- -1
    width <- max(8, min(16, 0.7 * ncol(plot_mat) + 4))
    height <- max(4, min(16, 0.28 * nrow(plot_mat) + 2))
    draw <- function() {
      cols <- if (diverging) {
        grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
      } else {
        grDevices::colorRampPalette(c("white", "#F4A582", "#B2182B"))(100)
      }
      args <- list(
        plot_mat,
        cluster_rows = nrow(plot_mat) >= 2 && !anyNA(plot_mat),
        cluster_cols = ncol(plot_mat) >= 2 && !anyNA(plot_mat),
        na_col = "grey92",
        color = cols,
        main = title,
        fontsize_row = 7,
        fontsize_col = 8,
        border_color = NA
      )
      if (diverging) args$breaks <- seq(-1, 1, length.out = 101)
      do.call(pheatmap::pheatmap, args)
    }
    grDevices::pdf(paste0(stub, ".pdf"), width = width, height = height)
    on.exit(while (grDevices::dev.cur() > 1) grDevices::dev.off(), add = TRUE)
    tryCatch(draw(), error = function(e) log_msg("囊泡汇总热图失败: ", e$message))
    grDevices::dev.off()
    grDevices::png(paste0(stub, ".png"), width = width, height = height, units = "in", res = 200)
    tryCatch(draw(), error = function(e) log_msg("囊泡汇总热图失败: ", e$message))
    grDevices::dev.off()
  }
  for (fc in names(fc_levels)) {
    sub <- all[all$fc_label == fc, , drop = FALSE]
    fc_dir <- file.path(outdir, fc)
    dir.create(fc_dir, recursive = TRUE, showWarnings = FALSE)
    if (nrow(sub) == 0) {
      note_empty(file.path(fc_dir, "EMPTY.txt"), "no rows for this FC")
      next
    }
    terms <- unique(sub$ID)
    row_lab <- vapply(terms, function(gid) {
      hit <- sub[sub$ID == gid, , drop = FALSE][1, ]
      paste(hit$class_label, hit$Description, sep = " | ")
    }, character(1))
    mat <- matrix(NA_real_, nrow = length(terms), ncol = length(comps), dimnames = list(row_lab, comps))
    for (i in seq_len(nrow(sub))) {
      rid <- match(sub$ID[i], terms)
      cid <- match(sub$comparison[i], comps)
      if (is.na(rid) || is.na(cid)) next
      padj <- sub$p.adjust[i]
      if (is.na(padj)) padj <- sub$pvalue[i]
      if (!is.na(padj) && padj > 0) mat[rid, cid] <- -log10(padj)
    }
    write_csv(
      data.frame(term = row_lab, go_id = terms, mat, check.names = FALSE),
      file.path(fc_dir, "vesicle_GO_neglogp_by_comparison.csv")
    )
    draw_heat(mat, paste0(fc, " 囊泡 GO 在六类比较中的 -log10(p.adjust)"), file.path(fc_dir, "vesicle_GO_neglogp_heatmap"))
    class_rows <- lapply(unique(sub$vesicle_class), function(cl) {
      piece <- sub[sub$vesicle_class == cl, , drop = FALSE]
      lapply(comps, function(comp) {
        x <- piece[piece$comparison == comp, , drop = FALSE]
        data.frame(
          fc_label = fc,
          comparison = comp,
          vesicle_class = cl,
          class_label = vesicle_class_labels[[cl]],
          n_terms = nrow(x),
          n_p_lt_cutoff = sum(!is.na(x$pvalue) & x$pvalue < p_cutoff),
          best_p = if (nrow(x) == 0 || all(is.na(x$pvalue))) NA_real_ else min(x$pvalue, na.rm = TRUE),
          best_padj = if (nrow(x) == 0 || all(is.na(x$p.adjust))) NA_real_ else min(x$p.adjust, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      })
    })
    class_df <- do.call(rbind, unlist(class_rows, recursive = FALSE))
    write_csv(class_df, file.path(fc_dir, "vesicle_class_summary.csv"))
    class_mat <- matrix(NA_real_, nrow = length(vesicle_class_labels), ncol = length(comps), dimnames = list(vesicle_class_labels, comps))
    for (i in seq_len(nrow(class_df))) {
      class_mat[class_df$class_label[i], class_df$comparison[i]] <- -log10(pmax(class_df$best_padj[i], 1e-300))
    }
    class_mat[!is.finite(class_mat)] <- NA_real_
    draw_heat(class_mat, paste0(fc, " 三类囊泡通路在六类比较中的 -log10(最佳 p.adjust)"), file.path(fc_dir, "vesicle_class_neglogp_heatmap"))
    usable <- mat[, colSums(is.finite(mat)) >= 3, drop = FALSE]
    if (ncol(usable) >= 2 && nrow(usable) >= 3) {
      corm <- stats::cor(usable, use = "pairwise.complete.obs")
      write_csv(
        data.frame(comparison = rownames(corm), corm, check.names = FALSE),
        file.path(fc_dir, "vesicle_comparison_correlation.csv")
      )
      draw_heat(corm, paste0(fc, " 六类比较的囊泡 GO 富集相关性"), file.path(fc_dir, "vesicle_comparison_correlation_heatmap"))
    } else {
      note_empty(file.path(fc_dir, "vesicle_comparison_correlation_EMPTY.txt"), "fewer than 3 shared vesicle terms")
    }
  }
  log_msg("囊泡专项六类比较汇总: ", outdir)
  invisible(outdir)
}

run_downstream <- function(sub, de, heat, sample_info, comparison, outdir, id_map, fc_label, comparison_id = comparison) {
  up <- sub
  write_csv(up, file.path(outdir, paste0(fc_label, "_genes.csv")))
  log_msg(comparison, " ", fc_label, " 上调基因: ", nrow(up))
  tryCatch(
    plot_expression_heatmap(
      heat, sample_info, up$gene,
      paste0(comparison, " ", fc_label, " 上调基因"),
      file.path(outdir, paste0(fc_label, "_expression_heatmap"))
    ),
    error = function(e) log_msg("表达热图失败: ", e$message)
  )
  if (nrow(up) < 3) {
    note_empty(file.path(outdir, "enrichment_EMPTY.txt"), "upregulated genes < 3; skip GO/KEGG/pathway")
    log_msg(comparison, " 上调基因少于 3 个，跳过 GO/KEGG/通路")
    return(invisible(NULL))
  }

  all_genes <- de$gene[!grepl("^(XLOC|TCONS|CUFF)_", de$gene)]
  mp_all <- id_map[id_map$gene %in% all_genes, , drop = FALSE]
  mp_up <- id_map[id_map$gene %in% up$gene, , drop = FALSE]
  universe <- unique(mp_all$entrez)
  entrez <- unique(mp_up$entrez)
  entrez <- intersect(entrez, universe)
  log_msg(comparison, " 上调基因映射到 Entrez: ", length(entrez), " / ", nrow(up))
  if (length(entrez) < 3) {
    note_empty(file.path(outdir, "enrichment_EMPTY.txt"), "fewer than 3 genes mapped to Entrez")
    return(invisible(NULL))
  }

  fc_symbol <- setNames(up$log2FC, up$gene)
  fc_entrez <- setNames(up$log2FC[match(mp_up$gene, up$gene)], mp_up$entrez)
  fc_entrez <- fc_entrez[!duplicated(names(fc_entrez)) & nzchar(names(fc_entrez))]

  go_dir <- file.path(outdir, "GO")
  kegg_dir <- file.path(outdir, "KEGG")
  pw_dir <- file.path(outdir, "Pathway")
  dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(kegg_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pw_dir, recursive = TRUE, showWarnings = FALSE)

  vesicle_from_go <- list()
  for (ont in c("BP", "CC", "MF")) {
    ont_name <- c(BP = "Biological Process", CC = "Cellular Component", MF = "Molecular Function")[[ont]]
    obj <- tryCatch(
      clusterProfiler::enrichGO(
        gene = entrez, universe = universe, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
        ont = ont, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
        minGSSize = 5, maxGSSize = 500, readable = TRUE
      ),
      error = function(e) {
        log_msg("GO ", ont, " 失败: ", e$message)
        NULL
      }
    )
    full_df <- enrich_df_ranked(obj)
    if (is.null(full_df)) {
      note_empty(file.path(go_dir, paste0("ORA_GO_", ont, "_EMPTY.txt")), "no GO terms")
      log_msg("没有富集条目: ", comparison, " ", fc_label, " GO ", ont_name)
    } else {
      sig <- full_df[!is.na(full_df$pvalue) & full_df$pvalue < p_cutoff, , drop = FALSE]
      relaxed <- nrow(sig) == 0
      show <- if (relaxed) full_df else sig
      show$genome_wide_rank <- NULL
      show$vesicle_class <- NULL
      show$class_label <- NULL
      save_enrichment(
        show, relaxed,
        file.path(go_dir, paste0("ORA_GO_", ont)),
        paste0(comparison, " ", fc_label, " 上调基因 GO ", ont_name),
        fc_symbol, fc_entrez
      )
      vesicle_from_go[[ont]] <- export_go_vesicle_focus(
        full_df, file.path(go_dir, paste0("ORA_GO_", ont, "_FOCUS_vesicle_transport"))
      )
    }
  }
  tryCatch(
    run_vesicle_focus(
      entrez, universe, up, id_map, heat, sample_info, comparison, comparison_id,
      outdir, fc_label, fc_symbol, fc_entrez, vesicle_from_go
    ),
    error = function(e) log_msg("囊泡专项失败: ", e$message)
  )

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
    file.path(kegg_dir, "ORA_KEGG"),
    paste0(comparison, " ", fc_label, " 上调基因 KEGG"),
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
    file.path(pw_dir, "ORA_Reactome"),
    paste0(comparison, " ", fc_label, " 上调基因 Reactome 通路"),
    fc_symbol, fc_entrez
  )
}

write_readme <- function(result_dir) {
  txt <- c(
    "PCY_RNA 结果说明",
    "四个样品：NTC_rep0、NTC_rep1、PCY_sh1、PCY_sh4。1-vs-1 不合并两个 NTC。",
    "log2FC = log2(PCY / NTC)。正值表示 PCY 高于 NTC（上调）。",
    "上调基因分层：有 P 值时要求 P < 0.05，再分别取 FC > 1 和 FC > 1.25。",
    "1-vs-1 没有生物学重复时不算 P 值；同目录若有 gene_exp.diff，用它的 P 值。",
    "两个 knockdown 对两个 NTC、或对其中一个 NTC 的比较，在能估计时用 limma 的 P 值。",
    "",
    "六类比较目录：",
    "  PCY_sh1_vs_NTC_rep0/  PCY_sh4_vs_NTC_rep0/",
    "  PCY_sh1_vs_NTC_rep1/  PCY_sh4_vs_NTC_rep1/",
    "  PCYsh_mean_vs_NTC/",
    "  common_up_vs_NTC_rep0/  common_up_vs_NTC_rep1/",
    "  PCYsh_mean_vs_NTC_rep0/  PCYsh_mean_vs_NTC_rep1/",
    "",
    "每个比较：",
    "  DEG_all.csv、volcano",
    "  FoldChange/FC_1/ 与 FoldChange/FC_1.25/",
    "    *_genes.csv、表达热图",
    "    GO/  KEGG/  Pathway/   文件名以 ORA_ 开头，含气泡图和富集热图",
    "    GO/ORA_GO_*_FOCUS_vesicle_transport.csv  全基因组结果里的囊泡 GO，保留原始 p 和 genome_wide_rank",
    "    Focused_vesicle_transport/  内体运输与内体运动、囊泡运输、外泌体分泌",
    "",
    "00_vesicle_transport_across_6/ 把六类比较的囊泡 GO 放在一起，",
    "热图是 -log10(p.adjust)，相关性是这些数值在比较之间的 Pearson 相关。",
    "",
    "物种注释是人 org.Hs.eg.db。没有官方基因符号的 XLOC 会留在差异表，但通常进不了 GO/KEGG。"
  )
  writeLines(txt, file.path(result_dir, "00_README.txt"))
}

emit_comparison <- function(de, name, title, method, heat, sample_info, result_dir, id_map, summary_rows) {
  log_msg("开始 ", title)
  component_p <- all(c("pvalue_sh1", "pvalue_sh4") %in% names(de)) &&
    (has_real_p(de$pvalue_sh1) || has_real_p(de$pvalue_sh4))
  if (!has_real_p(de$pvalue) && !component_p) {
    log_msg(title, " 没有 P 值，FC 分层不再按 P 过滤")
  }
  outdir <- file.path(result_dir, name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  de <- de_export(de, name, method)
  write_csv(de, file.path(outdir, "DEG_all.csv"))
  tryCatch(plot_volcano(de, title, file.path(outdir, "volcano")), error = function(e) {
    log_msg("火山图失败: ", e$message)
  })
  for (tag in names(fc_levels)) {
    sub <- select_up(de, fc_levels[[tag]])
    log_msg(title, " ", tag, "（FC>", fc_levels[[tag]], "）: ", nrow(sub), " 个上调基因")
    fc_dir <- file.path(outdir, "FoldChange", tag)
    dir.create(fc_dir, recursive = TRUE, showWarnings = FALSE)
    tryCatch(
      run_downstream(sub, de, heat, sample_info, title, fc_dir, id_map, tag, name),
      error = function(e) log_msg(name, " ", tag, " 下游分析失败: ", e$message)
    )
  }
  summary_rows[[name]] <- data.frame(
    comparison = name,
    method = method,
    genes = nrow(de),
    up_FC_1 = nrow(select_up(de, 1)),
    up_FC_1.25 = nrow(select_up(de, 1.25)),
    has_pvalue = has_real_p(de$pvalue) || ("pvalue_sh1" %in% names(de) && has_real_p(de$pvalue_sh1)),
    stringsAsFactors = FALSE
  )
  summary_rows
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
  log_msg("上调基因：P < ", p_cutoff, "（能估计时），以及 FC > 1、FC > 1.25")

  expr <- load_expression(project_dir)
  mat <- prepare_matrix(expr$mat)
  sample_info <- expr$sample_info
  sample_info <- sample_info[sample_info$sample %in% colnames(mat), , drop = FALSE]
  sample_info <- align_samples(mat, sample_info)
  log_msg("样品: ", paste(paste0(sample_info$sample, "=", sample_info$group), collapse = "; "))
  present <- as.character(sample_info$group)
  log_msg("每组样品数: ", paste(paste0(sample_levels, "=", vapply(sample_levels, function(id) sum(present == id), integer(1))), collapse = ", "))
  miss <- setdiff(sample_levels, present)
  if (length(miss) > 0) {
    stop(
      "缺少组别 ", paste(miss, collapse = ", "),
      "。需要 NTC_rep0、NTC_rep1、PCY_sh1、PCY_sh4。当前样品: ",
      paste(paste0(sample_info$sample, "=", sample_info$group), collapse = "; "),
      call. = FALSE
    )
  }
  write_csv(data.frame(sample_info, stringsAsFactors = FALSE), file.path(result_dir, "00_sample_info.csv"))
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
  norm <- normalize_expression(mat, value_type)
  log_mat <- norm$log_mat
  heat <- norm$heat
  cols <- stats::setNames(vapply(sample_levels, function(id) sample_column(sample_info, id), character(1)), sample_levels)

  one_vs_one <- function(treat, control, method_base) {
    de <- fc_from_columns(log_mat, cols[[treat]], cols[[control]])
    de <- attach_cuffdiff_p(de, project_dir, treat, control)
    if (!has_real_p(de$pvalue)) {
      log_msg(treat, " vs ", control, " 是 1 对 1，没有重复，不计算 P 值")
    }
    list(de = de, method = paste0(method_base, if (has_real_p(de$pvalue)) "+Cuffdiff_P" else "+no_P"))
  }
  with_p <- function(treat_cols, control_cols, method_base) {
    de <- fc_from_columns(log_mat, treat_cols, control_cols)
    de <- attach_limma_p(de, log_mat, treat_cols, control_cols)
    if (has_real_p(de$pvalue)) {
      log_msg("limma P 值：", paste(treat_cols, collapse = ","), " vs ", paste(control_cols, collapse = ","))
    } else {
      log_msg("该比较无法估计 P 值：", paste(treat_cols, collapse = ","), " vs ", paste(control_cols, collapse = ","))
    }
    list(de = de, method = paste0(method_base, if (has_real_p(de$pvalue)) "+limma_P" else "+no_P"))
  }

  base_method <- norm$method
  sets <- list(
    one_vs_one("PCY_sh1", "NTC_rep0", base_method),
    one_vs_one("PCY_sh4", "NTC_rep0", base_method),
    one_vs_one("PCY_sh1", "NTC_rep1", base_method),
    one_vs_one("PCY_sh4", "NTC_rep1", base_method)
  )
  names(sets) <- c("PCY_sh1_vs_NTC_rep0", "PCY_sh4_vs_NTC_rep0", "PCY_sh1_vs_NTC_rep1", "PCY_sh4_vs_NTC_rep1")
  sets$PCYsh_mean_vs_NTC <- with_p(
    c(cols[["PCY_sh1"]], cols[["PCY_sh4"]]),
    c(cols[["NTC_rep0"]], cols[["NTC_rep1"]]),
    base_method
  )
  sets$PCYsh_mean_vs_NTC_rep0 <- with_p(
    c(cols[["PCY_sh1"]], cols[["PCY_sh4"]]),
    cols[["NTC_rep0"]],
    base_method
  )
  sets$PCYsh_mean_vs_NTC_rep1 <- with_p(
    c(cols[["PCY_sh1"]], cols[["PCY_sh4"]]),
    cols[["NTC_rep1"]],
    base_method
  )
  sets$common_up_vs_NTC_rep0 <- list(
    de = build_common_up(sets$PCY_sh1_vs_NTC_rep0$de, sets$PCY_sh4_vs_NTC_rep0$de),
    method = paste0(base_method, "+common_up")
  )
  sets$common_up_vs_NTC_rep1 <- list(
    de = build_common_up(sets$PCY_sh1_vs_NTC_rep1$de, sets$PCY_sh4_vs_NTC_rep1$de),
    method = paste0(base_method, "+common_up")
  )
  titles <- c(
    PCY_sh1_vs_NTC_rep0 = "PCY_sh1 vs NTC_rep0",
    PCY_sh4_vs_NTC_rep0 = "PCY_sh4 vs NTC_rep0",
    PCY_sh1_vs_NTC_rep1 = "PCY_sh1 vs NTC_rep1",
    PCY_sh4_vs_NTC_rep1 = "PCY_sh4 vs NTC_rep1",
    PCYsh_mean_vs_NTC = "mean(PCY_sh1, PCY_sh4) vs mean(NTC_rep0, NTC_rep1)",
    common_up_vs_NTC_rep0 = "PCY_sh1 与 PCY_sh4 相对 NTC_rep0 的共同上调",
    common_up_vs_NTC_rep1 = "PCY_sh1 与 PCY_sh4 相对 NTC_rep1 的共同上调",
    PCYsh_mean_vs_NTC_rep0 = "mean(PCY_sh1, PCY_sh4) vs NTC_rep0",
    PCYsh_mean_vs_NTC_rep1 = "mean(PCY_sh1, PCY_sh4) vs NTC_rep1"
  )

  write_readme(result_dir)
  .vesicle_env$focused <- list()
  .vesicle_env$comparison_order <- names(titles)
  id_map <- map_to_entrez(rownames(mat))
  log_msg("全基因映射到 Entrez: ", nrow(id_map), " / ", nrow(mat))
  summary_rows <- list()
  for (nm in names(titles)) {
    summary_rows <- tryCatch(
      emit_comparison(sets[[nm]]$de, nm, titles[[nm]], sets[[nm]]$method, heat, sample_info, result_dir, id_map, summary_rows),
      error = function(e) {
        log_msg(nm, " 失败: ", e$message)
        summary_rows
      }
    )
  }
  if (length(summary_rows) > 0) {
    write_csv(do.call(rbind, summary_rows), file.path(result_dir, "00_summary.csv"))
  }
  tryCatch(
    write_vesicle_across(result_dir),
    error = function(e) log_msg("囊泡六类比较汇总失败: ", e$message)
  )
  log_msg("完成。结果在: ", result_dir)
  invisible(result_dir)
}

main()
