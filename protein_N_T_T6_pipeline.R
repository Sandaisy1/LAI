#!/usr/bin/env Rscript
# =============================================================================
# N / T / T6 蛋白质组分析（DIA-NN report.pg_matrix）
#
# 用法（在数据目录或仓库根运行均可）：
#   setwd("E:/R/protein T T6")
#   source("protein_N_T_T6_pipeline.R")
# 或：
#   Rscript protein_N_T_T6_pipeline.R
#
# 只读 report.pg_matrix（可带 .tsv/.txt）。不要用 report.pr_matrix。
# 比较：
#   1) N vs T   （N1,N2,N3,N5,N7,N10 vs T1–T6；参考为 T）
#   2) T6 vs T  （T6-1–T6-7 vs T1–T6；样品 T6 属于 T 组，不属于 T6 组）
# 每个比较：全表 + 火山图 + 仅上调蛋白的 GO / KEGG。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")
options(clusterProfiler.download.method = "auto")

# -----------------------------------------------------------------------------
# 0. 依赖
# -----------------------------------------------------------------------------
cran_required <- c(
  "dplyr", "tidyr", "tibble", "stringr", "ggplot2", "ggrepel", "matrixStats"
)
bioc_required <- c("limma", "AnnotationDbi", "org.Hs.eg.db")
bioc_ora <- c("clusterProfiler", "enrichplot", "GO.db")

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
  if (length(still) > 0 && required) stop("缺少必需 R 包: ", paste(still, collapse = ", "))
  if (length(still) > 0) message("可选包未安装，相关分析将跳过: ", paste(still, collapse = ", "))
  invisible(TRUE)
}

if (!isTRUE(as.logical(Sys.getenv("PROTEIN_NT_T6_SKIP_INSTALL", "false")))) {
  install_if_missing(cran_required, bioc = FALSE, required = TRUE)
  install_if_missing(bioc_required, bioc = TRUE, required = TRUE)
  install_if_missing(bioc_ora, bioc = TRUE, required = FALSE)
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
    !requireNamespace("AnnotationDbi", quietly = TRUE)) {
  stop("GO/KEGG 需要 org.Hs.eg.db 与 AnnotationDbi，请先安装后再跑脚本，不要跳过富集。")
}
for (p in c(cran_required, bioc_required, bioc_ora)) {
  if (requireNamespace(p, quietly = TRUE)) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}
has_clusterProfiler <- requireNamespace("clusterProfiler", quietly = TRUE)
has_enrichplot <- requireNamespace("enrichplot", quietly = TRUE)
has_godb <- requireNamespace("GO.db", quietly = TRUE)

# -----------------------------------------------------------------------------
# 1. 路径与参数
# -----------------------------------------------------------------------------
N_SAMPLES  <- c("N1", "N2", "N3", "N5", "N7", "N10")
T_SAMPLES  <- c("T1", "T2", "T3", "T4", "T5", "T6")
T6_SAMPLES <- c("T6-1", "T6-2", "T6-3", "T6-4", "T6-5", "T6-6", "T6-7")
ALL_SAMPLES <- c(N_SAMPLES, T_SAMPLES, T6_SAMPLES)

P_CUTOFF <- 0.01
FC_CUTOFF <- 1.5
MIN_DETECTED_FRAC <- 0.5

pg_matrix_names <- c(
  "report.pg_matrix", "report.pg_matrix.tsv", "report.pg_matrix.txt",
  "report.pg_matrix.tsv.txt"
)

has_pg_matrix <- function(d) {
  any(file.exists(file.path(d, pg_matrix_names)))
}

resolve_project_dir <- function() {
  env_dir <- Sys.getenv("PROTEIN_NT_T6_DIR", unset = "")
  candidates <- c(
    env_dir,
    "E:/R/protein T T6",
    "E:\\R\\protein T T6",
    getwd()
  )
  script <- tryCatch(
    dirname(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE)),
    error = function(e) NA_character_
  )
  if (!is.na(script)) candidates <- c(candidates, script)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    candidates <- c(candidates, dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)))
  }
  candidates <- unique(candidates[nzchar(candidates)])
  for (d in candidates) {
    if (dir.exists(d) && has_pg_matrix(d)) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

find_pg_matrix <- function(project_dir) {
  hits <- file.path(project_dir, pg_matrix_names)
  hits <- hits[file.exists(hits)]
  if (length(hits) == 0) {
    extra <- list.files(project_dir, pattern = "^report\\.pg_matrix", full.names = TRUE)
    extra <- extra[!grepl("pr_matrix", basename(extra), ignore.case = TRUE)]
    hits <- extra
  }
  if (length(hits) == 0) {
    stop(
      "找不到 report.pg_matrix。请把 DIA-NN 蛋白组矩阵放到：", project_dir,
      "\n不要使用 report.pr_matrix。"
    )
  }
  hits[1]
}

project_dir <- resolve_project_dir()
result_dir  <- file.path(project_dir, "results")
log_dir     <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("protein_pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_msg("项目目录: ", project_dir)

# -----------------------------------------------------------------------------
# 2. 样品列匹配：先 T6-1…T6-7，再 T1–T6，再 N*
# -----------------------------------------------------------------------------
basename_token <- function(col) {
  x <- gsub("\\\\", "/", as.character(col))
  x <- sub("^.*/", "", x)
  sub("\\.(raw|dia|mzml|wiff|d)$", "", x, ignore.case = TRUE)
}

is_annotation_col <- function(nm) {
  low <- tolower(gsub("[^A-Za-z0-9]", "", nm))
  low %in% c(
    "proteingroup", "proteinids", "proteinnames", "genes",
    "firstproteindescription", "proteinid", "majorityproteinids",
    "fastaheaders", "nproteins", "nsequences", "nproteotypicsequences",
    "firstproteinids", "organisms", "organism", "proteindescription",
    "pgmaxlfq", "nfiles"
  ) || grepl("^(protein|gene|fasta|organism|description|nproteins|nsequences)", low)
}

match_t6_id <- function(token) {
  m <- regexpr("(?:^|[^A-Za-z0-9])T6[-_.]([1-7])(?:[^0-9]|$)", token, ignore.case = TRUE, perl = TRUE)
  if (m[1] < 0) return(NA_character_)
  cap <- regmatches(token, m)
  d <- sub(".*T6[-_.]([1-7]).*", "\\1", cap, ignore.case = TRUE)
  paste0("T6-", d)
}

match_t_id <- function(token) {
  if (!is.na(match_t6_id(token))) return(NA_character_)
  m <- regexpr("(?:^|[^A-Za-z0-9])T([1-6])(?:[^0-9]|$)", token, ignore.case = TRUE, perl = TRUE)
  if (m[1] < 0) return(NA_character_)
  cap <- regmatches(token, m)
  d <- sub(".*T([1-6]).*", "\\1", cap, ignore.case = TRUE)
  paste0("T", d)
}

match_n_id <- function(token) {
  if (grepl("(?:^|[^A-Za-z0-9])N10(?:[^0-9]|$)", token, ignore.case = TRUE, perl = TRUE)) {
    return("N10")
  }
  m <- regexpr("(?:^|[^A-Za-z0-9])N([12357])(?:[^0-9]|$)", token, ignore.case = TRUE, perl = TRUE)
  if (m[1] < 0) return(NA_character_)
  cap <- regmatches(token, m)
  d <- sub(".*N([12357]).*", "\\1", cap, ignore.case = TRUE)
  paste0("N", d)
}

classify_sample_col <- function(col) {
  token <- basename_token(col)
  t6 <- match_t6_id(token)
  if (!is.na(t6)) return(list(id = t6, group = "T6"))
  tt <- match_t_id(token)
  if (!is.na(tt)) return(list(id = tt, group = "T"))
  nn <- match_n_id(token)
  if (!is.na(nn)) return(list(id = nn, group = "N"))
  list(id = NA_character_, group = NA_character_)
}

pick_id_column <- function(nms) {
  low <- tolower(nms)
  for (k in c("protein.group", "protein.ids", "protein.ids", "protein.id")) {
    hit <- which(low == k)
    if (length(hit) > 0) return(nms[hit[1]])
  }
  hit <- which(grepl("protein\\.(group|ids|id)$", low))
  if (length(hit) > 0) return(nms[hit[1]])
  nms[1]
}

pick_gene_column <- function(nms) {
  low <- tolower(nms)
  for (k in c("genes", "gene", "gene.names", "gene.name", "gene.symbol")) {
    hit <- which(low == k)
    if (length(hit) > 0) return(nms[hit[1]])
  }
  NA_character_
}

pick_official_symbol <- function(x) {
  x <- trimws(as.character(x)[1])
  if (length(x) != 1 || is.na(x) || x %in% c("", "-", ".", "NA")) return(NA_character_)
  parts <- unlist(strsplit(x, "[,;|/]+"))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts) & !parts %in% c("-", ".", "NA")]
  if (length(parts) == 0) return(NA_character_)
  parts[1]
}

first_accession <- function(x) {
  x <- trimws(as.character(x)[1])
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  parts <- unlist(strsplit(x, "[;|,]+"))
  trimws(parts[1])
}

# -----------------------------------------------------------------------------
# 3. 读入 pg_matrix
# -----------------------------------------------------------------------------
read_pg_matrix <- function(path) {
  log_msg("读取蛋白组矩阵: ", path)
  if (grepl("pr_matrix", basename(path), ignore.case = TRUE)) {
    stop("拒绝读取 report.pr_matrix；请使用 report.pg_matrix")
  }
  raw <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(raw) < 3) stop("矩阵列数过少，请确认文件是制表符分隔的 report.pg_matrix")
  nms <- names(raw)
  id_col <- pick_id_column(nms)
  gene_col <- pick_gene_column(nms)
  log_msg("蛋白 ID 列: ", id_col, if (is.na(gene_col)) "" else paste0("; 基因列: ", gene_col))

  sample_info <- list()
  used_ids <- character()
  for (col in nms) {
    if (identical(col, id_col) || identical(col, gene_col) || is_annotation_col(col)) next
    if (!is.numeric(raw[[col]])) {
      suppressWarnings(num <- as.numeric(raw[[col]]))
      if (mean(is.na(num)) > 0.5) next
      raw[[col]] <- num
    }
    hit <- classify_sample_col(col)
    if (is.na(hit$id)) {
      log_msg("未匹配样品列（忽略）: ", col)
      next
    }
    if (hit$id %in% used_ids) {
      log_msg("样品 ", hit$id, " 重复匹配，保留先出现的列，忽略: ", col)
      next
    }
    used_ids <- c(used_ids, hit$id)
    sample_info[[length(sample_info) + 1]] <- data.frame(
      sample = hit$id, group = hit$group, column = col, stringsAsFactors = FALSE
    )
  }
  if (length(sample_info) == 0) stop("没有匹配到 N / T / T6 样品列")
  sample_info <- do.call(rbind, sample_info)
  sample_info$group <- factor(sample_info$group, levels = c("N", "T", "T6"))
  sample_info <- sample_info[order(sample_info$group, sample_info$sample), , drop = FALSE]

  missing <- setdiff(ALL_SAMPLES, sample_info$sample)
  if (length(missing) > 0) {
    log_msg("警告：未找到这些样品列: ", paste(missing, collapse = ", "))
  }
  extra_n <- setdiff(sample_info$sample[sample_info$group == "N"], N_SAMPLES)
  extra_t <- setdiff(sample_info$sample[sample_info$group == "T"], T_SAMPLES)
  extra_t6 <- setdiff(sample_info$sample[sample_info$group == "T6"], T6_SAMPLES)
  if (length(c(extra_n, extra_t, extra_t6)) > 0) {
    log_msg("警告：多出的样品（将忽略）: ", paste(c(extra_n, extra_t, extra_t6), collapse = ", "))
    sample_info <- sample_info[sample_info$sample %in% ALL_SAMPLES, , drop = FALSE]
  }

  t6_in_t <- intersect(sample_info$sample[sample_info$group == "T"], T6_SAMPLES)
  if (length(t6_in_t) > 0) stop("样品分组错误：T6 组样品被分进了 T")
  if ("T6" %in% sample_info$sample[sample_info$group == "T6"]) {
    stop("样品分组错误：T 组样品 T6 被分进了 T6 组")
  }

  protein_id <- vapply(raw[[id_col]], first_accession, character(1), USE.NAMES = FALSE)
  gene <- if (!is.na(gene_col)) {
    vapply(raw[[gene_col]], pick_official_symbol, character(1), USE.NAMES = FALSE)
  } else {
    rep(NA_character_, nrow(raw))
  }
  gene[is.na(gene) | !nzchar(gene)] <- protein_id[is.na(gene) | !nzchar(gene)]

  mat <- as.matrix(raw[, sample_info$column, drop = FALSE])
  storage.mode(mat) <- "double"
  colnames(mat) <- sample_info$sample
  mat[!is.finite(mat) | mat <= 0] <- NA_real_

  keep <- !is.na(protein_id) & nzchar(protein_id)
  if (any(duplicated(protein_id[keep]))) {
    log_msg("存在重复 Protein.Group，按行均值保留较高的一条")
  }
  row_mean <- rowMeans(mat, na.rm = TRUE)
  ord <- order(keep, row_mean, decreasing = TRUE, na.last = TRUE)
  mat <- mat[ord, , drop = FALSE]
  protein_id <- protein_id[ord]
  gene <- gene[ord]
  keep <- keep[ord]
  mat <- mat[keep, , drop = FALSE]
  protein_id <- protein_id[keep]
  gene <- gene[keep]
  dup <- duplicated(protein_id)
  if (any(dup)) {
    mat <- mat[!dup, , drop = FALSE]
    protein_id <- protein_id[!dup]
    gene <- gene[!dup]
  }
  rownames(mat) <- protein_id

  list(mat = mat, gene = gene, protein_id = protein_id, sample_info = sample_info)
}

# -----------------------------------------------------------------------------
# 4. 过滤 + log2 + 中位数标准化（三组一起做一次）
# -----------------------------------------------------------------------------
preprocess_intensity <- function(mat, sample_info) {
  det_frac <- rowMeans(is.finite(mat))
  keep <- det_frac >= MIN_DETECTED_FRAC
  log_msg(
    "过滤前蛋白 ", nrow(mat), "；检测率 >= ", MIN_DETECTED_FRAC,
    " 后保留 ", sum(keep)
  )
  if (sum(keep) < 10) stop("过滤后蛋白太少，请检查矩阵或放宽缺失阈值")
  mat <- mat[keep, , drop = FALSE]

  log_mat <- log2(mat)
  col_med <- apply(log_mat, 2, stats::median, na.rm = TRUE)
  if (any(!is.finite(col_med))) {
    bad <- names(col_med)[!is.finite(col_med)]
    stop("这些样品在过滤后全是缺失，未丢弃整组，请检查定量: ", paste(bad, collapse = ", "))
  }
  log_mat <- sweep(log_mat, 2, col_med, "-")
  log_msg("已做 log2 + 中位数标准化；样品中位数已写出日志")
  utils::write.csv(
    data.frame(sample = names(col_med), group = sample_info$group[match(names(col_med), sample_info$sample)],
               log2_median_before_norm = as.numeric(col_med)),
    file.path(log_dir, "sample_log2_medians.csv"),
    row.names = FALSE
  )
  list(log_mat = log_mat, keep_index = which(keep))
}

# -----------------------------------------------------------------------------
# 5. limma：A vs T，正 logFC = A 更高
# -----------------------------------------------------------------------------
run_limma <- function(log_mat, sample_info, test_group, comp_name) {
  ref <- "T"
  keep_s <- sample_info$sample[sample_info$group %in% c(test_group, ref)]
  if (test_group == "N") keep_s <- intersect(keep_s, c(N_SAMPLES, T_SAMPLES))
  if (test_group == "T6") keep_s <- intersect(keep_s, c(T6_SAMPLES, T_SAMPLES))
  if ("T6" %in% keep_s && test_group == "T6") {
    # T6 样品名属于 T 组；T6 组样品是 T6-1…。上面 already intersected.
  }
  if (any(keep_s %in% T6_SAMPLES) && test_group == "N") {
    keep_s <- setdiff(keep_s, T6_SAMPLES)
  }
  if (test_group == "T6" && any(keep_s %in% N_SAMPLES)) {
    keep_s <- setdiff(keep_s, N_SAMPLES)
  }

  grp <- droplevels(factor(
    sample_info$group[match(keep_s, sample_info$sample)],
    levels = c(ref, test_group)
  ))
  n_test <- sum(grp == test_group)
  n_ref <- sum(grp == ref)
  log_msg(comp_name, " 样品: ", paste(keep_s, collapse = ", "))
  log_msg(comp_name, " 设计: ", test_group, " n=", n_test, " vs ", ref, " n=", n_ref)
  if (n_test < 2 || n_ref < 2) stop(comp_name, " 每组至少需要 2 个样品")

  sub <- log_mat[, keep_s, drop = FALSE]
  ok <- rowSums(is.finite(sub[, grp == test_group, drop = FALSE])) >= 2 &
    rowSums(is.finite(sub[, grp == ref, drop = FALSE])) >= 2
  log_msg(comp_name, " 两组各至少 2 个非缺失后保留 ", sum(ok), " / ", length(ok))
  sub <- sub[ok, , drop = FALSE]

  design <- stats::model.matrix(~ grp)
  colnames(design) <- c("Intercept", test_group)
  fit <- limma::lmFit(sub, design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit, coef = test_group, number = Inf, sort.by = "none")
  tt$protein_id <- rownames(tt)
  tt$log2FC <- tt$logFC
  tt$pvalue <- tt$P.Value
  tt$padj <- tt$adj.P.Val
  tt$FC <- 2^tt$log2FC
  tt
}

# -----------------------------------------------------------------------------
# 6. 绘图与 ORA
# -----------------------------------------------------------------------------
save_gg <- function(plot, path_stub, width = 8, height = 6) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".pdf"), plot, width = width, height = height),
           error = function(e) log_msg("pdf ggsave failed: ", e$message))
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".png"), plot, width = width, height = height, dpi = 300),
           error = function(e) log_msg("png ggsave failed: ", e$message))
}

plot_volcano <- function(de, title, outfile) {
  df <- de
  df$y <- -log10(pmax(df$pvalue, 1e-300))
  lfc <- log2(FC_CUTOFF)
  df$set <- "NS"
  df$set[df$pvalue < P_CUTOFF & df$log2FC <= -lfc] <- "Down"
  df$set[df$pvalue < P_CUTOFF & df$log2FC >= lfc] <- "Up"
  df$set <- factor(df$set, levels = c("NS", "Down", "Up"))
  lab_n <- min(12, sum(df$set == "Up", na.rm = TRUE))
  up <- df[df$set == "Up", , drop = FALSE]
  up <- up[order(up$pvalue, -up$log2FC), , drop = FALSE]
  df$label <- NA_character_
  if (nrow(up) > 0 && lab_n > 0) {
    df$label[match(utils::head(up$gene, lab_n), df$gene)] <- utils::head(up$gene, lab_n)
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = y, color = set)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.4) +
    ggplot2::scale_color_manual(values = c(NS = "grey70", Down = "#4C78A8", Up = "#D62828")) +
    ggplot2::geom_vline(xintercept = c(-lfc, lfc), linetype = 2, color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(P_CUTOFF), linetype = 2, color = "grey40") +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, x = "log2 Fold Change", y = "-log10(p value)", color = NULL)
  save_gg(p, outfile, width = 8, height = 6)
}

map_ids_orgdb <- function(keys, keytype, column) {
  keys <- unique(keys[!is.na(keys) & nzchar(keys)])
  if (length(keys) == 0) return(setNames(character(), character()))
  mapped <- tryCatch(
    AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = keys, column = column, keytype = keytype, multiVals = "first"
    ),
    error = function(e) setNames(rep(NA_character_, length(keys)), keys)
  )
  mapped <- mapped[!is.na(mapped) & nzchar(mapped)]
  mapped
}

map_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (length(symbols) == 0) return(data.frame(gene = character(), entrez = character()))
  if (isTRUE(has_clusterProfiler)) {
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
    if (nrow(m) > 0) {
      m <- m[!duplicated(m[[1]]), ]
      return(data.frame(gene = m[[1]], entrez = as.character(m[[2]]), stringsAsFactors = FALSE))
    }
  }
  sym <- map_ids_orgdb(symbols, "SYMBOL", "ENTREZID")
  left <- setdiff(symbols, names(sym))
  uni <- if (length(left) > 0) map_ids_orgdb(left, "UNIPROT", "ENTREZID") else character()
  all_map <- c(sym, uni)
  if (length(all_map) == 0) return(data.frame(gene = character(), entrez = character()))
  data.frame(gene = names(all_map), entrez = as.character(unname(all_map)), stringsAsFactors = FALSE)
}

go_term_name <- function(ids) {
  ids <- as.character(ids)
  out <- setNames(ids, ids)
  if (isTRUE(has_godb)) {
    got <- tryCatch(AnnotationDbi::Term(GO.db::GOTERM[ids[ids %in% names(as.list(GO.db::GOTERM))]]), error = function(e) NULL)
    if (!is.null(got)) out[names(got)] <- unname(got)
  }
  out
}

kegg_pathway_names <- function(ids) {
  ids <- unique(as.character(ids))
  out <- setNames(ids, ids)
  cache <- file.path(log_dir, "kegg_hsa_pathway_names.tsv")
  tab <- NULL
  if (file.exists(cache)) {
    tab <- tryCatch(utils::read.delim(cache, header = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  }
  if (is.null(tab)) {
    tmp <- tempfile()
    ok <- tryCatch({
      utils::download.file("https://rest.kegg.jp/list/pathway/hsa", tmp, quiet = TRUE, mode = "wb")
      TRUE
    }, error = function(e) FALSE)
    if (ok && file.exists(tmp) && file.info(tmp)$size > 0) {
      tab <- tryCatch(utils::read.delim(tmp, header = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(tab)) utils::write.table(tab, cache, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }
  if (!is.null(tab) && ncol(tab) >= 2) {
    kid <- sub("^path:", "", tab[[1]])
    kid <- sub("^hsa", "", kid)
    names_k <- setNames(as.character(tab[[2]]), kid)
    hit <- intersect(ids, names(names_k))
    out[hit] <- unname(names_k[hit])
    hit2 <- intersect(ids, paste0("hsa", names(names_k)))
    out[hit2] <- unname(names_k[sub("^hsa", "", hit2)])
  }
  out
}

ora_hyper <- function(query, universe, term2gene, term2name = NULL, min_term = 2) {
  query <- unique(intersect(as.character(query), as.character(universe)))
  universe <- unique(as.character(universe))
  term2gene <- unique(term2gene[term2gene$gene %in% universe, , drop = FALSE])
  if (nrow(term2gene) == 0 || length(query) == 0) return(data.frame())
  split_g <- split(term2gene$gene, term2gene$term)
  n_u <- length(universe)
  n_q <- length(query)
  rows <- lapply(names(split_g), function(tm) {
    in_term <- unique(split_g[[tm]])
    hit <- intersect(query, in_term)
    k <- length(hit)
    m <- length(in_term)
    if (k == 0 || m < min_term) return(NULL)
    p <- stats::phyper(k - 1, m, n_u - m, n_q, lower.tail = FALSE)
    desc <- tm
    if (!is.null(term2name) && tm %in% names(term2name)) desc <- unname(term2name[[tm]])
    data.frame(
      ID = tm,
      Description = desc,
      GeneRatio = paste0(k, "/", n_q),
      BgRatio = paste0(m, "/", n_u),
      pvalue = p,
      Count = k,
      geneID = paste(hit, collapse = "/"),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(data.frame())
  df <- do.call(rbind, rows)
  df$p.adjust <- stats::p.adjust(df$pvalue, method = "BH")
  df$qvalue <- df$p.adjust
  df[order(df$pvalue, -df$Count), , drop = FALSE]
}

build_go_term2gene <- function(universe, ont) {
  tab <- tryCatch(
    suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db, keys = universe, keytype = "ENTREZID",
      columns = c("GOALL", "ONTOLOGYALL")
    )),
    error = function(e) {
      suppressMessages(AnnotationDbi::select(
        org.Hs.eg.db, keys = universe, keytype = "ENTREZID",
        columns = c("GO", "ONTOLOGY")
      ))
    }
  )
  if (is.null(tab) || nrow(tab) == 0) return(data.frame(term = character(), gene = character()))
  ont_col <- if ("ONTOLOGYALL" %in% names(tab)) "ONTOLOGYALL" else "ONTOLOGY"
  go_col <- if ("GOALL" %in% names(tab)) "GOALL" else "GO"
  tab <- tab[!is.na(tab[[go_col]]) & tab[[ont_col]] == ont, , drop = FALSE]
  data.frame(term = as.character(tab[[go_col]]), gene = as.character(tab$ENTREZID), stringsAsFactors = FALSE)
}

build_kegg_term2gene <- function(universe) {
  tab <- tryCatch(
    suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db, keys = universe, keytype = "ENTREZID", columns = "PATH"
    )),
    error = function(e) NULL
  )
  if (!is.null(tab) && "PATH" %in% names(tab)) {
    tab <- tab[!is.na(tab$PATH) & nzchar(tab$PATH), , drop = FALSE]
    if (nrow(tab) > 0) {
      return(data.frame(
        term = paste0("hsa", tab$PATH),
        gene = as.character(tab$ENTREZID),
        stringsAsFactors = FALSE
      ))
    }
  }
  tmp <- tempfile()
  ok <- tryCatch({
    utils::download.file("https://rest.kegg.jp/link/pathway/hsa", tmp, quiet = TRUE, mode = "wb")
    TRUE
  }, error = function(e) FALSE)
  if (!ok || !file.exists(tmp) || file.info(tmp)$size == 0) {
    return(data.frame(term = character(), gene = character()))
  }
  raw <- utils::read.delim(tmp, header = FALSE, stringsAsFactors = FALSE)
  term <- sub("^path:", "", raw[[1]])
  gene <- sub("^hsa:", "", raw[[2]])
  keep <- gene %in% universe
  data.frame(term = term[keep], gene = gene[keep], stringsAsFactors = FALSE)
}

note_empty <- function(stub, msg) {
  dir.create(dirname(stub), recursive = TRUE, showWarnings = FALSE)
  writeLines(msg, paste0(stub, "_EMPTY.txt"))
}

plot_ora_df <- function(df, stub, title) {
  dir.create(dirname(stub), recursive = TRUE, showWarnings = FALSE)
  if (is.null(df) || nrow(df) == 0) {
    note_empty(stub, "no enrichment terms")
    return(invisible(NULL))
  }
  utils::write.csv(df, paste0(stub, ".csv"), row.names = FALSE)
  show <- utils::head(df, 15)
  show$Description <- factor(show$Description, levels = rev(unique(show$Description)))
  show$neglogp <- -log10(pmax(show$p.adjust, 1e-300))
  gr <- as.numeric(sub("/.*", "", show$GeneRatio)) / as.numeric(sub(".*/", "", show$GeneRatio))
  show$GeneRatioNum <- gr
  bar <- ggplot2::ggplot(show, ggplot2::aes(x = Description, y = neglogp)) +
    ggplot2::geom_col(fill = "#D62828", width = 0.75) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, x = NULL, y = "-log10(adjusted p)")
  save_gg(bar, paste0(stub, "_barplot"), width = 9, height = max(5, min(10, 0.35 * nrow(show) + 3)))
  dot <- ggplot2::ggplot(show, ggplot2::aes(x = GeneRatioNum, y = Description, size = Count, color = p.adjust)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "#D62828", high = "#4C78A8") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, x = "GeneRatio", y = NULL)
  save_gg(dot, paste0(stub, "_dotplot"), width = 9, height = max(5, min(10, 0.35 * nrow(show) + 3)))
}

plot_ora_object <- function(x, stub, title) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) {
    note_empty(stub, "no enrichment terms")
    return(invisible(NULL))
  }
  df <- as.data.frame(x)
  plot_ora_df(df, stub, title)
  if (isTRUE(has_enrichplot)) {
    nshow <- min(15, nrow(df))
    tryCatch({
      save_gg(enrichplot::dotplot(x, showCategory = nshow) + ggplot2::ggtitle(title),
              paste0(stub, "_enrichplot_dotplot"), 9, 7)
    }, error = function(e) log_msg("enrichplot dotplot failed: ", e$message))
  }
}

run_local_go_kegg <- function(query, universe, go_dir, kg_dir, label) {
  log_msg(label, " 使用 org.Hs.eg.db 本地超几何检验做 GO / KEGG（不上 clusterProfiler 也出图）")
  for (ont in c("BP", "MF", "CC")) {
    t2g <- build_go_term2gene(universe, ont)
    names_go <- go_term_name(unique(t2g$term))
    df <- ora_hyper(query, universe, t2g, names_go)
    plot_ora_df(df, file.path(go_dir, paste0("ORA_GO_", ont)),
                paste(label, "| ORA GO", ont, "| upregulated"))
  }
  t2k <- build_kegg_term2gene(universe)
  names_k <- kegg_pathway_names(unique(t2k$term))
  dfk <- ora_hyper(query, universe, t2k, names_k)
  plot_ora_df(dfk, file.path(kg_dir, "ORA_KEGG"),
              paste(label, "| ORA KEGG | upregulated"))
}

run_up_ora <- function(up_genes, universe_genes, outdir, label) {
  go_dir <- file.path(outdir, "GO")
  kg_dir <- file.path(outdir, "KEGG")
  dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(kg_dir, recursive = TRUE, showWarnings = FALSE)

  mp_up <- map_to_entrez(up_genes)
  mp_bg <- map_to_entrez(universe_genes)
  query <- unique(mp_up$entrez)
  universe <- unique(mp_bg$entrez)
  log_msg(label, " 上调蛋白 ", length(up_genes), " → Entrez ", length(query),
          "；背景 ", length(universe_genes), " → Entrez ", length(universe))
  if (length(query) < 3) {
    writeLines(paste("mapped_entrez", length(query)), file.path(outdir, "ORA_too_few_genes.txt"))
    log_msg(label, " 上调映射基因 < 3，仍写出空表说明，但会继续尝试富集")
  }

  used_cp <- FALSE
  if (isTRUE(has_clusterProfiler) && length(query) >= 3) {
    for (ont in c("BP", "MF", "CC")) {
      ego <- tryCatch(
        clusterProfiler::enrichGO(
          gene = query, universe = universe, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
          ont = ont, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
        ),
        error = function(e) {
          log_msg("enrichGO ", ont, " failed: ", e$message)
          NULL
        }
      )
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        plot_ora_object(ego, file.path(go_dir, paste0("ORA_GO_", ont)),
                        paste(label, "| ORA GO", ont, "| upregulated"))
        used_cp <- TRUE
      }
    }
    ek <- tryCatch(
      clusterProfiler::enrichKEGG(
        gene = query, universe = universe, organism = "hsa",
        pvalueCutoff = 1, qvalueCutoff = 1
      ),
      error = function(e) {
        log_msg("enrichKEGG failed: ", e$message)
        NULL
      }
    )
    if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) {
      ek <- tryCatch(
        clusterProfiler::setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
        error = function(e) ek
      )
      plot_ora_object(ek, file.path(kg_dir, "ORA_KEGG"),
                      paste(label, "| ORA KEGG | upregulated"))
      used_cp <- TRUE
    }
  }

  go_ok <- length(list.files(go_dir, pattern = "^ORA_GO_.*\\.csv$")) >= 3
  kg_ok <- length(list.files(kg_dir, pattern = "^ORA_KEGG\\.csv$")) >= 1
  if (!go_ok || !kg_ok) {
    run_local_go_kegg(query, universe, go_dir, kg_dir, label)
  } else if (isTRUE(used_cp)) {
    log_msg(label, " GO / KEGG 已由 clusterProfiler 写出")
  }
}

# -----------------------------------------------------------------------------
# 7. 单个比较出表出图
# -----------------------------------------------------------------------------
analyze_comparison <- function(comp_name, test_group, log_mat, gene_map, sample_info) {
  outdir <- file.path(result_dir, comp_name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  log_msg("==== ", comp_name, " ====")

  tt <- run_limma(log_mat, sample_info, test_group, comp_name)
  tt$gene <- gene_map[tt$protein_id]
  tt$gene[is.na(tt$gene) | !nzchar(tt$gene)] <- tt$protein_id
  tt$significant <- tt$pvalue < P_CUTOFF
  tt$up <- tt$significant & tt$log2FC > 0 & tt$FC >= FC_CUTOFF
  tt$down <- tt$significant & tt$log2FC < 0 & tt$FC <= 1 / FC_CUTOFF
  tt <- tt[order(tt$pvalue, -abs(tt$log2FC)), , drop = FALSE]

  de_out <- tt[, c("protein_id", "gene", "log2FC", "FC", "pvalue", "padj",
                   "AveExpr", "t", "B", "significant", "up", "down")]
  utils::write.csv(de_out, file.path(outdir, paste0(comp_name, "_DE_all.csv")), row.names = FALSE)
  up <- tt[tt$up, , drop = FALSE]
  utils::write.csv(
    up[, c("protein_id", "gene", "log2FC", "FC", "pvalue", "padj")],
    file.path(outdir, paste0(comp_name, "_DE_up_p", P_CUTOFF, "_FC", FC_CUTOFF, ".csv")),
    row.names = FALSE
  )
  log_msg(comp_name, " 全表 ", nrow(tt), "；上调 ", nrow(up),
          " (p < ", P_CUTOFF, " 且 FC >= ", FC_CUTOFF, ")")

  plot_volcano(
    tt,
    paste0(comp_name, " | volcano | p < ", P_CUTOFF, ", FC >= ", FC_CUTOFF),
    file.path(outdir, paste0(comp_name, "_volcano"))
  )
  run_up_ora(up$gene, tt$gene, outdir, paste(comp_name, "upregulated"))
  invisible(tt)
}

# -----------------------------------------------------------------------------
# 8. 主流程（PROTEIN_NT_T6_FUNCTIONS_ONLY=true 时只加载函数）
# -----------------------------------------------------------------------------
if (isTRUE(as.logical(Sys.getenv("PROTEIN_NT_T6_FUNCTIONS_ONLY", "false")))) {
  log_msg("仅加载函数，不跑主流程")
} else {
pg_path <- find_pg_matrix(project_dir)
dat <- read_pg_matrix(pg_path)
utils::write.csv(dat$sample_info, file.path(log_dir, "sample_map.csv"), row.names = FALSE)
log_msg(
  "样品计数 N=", sum(dat$sample_info$group == "N"),
  " T=", sum(dat$sample_info$group == "T"),
  " T6=", sum(dat$sample_info$group == "T6")
)

prep <- preprocess_intensity(dat$mat, dat$sample_info)
gene_map <- setNames(dat$gene, dat$protein_id)
gene_map <- gene_map[rownames(prep$log_mat)]

analyze_comparison("N_vs_T", "N", prep$log_mat, gene_map, dat$sample_info)
analyze_comparison("T6_vs_T", "T6", prep$log_mat, gene_map, dat$sample_info)

log_msg("完成。结果在: ", result_dir)
log_msg("  results/N_vs_T/")
log_msg("  results/T6_vs_T/")
}
