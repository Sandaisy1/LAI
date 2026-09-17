#!/usr/bin/env Rscript
# =============================================================================
# E:\R\Nerve RNA  —— 乳腺癌转移 bulk（GSE175692）独立分析
#
# 问题1：肿瘤活检里哪些高表达基因可能促进神经浸润
# 问题2：神经浸润更倾向哪个继发部位（肺 / 脑 / 骨 …）
#
# 数据目录默认：E:/R/Nerve RNA
# 需要放入（二选一或都放）：
#   GSE175692_series_matrix.txt.gz   （推荐，含 organ 注释）
#   GSE175692_raw_data.txt.gz
# 若本地没有且能上网，会尝试 GEOquery 下载。
#
# 用法（Windows R / RStudio）：
#   setwd("E:/R/Nerve RNA")
#   source("Nerve_RNA_nerve_infiltration.R")
#
# 独立于 TG_RNAseq_*.R / ST / scRNA，不修改原流程。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")

need_cran <- c("ggplot2", "matrixStats")
opt_cran  <- c("writexl", "pheatmap", "preprocessCore")
need_bioc <- c("limma")
opt_bioc  <- c("GEOquery", "Biobase", "org.Hs.eg.db", "clusterProfiler",
               "enrichplot", "fgsea", "msigdbr")

install_if_missing <- function(pkgs, bioc = FALSE, required = FALSE) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) == 0) return(invisible(TRUE))
  # 优先写到可写的用户库，避免系统库无权限时直接失败
  writable <- Filter(function(p) {
    dir.exists(p) && file.access(p, 2) == 0
  }, .libPaths())
  if (length(writable) == 0) {
    user_lib <- path.expand(file.path("~", "R", "library"))
    dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
    .libPaths(c(user_lib, .libPaths()))
    writable <- user_lib
  }
  lib <- writable[1]
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      tryCatch(install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = lib),
               error = function(e) NULL)
    }
    tryCatch(BiocManager::install(miss, update = FALSE, ask = FALSE, lib = lib),
             error = function(e) message("Bioc install failed: ", e$message))
  } else {
    tryCatch(install.packages(miss, repos = "https://cloud.r-project.org", lib = lib),
             error = function(e) message("CRAN install failed: ", e$message))
  }
  still <- miss[!vapply(miss, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still) > 0 && required) stop("缺少 R 包: ", paste(still, collapse = ", "))
  if (length(still) > 0) message("可选包未装，相关步骤跳过: ", paste(still, collapse = ", "))
  invisible(TRUE)
}
install_if_missing(need_cran, FALSE, TRUE)
install_if_missing(opt_cran, FALSE, FALSE)
install_if_missing(need_bioc, TRUE, TRUE)
install_if_missing(opt_bioc, TRUE, FALSE)
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
for (p in c(need_cran, opt_cran, need_bioc, opt_bioc)) {
  if (has_pkg(p)) suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# -----------------------------------------------------------------------------
# 0. 路径：E:/R/Nerve RNA
# -----------------------------------------------------------------------------
resolve_nerve_rna_dir <- function() {
  env_dir <- Sys.getenv("NERVE_RNA_DIR", unset = "")
  cands <- unique(c(
    env_dir,
    "E:/R/Nerve RNA", "E:\\R\\Nerve RNA",
    "E:/R/Nerve_RNA", "E:\\R\\Nerve_RNA",
    getwd()
  ))
  looks <- function(d) {
    dir.exists(d) && (
      length(list.files(d, pattern = "GSE175692", ignore.case = TRUE)) > 0 ||
        length(list.files(d, pattern = "series_matrix|raw_data", ignore.case = TRUE)) > 0 ||
        identical(normalizePath(d, winslash = "/", mustWork = FALSE),
                  normalizePath(getwd(), winslash = "/", mustWork = FALSE))
    )
  }
  for (d in cands[nzchar(cands)]) {
    if (looks(d) || dir.exists(d)) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  stop("找不到数据目录。请创建 E:/R/Nerve RNA 并放入 GSE175692 文件，或 setwd 到该目录。")
}

data_dir   <- resolve_nerve_rna_dir()
result_dir <- file.path(data_dir, "results")
log_dir    <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("Nerve_RNA_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
log_msg("Data directory: ", data_dir)

p_cutoff <- 0.01
fc_cutoffs <- c("FC_1.25" = 1.25, "FC_1.5" = 1.5, "FC_2" = 2)

# GSE175692 NanoString 面板上实际存在的神经亲和 / 侵袭相关基因
# （完整 NGF/ARTN/SEMA3 等多数不在 771 基因里，见日志）
ligand_priority <- c(
  "BDNF", "NGFR", "NTRK2", "CXCL12", "VEGFA", "VEGFD", "TGFB1", "TGFB2", "TGFB3",
  "WNT5A", "WNT5B", "WNT11", "NCAM1", "HGF", "MET", "ERBB2", "ERBB4",
  "CCL2", "EFNA5", "ROBO4", "CDH2", "VIM", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2",
  "MMP9", "MMP14", "ITGB1", "ITGAV"
)
stroma_genes <- c("SOX10", "S100A7", "S100A14")  # 面板几乎无 GFAP/MBP/MPZ
nerve_term_pat <- paste(
  "axon", "neur", "nerve", "schwann", "synap", "neurotroph", "semaphorin",
  "ephrin", "glia", "wnt", "emt", "invas", "migrat",
  sep = "|"
)

# -----------------------------------------------------------------------------
# 1. 工具函数
# -----------------------------------------------------------------------------
norm_organ <- function(x) {
  x <- tolower(trimws(as.character(x)))
  if (grepl("brain|cns|cerebr", x)) return("brain")
  if (grepl("lung|pulmon", x)) return("lung")
  if (grepl("bone|osseo|marrow", x)) return("bone")
  if (grepl("liver|hepat", x)) return("liver")
  if (grepl("skin|derm|cutan", x)) return("skin")
  if (grepl("lymph|node|axilla", x)) return("lymph_node")
  if (grepl("breast", x)) return("breast")
  if (grepl("pleura", x)) return("pleura")
  if (grepl("ovar", x)) return("ovary")
  if (grepl("periton", x)) return("peritoneum")
  if (grepl("muscle", x)) return("muscle")
  "other"
}

present <- function(genes, universe) intersect(unique(genes), universe)

zmean_score <- function(mat, genes) {
  g <- present(genes, rownames(mat))
  if (length(g) == 0) return(rep(NA_real_, ncol(mat)))
  if (length(g) == 1) {
    z <- scale(as.numeric(mat[g, ]))[, 1]
    z[is.na(z)] <- 0
    return(as.numeric(z))
  }
  z <- t(scale(t(mat[g, , drop = FALSE])))
  z[is.na(z)] <- 0
  colMeans(z)
}

quantile_norm <- function(mat) {
  if (has_pkg("preprocessCore")) {
    out <- preprocessCore::normalize.quantiles(as.matrix(mat))
    dimnames(out) <- dimnames(mat)
    return(out)
  }
  apply(mat, 2, function(v) {
    r <- rank(v, ties.method = "average", na.last = "keep")
    stats::qnorm((r - 0.5) / sum(is.finite(r)))
  })
}

empty_de <- function() {
  data.frame(gene = character(), log2FC = numeric(), AveExpr = numeric(),
             pvalue = numeric(), padj = numeric(), padj_BH = numeric(),
             stringsAsFactors = FALSE)
}

limma_near_far <- function(mat, group, label = "") {
  group <- factor(as.character(group), levels = c("far", "near"))
  if (nlevels(droplevels(group)) < 2 || ncol(mat) < 4) return(empty_de())
  design <- stats::model.matrix(~ group)
  fit <- tryCatch(
    limma::eBayes(limma::lmFit(mat, design), trend = TRUE, robust = TRUE),
    error = function(e) { log_msg("limma failed ", label, ": ", e$message); NULL }
  )
  if (is.null(fit)) return(empty_de())
  coefn <- grep("groupnear", colnames(design), value = TRUE)
  tt <- limma::topTable(fit, coef = coefn, number = Inf, sort.by = "none")
  data.frame(
    gene = rownames(tt), log2FC = tt$logFC, AveExpr = tt$AveExpr,
    pvalue = tt$P.Value, padj = tt$P.Value, padj_BH = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
}

save_gg <- function(p, stub, w = 7, h = 5) {
  dir.create(dirname(stub), recursive = TRUE, showWarnings = FALSE)
  if (!has_pkg("ggplot2") || is.null(p)) return(invisible(NULL))
  ggplot2::ggsave(paste0(stub, ".pdf"), p, width = w, height = h)
  ggplot2::ggsave(paste0(stub, ".png"), p, width = w, height = h, dpi = 140)
}

basic_volcano <- function(de, title, stub, fc_line = 1.25) {
  if (nrow(de) == 0) return(invisible(NULL))
  df <- de
  df$y <- -log10(pmax(df$pvalue, 1e-300))
  df$col <- "ns"
  df$col[!is.na(df$pvalue) & df$pvalue < p_cutoff & df$log2FC >= log2(fc_line)] <- "up"
  p <- ggplot2::ggplot(df, ggplot2::aes(log2FC, y, color = col)) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::scale_color_manual(values = c(ns = "grey70", up = "#D62828")) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, y = "-log10(p)", x = "log2FC")
  save_gg(p, stub, 7, 6)
}

# -----------------------------------------------------------------------------
# 2. 读入 GSE175692（本地优先）
# -----------------------------------------------------------------------------
find_file <- function(patterns) {
  hits <- character()
  for (pat in patterns) {
    hits <- c(hits, list.files(data_dir, pattern = pat, recursive = TRUE,
                               full.names = TRUE, ignore.case = TRUE))
  }
  unique(hits[file.exists(hits)])
}

parse_series_matrix <- function(path) {
  log_msg("Reading series matrix: ", path)
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  char_lines <- lines[grepl("^!Sample_characteristics_ch1", lines)]
  organ_line <- char_lines[grepl("organ:", char_lines, ignore.case = TRUE)][1]
  geo_line <- lines[grepl("^!Sample_geo_accession", lines)][1]
  title_line <- lines[grepl("^!Sample_title", lines)][1]
  if (is.na(organ_line) || !nzchar(organ_line)) stop("series matrix 里没有 organ 注释行")

  split_quoted <- function(line) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    parts <- parts[-1]
    gsub('^"|"$', "", parts)
  }
  organs_raw <- split_quoted(organ_line)
  organs_raw <- sub("^organ:\\s*", "", organs_raw, ignore.case = TRUE)
  geos <- if (!is.na(geo_line)) split_quoted(geo_line) else paste0("S", seq_along(organs_raw))
  titles <- if (!is.na(title_line)) split_quoted(title_line) else geos

  begin <- which(grepl("series_matrix_table_begin", lines))[1]
  end <- which(grepl("series_matrix_table_end", lines))[1]
  if (is.na(begin) || is.na(end) || end <= begin + 1) stop("找不到 series matrix 表达表")
  tab <- utils::read.delim(
    textConnection(lines[(begin + 1):(end - 1)]),
    check.names = FALSE, stringsAsFactors = FALSE, quote = "\""
  )
  gene <- as.character(tab[[1]])
  mat <- as.matrix(tab[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- gene
  # 列名常是 GSM；与 title 对齐
  if (ncol(mat) == length(geos)) colnames(mat) <- geos
  si <- data.frame(
    sample = colnames(mat),
    gsm = geos[seq_len(ncol(mat))],
    title = titles[seq_len(ncol(mat))],
    organ_raw = organs_raw[seq_len(ncol(mat))],
    organ = vapply(organs_raw[seq_len(ncol(mat))], norm_organ, character(1)),
    stringsAsFactors = FALSE
  )
  rownames(si) <- si$sample
  list(mat = mat, si = si, source = basename(path), is_log = TRUE)
}

parse_raw_data <- function(path, si_ref = NULL) {
  log_msg("Reading raw data: ", path)
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  # 跳过开头注释，定位 Gene Name 行
  lines <- readLines(con, warn = FALSE)
  hdr <- which(grepl("^Gene Name\\t|^Gene Name,", lines))[1]
  if (is.na(hdr)) stop("raw_data 找不到 Gene Name 表头")
  tab <- utils::read.delim(
    textConnection(lines[hdr:length(lines)]),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  gene <- as.character(tab[["Gene Name"]])
  # 跳过阳性/阴性对照探针
  cls <- if ("Class Name" %in% names(tab)) as.character(tab[["Class Name"]]) else rep("Endogenous", nrow(tab))
  keep_row <- is.na(cls) | grepl("Endogenous", cls, ignore.case = TRUE) | cls == ""
  drop_cols <- intersect(c("Gene Name", "Accession #", "Class Name"), names(tab))
  mat <- as.matrix(tab[, setdiff(names(tab), drop_cols), drop = FALSE])
  mat[mat %in% c("null", "NULL", "", "NA")] <- NA
  storage.mode(mat) <- "double"
  mat[is.na(mat)] <- 0
  mat <- mat[keep_row, , drop = FALSE]
  gene <- gene[keep_row]
  rownames(mat) <- make.unique(gene)
  # 样本名是 TGTT_xxxx；用 si_ref$title 映射 organ
  si <- data.frame(sample = colnames(mat), title = colnames(mat),
                   organ = "other", organ_raw = NA_character_, stringsAsFactors = FALSE)
  if (!is.null(si_ref) && "title" %in% names(si_ref)) {
    m <- match(si$title, si_ref$title)
    hit <- !is.na(m)
    si$organ[hit] <- si_ref$organ[m[hit]]
    si$organ_raw[hit] <- si_ref$organ_raw[m[hit]]
    si$gsm <- NA_character_
    si$gsm[hit] <- si_ref$gsm[m[hit]]
  }
  rownames(si) <- si$sample
  list(mat = mat, si = si, source = basename(path), is_log = FALSE)
}

try_geoquery <- function() {
  if (!has_pkg("GEOquery")) return(NULL)
  dest <- file.path(data_dir, "GEO")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  log_msg("本地文件不足，尝试 GEOquery::getGEO(GSE175692)")
  g <- tryCatch(GEOquery::getGEO("GSE175692", destdir = dest, getGPL = TRUE),
                error = function(e) { log_msg(e$message); NULL })
  if (is.null(g) || length(g) == 0) return(NULL)
  eset <- g[[1]]
  mat <- as.matrix(Biobase::exprs(eset))
  pd <- as.data.frame(Biobase::pData(eset))
  blob <- apply(pd, 1, function(r) paste(r, collapse = " | "))
  organ_raw <- vapply(blob, function(x) {
    m <- regmatches(x, regexpr("organ:\\s*[^|;]+", x, ignore.case = TRUE))
    if (length(m) == 0) return(NA_character_)
    sub("(?i)organ:\\s*", "", m)
  }, character(1))
  si <- data.frame(
    sample = colnames(mat),
    gsm = colnames(mat),
    title = if ("title" %in% names(pd)) as.character(pd$title) else colnames(mat),
    organ_raw = organ_raw,
    organ = vapply(organ_raw, norm_organ, character(1)),
    stringsAsFactors = FALSE
  )
  rownames(si) <- si$sample
  list(mat = mat, si = si, source = "GEOquery_GSE175692", is_log = TRUE)
}

matrix_files <- find_file(c("GSE175692_series_matrix\\.txt(\\.gz)?$", "series_matrix\\.txt(\\.gz)?$"))
raw_files    <- find_file(c("GSE175692_raw_data\\.txt(\\.gz)?$", "raw_data\\.txt(\\.gz)?$"))

ds_meta <- NULL
ds_expr <- NULL
if (length(matrix_files) > 0) {
  ds_meta <- tryCatch(parse_series_matrix(matrix_files[1]), error = function(e) {
    log_msg("series matrix 解析失败: ", e$message); NULL
  })
}
if (length(raw_files) > 0) {
  ds_expr <- tryCatch(parse_raw_data(raw_files[1], si_ref = if (!is.null(ds_meta)) ds_meta$si else NULL),
                      error = function(e) { log_msg("raw 解析失败: ", e$message); NULL })
}

# 优先：raw 计数 + matrix 的 organ；否则只用 matrix
if (!is.null(ds_expr) && !is.null(ds_meta)) {
  # 用 title 对齐 organ
  m <- match(ds_expr$si$title, ds_meta$si$title)
  if (all(is.na(m))) m <- match(ds_expr$si$sample, ds_meta$si$title)
  hit <- !is.na(m)
  ds_expr$si$organ[hit] <- ds_meta$si$organ[m[hit]]
  ds_expr$si$organ_raw[hit] <- ds_meta$si$organ_raw[m[hit]]
  # 分析用 log2(count+1) 后再分位数标准化
  mat_use <- log2(pmax(ds_expr$mat, 0) + 1)
  si_use <- ds_expr$si
  src_lab <- paste0(ds_expr$source, "+", ds_meta$source)
} else if (!is.null(ds_meta)) {
  mat_use <- ds_meta$mat
  si_use <- ds_meta$si
  src_lab <- ds_meta$source
} else if (!is.null(ds_expr) && any(si_use_ok <- ds_expr$si$organ != "other")) {
  mat_use <- log2(pmax(ds_expr$mat, 0) + 1)
  si_use <- ds_expr$si
  src_lab <- ds_expr$source
} else {
  gq <- try_geoquery()
  if (is.null(gq)) {
    stop(
      "E:/R/Nerve RNA 下没有可读的 GSE175692。\n",
      "请下载并放入：\n",
      "  GSE175692_series_matrix.txt.gz\n",
      "  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175692/matrix/\n",
      "可选：GSE175692_raw_data.txt.gz（suppl/）"
    )
  }
  mat_use <- gq$mat
  si_use <- gq$si
  src_lab <- gq$source
}

# 过滤 + 标准化
keep_gene <- rowSums(is.finite(mat_use)) >= max(10, floor(0.5 * ncol(mat_use)))
mat_use <- mat_use[keep_gene, , drop = FALSE]
mat_use[!is.finite(mat_use)] <- min(mat_use[is.finite(mat_use)], na.rm = TRUE)
mat_n <- quantile_norm(mat_use)
si_use <- si_use[match(colnames(mat_n), si_use$sample), , drop = FALSE]
rownames(si_use) <- si_use$sample

utils::write.csv(si_use, file.path(log_dir, "sample_organ_annotation.csv"), row.names = FALSE)
org_tab <- sort(table(si_use$organ), decreasing = TRUE)
log_msg("Source: ", src_lab)
log_msg("Genes: ", nrow(mat_n), "  Samples: ", ncol(mat_n))
log_msg("Organs: ", paste(sprintf("%s=%d", names(org_tab), as.integer(org_tab)), collapse = ", "))

lig_in <- present(ligand_priority, rownames(mat_n))
str_in <- present(stroma_genes, rownames(mat_n))
log_msg("神经亲和候选基因（面板命中 ", length(lig_in), "/", length(ligand_priority), "）: ",
        paste(lig_in, collapse = ", "))
log_msg("神经基质基因（面板命中 ", length(str_in), "/", length(stroma_genes), "）: ",
        paste(str_in, collapse = ", "))
if (length(lig_in) < 3) log_msg("WARNING: 面板神经亲和基因很少，问题1结果会偏窄")

# -----------------------------------------------------------------------------
# 3. 问题2：继发部位倾向（先算分数，后面写总表）
# -----------------------------------------------------------------------------
si_use$score_ligand <- zmean_score(mat_n, lig_in)
si_use$score_stroma <- zmean_score(mat_n, str_in)
ok <- is.finite(si_use$score_ligand)
si_use$score_ligand_resid <- NA_real_
if (sum(ok) >= 8 && any(is.finite(si_use$score_stroma))) {
  ok2 <- ok & is.finite(si_use$score_stroma)
  if (sum(ok2) >= 8) {
    fit <- stats::lm(score_ligand ~ score_stroma, data = si_use[ok2, ])
    si_use$score_ligand_resid[ok2] <- stats::resid(fit)
  }
} else if (sum(ok) >= 8) {
  # 面板几乎无神经基质基因时，残差 = 配体分数本身
  si_use$score_ligand_resid[ok] <- si_use$score_ligand[ok]
  log_msg("面板神经基质基因不足，部位排名直接用配体分数（并在表中注明）")
}

# -----------------------------------------------------------------------------
# 4. 问题1：各器官 vs 其余；脑/肺/骨两两比较（只看上调）
# -----------------------------------------------------------------------------
emit_one <- function(comp_name, de) {
  base <- file.path(result_dir, "Q1_genes_promote_nerve", comp_name)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  de <- de[!is.na(de$log2FC), , drop = FALSE]
  utils::write.csv(de, file.path(base, "DE_full.csv"), row.names = FALSE)
  writeLines(
    c("问题1：该比较中上调基因（先 p < 0.01，再 FC >= 1.25 / 1.5 / 2）",
      "只分析上调。关联 ≠ 已证明促进浸润。",
      "优先看候选表里的 BDNF/NGFR/NTRK2/CXCL12/VEGFA/TGFB1/WNT5A/NCAM1 等。"),
    file.path(base, "00_READ_ME.txt")
  )
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    od <- file.path(base, "FoldChange", nm)
    dir.create(od, recursive = TRUE, showWarnings = FALSE)
    keep <- !is.na(de$pvalue) & de$pvalue < p_cutoff & de$log2FC > 0 & (2^de$log2FC >= fc)
    sub <- de[keep, , drop = FALSE]
    if (nrow(sub) > 0) sub <- sub[order(sub$log2FC, decreasing = TRUE), , drop = FALSE]
    utils::write.csv(sub, file.path(od, paste0(nm, "_up_genes.csv")), row.names = FALSE)
    basic_volcano(de, paste(comp_name, nm), file.path(od, paste0(nm, "_volcano")), fc)
    # 神经相关基因子集
    focus <- sub[sub$gene %in% lig_in, , drop = FALSE]
    utils::write.csv(focus, file.path(od, paste0(nm, "_FOCUS_nerve_ligands.csv")), row.names = FALSE)
  }
  invisible(de)
}

de_store <- list()
focus_organs <- c("brain", "lung", "bone", "liver", "skin", "lymph_node", "breast", "pleura")
for (org in intersect(focus_organs, unique(si_use$organ))) {
  n1 <- sum(si_use$organ == org)
  n0 <- sum(si_use$organ != org)
  if (n1 < 3 || n0 < 3) {
    log_msg("跳过 ", org, " vs rest（n=", n1, "/", n0, "）")
    next
  }
  grp <- ifelse(si_use$organ == org, "near", "far")
  de <- limma_near_far(mat_n, grp, paste(org, "vs rest"))
  de_store[[paste0(org, "_vs_other")]] <- emit_one(paste0(org, "_vs_other_mets"), de)
}

pairs <- list(c("brain", "lung"), c("brain", "bone"), c("lung", "bone"),
              c("brain", "liver"), c("lung", "liver"), c("bone", "liver"))
for (pr in pairs) {
  a <- pr[1]; b <- pr[2]
  if (sum(si_use$organ == a) < 3 || sum(si_use$organ == b) < 3) next
  keep <- si_use$organ %in% c(a, b)
  grp <- ifelse(si_use$organ[keep] == a, "near", "far")
  de <- limma_near_far(mat_n[, keep, drop = FALSE], grp, paste(a, "vs", b))
  de_store[[paste0(a, "_vs_", b)]] <- emit_one(paste0(a, "_vs_", b), de)
}

# -----------------------------------------------------------------------------
# 5. 候选分子总表（问题1）
# -----------------------------------------------------------------------------
cand <- data.frame(
  gene = lig_in,
  on_panel = TRUE,
  mean_expr = rowMeans(mat_n[lig_in, , drop = FALSE]),
  stringsAsFactors = FALSE
)
# 各器官 vs 其余的 log2FC / p
for (org in c("brain", "lung", "bone")) {
  key <- paste0(org, "_vs_other")
  if (!key %in% names(de_store)) next
  d <- de_store[[key]]
  cand[[paste0(org, "_vs_rest_log2FC")]] <- d$log2FC[match(cand$gene, d$gene)]
  cand[[paste0(org, "_vs_rest_pvalue")]] <- d$pvalue[match(cand$gene, d$gene)]
}
fc_cols <- grep("_log2FC$", names(cand), value = TRUE)
pv_cols <- grep("_pvalue$", names(cand), value = TRUE)
cand$n_up_p01 <- 0L
if (length(fc_cols) > 0) {
  for (i in seq_along(fc_cols)) {
    fc <- cand[[fc_cols[i]]]
    pv <- cand[[pv_cols[i]]]
    cand$n_up_p01 <- cand$n_up_p01 + as.integer(!is.na(fc) & fc > 0 & !is.na(pv) & pv < p_cutoff)
  }
}
# 与「含脑/肺/骨」配体分数的相关（样本级：基因在各样本表达 vs 该样本是否脑）
for (org in c("brain", "lung", "bone")) {
  if (sum(si_use$organ == org) < 3) next
  y <- as.integer(si_use$organ == org)
  cors <- apply(mat_n[cand$gene, , drop = FALSE], 1, function(v) {
    suppressWarnings(stats::cor(v, y, method = "spearman"))
  })
  cand[[paste0("spearman_vs_", org)]] <- as.numeric(cors)
}
cand$priority <- cand$n_up_p01 +
  as.integer(!is.na(cand$mean_expr) & cand$mean_expr > median(cand$mean_expr, na.rm = TRUE))
cand <- cand[order(-cand$priority, -cand$mean_expr), ]
utils::write.csv(cand, file.path(result_dir, "01_CANDIDATE_MOLECULES_tumor_to_nerve.csv"),
                 row.names = FALSE)
if (has_pkg("writexl")) {
  tryCatch(writexl::write_xlsx(cand, file.path(result_dir, "01_CANDIDATE_MOLECULES_tumor_to_nerve.xlsx")),
           error = function(e) NULL)
}

# 未在面板上的经典配体 → 记入说明，避免误以为分析过
missing_classic <- setdiff(
  c("NGF", "ARTN", "GDNF", "NTN1", "SLIT2", "SEMA3A", "SEMA3C", "L1CAM", "NTF3", "NRG1"),
  rownames(mat_n)
)
writeLines(
  c("下列经典神经浸润配体不在 GSE175692 的 771 基因面板上，本分析无法评估：",
    paste(missing_classic, collapse = ", "),
    "",
    "面板上可用于问题1的基因见 01_CANDIDATE_MOLECULES_tumor_to_nerve.csv"),
  file.path(result_dir, "01_NOTE_genes_not_on_NanoString_panel.txt")
)

# -----------------------------------------------------------------------------
# 6. 部位排名总表（问题2）
# -----------------------------------------------------------------------------
organs_ok <- names(which(table(si_use$organ) >= 3))
site_rows <- lapply(organs_ok, function(org) {
  in_o <- si_use$organ == org
  rest <- !in_o
  wl <- tryCatch(stats::wilcox.test(si_use$score_ligand[in_o], si_use$score_ligand[rest])$p.value,
                 error = function(e) NA_real_)
  wr <- tryCatch(stats::wilcox.test(si_use$score_ligand_resid[in_o], si_use$score_ligand_resid[rest])$p.value,
                 error = function(e) NA_real_)
  ws <- tryCatch(stats::wilcox.test(si_use$score_stroma[in_o], si_use$score_stroma[rest])$p.value,
                 error = function(e) NA_real_)
  data.frame(
    organ = org,
    n = sum(in_o),
    mean_ligand = mean(si_use$score_ligand[in_o], na.rm = TRUE),
    mean_ligand_resid = mean(si_use$score_ligand_resid[in_o], na.rm = TRUE),
    mean_stroma = mean(si_use$score_stroma[in_o], na.rm = TRUE),
    p_ligand_vs_rest = wl,
    p_ligand_resid_vs_rest = wr,
    p_stroma_vs_rest = ws,
    stringsAsFactors = FALSE
  )
})
site <- do.call(rbind, site_rows)
site <- site[order(-site$mean_ligand_resid, -site$mean_ligand), ]
site$rank_by_neurotropic_ligand <- seq_len(nrow(site))
site$note <- ifelse(
  site$organ == "brain",
  "脑实质转移不等于外周神经浸润(PNI)；本面板几乎无GFAP/MBP，勿把脑基质信号当PNI",
  ""
)
utils::write.csv(site, file.path(result_dir, "02_SITE_RANK_neural_invasion.csv"), row.names = FALSE)
utils::write.csv(si_use, file.path(result_dir, "02_sample_scores_by_organ.csv"), row.names = FALSE)

# 箱线图：肺脑骨等
plot_df <- si_use[si_use$organ %in% organs_ok, ]
if (nrow(plot_df) > 0 && has_pkg("ggplot2")) {
  plot_df$organ <- factor(plot_df$organ, levels = site$organ)
  p1 <- ggplot2::ggplot(plot_df, ggplot2::aes(organ, score_ligand, fill = organ)) +
    ggplot2::geom_boxplot(outlier.size = 0.7, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1), legend.position = "none") +
    ggplot2::labs(
      title = "Q2: Neurotropic ligand score by metastatic site (GSE175692)",
      subtitle = "Higher = stronger tumor-side neurotrophic / invasion program on this panel",
      y = "score_ligand", x = NULL
    )
  save_gg(p1, file.path(result_dir, "02_SITE_RANK_ligand_boxplot"), 9, 5)

  p2 <- ggplot2::ggplot(plot_df, ggplot2::aes(organ, score_ligand_resid, fill = organ)) +
    ggplot2::geom_boxplot(outlier.size = 0.7, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1), legend.position = "none") +
    ggplot2::labs(
      title = "Q2: Ligand score residualized on stroma score",
      subtitle = "Use this ranking for site tropism; brain stroma is NOT PNI",
      y = "score_ligand | stroma", x = NULL
    )
  save_gg(p2, file.path(result_dir, "02_SITE_RANK_ligand_resid_boxplot"), 9, 5)
}

top_site <- site$organ[1]
top_msg <- paste0(
  "按神经亲和配体残差排名第1的继发部位: ", top_site,
  " (n=", site$n[1], ", mean_ligand_resid=",
  signif(site$mean_ligand_resid[1], 3), ")"
)
log_msg(top_msg)

# -----------------------------------------------------------------------------
# 7. 可选 ORA（上调候选基因）
# -----------------------------------------------------------------------------
if (has_pkg("clusterProfiler") && has_pkg("org.Hs.eg.db") && nrow(cand) >= 5) {
  up_genes <- cand$gene[cand$n_up_p01 > 0]
  if (length(up_genes) < 3) up_genes <- utils::head(cand$gene, 15)
  eg <- tryCatch(
    clusterProfiler::bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) NULL
  )
  if (!is.null(eg) && nrow(eg) >= 3) {
    ora_dir <- file.path(result_dir, "Q1_genes_promote_nerve", "ORA_candidate_ligands")
    dir.create(ora_dir, recursive = TRUE, showWarnings = FALSE)
    ego <- tryCatch(
      clusterProfiler::enrichGO(
        eg$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
        pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
      ),
      error = function(e) NULL
    )
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      df <- as.data.frame(ego)
      utils::write.csv(df, file.path(ora_dir, "ORA_GO_BP.csv"), row.names = FALSE)
      hit <- grepl(nerve_term_pat, df$Description, ignore.case = TRUE)
      if (any(hit)) {
        foc <- df[hit, , drop = FALSE]
        foc$genome_wide_rank <- match(which(hit), order(df$pvalue))
        utils::write.csv(foc, file.path(ora_dir, "FOCUS_nerve_invasion.csv"), row.names = FALSE)
      }
    }
  }
}

# -----------------------------------------------------------------------------
# 8. 方案说明
# -----------------------------------------------------------------------------
protocol <- c(
  "============================================================",
  "E:/R/Nerve RNA  |  GSE175692 独立分析",
  "数据源: ", src_lab,
  "============================================================",
  "",
  "【问题1】肿瘤高表达哪些基因可能促进神经浸润？",
  "  主表: results/01_CANDIDATE_MOLECULES_tumor_to_nerve.csv",
  "  详细: results/Q1_genes_promote_nerve/<器官>_vs_other_mets/FoldChange/",
  "  规则: 先 p < 0.01，再上调 FC >= 1.25 / 1.5 / 2",
  "  本队列是转移灶活检 bulk（NanoString 771 基因），不是空间邻域。",
  "  面板命中的关键候选: ", paste(lig_in, collapse = ", "),
  "  不在面板、无法评估: ", paste(missing_classic, collapse = ", "),
  "  关联 ≠ 因果；促浸润需后续功能实验验证。",
  "",
  "【问题2】神经浸润更倾向哪个继发部位？",
  "  主表: results/02_SITE_RANK_neural_invasion.csv",
  "  排序列: rank_by_neurotropic_ligand（按 mean_ligand_resid）",
  "  ", top_msg,
  "  箱线图: 02_SITE_RANK_ligand_boxplot.png / _ligand_resid_boxplot.png",
  "  重要: 脑实质转移 ≠ 外周神经浸润(PNI)。",
  "        本面板几乎没有 GFAP/MBP/MPZ，不要把脑组织信号当成 PNI。",
  "",
  "器官样本数: ", paste(sprintf("%s=%d", names(org_tab), as.integer(org_tab)), collapse = ", "),
  "",
  "运行方式:",
  "  setwd(\"E:/R/Nerve RNA\")",
  "  source(\"Nerve_RNA_nerve_infiltration.R\")"
)
writeLines(protocol, file.path(result_dir, "00_PROTOCOL_两个问题怎么看.txt"))
log_msg("Done.")
log_msg("Q1 -> ", file.path(result_dir, "01_CANDIDATE_MOLECULES_tumor_to_nerve.csv"))
log_msg("Q2 -> ", file.path(result_dir, "02_SITE_RANK_neural_invasion.csv"))
log_msg("Read -> ", file.path(result_dir, "00_PROTOCOL_两个问题怎么看.txt"))
