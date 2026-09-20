#!/usr/bin/env Rscript
# =============================================================================
# GSE273439_Mouse_breast.R
# 4T1 Visium FFPE：配对原发乳腺瘤 vs 肺（GSE273439）
#
# 用法（Windows）：
#   setwd("E:/R/Mouse Breast/GSE273439")
#   source("GSE273439_Mouse_breast.R")
# 或把本脚本放在 E:/R/Mouse Breast，数据在其子目录 GSE273439。
#
# 四个问题（「低表达促进转移」= 配对肺相对该鼠原位下调）：
#   1) 原位内低表达、与肺转移相关的下调基因（627 原位对 627 肺；628 对 628）
#      骨转移：本 GEO 没有骨切片，目录 03_ 里写明，不编造
#   2) 肺转移灶特异下调（相对未受累肺实质），不能做骨特异
#   3) 原位瘤按施旺细胞 / 神经营养因子 / 轴突导向打分，找负相关基因
#   4) 三种神经分数分别与配对肺转移的关系
#
# 入选：p < 0.05，下调倍数 down_FC = 原位/肺 >= 1 与 >= 1.25
# （即肺相对原位的 FC < 1 与 < 1/1.25）
#
# 数据注意：
#   - 两只 BALB/c，各一块 TM + 一块 LUNG，Visium FFPE，mm10
#   - 不是 bulk TPM；spot 不是独立生物学重复，p 值偏乐观，故以「两鼠都过阈值」为主结论
#   - 没有病理 PNI；神经浸润用三套基因集在原位 spot 上打分
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")
options(clusterProfiler.download.method = "auto")

# -----------------------------------------------------------------------------
# 0. 依赖
# -----------------------------------------------------------------------------
cran_required <- c(
  "Matrix", "dplyr", "tidyr", "tibble", "stringr", "ggplot2", "ggrepel",
  "pheatmap", "RColorBrewer", "matrixStats", "cowplot", "writexl"
)
cran_optional <- c("ggvenn", "GSVA", "jsonlite")
skip_enrich <- nzchar(Sys.getenv("GSE273439_SKIP_ENRICHMENT"))
bioc_required <- c("limma")
bioc_optional <- c(
  "AnnotationDbi", "org.Mm.eg.db", "clusterProfiler", "enrichplot",
  "DOSE", "fgsea", "msigdbr", "ReactomePA"
)

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
if (!isTRUE(skip_enrich)) {
  install_if_missing(cran_optional, bioc = FALSE, required = FALSE)
  install_if_missing(bioc_optional, bioc = TRUE, required = FALSE)
} else {
  miss_opt <- c(cran_optional, bioc_optional)
  miss_opt <- miss_opt[!vapply(miss_opt, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_opt) > 0) {
    message("GSE273439_SKIP_ENRICHMENT: skip installing ", paste(miss_opt, collapse = ", "))
  }
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
resolve_project_dir <- function() {
  env_dir <- Sys.getenv("GSE273439_DIR", unset = "")
  candidates <- c(
    env_dir,
    "E:/R/Mouse Breast/GSE273439",
    "E:\\R\\Mouse Breast\\GSE273439",
    "E:/R/Mouse Breast",
    file.path(getwd(), "GSE273439"),
    getwd()
  )
  candidates <- unique(candidates[nzchar(candidates)])
  looks_like <- function(d) {
    if (!dir.exists(d)) return(FALSE)
    mtx <- list.files(d, pattern = "matrix\\.mtx(\\.gz)?$", recursive = TRUE, ignore.case = TRUE)
    tar <- list.files(d, pattern = "GSE273439_RAW\\.tar$", recursive = TRUE, ignore.case = TRUE)
    gsm <- list.files(d, pattern = "GSM84284(08|09|10|11)", recursive = TRUE, ignore.case = TRUE)
    length(mtx) > 0 || length(tar) > 0 || length(gsm) > 0
  }
  for (d in candidates) {
    if (looks_like(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

project_dir <- resolve_project_dir()
result_dir  <- file.path(project_dir, "results_GSE273439")
log_dir     <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("GSE273439_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

p_cutoff   <- 0.05
fc_cutoffs <- c("FC_1" = 1, "FC_1.25" = 1.25)

log_msg("Project dir: ", project_dir)
log_msg("Results dir: ", result_dir)
log_msg("Gate: p < ", p_cutoff, " ; down-FC (primary/lung) >= 1 and >= 1.25")
log_msg("Paired only: Mouse627 TM vs Mouse627 LUNG; Mouse628 TM vs Mouse628 LUNG")

# -----------------------------------------------------------------------------
# 2. 发现并整理 4 个 Visium 样品
# -----------------------------------------------------------------------------
maybe_untar <- function(root) {
  mtx <- list.files(root, pattern = "matrix\\.mtx(\\.gz)?$", recursive = TRUE)
  if (length(mtx) >= 4) return(invisible(TRUE))
  tars <- list.files(root, pattern = "GSE273439_RAW\\.tar$", recursive = TRUE, full.names = TRUE)
  if (length(tars) == 0) return(invisible(FALSE))
  dest <- file.path(root, "_RAW_extracted")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  log_msg("Untarring ", tars[[1]], " -> ", dest)
  utils::untar(tars[[1]], exdir = dest)
  invisible(TRUE)
}
maybe_untar(project_dir)

infer_sample <- function(path) {
  b <- basename(path)
  d <- basename(dirname(path))
  txt <- paste(b, d, path)
  gsm_map <- c(
    GSM8428408 = "Mouse627_TM",
    GSM8428409 = "Mouse627_LUNG",
    GSM8428410 = "Mouse628_TM",
    GSM8428411 = "Mouse628_LUNG"
  )
  for (g in names(gsm_map)) {
    if (grepl(g, txt, ignore.case = TRUE)) return(unname(gsm_map[[g]]))
  }
  if (grepl("627.*LUNG|LUNG.*627", txt, ignore.case = TRUE)) return("Mouse627_LUNG")
  if (grepl("628.*LUNG|LUNG.*628", txt, ignore.case = TRUE)) return("Mouse628_LUNG")
  if (grepl("627.*(TM|tumor|primary)|TM.*627", txt, ignore.case = TRUE)) return("Mouse627_TM")
  if (grepl("628.*(TM|tumor|primary)|TM.*628", txt, ignore.case = TRUE)) return("Mouse628_TM")
  NA_character_
}

pick_one <- function(files, pattern) {
  hit <- files[grepl(pattern, basename(files), ignore.case = TRUE)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

discover_visium <- function(root) {
  allf <- list.files(root, recursive = TRUE, full.names = TRUE)
  allf <- allf[!grepl("results_GSE273439|_seurat_input", allf)]
  mtx <- allf[grepl("matrix\\.mtx(\\.gz)?$", basename(allf), ignore.case = TRUE)]
  if (length(mtx) == 0) {
    stop(
      "未找到 Visium matrix.mtx.gz。请把 GSE273439_RAW.tar 解压到: ", root,
      "\n需要四个样品：627_TM / 627_LUNG / 628_TM / 628_LUNG"
    )
  }
  rows <- list()
  for (m in mtx) {
    sid <- infer_sample(m)
    if (is.na(sid)) {
      log_msg("Skip unmatched mtx: ", m)
      next
    }
    gsm <- regmatches(basename(m), regexpr("GSM[0-9]+", basename(m), ignore.case = TRUE))
    if (length(gsm) == 1 && nzchar(gsm)) {
      pool <- allf[grepl(gsm, basename(allf), ignore.case = TRUE)]
    } else {
      pool <- list.files(dirname(m), full.names = TRUE)
    }
    rows[[sid]] <- data.frame(
      sample = sid,
      mouse = sub("_.*", "", sid),
      tissue = sub(".*_", "", sid),
      mtx = m,
      features = pick_one(pool, "features\\.tsv"),
      barcodes = pick_one(pool, "barcodes\\.tsv"),
      positions = pick_one(pool, "tissue_positions"),
      scalefactors = pick_one(pool, "scalefactors_json"),
      hires = pick_one(pool, "tissue_hires_image"),
      lowres = pick_one(pool, "tissue_lowres_image"),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) stop("找到 mtx 但无法识别 627/628 TM/LUNG 文件名")
  tab <- do.call(rbind, rows)
  rownames(tab) <- tab$sample
  need <- c("Mouse627_TM", "Mouse627_LUNG", "Mouse628_TM", "Mouse628_LUNG")
  miss <- setdiff(need, tab$sample)
  if (length(miss) > 0) {
    log_msg("WARNING missing samples: ", paste(miss, collapse = ", "))
  }
  tab
}

sample_files <- discover_visium(project_dir)
log_msg("Discovered samples: ", paste(sample_files$sample, collapse = ", "))
print(sample_files[, c("sample", "mouse", "tissue")])

open_maybe_gz <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
}

read_tsv_nogz <- function(path, header = FALSE) {
  con <- open_maybe_gz(path)
  on.exit(close(con), add = TRUE)
  utils::read.delim(con, header = header, stringsAsFactors = FALSE, check.names = FALSE)
}

read_positions <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  con <- open_maybe_gz(path)
  first <- readLines(con, n = 1)
  close(con)
  has_header <- grepl("barcode", first, ignore.case = TRUE)
  sep <- if (grepl(",", first, fixed = TRUE)) "," else "\t"
  con <- open_maybe_gz(path)
  on.exit(close(con), add = TRUE)
  df <- utils::read.table(
    con, header = has_header, sep = sep, stringsAsFactors = FALSE,
    check.names = FALSE, comment.char = "", quote = ""
  )
  if (!has_header) {
    nms <- c(
      "barcode", "in_tissue", "array_row", "array_col",
      "pxl_row_in_fullres", "pxl_col_in_fullres"
    )
    names(df)[seq_len(min(length(nms), ncol(df)))] <- nms[seq_len(min(length(nms), ncol(df)))]
  }
  names(df) <- tolower(gsub("[^A-Za-z0-9]+", "_", names(df)))
  names(df) <- gsub("_+$", "", names(df))
  if (!"barcode" %in% names(df)) names(df)[1] <- "barcode"
  if (!"in_tissue" %in% names(df) && ncol(df) >= 2) names(df)[2] <- "in_tissue"
  df$barcode <- as.character(df[[1]])
  if ("in_tissue" %in% names(df)) {
    df$in_tissue <- as.integer(df$in_tissue)
  } else {
    df$in_tissue <- 1L
  }
  df
}

collapse_duplicate_symbols <- function(mat, symbols) {
  symbols <- as.character(symbols)
  empty <- is.na(symbols) | !nzchar(symbols)
  if (any(empty)) symbols[empty] <- paste0("NOVEL_", seq_len(sum(empty)))
  if (!any(duplicated(symbols))) {
    rownames(mat) <- symbols
    return(mat)
  }
  log_msg("Collapsing duplicate symbols by sum")
  uniq <- unique(symbols)
  idx <- split(seq_along(symbols), factor(symbols, levels = uniq))
  parts <- lapply(uniq, function(s) {
    rows <- idx[[s]]
    if (length(rows) == 1L) {
      mat[rows, , drop = FALSE]
    } else {
      cs <- Matrix::colSums(mat[rows, , drop = FALSE])
      out <- Matrix::Matrix(cs, nrow = 1L, sparse = TRUE)
      rownames(out) <- s
      colnames(out) <- colnames(mat)
      out
    }
  })
  out <- do.call(rbind, parts)
  rownames(out) <- uniq
  methods::as(out, "dgCMatrix")
}

load_one_sample <- function(row) {
  if (is.na(row$features) || is.na(row$barcodes)) {
    stop("样品 ", row$sample, " 缺少 features.tsv 或 barcodes.tsv")
  }
  log_msg("Loading ", row$sample)
  feat <- read_tsv_nogz(row$features, header = FALSE)
  bc <- read_tsv_nogz(row$barcodes, header = FALSE)[[1]]
  mat <- Matrix::readMM(row$mtx)
  mat <- tryCatch(
    methods::as(mat, "CsparseMatrix"),
    error = function(e) methods::as(mat, "dgCMatrix")
  )
  if (nrow(mat) != nrow(feat)) {
    stop(row$sample, ": matrix rows (", nrow(mat), ") != features (", nrow(feat), ")")
  }
  if (ncol(mat) != length(bc)) {
    stop(row$sample, ": matrix cols (", ncol(mat), ") != barcodes (", length(bc), ")")
  }
  ens <- as.character(feat[[1]])
  sym <- if (ncol(feat) >= 2) as.character(feat[[2]]) else ens
  sym <- ifelse(!is.na(sym) & nzchar(sym), sym, ens)
  colnames(mat) <- paste(row$sample, bc, sep = "|")
  mat <- collapse_duplicate_symbols(mat, sym)
  ens_map <- ens[!duplicated(sym)]
  names(ens_map) <- unique(sym)
  pos <- read_positions(row$positions)
  md <- data.frame(
    cell = colnames(mat),
    barcode = bc,
    sample = row$sample,
    mouse = row$mouse,
    tissue = row$tissue,
    stringsAsFactors = FALSE
  )
  if (!is.null(pos)) {
    pos$barcode <- as.character(pos$barcode)
    md <- merge(md, pos, by = "barcode", all.x = TRUE, sort = FALSE)
    md <- md[match(colnames(mat), md$cell), ]
  } else {
    md$in_tissue <- 1L
    md$array_row <- NA_integer_
    md$array_col <- NA_integer_
  }
  if (!"in_tissue" %in% names(md) || all(is.na(md$in_tissue))) md$in_tissue <- 1L
  md$in_tissue[is.na(md$in_tissue)] <- 1L
  list(counts = mat, meta = md, ensembl = ens_map)
}

loaded <- lapply(seq_len(nrow(sample_files)), function(i) load_one_sample(sample_files[i, ]))
names(loaded) <- sample_files$sample

# 基因并集
all_genes <- Reduce(union, lapply(loaded, function(x) rownames(x$counts)))
align_mat <- function(mat, genes) {
  miss <- setdiff(genes, rownames(mat))
  if (length(miss) > 0) {
    extra <- Matrix::Matrix(0, nrow = length(miss), ncol = ncol(mat), sparse = TRUE)
    rownames(extra) <- miss
    colnames(extra) <- colnames(mat)
    mat <- rbind(mat, extra)
  }
  mat[genes, , drop = FALSE]
}
counts <- do.call(cbind, lapply(loaded, function(x) align_mat(x$counts, all_genes)))
meta <- do.call(rbind, lapply(loaded, function(x) x$meta))
rownames(meta) <- meta$cell
meta <- meta[colnames(counts), , drop = FALSE]
meta$nCount <- Matrix::colSums(counts)
meta$nFeature <- Matrix::colSums(counts > 0)
mt_genes <- grepl("^mt-", rownames(counts), ignore.case = TRUE)
meta$percent_mt <- if (any(mt_genes)) {
  as.numeric(Matrix::colSums(counts[mt_genes, , drop = FALSE]) / pmax(meta$nCount, 1) * 100)
} else {
  0
}

log_msg("Combined: ", nrow(counts), " genes x ", ncol(counts), " spots")

# 基因注释
map_symbols <- function(symbols) {
  uniq <- unique(symbols)
  if (!has_pkg("org.Mm.eg.db") || !has_pkg("AnnotationDbi")) {
    log_msg("org.Mm.eg.db not available; keep gene symbols without Entrez map")
    return(data.frame(gene = uniq, ensembl = NA_character_, entrez = NA_character_,
                      stringsAsFactors = FALSE))
  }
  ens <- tryCatch(
    AnnotationDbi::mapIds(org.Mm.eg.db, keys = uniq, column = "ENSEMBL",
                          keytype = "SYMBOL", multiVals = "first"),
    error = function(e) setNames(rep(NA_character_, length(uniq)), uniq)
  )
  ent <- tryCatch(
    AnnotationDbi::mapIds(org.Mm.eg.db, keys = uniq, column = "ENTREZID",
                          keytype = "SYMBOL", multiVals = "first"),
    error = function(e) setNames(rep(NA_character_, length(uniq)), uniq)
  )
  data.frame(
    gene = uniq,
    ensembl = unname(ens[uniq]),
    entrez = unname(ent[uniq]),
    stringsAsFactors = FALSE
  )
}
id_map <- map_symbols(rownames(counts))

# -----------------------------------------------------------------------------
# 3. QC、过滤、标准化
# -----------------------------------------------------------------------------
qc_dir <- file.path(result_dir, "00_QC")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(sample_files, file.path(qc_dir, "discovered_files.csv"), row.names = FALSE)

save_gg <- function(plot, path_stub, width = 8, height = 6) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".pdf"), plot, width = width, height = height),
           error = function(e) log_msg("pdf ggsave failed: ", e$message))
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".png"), plot, width = width, height = height, dpi = 300),
           error = function(e) log_msg("png ggsave failed: ", e$message))
}

p_qc <- ggplot2::ggplot(meta, ggplot2::aes(x = sample, y = nCount, fill = tissue)) +
  ggplot2::geom_violin(scale = "width") +
  ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.3) +
  ggplot2::scale_y_log10() +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
  ggplot2::labs(title = "UMI per spot", x = NULL)
save_gg(p_qc, file.path(qc_dir, "nCount_violin"), 9, 5)

keep_spot <- meta$in_tissue == 1 & meta$nCount >= 50 & meta$nFeature >= 20
if (sum(keep_spot) < 50) {
  log_msg("QC too strict; keep all in-tissue spots")
  keep_spot <- meta$in_tissue == 1
}
log_msg("Keep spots: ", sum(keep_spot), " / ", ncol(counts))
counts <- counts[, keep_spot, drop = FALSE]
meta <- meta[keep_spot, , drop = FALSE]

keep_gene <- Matrix::rowSums(counts > 0) >= 10
if (sum(keep_gene) < 200) keep_gene <- Matrix::rowSums(counts) > 0
log_msg("Keep genes: ", sum(keep_gene), " / ", nrow(counts))
counts <- counts[keep_gene, , drop = FALSE]
id_map <- id_map[match(rownames(counts), id_map$gene), ]
id_map$gene <- rownames(counts)

lib <- pmax(meta$nCount, 1)
log_mat <- as.matrix(Matrix::t(Matrix::t(counts) / lib * 1e4))
log_mat <- log2(log_mat + 1)
storage.mode(log_mat) <- "double"
log_msg("Normalization: log2(CP10k + 1) on in-tissue spots")

# 肿瘤上皮分数，用来圈肺转移灶
tumor_markers <- c(
  "Epcam", "Krt8", "Krt18", "Krt19", "Krt14", "Krt5", "Cdh1", "Mki67", "Krt7"
)
tumor_use <- intersect(tumor_markers, rownames(log_mat))
log_msg("Tumor markers present: ", paste(tumor_use, collapse = ", "))
if (length(tumor_use) >= 2) {
  z <- t(scale(t(log_mat[tumor_use, , drop = FALSE])))
  meta$tumor_score <- colMeans(z, na.rm = TRUE)
} else if (length(tumor_use) == 1) {
  meta$tumor_score <- as.numeric(scale(log_mat[tumor_use, ]))
} else {
  meta$tumor_score <- 0
  log_msg("WARNING: no epithelial markers; cannot isolate lung met foci")
}

meta$compartment <- "other"
for (ms in unique(meta$mouse)) {
  tm <- meta$mouse == ms & meta$tissue == "TM"
  lu <- meta$mouse == ms & meta$tissue == "LUNG"
  tm_cut <- stats::median(meta$tumor_score[tm], na.rm = TRUE)
  lu_cut <- stats::quantile(meta$tumor_score[lu], 0.75, na.rm = TRUE)
  # 肺灶：高于该肺 75% 分位，且不低于该鼠原位中位数的 50%
  foci <- lu & meta$tumor_score >= lu_cut & meta$tumor_score >= (tm_cut * 0.5)
  uninv <- lu & !foci
  tm_tumor <- tm & meta$tumor_score >= tm_cut
  meta$compartment[tm] <- "TM_all"
  meta$compartment[tm_tumor] <- "TM_tumor"
  meta$compartment[uninv] <- "LUNG_uninvolved"
  meta$compartment[foci] <- "LUNG_foci"
  log_msg(
    ms, " TM_tumor=", sum(meta$compartment == "TM_tumor" & meta$mouse == ms),
    " LUNG_foci=", sum(meta$compartment == "LUNG_foci" & meta$mouse == ms),
    " LUNG_uninvolved=", sum(meta$compartment == "LUNG_uninvolved" & meta$mouse == ms)
  )
}
utils::write.csv(meta, file.path(qc_dir, "spot_metadata.csv"), row.names = FALSE)

comp_tab <- as.data.frame(table(mouse = meta$mouse, tissue = meta$tissue, compartment = meta$compartment))
comp_tab <- comp_tab[comp_tab$Freq > 0, ]
utils::write.csv(comp_tab, file.path(qc_dir, "compartment_counts.csv"), row.names = FALSE)

plot_spatial <- function(md, color_col, title, outfile, discrete = FALSE) {
  if (!"array_col" %in% names(md) || all(is.na(md$array_col))) {
    log_msg("No array coordinates, skip spatial: ", title)
    return(invisible(NULL))
  }
  df <- md[is.finite(md$array_col) & is.finite(md$array_row), ]
  if (nrow(df) < 10) return(invisible(NULL))
  df$color_value <- df[[color_col]]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = array_col, y = -array_row, color = color_value)) +
    ggplot2::geom_point(size = 0.6) +
    ggplot2::coord_fixed() +
    ggplot2::facet_wrap(~ sample) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::labs(title = title, color = color_col)
  if (!discrete) {
    p <- p + ggplot2::scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(9, "RdBu")))
  }
  save_gg(p, outfile, 10, 8)
}
plot_spatial(meta, "nCount", "UMI (in tissue)", file.path(qc_dir, "spatial_nCount"))
plot_spatial(meta, "tumor_score", "Epithelial / tumor module score", file.path(qc_dir, "spatial_tumor_score"))
plot_spatial(meta, "compartment", "Spot compartments", file.path(qc_dir, "spatial_compartment"), discrete = TRUE)

# 四样品伪 bulk（热图用）
sample_levels <- c("Mouse627_TM", "Mouse627_LUNG", "Mouse628_TM", "Mouse628_LUNG")
sample_levels <- intersect(sample_levels, unique(meta$sample))
pb <- sapply(sample_levels, function(s) {
  idx <- which(meta$sample == s)
  matrixStats::rowMeans2(log_mat[, idx, drop = FALSE])
})
rownames(pb) <- rownames(log_mat)
sample_info <- data.frame(
  sample = sample_levels,
  mouse = sub("_.*", "", sample_levels),
  tissue = sub(".*_", "", sample_levels),
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_info$sample
pb_tm <- rowMeans(pb[, grepl("_TM$", colnames(pb)), drop = FALSE])
tm_pct <- rank(pb_tm, ties.method = "average") / length(pb_tm) * 100

# -----------------------------------------------------------------------------
# 4. 作图 / 富集工具
# -----------------------------------------------------------------------------
note_empty <- function(stub, msg) {
  dir.create(dirname(stub), recursive = TRUE, showWarnings = FALSE)
  writeLines(msg, paste0(stub, "_EMPTY.txt"))
}

write_table <- function(df, stub) {
  dir.create(dirname(stub), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, paste0(stub, ".csv"), row.names = FALSE)
  tryCatch(writexl::write_xlsx(df, paste0(stub, ".xlsx")),
           error = function(e) log_msg("xlsx failed: ", e$message))
}

map_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (length(symbols) == 0) return(data.frame(gene = character(), entrez = character()))
  from_map <- id_map[id_map$gene %in% symbols & !is.na(id_map$entrez) & nzchar(id_map$entrez),
                     c("gene", "entrez")]
  leftover <- setdiff(symbols, from_map$gene)
  extra <- data.frame(gene = character(), entrez = character())
  if (length(leftover) > 0 && has_pkg("clusterProfiler") && has_pkg("org.Mm.eg.db")) {
    m <- tryCatch(
      clusterProfiler::bitr(leftover, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db),
      error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
    )
    if (nrow(m) > 0) extra <- data.frame(gene = m$SYMBOL, entrez = m$ENTREZID)
  }
  out <- rbind(from_map, extra)
  out[!duplicated(out$gene), ]
}

ranked_entrez <- function(de, stat_col = "log2FC") {
  mp <- map_to_entrez(de$gene)
  de2 <- merge(de, mp, by = "gene")
  de2 <- de2[!is.na(de2$entrez) & !is.na(de2[[stat_col]]), ]
  de2 <- de2[order(abs(de2[[stat_col]]), decreasing = TRUE), ]
  de2 <- de2[!duplicated(de2$entrez), ]
  stats <- de2[[stat_col]]
  names(stats) <- de2$entrez
  sort(stats, decreasing = TRUE)
}

plot_volcano <- function(de, highlight, title, outfile, fc_line = 1, lfc_col = "log2FC", p_col = "pvalue") {
  df <- de
  df$log2FC <- df[[lfc_col]]
  if (p_col %in% names(df) && any(!is.na(df[[p_col]]))) {
    df$y <- -log10(pmax(df[[p_col]], 1e-300))
    ylab <- paste0("-log10(", p_col, ")")
    hline <- -log10(p_cutoff)
  } else {
    df$y <- df$AveExpr
    ylab <- "Average expression"
    hline <- NULL
  }
  df$set <- ifelse(df$gene %in% highlight, "selected", "other")
  df$label <- ifelse(df$gene %in% utils::head(highlight, 15), df$gene, NA)
  lfc_line <- log2(fc_line)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = y, color = set)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.2) +
    ggplot2::scale_color_manual(values = c(other = "grey70", selected = "#1F77B4")) +
    ggplot2::geom_vline(xintercept = c(-lfc_line, lfc_line), linetype = 2, color = "grey40") +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, x = "log2FC (lung - matched primary)", y = ylab, color = NULL)
  if (!is.null(hline)) p <- p + ggplot2::geom_hline(yintercept = hline, linetype = 2, color = "grey40")
  save_gg(p, outfile)
}

plot_heatmap <- function(heat_mat, sample_info, genes, title, outfile, extra_ann = NULL) {
  genes <- intersect(genes, rownames(heat_mat))
  if (length(genes) > 200) {
    genes <- genes[seq_len(200)]
  }
  if (length(genes) < 2) {
    note_empty(outfile, "heatmap skipped (<2 genes)")
    return(invisible(NULL))
  }
  sub <- heat_mat[genes, , drop = FALSE]
  ann <- data.frame(
    Mouse = sample_info$mouse[match(colnames(sub), sample_info$sample)],
    Tissue = sample_info$tissue[match(colnames(sub), sample_info$sample)],
    row.names = colnames(sub)
  )
  if (!is.null(extra_ann)) {
    extra_ann <- extra_ann[match(colnames(sub), rownames(extra_ann)), , drop = FALSE]
    ann <- cbind(ann, extra_ann)
  }
  pal <- list(
    Tissue = c(TM = "#4C78A8", LUNG = "#F58518"),
    Mouse = c(Mouse627 = "#54A24B", Mouse628 = "#E45756")
  )
  draw_hm <- function() {
    args <- list(
      mat = sub, scale = "row", annotation_col = ann, annotation_colors = pal,
      show_rownames = nrow(sub) <= 80, fontsize_row = 6, main = title,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100)
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
  grDevices::png(paste0(outfile, ".png"), width = 2400, height = max(1600, 40 * nrow(sub) + 400), res = 300)
  draw_hm()
  grDevices::dev.off()
}

plot_de_bar <- function(sub, title, outfile, lfc_col = "log2FC") {
  if (nrow(sub) == 0) return(invisible(NULL))
  df <- sub
  df$log2FC <- df[[lfc_col]]
  df <- df[order(df$log2FC, decreasing = FALSE), , drop = FALSE]
  if (nrow(df) > 60) df <- rbind(utils::head(df, 30), utils::tail(df, 30))
  df$gene <- factor(df$gene, levels = rev(unique(df$gene)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = gene, y = log2FC)) +
    ggplot2::geom_col(fill = "#1F77B4", width = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title, x = NULL, y = "log2FC (lung - primary); negative = down in lung")
  save_gg(p, outfile, width = 8, height = max(5, min(16, 0.22 * nrow(df) + 2)))
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
    msigdbr::msigdbr(species = "Mus musculus", collection = "H"),
    error = function(e) msigdbr::msigdbr(species = "Mus musculus", category = "H")
  )
  gene_col <- intersect(c("ncbi_gene", "entrez_gene"), names(msig))[1]
  msig[, c("gs_name", gene_col)]
}

neural_keyword <- function(id, desc) {
  txt <- tolower(paste(id, desc))
  if (grepl("axon|neur|nerve|pni|perineur|schwann|myelin|synap|neural crest|gangli|glia", txt)) {
    return("neural_pni")
  }
  NA_character_
}

export_neural_terms <- function(x, stub, title) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) return(invisible(NULL))
  df <- as.data.frame(x)
  desc <- if ("Description" %in% names(df)) df$Description else df$ID
  df$genome_wide_rank <- seq_len(nrow(df))
  df$focus_class <- vapply(seq_len(nrow(df)), function(i) neural_keyword(df$ID[i], desc[i]), character(1))
  hit <- df[!is.na(df$focus_class), , drop = FALSE]
  if (nrow(hit) == 0) {
    note_empty(paste0(stub, "_FOCUS_neural_pni"), "no neural/PNI terms in this result")
    return(invisible(NULL))
  }
  utils::write.csv(hit, paste0(stub, "_FOCUS_neural_pni.csv"), row.names = FALSE)
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
  export_neural_terms(x, stub, title)
}

plot_gsea_object <- function(x, stub, title) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) {
    note_empty(stub, "no GSEA terms")
    return(invisible(NULL))
  }
  df <- as.data.frame(x)
  utils::write.csv(df, paste0(stub, ".csv"), row.names = FALSE)
  nshow <- min(15, nrow(df))
  try_save_plot(function() {
    p <- enrichplot::dotplot(x, showCategory = nshow, split = ".sign")
    tryCatch(p + ggplot2::facet_grid(. ~ .sign) + ggplot2::ggtitle(title),
             error = function(e) p + ggplot2::ggtitle(title))
  }, paste0(stub, "_dotplot"), 10, 7)
  ncurve <- min(5, nrow(df))
  try_save_plot(function() enrichplot::gseaplot2(x, geneSetID = seq_len(ncurve), pvalue_table = TRUE, title = title),
                paste0(stub, "_gseaplot"), 10, 8)
  export_neural_terms(x, stub, title)
}

run_ora_plots <- function(genes, de_sub, outdir, label, tag, lfc_col = "log2FC") {
  if (!has_pkg("clusterProfiler")) {
    note_empty(file.path(outdir, paste0(tag, "_ORA")), "clusterProfiler not installed")
    return(invisible(NULL))
  }
  go_dir <- file.path(outdir, "GO")
  pw_dir <- file.path(outdir, "Pathway")
  kg_dir <- file.path(outdir, "KEGG")
  dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(kg_dir, recursive = TRUE, showWarnings = FALSE)
  pref <- paste0(tag, "_")
  mp <- map_to_entrez(genes)
  entrez <- unique(mp$entrez)
  fc_sym <- setNames(de_sub[[lfc_col]], de_sub$gene)
  if (length(entrez) < 3) {
    log_msg("ORA skipped, mapped genes < 3: ", outdir)
    note_empty(file.path(go_dir, paste0(pref, "ORA_GO")), "too few mapped genes")
    note_empty(file.path(pw_dir, paste0(pref, "ORA_Pathway")), "too few mapped genes")
    note_empty(file.path(kg_dir, paste0(pref, "ORA_KEGG")), "too few mapped genes")
    return(invisible(NULL))
  }
  for (ont in c("BP", "MF", "CC")) {
    ego <- enrich_or_relax(
      function() clusterProfiler::enrichGO(
        gene = entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = ont,
        pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
      ),
      function() clusterProfiler::enrichGO(
        gene = entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = ont,
        pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
      ),
      paste("enrichGO", ont)
    )
    plot_ora_object(ego, file.path(go_dir, paste0(pref, "ORA_GO_", ont)),
                    title_maybe_relaxed(ego, paste(label, "| ORA GO", ont)), fold_change = fc_sym)
  }
  ek <- enrich_or_relax(
    function() clusterProfiler::enrichKEGG(
      gene = entrez, organism = "mmu", pvalueCutoff = 0.05, qvalueCutoff = 0.2
    ),
    function() clusterProfiler::enrichKEGG(
      gene = entrez, organism = "mmu", pvalueCutoff = 1, qvalueCutoff = 1
    ),
    "enrichKEGG"
  )
  if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) {
    ek <- tryCatch(clusterProfiler::setReadable(ek, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"),
                   error = function(e) ek)
  }
  plot_ora_object(ek, file.path(kg_dir, paste0(pref, "ORA_KEGG")),
                  title_maybe_relaxed(ek, paste(label, "| ORA KEGG")), fold_change = fc_sym)
  if (has_pkg("ReactomePA")) {
    er <- enrich_or_relax(
      function() ReactomePA::enrichPathway(
        gene = entrez, organism = "mouse", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
      ),
      function() ReactomePA::enrichPathway(
        gene = entrez, organism = "mouse", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
      ),
      "enrichPathway"
    )
    plot_ora_object(er, file.path(pw_dir, paste0(pref, "ORA_Reactome_pathway")),
                    title_maybe_relaxed(er, paste(label, "| ORA Reactome")), fold_change = fc_sym)
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
    hm <- tryCatch(clusterProfiler::setReadable(hm, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"),
                   error = function(e) hm)
  }
  plot_ora_object(hm, file.path(pw_dir, paste0(pref, "ORA_MSigDB_Hallmark_pathway")),
                  title_maybe_relaxed(hm, paste(label, "| ORA Hallmark")), fold_change = fc_sym)
  writeLines(
    c("This GO/Pathway/KEGG folder is ORA, NOT GSEA.",
      "GSEA files are in GSEA/ and start with GSEA_."),
    file.path(outdir, paste0(pref, "00_ORA_is_not_GSEA.txt"))
  )
}

run_gsea_full <- function(de, outdir, label, stat_col = "log2FC") {
  if (!has_pkg("clusterProfiler")) {
    note_empty(file.path(outdir, "GSEA_all"), "clusterProfiler not installed")
    return(invisible(NULL))
  }
  gsea_dir <- file.path(outdir, "GSEA")
  dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
  stats <- ranked_entrez(de, stat_col = stat_col)
  if (length(stats) < 10) {
    note_empty(file.path(gsea_dir, "GSEA_all"), "too few ranked genes")
    return(invisible(NULL))
  }
  gsea_one <- function(fun, lab) {
    enrich_or_relax(
      function() fun(pvalueCutoff = 0.05, minGSSize = 10),
      function() fun(pvalueCutoff = 1, minGSSize = 5),
      lab
    )
  }
  go_bp <- gsea_one(function(pvalueCutoff, minGSSize) {
    clusterProfiler::gseGO(
      geneList = stats, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "ENTREZID",
      minGSSize = minGSSize, maxGSSize = 500, pvalueCutoff = pvalueCutoff,
      verbose = FALSE, eps = 0
    )
  }, "gseGO_BP")
  plot_gsea_object(go_bp, file.path(gsea_dir, "GSEA_GO_BP"), paste(label, "| GSEA GO BP"))
  kegg <- gsea_one(function(pvalueCutoff, minGSSize) {
    clusterProfiler::gseKEGG(
      geneList = stats, organism = "mmu", minGSSize = minGSSize, maxGSSize = 500,
      pvalueCutoff = pvalueCutoff, verbose = FALSE, eps = 0
    )
  }, "gseKEGG")
  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
    kegg <- tryCatch(clusterProfiler::setReadable(kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"),
                     error = function(e) kegg)
  }
  plot_gsea_object(kegg, file.path(gsea_dir, "GSEA_KEGG"), paste(label, "| GSEA KEGG"))
  term2gene <- tryCatch(msig_hallmark_map(), error = function(e) NULL)
  if (!is.null(term2gene)) {
    hm <- gsea_one(function(pvalueCutoff, minGSSize) {
      clusterProfiler::GSEA(
        geneList = stats, TERM2GENE = term2gene, minGSSize = minGSSize,
        maxGSSize = 500, pvalueCutoff = pvalueCutoff, eps = 0, verbose = FALSE
      )
    }, "GSEA_Hallmark")
    plot_gsea_object(hm, file.path(gsea_dir, "GSEA_Hallmark"), paste(label, "| GSEA Hallmark"))
  }
}

emit_subset <- function(comp_name, sub, tag, title, outdir, full_de, heat_mat, sample_info,
                        lfc_col = "log2FC", p_col = "pvalue", fc_line = 1, extra_ann = NULL,
                        do_ora = FALSE) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(paste("comparison:", comp_name), paste("subset:", tag), paste("n_genes:", nrow(sub))),
    file.path(outdir, paste0("00_", tag, "_THIS_FOLDER.txt"))
  )
  write_table(sub, file.path(outdir, paste0(tag, "_selected_genes")))
  log_msg(comp_name, " ", tag, ": n = ", nrow(sub))
  if (nrow(sub) == 0) {
    writeLines("no genes", file.path(outdir, paste0(tag, "_EMPTY.txt")))
    return(invisible(NULL))
  }
  tryCatch(plot_de_bar(sub, paste(title, "| selected genes"), file.path(outdir, paste0(tag, "_log2FC_barplot")), lfc_col),
           error = function(e) log_msg("bar failed: ", e$message))
  tryCatch(plot_volcano(full_de, sub$gene, title, file.path(outdir, paste0(tag, "_volcano")),
                        fc_line = fc_line, lfc_col = lfc_col, p_col = p_col),
           error = function(e) log_msg("volcano failed: ", e$message))
  tryCatch(plot_heatmap(heat_mat, sample_info, sub$gene, title, file.path(outdir, paste0(tag, "_heatmap")), extra_ann),
           error = function(e) {
             while (grDevices::dev.cur() > 1) grDevices::dev.off()
             log_msg("heatmap failed: ", e$message)
           })
  if (isTRUE(do_ora) && !isTRUE(skip_enrich)) {
    tryCatch(run_ora_plots(sub$gene, sub, outdir, title, tag, lfc_col = lfc_col),
             error = function(e) log_msg("ORA failed: ", e$message))
  }
}

analyze_down_list <- function(comp_name, de, heat_mat, sample_info, lfc_col, p_col, do_gsea = TRUE) {
  base <- file.path(result_dir, comp_name)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  write_table(de, file.path(base, "DE_full"))
  writeLines(
    c("先看本目录 DE_full，再看 FoldChange/。",
      "筛选：p < 0.05 且配对肺相对该鼠原位下调（down_FC = 原位/肺 >= 1 与 >= 1.25）。",
      "627 只和 627 比，628 只和 628 比，两鼠不混成一组。",
      "GO/Pathway/KEGG 是 ORA；GSEA 在 00_GSEA_all_genes_NOT_FC_or_topN。"),
    file.path(base, "00_READ_ME.txt")
  )
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    keep <- !is.na(de[[lfc_col]]) & de[[lfc_col]] < 0 &
      (2^(-de[[lfc_col]]) >= fc) & !is.na(de[[p_col]]) & de[[p_col]] < p_cutoff
    sub <- de[keep, , drop = FALSE]
    sub <- sub[order(sub[[lfc_col]], decreasing = FALSE), , drop = FALSE]
    emit_subset(comp_name, sub, nm, paste0(comp_name, " | down FC >= ", fc),
                file.path(base, "FoldChange", nm), de, heat_mat, sample_info,
                lfc_col, p_col, fc_line = fc,
                do_ora = identical(nm, "FC_1.25") && nrow(sub) <= 800)
  }
  if (isTRUE(do_gsea) && !isTRUE(skip_enrich)) {
    tryCatch(run_gsea_full(de, file.path(base, "00_GSEA_all_genes_NOT_FC_or_topN"), comp_name, lfc_col),
             error = function(e) log_msg("GSEA failed: ", e$message))
  }
}

is_down <- function(lfc, p, fc) {
  !is.na(lfc) & !is.na(p) & p < p_cutoff & lfc < 0 & (2^(-lfc) >= fc)
}

# -----------------------------------------------------------------------------
# 5. 配对 limma：每只鼠 原位 vs 自己的肺
# -----------------------------------------------------------------------------
limma_pair <- function(group_a, group_b, label) {
  a <- which(group_a)
  b <- which(group_b)
  if (length(a) < 8 || length(b) < 8) {
    log_msg("Too few spots for ", label, " (n_a=", length(a), " n_b=", length(b), ")")
    return(NULL)
  }
  use <- c(a, b)
  y <- log_mat[, use, drop = FALSE]
  grp <- factor(c(rep("A", length(a)), rep("B", length(b))), levels = c("A", "B"))
  design <- stats::model.matrix(~ 0 + grp)
  colnames(design) <- c("A", "B")
  fit <- limma::eBayes(limma::lmFit(y, design), trend = TRUE, robust = TRUE)
  cont <- limma::makeContrasts(B_vs_A = B - A, levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont), trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "none")
  data.frame(
    gene = rownames(tt),
    log2FC = tt$logFC,
    AveExpr = tt$AveExpr,
    t = tt$t,
    pvalue = tt$P.Value,
    padj = tt$adj.P.Val,
    n_primary = length(a),
    n_lung = length(b),
    down_FC = 2^(-tt$logFC),
    stringsAsFactors = FALSE
  )
}

anno_primary <- function(de) {
  if (is.null(de) || nrow(de) == 0) return(de)
  de$ensembl <- id_map$ensembl[match(de$gene, id_map$gene)]
  de$entrez <- id_map$entrez[match(de$gene, id_map$gene)]
  de$TM_mean <- pb_tm[de$gene]
  de$TM_expr_percentile <- tm_pct[de$gene]
  de$TM_low_quartile <- !is.na(de$TM_expr_percentile) & de$TM_expr_percentile <= 25
  de
}

mice <- intersect(c("Mouse627", "Mouse628"), unique(meta$mouse))
if (length(mice) == 0) stop("未识别到 Mouse627 / Mouse628")

de_all_by_mouse <- list()
de_foci_by_mouse <- list()
de_uninv_by_mouse <- list()

for (ms in mice) {
  log_msg("Paired DE: ", ms, " LUNG vs matched TM (all in-tissue)")
  de_all_by_mouse[[ms]] <- anno_primary(limma_pair(
    meta$mouse == ms & meta$tissue == "TM",
    meta$mouse == ms & meta$tissue == "LUNG",
    paste(ms, "LUNG_all vs TM_all")
  ))
  n_foci <- sum(meta$mouse == ms & meta$compartment == "LUNG_foci")
  n_tm_t <- sum(meta$mouse == ms & meta$compartment == "TM_tumor")
  if (n_foci >= 8 && n_tm_t >= 8) {
    log_msg("Paired DE: ", ms, " LUNG_foci vs TM_tumor")
    de_foci_by_mouse[[ms]] <- anno_primary(limma_pair(
      meta$mouse == ms & meta$compartment == "TM_tumor",
      meta$mouse == ms & meta$compartment == "LUNG_foci",
      paste(ms, "foci vs TM_tumor")
    ))
  } else {
    log_msg(ms, " too few foci (", n_foci, "); reuse all-lung DE for foci table")
    de_foci_by_mouse[[ms]] <- de_all_by_mouse[[ms]]
  }
  n_un <- sum(meta$mouse == ms & meta$compartment == "LUNG_uninvolved")
  if (n_un >= 8 && n_tm_t >= 8) {
    de_uninv_by_mouse[[ms]] <- anno_primary(limma_pair(
      meta$mouse == ms & meta$compartment == "TM_tumor",
      meta$mouse == ms & meta$compartment == "LUNG_uninvolved",
      paste(ms, "uninvolved vs TM_tumor")
    ))
  } else {
    de_uninv_by_mouse[[ms]] <- NULL
  }
}

merge_paired <- function(de_list, col_prefix = "") {
  stopifnot(length(de_list) >= 1)
  base <- de_list[[1]][, c("gene", "ensembl", "entrez", "TM_mean", "TM_expr_percentile", "TM_low_quartile")]
  for (ms in names(de_list)) {
    d <- de_list[[ms]]
    if (is.null(d)) next
    stub <- paste0(col_prefix, ms)
    add <- data.frame(
      gene = d$gene,
      log2FC = d$log2FC,
      pvalue = d$pvalue,
      padj = d$padj,
      down_FC = d$down_FC,
      AveExpr = d$AveExpr
    )
    names(add)[-1] <- paste0(names(add)[-1], "_", stub)
    base <- merge(base, add, by = "gene", all = TRUE)
  }
  ms_ok <- names(de_list)[!vapply(de_list, is.null, logical(1))]
  lfc_cols <- paste0("log2FC_", col_prefix, ms_ok)
  p_cols <- paste0("pvalue_", col_prefix, ms_ok)
  lfc_cols <- intersect(lfc_cols, names(base))
  p_cols <- intersect(p_cols, names(base))
  if (length(lfc_cols) == 2) {
    base$mean_log2FC <- rowMeans(base[, lfc_cols, drop = FALSE], na.rm = TRUE)
    base$pvalue_both <- pmax(base[[p_cols[1]]], base[[p_cols[2]]], na.rm = FALSE)
  } else {
    base$mean_log2FC <- base[[lfc_cols[1]]]
    base$pvalue_both <- base[[p_cols[1]]]
  }
  base$down_FC_mean <- 2^(-base$mean_log2FC)
  base
}

merged_all <- merge_paired(de_all_by_mouse)
merged_foci <- merge_paired(de_foci_by_mouse, col_prefix = "foci_")
write_table(merged_all, file.path(qc_dir, "paired_limma_LUNG_vs_TM_all_spots"))
write_table(merged_foci, file.path(qc_dir, "paired_limma_LUNG_foci_vs_TM_tumor"))

# -----------------------------------------------------------------------------
# 6. 问题 1：配对肺下调（每鼠单独 + 两鼠都下调）
# -----------------------------------------------------------------------------
log_msg("Q1: paired lung down vs matched primary")

for (ms in names(de_all_by_mouse)) {
  if (is.null(de_all_by_mouse[[ms]])) next
  analyze_down_list(
    paste0("01_lung_down_", ms),
    de_all_by_mouse[[ms]], pb, sample_info,
    lfc_col = "log2FC", p_col = "pvalue",
    do_gsea = identical(ms, names(de_all_by_mouse)[1])
  )
}

both_de <- merged_all
both_de$log2FC <- both_de$mean_log2FC
both_de$pvalue <- both_de$pvalue_both
keep_both <- TRUE
if (length(mice) == 2) {
  c1 <- paste0("log2FC_", mice[1])
  c2 <- paste0("log2FC_", mice[2])
  p1 <- paste0("pvalue_", mice[1])
  p2 <- paste0("pvalue_", mice[2])
  keep_both <- is_down(both_de[[c1]], both_de[[p1]], 1) &
    is_down(both_de[[c2]], both_de[[p2]], 1)
  both_de$pvalue[!keep_both] <- 1
}
analyze_down_list("01_lung_down_both_mice", both_de, pb, sample_info,
                  lfc_col = "mean_log2FC", p_col = "pvalue")

base_both <- file.path(result_dir, "01_lung_down_both_mice")
if (length(mice) == 2) {
  c1 <- paste0("log2FC_", mice[1])
  c2 <- paste0("log2FC_", mice[2])
  p1 <- paste0("pvalue_", mice[1])
  p2 <- paste0("pvalue_", mice[2])
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    sub <- merged_all[is_down(merged_all[[c1]], merged_all[[p1]], fc) &
                        is_down(merged_all[[c2]], merged_all[[p2]], fc), ]
    sub <- sub[order(sub$mean_log2FC), ]
    sub$log2FC <- sub$mean_log2FC
    sub$pvalue <- pmax(sub[[p1]], sub[[p2]])
    emit_subset(
      "01_lung_down_both_mice", sub, paste0(nm, "_BOTH_mice"),
      paste0("Q1 both mice down | FC >= ", fc),
      file.path(base_both, "FoldChange", paste0(nm, "_BOTH_mice")),
      transform(merged_all, log2FC = mean_log2FC, pvalue = pvalue_both),
      pb, sample_info, "mean_log2FC", "pvalue_both", fc,
      do_ora = identical(nm, "FC_1.25") && nrow(sub) <= 800
    )
  }
}

bone_dir <- file.path(result_dir, "03_bone_not_in_GSE273439")
dir.create(bone_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  c("GSE273439 没有骨 / 骨髓切片。",
    "四个 visium 样品只有：Mouse627_TM、Mouse627_LUNG、Mouse628_TM、Mouse628_LUNG。",
    "因此不能从本数据寻找「原位低表达促进骨转移」，也不能做肺 vs 骨的器官特异交集。",
    "骨转移请用 GSE165393（骨髓来源细胞系）或 GSE37975（4T1.2 配对芯片）。",
    "不要把两只鼠的肺混在一起冒充多器官。"),
  file.path(bone_dir, "00_NOTE.txt")
)
log_msg("Q1 bone: skipped (no bone samples in GSE273439)")

# -----------------------------------------------------------------------------
# 7. 问题 2：肺转移灶特异（相对未受累肺），不是骨特异
# -----------------------------------------------------------------------------
log_msg("Q2: lung-met-foci specific down (not uninvolved lung)")
foci_de <- merged_foci
foci_de$log2FC <- foci_de$mean_log2FC
foci_de$pvalue <- foci_de$pvalue_both

if (length(mice) == 2 && all(paste0("log2FC_foci_", mice) %in% names(merged_foci))) {
  f1 <- paste0("log2FC_foci_", mice[1])
  f2 <- paste0("log2FC_foci_", mice[2])
  pf1 <- paste0("pvalue_foci_", mice[1])
  pf2 <- paste0("pvalue_foci_", mice[2])
  mark_foci <- is_down(merged_foci[[f1]], merged_foci[[pf1]], 1) &
    is_down(merged_foci[[f2]], merged_foci[[pf2]], 1)
} else {
  mark_foci <- is_down(foci_de$mean_log2FC, foci_de$pvalue_both, 1)
}

# 未受累肺也下调的基因 = 肺组织差异，不算转移灶特异
uninv_down <- character()
if (length(de_uninv_by_mouse) > 0) {
  un_ok <- names(de_uninv_by_mouse)[!vapply(de_uninv_by_mouse, is.null, logical(1))]
  if (length(un_ok) > 0) {
    merged_un <- merge_paired(de_uninv_by_mouse[un_ok], col_prefix = "uninv_")
    write_table(merged_un, file.path(qc_dir, "paired_limma_uninvolved_lung_vs_TM"))
    if (length(un_ok) == 2 && all(paste0("log2FC_uninv_", un_ok) %in% names(merged_un))) {
      uninv_down <- merged_un$gene[
        is_down(merged_un[[paste0("log2FC_uninv_", un_ok[1])]],
                merged_un[[paste0("pvalue_uninv_", un_ok[1])]], 1) &
          is_down(merged_un[[paste0("log2FC_uninv_", un_ok[2])]],
                  merged_un[[paste0("pvalue_uninv_", un_ok[2])]], 1)
      ]
    } else {
      uninv_down <- merged_un$gene[is_down(merged_un$mean_log2FC, merged_un$pvalue_both, 1)]
    }
  }
}

foci_spec <- foci_de
mark_spec <- mark_foci & !(foci_spec$gene %in% uninv_down)
foci_spec$pvalue[!mark_spec] <- 1
analyze_down_list("02_lung_foci_specific_down", foci_spec, pb, sample_info,
                  lfc_col = "mean_log2FC", p_col = "pvalue")

if (has_pkg("ggvenn") && length(mice) == 2) {
  vdf <- list()
  vdf[[paste0(mice[1], "_lung_down")]] <- merged_all$gene[is_down(merged_all[[paste0("log2FC_", mice[1])]], merged_all[[paste0("pvalue_", mice[1])]], 1)]
  vdf[[paste0(mice[2], "_lung_down")]] <- merged_all$gene[is_down(merged_all[[paste0("log2FC_", mice[2])]], merged_all[[paste0("pvalue_", mice[2])]], 1)]
  p <- ggvenn::ggvenn(vdf, fill_color = c("#F58518", "#54A24B")) +
    ggplot2::labs(title = "Paired lung-down genes (p<0.05, FC>=1) per mouse")
  save_gg(p, file.path(qc_dir, "venn_two_mice_lung_down"))
}

# -----------------------------------------------------------------------------
# 8. 问题 3：原位瘤三套神经浸润分数 + 负相关基因
# -----------------------------------------------------------------------------
log_msg("Q3: three neural signatures on primary (TM) spots")
pni_dir <- file.path(result_dir, "06_neural_invasion")
dir.create(pni_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  c("没有病理神经浸润标签。只在原位 TM spots 上打三种分数：",
    "  1) 施旺细胞 Schwann",
    "  2) 神经营养因子 Neurotrophin",
    "  3) 轴突导向 Axon_guidance",
    "负相关：该鼠原位 spots 上 Spearman rho<0 且 p<0.05。",
    "两鼠各自计算，主结论取两鼠都负相关的基因。"),
  file.path(pni_dir, "00_READ_ME.txt")
)

msig_bp <- tryCatch({
  tryCatch(
    msigdbr::msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP"),
    error = function(e) msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
  )
}, error = function(e) {
  log_msg("msigdbr GO:BP failed: ", e$message)
  NULL
})

msig_symbols <- function(pattern) {
  if (is.null(msig_bp) || nrow(msig_bp) == 0) return(character())
  col <- intersect(c("gene_symbol", "gs_exact_source"), names(msig_bp))[1]
  if (!"gene_symbol" %in% names(msig_bp)) return(character())
  unique(msig_bp$gene_symbol[grepl(pattern, msig_bp$gs_name)])
}

curated <- list(
  Schwann = c(
    "Sox10", "S100b", "Mag", "Mpz", "Mbp", "Pmp22", "Plp1", "Gfap", "Egr2",
    "Pou3f1", "Pou3f2", "Ngfr", "Ncam1", "L1cam", "Cdh19", "Dhh", "Mpzl1",
    "Prx", "Gjb1", "Mal", "Cnp", "Erbb3", "Erbb2", "Nrg1"
  ),
  Neurotrophin = c(
    "Ngf", "Bdnf", "Ntf3", "Ntf5", "Ntrk1", "Ntrk2", "Ntrk3", "Ngfr",
    "Artn", "Gdnf", "Nrtn", "Pspn", "Gfra1", "Gfra2", "Gfra3", "Ret",
    "Gal", "Tac1", "Calca", "Vgf", "Cntf", "Lif", "Igf1"
  ),
  Axon_guidance = c(
    "Sema3a", "Sema3d", "Sema3e", "Sema3f", "Sema4d", "Sema5a", "Sema6a",
    "Nrp1", "Nrp2", "Plxna1", "Plxna3", "Plxnb1", "Robo1", "Robo2",
    "Slit1", "Slit2", "Slit3", "Ephb2", "Efnb1", "Epha4", "Efna1",
    "Unc5b", "Unc5c", "Dcc", "Ntn1", "L1cam", "Ncam1", "Efnb2", "Ephb1"
  )
)
go_pat <- list(
  Schwann = "SCHWANN|MYELIN|PERIPHERAL_NERVOUS",
  Neurotrophin = "NEUROTROPHIN|NERVE_GROWTH_FACTOR|GDNF|NEUROTROPHIC",
  Axon_guidance = "AXON_GUIDANCE|NEURON_PROJECTION_GUIDANCE|SEMASPHORIN|NETRIN|SLIT_ROBO"
)

sig_genes <- lapply(names(curated), function(nm) {
  g <- unique(c(curated[[nm]], msig_symbols(go_pat[[nm]])))
  intersect(g, rownames(log_mat))
})
names(sig_genes) <- names(curated)
for (nm in names(sig_genes)) {
  log_msg("Signature ", nm, ": ", length(sig_genes[[nm]]), " genes")
  utils::write.csv(
    data.frame(gene = sig_genes[[nm]]),
    file.path(pni_dir, paste0("signature_genes_", nm, ".csv")),
    row.names = FALSE
  )
}

mean_z_spots <- function(genes, cols) {
  genes <- intersect(genes, rownames(log_mat))
  if (length(genes) < 2) return(setNames(rep(NA_real_, length(cols)), cols))
  z <- t(scale(t(log_mat[genes, cols, drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

for (nm in names(sig_genes)) {
  meta[[nm]] <- mean_z_spots(sig_genes[[nm]], colnames(log_mat))
}

spot_score_tab <- meta[, c("cell", "sample", "mouse", "tissue", "compartment",
                           "tumor_score", names(curated)), drop = FALSE]
write_table(spot_score_tab, file.path(pni_dir, "spot_three_neural_scores"))

tm_mark <- dplyr::bind_rows(lapply(mice, function(ms) {
  sub <- meta[meta$mouse == ms & meta$tissue == "TM", ]
  data.frame(
    mouse = ms, n_TM_spots = nrow(sub),
    Schwann_mean = mean(sub$Schwann, na.rm = TRUE),
    Neurotrophin_mean = mean(sub$Neurotrophin, na.rm = TRUE),
    Axon_guidance_mean = mean(sub$Axon_guidance, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write_table(tm_mark, file.path(pni_dir, "TM_neural_invasion_marker_by_mouse"))

score_long <- tidyr::pivot_longer(
  meta[, c("sample", "mouse", "tissue", names(curated))],
  cols = dplyr::all_of(names(curated)),
  names_to = "signature", values_to = "score"
)
p_box <- ggplot2::ggplot(score_long, ggplot2::aes(x = tissue, y = score, fill = tissue)) +
  ggplot2::geom_violin(scale = "width", alpha = 0.7) +
  ggplot2::facet_grid(signature ~ mouse, scales = "free_y") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(title = "Neural scores: matched TM vs LUNG (spots)", x = NULL, y = "score")
save_gg(p_box, file.path(pni_dir, "three_scores_TM_vs_LUNG_by_mouse"), 9, 8)

for (nm in names(curated)) {
  plot_spatial(meta, nm, paste(nm, "score"), file.path(pni_dir, paste0("spatial_", nm)))
}

spearman_vs_score_fast <- function(score_vec, cols) {
  sc <- as.numeric(score_vec[cols])
  ok <- is.finite(sc)
  cols <- cols[ok]
  sc <- sc[ok]
  if (length(sc) < 8) {
    return(data.frame(gene = rownames(log_mat), rho = NA_real_, pvalue = NA_real_))
  }
  ranked_sc <- rank(sc)
  n <- length(sc)
  mat <- log_mat[, cols, drop = FALSE]
  rho <- apply(mat, 1, function(x) {
    if (sum(is.finite(x)) < 8) return(NA_real_)
    stats::cor(rank(x), ranked_sc, method = "pearson", use = "complete.obs")
  })
  tstat <- rho * sqrt((n - 2) / pmax(1 - rho^2, 1e-8))
  pval <- 2 * stats::pt(-abs(tstat), df = n - 2)
  data.frame(gene = names(rho), rho = as.numeric(rho), pvalue = as.numeric(pval), stringsAsFactors = FALSE)
}

neg_by_sig_mouse <- list()
neg_both_by_sig <- list()

for (nm in names(curated)) {
  log_msg("Neural axis: ", nm)
  ax_dir <- file.path(pni_dir, nm)
  dir.create(ax_dir, recursive = TRUE, showWarnings = FALSE)
  plot_heatmap(
    pb, sample_info, sig_genes[[nm]],
    paste(nm, "signature genes (pseudobulk)"),
    file.path(ax_dir, paste0("heatmap_", nm, "_signature"))
  )

  mouse_neg <- list()
  mouse_cor <- list()
  for (ms in mice) {
    cols <- meta$cell[meta$mouse == ms & meta$tissue == "TM"]
    cor_df <- spearman_vs_score_fast(setNames(meta[[nm]], meta$cell), cols)
    cor_df$mouse <- ms
    cor_df$padj <- stats::p.adjust(cor_df$pvalue, method = "BH")
    cor_df$ensembl <- id_map$ensembl[match(cor_df$gene, id_map$gene)]
    cor_df$TM_mean <- pb_tm[cor_df$gene]
    mouse_cor[[ms]] <- cor_df
    neg <- cor_df[is.finite(cor_df$rho) & cor_df$rho < 0 & cor_df$pvalue < p_cutoff, ]
    neg <- neg[order(neg$rho), ]
    neg$log2FC <- neg$rho
    mouse_neg[[ms]] <- neg
    write_table(cor_df, file.path(ax_dir, paste0(nm, "_", ms, "_gene_vs_TM_score_spearman")))
    write_table(neg, file.path(ax_dir, paste0(nm, "_", ms, "_neg_cor_p0.05")))
    emit_subset(
      paste0("06_neural_invasion/", nm, "/", ms), neg,
      paste0(nm, "_", ms, "_neg_cor_p0.05"),
      paste0(nm, " | ", ms, " TM spots | rho<0 p<0.05"),
      file.path(ax_dir, ms, "neg_cor_p0.05"),
      transform(cor_df, log2FC = rho), pb, sample_info,
      lfc_col = "rho", p_col = "pvalue",
      do_ora = nrow(neg) <= 800 && !isTRUE(skip_enrich)
    )
  }
  neg_by_sig_mouse[[nm]] <- mouse_neg

  if (length(mice) == 2) {
    g1 <- mouse_neg[[mice[1]]]$gene
    g2 <- mouse_neg[[mice[2]]]$gene
    both_g <- intersect(g1, g2)
    both <- mouse_neg[[mice[1]]][mouse_neg[[mice[1]]]$gene %in% both_g, ]
    both$rho_other <- mouse_neg[[mice[2]]]$rho[match(both$gene, mouse_neg[[mice[2]]]$gene)]
    both$pvalue_other <- mouse_neg[[mice[2]]]$pvalue[match(both$gene, mouse_neg[[mice[2]]]$gene)]
    both$rho_mean <- (both$rho + both$rho_other) / 2
    both$log2FC <- both$rho_mean
    both$pvalue <- pmax(both$pvalue, both$pvalue_other)
    both <- both[order(both$rho_mean), ]
  } else {
    both <- mouse_neg[[mice[1]]]
    both$rho_mean <- both$rho
  }
  neg_both_by_sig[[nm]] <- both
  write_table(both, file.path(ax_dir, paste0(nm, "_BOTH_mice_neg_cor_p0.05")))
  emit_subset(
    paste0("06_neural_invasion/", nm), both,
    paste0(nm, "_BOTH_mice_neg_cor"),
    paste0(nm, " | both mice TM negative correlation"),
    file.path(ax_dir, "BOTH_mice_neg_cor"),
    transform(mouse_cor[[mice[1]]], log2FC = rho), pb, sample_info,
    lfc_col = "rho", p_col = "pvalue",
    do_ora = nrow(both) <= 800 && !isTRUE(skip_enrich)
  )
  if (!isTRUE(skip_enrich)) {
    tryCatch(
      run_gsea_full(
        transform(mouse_cor[[mice[1]]], log2FC = rho),
        file.path(ax_dir, "00_GSEA_all_genes_NOT_FC_or_topN"),
        paste(nm, "TM correlation"), stat_col = "rho"
      ),
      error = function(e) log_msg(nm, " GSEA failed: ", e$message)
    )
  }
}

if (has_pkg("ggvenn")) {
  v3 <- lapply(neg_both_by_sig, function(x) x$gene)
  p <- ggvenn::ggvenn(v3, fill_color = c("#4C78A8", "#F58518", "#54A24B")) +
    ggplot2::labs(title = "Both-mice TM negative-correlation genes")
  save_gg(p, file.path(pni_dir, "venn_three_neg_cor_both_mice"), 8, 6)
}

# -----------------------------------------------------------------------------
# 9. 问题 4：三种神经分数分别 vs 配对肺转移
# -----------------------------------------------------------------------------
log_msg("Q4: each neural score vs matched lung metastasis")
q4_dir <- file.path(result_dir, "07_neural_vs_lung")
dir.create(q4_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  c("每种神经分数单独对配对肺：",
    "  - 该鼠 TM spots vs 该鼠 LUNG spots 的 Wilcoxon（分数高低）",
    "  - 该鼠 TM spots vs 该鼠 LUNG_foci",
    "  - 负相关基因 ∩ 该鼠肺下调基因",
    "  - 两鼠都成立的交集"),
  file.path(q4_dir, "00_READ_ME.txt")
)

wilcox_score <- function(a, b) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) < 5 || length(b) < 5) {
    return(data.frame(n_a = length(a), n_b = length(b), mean_a = mean(a), mean_b = mean(b),
                      delta = mean(b) - mean(a), pvalue = NA_real_))
  }
  wt <- suppressWarnings(stats::wilcox.test(b, a))
  data.frame(
    n_a = length(a), n_b = length(b),
    mean_a = mean(a), mean_b = mean(b),
    delta = mean(b) - mean(a),
    pvalue = wt$p.value
  )
}

lung_down_mouse <- list()
if (length(mice) >= 1) {
  for (ms in mice) {
    lung_down_mouse[[ms]] <- merged_all$gene[is_down(
      merged_all[[paste0("log2FC_", ms)]],
      merged_all[[paste0("pvalue_", ms)]], 1
    )]
  }
}
lung_down_both <- if (length(mice) == 2) {
  intersect(lung_down_mouse[[mice[1]]], lung_down_mouse[[mice[2]]])
} else {
  lung_down_mouse[[mice[1]]]
}

assoc_rows <- list()
overlap_rows <- list()
prog_rows <- list()

for (nm in names(curated)) {
  ax_dir <- file.path(q4_dir, nm)
  dir.create(ax_dir, recursive = TRUE, showWarnings = FALSE)
  for (ms in mice) {
    tm_sc <- meta[[nm]][meta$mouse == ms & meta$tissue == "TM"]
    lu_sc <- meta[[nm]][meta$mouse == ms & meta$tissue == "LUNG"]
    foci_sc <- meta[[nm]][meta$mouse == ms & meta$compartment == "LUNG_foci"]
    tt_lu <- wilcox_score(tm_sc, lu_sc)
    tt_foci <- wilcox_score(tm_sc, foci_sc)
    assoc_rows[[paste(nm, ms, "LUNG")]] <- data.frame(
      signature = nm, mouse = ms, contrast = "LUNG_vs_matched_TM", tt_lu
    )
    assoc_rows[[paste(nm, ms, "foci")]] <- data.frame(
      signature = nm, mouse = ms, contrast = "LUNG_foci_vs_matched_TM", tt_foci
    )

    neg <- neg_by_sig_mouse[[nm]][[ms]]
    ov <- neg[neg$gene %in% lung_down_mouse[[ms]], ]
    ov$log2FC_lung <- merged_all[[paste0("log2FC_", ms)]][match(ov$gene, merged_all$gene)]
    ov$pvalue_lung <- merged_all[[paste0("pvalue_", ms)]][match(ov$gene, merged_all$gene)]
    write_table(ov, file.path(ax_dir, paste0(nm, "_", ms, "_neg_cor_AND_paired_lung_down")))
    overlap_rows[[paste(nm, ms)]] <- data.frame(
      signature = nm, mouse = ms,
      n_neg_cor = nrow(neg),
      n_lung_down = length(lung_down_mouse[[ms]]),
      n_overlap = nrow(ov),
      stringsAsFactors = FALSE
    )
  }
  both_neg <- neg_both_by_sig[[nm]]
  both_ov <- both_neg[both_neg$gene %in% lung_down_both, ]
  both_ov$mean_log2FC_lung <- merged_all$mean_log2FC[match(both_ov$gene, merged_all$gene)]
  write_table(both_ov, file.path(ax_dir, paste0(nm, "_BOTH_mice_neg_cor_AND_both_lung_down")))

  # 原位高神经分数 spots 上，肺下调基因是否更低
  for (ms in mice) {
    tm <- meta$mouse == ms & meta$tissue == "TM"
    sc <- meta[[nm]][tm]
    hi <- meta$cell[tm][sc >= stats::median(sc, na.rm = TRUE)]
    lo <- meta$cell[tm][sc < stats::median(sc, na.rm = TRUE)]
    guse <- intersect(lung_down_mouse[[ms]], rownames(log_mat))
    if (length(guse) >= 5 && length(hi) >= 10 && length(lo) >= 10) {
      prog_hi <- colMeans(log_mat[guse, hi, drop = FALSE], na.rm = TRUE)
      prog_lo <- colMeans(log_mat[guse, lo, drop = FALSE], na.rm = TRUE)
      wt <- suppressWarnings(stats::wilcox.test(prog_hi, prog_lo))
      prog_rows[[paste(nm, ms)]] <- data.frame(
        signature = nm, mouse = ms,
        n_lung_down_genes = length(guse),
        lung_down_program_highNeural_minus_low = mean(prog_hi) - mean(prog_lo),
        p_program = wt$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
}

assoc_all <- do.call(rbind, assoc_rows)
write_table(assoc_all, file.path(q4_dir, "three_signatures_score_vs_matched_lung"))
ov_tab <- dplyr::bind_rows(overlap_rows)
write_table(ov_tab, file.path(q4_dir, "neg_cor_overlap_counts"))
if (length(prog_rows) > 0) {
  write_table(dplyr::bind_rows(prog_rows), file.path(q4_dir, "lung_down_program_in_high_vs_low_neural_TM"))
}

p_assoc <- ggplot2::ggplot(
  assoc_all[assoc_all$contrast == "LUNG_vs_matched_TM", ],
  ggplot2::aes(x = signature, y = delta, fill = mouse)
) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.65) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    title = "Neural score: matched lung minus primary (spots)",
    y = "mean(LUNG) - mean(TM)", x = NULL
  )
save_gg(p_assoc, file.path(q4_dir, "delta_neural_score_lung_vs_TM"), 8, 5)

p_foci <- ggplot2::ggplot(
  assoc_all[assoc_all$contrast == "LUNG_foci_vs_matched_TM", ],
  ggplot2::aes(x = signature, y = delta, fill = mouse)
) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.65) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    title = "Neural score: lung met foci minus matched primary",
    y = "mean(foci) - mean(TM)", x = NULL
  )
save_gg(p_foci, file.path(q4_dir, "delta_neural_score_foci_vs_TM"), 8, 5)

# -----------------------------------------------------------------------------
# 10. 总表
# -----------------------------------------------------------------------------
sum_dir <- file.path(result_dir, "08_summary")
dir.create(sum_dir, recursive = TRUE, showWarnings = FALSE)

q1_list <- list()
if (length(mice) == 2) {
  c1 <- paste0("log2FC_", mice[1])
  c2 <- paste0("log2FC_", mice[2])
  p1 <- paste0("pvalue_", mice[1])
  p2 <- paste0("pvalue_", mice[2])
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    sub <- merged_all[is_down(merged_all[[c1]], merged_all[[p1]], fc) &
                        is_down(merged_all[[c2]], merged_all[[p2]], fc), ]
    sub <- sub[order(sub$mean_log2FC), ]
    q1_list[[nm]] <- sub
    write_table(sub, file.path(sum_dir, paste0("Q1_both_mice_lung_down_p0.05_", nm)))
  }
  for (ms in mice) {
    sub <- merged_all[is_down(merged_all[[paste0("log2FC_", ms)]],
                              merged_all[[paste0("pvalue_", ms)]], 1.25), ]
    sub <- sub[order(sub[[paste0("log2FC_", ms)]]), ]
    write_table(sub, file.path(sum_dir, paste0("Q1_", ms, "_lung_down_p0.05_FC1.25")))
  }
} else {
  ms <- mice[1]
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    sub <- de_all_by_mouse[[ms]][is_down(de_all_by_mouse[[ms]]$log2FC, de_all_by_mouse[[ms]]$pvalue, fc), ]
    q1_list[[nm]] <- sub
    write_table(sub, file.path(sum_dir, paste0("Q1_", ms, "_lung_down_p0.05_", nm)))
  }
}

q2 <- foci_spec[is_down(foci_spec$mean_log2FC, foci_spec$pvalue, 1.25), ]
q2 <- q2[order(q2$mean_log2FC), ]
write_table(q2, file.path(sum_dir, "Q2_lung_foci_specific_down_p0.05_FC1.25"))
writeLines(
  "GSE273439 has no bone. Q2 bone-specific list does not exist.",
  file.path(sum_dir, "Q2_bone_specific_NOT_AVAILABLE.txt")
)

write_table(tm_mark, file.path(sum_dir, "Q3_TM_three_neural_scores"))
write_table(assoc_all, file.path(sum_dir, "Q4_neural_scores_vs_matched_lung"))
write_table(ov_tab, file.path(sum_dir, "Q4_neg_cor_overlap_with_lung_down"))

for (nm in names(neg_both_by_sig)) {
  write_table(utils::head(neg_both_by_sig[[nm]], 500),
              file.path(sum_dir, paste0("Q3_", nm, "_BOTH_mice_neg_cor_p0.05_top500")))
}

n_tab <- data.frame(
  question = c(
    paste0("Q1 both-mice lung down p<0.05 ", names(fc_cutoffs)),
    "Q2 lung-foci-specific down p<0.05 FC>=1.25",
    "Q2 bone-specific",
    paste0("Q3 ", names(neg_both_by_sig), " both-mice TM neg-cor p<0.05")
  ),
  n_genes = c(
    vapply(q1_list, nrow, integer(1)),
    nrow(q2),
    0,
    vapply(neg_both_by_sig, nrow, integer(1))
  ),
  stringsAsFactors = FALSE
)
write_table(n_tab, file.path(sum_dir, "gene_counts"))
log_msg("Done. Summary counts:")
print(n_tab)
log_msg("Open: ", sum_dir)
log_msg("Bone metastasis cannot be analysed in GSE273439 (no bone slides).")
