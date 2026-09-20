#!/usr/bin/env Rscript
# =============================================================================
# GSE146012_Mouse_breast.R
# 4T1 GFP+ 肿瘤细胞 bulk RNA-seq：原发瘤 vs 肺转移（GSE146012）
#
# 用法（Windows）：
#   setwd("E:/R/Mouse Breast/GSE146012")
#   source("GSE146012_Mouse_breast.R")
#
# 四个问题（「低表达促进转移」= 配对肺相对该号原位下调）：
#   1) 原位低表达与肺转移相关的下调基因（1 号原位对 1 号肺，2 号对 2 号）
#      骨转移：本 GEO 没有骨，目录 03_ 写明，不编造
#   2) 肺路径特异下调（MFP 自发肺转移下调、TVI 实验性肺定植不下调）；不能做骨特异
#   3) 原位瘤按施旺 / 神经营养 / 轴突导向打分，找负相关基因
#   4) 三种神经分数分别与配对肺转移的关系
#
# 入选：p < 0.05，下调倍数 down_FC = 原位/肺 >= 1 与 >= 1.25
# 不做 Top50–300。
#
# 数据注意：
#   - 已是 DESeq2 log 矩阵，不要再跑 DESeq2
#   - 原发只有 Tumor-1/2；MFP-Met-3、TVI-Met-3 没有对应原发 3 号
#   - 1-vs-1 无法估计 p，不伪造；p 值来自两对配对 limma
#   - MFP-Met 与 TVI-Met 不合并
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")
options(clusterProfiler.download.method = "auto")
if (capabilities("cairo")) options(bitmapType = "cairo")

cran_required <- c(
  "dplyr", "tidyr", "tibble", "stringr", "ggplot2", "ggrepel", "pheatmap",
  "RColorBrewer", "matrixStats", "cowplot", "writexl"
)
cran_optional <- c("ggvenn", "GSVA")
skip_enrich <- nzchar(Sys.getenv("GSE146012_SKIP_ENRICHMENT"))
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
    message("GSE146012_SKIP_ENRICHMENT: skip installing ", paste(miss_opt, collapse = ", "))
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
  env_dir <- Sys.getenv("GSE146012_DIR", unset = "")
  candidates <- c(
    env_dir,
    "E:/R/Mouse Breast/GSE146012",
    "E:\\R\\Mouse Breast\\GSE146012",
    "E:/R/Mouse Breast",
    file.path(getwd(), "GSE146012"),
    getwd()
  )
  candidates <- unique(candidates[nzchar(candidates)])
  looks_like <- function(d) {
    if (!dir.exists(d)) return(FALSE)
    hits <- list.files(d, pattern = "GSE146012.*RNA_seq_logtransformed_count",
                       recursive = TRUE, ignore.case = TRUE)
    length(hits) > 0
  }
  for (d in candidates) {
    if (looks_like(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

find_count_file <- function(root) {
  hits <- list.files(
    root, pattern = "GSE146012.*RNA_seq_logtransformed_count.*\\.(txt|csv)(\\.gz)?$",
    recursive = TRUE, full.names = TRUE, ignore.case = TRUE
  )
  hits <- hits[!grepl("results_GSE146012", hits)]
  if (length(hits) == 0) {
    stop("未找到 GSE146012_RNA_seq_logtransformed_count.txt.gz。请放到: ", root)
  }
  hits[[1]]
}

project_dir <- resolve_project_dir()
result_dir  <- file.path(project_dir, "results_GSE146012")
log_dir     <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("GSE146012_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

p_cutoff   <- 0.05
# 只这两档；不要 Top50–300
fc_cutoffs <- c("FC_1" = 1, "FC_1.25" = 1.25)

log_msg("Project dir: ", project_dir)
log_msg("Results dir: ", result_dir)
log_msg("Gate: p < ", p_cutoff, " ; down-FC (primary/lung) >= 1 and >= 1.25 ; no Top50-300")
log_msg("Paired by replicate number: Tumor-1 vs Met-1, Tumor-2 vs Met-2")

# -----------------------------------------------------------------------------
# 2. 读入 log 矩阵
# -----------------------------------------------------------------------------
count_file <- find_count_file(project_dir)
log_msg("Count file: ", count_file)
raw <- utils::read.delim(count_file, check.names = FALSE, stringsAsFactors = FALSE)
names(raw) <- gsub('^\\"|\\"$', "", names(raw))
ens_col <- intersect(c("EnsembleID", "EnsemblID", "ensembl", "gene_id"), names(raw))
sym_col <- intersect(c("GeneSymbol", "gene", "symbol", "Gene"), names(raw))
if (length(ens_col) == 0) ens_col <- names(raw)[1]
if (length(sym_col) == 0) {
  raw$GeneSymbol <- as.character(raw[[ens_col[1]]])
  sym_col <- "GeneSymbol"
}
ens_col <- ens_col[[1]]
sym_col <- sym_col[[1]]
sample_cols <- setdiff(names(raw), c(ens_col, sym_col))
mat <- as.matrix(raw[, sample_cols, drop = FALSE])
storage.mode(mat) <- "numeric"
ensembl <- sub("\\.\\d+$", "", as.character(raw[[ens_col]]))
symbols <- as.character(raw[[sym_col]])
symbols[is.na(symbols) | !nzchar(symbols)] <- ensembl[is.na(symbols) | !nzchar(symbols)]

# 重复符号取均值
if (any(duplicated(symbols))) {
  log_msg("Collapsing duplicate symbols by mean")
  uniq <- unique(symbols)
  idx <- split(seq_along(symbols), factor(symbols, levels = uniq))
  mat <- do.call(rbind, lapply(uniq, function(s) {
    rows <- idx[[s]]
    if (length(rows) == 1L) mat[rows, , drop = FALSE] else {
      m <- matrix(colMeans(mat[rows, , drop = FALSE], na.rm = TRUE), nrow = 1)
      colnames(m) <- colnames(mat)
      m
    }
  }))
  rownames(mat) <- uniq
  ensembl <- vapply(uniq, function(s) ensembl[match(s, symbols)], character(1))
} else {
  rownames(mat) <- symbols
}
colnames(mat) <- gsub("^4T1-", "", colnames(mat))

parse_sample <- function(nm) {
  tissue <- dplyr::case_when(
    grepl("Tumor", nm, ignore.case = TRUE) ~ "Tumor",
    grepl("MFP", nm, ignore.case = TRUE) ~ "MFP_Met",
    grepl("TVI", nm, ignore.case = TRUE) ~ "TVI_Met",
    TRUE ~ "other"
  )
  pair <- as.integer(sub(".*-(\\d+)$", "\\1", nm))
  if (!is.finite(pair)) pair <- NA_integer_
  data.frame(sample = nm, tissue = tissue, pair = pair,
             route = ifelse(tissue == "TVI_Met", "TVI", "MFP"),
             stringsAsFactors = FALSE)
}
sample_info <- do.call(rbind, lapply(colnames(mat), parse_sample))
rownames(sample_info) <- sample_info$sample
log_msg("Samples: ", paste(sample_info$sample, collapse = ", "))
print(sample_info)
n_tumor <- sum(sample_info$tissue == "Tumor")
n_mfp <- sum(sample_info$tissue == "MFP_Met")
if (n_tumor < 1 || n_mfp < 1) {
  stop("未能识别 Tumor / MFP-Met 列。当前样品: ",
       paste(sample_info$sample, collapse = ", "))
}

qc_dir <- file.path(result_dir, "00_QC")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(sample_info, file.path(qc_dir, "sample_info.csv"), row.names = FALSE)

# 低表达过滤（log 矩阵）
keep_gene <- rowSums(is.finite(mat) & mat > 1) >= 2
if (sum(keep_gene) < 200) keep_gene <- rowSums(is.finite(mat) & mat > 0) >= 2
log_msg("Low-expression filter: keep ", sum(keep_gene), " / ", nrow(mat))
mat <- mat[keep_gene, , drop = FALSE]
ensembl <- ensembl[keep_gene]

map_symbols <- function(symbols) {
  uniq <- unique(symbols)
  if (!has_pkg("org.Mm.eg.db") || !has_pkg("AnnotationDbi")) {
    log_msg("org.Mm.eg.db not available; keep gene symbols without Entrez map")
    return(data.frame(gene = uniq, ensembl = NA_character_, entrez = NA_character_,
                      stringsAsFactors = FALSE))
  }
  ent <- tryCatch(
    AnnotationDbi::mapIds(org.Mm.eg.db, keys = uniq, column = "ENTREZID",
                          keytype = "SYMBOL", multiVals = "first"),
    error = function(e) setNames(rep(NA_character_, length(uniq)), uniq)
  )
  ens <- tryCatch(
    AnnotationDbi::mapIds(org.Mm.eg.db, keys = uniq, column = "ENSEMBL",
                          keytype = "SYMBOL", multiVals = "first"),
    error = function(e) setNames(rep(NA_character_, length(uniq)), uniq)
  )
  data.frame(gene = uniq, ensembl = unname(ens[uniq]), entrez = unname(ent[uniq]),
             stringsAsFactors = FALSE)
}
id_map <- map_symbols(rownames(mat))
id_map$ensembl <- ifelse(!is.na(id_map$ensembl) & nzchar(id_map$ensembl),
                         id_map$ensembl, ensembl[match(id_map$gene, rownames(mat))])
id_map <- id_map[match(rownames(mat), id_map$gene), ]
id_map$gene <- rownames(mat)

# 已是 DESeq2 log，不再 DESeq2；热图用原矩阵
heat_mat <- mat
tumor_cols <- sample_info$sample[sample_info$tissue == "Tumor"]
tm_mean <- rowMeans(heat_mat[, tumor_cols, drop = FALSE], na.rm = TRUE)
tm_pct <- rank(tm_mean, ties.method = "average") / length(tm_mean) * 100

save_gg <- function(plot, path_stub, width = 8, height = 6) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".pdf"), plot, width = width, height = height),
           error = function(e) log_msg("pdf ggsave failed: ", e$message))
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".png"), plot, width = width, height = height, dpi = 300),
           error = function(e) log_msg("png ggsave failed: ", e$message))
}

pca <- tryCatch(stats::prcomp(t(heat_mat), scale. = TRUE), error = function(e) NULL)
if (!is.null(pca)) {
  df <- data.frame(
    pca$x[, 1:2, drop = FALSE],
    tissue = sample_info$tissue[match(rownames(pca$x), sample_info$sample)],
    pair = factor(sample_info$pair[match(rownames(pca$x), sample_info$sample)]),
    sample = rownames(pca$x)
  )
  varp <- summary(pca)$importance[2, 1:2] * 100
  p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, color = tissue, shape = pair, label = sample)) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = "PCA (DESeq2 log matrix)",
                  x = sprintf("PC1 (%.1f%%)", varp[1]), y = sprintf("PC2 (%.1f%%)", varp[2]))
  save_gg(p, file.path(qc_dir, "PCA_all_samples"))
}

# -----------------------------------------------------------------------------
# 3. 作图 / 富集
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
    df$y <- if ("AveExpr" %in% names(df)) df$AveExpr else abs(df$log2FC)
    ylab <- "Average expression"
    hline <- NULL
  }
  df$set <- ifelse(df$gene %in% highlight, "selected", "other")
  df$label <- ifelse(df$gene %in% utils::head(highlight, 15), df$gene, NA)
  lfc_line <- log2(fc_line)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = y, color = set)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.4) +
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
  if (length(genes) > 200) genes <- genes[seq_len(200)]
  if (length(genes) < 2) {
    note_empty(outfile, "heatmap skipped (<2 genes)")
    return(invisible(NULL))
  }
  sub <- heat_mat[genes, , drop = FALSE]
  ann <- data.frame(
    Tissue = sample_info$tissue[match(colnames(sub), sample_info$sample)],
    Pair = factor(sample_info$pair[match(colnames(sub), sample_info$sample)]),
    row.names = colnames(sub)
  )
  if (!is.null(extra_ann)) {
    extra_ann <- extra_ann[match(colnames(sub), rownames(extra_ann)), , drop = FALSE]
    ann <- cbind(ann, extra_ann)
  }
  pal <- list(
    Tissue = c(Tumor = "#4C78A8", MFP_Met = "#F58518", TVI_Met = "#54A24B"),
    Pair = c(`1` = "#E45756", `2` = "#72B7B2", `3` = "grey70")
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
  if (length(entrez) < 3) {
    note_empty(file.path(go_dir, paste0(pref, "ORA_GO")), "too few mapped genes")
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
                    title_maybe_relaxed(ego, paste(label, "| ORA GO", ont)))
  }
  ek <- enrich_or_relax(
    function() clusterProfiler::enrichKEGG(gene = entrez, organism = "mmu", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
    function() clusterProfiler::enrichKEGG(gene = entrez, organism = "mmu", pvalueCutoff = 1, qvalueCutoff = 1),
    "enrichKEGG"
  )
  if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) {
    ek <- tryCatch(clusterProfiler::setReadable(ek, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"), error = function(e) ek)
  }
  plot_ora_object(ek, file.path(kg_dir, paste0(pref, "ORA_KEGG")),
                  title_maybe_relaxed(ek, paste(label, "| ORA KEGG")))
  if (has_pkg("ReactomePA")) {
    er <- enrich_or_relax(
      function() ReactomePA::enrichPathway(gene = entrez, organism = "mouse", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
      function() ReactomePA::enrichPathway(gene = entrez, organism = "mouse", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE),
      "enrichPathway"
    )
    plot_ora_object(er, file.path(pw_dir, paste0(pref, "ORA_Reactome_pathway")),
                    title_maybe_relaxed(er, paste(label, "| ORA Reactome")))
  }
  writeLines("This GO/Pathway/KEGG folder is ORA, NOT GSEA.",
             file.path(outdir, paste0(pref, "00_ORA_is_not_GSEA.txt")))
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
    enrich_or_relax(function() fun(pvalueCutoff = 0.05, minGSSize = 10),
                    function() fun(pvalueCutoff = 1, minGSSize = 5), lab)
  }
  go_bp <- gsea_one(function(pvalueCutoff, minGSSize) {
    clusterProfiler::gseGO(geneList = stats, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "ENTREZID",
                           minGSSize = minGSSize, maxGSSize = 500, pvalueCutoff = pvalueCutoff,
                           verbose = FALSE, eps = 0)
  }, "gseGO_BP")
  plot_gsea_object(go_bp, file.path(gsea_dir, "GSEA_GO_BP"), paste(label, "| GSEA GO BP"))
  kegg <- gsea_one(function(pvalueCutoff, minGSSize) {
    clusterProfiler::gseKEGG(geneList = stats, organism = "mmu", minGSSize = minGSSize, maxGSSize = 500,
                             pvalueCutoff = pvalueCutoff, verbose = FALSE, eps = 0)
  }, "gseKEGG")
  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
    kegg <- tryCatch(clusterProfiler::setReadable(kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"), error = function(e) kegg)
  }
  plot_gsea_object(kegg, file.path(gsea_dir, "GSEA_KEGG"), paste(label, "| GSEA KEGG"))
}

emit_subset <- function(comp_name, sub, tag, title, outdir, full_de, heat_mat, sample_info,
                        lfc_col = "log2FC", p_col = "pvalue", fc_line = 1, do_ora = FALSE) {
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
  tryCatch(plot_heatmap(heat_mat, sample_info, sub$gene, title, file.path(outdir, paste0(tag, "_heatmap"))),
           error = function(e) {
             while (grDevices::dev.cur() > 1) grDevices::dev.off()
             log_msg("heatmap failed: ", e$message)
           })
  if (isTRUE(do_ora) && !isTRUE(skip_enrich)) {
    tryCatch(run_ora_plots(sub$gene, sub, outdir, title, tag, lfc_col = lfc_col),
             error = function(e) log_msg("ORA failed: ", e$message))
  }
}

analyze_down_list <- function(comp_name, de, heat_mat, sample_info, lfc_col, p_col,
                              require_p = TRUE, do_gsea = TRUE, note = NULL) {
  base <- file.path(result_dir, comp_name)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  write_table(de, file.path(base, "DE_full"))
  lines <- c(
    "先看本目录 DE_full，再看 FoldChange/FC_1 与 FoldChange/FC_1.25。",
    "筛选：下调倍数 down_FC = 原位/肺 >= 1 与 >= 1.25。",
    "不做 Top50–300，没有 TopRank/。",
    "1 号原位只对 1 号肺，2 号只对 2 号肺；MFP 与 TVI 不合并。"
  )
  if (isTRUE(require_p)) {
    lines <- c(lines, "本目录还要求配对 limma p < 0.05。")
  } else {
    lines <- c(lines, "本目录是 1-vs-1，无法估计 p，只按 FC，不伪造 p 值。")
  }
  if (!is.null(note)) lines <- c(lines, note)
  writeLines(lines, file.path(base, "00_READ_ME.txt"))
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    keep <- !is.na(de[[lfc_col]]) & de[[lfc_col]] < 0 & (2^(-de[[lfc_col]]) >= fc)
    if (isTRUE(require_p)) {
      keep <- keep & !is.na(de[[p_col]]) & de[[p_col]] < p_cutoff
    }
    sub <- de[keep, , drop = FALSE]
    sub <- sub[order(sub[[lfc_col]], decreasing = FALSE), , drop = FALSE]
    emit_subset(comp_name, sub, nm, paste0(comp_name, " | down FC >= ", fc),
                file.path(base, "FoldChange", nm), de, heat_mat, sample_info,
                lfc_col, p_col, fc_line = fc,
                do_ora = identical(nm, "FC_1.25") && nrow(sub) <= 800)
  }
  if (isTRUE(do_gsea) && !isTRUE(skip_enrich)) {
    tryCatch(run_gsea_full(de, file.path(base, "00_GSEA_all_genes"), comp_name, lfc_col),
             error = function(e) log_msg("GSEA failed: ", e$message))
  }
}

is_down_fc <- function(lfc, fc) {
  !is.na(lfc) & lfc < 0 & (2^(-lfc) >= fc)
}
is_down <- function(lfc, p, fc) {
  is_down_fc(lfc, fc) & !is.na(p) & p < p_cutoff
}

# -----------------------------------------------------------------------------
# 4. 一一对应：Tumor-i vs Met-i
# -----------------------------------------------------------------------------
find_sample <- function(tissue, pair) {
  hit <- sample_info$sample[sample_info$tissue == tissue & sample_info$pair == pair]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

pair_lfc <- function(met_tissue, pair_id) {
  tm <- find_sample("Tumor", pair_id)
  met <- find_sample(met_tissue, pair_id)
  if (is.na(tm) || is.na(met)) return(NULL)
  lfc <- heat_mat[, met] - heat_mat[, tm]
  data.frame(
    gene = rownames(heat_mat),
    log2FC = as.numeric(lfc),
    AveExpr = rowMeans(heat_mat[, c(tm, met), drop = FALSE]),
    down_FC = 2^(-as.numeric(lfc)),
    primary = tm,
    lung = met,
    pair = pair_id,
    pvalue = NA_real_,
    stringsAsFactors = FALSE
  )
}

anno_de <- function(de) {
  if (is.null(de) || !is.data.frame(de) || !("gene" %in% names(de))) return(NULL)
  de$ensembl <- id_map$ensembl[match(de$gene, id_map$gene)]
  de$entrez <- id_map$entrez[match(de$gene, id_map$gene)]
  de$TM_mean <- tm_mean[de$gene]
  de$TM_expr_percentile <- tm_pct[de$gene]
  de$TM_low_quartile <- !is.na(de$TM_expr_percentile) & de$TM_expr_percentile <= 25
  de
}

paired_limma <- function(met_tissue, label) {
  pairs <- intersect(sample_info$pair[sample_info$tissue == "Tumor"],
                     sample_info$pair[sample_info$tissue == met_tissue])
  pairs <- sort(pairs[is.finite(pairs)])
  if (length(pairs) < 2) {
    log_msg("Need >=2 complete pairs for limma p: ", label)
    return(NULL)
  }
  tm <- vapply(pairs, function(p) find_sample("Tumor", p), character(1))
  mt <- vapply(pairs, function(p) find_sample(met_tissue, p), character(1))
  use <- c(tm, mt)
  y <- heat_mat[, use, drop = FALSE]
  pair <- factor(c(pairs, pairs))
  tissue <- factor(c(rep("Tumor", length(pairs)), rep("Met", length(pairs))), levels = c("Tumor", "Met"))
  design <- stats::model.matrix(~ pair + tissue)
  fit <- tryCatch(
    limma::eBayes(limma::lmFit(y, design), trend = TRUE, robust = TRUE),
    error = function(e) {
      log_msg("Paired limma failed (", label, "): ", e$message)
      NULL
    }
  )
  if (is.null(fit)) return(NULL)
  coef <- grep("tissueMet", colnames(design), value = TRUE)
  tt <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
  data.frame(
    gene = rownames(tt),
    log2FC = tt$logFC,
    AveExpr = tt$AveExpr,
    t = tt$t,
    pvalue = tt$P.Value,
    padj = tt$adj.P.Val,
    down_FC = 2^(-tt$logFC),
    n_pairs = length(pairs),
    stringsAsFactors = FALSE
  )
}

# 主分析：MFP 脂肪垫路径（原位 → 肺）
log_msg("Q1 lung: paired MFP-Met vs Tumor by replicate number")
de_mfp_p1 <- anno_de(pair_lfc("MFP_Met", 1))
de_mfp_p2 <- anno_de(pair_lfc("MFP_Met", 2))
de_mfp_p3 <- anno_de(pair_lfc("MFP_Met", 3))
mfp_p3_true_pair <- !is.null(de_mfp_p3)
if (is.null(de_mfp_p3)) {
  log_msg("No Tumor-3; MFP-Met-3 FC vs mean of Tumor-1/2 (not a true pair)")
  met3 <- find_sample("MFP_Met", 3)
  if (!is.na(met3)) {
    lfc <- heat_mat[, met3] - tm_mean
    de_mfp_p3 <- anno_de(data.frame(
      gene = rownames(heat_mat), log2FC = as.numeric(lfc),
      AveExpr = rowMeans(cbind(tm_mean, heat_mat[, met3])),
      down_FC = 2^(-as.numeric(lfc)), primary = "Tumor_mean", lung = met3,
      pair = 3, pvalue = NA_real_, stringsAsFactors = FALSE
    ))
  }
}
de_mfp_paired <- anno_de(paired_limma("MFP_Met", "MFP paired Tumor vs Met"))

if (!is.null(de_mfp_p1)) {
  analyze_down_list("01_lung_down_pair1_Tumor1_vs_MFPMet1", de_mfp_p1, heat_mat, sample_info,
                    "log2FC", "pvalue", require_p = FALSE, do_gsea = FALSE)
}
if (!is.null(de_mfp_p2)) {
  analyze_down_list("01_lung_down_pair2_Tumor2_vs_MFPMet2", de_mfp_p2, heat_mat, sample_info,
                    "log2FC", "pvalue", require_p = FALSE, do_gsea = FALSE)
}
if (!is.null(de_mfp_p3)) {
  if (isTRUE(mfp_p3_true_pair)) {
    analyze_down_list("01_lung_down_pair3_Tumor3_vs_MFPMet3", de_mfp_p3, heat_mat, sample_info,
                      "log2FC", "pvalue", require_p = FALSE, do_gsea = FALSE)
  } else {
    analyze_down_list(
      "01_lung_down_MFPMet3_vs_Tumor_mean", de_mfp_p3, heat_mat, sample_info,
      "log2FC", "pvalue", require_p = FALSE, do_gsea = FALSE,
      note = "没有 Tumor-3。MFP-Met-3 相对 Tumor-1/2 均值，不是 1 对 1。"
    )
  }
}

# 两对都下调，并且配对 limma p < 0.05
merged <- data.frame(
  gene = rownames(heat_mat),
  ensembl = id_map$ensembl,
  entrez = id_map$entrez,
  TM_mean = tm_mean,
  TM_expr_percentile = tm_pct,
  stringsAsFactors = FALSE
)
if (!is.null(de_mfp_p1)) {
  merged$log2FC_pair1 <- de_mfp_p1$log2FC[match(merged$gene, de_mfp_p1$gene)]
  merged$down_FC_pair1 <- de_mfp_p1$down_FC[match(merged$gene, de_mfp_p1$gene)]
}
if (!is.null(de_mfp_p2)) {
  merged$log2FC_pair2 <- de_mfp_p2$log2FC[match(merged$gene, de_mfp_p2$gene)]
  merged$down_FC_pair2 <- de_mfp_p2$down_FC[match(merged$gene, de_mfp_p2$gene)]
}
if (!is.null(de_mfp_paired)) {
  merged$log2FC <- de_mfp_paired$log2FC[match(merged$gene, de_mfp_paired$gene)]
  merged$pvalue <- de_mfp_paired$pvalue[match(merged$gene, de_mfp_paired$gene)]
  merged$padj <- de_mfp_paired$padj[match(merged$gene, de_mfp_paired$gene)]
  merged$down_FC <- de_mfp_paired$down_FC[match(merged$gene, de_mfp_paired$gene)]
} else {
  merged$log2FC <- rowMeans(cbind(merged$log2FC_pair1, merged$log2FC_pair2), na.rm = TRUE)
  merged$pvalue <- NA_real_
  merged$down_FC <- 2^(-merged$log2FC)
}
merged$mean_pair_log2FC <- rowMeans(cbind(merged$log2FC_pair1, merged$log2FC_pair2), na.rm = TRUE)
write_table(merged, file.path(qc_dir, "paired_MFP_Met_vs_Tumor_full"))

both_de <- merged
both_de$log2FC <- both_de$mean_pair_log2FC
keep_both <- is_down_fc(merged$log2FC_pair1, 1) & is_down_fc(merged$log2FC_pair2, 1)
if ("pvalue" %in% names(merged) && any(!is.na(merged$pvalue))) {
  both_de$pvalue[!keep_both] <- 1
  analyze_down_list(
    "01_lung_down_both_pairs_MFP", both_de, heat_mat, sample_info,
    "mean_pair_log2FC", "pvalue", require_p = TRUE, do_gsea = TRUE,
    note = "两对都下调，且配对 limma p<0.05。MFP-Met-3 不进入这一步。"
  )
} else {
  both_de$pvalue <- ifelse(keep_both, 0, 1)
  analyze_down_list(
    "01_lung_down_both_pairs_MFP", both_de, heat_mat, sample_info,
    "mean_pair_log2FC", "pvalue", require_p = FALSE, do_gsea = TRUE,
    note = "配对 limma 失败时只要求两对 FC 都下调。"
  )
}

base_both <- file.path(result_dir, "01_lung_down_both_pairs_MFP")
for (nm in names(fc_cutoffs)) {
  fc <- unname(fc_cutoffs[[nm]])
  sub <- merged[is_down_fc(merged$log2FC_pair1, fc) & is_down_fc(merged$log2FC_pair2, fc), ]
  if ("pvalue" %in% names(merged) && any(!is.na(merged$pvalue))) {
    sub <- sub[!is.na(sub$pvalue) & sub$pvalue < p_cutoff, ]
  }
  sub <- sub[order(sub$mean_pair_log2FC), ]
  sub$log2FC <- sub$mean_pair_log2FC
  emit_subset(
    "01_lung_down_both_pairs_MFP", sub, paste0(nm, "_BOTH_pairs"),
    paste0("Q1 both pairs MFP down | FC >= ", fc),
    file.path(base_both, "FoldChange", paste0(nm, "_BOTH_pairs")),
    transform(merged, log2FC = mean_pair_log2FC), heat_mat, sample_info,
    "mean_pair_log2FC", "pvalue", fc,
    do_ora = identical(nm, "FC_1.25") && nrow(sub) <= 800
  )
}

# TVI 路径（实验性肺定植，单独出表，不和 MFP 合并）
log_msg("TVI path: paired TVI-Met vs Tumor (not merged with MFP)")
de_tvi_p1 <- anno_de(pair_lfc("TVI_Met", 1))
de_tvi_p2 <- anno_de(pair_lfc("TVI_Met", 2))
de_tvi_paired <- anno_de(paired_limma("TVI_Met", "TVI paired"))
merged_tvi <- merged[, c("gene", "ensembl", "entrez", "TM_mean"), drop = FALSE]
if (!is.null(de_tvi_p1)) merged_tvi$log2FC_pair1 <- de_tvi_p1$log2FC[match(merged_tvi$gene, de_tvi_p1$gene)]
if (!is.null(de_tvi_p2)) merged_tvi$log2FC_pair2 <- de_tvi_p2$log2FC[match(merged_tvi$gene, de_tvi_p2$gene)]
if (!is.null(de_tvi_paired)) {
  merged_tvi$pvalue <- de_tvi_paired$pvalue[match(merged_tvi$gene, de_tvi_paired$gene)]
  merged_tvi$log2FC <- de_tvi_paired$log2FC[match(merged_tvi$gene, de_tvi_paired$gene)]
} else {
  merged_tvi$log2FC <- rowMeans(cbind(merged_tvi$log2FC_pair1, merged_tvi$log2FC_pair2), na.rm = TRUE)
  merged_tvi$pvalue <- NA_real_
}
merged_tvi$mean_pair_log2FC <- rowMeans(cbind(merged_tvi$log2FC_pair1, merged_tvi$log2FC_pair2), na.rm = TRUE)
write_table(merged_tvi, file.path(qc_dir, "paired_TVI_Met_vs_Tumor_full"))
if (!is.null(de_tvi_p1)) {
  analyze_down_list("01b_TVI_lung_down_pair1", de_tvi_p1, heat_mat, sample_info,
                    "log2FC", "pvalue", require_p = FALSE, do_gsea = FALSE)
}
if (!is.null(de_tvi_p2)) {
  analyze_down_list("01b_TVI_lung_down_pair2", de_tvi_p2, heat_mat, sample_info,
                    "log2FC", "pvalue", require_p = FALSE, do_gsea = FALSE)
}
if (!is.null(de_tvi_p1) && !is.null(de_tvi_p2)) {
  tvi_both <- merged_tvi
  tvi_both$pvalue[!is_down_fc(merged_tvi$log2FC_pair1, 1) | !is_down_fc(merged_tvi$log2FC_pair2, 1)] <- 1
  analyze_down_list(
    "01b_TVI_lung_down_both_pairs", tvi_both, heat_mat, sample_info,
    "mean_pair_log2FC", "pvalue",
    require_p = any(!is.na(merged_tvi$pvalue)), do_gsea = FALSE,
    note = "TVI 两对都下调；与 MFP 分开，不合并。"
  )
}

bone_dir <- file.path(result_dir, "03_bone_not_in_GSE146012")
dir.create(bone_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  c("GSE146012 没有骨 / 骨髓样品。",
    "8 个样品只有：Tumor-1/2、MFP-Met-1/2/3、TVI-Met-1/2/3。",
    "不能从本数据寻找「原位低表达促进骨转移」，也不能做肺 vs 骨器官特异。",
    "骨转移请用 GSE165393 或 GSE37975。",
    "不要把 MFP 肺和 TVI 肺并成一组冒充多器官。"),
  file.path(bone_dir, "00_NOTE.txt")
)
log_msg("Q1 bone: skipped (no bone samples in GSE146012)")

# -----------------------------------------------------------------------------
# 5. Q2：MFP 肺特异（相对 TVI），不是骨特异
# -----------------------------------------------------------------------------
log_msg("Q2: MFP-lung-specific down (not TVI); bone not available")
q2 <- merged
q2$log2FC <- q2$mean_pair_log2FC
mark_mfp <- is_down_fc(merged$log2FC_pair1, 1) & is_down_fc(merged$log2FC_pair2, 1)
mark_tvi <- is_down_fc(merged_tvi$log2FC_pair1, 1) & is_down_fc(merged_tvi$log2FC_pair2, 1)
mark_spec <- mark_mfp & !mark_tvi
q2$pvalue <- ifelse(mark_spec, q2$pvalue, 1)
if (all(is.na(q2$pvalue[mark_spec]))) q2$pvalue[mark_spec] <- 0
analyze_down_list(
  "02_MFP_lung_specific_vs_TVI", q2, heat_mat, sample_info,
  "mean_pair_log2FC", "pvalue",
  require_p = any(!is.na(merged$pvalue)),
  note = "MFP 两对都下调，且 TVI 两对未同时下调。这是接种路径特异，不是骨器官特异。"
)

if (has_pkg("ggvenn")) {
  vdf <- list(
    MFP_pair_down = merged$gene[mark_mfp],
    TVI_pair_down = merged_tvi$gene[mark_tvi]
  )
  p <- ggvenn::ggvenn(vdf, fill_color = c("#F58518", "#54A24B")) +
    ggplot2::labs(title = "Paired down genes: MFP lung vs TVI lung")
  save_gg(p, file.path(qc_dir, "venn_MFP_vs_TVI_lung_down"))
}

# -----------------------------------------------------------------------------
# 6. Q3 三套神经分数（原位标记）+ 负相关
# -----------------------------------------------------------------------------
log_msg("Q3: three neural signatures")
pni_dir <- file.path(result_dir, "06_neural_invasion")
dir.create(pni_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  c("没有病理 PNI。三套分数：施旺 / 神经营养 / 轴突导向。",
    "原位只有 Tumor-1/2，不能在 2 个点上做 Spearman p。",
    "负相关在 MFP 路径 5 个样品上计算（2 原位 + 3 肺），p<0.05 且 rho<0。",
    "同时给出两对原位各自的分数。"),
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
  if (is.null(msig_bp) || !"gene_symbol" %in% names(msig_bp)) return(character())
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
  Axon_guidance = "AXON_GUIDANCE|NEURON_PROJECTION_GUIDANCE|SEMAPHORIN|NETRIN|SLIT_ROBO"
)
sig_genes <- lapply(names(curated), function(nm) {
  unique(intersect(c(curated[[nm]], msig_symbols(go_pat[[nm]])), rownames(heat_mat)))
})
names(sig_genes) <- names(curated)
for (nm in names(sig_genes)) {
  log_msg("Signature ", nm, ": ", length(sig_genes[[nm]]), " genes")
  utils::write.csv(data.frame(gene = sig_genes[[nm]]),
                   file.path(pni_dir, paste0("signature_genes_", nm, ".csv")), row.names = FALSE)
}

mean_z_score <- function(genes) {
  genes <- intersect(genes, rownames(heat_mat))
  if (length(genes) < 2) return(setNames(rep(NA_real_, ncol(heat_mat)), colnames(heat_mat)))
  z <- t(scale(t(heat_mat[genes, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}
score_mat <- sapply(sig_genes, mean_z_score)
if (has_pkg("GSVA")) {
  gs_ok <- sig_genes[vapply(sig_genes, length, integer(1)) >= 8]
  if (length(gs_ok) > 0) {
    gsva_mat <- tryCatch({
      if (utils::packageVersion("GSVA") >= "1.50.0") {
        GSVA::gsva(GSVA::ssgseaParam(as.matrix(heat_mat), gs_ok), verbose = FALSE)
      } else {
        GSVA::gsva(as.matrix(heat_mat), gs_ok, method = "ssgsea", verbose = FALSE)
      }
    }, error = function(e) {
      log_msg("GSVA failed, keep mean z-score: ", e$message)
      NULL
    })
    if (!is.null(gsva_mat)) {
      for (nm in rownames(gsva_mat)) score_mat[, nm] <- as.numeric(gsva_mat[nm, colnames(heat_mat)])
    }
  }
}
score_mat <- as.data.frame(score_mat)
score_mat$sample <- colnames(heat_mat)
score_tab <- merge(sample_info, score_mat, by = "sample", all.x = TRUE)
write_table(score_tab, file.path(pni_dir, "sample_three_neural_scores"))
tm_mark <- score_tab[score_tab$tissue == "Tumor", ]
write_table(tm_mark, file.path(pni_dir, "Tumor_neural_invasion_marker"))

score_long <- tidyr::pivot_longer(
  score_tab, cols = dplyr::all_of(names(curated)),
  names_to = "signature", values_to = "score"
)
p_bar <- ggplot2::ggplot(score_long, ggplot2::aes(x = sample, y = score, fill = tissue)) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(~ signature, ncol = 1, scales = "free_y") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(title = "Three neural invasion signature scores", x = NULL, y = "score")
save_gg(p_bar, file.path(pni_dir, "three_scores_by_sample"), 9, 10)

mfp_cols <- sample_info$sample[sample_info$tissue %in% c("Tumor", "MFP_Met")]
spearman_fast <- function(score_vec, cols) {
  sc <- as.numeric(score_vec[cols])
  ok <- is.finite(sc)
  cols <- cols[ok]
  sc <- sc[ok]
  n <- length(sc)
  if (n < 4) {
    return(data.frame(gene = rownames(heat_mat), rho = NA_real_, pvalue = NA_real_))
  }
  ranked_sc <- rank(sc)
  rho <- apply(heat_mat[, cols, drop = FALSE], 1, function(x) {
    stats::cor(rank(x), ranked_sc, method = "pearson", use = "complete.obs")
  })
  tstat <- rho * sqrt((n - 2) / pmax(1 - rho^2, 1e-8))
  pval <- 2 * stats::pt(-abs(tstat), df = n - 2)
  data.frame(gene = names(rho), rho = as.numeric(rho), pvalue = as.numeric(pval), stringsAsFactors = FALSE)
}

neg_both <- list()
lung_down_genes <- merged$gene[is_down_fc(merged$log2FC_pair1, 1) & is_down_fc(merged$log2FC_pair2, 1)]

for (nm in names(curated)) {
  log_msg("Neural axis: ", nm)
  ax_dir <- file.path(pni_dir, nm)
  dir.create(ax_dir, recursive = TRUE, showWarnings = FALSE)
  sc <- setNames(score_tab[[nm]], score_tab$sample)
  plot_heatmap(heat_mat, sample_info, sig_genes[[nm]], paste(nm, "signature genes"),
               file.path(ax_dir, paste0("heatmap_", nm, "_signature")))
  cor_df <- spearman_fast(sc, mfp_cols)
  cor_df$padj <- stats::p.adjust(cor_df$pvalue, method = "BH")
  cor_df$ensembl <- id_map$ensembl[match(cor_df$gene, id_map$gene)]
  cor_df$log2FC_pair1 <- merged$log2FC_pair1[match(cor_df$gene, merged$gene)]
  cor_df$log2FC_pair2 <- merged$log2FC_pair2[match(cor_df$gene, merged$gene)]
  cor_df$down_in_MFP_both_pairs <- cor_df$gene %in% lung_down_genes
  write_table(cor_df, file.path(ax_dir, paste0(nm, "_gene_vs_score_spearman_MFP_path")))
  neg <- cor_df[is.finite(cor_df$rho) & cor_df$rho < 0 & cor_df$pvalue < p_cutoff, ]
  neg <- neg[order(neg$rho), ]
  neg$log2FC <- neg$rho
  neg_both[[nm]] <- neg
  write_table(neg, file.path(ax_dir, paste0(nm, "_neg_cor_p0.05")))
  emit_subset(
    paste0("06_neural_invasion/", nm), neg, paste0(nm, "_neg_cor_p0.05"),
    paste0(nm, " | negative correlation on MFP-path samples"),
    file.path(ax_dir, "neg_cor_p0.05"),
    transform(cor_df, log2FC = rho), heat_mat, sample_info,
    lfc_col = "rho", p_col = "pvalue",
    do_ora = nrow(neg) <= 800 && !isTRUE(skip_enrich)
  )
  write_table(neg[neg$gene %in% lung_down_genes, ],
              file.path(ax_dir, paste0(nm, "_neg_cor_AND_paired_MFP_lung_down")))
  if (!isTRUE(skip_enrich)) {
    tryCatch(
      run_gsea_full(transform(cor_df, log2FC = rho), file.path(ax_dir, "00_GSEA_all_genes"),
                    paste(nm, "correlation"), stat_col = "rho"),
      error = function(e) log_msg(nm, " GSEA failed: ", e$message)
    )
  }
}

if (has_pkg("ggvenn")) {
  v3 <- lapply(neg_both, function(x) x$gene)
  p <- ggvenn::ggvenn(v3, fill_color = c("#4C78A8", "#F58518", "#54A24B")) +
    ggplot2::labs(title = "Negatively correlated genes (p<0.05) across three neural scores")
  save_gg(p, file.path(pni_dir, "venn_three_neg_cor"), 8, 6)
}

# -----------------------------------------------------------------------------
# 7. Q4：三种神经分数 vs 配对肺
# -----------------------------------------------------------------------------
log_msg("Q4: neural scores vs matched MFP lung")
q4_dir <- file.path(result_dir, "07_neural_vs_lung")
dir.create(q4_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  c("每种神经分数单独对配对 MFP 肺：",
    "  - 1 号：Tumor-1 vs MFP-Met-1 分数差",
    "  - 2 号：Tumor-2 vs MFP-Met-2 分数差",
    "  - 原位组 vs MFP 肺组的 Wilcoxon（n 很小）",
    "  - 负相关基因 ∩ 两对都下调的肺基因"),
  file.path(q4_dir, "00_READ_ME.txt")
)

assoc_rows <- list()
ov_rows <- list()
for (nm in names(curated)) {
  ax_dir <- file.path(q4_dir, nm)
  dir.create(ax_dir, recursive = TRUE, showWarnings = FALSE)
  sc <- setNames(score_tab[[nm]], score_tab$sample)
  for (pr in c(1, 2)) {
    tm <- find_sample("Tumor", pr)
    mt <- find_sample("MFP_Met", pr)
    if (is.na(tm) || is.na(mt)) next
    assoc_rows[[paste(nm, pr)]] <- data.frame(
      signature = nm, pair = pr, primary = tm, lung = mt,
      score_primary = sc[[tm]], score_lung = sc[[mt]],
      delta = sc[[mt]] - sc[[tm]],
      stringsAsFactors = FALSE
    )
  }
  tm_sc <- sc[sample_info$sample[sample_info$tissue == "Tumor"]]
  lu_sc <- sc[sample_info$sample[sample_info$tissue == "MFP_Met"]]
  tm_sc <- tm_sc[is.finite(tm_sc)]
  lu_sc <- lu_sc[is.finite(lu_sc)]
  p_w <- if (length(tm_sc) >= 2 && length(lu_sc) >= 2) {
    suppressWarnings(stats::wilcox.test(lu_sc, tm_sc)$p.value)
  } else NA_real_
  write_table(
    data.frame(signature = nm, n_tumor = length(tm_sc), n_lung = length(lu_sc),
               mean_tumor = mean(tm_sc), mean_lung = mean(lu_sc),
               delta = mean(lu_sc) - mean(tm_sc), wilcox_p = p_w),
    file.path(ax_dir, paste0(nm, "_score_Tumor_vs_MFP_Met"))
  )
  ov <- neg_both[[nm]][neg_both[[nm]]$gene %in% lung_down_genes, ]
  write_table(ov, file.path(ax_dir, paste0(nm, "_neg_cor_AND_both_pairs_MFP_down")))
  ov_rows[[nm]] <- data.frame(
    signature = nm, n_neg_cor = nrow(neg_both[[nm]]),
    n_lung_down = length(lung_down_genes), n_overlap = nrow(ov)
  )
}
assoc_all <- do.call(rbind, assoc_rows)
write_table(assoc_all, file.path(q4_dir, "paired_neural_score_delta_MFP"))
write_table(dplyr::bind_rows(ov_rows), file.path(q4_dir, "neg_cor_overlap_counts"))

p_assoc <- ggplot2::ggplot(assoc_all, ggplot2::aes(x = signature, y = delta, fill = factor(pair))) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.65) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(title = "Neural score: matched MFP lung minus primary",
                y = "score(MFP-Met-i) - score(Tumor-i)", x = NULL, fill = "pair")
save_gg(p_assoc, file.path(q4_dir, "delta_neural_score_paired_MFP"), 8, 5)

# -----------------------------------------------------------------------------
# 8. 总表
# -----------------------------------------------------------------------------
sum_dir <- file.path(result_dir, "08_summary")
dir.create(sum_dir, recursive = TRUE, showWarnings = FALSE)

q1_list <- list()
for (nm in names(fc_cutoffs)) {
  fc <- unname(fc_cutoffs[[nm]])
  sub <- merged[is_down_fc(merged$log2FC_pair1, fc) & is_down_fc(merged$log2FC_pair2, fc), ]
  if ("pvalue" %in% names(merged) && any(!is.na(merged$pvalue))) {
    sub <- sub[!is.na(sub$pvalue) & sub$pvalue < p_cutoff, ]
  }
  sub <- sub[order(sub$mean_pair_log2FC), ]
  q1_list[[nm]] <- sub
  write_table(sub, file.path(sum_dir, paste0("Q1_both_pairs_MFP_lung_down_p0.05_", nm)))
}
if (!is.null(de_mfp_p1)) {
  write_table(de_mfp_p1[is_down_fc(de_mfp_p1$log2FC, 1.25), ],
              file.path(sum_dir, "Q1_pair1_Tumor1_vs_MFPMet1_FC1.25"))
}
if (!is.null(de_mfp_p2)) {
  write_table(de_mfp_p2[is_down_fc(de_mfp_p2$log2FC, 1.25), ],
              file.path(sum_dir, "Q1_pair2_Tumor2_vs_MFPMet2_FC1.25"))
}
q2s <- merged[mark_spec, , drop = FALSE]
if ("pvalue" %in% names(merged) && any(!is.na(merged$pvalue))) {
  q2s <- q2s[!is.na(q2s$pvalue) & q2s$pvalue < p_cutoff, , drop = FALSE]
}
q2s <- q2s[order(q2s$mean_pair_log2FC), , drop = FALSE]
write_table(q2s, file.path(sum_dir, "Q2_MFP_lung_specific_vs_TVI"))
writeLines("GSE146012 has no bone. Q2 bone-specific list does not exist.",
           file.path(sum_dir, "Q2_bone_specific_NOT_AVAILABLE.txt"))
write_table(tm_mark, file.path(sum_dir, "Q3_Tumor_three_neural_scores"))
write_table(assoc_all, file.path(sum_dir, "Q4_neural_scores_vs_matched_MFP_lung"))
for (nm in names(neg_both)) {
  write_table(neg_both[[nm]], file.path(sum_dir, paste0("Q3_", nm, "_neg_cor_p0.05")))
}

n_tab <- data.frame(
  question = c(
    paste0("Q1 both-pairs MFP lung down ", names(fc_cutoffs)),
    "Q2 MFP-specific vs TVI",
    "Q2 bone-specific",
    paste0("Q3 ", names(neg_both), " neg-cor p<0.05")
  ),
  n_genes = c(
    vapply(q1_list, nrow, integer(1)),
    nrow(q2s),
    0,
    vapply(neg_both, nrow, integer(1))
  ),
  stringsAsFactors = FALSE
)
write_table(n_tab, file.path(sum_dir, "gene_counts"))
log_msg("Done. Summary counts:")
print(n_tab)
log_msg("Open: ", sum_dir)
log_msg("Bone metastasis cannot be analysed in GSE146012 (no bone samples).")
