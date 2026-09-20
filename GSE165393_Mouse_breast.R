#!/usr/bin/env Rscript
# =============================================================================
# GSE165393_Mouse_breast.R
# MMTV-PyMT 器官来源细胞系（GSE165393）
#
# 用法（Windows）：
#   setwd("E:/R/Mouse Breast")
#   source("GSE165393_Mouse_breast.R")
#
# 三个问题（「低表达促进转移」= 转移抑制候选：转移灶相对原位下调）：
#   1) 原位内低表达、同时与肺转移和骨（骨髓）转移相关的共享下调基因
#   2) 只与肺、或只与骨相关的器官特异下调基因
#   3) 按施旺细胞、神经营养因子、轴突导向分别打分，找负相关基因，
#      并分别检验这三种分数与肺转移、骨转移的关系
#
# 数据注意：
#   - 14 个 bulk 是 FACS 后培养的细胞系，不是新鲜肿瘤块
#   - BM 是骨髓来源，不是皮质骨转移灶
#   - 本套数据没有病理 PNI 标签，第 3 问用神经/PNI 基因集打分
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")
options(clusterProfiler.download.method = "auto")

# -----------------------------------------------------------------------------
# 0. 依赖
# -----------------------------------------------------------------------------
cran_required <- c(
  "dplyr", "tidyr", "tibble", "stringr", "ggplot2", "ggrepel", "pheatmap",
  "RColorBrewer", "matrixStats", "cowplot", "writexl"
)
cran_optional <- c("ggvenn", "GSVA")
bioc_required <- c(
  "limma", "edgeR", "clusterProfiler", "org.Mm.eg.db", "enrichplot",
  "DOSE", "AnnotationDbi", "fgsea", "msigdbr"
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
  if (length(still) > 0) {
    message("可选包未安装，相关分析将跳过: ", paste(still, collapse = ", "))
  }
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
# 1. 路径
# -----------------------------------------------------------------------------
resolve_project_dir <- function() {
  env_dir <- Sys.getenv("GSE165393_DIR", unset = "")
  candidates <- c(
    env_dir,
    "E:/R/Mouse Breast",
    "E:\\R\\Mouse Breast",
    file.path(getwd()),
    file.path(getwd(), "Mouse Breast"),
    file.path(getwd(), "data")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  is_data_dir <- function(d) {
    if (!dir.exists(d)) return(FALSE)
    hits <- list.files(
      d, pattern = "GSE165393.*AllTissues_TPM.*", recursive = TRUE,
      full.names = FALSE, ignore.case = TRUE
    )
    length(hits) > 0
  }
  for (d in candidates) {
    if (is_data_dir(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

find_tpm_file <- function(root) {
  hits <- list.files(
    root, pattern = "GSE165393.*AllTissues_TPM.*\\.(csv|txt)(\\.gz)?$",
    recursive = TRUE, full.names = TRUE, ignore.case = TRUE
  )
  if (length(hits) == 0) {
    hits <- list.files(
      root, pattern = "AllTissues_TPM", recursive = TRUE,
      full.names = TRUE, ignore.case = TRUE
    )
  }
  if (length(hits) == 0) {
    stop(
      "未找到 GSE165393_AllTissues_TPM.csv.gz。请把 GEO 补充文件放到: ",
      root
    )
  }
  hits[[1]]
}

project_dir <- resolve_project_dir()
result_dir  <- file.path(project_dir, "results_GSE165393")
log_dir     <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("GSE165393_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

p_cutoff   <- 0.05
fc_cutoffs <- c("FC_1" = 1, "FC_1.25" = 1.25)
top_ns     <- c(50, 75, 100, 150, 200, 250, 300)
skip_enrich <- nzchar(Sys.getenv("GSE165393_SKIP_ENRICHMENT"))

log_msg("Project dir: ", project_dir)
log_msg("Results dir: ", result_dir)
log_msg("Significance gate: p < ", p_cutoff, " ; down-FC strata: met/MFP FC < 1 and < 1/1.25")

# -----------------------------------------------------------------------------
# 2. 读入 TPM 与样品表
# -----------------------------------------------------------------------------
tpm_file <- find_tpm_file(project_dir)
log_msg("TPM file: ", tpm_file)
raw <- utils::read.csv(tpm_file, check.names = FALSE)
names(raw) <- gsub('^\\"|\\"$', "", names(raw))
if (!"gene_id" %in% names(raw)) {
  names(raw)[1] <- "gene_id"
}
gene_id <- as.character(raw$gene_id)
ensembl <- sub("\\.\\d+$", "", gene_id)
mat_tpm <- as.matrix(raw[, setdiff(names(raw), "gene_id"), drop = FALSE])
storage.mode(mat_tpm) <- "numeric"
colnames(mat_tpm) <- gsub("\\.TPM$", "", colnames(mat_tpm))
rownames(mat_tpm) <- gene_id

parse_sample <- function(nm) {
  tissue <- dplyr::case_when(
    grepl("^MFP", nm, ignore.case = TRUE) ~ "MFP",
    grepl("^LU|^Lun", nm, ignore.case = TRUE) ~ "LU",
    grepl("^LN", nm, ignore.case = TRUE) ~ "LN",
    grepl("^BM", nm, ignore.case = TRUE) ~ "BM",
    TRUE ~ "other"
  )
  cd44 <- dplyr::case_when(
    grepl("high", nm, ignore.case = TRUE) ~ "high",
    grepl("low", nm, ignore.case = TRUE) ~ "low",
    TRUE ~ "unsplit"
  )
  data.frame(sample = nm, tissue = tissue, cd44 = cd44, stringsAsFactors = FALSE)
}

sample_info <- do.call(rbind, lapply(colnames(mat_tpm), parse_sample))
rownames(sample_info) <- sample_info$sample
sample_info$group <- paste(sample_info$tissue, sample_info$cd44, sep = "_")
log_msg("Samples: ", paste(sample_info$sample, collapse = ", "))
print(sample_info)

# -----------------------------------------------------------------------------
# 3. 基因符号
# -----------------------------------------------------------------------------
map_ensembl <- function(ens) {
  uniq <- unique(ens)
  sym <- tryCatch(
    AnnotationDbi::mapIds(org.Mm.eg.db, keys = uniq, column = "SYMBOL",
                          keytype = "ENSEMBL", multiVals = "first"),
    error = function(e) {
      log_msg("SYMBOL map failed: ", e$message)
      setNames(rep(NA_character_, length(uniq)), uniq)
    }
  )
  ent <- tryCatch(
    AnnotationDbi::mapIds(org.Mm.eg.db, keys = uniq, column = "ENTREZID",
                          keytype = "ENSEMBL", multiVals = "first"),
    error = function(e) setNames(rep(NA_character_, length(uniq)), uniq)
  )
  data.frame(
    ensembl = uniq,
    gene = ifelse(!is.na(sym[uniq]) & nzchar(sym[uniq]), unname(sym[uniq]), uniq),
    entrez = unname(ent[uniq]),
    stringsAsFactors = FALSE
  )
}

id_map <- map_ensembl(ensembl)
id_map <- id_map[match(ensembl, id_map$ensembl), ]
id_map$gene_id <- gene_id
mean_tpm <- rowMeans(mat_tpm, na.rm = TRUE)
ord <- order(mean_tpm, decreasing = TRUE)
keep_row <- !duplicated(id_map$gene[ord])
keep_idx <- sort(ord[keep_row])
mat_tpm <- mat_tpm[keep_idx, , drop = FALSE]
id_map <- id_map[keep_idx, , drop = FALSE]
rownames(mat_tpm) <- id_map$gene
log_msg("Unique symbols/IDs after collapse: ", nrow(mat_tpm))

# -----------------------------------------------------------------------------
# 4. 过滤 + 分位数标准化（TPM，不用 DESeq2）
# -----------------------------------------------------------------------------
filter_low_expression <- function(mat) {
  min_n <- 2
  keep <- rowSums(mat > 1, na.rm = TRUE) >= min_n
  if (sum(keep) < 200) {
    keep <- rowSums(mat > 0, na.rm = TRUE) >= min_n
    log_msg("Strict TPM filter too small; fallback to >0 in >=2 samples")
  }
  log_msg("Low-expression filter: keep ", sum(keep), " / ", nrow(mat))
  mat[keep, , drop = FALSE]
}

mat_f <- filter_low_expression(mat_tpm)
id_map <- id_map[match(rownames(mat_f), id_map$gene), ]
log_mat <- limma::normalizeBetweenArrays(log2(pmax(mat_f, 0) + 1), method = "quantile")
heat_mat <- log_mat
log_msg("Normalization: quantile on log2(TPM+1)")

mfp_mean_tpm <- rowMeans(mat_f[, sample_info$sample[sample_info$tissue == "MFP"], drop = FALSE], na.rm = TRUE)
mfp_pct <- rank(mfp_mean_tpm, ties.method = "average") / length(mfp_mean_tpm) * 100

# -----------------------------------------------------------------------------
# 5. 工具函数
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
  if (length(leftover) > 0) {
    m <- tryCatch(
      clusterProfiler::bitr(leftover, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db),
      error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
    )
    if (nrow(m) > 0) extra <- data.frame(gene = m$SYMBOL, entrez = m$ENTREZID)
  }
  out <- rbind(from_map, extra)
  out <- out[!duplicated(out$gene), ]
  out
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
    ggplot2::geom_point(alpha = 0.7, size = 1.4) +
    ggplot2::scale_color_manual(values = c(other = "grey70", selected = "#1F77B4")) +
    ggplot2::geom_vline(xintercept = c(-lfc_line, lfc_line), linetype = 2, color = "grey40") +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, x = "log2 Fold Change (met - MFP)", y = ylab, color = NULL)
  if (!is.null(hline)) p <- p + ggplot2::geom_hline(yintercept = hline, linetype = 2, color = "grey40")
  save_gg(p, outfile)
}

plot_heatmap <- function(heat_mat, sample_info, genes, title, outfile, extra_ann = NULL) {
  genes <- intersect(genes, rownames(heat_mat))
  if (length(genes) > 200) {
    log_msg("Heatmap truncated to 200 genes: ", title)
    genes <- genes[seq_len(200)]
  }
  if (length(genes) < 2) {
    note_empty(outfile, "heatmap skipped (<2 genes)")
    return(invisible(NULL))
  }
  sub <- heat_mat[genes, , drop = FALSE]
  ann <- data.frame(
    Tissue = sample_info$tissue[match(colnames(sub), sample_info$sample)],
    CD44 = sample_info$cd44[match(colnames(sub), sample_info$sample)],
    row.names = colnames(sub)
  )
  if (!is.null(extra_ann)) {
    extra_ann <- extra_ann[match(colnames(sub), rownames(extra_ann)), , drop = FALSE]
    ann <- cbind(ann, extra_ann)
  }
  pal <- list(
    Tissue = c(MFP = "#4C78A8", LU = "#F58518", BM = "#54A24B", LN = "#E45756"),
    CD44 = c(high = "#B279A2", low = "#72B7B2", unsplit = "grey70")
  )
  draw_hm <- function() {
    args <- list(
      mat = sub, scale = "row", annotation_col = ann, annotation_colors = pal,
      show_rownames = nrow(sub) <= 80, fontsize_row = 6, main = title,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100),
      clustering_distance_cols = "euclidean"
    )
    tryCatch(
      do.call(pheatmap::pheatmap, c(args, list(clustering_distance_rows = "correlation"))),
      error = function(e) do.call(pheatmap::pheatmap, c(args, list(clustering_distance_rows = "euclidean")))
    )
  }
  grDevices::pdf(paste0(outfile, ".pdf"), width = 9, height = max(6, min(18, 0.18 * nrow(sub) + 3)))
  on.exit({
    while (grDevices::dev.cur() > 1) grDevices::dev.off()
  }, add = TRUE)
  draw_hm()
  grDevices::dev.off()
  grDevices::png(paste0(outfile, ".png"), width = 2600, height = max(1800, 40 * nrow(sub) + 400), res = 300)
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
    ggplot2::labs(title = title, x = NULL, y = "log2FC (met - MFP); negative = down in met")
  save_gg(p, outfile, width = 8, height = max(5, min(16, 0.22 * nrow(df) + 2)))
}

plot_pca <- function(heat_mat, sample_info, outfile) {
  pca <- tryCatch(stats::prcomp(t(heat_mat), scale. = TRUE), error = function(e) {
    log_msg("PCA failed: ", e$message)
    NULL
  })
  if (is.null(pca)) return(invisible(NULL))
  df <- data.frame(
    pca$x[, 1:2, drop = FALSE],
    tissue = sample_info$tissue[match(rownames(pca$x), sample_info$sample)],
    cd44 = sample_info$cd44[match(rownames(pca$x), sample_info$sample)],
    sample = rownames(pca$x)
  )
  varp <- summary(pca)$importance[2, 1:2] * 100
  p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, color = tissue, shape = cd44, label = sample)) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = "PCA after quantile log2(TPM+1)",
      x = sprintf("PC1 (%.1f%%)", varp[1]),
      y = sprintf("PC2 (%.1f%%)", varp[2])
    )
  save_gg(p, outfile)
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

analyze_down_list <- function(comp_name, de, heat_mat, sample_info, lfc_col, p_col,
                              extra_ann = NULL, do_gsea = TRUE) {
  base <- file.path(result_dir, comp_name)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  write_table(de, file.path(base, "DE_full"))
  writeLines(
    c("先看本目录 DE_full，再看 FoldChange/ 与 TopRank/。",
      "筛选：p < 0.05 且相对原位（MFP）下调（FC < 1 即入选；再加 FC < 1/1.25）。",
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
                lfc_col, p_col, fc_line = fc, extra_ann = extra_ann,
                do_ora = identical(nm, "FC_1.25") && nrow(sub) <= 800)
  }
  for (n in top_ns) {
    tag <- paste0("top", n)
    x <- de[!is.na(de[[lfc_col]]) & de[[lfc_col]] < 0, , drop = FALSE]
    x <- x[!is.na(x[[p_col]]) & x[[p_col]] < p_cutoff, , drop = FALSE]
    x <- x[order(x[[lfc_col]], decreasing = FALSE), , drop = FALSE]
    sub <- utils::head(x, n)
    emit_subset(comp_name, sub, tag, paste0(comp_name, " | downregulated top ", n),
                file.path(base, "TopRank", tag), de, heat_mat, sample_info,
                lfc_col, p_col, fc_line = 1, extra_ann = extra_ann,
                do_ora = n == 100)
  }
  if (isTRUE(do_gsea) && !isTRUE(skip_enrich)) {
    tryCatch(run_gsea_full(de, file.path(base, "00_GSEA_all_genes_NOT_FC_or_topN"), comp_name, lfc_col),
             error = function(e) log_msg("GSEA failed: ", e$message))
  }
}

# -----------------------------------------------------------------------------
# 6. limma：肺 vs 原位、骨髓 vs 原位（不把 CD44 high/low 先平均）
# -----------------------------------------------------------------------------
qc_dir <- file.path(result_dir, "00_QC")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
plot_pca(heat_mat, sample_info, file.path(qc_dir, "PCA_all_samples"))
utils::write.csv(sample_info, file.path(qc_dir, "sample_info.csv"), row.names = FALSE)

used <- sample_info$tissue %in% c("MFP", "LU", "BM")
si3 <- sample_info[used, , drop = FALSE]
log3 <- heat_mat[, si3$sample, drop = FALSE]
si3$tissue <- factor(si3$tissue, levels = c("MFP", "LU", "BM"))
design <- stats::model.matrix(~ 0 + tissue, data = si3)
colnames(design) <- levels(si3$tissue)
cont <- limma::makeContrasts(
  LU_vs_MFP = LU - MFP,
  BM_vs_MFP = BM - MFP,
  LU_vs_BM  = LU - BM,
  levels = design
)
fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(log3, design), cont))

tt_from <- function(coef) {
  tt <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
  data.frame(
    gene = rownames(tt),
    log2FC = tt$logFC,
    AveExpr = tt$AveExpr,
    t = tt$t,
    pvalue = tt$P.Value,
    padj = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
}

de_lu <- tt_from("LU_vs_MFP")
de_bm <- tt_from("BM_vs_MFP")
de_lu_bm <- tt_from("LU_vs_BM")

anno_primary <- function(de) {
  de$ensembl <- id_map$ensembl[match(de$gene, id_map$gene)]
  de$entrez <- id_map$entrez[match(de$gene, id_map$gene)]
  de$MFP_mean_tpm <- mfp_mean_tpm[de$gene]
  de$MFP_expr_percentile <- mfp_pct[de$gene]
  de$MFP_low_quartile <- !is.na(de$MFP_expr_percentile) & de$MFP_expr_percentile <= 25
  de$down_FC <- 2^(-de$log2FC)
  de
}
de_lu <- anno_primary(de_lu)
de_bm <- anno_primary(de_bm)
de_lu_bm <- anno_primary(de_lu_bm)

merged <- data.frame(
  gene = de_lu$gene,
  ensembl = de_lu$ensembl,
  entrez = de_lu$entrez,
  MFP_mean_tpm = de_lu$MFP_mean_tpm,
  MFP_expr_percentile = de_lu$MFP_expr_percentile,
  MFP_low_quartile = de_lu$MFP_low_quartile,
  log2FC_LU = de_lu$log2FC,
  pvalue_LU = de_lu$pvalue,
  padj_LU = de_lu$padj,
  down_FC_LU = de_lu$down_FC,
  log2FC_BM = de_bm$log2FC[match(de_lu$gene, de_bm$gene)],
  pvalue_BM = de_bm$pvalue[match(de_lu$gene, de_bm$gene)],
  padj_BM = de_bm$padj[match(de_lu$gene, de_bm$gene)],
  down_FC_BM = de_bm$down_FC[match(de_lu$gene, de_bm$gene)],
  log2FC_LU_vs_BM = de_lu_bm$log2FC[match(de_lu$gene, de_lu_bm$gene)],
  pvalue_LU_vs_BM = de_lu_bm$pvalue[match(de_lu$gene, de_lu_bm$gene)],
  stringsAsFactors = FALSE
)
merged$mean_met_log2FC <- (merged$log2FC_LU + merged$log2FC_BM) / 2
merged$min_down_FC <- pmin(merged$down_FC_LU, merged$down_FC_BM)
write_table(merged, file.path(result_dir, "00_QC", "limma_LU_and_BM_vs_MFP_full"))

is_down <- function(lfc, p, fc) {
  !is.na(lfc) & !is.na(p) & p < p_cutoff & lfc < 0 & (2^(-lfc) >= fc)
}

# -----------------------------------------------------------------------------
# 7. 问题 1：共享下调（肺 + 骨）
# -----------------------------------------------------------------------------
log_msg("Q1: shared down in LU and BM vs MFP")
merged$log2FC <- merged$mean_met_log2FC
merged$pvalue <- pmax(merged$pvalue_LU, merged$pvalue_BM)

shared_de <- merged
shared_keep_full <- is_down(merged$log2FC_LU, merged$pvalue_LU, 1) &
  is_down(merged$log2FC_BM, merged$pvalue_BM, 1)
shared_de$pvalue[!shared_keep_full] <- 1
analyze_down_list(
  "01_shared_lung_and_bone_down",
  shared_de, heat_mat[, si3$sample, drop = FALSE], si3,
  lfc_col = "mean_met_log2FC", p_col = "pvalue"
)

# 覆盖 FC 分层：共享必须两个器官都过同一 FC
base_shared <- file.path(result_dir, "01_shared_lung_and_bone_down")
for (nm in names(fc_cutoffs)) {
  fc <- unname(fc_cutoffs[[nm]])
  sub <- merged[is_down(merged$log2FC_LU, merged$pvalue_LU, fc) &
                  is_down(merged$log2FC_BM, merged$pvalue_BM, fc), ]
  sub <- sub[order(sub$mean_met_log2FC, decreasing = FALSE), ]
  sub$log2FC <- sub$mean_met_log2FC
  sub$pvalue <- pmax(sub$pvalue_LU, sub$pvalue_BM)
  emit_subset(
    "01_shared_lung_and_bone_down", sub, paste0(nm, "_BOTH_organs"),
    paste0("Q1 shared down | both organs FC >= ", fc),
    file.path(base_shared, "FoldChange", paste0(nm, "_BOTH_organs")),
    transform(merged, log2FC = mean_met_log2FC, pvalue = pmax(pvalue_LU, pvalue_BM)),
    heat_mat[, si3$sample, drop = FALSE], si3,
    lfc_col = "mean_met_log2FC", p_col = "pvalue", fc_line = fc,
    do_ora = identical(nm, "FC_1.25") && nrow(sub) <= 800
  )
}

# -----------------------------------------------------------------------------
# 8. 问题 2：器官特异下调
# -----------------------------------------------------------------------------
log_msg("Q2: organ-specific down")
lung_spec <- merged
lung_spec$log2FC <- lung_spec$log2FC_LU
lung_spec$pvalue <- lung_spec$pvalue_LU
# 肺特异：肺下调，且未达到骨下调标准；并要求肺比骨更低
mark_lung <- is_down(merged$log2FC_LU, merged$pvalue_LU, 1) &
  !is_down(merged$log2FC_BM, merged$pvalue_BM, 1) &
  merged$log2FC_LU < merged$log2FC_BM
lung_spec$pvalue[!mark_lung] <- 1
analyze_down_list(
  "02_lung_specific_down",
  lung_spec, heat_mat[, si3$sample, drop = FALSE], si3,
  lfc_col = "log2FC_LU", p_col = "pvalue"
)

bone_spec <- merged
bone_spec$log2FC <- bone_spec$log2FC_BM
bone_spec$pvalue <- bone_spec$pvalue_BM
mark_bone <- is_down(merged$log2FC_BM, merged$pvalue_BM, 1) &
  !is_down(merged$log2FC_LU, merged$pvalue_LU, 1) &
  merged$log2FC_BM < merged$log2FC_LU
bone_spec$pvalue[!mark_bone] <- 1
analyze_down_list(
  "03_bone_specific_down",
  bone_spec, heat_mat[, si3$sample, drop = FALSE], si3,
  lfc_col = "log2FC_BM", p_col = "pvalue"
)

base_lu <- file.path(result_dir, "02_lung_specific_down")
base_bm <- file.path(result_dir, "03_bone_specific_down")
for (nm in names(fc_cutoffs)) {
  fc <- unname(fc_cutoffs[[nm]])
  lu <- merged[is_down(merged$log2FC_LU, merged$pvalue_LU, fc) &
                 !is_down(merged$log2FC_BM, merged$pvalue_BM, fc) &
                 merged$log2FC_LU < merged$log2FC_BM, ]
  lu <- lu[order(lu$log2FC_LU), ]
  lu$log2FC <- lu$log2FC_LU
  lu$pvalue <- lu$pvalue_LU
  emit_subset(
    "02_lung_specific_down", lu, paste0(nm, "_lung_only"),
    paste0("Q2 lung-specific down | FC >= ", fc),
    file.path(base_lu, "FoldChange", paste0(nm, "_lung_only")),
    transform(merged, log2FC = log2FC_LU, pvalue = pvalue_LU),
    heat_mat[, si3$sample, drop = FALSE], si3,
    "log2FC_LU", "pvalue_LU", fc, do_ora = identical(nm, "FC_1.25") && nrow(lu) <= 800
  )
  bm <- merged[is_down(merged$log2FC_BM, merged$pvalue_BM, fc) &
                 !is_down(merged$log2FC_LU, merged$pvalue_LU, fc) &
                 merged$log2FC_BM < merged$log2FC_LU, ]
  bm <- bm[order(bm$log2FC_BM), ]
  bm$log2FC <- bm$log2FC_BM
  bm$pvalue <- bm$pvalue_BM
  emit_subset(
    "03_bone_specific_down", bm, paste0(nm, "_bone_only"),
    paste0("Q2 bone-specific down | FC >= ", fc),
    file.path(base_bm, "FoldChange", paste0(nm, "_bone_only")),
    transform(merged, log2FC = log2FC_BM, pvalue = pvalue_BM),
    heat_mat[, si3$sample, drop = FALSE], si3,
    "log2FC_BM", "pvalue_BM", fc, do_ora = identical(nm, "FC_1.25") && nrow(bm) <= 800
  )
}

if (has_pkg("ggvenn")) {
  vdf <- list(
    Lung_down = merged$gene[is_down(merged$log2FC_LU, merged$pvalue_LU, 1)],
    Bone_down = merged$gene[is_down(merged$log2FC_BM, merged$pvalue_BM, 1)]
  )
  p <- ggvenn::ggvenn(vdf, fill_color = c("#F58518", "#54A24B")) +
    ggplot2::labs(title = "p<0.05 downregulated vs MFP")
  save_gg(p, file.path(result_dir, "00_QC", "venn_lung_vs_bone_down"))
}

# 单独的肺 vs MFP、骨 vs MFP 全表（便于核对）
analyze_down_list("04_lung_vs_MFP_all_down", de_lu, heat_mat[, si3$sample, drop = FALSE], si3,
                  "log2FC", "pvalue", do_gsea = FALSE)
analyze_down_list("05_bone_vs_MFP_all_down", de_bm, heat_mat[, si3$sample, drop = FALSE], si3,
                  "log2FC", "pvalue", do_gsea = FALSE)


# -----------------------------------------------------------------------------
# 9. 问题 3：施旺细胞 / 神经营养因子 / 轴突导向 分别打分
# -----------------------------------------------------------------------------
log_msg("Q3: three neural signatures (Schwann, neurotrophin, axon guidance)")
pni_dir <- file.path(result_dir, "06_neural_invasion")
dir.create(pni_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  c("没有病理神经浸润标签。三套 marker 分开打分：",
    "  1) 施旺细胞 Schwann",
    "  2) 神经营养因子 Neurotrophin",
    "  3) 轴突导向 Axon_guidance",
    "每套：样品分数、原位 MFP 标记、与分数负相关的基因（p<0.05, rho<0）。",
    "另外比较三套分数在原位 vs 肺、原位 vs 骨之间的差异，并与肺/骨下调基因取交集。"),
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
  unique(msig_bp$gene_symbol[grepl(pattern, msig_bp$gs_name)])
}

curated <- list(
  Schwann = c(
    "Sox10", "S100b", "Mag", "Mpz", "Mbp", "Pmp22", "Plp1", "Gfap", "Egr2",
    "Pou3f1", "Pou3f2", "Ngfr", "Ncam1", "L1cam", "Cdh19", "Dhh", "Mpzl1",
    "Prx", "Gjb1", "Mal", "Cnp", "Erbb3", "Erbb2", "Nrg1", "Mpz", "Oct6"
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
  intersect(g, rownames(heat_mat))
})
names(sig_genes) <- names(curated)
for (nm in names(sig_genes)) {
  log_msg("Signature ", nm, ": ", length(sig_genes[[nm]]), " genes in matrix")
  utils::write.csv(
    data.frame(gene = sig_genes[[nm]]),
    file.path(pni_dir, paste0("signature_genes_", nm, ".csv")),
    row.names = FALSE
  )
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
      for (nm in rownames(gsva_mat)) {
        score_mat[, nm] <- as.numeric(gsva_mat[nm, colnames(heat_mat)])
      }
    }
  }
}
score_mat <- as.data.frame(score_mat)
score_mat$sample <- colnames(heat_mat)
score_tab <- merge(sample_info, score_mat, by = "sample", all.x = TRUE)
write_table(score_tab, file.path(pni_dir, "sample_three_neural_scores"))
mfp_mark <- score_tab[score_tab$tissue == "MFP", ]
write_table(mfp_mark, file.path(pni_dir, "MFP_neural_invasion_marker"))

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

p_box <- ggplot2::ggplot(score_long, ggplot2::aes(x = tissue, y = score, fill = tissue)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  ggplot2::geom_jitter(width = 0.15, size = 2) +
  ggplot2::facet_wrap(~ signature, scales = "free_y") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(title = "Neural scores by tissue (MFP / LU / BM / LN)", x = NULL, y = "score")
save_gg(p_box, file.path(pni_dir, "three_scores_by_tissue"), 9, 5)

annotate_cor <- function(df) {
  df$padj <- stats::p.adjust(df$pvalue, method = "BH")
  df$ensembl <- id_map$ensembl[match(df$gene, id_map$gene)]
  df$MFP_mean_tpm <- mfp_mean_tpm[df$gene]
  df$MFP_expr_percentile <- mfp_pct[df$gene]
  df$log2FC_LU <- merged$log2FC_LU[match(df$gene, merged$gene)]
  df$pvalue_LU <- merged$pvalue_LU[match(df$gene, merged$gene)]
  df$log2FC_BM <- merged$log2FC_BM[match(df$gene, merged$gene)]
  df$pvalue_BM <- merged$pvalue_BM[match(df$gene, merged$gene)]
  df$down_in_lung <- is_down(df$log2FC_LU, df$pvalue_LU, 1)
  df$down_in_bone <- is_down(df$log2FC_BM, df$pvalue_BM, 1)
  df[order(df$rho), ]
}

spearman_vs_score <- function(score_vec) {
  do.call(rbind, lapply(rownames(heat_mat), function(g) {
    ct <- suppressWarnings(stats::cor.test(
      as.numeric(heat_mat[g, ]), score_vec[colnames(heat_mat)],
      method = "spearman", exact = FALSE
    ))
    data.frame(gene = g, rho = unname(ct$estimate), pvalue = ct$p.value, stringsAsFactors = FALSE)
  }))
}

tissue_ttest <- function(score_vec, grp_a, grp_b) {
  a <- score_vec[sample_info$sample[sample_info$tissue == grp_a]]
  b <- score_vec[sample_info$sample[sample_info$tissue == grp_b]]
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) < 2 || length(b) < 2) {
    return(data.frame(group_a = grp_a, group_b = grp_b, n_a = length(a), n_b = length(b),
                      mean_a = mean(a), mean_b = mean(b), delta = mean(b) - mean(a),
                      pvalue = NA_real_))
  }
  tt <- stats::t.test(b, a)
  data.frame(
    group_a = grp_a, group_b = grp_b, n_a = length(a), n_b = length(b),
    mean_a = mean(a), mean_b = mean(b), delta = unname(diff(tt$estimate)),
    pvalue = tt$p.value
  )
}

neg_lists <- list()
assoc_rows <- list()
lung_down_genes <- merged$gene[is_down(merged$log2FC_LU, merged$pvalue_LU, 1)]
bone_down_genes <- merged$gene[is_down(merged$log2FC_BM, merged$pvalue_BM, 1)]

for (nm in names(curated)) {
  log_msg("Neural axis: ", nm)
  ax_dir <- file.path(pni_dir, nm)
  dir.create(ax_dir, recursive = TRUE, showWarnings = FALSE)
  sc <- setNames(score_tab[[nm]], score_tab$sample)
  plot_heatmap(
    heat_mat, sample_info, sig_genes[[nm]],
    paste(nm, "signature genes"),
    file.path(ax_dir, paste0("heatmap_", nm, "_signature")),
    extra_ann = data.frame(score = sc[sample_info$sample], row.names = sample_info$sample)
  )

  cor_df <- annotate_cor(spearman_vs_score(sc))
  write_table(cor_df, file.path(ax_dir, paste0(nm, "_gene_vs_score_spearman")))
  neg <- cor_df[cor_df$rho < 0 & cor_df$pvalue < p_cutoff, ]
  neg$log2FC <- neg$rho
  neg_lists[[nm]] <- neg
  write_table(neg, file.path(ax_dir, paste0(nm, "_neg_cor_p0.05")))
  emit_subset(
    paste0("06_neural_invasion/", nm), neg, paste0(nm, "_neg_cor_p0.05"),
    paste0(nm, " | negative correlation with score (p<0.05)"),
    file.path(ax_dir, "neg_cor_p0.05"),
    transform(cor_df, log2FC = rho), heat_mat, sample_info,
    lfc_col = "rho", p_col = "pvalue", do_ora = nrow(neg) <= 800
  )
  emit_subset(
    paste0("06_neural_invasion/", nm), utils::head(neg, 100), paste0(nm, "_top100"),
    paste0(nm, " | top100 negative correlation"),
    file.path(ax_dir, "TopRank", "top100"),
    transform(cor_df, log2FC = rho), heat_mat, sample_info,
    lfc_col = "rho", p_col = "pvalue", do_ora = TRUE
  )
  if (!isTRUE(skip_enrich)) {
    tryCatch(
      run_gsea_full(transform(cor_df, log2FC = rho), file.path(ax_dir, "00_GSEA_all_genes_NOT_FC_or_topN"),
                    paste(nm, "negative correlation"), stat_col = "rho"),
      error = function(e) log_msg(nm, " GSEA failed: ", e$message)
    )
  }

  lu_ov <- neg[neg$gene %in% lung_down_genes, ]
  bm_ov <- neg[neg$gene %in% bone_down_genes, ]
  both_ov <- neg[neg$gene %in% intersect(lung_down_genes, bone_down_genes), ]
  write_table(lu_ov, file.path(ax_dir, paste0(nm, "_neg_cor_AND_lung_down")))
  write_table(bm_ov, file.path(ax_dir, paste0(nm, "_neg_cor_AND_bone_down")))
  write_table(both_ov, file.path(ax_dir, paste0(nm, "_neg_cor_AND_shared_lung_bone_down")))

  tt_lu <- tissue_ttest(sc, "MFP", "LU")
  tt_bm <- tissue_ttest(sc, "MFP", "BM")
  tt_lu_bm <- tissue_ttest(sc, "LU", "BM")
  assoc <- rbind(
    data.frame(signature = nm, contrast = "LU_vs_MFP", tt_lu, stringsAsFactors = FALSE),
    data.frame(signature = nm, contrast = "BM_vs_MFP", tt_bm, stringsAsFactors = FALSE),
    data.frame(signature = nm, contrast = "BM_vs_LU", tt_lu_bm, stringsAsFactors = FALSE)
  )
  write_table(assoc, file.path(ax_dir, paste0(nm, "_score_vs_lung_bone_ttest")))
  assoc_rows[[nm]] <- assoc

  # 签名基因本身在肺/骨 vs 原位的 limma
  sg <- intersect(sig_genes[[nm]], merged$gene)
  sig_de <- merged[merged$gene %in% sg, ]
  write_table(sig_de, file.path(ax_dir, paste0(nm, "_signature_genes_limma_vs_MFP")))
}

assoc_all <- do.call(rbind, assoc_rows)
write_table(assoc_all, file.path(pni_dir, "three_signatures_score_vs_lung_bone"))

p_assoc <- ggplot2::ggplot(assoc_all[assoc_all$contrast %in% c("LU_vs_MFP", "BM_vs_MFP"), ],
                           ggplot2::aes(x = signature, y = delta, fill = contrast)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.65) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    title = "Neural score change vs primary (positive = higher in metastasis)",
    y = "mean(met) - mean(MFP)", x = NULL
  )
save_gg(p_assoc, file.path(pni_dir, "three_signatures_delta_lung_bone"), 8, 5)

if (has_pkg("ggvenn")) {
  v3 <- lapply(neg_lists, function(x) x$gene)
  names(v3) <- names(neg_lists)
  p <- ggvenn::ggvenn(v3, fill_color = c("#4C78A8", "#F58518", "#54A24B")) +
    ggplot2::labs(title = "Negatively correlated genes (p<0.05) across three neural scores")
  save_gg(p, file.path(pni_dir, "venn_three_neg_cor"), 8, 6)
}

# -----------------------------------------------------------------------------
# 10. 三问总表
# -----------------------------------------------------------------------------
sum_dir <- file.path(result_dir, "07_summary")
dir.create(sum_dir, recursive = TRUE, showWarnings = FALSE)

q1_fc1 <- merged[is_down(merged$log2FC_LU, merged$pvalue_LU, 1) &
                   is_down(merged$log2FC_BM, merged$pvalue_BM, 1), ]
q1_fc1 <- q1_fc1[order(q1_fc1$mean_met_log2FC), ]
q1 <- merged[is_down(merged$log2FC_LU, merged$pvalue_LU, 1.25) &
               is_down(merged$log2FC_BM, merged$pvalue_BM, 1.25), ]
q1 <- q1[order(q1$mean_met_log2FC), ]
q2l <- merged[is_down(merged$log2FC_LU, merged$pvalue_LU, 1.25) &
                !is_down(merged$log2FC_BM, merged$pvalue_BM, 1.25) &
                merged$log2FC_LU < merged$log2FC_BM, ]
q2l <- q2l[order(q2l$log2FC_LU), ]
q2b <- merged[is_down(merged$log2FC_BM, merged$pvalue_BM, 1.25) &
                !is_down(merged$log2FC_LU, merged$pvalue_LU, 1.25) &
                merged$log2FC_BM < merged$log2FC_LU, ]
q2b <- q2b[order(q2b$log2FC_BM), ]

write_table(q1_fc1, file.path(sum_dir, "Q1_shared_lung_and_bone_down_p0.05_FC1"))
write_table(q1, file.path(sum_dir, "Q1_shared_lung_and_bone_down_p0.05_FC1.25"))
write_table(q2l, file.path(sum_dir, "Q2_lung_specific_down_p0.05_FC1.25"))
write_table(q2b, file.path(sum_dir, "Q2_bone_specific_down_p0.05_FC1.25"))
write_table(mfp_mark, file.path(sum_dir, "Q3_MFP_three_neural_scores"))
write_table(assoc_all, file.path(sum_dir, "Q3_neural_scores_vs_lung_bone"))

q3_counts <- data.frame(
  signature = names(neg_lists),
  n_neg_cor_p0.05 = vapply(neg_lists, nrow, integer(1)),
  n_neg_and_lung_down = vapply(neg_lists, function(x) sum(x$gene %in% lung_down_genes), integer(1)),
  n_neg_and_bone_down = vapply(neg_lists, function(x) sum(x$gene %in% bone_down_genes), integer(1)),
  stringsAsFactors = FALSE
)
write_table(q3_counts, file.path(sum_dir, "Q3_neg_cor_overlap_with_lung_bone"))
for (nm in names(neg_lists)) {
  write_table(utils::head(neg_lists[[nm]], 500),
              file.path(sum_dir, paste0("Q3_", nm, "_neg_cor_p0.05_top500")))
}

n_tab <- data.frame(
  question = c(
    "Q1 shared down p<0.05 FC<1",
    "Q1 shared down p<0.05 FC<1/1.25",
    "Q2 lung-specific p<0.05 FC<1/1.25",
    "Q2 bone-specific p<0.05 FC<1/1.25",
    paste0("Q3 ", q3_counts$signature, " neg-cor p<0.05")
  ),
  n_genes = c(nrow(q1_fc1), nrow(q1), nrow(q2l), nrow(q2b), q3_counts$n_neg_cor_p0.05)
)
write_table(n_tab, file.path(sum_dir, "gene_counts"))
log_msg("Done. Summary counts:")
print(n_tab)
log_msg("Open: ", sum_dir)
