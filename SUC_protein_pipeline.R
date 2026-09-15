#!/usr/bin/env Rscript
# =============================================================================
# SUC 蛋白质组分析流程（DIA-NN）
# 样品 1–16 独立保留，禁止合并。
# 六类蛋白（每类单独出表、单独作图）：
#   A) 2 vs 1 上调，且 3 vs 4 无变化
#   B) 5 vs 6 下调，且 7 vs 8 无变化
#   C) 9 vs 10 上调，且 11 vs 12 无变化
#   D) 13 vs 14 下调，且 15 vs 16 无变化
#   E) 同一档位下 A ∩ B
#   F) 同一档位下 C ∩ D
# 预处理：16 个样品一起过滤低丰度 + log2 分位数标准化，再用标准化值算 FC。
# 1-vs-1 无重复：不伪造 p 值。
# 分层：FC >= 1 / 1.25 / 1.5 / 2（下调为倒数）以及 top 50–300。
# 每个非空子集：差异表、火山图、热图、ORA GO / 通路 / KEGG、GSEA，
# 以及 mitochondria 文本内通路的专项 GO / KEGG / 通路图。
# 「GWAS」按 GSEA（基因集富集）实现，不是 SNP 全基因组关联。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")
options(clusterProfiler.download.method = "auto")

# -----------------------------------------------------------------------------
# 0. 依赖包
# -----------------------------------------------------------------------------
cran_required <- c(
  "dplyr", "tidyr", "tibble", "stringr", "ggplot2",
  "ggrepel", "pheatmap", "RColorBrewer", "matrixStats", "cowplot",
  "ggridges", "ggnewscale", "igraph"
)
cran_optional <- c("ggvenn", "writexl")
bioc_required <- c(
  "limma", "clusterProfiler", "org.Hs.eg.db",
  "enrichplot", "DOSE", "AnnotationDbi", "fgsea", "msigdbr"
)
bioc_optional <- c("ReactomePA", "pathview", "GO.db")

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
  env_dir <- Sys.getenv("SUC_PROTEIN_DIR", unset = "")
  candidates <- c(
    env_dir,
    "E:/R/SUC-protein",
    "E:\\R\\SUC-protein",
    file.path(getwd(), "SUC-protein"),
    getwd()
  )
  candidates <- unique(candidates[nzchar(candidates)])
  wanted <- c("report.pg_matrix", "report.pg_matrix.tsv",
              "report.pr_matrix", "report.pr_matrix.tsv")
  for (d in candidates) {
    if (!dir.exists(d)) next
    hits <- file.exists(file.path(d, wanted))
    if (any(hits)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
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

fc_cutoffs <- c("FC_1" = 1, "FC_1.25" = 1.25, "FC_1.5" = 1.5, "FC_2" = 2)
top_ns <- c(50, 75, 100, 150, 200, 250, 300)
unchanged_max_fc <- 1.25
pseudo <- 1

# -----------------------------------------------------------------------------
# 2. 读入 DIA-NN 矩阵与 mitochondria 通路表
# -----------------------------------------------------------------------------
first_existing <- function(dir, names) {
  for (nm in names) {
    p <- file.path(dir, nm)
    if (file.exists(p)) return(p)
  }
  NA_character_
}

clean_symbol <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("\\s+", "", x)
  vapply(x, function(z) {
    if (!nzchar(z) || toupper(z) %in% c("NA", "NULL", "N/A")) return(NA_character_)
    parts <- unlist(strsplit(z, "[;,/|]"))
    parts <- parts[nzchar(parts)]
    if (length(parts) == 0) return(NA_character_)
    parts[[1]]
  }, character(1), USE.NAMES = FALSE)
}

extract_sample_number <- function(colname) {
  b <- basename(as.character(colname))
  b <- sub("\\.(raw|mzML|mzml|dia|parquet|tsv|txt)$", "", b, ignore.case = TRUE)
  if (grepl("^(sample[_-]?)?(1[0-6]|[1-9])$", b, ignore.case = TRUE)) {
    return(as.integer(sub("^sample[_-]?", "", b, ignore.case = TRUE)))
  }
  m <- gregexpr("(?<![0-9])(1[0-6]|[1-9])(?![0-9])", b, perl = TRUE)[[1]]
  if (length(m) == 0 || m[1] == -1) return(NA_integer_)
  last <- m[length(m)]
  len <- attr(m, "match.length")[length(m)]
  as.integer(substr(b, last, last + len - 1))
}

is_meta_column <- function(nm) {
  grepl(
    paste0(
      "^(protein|gene|first\\.|description|pg\\.|pr\\.|modified|stripped|",
      "precursor|peptide|proteotypic|q\\.value|quantity|ids?)$"
    ),
    gsub("[^A-Za-z0-9]+", ".", nm),
    ignore.case = TRUE
  ) || grepl(
    "protein\\.(group|ids|names)|gene|description|protein\\.names|first\\.protein",
    nm, ignore.case = TRUE
  )
}

pick_id_column <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

read_diann_matrix <- function(path, aggregate_precursor = FALSE) {
  log_msg("Reading: ", path)
  df <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
  if (ncol(df) < 3) stop("矩阵列数过少: ", path)
  gene_col <- pick_id_column(df, c("Genes", "Gene", "Gene.Names", "Gene.Name", "genes"))
  prot_col <- pick_id_column(df, c("Protein.Group", "Protein.Ids", "Protein.Id", "Protein.Names", "Protein.Group.Ids"))
  if (is.na(prot_col)) prot_col <- names(df)[1]
  sample_cols <- names(df)[!vapply(names(df), is_meta_column, logical(1))]
  if (length(sample_cols) == 0) {
    sample_cols <- names(df)[!names(df) %in% c(gene_col, prot_col)]
  }
  nums <- vapply(sample_cols, extract_sample_number, integer(1))
  mapped <- !is.na(nums) & nums >= 1 & nums <= 16
  if (sum(unique(nums[mapped])) < 16 && length(sample_cols) == 16 && all(is.na(nums))) {
    log_msg("WARNING: 未能从列名解析 1–16，按矩阵中 16 个定量列的从左到右顺序编号")
    nums <- seq_along(sample_cols)
    mapped <- rep(TRUE, length(sample_cols))
  }
  if (sum(mapped) < 16) {
    log_msg("WARNING: 只解析到 ", sum(mapped), " 个样品列: ",
            paste(paste0(sample_cols[mapped], "=", nums[mapped]), collapse = "; "))
  }
  keep_cols <- sample_cols[mapped]
  keep_nums <- as.integer(nums[mapped])
  if (anyDuplicated(keep_nums)) {
    log_msg("WARNING: 样品编号重复，保留每个编号的第一列")
    ok <- !duplicated(keep_nums)
    keep_cols <- keep_cols[ok]
    keep_nums <- keep_nums[ok]
  }
  mat <- as.matrix(df[, keep_cols, drop = FALSE])
  storage.mode(mat) <- "double"
  colnames(mat) <- as.character(keep_nums)
  symbols <- if (!is.na(gene_col)) clean_symbol(df[[gene_col]]) else rep(NA_character_, nrow(df))
  prot <- as.character(df[[prot_col]])
  uniprot <- NA_character_
  uid_col <- pick_id_column(df, c("Protein.Ids", "Protein.Id", "Protein.Group"))
  if (!is.na(uid_col)) uniprot <- clean_symbol(df[[uid_col]])
  gene <- symbols
  gene[is.na(gene) | !nzchar(gene)] <- prot[is.na(gene) | !nzchar(gene)]
  gene[is.na(gene) | !nzchar(gene)] <- paste0("PROT_", seq_len(nrow(df)))[is.na(gene) | !nzchar(gene)]
  if (isTRUE(aggregate_precursor)) {
    log_msg("Aggregating precursor rows by protein/gene (max intensity)")
    tmp <- as.data.frame(mat, check.names = FALSE)
    tmp$gene <- gene
    tmp$symbol <- symbols
    tmp$uniprot <- uniprot
    tmp$prot <- prot
    first_nz <- function(x) {
      x <- x[!is.na(x) & nzchar(x)]
      if (length(x) == 0) NA_character_ else x[[1]]
    }
    agg <- dplyr::summarise(
      dplyr::group_by(tmp, gene),
      dplyr::across(dplyr::all_of(colnames(mat)), function(x) max(x, na.rm = TRUE)),
      symbol = first_nz(symbol),
      uniprot = first_nz(uniprot),
      prot = prot[[1]],
      .groups = "drop"
    )
    mat <- as.matrix(agg[, colnames(mat), drop = FALSE])
    storage.mode(mat) <- "double"
    mat[!is.finite(mat)] <- NA_real_
    rownames(mat) <- agg$gene
    meta <- data.frame(
      gene = agg$gene, symbol = agg$symbol, uniprot = agg$uniprot,
      protein_id = agg$prot, stringsAsFactors = FALSE
    )
  } else {
    if (anyDuplicated(gene)) {
      log_msg("Duplicate protein labels: keep the row with highest mean intensity")
      ord <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
      mat <- mat[ord, , drop = FALSE]
      gene <- gene[ord]
      symbols <- symbols[ord]
      uniprot <- uniprot[ord]
      prot <- prot[ord]
      keep <- !duplicated(gene)
      mat <- mat[keep, , drop = FALSE]
      gene <- gene[keep]
      symbols <- symbols[keep]
      uniprot <- uniprot[keep]
      prot <- prot[keep]
    }
    rownames(mat) <- gene
    meta <- data.frame(
      gene = gene, symbol = symbols, uniprot = uniprot,
      protein_id = prot, stringsAsFactors = FALSE
    )
  }
  miss <- setdiff(as.character(1:16), colnames(mat))
  if (length(miss) > 0) {
    log_msg("WARNING: 缺少样品列 ", paste(miss, collapse = ", "))
  }
  mat <- mat[, intersect(as.character(1:16), colnames(mat)), drop = FALSE]
  sample_info <- data.frame(
    sample = colnames(mat),
    group = paste0("S", colnames(mat)),
    stringsAsFactors = FALSE
  )
  list(mat = mat, meta = meta, sample_info = sample_info, source = basename(path))
}

load_protein_matrix <- function(project_dir) {
  pg <- first_existing(project_dir, c("report.pg_matrix", "report.pg_matrix.tsv"))
  if (!is.na(pg)) {
    log_msg("Using protein-group matrix: ", basename(pg))
    return(read_diann_matrix(pg, aggregate_precursor = FALSE))
  }
  pr <- first_existing(project_dir, c("report.pr_matrix", "report.pr_matrix.tsv"))
  if (!is.na(pr)) {
    log_msg("Protein-group matrix missing; aggregating precursor matrix")
    return(read_diann_matrix(pr, aggregate_precursor = TRUE))
  }
  stop("未找到 report.pg_matrix / report.pr_matrix: ", project_dir)
}

parse_mitochondria_file <- function(project_dir) {
  path <- first_existing(project_dir, c("mitochondria", "mitochondria.txt", "mitochondria.tsv", "mitochondria.csv"))
  if (is.na(path)) {
    log_msg("WARNING: 未找到 mitochondria 文本，专项分析将用线粒体关键词回退")
    return(list(path = NA_character_, go_ids = character(), names = character(), raw = character()))
  }
  raw <- readLines(path, warn = FALSE, encoding = "UTF-8")
  raw <- trimws(raw)
  raw <- raw[nzchar(raw) & !startsWith(raw, "#")]
  ids <- unique(unlist(regmatches(raw, gregexpr("GO:[0-9]{7}", raw, ignore.case = TRUE))))
  ids <- toupper(ids)
  names_only <- gsub("GO:[0-9]{7}", "", raw, ignore.case = TRUE)
  names_only <- trimws(gsub("[,;|/]+$", "", names_only))
  names_only <- names_only[nzchar(names_only)]
  log_msg("mitochondria file: ", basename(path), " | lines=", length(raw),
          " GO IDs=", length(ids), " names=", length(names_only))
  list(path = path, go_ids = ids, names = unique(names_only), raw = raw)
}

# -----------------------------------------------------------------------------
# 3. 过滤 + 标准化（全部 16 个样品一次完成）
# -----------------------------------------------------------------------------
filter_low_abundance <- function(mat) {
  n <- ncol(mat)
  detected <- rowSums(is.finite(mat) & mat > 0, na.rm = TRUE)
  row_sum <- rowSums(pmax(mat, 0), na.rm = TRUE)
  keep <- detected >= max(2, floor(n / 4)) & row_sum > 0
  if (sum(keep) < 200) {
    keep <- detected >= 1 & row_sum > 0
    log_msg("Strict abundance filter left too few proteins; fallback to detected-in >=1 sample")
  }
  log_msg("Low-abundance filter: keep ", sum(keep), " / ", nrow(mat), " proteins")
  mat[keep, , drop = FALSE]
}

normalize_intensity <- function(mat) {
  mat[!is.finite(mat)] <- NA_real_
  log_mat <- log2(pmax(mat, 0) + pseudo)
  log_msg("Quantile normalize log2(intensity + ", pseudo, ")")
  # limma quantile 不能有全 NA 列
  for (j in seq_len(ncol(log_mat))) {
    if (all(!is.finite(log_mat[, j]))) stop("样品 ", colnames(log_mat)[j], " 全部缺失")
  }
  med <- apply(log_mat, 2, stats::median, na.rm = TRUE)
  log_imp <- log_mat
  for (j in seq_len(ncol(log_imp))) {
    miss <- !is.finite(log_imp[, j])
    if (any(miss)) log_imp[miss, j] <- med[j]
  }
  qn <- limma::normalizeBetweenArrays(log_imp, method = "quantile")
  qn[!is.finite(log_mat)] <- NA_real_
  list(log_mat = qn, heat_mat = qn, method = "quantile_log2")
}

# -----------------------------------------------------------------------------
# 4. 1-vs-1 FC（不伪造 p）
# -----------------------------------------------------------------------------
pairwise_fc <- function(log_mat, treat, ctrl, comp_name) {
  treat <- as.character(treat)
  ctrl <- as.character(ctrl)
  if (!all(c(treat, ctrl) %in% colnames(log_mat))) {
    log_msg("SKIP ", comp_name, ": missing samples ", treat, " or ", ctrl)
    return(NULL)
  }
  log_msg(comp_name, " : sample ", treat, " vs ", ctrl, " (1-vs-1, FC only, no p-value)")
  log2FC <- log_mat[, treat] - log_mat[, ctrl]
  ave <- (log_mat[, treat] + log_mat[, ctrl]) / 2
  data.frame(
    gene = rownames(log_mat),
    log2FC = as.numeric(log2FC),
    FoldChange = 2^as.numeric(log2FC),
    AveExpr = as.numeric(ave),
    pvalue = NA_real_,
    padj = NA_real_,
    treat_sample = treat,
    ctrl_sample = ctrl,
    stringsAsFactors = FALSE
  )
}

is_unchanged <- function(log2fc, max_fc = unchanged_max_fc) {
  fc <- 2^log2fc
  is.finite(fc) & fc > (1 / max_fc) & fc < max_fc
}

# -----------------------------------------------------------------------------
# 5. 基因 ID 转换
# -----------------------------------------------------------------------------
map_to_entrez <- function(symbols, uniprot = NULL) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  out <- data.frame(gene = character(), entrez = character(), stringsAsFactors = FALSE)
  if (length(symbols) == 0 && (is.null(uniprot) || all(is.na(uniprot)))) return(out)
  m <- tryCatch(
    clusterProfiler::bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
  )
  if (nrow(m) > 0) {
    m <- m[!duplicated(m[[1]]), ]
    out <- data.frame(gene = m[[1]], entrez = m[[2]], stringsAsFactors = FALSE)
  }
  leftover <- setdiff(symbols, out$gene)
  if (length(leftover) > 0) {
    m2 <- tryCatch(
      clusterProfiler::bitr(leftover, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
      error = function(e) data.frame()
    )
    if (nrow(m2) > 0) {
      out <- rbind(out, data.frame(gene = m2[[1]], entrez = m2[[2]], stringsAsFactors = FALSE))
    }
  }
  if (!is.null(uniprot)) {
    uni <- unique(as.character(uniprot))
    uni <- uni[!is.na(uni) & nzchar(uni)]
    if (length(uni) > 0) {
      m3 <- tryCatch(
        clusterProfiler::bitr(uni, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
        error = function(e) data.frame()
      )
      if (nrow(m3) > 0) {
        out <- rbind(out, data.frame(gene = m3[[1]], entrez = m3[[2]], stringsAsFactors = FALSE))
      }
    }
  }
  out[!duplicated(out$gene), , drop = FALSE]
}

ranked_entrez <- function(de) {
  mp <- map_to_entrez(de$gene, if ("uniprot" %in% names(de)) de$uniprot else NULL)
  de2 <- merge(de, mp, by = "gene")
  de2 <- de2[!is.na(de2$entrez) & !is.na(de2$log2FC), ]
  de2 <- de2[order(abs(de2$log2FC), decreasing = TRUE), ]
  de2 <- de2[!duplicated(de2$entrez), ]
  stats <- de2$log2FC
  names(stats) <- de2$entrez
  sort(stats, decreasing = TRUE)
}

# -----------------------------------------------------------------------------
# 6. 绘图工具
# -----------------------------------------------------------------------------
save_gg <- function(plot, path_stub, width = 8, height = 6) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".pdf"), plot, width = width, height = height),
           error = function(e) log_msg("pdf ggsave failed: ", e$message))
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".png"), plot, width = width, height = height, dpi = 300),
           error = function(e) log_msg("png ggsave failed: ", e$message))
}

try_save_plot <- function(fun, stub, width = 9, height = 7) {
  p <- tryCatch(fun(), error = function(e) {
    log_msg("Plot failed (", basename(stub), "): ", e$message)
    NULL
  })
  if (is.null(p)) return(invisible(FALSE))
  ok <- tryCatch({
    save_gg(p, stub, width = width, height = height)
    TRUE
  }, error = function(e) {
    log_msg("ggsave failed (", basename(stub), "): ", e$message)
    FALSE
  })
  ok
}

note_empty <- function(stub, msg) {
  dir.create(dirname(stub), recursive = TRUE, showWarnings = FALSE)
  writeLines(msg, paste0(stub, "_EMPTY.txt"))
}

plot_volcano <- function(de, highlight, title, outfile, fc_line = 1) {
  df <- de
  df$y <- df$AveExpr
  ylab <- "Average log2 intensity (1-vs-1, no p-value)"
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
  save_gg(p, outfile)
}

plot_contrast_scatter <- function(change_de, stable_de, highlight, title, outfile,
                                  unchanged_fc = unchanged_max_fc) {
  if (is.null(change_de) || is.null(stable_de)) return(invisible(NULL))
  df <- merge(
    change_de[, c("gene", "log2FC")],
    stable_de[, c("gene", "log2FC")],
    by = "gene", suffixes = c("_change", "_stable")
  )
  df$set <- ifelse(df$gene %in% highlight, "selected", "other")
  df$label <- ifelse(df$gene %in% utils::head(highlight, 15), df$gene, NA)
  band <- log2(unchanged_fc)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC_change, y = log2FC_stable, color = set)) +
    ggplot2::geom_point(alpha = 0.75, size = 1.5) +
    ggplot2::geom_hline(yintercept = c(-band, band), linetype = 2, color = "grey40") +
    ggplot2::geom_vline(xintercept = 0, linetype = 3, color = "grey60") +
    ggplot2::scale_color_manual(values = c(other = "grey70", selected = "#D62828")) +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = title,
      x = "log2FC (change contrast)",
      y = paste0("log2FC (stable contrast; unchanged |FC| < ", unchanged_fc, ")"),
      color = NULL
    )
  save_gg(p, outfile)
}

plot_heatmap <- function(heat_mat, sample_info, genes, title, outfile, samples = NULL) {
  genes <- intersect(genes, rownames(heat_mat))
  if (length(genes) > 200) {
    log_msg("Heatmap truncated to 200 proteins from ", length(genes), ": ", title)
    genes <- genes[seq_len(200)]
  }
  if (length(genes) < 2) {
    log_msg("Heatmap skipped (<2 proteins): ", title)
    note_empty(outfile, "fewer than 2 proteins")
    return(invisible(NULL))
  }
  cols <- if (is.null(samples)) colnames(heat_mat) else intersect(as.character(samples), colnames(heat_mat))
  if (length(cols) < 2) {
    note_empty(outfile, "fewer than 2 samples")
    return(invisible(NULL))
  }
  sub <- heat_mat[genes, cols, drop = FALSE]
  sub[!is.finite(sub)] <- stats::median(sub[is.finite(sub)])
  ann <- data.frame(Sample = sample_info$group, row.names = sample_info$sample)
  ann <- ann[colnames(sub), , drop = FALSE]
  draw_hm <- function() {
    args <- list(
      mat = sub, scale = "row", annotation_col = ann,
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
  x <- heat_mat
  x[!is.finite(x)] <- stats::median(x[is.finite(x)])
  pca <- tryCatch(stats::prcomp(t(x), scale. = TRUE), error = function(e) {
    log_msg("PCA failed: ", e$message)
    NULL
  })
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
      title = "PCA after log2 quantile normalization",
      x = sprintf("PC1 (%.1f%%)", varp[1]),
      y = sprintf("PC2 (%.1f%%)", varp[2])
    )
  save_gg(p, outfile)
}

plot_de_bar <- function(sub, title, outfile) {
  if (nrow(sub) == 0) return(invisible(NULL))
  df <- sub[order(sub$log2FC, decreasing = TRUE), , drop = FALSE]
  if (nrow(df) > 60) df <- rbind(utils::head(df, 30), utils::tail(df, 30))
  df$gene <- factor(df$gene, levels = rev(unique(df$gene)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = gene, y = log2FC)) +
    ggplot2::geom_col(fill = "#D62828", width = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title, x = NULL, y = "log2 Fold Change")
  save_gg(p, outfile, width = 8, height = max(5, min(16, 0.22 * nrow(df) + 2)))
}

# -----------------------------------------------------------------------------
# 7. mitochondria 专项匹配（不改全库 p 值）
# -----------------------------------------------------------------------------
.mito_env <- new.env(parent = emptyenv())

term_matches_mito <- function(id, desc, mito) {
  id <- as.character(id)[1]
  desc <- as.character(desc)[1]
  if (is.na(id)) id <- ""
  if (is.na(desc)) desc <- ""
  id_u <- toupper(id)
  txt <- tolower(paste(id, desc))
  if (length(mito$go_ids) > 0 && id_u %in% toupper(mito$go_ids)) return(TRUE)
  needles <- unique(c(mito$names, mito$raw))
  needles <- needles[nzchar(needles)]
  if (length(needles) > 0) {
    hit <- vapply(needles, function(x) {
      x <- tolower(x)
      if (!nzchar(x) || nchar(x) < 4) return(FALSE)
      grepl(x, txt, fixed = TRUE) || grepl(txt, x, fixed = TRUE)
    }, logical(1))
    if (any(hit)) return(TRUE)
  }
  if (length(mito$go_ids) == 0 && length(mito$names) == 0) {
    return(grepl("mitochondr|oxidative phosphorylat|respiratory chain|tca cycle|electron transport|oxphos", txt))
  }
  FALSE
}

entrez_for_go <- function(go_id) {
  ids <- tryCatch(
    AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = go_id, column = "ENTREZID",
      keytype = "GOALL", multiVals = "list"
    )[[go_id]],
    error = function(e) {
      tryCatch(
        AnnotationDbi::mapIds(
          org.Hs.eg.db, keys = go_id, column = "ENTREZID",
          keytype = "GO", multiVals = "list"
        )[[go_id]],
        error = function(e2) character()
      )
    }
  )
  unique(as.character(ids[!is.na(ids)]))
}

append_term2gene <- function(lst, name, entrez) {
  entrez <- unique(as.character(entrez))
  entrez <- entrez[nzchar(entrez) & !is.na(entrez)]
  if (length(entrez) < 5) return(lst)
  lst[[length(lst) + 1]] <- data.frame(
    gs_name = name, entrez = entrez, stringsAsFactors = FALSE
  )
  lst
}

get_mito_term2gene <- function(mito) {
  if (!is.null(.mito_env$term2gene)) return(.mito_env$term2gene)
  rows <- list()
  for (gid in unique(mito$go_ids)) {
    rows <- append_term2gene(rows, gid, entrez_for_go(gid))
  }
  if (length(rows) == 0 && has_pkg("GO.db") && length(mito$names) > 0) {
    terms <- tryCatch(AnnotationDbi::Term(GO.db::GOTERM), error = function(e) NULL)
    if (!is.null(terms)) {
      for (nm in mito$names) {
        hit <- names(terms)[tolower(terms) == tolower(nm)]
        if (length(hit) == 0) hit <- names(terms)[grepl(tolower(nm), tolower(terms), fixed = TRUE)]
        for (gid in utils::head(hit, 5)) {
          rows <- append_term2gene(rows, paste(gid, terms[[gid]]), entrez_for_go(gid))
        }
      }
    }
  }
  if (length(rows) == 0) {
    log_msg("WARNING: mitochondria TERM2GENE empty; focused enricher will skip")
    .mito_env$term2gene <- data.frame(gs_name = character(), entrez = character())
    return(.mito_env$term2gene)
  }
  t2g <- unique(do.call(rbind, rows))
  log_msg("mitochondria gene sets: ", length(unique(t2g$gs_name)),
          " terms, ", length(unique(t2g$entrez)), " Entrez genes")
  .mito_env$term2gene <- t2g
  t2g
}

plot_focus_term_bar <- function(df, stub, title) {
  if (nrow(df) == 0) return(invisible(NULL))
  lab <- if ("Description" %in% names(df)) as.character(df$Description) else as.character(df$ID)
  df$lab <- lab
  if ("NES" %in% names(df) && any(is.finite(df$NES))) {
    df <- df[order(abs(df$NES), decreasing = TRUE), , drop = FALSE]
    df <- utils::head(df, 20)
    df$lab <- factor(df$lab, levels = rev(unique(df$lab)))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = NES, y = lab)) +
      ggplot2::geom_col(fill = "#2A9D8F") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::labs(title = title, y = NULL)
  } else {
    yv <- if ("p.adjust" %in% names(df)) df$p.adjust else df$pvalue
    df$neglog <- -log10(pmax(as.numeric(yv), 1e-300))
    df <- df[order(yv), , drop = FALSE]
    df <- utils::head(df, 20)
    df$lab <- factor(df$lab, levels = rev(unique(df$lab)))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = neglog, y = lab)) +
      ggplot2::geom_col(fill = "#2A9D8F") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::labs(title = title, x = "-log10(p.adjust)", y = NULL)
  }
  save_gg(p, stub, width = 12, height = max(5, min(12, 0.38 * nrow(df) + 2)))
}

export_mito_focus_terms <- function(x, stub, title, mito) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) return(invisible(NULL))
  df <- as.data.frame(x)
  desc <- if ("Description" %in% names(df)) df$Description else df$ID
  df$genome_wide_rank <- seq_len(nrow(df))
  keep <- vapply(seq_len(nrow(df)), function(i) {
    term_matches_mito(df$ID[i], desc[i], mito)
  }, logical(1))
  hit <- df[keep, , drop = FALSE]
  if (nrow(hit) == 0) {
    note_empty(paste0(stub, "_FOCUS_mitochondria"), "no mitochondria-file terms in this genome-wide result")
    return(invisible(NULL))
  }
  utils::write.csv(hit, paste0(stub, "_FOCUS_mitochondria.csv"), row.names = FALSE)
  plot_focus_term_bar(
    hit, paste0(stub, "_FOCUS_mitochondria_barplot"),
    paste(title, "| mitochondria terms (original p-values)")
  )
}

# -----------------------------------------------------------------------------
# 8. 富集分析（ORA 与 GSEA 分开目录）
# -----------------------------------------------------------------------------
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
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  )
  gene_col <- intersect(c("ncbi_gene", "entrez_gene"), names(msig))[1]
  msig[, c("gs_name", gene_col)]
}

plot_ora_object <- function(x, stub, title, fold_change = NULL, mito = NULL) {
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
  x2 <- tryCatch(enrichplot::pairwise_termsim(x), error = function(e) NULL)
  if (!is.null(x2) && nrow(df) >= 2) {
    try_save_plot(function() enrichplot::emapplot(x2, showCategory = min(20, nrow(df))) + ggplot2::ggtitle(title),
                  paste0(stub, "_emapplot"), 10, 8)
    try_save_plot(function() enrichplot::treeplot(x2, showCategory = min(20, nrow(df))) + ggplot2::ggtitle(title),
                  paste0(stub, "_treeplot"), 11, 8)
  }
  try_save_plot(function() enrichplot::cnetplot(
    x, showCategory = min(8, nshow), foldChange = fold_change, circular = FALSE
  ) + ggplot2::ggtitle(title), paste0(stub, "_cnetplot"), 10, 8)
  try_save_plot(function() enrichplot::heatplot(x, showCategory = nshow, foldChange = fold_change) +
                  ggplot2::ggtitle(title), paste0(stub, "_heatplot"), 11, 6)
  if (!is.null(mito)) export_mito_focus_terms(x, stub, title, mito)
}

plot_gsea_object <- function(x, stub, title, mito = NULL) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) {
    note_empty(stub, "no GSEA terms")
    return(invisible(NULL))
  }
  df <- as.data.frame(x)
  utils::write.csv(df, paste0(stub, ".csv"), row.names = FALSE)
  nshow <- min(15, nrow(df))
  try_save_plot(function() {
    p <- enrichplot::dotplot(x, showCategory = nshow, split = ".sign")
    p <- tryCatch(p + ggplot2::facet_grid(. ~ .sign) + ggplot2::ggtitle(title),
                  error = function(e) p + ggplot2::ggtitle(title))
    p
  }, paste0(stub, "_dotplot"), 10, 7)
  try_save_plot(function() enrichplot::ridgeplot(x, showCategory = nshow) + ggplot2::ggtitle(title),
                paste0(stub, "_ridgeplot"), 10, 8)
  ncurve <- min(5, nrow(df))
  try_save_plot(function() enrichplot::gseaplot2(x, geneSetID = seq_len(ncurve), pvalue_table = TRUE, title = title),
                paste0(stub, "_gseaplot"), 10, 8)
  for (i in seq_len(min(3, nrow(df)))) {
    desc <- as.character(df$Description[i])
    try_save_plot(function() enrichplot::gseaplot2(x, geneSetID = i, title = desc),
                  paste0(stub, "_gseaplot_top", i), 8, 6)
  }
  x2 <- tryCatch(enrichplot::pairwise_termsim(x), error = function(e) NULL)
  if (!is.null(x2) && nrow(df) >= 2) {
    try_save_plot(function() enrichplot::emapplot(x2, showCategory = min(20, nrow(df))) + ggplot2::ggtitle(title),
                  paste0(stub, "_emapplot"), 10, 8)
    try_save_plot(function() enrichplot::cnetplot(x, showCategory = min(8, nshow)) + ggplot2::ggtitle(title),
                  paste0(stub, "_cnetplot"), 10, 8)
  }
  if (!is.null(mito)) export_mito_focus_terms(x, stub, title, mito)
}

plot_kegg_pathview <- function(kegg_obj, stats, outdir) {
  if (is.null(kegg_obj) || nrow(as.data.frame(kegg_obj)) == 0) return(invisible(NULL))
  if (!has_pkg("pathview")) {
    log_msg("pathview not installed; skip KEGG pathway maps")
    return(invisible(NULL))
  }
  ids <- utils::head(as.character(as.data.frame(kegg_obj)$ID), 3)
  old <- getwd()
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  setwd(outdir)
  on.exit(setwd(old), add = TRUE)
  for (id in ids) {
    pid <- sub("^hsa", "", id)
    tryCatch(
      pathview::pathview(gene.data = stats, pathway.id = pid, species = "hsa", kegg.native = TRUE),
      error = function(e) log_msg("pathview failed for ", id, ": ", e$message)
    )
  }
}

run_focused_mito_ora <- function(entrez, outdir, label, tag, fc_sym, mito) {
  t2g <- get_mito_term2gene(mito)
  fdir <- file.path(outdir, "Focused_mitochondria")
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  pref <- paste0(tag, "_")
  writeLines(
    c("这不是改全库 GO/KEGG 的 p 值或排名。",
      "全库结果旁边的 *_FOCUS_mitochondria.csv 保留原始 p 值和 genome_wide_rank。",
      "本文件夹只检验 mitochondria 文本中的通路。"),
    file.path(fdir, "00_README.txt")
  )
  if (is.null(t2g) || nrow(t2g) < 5) {
    note_empty(file.path(fdir, paste0(pref, "ORA_focused_mitochondria")), "empty mitochondria gene sets")
    return(invisible(NULL))
  }
  obj <- enrich_or_relax(
    function() clusterProfiler::enricher(
      entrez, TERM2GENE = t2g, minGSSize = 5, maxGSSize = 2500,
      pvalueCutoff = 0.05, qvalueCutoff = 0.2
    ),
    function() clusterProfiler::enricher(
      entrez, TERM2GENE = t2g, minGSSize = 5, maxGSSize = 2500,
      pvalueCutoff = 1, qvalueCutoff = 1
    ),
    "focused ORA mitochondria"
  )
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    obj <- tryCatch(
      clusterProfiler::setReadable(obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
      error = function(e) obj
    )
  }
  plot_ora_object(
    obj, file.path(fdir, paste0(pref, "ORA_focused_mitochondria")),
    title_maybe_relaxed(obj, paste(label, "| ORA focused mitochondria file")),
    fold_change = fc_sym, mito = NULL
  )
}

run_ora_plots <- function(genes, de_sub, outdir, label, tag, mito) {
  go_dir <- file.path(outdir, "GO")
  pw_dir <- file.path(outdir, "Pathway")
  kg_dir <- file.path(outdir, "KEGG")
  dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(kg_dir, recursive = TRUE, showWarnings = FALSE)
  pref <- paste0(tag, "_")

  mp <- map_to_entrez(genes)
  entrez <- unique(mp$entrez)
  fc_sym <- setNames(de_sub$log2FC, de_sub$gene)
  fc_entrez <- setNames(de_sub$log2FC[match(mp$gene, de_sub$gene)], mp$entrez)
  if (length(entrez) < 3) {
    log_msg("ORA skipped, mapped genes < 3: ", outdir)
    writeLines(paste("mapped_entrez", length(entrez)), file.path(outdir, paste0(pref, "ORA_skipped.txt")))
    note_empty(file.path(go_dir, paste0(pref, "ORA_GO")), "too few mapped genes")
    note_empty(file.path(pw_dir, paste0(pref, "ORA_Pathway")), "too few mapped genes")
    note_empty(file.path(kg_dir, paste0(pref, "ORA_KEGG")), "too few mapped genes")
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
    plot_ora_object(ego, file.path(go_dir, paste0(pref, "ORA_GO_", ont)),
                    title_maybe_relaxed(ego, paste(label, "| ORA GO", ont, "(not GSEA)")),
                    fold_change = fc_sym, mito = mito)
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
    ek <- tryCatch(clusterProfiler::setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"), error = function(e) ek)
  }
  plot_ora_object(ek, file.path(kg_dir, paste0(pref, "ORA_KEGG")),
                  title_maybe_relaxed(ek, paste(label, "| ORA KEGG (not GSEA)")),
                  fold_change = fc_sym, mito = mito)
  plot_kegg_pathview(ek, fc_entrez, kg_dir)

  if (has_pkg("ReactomePA")) {
    er <- enrich_or_relax(
      function() ReactomePA::enrichPathway(
        gene = entrez, organism = "human", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
      ),
      function() ReactomePA::enrichPathway(
        gene = entrez, organism = "human", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
      ),
      "enrichPathway"
    )
    plot_ora_object(er, file.path(pw_dir, paste0(pref, "ORA_Reactome_pathway")),
                    title_maybe_relaxed(er, paste(label, "| ORA Reactome pathway (not GSEA)")),
                    fold_change = fc_sym, mito = mito)
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
    hm <- tryCatch(clusterProfiler::setReadable(hm, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"), error = function(e) hm)
  }
  plot_ora_object(hm, file.path(pw_dir, paste0(pref, "ORA_MSigDB_Hallmark_pathway")),
                  title_maybe_relaxed(hm, paste(label, "| ORA Hallmark pathway (not GSEA)")),
                  fold_change = fc_sym, mito = mito)
  tryCatch(
    run_focused_mito_ora(entrez, outdir, label, tag, fc_sym, mito),
    error = function(e) log_msg("focused mito ORA failed: ", e$message)
  )
  writeLines(
    c("This GO/Pathway/KEGG folder is ORA (over-representation), NOT GSEA.",
      "GSEA files are in the sibling folder named GSEA/ and start with GSEA_.",
      "User-facing 'GWAS' in the request is implemented as GSEA."),
    file.path(outdir, paste0(pref, "00_ORA_is_not_GSEA.txt"))
  )
}

plot_fgsea_hallmark <- function(stats, outdir, title, prefix = "") {
  stub <- function(x) file.path(outdir, paste0(prefix, x))
  if (length(stats) < 5) {
    note_empty(stub("GSEA_Hallmark_fgsea"), "too few ranked genes")
    return(invisible(NULL))
  }
  term2gene <- tryCatch(msig_hallmark_map(), error = function(e) NULL)
  if (is.null(term2gene)) return(invisible(NULL))
  pathways <- split(as.character(term2gene[[2]]), term2gene[[1]])
  fg <- tryCatch(fgsea::fgsea(pathways = pathways, stats = stats, minSize = 5, maxSize = 500), error = function(e) {
    log_msg("fgsea Hallmark failed: ", e$message)
    NULL
  })
  if (is.null(fg) || nrow(fg) == 0) {
    note_empty(stub("GSEA_Hallmark_fgsea"), "no fgsea terms")
    return(invisible(NULL))
  }
  fg <- as.data.frame(fg)
  fg <- fg[order(fg$pval), ]
  utils::write.csv(fg, stub("GSEA_Hallmark_fgsea.csv"), row.names = FALSE)
  plot_df <- utils::head(fg, 15)
  plot_df$pathway <- factor(plot_df$pathway, levels = rev(plot_df$pathway))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = NES, y = pathway, fill = padj < 0.05)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#D62828", "FALSE" = "grey70")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title, y = NULL, fill = "padj < 0.05")
  save_gg(p, stub("GSEA_Hallmark_fgsea_barplot"), 10, 7)
}

gsea_ids_overlapping_genes <- function(gsea_obj, symbols, entrez) {
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) return(character())
  df <- as.data.frame(gsea_obj)
  if (!"core_enrichment" %in% names(df)) return(character())
  keep <- vapply(df$core_enrichment, function(s) {
    gs <- unlist(strsplit(as.character(s), "/"))
    any(gs %in% symbols) || any(gs %in% entrez)
  }, logical(1))
  df$ID[keep]
}

plot_gsea_selected_ids <- function(x, ids, stub, title) {
  if (is.null(x) || length(ids) == 0) {
    note_empty(stub, "no overlapping GSEA terms")
    return(invisible(NULL))
  }
  ids <- ids[ids %in% as.data.frame(x)$ID]
  if (length(ids) == 0) {
    note_empty(stub, "no overlapping GSEA terms")
    return(invisible(NULL))
  }
  ids <- utils::head(ids, 5)
  try_save_plot(function() enrichplot::gseaplot2(x, geneSetID = ids, pvalue_table = TRUE, title = title),
                stub, 10, 8)
}

build_gsea_cache <- function(de) {
  stats <- ranked_entrez(de)
  out <- list(stats = stats)
  if (length(stats) < 10) return(out)
  gsea_one <- function(fun, label) {
    enrich_or_relax(
      function() fun(pvalueCutoff = 0.05, minGSSize = 10),
      function() fun(pvalueCutoff = 1, minGSSize = 5),
      label
    )
  }
  out$GO_BP <- gsea_one(function(pvalueCutoff, minGSSize) {
    clusterProfiler::gseGO(
      geneList = stats, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID",
      minGSSize = minGSSize, maxGSSize = 500, pvalueCutoff = pvalueCutoff,
      verbose = FALSE, eps = 0
    )
  }, "gseGO_BP")
  out$GO_MF <- gsea_one(function(pvalueCutoff, minGSSize) {
    clusterProfiler::gseGO(
      geneList = stats, OrgDb = org.Hs.eg.db, ont = "MF", keyType = "ENTREZID",
      minGSSize = minGSSize, maxGSSize = 500, pvalueCutoff = pvalueCutoff,
      verbose = FALSE, eps = 0
    )
  }, "gseGO_MF")
  out$GO_CC <- gsea_one(function(pvalueCutoff, minGSSize) {
    clusterProfiler::gseGO(
      geneList = stats, OrgDb = org.Hs.eg.db, ont = "CC", keyType = "ENTREZID",
      minGSSize = minGSSize, maxGSSize = 500, pvalueCutoff = pvalueCutoff,
      verbose = FALSE, eps = 0
    )
  }, "gseGO_CC")
  out$KEGG <- gsea_one(function(pvalueCutoff, minGSSize) {
    clusterProfiler::gseKEGG(
      geneList = stats, organism = "hsa", minGSSize = minGSSize, maxGSSize = 500,
      pvalueCutoff = pvalueCutoff, verbose = FALSE, eps = 0
    )
  }, "gseKEGG")
  if (has_pkg("ReactomePA")) {
    out$Reactome <- gsea_one(function(pvalueCutoff, minGSSize) {
      ReactomePA::gsePathway(
        stats, organism = "human", minGSSize = minGSSize, maxGSSize = 500,
        pvalueCutoff = pvalueCutoff, verbose = FALSE, eps = 0
      )
    }, "gsePathway")
  }
  term2gene <- tryCatch(msig_hallmark_map(), error = function(e) NULL)
  if (!is.null(term2gene)) {
    out$Hallmark <- gsea_one(function(pvalueCutoff, minGSSize) {
      clusterProfiler::GSEA(
        geneList = stats, TERM2GENE = term2gene, minGSSize = minGSSize,
        maxGSSize = 500, pvalueCutoff = pvalueCutoff, eps = 0, verbose = FALSE
      )
    }, "GSEA_Hallmark")
  }
  for (nm in c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Hallmark")) {
    if (!is.null(out[[nm]]) && nrow(as.data.frame(out[[nm]])) > 0) {
      out[[nm]] <- tryCatch(
        clusterProfiler::setReadable(out[[nm]], OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
        error = function(e) out[[nm]]
      )
    }
  }
  out
}

run_gsea_plots <- function(sub, gsea_cache, outdir, tag, label, mito) {
  gsea_dir <- file.path(outdir, "GSEA")
  dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
  pref <- paste0(tag, "_")
  sub_stats <- ranked_entrez(sub)
  plot_fgsea_hallmark(sub_stats, gsea_dir, paste(label, "| GSEA Hallmark"), prefix = pref)

  if (length(sub_stats) >= 8) {
    term2gene <- tryCatch(msig_hallmark_map(), error = function(e) NULL)
    if (!is.null(term2gene)) {
      hm <- enrich_or_relax(
        function() clusterProfiler::GSEA(
          geneList = sub_stats, TERM2GENE = term2gene, minGSSize = 5,
          maxGSSize = 500, pvalueCutoff = 0.05, eps = 0, verbose = FALSE
        ),
        function() clusterProfiler::GSEA(
          geneList = sub_stats, TERM2GENE = term2gene, minGSSize = 3,
          maxGSSize = 500, pvalueCutoff = 1, eps = 0, verbose = FALSE
        ),
        paste("subset Hallmark GSEA", tag)
      )
      plot_gsea_object(hm, file.path(gsea_dir, paste0(pref, "GSEA_Hallmark")),
                       paste(label, "| GSEA Hallmark (subset ranked)"), mito = mito)
    }
    kegg <- enrich_or_relax(
      function() clusterProfiler::gseKEGG(
        geneList = sub_stats, organism = "hsa", minGSSize = 5, maxGSSize = 500,
        pvalueCutoff = 0.05, verbose = FALSE, eps = 0
      ),
      function() clusterProfiler::gseKEGG(
        geneList = sub_stats, organism = "hsa", minGSSize = 3, maxGSSize = 500,
        pvalueCutoff = 1, verbose = FALSE, eps = 0
      ),
      paste("subset KEGG GSEA", tag)
    )
    if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
      kegg <- tryCatch(clusterProfiler::setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"), error = function(e) kegg)
    }
    plot_gsea_object(kegg, file.path(gsea_dir, paste0(pref, "GSEA_KEGG")),
                     paste(label, "| GSEA KEGG (subset ranked)"), mito = mito)
  }

  mp <- map_to_entrez(sub$gene)
  for (nm in c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Hallmark")) {
    ids <- gsea_ids_overlapping_genes(gsea_cache[[nm]], sub$gene, mp$entrez)
    plot_gsea_selected_ids(
      gsea_cache[[nm]], ids,
      file.path(gsea_dir, paste0(pref, "GSEA_fullrank_overlap_", nm)),
      paste(label, "| GSEA", nm, "(full-rank overlap)")
    )
  }
}

run_focused_mito_gsea <- function(stats, outdir, label, mito) {
  t2g <- get_mito_term2gene(mito)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c("这不是改全库 GO/KEGG 的 p 值或排名。",
      "本文件夹只检验 mitochondria 文本中的通路。",
      "全库结果旁边的 *_FOCUS_mitochondria.csv 保留原始 p 值和 genome_wide_rank。"),
    file.path(outdir, "00_README.txt")
  )
  if (is.null(t2g) || nrow(t2g) < 5 || length(stats) < 10) {
    note_empty(file.path(outdir, "GSEA_focused_mitochondria"), "too few genes or empty gene sets")
    return(invisible(NULL))
  }
  obj <- enrich_or_relax(
    function() clusterProfiler::GSEA(
      geneList = stats, TERM2GENE = t2g, minGSSize = 5, maxGSSize = 2500,
      pvalueCutoff = 0.05, eps = 0, verbose = FALSE
    ),
    function() clusterProfiler::GSEA(
      geneList = stats, TERM2GENE = t2g, minGSSize = 5, maxGSSize = 2500,
      pvalueCutoff = 1, eps = 0, verbose = FALSE
    ),
    "focused GSEA mitochondria"
  )
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    obj <- tryCatch(
      clusterProfiler::setReadable(obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
      error = function(e) obj
    )
  }
  plot_gsea_object(
    obj, file.path(outdir, "GSEA_focused_mitochondria"),
    paste(label, "| GSEA focused mitochondria file"),
    mito = NULL
  )
}

emit_subset_analysis <- function(class_name, sub, tag, title, outdir, full_de_for_volcano,
                                 heat_mat, sample_info, gsea_cache, mito, fc_line = 1,
                                 heat_samples = NULL, stable_de = NULL) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(paste("class:", class_name),
      paste("subset:", tag),
      paste("title:", title),
      paste("n_proteins:", nrow(sub)),
      "p-values were not estimated for 1-vs-1 samples."),
    file.path(outdir, paste0("00_", tag, "_THIS_FOLDER.txt"))
  )
  utils::write.csv(sub, file.path(outdir, paste0(tag, "_DE_selected_proteins.csv")), row.names = FALSE)
  if (has_pkg("writexl")) {
    tryCatch(writexl::write_xlsx(sub, file.path(outdir, paste0(tag, "_DE_selected_proteins.xlsx"))),
             error = function(e) log_msg("xlsx write failed: ", e$message))
  }
  log_msg(class_name, " ", tag, ": n = ", nrow(sub))
  if (nrow(sub) == 0) {
    writeLines("no proteins", file.path(outdir, paste0(tag, "_EMPTY.txt")))
    return(invisible(NULL))
  }
  tryCatch(plot_de_bar(sub, paste0(title, " | selected proteins"), file.path(outdir, paste0(tag, "_DE_log2FC_barplot"))),
           error = function(e) log_msg("DE barplot failed: ", e$message))
  tryCatch(plot_volcano(full_de_for_volcano, sub$gene, title, file.path(outdir, paste0(tag, "_volcano")), fc_line = fc_line),
           error = function(e) log_msg("volcano failed: ", e$message))
  if (!is.null(stable_de)) {
    tryCatch(
      plot_contrast_scatter(full_de_for_volcano, stable_de, sub$gene, title,
                            file.path(outdir, paste0(tag, "_scatter_change_vs_stable"))),
      error = function(e) log_msg("scatter failed: ", e$message)
    )
  }
  if (all(c("log2FC_A_or_C", "log2FC_B_or_D") %in% names(sub))) {
    tryCatch({
      df <- sub
      df$set <- "selected"
      df$label <- ifelse(df$gene %in% utils::head(df$gene, 15), df$gene, NA)
      p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC_A_or_C, y = log2FC_B_or_D, color = set)) +
        ggplot2::geom_point(alpha = 0.8, size = 1.8) +
        ggplot2::geom_vline(xintercept = 0, linetype = 3, color = "grey60") +
        ggplot2::geom_hline(yintercept = 0, linetype = 3, color = "grey60") +
        ggplot2::scale_color_manual(values = c(selected = "#D62828")) +
        ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::labs(title = title, x = "log2FC class A or C", y = "log2FC class B or D", color = NULL)
      save_gg(p, file.path(outdir, paste0(tag, "_scatter_intersect_contrasts")))
    }, error = function(e) log_msg("intersect scatter failed: ", e$message))
  }
  tryCatch(plot_heatmap(heat_mat, sample_info, sub$gene, title,
                        file.path(outdir, paste0(tag, "_heatmap")), samples = heat_samples),
           error = function(e) {
             while (grDevices::dev.cur() > 1) grDevices::dev.off()
             log_msg("heatmap failed: ", e$message)
           })
  tryCatch(run_ora_plots(sub$gene, sub, outdir, title, tag, mito),
           error = function(e) log_msg("ORA/plots failed: ", e$message))
  tryCatch(run_gsea_plots(sub, gsea_cache, outdir, tag, title, mito),
           error = function(e) log_msg("GSEA/plots failed: ", e$message))
}

analyze_one_class <- function(cls, heat_mat, sample_info, mito) {
  base <- file.path(result_dir, cls$dir)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(cls$de, file.path(base, "DE_full_change_contrast.csv"), row.names = FALSE)
  writeLines(
    c(paste("class:", cls$id),
      paste("label:", cls$label),
      paste("change:", cls$change_treat, "vs", cls$change_ctrl, cls$direction),
      paste("stable/no-change:", cls$stable_treat, "vs", cls$stable_ctrl,
            "| 1/", unchanged_max_fc, "< FC <", unchanged_max_fc),
      "1-vs-1: no p-value was estimated."),
    file.path(base, "00_READ_ME_先看这里.txt")
  )
  fc_dirs <- file.path(base, "FoldChange", names(fc_cutoffs))
  top_dirs <- file.path(base, "TopRank", paste0("top", top_ns))
  invisible(lapply(c(fc_dirs, top_dirs), dir.create, recursive = TRUE, showWarnings = FALSE))

  gsea_cache <- list()
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    sub <- cls$select_fc(fc)
    tryCatch(
      emit_subset_analysis(
        cls$id, sub, nm, paste0(cls$label, " | ", nm),
        file.path(base, "FoldChange", nm),
        cls$de, heat_mat, sample_info, gsea_cache, mito,
        fc_line = fc, heat_samples = cls$heat_samples, stable_de = cls$stable_de
      ),
      error = function(e) log_msg("ERROR subset ", cls$id, " ", nm, ": ", e$message)
    )
  }
  for (n in top_ns) {
    tag <- paste0("top", n)
    sub <- cls$select_top(n)
    tryCatch(
      emit_subset_analysis(
        cls$id, sub, tag, paste0(cls$label, " | ", tag),
        file.path(base, "TopRank", tag),
        cls$de, heat_mat, sample_info, gsea_cache, mito,
        fc_line = 1, heat_samples = cls$heat_samples, stable_de = cls$stable_de
      ),
      error = function(e) log_msg("ERROR subset ", cls$id, " ", tag, ": ", e$message)
    )
  }

  tryCatch({
    log_msg("Building full-list GSEA after subset plots: ", cls$id)
    gsea_cache <- build_gsea_cache(cls$gsea_de)
    full_gsea_dir <- file.path(base, "00_GSEA_all_genes_NOT_FC_or_topN")
    dir.create(full_gsea_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines("全蛋白 GSEA，不是 FC/topN 分层结果。分层图在 FoldChange/ 和 TopRank/。用户所称 GWAS 在此按 GSEA 输出。",
               file.path(full_gsea_dir, "00_README.txt"))
    for (nm in c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Hallmark")) {
      plot_gsea_object(gsea_cache[[nm]], file.path(full_gsea_dir, paste0("allGenes_GSEA_", nm)),
                       paste("GSEA", nm, "|", cls$label, "| ALL proteins, NOT FC/topN"), mito = mito)
    }
    plot_fgsea_hallmark(gsea_cache$stats, full_gsea_dir,
                        paste("GSEA Hallmark |", cls$label, "| ALL proteins"), prefix = "allGenes_")
  }, error = function(e) log_msg("full-list GSEA failed for ", cls$id, ": ", e$message))

  tryCatch({
    focus_stats <- if (exists("gsea_cache", inherits = FALSE) && !is.null(gsea_cache$stats)) {
      gsea_cache$stats
    } else {
      ranked_entrez(cls$gsea_de)
    }
    run_focused_mito_gsea(
      focus_stats,
      file.path(base, "Focused_mitochondria"),
      paste(cls$label, "| all proteins"),
      mito
    )
  }, error = function(e) log_msg("focused mito GSEA failed for ", cls$id, ": ", e$message))
}

# -----------------------------------------------------------------------------
# 9. 六类蛋白定义
# -----------------------------------------------------------------------------
annotate_de <- function(de, meta) {
  if (is.null(de) || nrow(de) == 0) return(de)
  hit <- match(de$gene, meta$gene)
  de$symbol <- meta$symbol[hit]
  de$uniprot <- meta$uniprot[hit]
  de$protein_id <- meta$protein_id[hit]
  de
}

select_up_stable <- function(change_de, stable_de, fc) {
  keep <- is.finite(change_de$log2FC) & (2^change_de$log2FC >= fc)
  stab <- is_unchanged(stable_de$log2FC[match(change_de$gene, stable_de$gene)])
  out <- change_de[keep & stab, , drop = FALSE]
  if (nrow(out) > 0) out <- out[order(out$log2FC, decreasing = TRUE), , drop = FALSE]
  out
}

select_down_stable <- function(change_de, stable_de, fc) {
  keep <- is.finite(change_de$log2FC) & (2^change_de$log2FC <= (1 / fc))
  stab <- is_unchanged(stable_de$log2FC[match(change_de$gene, stable_de$gene)])
  out <- change_de[keep & stab, , drop = FALSE]
  if (nrow(out) > 0) out <- out[order(out$log2FC, decreasing = FALSE), , drop = FALSE]
  out
}

select_up_stable_top <- function(change_de, stable_de, n) {
  pool <- select_up_stable(change_de, stable_de, fc = 1)
  utils::head(pool, n)
}

select_down_stable_top <- function(change_de, stable_de, n) {
  pool <- select_down_stable(change_de, stable_de, fc = 1)
  utils::head(pool, n)
}

intersect_tables <- function(a, b, log2fc_fun) {
  genes <- intersect(a$gene, b$gene)
  if (length(genes) == 0) {
    return(a[0, , drop = FALSE])
  }
  aa <- a[match(genes, a$gene), ]
  bb <- b[match(genes, b$gene), ]
  out <- aa
  out$log2FC_A_or_C <- aa$log2FC
  out$log2FC_B_or_D <- bb$log2FC
  out$log2FC <- log2fc_fun(aa$log2FC, bb$log2FC)
  out$FoldChange <- 2^out$log2FC
  out$AveExpr <- (aa$AveExpr + bb$AveExpr) / 2
  out[order(abs(out$log2FC), decreasing = TRUE), , drop = FALSE]
}

make_class <- function(id, dir, label, direction, change_de, stable_de,
                       change_treat, change_ctrl, stable_treat, stable_ctrl,
                       heat_samples, gsea_de, select_fc, select_top) {
  list(
    id = id, dir = dir, label = label, direction = direction,
    de = change_de, stable_de = stable_de, gsea_de = gsea_de,
    change_treat = change_treat, change_ctrl = change_ctrl,
    stable_treat = stable_treat, stable_ctrl = stable_ctrl,
    heat_samples = as.character(heat_samples),
    select_fc = select_fc, select_top = select_top
  )
}

save_two_set_venn <- function(lst, outfile, title, fill_color) {
  if (!has_pkg("ggvenn")) return(invisible(NULL))
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  p <- tryCatch(
    ggvenn::ggvenn(lst, fill_color = fill_color) + ggplot2::labs(title = title),
    error = function(e) {
      log_msg("venn failed: ", e$message)
      NULL
    }
  )
  if (!is.null(p)) save_gg(p, outfile, width = 7, height = 6)
}

# -----------------------------------------------------------------------------
# 10. 主流程
# -----------------------------------------------------------------------------
log_msg("Project dir: ", project_dir)
log_msg("Unchanged definition: 1/", unchanged_max_fc, " < FC < ", unchanged_max_fc)
log_msg("No p-values will be fabricated for 1-vs-1 sample pairs.")

obj <- load_protein_matrix(project_dir)
mito <- parse_mitochondria_file(project_dir)
log_msg("Loaded from ", obj$source, " | proteins=", nrow(obj$mat), " samples=", ncol(obj$mat))
utils::write.csv(obj$sample_info, file.path(log_dir, "sample_info.csv"), row.names = FALSE)

need <- as.character(1:16)
if (!all(need %in% colnames(obj$mat))) {
  stop("需要样品 1–16 全部存在。当前有: ", paste(colnames(obj$mat), collapse = ", "))
}

filt <- filter_low_abundance(obj$mat)
meta <- obj$meta[match(rownames(filt), obj$meta$gene), ]
norm <- normalize_intensity(filt)
log_mat <- norm$log_mat
si <- obj$sample_info[match(colnames(log_mat), obj$sample_info$sample), ]

utils::write.csv(
  cbind(gene = rownames(log_mat), as.data.frame(log_mat)),
  file.path(result_dir, "normalized_log_matrix.csv"),
  row.names = FALSE
)
plot_pca(norm$heat_mat, si, file.path(result_dir, "00_QC_PCA"))

fc_2v1  <- annotate_de(pairwise_fc(log_mat, 2, 1, "2_vs_1"), meta)
fc_3v4  <- annotate_de(pairwise_fc(log_mat, 3, 4, "3_vs_4"), meta)
fc_5v6  <- annotate_de(pairwise_fc(log_mat, 5, 6, "5_vs_6"), meta)
fc_7v8  <- annotate_de(pairwise_fc(log_mat, 7, 8, "7_vs_8"), meta)
fc_9v10 <- annotate_de(pairwise_fc(log_mat, 9, 10, "9_vs_10"), meta)
fc_11v12 <- annotate_de(pairwise_fc(log_mat, 11, 12, "11_vs_12"), meta)
fc_13v14 <- annotate_de(pairwise_fc(log_mat, 13, 14, "13_vs_14"), meta)
fc_15v16 <- annotate_de(pairwise_fc(log_mat, 15, 16, "15_vs_16"), meta)

classes <- list(
  make_class(
    "A", "classA_2vs1_up_3vs4_unchanged",
    "Class A: 2 vs 1 up AND 3 vs 4 unchanged", "up",
    fc_2v1, fc_3v4, 2, 1, 3, 4, c(1, 2, 3, 4), fc_2v1,
    function(fc) select_up_stable(fc_2v1, fc_3v4, fc),
    function(n) select_up_stable_top(fc_2v1, fc_3v4, n)
  ),
  make_class(
    "B", "classB_5vs6_down_7vs8_unchanged",
    "Class B: 5 vs 6 down AND 7 vs 8 unchanged", "down",
    fc_5v6, fc_7v8, 5, 6, 7, 8, c(5, 6, 7, 8), fc_5v6,
    function(fc) select_down_stable(fc_5v6, fc_7v8, fc),
    function(n) select_down_stable_top(fc_5v6, fc_7v8, n)
  ),
  make_class(
    "C", "classC_9vs10_up_11vs12_unchanged",
    "Class C: 9 vs 10 up AND 11 vs 12 unchanged", "up",
    fc_9v10, fc_11v12, 9, 10, 11, 12, c(9, 10, 11, 12), fc_9v10,
    function(fc) select_up_stable(fc_9v10, fc_11v12, fc),
    function(n) select_up_stable_top(fc_9v10, fc_11v12, n)
  ),
  make_class(
    "D", "classD_13vs14_down_15vs16_unchanged",
    "Class D: 13 vs 14 down AND 15 vs 16 unchanged", "down",
    fc_13v14, fc_15v16, 13, 14, 15, 16, c(13, 14, 15, 16), fc_13v14,
    function(fc) select_down_stable(fc_13v14, fc_15v16, fc),
    function(n) select_down_stable_top(fc_13v14, fc_15v16, n)
  )
)

# E = A ∩ B；GSEA 排序用「A 上调 + B 下调」的联合效应
gsea_E <- fc_2v1
gsea_E$log2FC <- fc_2v1$log2FC - fc_5v6$log2FC[match(fc_2v1$gene, fc_5v6$gene)]
gsea_E$FoldChange <- 2^gsea_E$log2FC
gsea_E$AveExpr <- (fc_2v1$AveExpr + fc_5v6$AveExpr[match(fc_2v1$gene, fc_5v6$gene)]) / 2

gsea_F <- fc_9v10
gsea_F$log2FC <- fc_9v10$log2FC - fc_13v14$log2FC[match(fc_9v10$gene, fc_13v14$gene)]
gsea_F$FoldChange <- 2^gsea_F$log2FC
gsea_F$AveExpr <- (fc_9v10$AveExpr + fc_13v14$AveExpr[match(fc_9v10$gene, fc_13v14$gene)]) / 2

classes[[5]] <- make_class(
  "E", "classE_A_intersect_B",
  "Class E: Class A ∩ Class B", "intersect",
  gsea_E, NULL, "A", "B", NA, NA, c(1, 2, 3, 4, 5, 6, 7, 8), gsea_E,
  function(fc) intersect_tables(
    select_up_stable(fc_2v1, fc_3v4, fc),
    select_down_stable(fc_5v6, fc_7v8, fc),
    function(a, b) (a - b) / 2
  ),
  function(n) intersect_tables(
    select_up_stable_top(fc_2v1, fc_3v4, n),
    select_down_stable_top(fc_5v6, fc_7v8, n),
    function(a, b) (a - b) / 2
  )
)
classes[[6]] <- make_class(
  "F", "classF_C_intersect_D",
  "Class F: Class C ∩ Class D", "intersect",
  gsea_F, NULL, "C", "D", NA, NA, c(9, 10, 11, 12, 13, 14, 15, 16), gsea_F,
  function(fc) intersect_tables(
    select_up_stable(fc_9v10, fc_11v12, fc),
    select_down_stable(fc_13v14, fc_15v16, fc),
    function(a, b) (a - b) / 2
  ),
  function(n) intersect_tables(
    select_up_stable_top(fc_9v10, fc_11v12, n),
    select_down_stable_top(fc_13v14, fc_15v16, n),
    function(a, b) (a - b) / 2
  )
)

for (cls in classes) {
  tryCatch(
    analyze_one_class(cls, norm$heat_mat, si, mito),
    error = function(e) log_msg("ERROR in class ", cls$id, ": ", e$message)
  )
}

tryCatch({
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    save_two_set_venn(
      list(A = select_up_stable(fc_2v1, fc_3v4, fc)$gene,
           B = select_down_stable(fc_5v6, fc_7v8, fc)$gene),
      file.path(result_dir, "classE_A_intersect_B", "FoldChange", nm, paste0("venn_", nm)),
      paste("Class E = A ∩ B |", nm), c("#F58518", "#54A24B")
    )
    save_two_set_venn(
      list(C = select_up_stable(fc_9v10, fc_11v12, fc)$gene,
           D = select_down_stable(fc_13v14, fc_15v16, fc)$gene),
      file.path(result_dir, "classF_C_intersect_D", "FoldChange", nm, paste0("venn_", nm)),
      paste("Class F = C ∩ D |", nm), c("#4C78A8", "#E45756")
    )
  }
  for (n in top_ns) {
    tag <- paste0("top", n)
    save_two_set_venn(
      list(A = select_up_stable_top(fc_2v1, fc_3v4, n)$gene,
           B = select_down_stable_top(fc_5v6, fc_7v8, n)$gene),
      file.path(result_dir, "classE_A_intersect_B", "TopRank", tag, paste0("venn_", tag)),
      paste("Class E = A ∩ B |", tag), c("#F58518", "#54A24B")
    )
    save_two_set_venn(
      list(C = select_up_stable_top(fc_9v10, fc_11v12, n)$gene,
           D = select_down_stable_top(fc_13v14, fc_15v16, n)$gene),
      file.path(result_dir, "classF_C_intersect_D", "TopRank", tag, paste0("venn_", tag)),
      paste("Class F = C ∩ D |", tag), c("#4C78A8", "#E45756")
    )
  }
}, error = function(e) log_msg("venn plots error: ", e$message))

base::writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
log_msg("All done. Results in: ", result_dir)
