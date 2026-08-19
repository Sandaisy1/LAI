#!/usr/bin/env Rscript
# =============================================================================
# TG BRCA 细胞 RNA-seq 分析流程
# 组别：NTC（对照）、TG_sh1、TG_sh5
# 三种比较：
#   1) TG_sh1 vs NTC ； TG_sh5 vs NTC
#   2) (TG_sh1 + TG_sh5)/2 vs NTC   （两个 knockdown 等权平均）
#   3) 两次单独比较中的共同上调基因
# 预处理：过滤低表达 + 标准化消除技术偏差
# 子集策略：
#   A) 上调 FC >= 1 / 1.25 / 1.5 / 2
#   B) 上调排名 top 50 / 75 / 100 / 150 / 200 / 250 / 300
# 每个比较 × 每个 FC 阈值 × 每个 topN 都必须出图：
#   差异基因表/柱状图、火山图、热图、GO图、通路富集图、KEGG图、GSEA图
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANGUAGE = "en")

# -----------------------------------------------------------------------------
# 0. 依赖包
# -----------------------------------------------------------------------------
cran_pkgs <- c(
  "readxl", "writexl", "dplyr", "tidyr", "tibble", "stringr", "ggplot2",
  "ggrepel", "pheatmap", "RColorBrewer", "matrixStats", "ggvenn",
  "cowplot", "ggridges", "ggnewscale"
)
bioc_pkgs <- c(
  "DESeq2", "edgeR", "limma", "clusterProfiler", "org.Hs.eg.db",
  "enrichplot", "DOSE", "ReactomePA", "AnnotationDbi", "fgsea", "msigdbr",
  "pathview"
)

install_if_missing <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  miss_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_cran) > 0) {
    install.packages(miss_cran, repos = "https://cloud.r-project.org")
  }
  miss_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_bioc) > 0) {
    BiocManager::install(miss_bioc, update = FALSE, ask = FALSE)
  }
}

install_if_missing()
invisible(lapply(c(cran_pkgs, bioc_pkgs), function(p) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}))

# -----------------------------------------------------------------------------
# 1. 路径与分析参数
# -----------------------------------------------------------------------------
resolve_project_dir <- function() {
  env_dir <- Sys.getenv("TG_RNASEQ_DIR", unset = "")
  candidates <- c(
    env_dir,
    "E:/R/TG_BRCA/TG",
    "E:\\R\\TG_BRCA\\TG",
    getwd()
  )
  candidates <- unique(candidates[nzchar(candidates)])
  for (d in candidates) {
    if (dir.exists(d) && (
      file.exists(file.path(d, "shTG(20220723).xlsx")) ||
      file.exists(file.path(d, "genes.read_group_tracking")) ||
      file.exists(file.path(d, "genes.fpkm_tracking")) ||
      file.exists(file.path(d, "genes.count_tracking"))
    )) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
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

padj_cutoff <- 0.05
fc_cutoffs <- c("FC_1" = 1, "FC_1.25" = 1.25, "FC_1.5" = 1.5, "FC_2" = 2)
top_ns     <- c(50, 75, 100, 150, 200, 250, 300)

# -----------------------------------------------------------------------------
# 2. 样本名识别
# -----------------------------------------------------------------------------
classify_sample <- function(name) {
  n <- toupper(gsub("[._\\- ]", "", name))
  if (grepl("SH5|SHRNA5|TGSH5", n)) return("TG_sh5")
  if (grepl("SH1|SHRNA1|TGSH1", n)) return("TG_sh1")
  if (grepl("NTC|SHNC|NEGCTRL|CTRL|CONTROL", n)) return("NTC")
  if (grepl("^NC[0-9]*$", n)) return("NTC")
  NA_character_
}

guess_gene_column <- function(df) {
  nms <- names(df)
  low <- tolower(nms)
  keys <- c("gene_short_name", "gene_name", "gene_symbol", "symbol",
            "geneid", "gene_id", "tracking_id", "id", "genes")
  for (k in keys) {
    hit <- which(low == k | grepl(paste0("^", k, "$"), low))
    if (length(hit) > 0) return(nms[hit[1]])
  }
  non_num <- nms[vapply(df, function(x) !is.numeric(x) && !is.integer(x), logical(1))]
  if (length(non_num) > 0) return(non_num[1])
  nms[1]
}

collapse_by_gene <- function(mat, genes) {
  genes[is.na(genes) | genes == "" | genes == "-"] <- NA
  keep <- !is.na(genes)
  mat <- mat[keep, , drop = FALSE]
  genes <- genes[keep]
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

# -----------------------------------------------------------------------------
# 3. 读入表达矩阵
# -----------------------------------------------------------------------------
read_excel_matrix <- function(xlsx) {
  sheets <- readxl::excel_sheets(xlsx)
  log_msg("Excel sheets: ", paste(sheets, collapse = ", "))
  sheet_use <- sheets[1]
  hit <- grep("matrix|fpkm|count|expr|tpm", sheets, ignore.case = TRUE)
  if (length(hit) > 0) sheet_use <- sheets[hit[1]]
  log_msg("Reading Excel sheet: ", sheet_use)
  df <- as.data.frame(readxl::read_excel(xlsx, sheet = sheet_use), stringsAsFactors = FALSE)
  gene_col <- guess_gene_column(df)
  maybe_num <- setdiff(names(df), gene_col)
  for (cc in maybe_num) {
    if (!is.numeric(df[[cc]])) {
      conv <- suppressWarnings(as.numeric(df[[cc]]))
      if (mean(is.finite(conv)) > 0.8) df[[cc]] <- conv
    }
  }
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  groups <- vapply(num_cols, classify_sample, character(1))
  keep <- !is.na(groups)
  if (sum(keep) < 2) {
    log_msg("Excel numeric columns could not be mapped to NTC/TG_sh1/TG_sh5: ",
            paste(num_cols, collapse = ", "))
    return(NULL)
  }
  mat <- as.matrix(df[, num_cols[keep], drop = FALSE])
  storage.mode(mat) <- "double"
  mat <- collapse_by_gene(mat, as.character(df[[gene_col]]))
  list(
    mat = mat,
    sample_info = data.frame(
      sample = colnames(mat),
      group = unname(groups[keep]),
      stringsAsFactors = FALSE
    ),
    source = basename(xlsx)
  )
}

read_read_group_tracking <- function(path) {
  rg <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  need <- c("tracking_id", "condition", "replicate")
  if (!all(need %in% names(rg))) return(NULL)
  value_col <- if ("raw_frags" %in% names(rg)) {
    "raw_frags"
  } else if ("external_scaled_frags" %in% names(rg)) {
    "external_scaled_frags"
  } else if ("FPKM" %in% names(rg)) {
    "FPKM"
  } else {
    return(NULL)
  }
  rg$sample <- paste(rg$condition, rg$replicate, sep = "_rep")
  rg$group <- vapply(as.character(rg$condition), classify_sample, character(1))
  if (all(is.na(rg$group))) {
    rg$group <- vapply(rg$sample, classify_sample, character(1))
  }
  rg <- rg[!is.na(rg$group), , drop = FALSE]
  if (nrow(rg) == 0) return(NULL)
  wide <- tidyr::pivot_wider(
    rg[, c("tracking_id", "sample", value_col)],
    names_from = "sample",
    values_from = value_col,
    values_fn = mean
  )
  gene_map <- NULL
  fpkm_file <- file.path(dirname(path), "genes.fpkm_tracking")
  if (file.exists(fpkm_file)) {
    fp <- utils::read.delim(fpkm_file, check.names = FALSE, stringsAsFactors = FALSE)
    if (all(c("tracking_id", "gene_short_name") %in% names(fp))) {
      gene_map <- fp[, c("tracking_id", "gene_short_name")]
    }
  }
  genes <- wide$tracking_id
  if (!is.null(gene_map)) {
    genes <- gene_map$gene_short_name[match(wide$tracking_id, gene_map$tracking_id)]
    genes[is.na(genes) | genes == "" | genes == "-"] <- wide$tracking_id[is.na(genes) | genes == "" | genes == "-"]
  }
  mat <- as.matrix(wide[, setdiff(names(wide), "tracking_id"), drop = FALSE])
  storage.mode(mat) <- "double"
  mat[is.na(mat)] <- 0
  mat <- collapse_by_gene(mat, genes)
  sample_info <- unique(rg[, c("sample", "group")])
  sample_info <- sample_info[match(colnames(mat), sample_info$sample), ]
  list(mat = mat, sample_info = sample_info, source = basename(path), value_col = value_col)
}

read_tracking_matrix <- function(path, value_pattern) {
  tr <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  gene_col <- if ("gene_short_name" %in% names(tr)) "gene_short_name" else if ("gene_id" %in% names(tr)) "gene_id" else names(tr)[1]
  val_cols <- grep(value_pattern, names(tr), value = TRUE, ignore.case = TRUE)
  val_cols <- val_cols[!grepl("variance|conf|status|dispersion|uncertainty", val_cols, ignore.case = TRUE)]
  if (length(val_cols) == 0) return(NULL)
  groups <- vapply(val_cols, classify_sample, character(1))
  if (all(is.na(groups))) {
    rg_info <- file.path(dirname(path), "read_groups.info")
    if (file.exists(rg_info)) {
      info <- utils::read.delim(rg_info, check.names = FALSE, stringsAsFactors = FALSE)
      cond <- if ("condition" %in% names(info)) info$condition else info[[2]]
      if (length(cond) == length(val_cols)) {
        groups <- vapply(as.character(cond), classify_sample, character(1))
        names(val_cols) <- paste(cond, seq_along(cond), sep = "_rep")
        colnames_keep <- names(val_cols)
      }
    }
  }
  keep <- !is.na(groups)
  if (sum(keep) < 2) return(NULL)
  mat <- as.matrix(tr[, val_cols[keep], drop = FALSE])
  storage.mode(mat) <- "double"
  mat[is.na(mat)] <- 0
  if (exists("colnames_keep")) colnames(mat) <- colnames_keep[keep]
  mat <- collapse_by_gene(mat, as.character(tr[[gene_col]]))
  list(
    mat = mat,
    sample_info = data.frame(sample = colnames(mat), group = unname(groups[keep]), stringsAsFactors = FALSE),
    source = basename(path)
  )
}

load_expression <- function(project_dir) {
  xlsx <- file.path(project_dir, "shTG(20220723).xlsx")
  if (file.exists(xlsx)) {
    obj <- tryCatch(read_excel_matrix(xlsx), error = function(e) {
      log_msg("Excel import failed: ", e$message)
      NULL
    })
    if (!is.null(obj)) return(obj)
  }
  rg <- file.path(project_dir, "genes.read_group_tracking")
  if (file.exists(rg)) {
    obj <- tryCatch(read_read_group_tracking(rg), error = function(e) {
      log_msg("read_group_tracking import failed: ", e$message)
      NULL
    })
    if (!is.null(obj)) return(obj)
  }
  ct <- file.path(project_dir, "genes.count_tracking")
  if (file.exists(ct)) {
    obj <- read_tracking_matrix(ct, "_count$|^q[0-9]+_count$")
    if (!is.null(obj)) return(obj)
  }
  fp <- file.path(project_dir, "genes.fpkm_tracking")
  if (file.exists(fp)) {
    obj <- read_tracking_matrix(fp, "_FPKM$|^q[0-9]+_FPKM$")
    if (!is.null(obj)) return(obj)
  }
  stop("未找到可用表达矩阵。请确认项目目录含 shTG(20220723).xlsx 或 Cuffdiff tracking 文件: ", project_dir)
}

detect_value_type <- function(mat) {
  x <- as.numeric(mat)
  x <- x[is.finite(x)]
  frac_int <- mean(abs(x - round(x)) < 1e-6)
  if (frac_int > 0.85 && stats::quantile(x, 0.95, na.rm = TRUE) > 50) {
    return("counts")
  }
  "fpkm"
}

# -----------------------------------------------------------------------------
# 4. 过滤低表达 + 标准化
# -----------------------------------------------------------------------------
filter_low_expression <- function(mat, sample_info, value_type) {
  min_n <- max(2, min(table(sample_info$group)))
  if (value_type == "counts") {
    keep <- rowSums(mat >= 10, na.rm = TRUE) >= min_n
  } else {
    keep <- rowSums(mat > 1, na.rm = TRUE) >= min_n
  }
  if (sum(keep) < 200) {
    keep <- rowSums(mat > 0, na.rm = TRUE) >= min_n
    log_msg("Strict filter left too few genes; fallback to expressed-in-", min_n, "-samples")
  }
  log_msg("Low-expression filter: keep ", sum(keep), " / ", nrow(mat), " genes")
  mat[keep, , drop = FALSE]
}

normalize_expression <- function(mat, sample_info, value_type) {
  group <- factor(sample_info$group, levels = c("NTC", "TG_sh1", "TG_sh5"))
  group <- droplevels(group)
  if (value_type == "counts") {
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(mat),
      colData = data.frame(row.names = colnames(mat), group = group),
      design = ~ group
    )
    dds <- DESeq2::estimateSizeFactors(dds)
    log_msg("DESeq2 size-factor normalization")
    norm_counts <- DESeq2::counts(dds, normalized = TRUE)
    log_mat <- log2(norm_counts + 1)
    vsd <- tryCatch({
      assay(DESeq2::vst(dds, blind = TRUE))
    }, error = function(e) {
      log_msg("vst failed, use log2(norm+1): ", e$message)
      log_mat
    })
    list(dds = dds, log_mat = log_mat, heat_mat = vsd, group = group, method = "DESeq2")
  } else {
    dge <- edgeR::DGEList(counts = pmax(mat, 0), group = group)
    keep <- edgeR::filterByExpr(dge, group = group)
    if (sum(keep) >= 200) {
      mat <- mat[keep, , drop = FALSE]
      dge <- dge[keep, , keep.lib.sizes = FALSE]
    }
    log_msg("Quantile normalize log2(FPKM+1)")
    log_mat <- limma::normalizeBetweenArrays(log2(pmax(mat, 0) + 1), method = "quantile")
    list(dds = NULL, log_mat = log_mat, heat_mat = log_mat, group = group, method = "quantile_logFPKM")
  }
}

# -----------------------------------------------------------------------------
# 5. 差异分析
# -----------------------------------------------------------------------------
safe_ebayes <- function(fit) {
  tryCatch(limma::eBayes(fit, trend = TRUE, robust = TRUE), error = function(e) {
    log_msg("eBayes(trend/robust) failed, fallback: ", e$message)
    tryCatch(limma::eBayes(fit), error = function(e2) {
      log_msg("eBayes failed: ", e2$message)
      fit$p.value <- matrix(NA_real_, nrow = nrow(fit$coefficients), ncol = ncol(fit$coefficients))
      fit$fdr <- fit$p.value
      class(fit) <- c("MArrayLM", class(fit))
      fit
    })
  })
}

run_limma_contrasts <- function(log_mat, sample_info) {
  group <- factor(sample_info$group, levels = c("NTC", "TG_sh1", "TG_sh5"))
  group <- droplevels(group)
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  fit <- limma::lmFit(log_mat, design)
  contrast_list <- list()
  if (all(c("TG_sh1", "NTC") %in% colnames(design))) {
    contrast_list$TG_sh1_vs_NTC <- "TG_sh1 - NTC"
  }
  if (all(c("TG_sh5", "NTC") %in% colnames(design))) {
    contrast_list$TG_sh5_vs_NTC <- "TG_sh5 - NTC"
  }
  if (all(c("TG_sh1", "TG_sh5", "NTC") %in% colnames(design))) {
    contrast_list$TGsh_mean_vs_NTC <- "(TG_sh1 + TG_sh5)/2 - NTC"
  }
  if (length(contrast_list) == 0) stop("无法构建任何对照：请检查样本是否包含 NTC 与 TG_sh1/TG_sh5")
  cons <- unlist(contrast_list, use.names = TRUE)
  cm <- limma::makeContrasts(contrasts = cons, levels = design)
  colnames(cm) <- names(cons)
  fit2 <- safe_ebayes(limma::contrasts.fit(fit, cm))
  out <- lapply(colnames(cm), function(cn) {
    tt <- tryCatch(
      limma::topTable(fit2, coef = cn, number = Inf, sort.by = "none"),
      error = function(e) {
        log_msg("topTable failed for ", cn, ": ", e$message, " ; use coefficient only")
        data.frame(
          logFC = as.numeric(fit2$coefficients[, cn]),
          AveExpr = rowMeans(log_mat),
          t = NA_real_,
          P.Value = NA_real_,
          adj.P.Val = NA_real_,
          B = NA_real_,
          row.names = rownames(log_mat)
        )
      }
    )
    data.frame(
      gene = rownames(tt),
      log2FC = tt$logFC,
      AveExpr = tt$AveExpr,
      pvalue = tt$P.Value,
      padj = tt$adj.P.Val,
      stringsAsFactors = FALSE
    )
  })
  names(out) <- colnames(cm)
  out
}

run_deseq2_contrasts <- function(dds) {
  dds$group <- stats::relevel(factor(dds$group), ref = "NTC")
  DESeq2::design(dds) <- ~ group
  dds <- DESeq2::DESeq(dds)
  out <- list()
  to_df <- function(res) {
    data.frame(
      gene = rownames(res),
      log2FC = as.numeric(res$log2FoldChange),
      AveExpr = as.numeric(res$baseMean),
      pvalue = as.numeric(res$pvalue),
      padj = as.numeric(res$padj),
      stringsAsFactors = FALSE
    )
  }
  levs <- levels(dds$group)
  if (all(c("TG_sh1", "NTC") %in% levs)) {
    out$TG_sh1_vs_NTC <- to_df(DESeq2::results(dds, contrast = c("group", "TG_sh1", "NTC")))
  }
  if (all(c("TG_sh5", "NTC") %in% levs)) {
    out$TG_sh5_vs_NTC <- to_df(DESeq2::results(dds, contrast = c("group", "TG_sh5", "NTC")))
  }
  if (all(c("TG_sh1", "TG_sh5", "NTC") %in% levs)) {
    rn <- DESeq2::resultsNames(dds)
    vec <- setNames(rep(0, length(rn)), rn)
    vec["groupTG_sh1"] <- 0.5
    vec["groupTG_sh5"] <- 0.5
    out$TGsh_mean_vs_NTC <- to_df(DESeq2::results(dds, contrast = vec))
  }
  out
}

build_common_up <- function(de_list) {
  a <- de_list[["TG_sh1_vs_NTC"]]
  b <- de_list[["TG_sh5_vs_NTC"]]
  if (is.null(a) || is.null(b)) return(NULL)
  a <- a[!is.na(a$log2FC), ]
  b <- b[!is.na(b$log2FC), ]
  a_up <- a$gene[a$log2FC > 0]
  b_up <- b$gene[b$log2FC > 0]
  common <- intersect(a_up, b_up)
  aa <- a[match(common, a$gene), ]
  bb <- b[match(common, b$gene), ]
  padj_both <- pmax(aa$padj, bb$padj, na.rm = FALSE)
  data.frame(
    gene = common,
    log2FC = (aa$log2FC + bb$log2FC) / 2,
    log2FC_sh1 = aa$log2FC,
    log2FC_sh5 = bb$log2FC,
    AveExpr = (aa$AveExpr + bb$AveExpr) / 2,
    pvalue = pmax(aa$pvalue, bb$pvalue, na.rm = FALSE),
    padj = padj_both,
    stringsAsFactors = FALSE
  )
}

passes_padj <- function(padj, have_pvalue) {
  if (!have_pvalue) return(rep(TRUE, length(padj)))
  !is.na(padj) & padj < padj_cutoff
}

select_by_fc <- function(de, fc, have_pvalue) {
  keep <- !is.na(de$log2FC) & (2^de$log2FC >= fc) & passes_padj(de$padj, have_pvalue)
  if ("log2FC_sh1" %in% names(de)) {
    keep <- keep & (2^de$log2FC_sh1 >= fc) & (2^de$log2FC_sh5 >= fc)
  }
  de[keep, , drop = FALSE]
}

select_by_topn <- function(de, n, have_pvalue) {
  x <- de[!is.na(de$log2FC) & de$log2FC > 0, , drop = FALSE]
  if ("log2FC_sh1" %in% names(x)) {
    x <- x[x$log2FC_sh1 > 0 & x$log2FC_sh5 > 0, , drop = FALSE]
  }
  sig <- x[passes_padj(x$padj, have_pvalue), , drop = FALSE]
  if (nrow(sig) == 0) sig <- x
  sig <- sig[order(sig$log2FC, decreasing = TRUE), , drop = FALSE]
  utils::head(sig, n)
}

# -----------------------------------------------------------------------------
# 6. 基因 ID 转换
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
      clusterProfiler::bitr(symbols, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
      error = function(e) data.frame(ENSEMBL = character(), ENTREZID = character())
    )
    if (nrow(m) > 0) names(m)[1] <- "SYMBOL"
  }
  if (nrow(m) == 0) return(data.frame(gene = character(), entrez = character()))
  m <- m[!duplicated(m[[1]]), ]
  data.frame(gene = m[[1]], entrez = m[[2]], stringsAsFactors = FALSE)
}

ranked_entrez <- function(de) {
  mp <- map_to_entrez(de$gene)
  de2 <- merge(de, mp, by = "gene")
  de2 <- de2[!is.na(de2$entrez) & !is.na(de2$log2FC), ]
  de2 <- de2[order(abs(de2$log2FC), decreasing = TRUE), ]
  de2 <- de2[!duplicated(de2$entrez), ]
  stats <- de2$log2FC
  names(stats) <- de2$entrez
  sort(stats, decreasing = TRUE)
}

# -----------------------------------------------------------------------------
# 7. 绘图
# -----------------------------------------------------------------------------
save_gg <- function(plot, path_stub, width = 8, height = 6) {
  ggplot2::ggsave(paste0(path_stub, ".pdf"), plot, width = width, height = height)
  ggplot2::ggsave(paste0(path_stub, ".png"), plot, width = width, height = height, dpi = 300)
}

plot_volcano <- function(de, highlight, title, outfile, fc_line = 1) {
  df <- de
  df$neglogp <- ifelse(is.na(df$pvalue), 0, -log10(pmax(df$pvalue, 1e-300)))
  df$set <- ifelse(df$gene %in% highlight, "selected", "other")
  df$label <- ifelse(df$gene %in% utils::head(highlight, 15), df$gene, NA)
  lfc_line <- log2(fc_line)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = neglogp, color = set)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.4) +
    ggplot2::scale_color_manual(values = c(other = "grey70", selected = "#D62828")) +
    ggplot2::geom_vline(xintercept = c(-lfc_line, lfc_line), linetype = 2, color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey40") +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, x = "log2 Fold Change", y = "-log10(p value)", color = NULL)
  save_gg(p, outfile)
}

plot_scatter_common <- function(de, highlight, title, outfile) {
  if (!all(c("log2FC_sh1", "log2FC_sh5") %in% names(de))) return(invisible(NULL))
  df <- de
  df$set <- ifelse(df$gene %in% highlight, "common_up", "other")
  df$label <- ifelse(df$gene %in% utils::head(highlight, 15), df$gene, NA)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC_sh1, y = log2FC_sh5, color = set)) +
    ggplot2::geom_point(alpha = 0.75, size = 1.6) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
    ggplot2::scale_color_manual(values = c(other = "grey70", common_up = "#D62828")) +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, x = "log2FC TG_sh1 vs NTC", y = "log2FC TG_sh5 vs NTC", color = NULL)
  save_gg(p, outfile)
}

plot_heatmap <- function(heat_mat, sample_info, genes, title, outfile) {
  genes <- intersect(genes, rownames(heat_mat))
  if (length(genes) < 2) {
    log_msg("Heatmap skipped (<2 genes): ", title)
    return(invisible(NULL))
  }
  sub <- heat_mat[genes, , drop = FALSE]
  ann <- data.frame(Group = sample_info$group, row.names = sample_info$sample)
  ann <- ann[colnames(sub), , drop = FALSE]
  pal <- c(NTC = "#4C78A8", TG_sh1 = "#F58518", TG_sh5 = "#54A24B")
  draw_hm <- function() {
    args <- list(
      mat = sub, scale = "row", annotation_col = ann,
      annotation_colors = list(Group = pal[names(pal) %in% unique(ann$Group)]),
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
  draw_hm()
  grDevices::dev.off()
  grDevices::png(paste0(outfile, ".png"), width = 2400, height = max(1800, 40 * nrow(sub) + 400), res = 300)
  draw_hm()
  grDevices::dev.off()
}

plot_pca <- function(heat_mat, sample_info, outfile) {
  if (ncol(heat_mat) < 2) {
    log_msg("PCA skipped: fewer than 2 samples")
    return(invisible(NULL))
  }
  pca <- tryCatch(stats::prcomp(t(heat_mat), scale. = TRUE), error = function(e) {
    log_msg("PCA failed: ", e$message)
    NULL
  })
  if (is.null(pca)) return(invisible(NULL))
  df <- data.frame(pca$x[, 1:2], group = sample_info$group[match(rownames(pca$x), sample_info$sample)],
                   sample = rownames(pca$x))
  varp <- summary(pca)$importance[2, 1:2] * 100
  p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, color = group, label = sample)) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = "PCA after normalization",
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
# 8. 富集分析与绘图（每个子集都必须出图）
# -----------------------------------------------------------------------------
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
  writeLines(msg, paste0(stub, "_EMPTY.txt"))
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
    p <- tryCatch(p + ggplot2::facet_grid(. ~ .sign) + ggplot2::ggtitle(title), error = function(e) p + ggplot2::ggtitle(title))
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
  out$Reactome <- gsea_one(function(pvalueCutoff, minGSSize) {
    ReactomePA::gsePathway(
      stats, organism = "human", minGSSize = minGSSize, maxGSSize = 500,
      pvalueCutoff = pvalueCutoff, verbose = FALSE, eps = 0
    )
  }, "gsePathway")
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

plot_fgsea_hallmark <- function(stats, outdir, title) {
  if (length(stats) < 5) {
    note_empty(file.path(outdir, "GSEA_Hallmark_fgsea"), "too few ranked genes")
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
    note_empty(file.path(outdir, "GSEA_Hallmark_fgsea"), "no fgsea terms")
    return(invisible(NULL))
  }
  fg <- as.data.frame(fg)
  fg <- fg[order(fg$pval), ]
  utils::write.csv(fg, file.path(outdir, "GSEA_Hallmark_fgsea.csv"), row.names = FALSE)
  plot_df <- utils::head(fg, 15)
  plot_df$pathway <- factor(plot_df$pathway, levels = rev(plot_df$pathway))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = NES, y = pathway, fill = padj < 0.05)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#D62828", "FALSE" = "grey70")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title, y = NULL, fill = "padj < 0.05")
  save_gg(p, file.path(outdir, "GSEA_Hallmark_fgsea_barplot"), 10, 7)
  top_ids <- utils::head(fg$pathway[is.finite(fg$NES)], 3)
  for (i in seq_along(top_ids)) {
    pid <- top_ids[i]
    pe <- tryCatch(fgsea::plotEnrichment(pathways[[pid]], stats) + ggplot2::labs(title = pid), error = function(e) NULL)
    if (!is.null(pe)) save_gg(pe, file.path(outdir, paste0("GSEA_Hallmark_enrichment_top", i)), 8, 5)
  }
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

plot_kegg_pathview <- function(kegg_obj, stats, outdir) {
  if (is.null(kegg_obj) || nrow(as.data.frame(kegg_obj)) == 0) return(invisible(NULL))
  if (!requireNamespace("pathview", quietly = TRUE)) return(invisible(NULL))
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

run_ora_plots <- function(genes, de_sub, outdir) {
  go_dir <- file.path(outdir, "GO")
  pw_dir <- file.path(outdir, "Pathway")
  kg_dir <- file.path(outdir, "KEGG")
  dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(kg_dir, recursive = TRUE, showWarnings = FALSE)

  mp <- map_to_entrez(genes)
  entrez <- unique(mp$entrez)
  fc_sym <- setNames(de_sub$log2FC, de_sub$gene)
  fc_entrez <- setNames(de_sub$log2FC[match(mp$gene, de_sub$gene)], mp$entrez)
  if (length(entrez) < 3) {
    log_msg("ORA skipped, mapped genes < 3: ", outdir)
    writeLines(paste("mapped_entrez", length(entrez)), file.path(outdir, "ORA_skipped.txt"))
    note_empty(file.path(go_dir, "GO"), "too few mapped genes")
    note_empty(file.path(pw_dir, "Pathway"), "too few mapped genes")
    note_empty(file.path(kg_dir, "KEGG"), "too few mapped genes")
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
    plot_ora_object(ego, file.path(go_dir, paste0("GO_", ont)),
                    title_maybe_relaxed(ego, paste("GO", ont)), fold_change = fc_sym)
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
  plot_ora_object(ek, file.path(kg_dir, "KEGG"), title_maybe_relaxed(ek, "KEGG"), fold_change = fc_sym)
  plot_kegg_pathview(ek, fc_entrez, kg_dir)

  er <- enrich_or_relax(
    function() ReactomePA::enrichPathway(
      gene = entrez, organism = "human", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
    ),
    function() ReactomePA::enrichPathway(
      gene = entrez, organism = "human", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
    ),
    "enrichPathway"
  )
  plot_ora_object(er, file.path(pw_dir, "Reactome"), title_maybe_relaxed(er, "Reactome"), fold_change = fc_sym)

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
  plot_ora_object(hm, file.path(pw_dir, "MSigDB_Hallmark"),
                  title_maybe_relaxed(hm, "MSigDB Hallmark"), fold_change = fc_sym)
}

run_gsea_plots <- function(sub, gsea_cache, outdir, tag) {
  gsea_dir <- file.path(outdir, "GSEA")
  dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
  sub_stats <- ranked_entrez(sub)
  plot_fgsea_hallmark(sub_stats, gsea_dir, paste("GSEA Hallmark |", tag))

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
      plot_gsea_object(hm, file.path(gsea_dir, "GSEA_subset_Hallmark"),
                       paste("GSEA Hallmark |", tag, "(subset ranked)"))
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
    plot_gsea_object(kegg, file.path(gsea_dir, "GSEA_subset_KEGG"),
                     paste("GSEA KEGG |", tag, "(subset ranked)"))
  }

  mp <- map_to_entrez(sub$gene)
  for (nm in c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Hallmark")) {
    ids <- gsea_ids_overlapping_genes(gsea_cache[[nm]], sub$gene, mp$entrez)
    plot_gsea_selected_ids(
      gsea_cache[[nm]], ids,
      file.path(gsea_dir, paste0("GSEA_fullrank_overlap_", nm)),
      paste("GSEA", nm, "|", tag, "(full-rank overlap)")
    )
  }
}

emit_subset_analysis <- function(comp_name, sub, tag, title, outdir, full_de_for_volcano,
                                 heat_mat, sample_info, gsea_cache, fc_line = 1) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(sub, file.path(outdir, "DE_selected_genes.csv"), row.names = FALSE)
  writexl::write_xlsx(sub, file.path(outdir, "DE_selected_genes.xlsx"))
  log_msg(comp_name, " ", tag, ": n = ", nrow(sub))
  if (nrow(sub) == 0) {
    writeLines("no genes", file.path(outdir, "EMPTY.txt"))
    return(invisible(NULL))
  }
  plot_de_bar(sub, paste0(title, " | DE genes"), file.path(outdir, "DE_log2FC_barplot"))
  plot_volcano(full_de_for_volcano, sub$gene, title, file.path(outdir, "volcano"), fc_line = fc_line)
  if ("log2FC_sh1" %in% names(full_de_for_volcano)) {
    plot_scatter_common(full_de_for_volcano, sub$gene, title, file.path(outdir, "scatter_sh1_sh5"))
  }
  plot_heatmap(heat_mat, sample_info, sub$gene, title, file.path(outdir, "heatmap"))
  run_ora_plots(sub$gene, sub, outdir)
  run_gsea_plots(sub, gsea_cache, outdir, tag)
}

# -----------------------------------------------------------------------------
# 9. 对单个比较执行全部子集分析（4 个 FC x 7 个 topN，全部出图）
# -----------------------------------------------------------------------------
analyze_one_comparison <- function(comp_name, de, full_de_for_volcano, heat_mat, sample_info, have_pvalue, gsea_de = NULL) {
  if (is.null(gsea_de)) gsea_de <- de
  base <- file.path(result_dir, comp_name)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(de, file.path(base, "DE_full.csv"), row.names = FALSE)
  writexl::write_xlsx(de, file.path(base, "DE_full.xlsx"))

  log_msg("Building GSEA cache for ", comp_name)
  gsea_cache <- build_gsea_cache(gsea_de)
  full_gsea_dir <- file.path(base, "GSEA_full_ranked")
  dir.create(full_gsea_dir, recursive = TRUE, showWarnings = FALSE)
  for (nm in c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Hallmark")) {
    plot_gsea_object(gsea_cache[[nm]], file.path(full_gsea_dir, paste0("GSEA_", nm)),
                     paste("GSEA", nm, "|", comp_name, "| full ranked list"))
  }
  plot_fgsea_hallmark(gsea_cache$stats, full_gsea_dir, paste("GSEA Hallmark |", comp_name))

  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    sub <- select_by_fc(de, fc, have_pvalue)
    emit_subset_analysis(
      comp_name, sub, nm, paste0(comp_name, " | up FC >= ", fc),
      file.path(base, "FoldChange", nm),
      full_de_for_volcano, heat_mat, sample_info, gsea_cache, fc_line = fc
    )
  }

  for (n in top_ns) {
    tag <- paste0("top", n)
    sub <- select_by_topn(de, n, have_pvalue)
    emit_subset_analysis(
      comp_name, sub, tag, paste0(comp_name, " | upregulated top ", n),
      file.path(base, "TopRank", tag),
      full_de_for_volcano, heat_mat, sample_info, gsea_cache, fc_line = 1
    )
  }
}

plot_venn_up <- function(de_list, have_pvalue) {
  if (is.null(de_list$TG_sh1_vs_NTC) || is.null(de_list$TG_sh5_vs_NTC)) return(invisible(NULL))
  outdir <- file.path(result_dir, "common_up")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    a <- select_by_fc(de_list$TG_sh1_vs_NTC, fc, have_pvalue)$gene
    b <- select_by_fc(de_list$TG_sh5_vs_NTC, fc, have_pvalue)$gene
    lst <- list(TG_sh1_up = a, TG_sh5_up = b)
    p <- tryCatch(ggvenn::ggvenn(lst, fill_color = c("#F58518", "#54A24B")) +
                    ggplot2::labs(title = paste0("Common up | FC >= ", fc)), error = function(e) NULL)
    vdir <- file.path(outdir, "FoldChange", nm)
    dir.create(vdir, recursive = TRUE, showWarnings = FALSE)
    if (!is.null(p)) save_gg(p, file.path(vdir, paste0("venn_", nm)), width = 7, height = 6)
  }
  for (n in top_ns) {
    tag <- paste0("top", n)
    a <- select_by_topn(de_list$TG_sh1_vs_NTC, n, have_pvalue)$gene
    b <- select_by_topn(de_list$TG_sh5_vs_NTC, n, have_pvalue)$gene
    lst <- list(TG_sh1_up = a, TG_sh5_up = b)
    p <- tryCatch(ggvenn::ggvenn(lst, fill_color = c("#F58518", "#54A24B")) +
                    ggplot2::labs(title = paste0("Common up | ", tag)), error = function(e) NULL)
    vdir <- file.path(outdir, "TopRank", tag)
    dir.create(vdir, recursive = TRUE, showWarnings = FALSE)
    if (!is.null(p)) save_gg(p, file.path(vdir, paste0("venn_", tag)), width = 7, height = 6)
  }
}

# -----------------------------------------------------------------------------
# 10. 主流程
# -----------------------------------------------------------------------------
log_msg("Project dir: ", project_dir)
expr <- load_expression(project_dir)
log_msg("Loaded from ", expr$source, " | genes=", nrow(expr$mat), " samples=", ncol(expr$mat))
print(expr$sample_info)
utils::write.csv(expr$sample_info, file.path(log_dir, "sample_info.csv"), row.names = FALSE)

present <- unique(expr$sample_info$group)
log_msg("Detected groups: ", paste(present, collapse = ", "))
if (!"NTC" %in% present) stop("未检测到 NTC 对照样本，请检查表达矩阵列名")
if (!"TG_sh1" %in% present) log_msg("WARNING: 未检测到 TG_sh1")
if (!"TG_sh5" %in% present) log_msg("WARNING: 未检测到 TG_sh5；将跳过 sh5 相关比较")

value_type <- detect_value_type(expr$mat)
log_msg("Value type inferred as: ", value_type)

filt <- filter_low_expression(expr$mat, expr$sample_info, value_type)
norm <- normalize_expression(filt, expr$sample_info, value_type)
expr$sample_info <- expr$sample_info[match(colnames(norm$log_mat), expr$sample_info$sample), ]
utils::write.csv(
  cbind(gene = rownames(norm$log_mat), as.data.frame(norm$log_mat)),
  file.path(result_dir, "normalized_log_matrix.csv"),
  row.names = FALSE
)
plot_pca(norm$heat_mat, expr$sample_info, file.path(result_dir, "00_QC_PCA"))

have_pvalue <- TRUE
if (!is.null(norm$dds) && value_type == "counts") {
  log_msg("Differential expression: DESeq2")
  de_list <- tryCatch(run_deseq2_contrasts(norm$dds), error = function(e) {
    log_msg("DESeq2 failed, fallback to limma: ", e$message)
    run_limma_contrasts(norm$log_mat, expr$sample_info)
  })
} else {
  log_msg("Differential expression: limma on normalized log matrix")
  de_list <- run_limma_contrasts(norm$log_mat, expr$sample_info)
}
have_pvalue <- any(vapply(de_list, function(x) any(!is.na(x$padj)), logical(1)))
if (!have_pvalue) log_msg("No p-values available (likely too few replicates); FC/rank filters only")

common <- build_common_up(de_list)
if (!is.null(common)) de_list$common_up <- common

common_full_rank <- NULL
if (!is.null(de_list$TG_sh1_vs_NTC) && !is.null(de_list$TG_sh5_vs_NTC)) {
  sh1 <- de_list$TG_sh1_vs_NTC
  sh5 <- de_list$TG_sh5_vs_NTC
  both <- merge(
    sh1[, c("gene", "log2FC", "AveExpr", "pvalue", "padj")],
    sh5[, c("gene", "log2FC", "AveExpr", "pvalue", "padj")],
    by = "gene", suffixes = c("_sh1", "_sh5")
  )
  both$log2FC <- (both$log2FC_sh1 + both$log2FC_sh5) / 2
  both$AveExpr <- (both$AveExpr_sh1 + both$AveExpr_sh5) / 2
  both$pvalue <- pmax(both$pvalue_sh1, both$pvalue_sh5, na.rm = FALSE)
  both$padj <- pmax(both$padj_sh1, both$padj_sh5, na.rm = FALSE)
  common_full_rank <- both
}

for (nm in names(de_list)) {
  volcano_df <- de_list[[nm]]
  gsea_de <- de_list[[nm]]
  if (nm == "common_up" && !is.null(common_full_rank)) {
    volcano_df <- common_full_rank
    gsea_de <- common_full_rank
  }
  analyze_one_comparison(
    nm, de_list[[nm]], volcano_df, norm$heat_mat, expr$sample_info, have_pvalue, gsea_de
  )
}
plot_venn_up(de_list, have_pvalue)

utils::writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
log_msg("All done. Results in: ", result_dir)
