#!/usr/bin/env Rscript
# NTC / TG_sh1 / TG_sh5 细胞 RNA-seq：
#   从 Cuffdiff tracking 重建表达矩阵（不直接用 Excel）
#   分析 1：TG_sh1 vs NTC、TG_sh5 vs NTC 各自 DEG + 火山图 + GO
#   分析 2：TG_sh1 与 TG_sh5 组均值再平均，相对 NTC 做 DEG + 火山图 + GO
#   GO + KEGG + Reactome/WikiPathways 仅上调基因，按线性 FC>=1 / 1.25 / 1.5 / 2 四组分别富集
#
# 用法:
#   Rscript scripts/tg_vs_ntc_deg_go.R [data_dir] [out_dir]
# 默认: E:/R/TG_BRCA/TG

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else "E:/R/TG_BRCA/TG"
out_dir  <- if (length(args) >= 2) args[[2]] else file.path(data_dir, "results")

control_group        <- "NTC"
treat_groups         <- c("TG_sh1", "TG_sh5")
pooled_treat_label   <- "TG_sh1_sh5_avg"
pooled_deseq_level   <- "TG"
padj_cutoff          <- 0.05
log2fc_cutoff        <- 1
# GO：仅上调；线性 FC 阈值（1 倍 = FC>=1，不是 log2FC>=1）
go_fc_cutoffs        <- c("up_FC1" = 1, "up_FC1.25" = 1.25, "up_FC1.5" = 1.5, "up_FC2" = 2)
go_ontologies        <- c("BP", "MF", "CC")
min_replicates       <- 2
species_orgdb        <- "org.Hs.eg.db"

dir.create(file.path(out_dir, "matrix"), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
need_pkg <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) return(invisible(TRUE))
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(missing, update = FALSE, ask = FALSE)
  } else {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
  invisible(TRUE)
}

find_input <- function(dir, pattern, required = TRUE) {
  hits <- list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (!length(hits)) {
    if (required) stop("找不到匹配 '", pattern, "' 的文件，目录: ", dir)
    return(NA_character_)
  }
  hits[[1]]
}

normalize_deg_names <- function(dt) {
  nms <- gsub("[()]", "", names(dt))
  nms <- gsub(" ", "_", nms)
  setnames(dt, nms)
  if ("log2fold_change" %in% names(dt)) setnames(dt, "log2fold_change", "log2FC")
  if ("log2_fold_change" %in% names(dt)) setnames(dt, "log2_fold_change", "log2FC")
  if ("q_value" %in% names(dt) && !"padj" %in% names(dt)) dt[, padj := q_value]
  if ("pvalue" %in% names(dt) && !"p_value" %in% names(dt)) dt[, p_value := pvalue]
  if ("gene_short_name" %in% names(dt) && !"gene" %in% names(dt)) {
    setnames(dt, "gene_short_name", "gene")
  }
  dt
}

read_cuffdiff_deg <- function(path, control, treat) {
  cuff <- fread(path)
  cuff <- normalize_deg_names(cuff)
  if ("status" %in% names(cuff)) cuff <- cuff[status == "OK"]
  if (all(c("sample_1", "sample_2") %in% names(cuff))) {
    cuff <- cuff[
      (sample_1 == control & sample_2 == treat) |
        (sample_1 == treat & sample_2 == control)
    ]
    cuff[sample_1 == treat & sample_2 == control, log2FC := -as.numeric(log2FC)]
  }
  if (!"log2FC" %in% names(cuff) || !nrow(cuff)) {
    return(NULL)
  }
  cuff
}

sanitize_gene_symbols <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x %in% c("", "-", ".", "NA")] <- NA_character_
  x
}

split_cuff_genes <- function(x) {
  x <- sanitize_gene_symbols(x)
  x <- x[!is.na(x)]
  unique(unlist(strsplit(x, ",\\s*"), use.names = FALSE))
}

write_tsv <- function(x, path) {
  fwrite(x, path, sep = "\t")
  message("Wrote ", path)
}

safe_rowmeans <- function(mat, cols) {
  cols <- intersect(cols, colnames(mat))
  if (!length(cols)) {
    return(rep(NA_real_, nrow(mat)))
  }
  rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
}

annotate_direction <- function(deg) {
  if (!"log2FC" %in% names(deg)) stop("deg 缺少 log2FC")
  deg[, log2FC := as.numeric(log2FC)]
  deg[, log2FC_plot := log2FC]
  finite_abs <- abs(deg$log2FC[is.finite(deg$log2FC)])
  max_finite <- if (length(finite_abs)) max(finite_abs, na.rm = TRUE) else NA_real_
  if (is.finite(max_finite)) {
    deg[!is.finite(log2FC_plot) & log2FC_plot > 0, log2FC_plot := max_finite]
    deg[!is.finite(log2FC_plot) & log2FC_plot < 0, log2FC_plot := -max_finite]
  }
  if (!"padj" %in% names(deg) && "q_value" %in% names(deg)) deg[, padj := q_value]
  if (!"p_value" %in% names(deg) && "pvalue" %in% names(deg)) deg[, p_value := pvalue]
  if (!"padj" %in% names(deg)) deg[, padj := NA_real_]
  deg[, direction := "NS"]
  deg[!is.na(padj) & padj < padj_cutoff & !is.na(log2FC) & log2FC >= log2fc_cutoff, direction := "Up"]
  deg[!is.na(padj) & padj < padj_cutoff & !is.na(log2FC) & log2FC <= -log2fc_cutoff, direction := "Down"]
  # 无 p 值时仍按 |log2FC| 标方向，便于平均后仅有倍数变化的情况
  deg[is.na(padj) & !is.na(log2FC) & log2FC >= log2fc_cutoff, direction := "Up"]
  deg[is.na(padj) & !is.na(log2FC) & log2FC <= -log2fc_cutoff, direction := "Down"]
  deg[, fold_change := 2^log2FC]
  deg
}

select_up_by_fc <- function(deg, fc_cutoff) {
  log2_cut <- log2(fc_cutoff)
  up <- deg[!is.na(log2FC) & is.finite(log2FC) & log2FC >= log2_cut]
  if ("padj" %in% names(up) && any(!is.na(up$padj))) {
    up <- up[!is.na(padj) & padj < padj_cutoff]
  }
  up
}

plot_volcano <- function(deg, sig, title, pdf_path) {
  need_pkg(c("ggplot2", "ggrepel"))
  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
  })
  volc <- deg[is.finite(log2FC_plot)]
  if (!nrow(volc)) {
    warning("火山图无有效点: ", title)
    return(invisible(NULL))
  }
  if ("padj" %in% names(volc) && any(!is.na(volc$padj))) {
    volc[, neglog10_padj := -log10(pmax(fifelse(is.na(padj), 1, padj), .Machine$double.xmin))]
    y_lab <- "-log10 adjusted p-value"
  } else {
    volc[, neglog10_padj := 0]
    y_lab <- "-log10 adjusted p-value (unavailable)"
  }

  label_dt <- rbindlist(list(
    sig[direction == "Up"][order(padj, -abs(log2FC))][seq_len(min(6L, .N))],
    sig[direction == "Down"][order(padj, -abs(log2FC))][seq_len(min(6L, .N))]
  ), fill = TRUE)
  if (nrow(label_dt)) {
    first_symbol <- function(z) {
      s <- split_cuff_genes(z)
      if (!length(s)) NA_character_ else s[[1]]
    }
    if ("gene" %in% names(label_dt)) {
      label_dt[, label := vapply(as.character(gene), first_symbol, character(1))]
    } else if ("tracking_id" %in% names(label_dt)) {
      label_dt[, label := tracking_id]
    } else {
      label_dt[, label := NA_character_]
    }
    label_dt <- label_dt[!is.na(label)]
  }

  p_volc <- ggplot(volc, aes(x = log2FC_plot, y = neglog10_padj, color = direction)) +
    geom_point(alpha = 0.7, size = 1.2) +
    geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = 2, linewidth = 0.3) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = 2, linewidth = 0.3) +
    scale_color_manual(values = c(Up = "#D62728", Down = "#1F77B4", NS = "grey70")) +
    labs(title = title, x = "log2 fold change", y = y_lab) +
    theme_bw(base_size = 12) +
    theme(legend.title = element_blank())

  if (nrow(label_dt)) {
    p_volc <- p_volc + geom_text_repel(
      data = label_dt,
      aes(
        x = log2FC_plot,
        y = -log10(pmax(fifelse(is.na(padj), 1, padj), .Machine$double.xmin)),
        label = label
      ),
      size = 3, max.overlaps = 20, show.legend = FALSE, color = "black"
    )
  }
  ggsave(pdf_path, p_volc, width = 7, height = 6)
  message("Wrote ", pdf_path)
}

run_go_analysis <- function(deg, sig, contrast_id, dest_dir) {
  # GO 只用上调基因；线性 FC 四组各自富集，不用 sig 的上+下调合并集
  invisible(sig)
  need_pkg("ggplot2")
  need_pkg(c("clusterProfiler", "enrichplot", "ReactomePA", species_orgdb), bioc = TRUE)
  suppressPackageStartupMessages({
    library(ggplot2)
    library(clusterProfiler)
    library(enrichplot)
    library(org.Hs.eg.db)
  })
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  uni_col <- if ("gene" %in% names(deg)) deg$gene else deg$tracking_id
  universe_symbols <- unique(split_cuff_genes(uni_col))

  map_ids <- function(symbols) {
    if (!length(symbols)) {
      return(data.frame(SYMBOL = character(), ENTREZID = character()))
    }
    suppressMessages(bitr(
      symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = species_orgdb
    ))
  }
  uni_map <- map_ids(universe_symbols)
  count_rows <- list()

  run_go <- function(entrez, ont, universe) {
    enrichGO(
      gene = unique(entrez),
      universe = unique(universe),
      OrgDb = get(species_orgdb),
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
  }

  kegg_helper <- c(
    "scripts/kegg_enrich_helpers.R",
    {
      fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
      if (length(fa)) file.path(dirname(normalizePath(sub("^--file=", "", fa[[1]]))), "kegg_enrich_helpers.R") else NA_character_
    }
  )
  kegg_helper <- kegg_helper[!is.na(kegg_helper) & file.exists(kegg_helper)][1]
  if (is.na(kegg_helper)) stop("找不到 kegg_enrich_helpers.R")
  if (!exists("save_kegg_and_pathways", mode = "function")) source(kegg_helper, local = FALSE)

  for (tag in names(go_fc_cutoffs)) {
    fc <- as.numeric(go_fc_cutoffs[[tag]])
    up <- select_up_by_fc(deg, fc)
    tag_dir <- file.path(dest_dir, tag)
    dir.create(tag_dir, recursive = TRUE, showWarnings = FALSE)

    gene_col <- if ("gene" %in% names(up)) up$gene else up$tracking_id
    up_symbols <- unique(split_cuff_genes(gene_col))
    write_tsv(up, file.path(tag_dir, paste0("up_genes_", tag, "_", contrast_id, ".tsv")))
    message(contrast_id, " ", tag, " (linear FC>=", fc, "): ", length(up_symbols), " up symbols")
    count_rows[[tag]] <- data.table(
      contrast = contrast_id,
      cutoff = tag,
      linear_FC_min = fc,
      n_up_rows = nrow(up),
      n_up_symbols = length(up_symbols)
    )

    if (length(up_symbols) < 5) {
      message(contrast_id, " skip GO/KEGG/pathway for ", tag, " (need >=5 unique symbols)")
      next
    }

    sig_map <- map_ids(up_symbols)
    message(contrast_id, " ", tag, ": mapped ", nrow(sig_map), "/", length(up_symbols), " symbols")
    if (!nrow(sig_map)) next

    for (ont in go_ontologies) {
      ego <- run_go(sig_map$ENTREZID, ont, uni_map$ENTREZID)
      if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
        message(contrast_id, " ", tag, " GO ", ont, ": no terms")
        next
      }
      write_tsv(
        as.data.table(as.data.frame(ego)),
        file.path(tag_dir, paste0("GO_", ont, "_", tag, "_", contrast_id, ".tsv"))
      )
      pdf(file.path(tag_dir, paste0("GO_", ont, "_", tag, "_dotplot.pdf")), width = 8, height = 6)
      print(dotplot(ego, showCategory = 15) + ggtitle(paste(contrast_id, tag, "GO", ont)))
      dev.off()
      pdf(file.path(tag_dir, paste0("GO_", ont, "_", tag, "_barplot.pdf")), width = 8, height = 6)
      print(barplot(ego, showCategory = 15) + ggtitle(paste(contrast_id, tag, "GO", ont)))
      dev.off()
    }

    save_kegg_and_pathways(
      entrez = sig_map$ENTREZID,
      universe = uni_map$ENTREZID,
      dest_dir = tag_dir,
      tag = tag,
      title_prefix = paste(contrast_id, tag)
    )
  }
  write_tsv(rbindlist(count_rows), file.path(dest_dir, paste0("up_FC_gene_counts_", contrast_id, ".tsv")))
}

write_deg_outputs <- function(deg, contrast_id, dest_dir, deg_source, treat, control) {
  dir.create(file.path(dest_dir, "deg"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dest_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dest_dir, "go"), recursive = TRUE, showWarnings = FALSE)

  deg <- annotate_direction(deg)
  deg[, `:=`(deg_source = deg_source, contrast = contrast_id)]
  write_tsv(deg, file.path(dest_dir, "deg", paste0("deg_", contrast_id, ".tsv")))
  sig <- deg[direction %in% c("Up", "Down")]
  write_tsv(sig, file.path(dest_dir, "deg", paste0("sig_deg_", contrast_id, ".tsv")))
  message(contrast_id, " significant: ", nrow(sig),
          " (Up ", sum(sig$direction == "Up"),
          ", Down ", sum(sig$direction == "Down"), ") via ", deg_source)

  plot_volcano(
    deg, sig,
    title = paste0(treat, " vs ", control, " (", deg_source, ")"),
    pdf_path = file.path(dest_dir, "plots", paste0("volcano_", contrast_id, ".pdf"))
  )
  run_go_analysis(deg, sig, contrast_id, file.path(dest_dir, "go"))
}

deseq2_contrast <- function(count_mat, coldata_df, treat, control) {
  need_pkg("DESeq2", bioc = TRUE)
  suppressPackageStartupMessages(library(DESeq2))
  keep <- coldata_df$condition %in% c(control, treat)
  mat <- count_mat[, coldata_df$sample[keep], drop = FALSE]
  cd <- droplevels(coldata_df[keep, , drop = FALSE])
  rownames(cd) <- cd$sample
  cd$condition <- relevel(factor(cd$condition), ref = control)
  dds <- DESeqDataSetFromMatrix(countData = mat, colData = cd, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", treat, control), alpha = padj_cutoff)
  res_dt <- as.data.table(as.data.frame(res), keep.rownames = "tracking_id")
  res_dt <- merge(res_dt, gene_map[, .(tracking_id, gene_short_name)], by = "tracking_id", all.x = TRUE)
  setnames(
    res_dt,
    old = c("log2FoldChange", "pvalue", "padj", "gene_short_name"),
    new = c("log2FC", "p_value", "padj", "gene"),
    skip_absent = TRUE
  )
  res_dt
}

fpkm_ttest_deg <- function(fpkm_mat, group_a_cols, group_b_cols, treat_mean = NULL) {
  # group_a = control samples; group_b = treat samples (or two shRNA-mean columns)
  a <- fpkm_mat[, group_a_cols, drop = FALSE]
  if (is.null(treat_mean)) {
    b <- fpkm_mat[, group_b_cols, drop = FALSE]
    treat_mean <- rowMeans(b, na.rm = TRUE)
  }
  ctrl_mean <- rowMeans(a, na.rm = TRUE)
  log2fc <- log2((treat_mean + 0.1) / (ctrl_mean + 0.1))
  pvals <- rep(NA_real_, nrow(fpkm_mat))
  if (ncol(a) >= 2 && length(group_b_cols) >= 2) {
    b <- fpkm_mat[, group_b_cols, drop = FALSE]
    for (i in seq_len(nrow(fpkm_mat))) {
      av <- as.numeric(a[i, ])
      bv <- as.numeric(b[i, ])
      if (sum(is.finite(av)) >= 2 && sum(is.finite(bv)) >= 2 &&
          (sd(av, na.rm = TRUE) > 0 || sd(bv, na.rm = TRUE) > 0)) {
        pvals[i] <- tryCatch(
          t.test(bv, av, var.equal = FALSE)$p.value,
          error = function(e) NA_real_
        )
      }
    }
  }
  data.table(
    tracking_id = rownames(fpkm_mat),
    gene = gene_map$gene_short_name[match(rownames(fpkm_mat), gene_map$tracking_id)],
    value_control = ctrl_mean,
    value_treat = treat_mean,
    log2FC = log2fc,
    p_value = pvals,
    padj = p.adjust(pvals, method = "BH")
  )
}

# ---------------------------------------------------------------------------
# Read Cuffdiff originals (not Excel)
# ---------------------------------------------------------------------------
tracking_path <- find_input(data_dir, "^genes\\.read_group_tracking$", required = FALSE)
fpkm_map_path <- find_input(data_dir, "^genes\\.fpkm_tracking$", required = FALSE)
diff_path     <- find_input(data_dir, "^gene_exp\\.diff$", required = FALSE)
rg_info_path  <- find_input(data_dir, "^read_groups\\.info$", required = FALSE)
excel_path    <- find_input(data_dir, "shTG.*\\.(xlsx|xls)$", required = FALSE)
if (is.na(excel_path)) excel_path <- find_input(data_dir, "^shTG", required = FALSE)

if (is.na(tracking_path)) {
  stop("缺少 genes.read_group_tracking。TG_sh5 与合并平均分析必须从该文件重建矩阵；Excel 只有 NTC vs TG_sh1。")
}

rg <- fread(tracking_path)
req_cols <- c("tracking_id", "condition", "replicate", "raw_frags", "FPKM")
missing_cols <- setdiff(req_cols, names(rg))
if (length(missing_cols)) {
  stop("genes.read_group_tracking 缺少列: ", paste(missing_cols, collapse = ", "))
}

rg[, `:=`(
  condition = as.character(condition),
  replicate = as.character(replicate),
  sample = paste(condition, replicate, sep = "_"),
  raw_frags = as.numeric(raw_frags),
  FPKM = as.numeric(FPKM)
)]

if (!is.na(rg_info_path) && file.exists(rg_info_path)) {
  message("Using sample metadata from ", basename(rg_info_path))
}

if (is.na(fpkm_map_path)) {
  gene_map <- unique(rg[, .(tracking_id, gene_short_name = tracking_id)])
} else {
  keep_cols <- intersect(
    c("tracking_id", "gene_id", "gene_short_name", "locus", "gene"),
    names(fread(fpkm_map_path, nrows = 0))
  )
  gene_map <- fread(fpkm_map_path, select = keep_cols)
  if (!"gene_short_name" %in% names(gene_map) && "gene" %in% names(gene_map)) {
    setnames(gene_map, "gene", "gene_short_name")
  }
  if (!"gene_short_name" %in% names(gene_map)) {
    gene_map[, gene_short_name := tracking_id]
  }
}
gene_map[, gene_short_name := sanitize_gene_symbols(gene_short_name)]

count_dt <- dcast(rg, tracking_id ~ sample, value.var = "raw_frags", fun.aggregate = sum)
fpkm_dt  <- dcast(rg, tracking_id ~ sample, value.var = "FPKM", fun.aggregate = mean)
count_dt <- merge(gene_map[, .(tracking_id, gene_short_name)], count_dt, by = "tracking_id", all.y = TRUE)
fpkm_dt  <- merge(gene_map[, .(tracking_id, gene_short_name)], fpkm_dt,  by = "tracking_id", all.y = TRUE)
write_tsv(count_dt, file.path(out_dir, "matrix", "counts_from_read_group_tracking.tsv"))
write_tsv(fpkm_dt,  file.path(out_dir, "matrix", "fpkm_from_read_group_tracking.tsv"))

coldata <- unique(rg[, .(sample, condition, replicate)])
coldata[, condition := factor(
  condition,
  levels = unique(c(control_group, treat_groups, as.character(unique(condition))))
)]
write_tsv(coldata, file.path(out_dir, "matrix", "sample_metadata.tsv"))

n_per <- coldata[, .N, by = condition]
message("Replicates per group:\n", paste(capture.output(n_per), collapse = "\n"))

present <- as.character(unique(coldata$condition))
if (!control_group %in% present) {
  stop("找不到对照组 ", control_group, "。实际分组: ", paste(present, collapse = ", "))
}
missing_treats <- setdiff(treat_groups, present)
if (length(missing_treats)) {
  warning("tracking 中缺少: ", paste(missing_treats, collapse = ", "),
          "。将跳过这些组的 pairwise / 合并分析。")
}

sample_cols <- setdiff(names(count_dt), c("tracking_id", "gene_short_name"))
count_mat <- as.matrix(count_dt[, ..sample_cols])
rownames(count_mat) <- count_dt$tracking_id
count_mat[is.na(count_mat)] <- 0
count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"

fpkm_cols <- setdiff(names(fpkm_dt), c("tracking_id", "gene_short_name"))
fpkm_mat <- as.matrix(fpkm_dt[, ..fpkm_cols])
rownames(fpkm_mat) <- fpkm_dt$tracking_id
storage.mode(fpkm_mat) <- "numeric"

coldata_df <- as.data.frame(coldata)
n_tbl <- setNames(n_per$N, as.character(n_per$condition))

can_deseq_pair <- function(treat) {
  isTRUE(n_tbl[[control_group]] >= min_replicates && n_tbl[[treat]] >= min_replicates)
}

# ---------------------------------------------------------------------------
# Analysis 1: TG_sh1 vs NTC, TG_sh5 vs NTC
# ---------------------------------------------------------------------------
for (treat in treat_groups) {
  contrast_id <- paste0(treat, "_vs_", control_group)
  dest <- file.path(out_dir, "pairwise", contrast_id)
  message("==== Analysis 1: ", contrast_id, " ====")

  if (!treat %in% present) {
    warning("跳过 ", contrast_id, "：矩阵中没有该组")
    next
  }

  deg <- NULL
  deg_source <- NULL
  if (can_deseq_pair(treat)) {
    deg <- deseq2_contrast(count_mat, coldata_df, treat, control_group)
    deg_source <- "DESeq2_from_read_group_tracking"
  } else {
    message(contrast_id, ": 重复不足，尝试 Cuffdiff / Excel / FPKM t-test")
    if (!is.na(diff_path)) {
      deg <- read_cuffdiff_deg(diff_path, control_group, treat)
      if (!is.null(deg)) deg_source <- "cuffdiff_gene_exp.diff"
    }
    if (is.null(deg) && treat == "TG_sh1" && !is.na(excel_path)) {
      need_pkg("readxl")
      deg <- as.data.table(readxl::read_excel(excel_path))
      deg <- normalize_deg_names(deg)
      if ("status" %in% names(deg)) deg <- deg[status == "OK"]
      deg_source <- "excel_fallback_NTC_vs_TG_sh1_only"
      warning("Excel 仅含 NTC vs TG_sh1，已用于 ", contrast_id)
    }
    if (is.null(deg)) {
      ctrl_s <- coldata[condition == control_group, sample]
      trt_s  <- coldata[condition == treat, sample]
      deg <- fpkm_ttest_deg(fpkm_mat, ctrl_s, trt_s)
      deg_source <- "fpkm_mean_ttest"
    }
  }
  write_deg_outputs(deg, contrast_id, dest, deg_source, treat, control_group)
}

# ---------------------------------------------------------------------------
# Analysis 2: average(TG_sh1 mean, TG_sh5 mean) vs NTC
# ---------------------------------------------------------------------------
message("==== Analysis 2: ", pooled_treat_label, " vs ", control_group, " ====")
if (!all(treat_groups %in% present)) {
  stop("合并平均分析需要 TG_sh1 与 TG_sh5 都在 tracking 中。Excel 无法做该分析。")
}

ntc_samples <- coldata[condition == control_group, sample]
sh1_samples <- coldata[condition == "TG_sh1", sample]
sh5_samples <- coldata[condition == "TG_sh5", sample]

sh1_mean <- safe_rowmeans(fpkm_mat, sh1_samples)
sh5_mean <- safe_rowmeans(fpkm_mat, sh5_samples)
ntc_mean <- safe_rowmeans(fpkm_mat, ntc_samples)
tg_avg   <- (sh1_mean + sh5_mean) / 2

avg_dt <- data.table(
  tracking_id = rownames(fpkm_mat),
  gene = gene_map$gene_short_name[match(rownames(fpkm_mat), gene_map$tracking_id)],
  NTC_mean = ntc_mean,
  TG_sh1_mean = sh1_mean,
  TG_sh5_mean = sh5_mean,
  TG_avg = tg_avg,
  log2FC_avg_vs_NTC = log2((tg_avg + 0.1) / (ntc_mean + 0.1))
)
write_tsv(avg_dt, file.path(out_dir, "matrix", "fpkm_TG_sh1_sh5_equal_weight_avg_vs_NTC.tsv"))

# counts: equal-weight shRNA means as two KD "replicates" (not a single averaged column)
sh1_count_mean <- round(safe_rowmeans(count_mat, sh1_samples))
sh5_count_mean <- round(safe_rowmeans(count_mat, sh5_samples))
avg_count_dt <- data.table(
  tracking_id = rownames(count_mat),
  gene = gene_map$gene_short_name[match(rownames(count_mat), gene_map$tracking_id)],
  TG_sh1_mean_count = sh1_count_mean,
  TG_sh5_mean_count = sh5_count_mean,
  TG_avg_count = round((sh1_count_mean + sh5_count_mean) / 2)
)
write_tsv(avg_count_dt, file.path(out_dir, "matrix", "counts_TG_sh1_sh5_equal_weight_avg.tsv"))

pooled_dest <- file.path(out_dir, "pooled_avg", paste0(pooled_treat_label, "_vs_", control_group))

# 敲低侧只用两个等权均值列（TG_sh1 组平均、TG_sh5 组平均），对照保留 NTC 各重复
sh1_count_mean[is.na(sh1_count_mean)] <- 0
sh5_count_mean[is.na(sh5_count_mean)] <- 0
collapsed_counts <- cbind(
  count_mat[, ntc_samples, drop = FALSE],
  TG_sh1_mean = pmax(0L, as.integer(round(sh1_count_mean))),
  TG_sh5_mean = pmax(0L, as.integer(round(sh5_count_mean)))
)
collapsed_cd <- data.frame(
  sample = colnames(collapsed_counts),
  condition = ifelse(colnames(collapsed_counts) %in% ntc_samples, control_group, pooled_deseq_level),
  replicate = colnames(collapsed_counts),
  stringsAsFactors = FALSE
)
can_deseq_pooled <- length(ntc_samples) >= min_replicates

if (can_deseq_pooled) {
  deg_pooled <- deseq2_contrast(
    collapsed_counts, collapsed_cd,
    treat = pooled_deseq_level, control = control_group
  )
  deg_source_pooled <- "DESeq2_equal_weight_shRNA_means_vs_NTC"
} else {
  message("NTC 重复不足，改用两个 shRNA 组均值对 NTC 的 t 检验 / 倍数变化")
  fake_mat <- cbind(fpkm_mat, TG_sh1_mean = sh1_mean, TG_sh5_mean = sh5_mean)
  deg_pooled <- fpkm_ttest_deg(
    fake_mat, ntc_samples, c("TG_sh1_mean", "TG_sh5_mean"),
    treat_mean = tg_avg
  )
  deg_source_pooled <- "fpkm_equal_weight_avg_ttest"
}

write_deg_outputs(
  deg_pooled,
  contrast_id = paste0(pooled_treat_label, "_vs_", control_group),
  dest_dir = pooled_dest,
  deg_source = deg_source_pooled,
  treat = pooled_treat_label,
  control = control_group
)

message("Done. Results in ", normalizePath(out_dir, mustWork = FALSE))
message("Pairwise: ", file.path(out_dir, "pairwise"))
message("Pooled average: ", file.path(out_dir, "pooled_avg"))
