#!/usr/bin/env Rscript
# TG_sh1 vs NTC 细胞 RNA-seq：
#   1) 从 Cuffdiff 原始 tracking 文件重建表达矩阵（不直接用 Excel）
#   2) 差异表达（DESeq2；重复不足时回退 gene_exp.diff）
#   3) 火山图 + GO 富集
#
# 用法:
#   Rscript scripts/tg_sh1_vs_ntc_deg_go.R [data_dir] [out_dir]
# 默认 data_dir / 分析目录: E:/R/TG_BRCA/TG

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else "E:/R/TG_BRCA/TG"
out_dir  <- if (length(args) >= 2) args[[2]] else file.path(data_dir, "results")

control_group   <- "NTC"
treat_group     <- "TG_sh1"
padj_cutoff     <- 0.05
log2fc_cutoff   <- 1
go_ontologies   <- c("BP", "MF", "CC")
min_replicates  <- 2
species_orgdb   <- "org.Hs.eg.db"

dir.create(file.path(out_dir, "matrix"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "deg"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "go"), recursive = TRUE, showWarnings = FALSE)

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
  nms <- names(dt)
  nms <- gsub("[()]", "", nms)
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
  if (!"log2FC" %in% names(cuff)) stop("差异表缺少 log2 fold change 列: ", path)
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
  stop("缺少 genes.read_group_tracking，无法从 Cuffdiff 原始文件重建表达矩阵。")
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

# Wide matrices: genes x samples
count_dt <- dcast(rg, tracking_id ~ sample, value.var = "raw_frags", fun.aggregate = sum)
fpkm_dt  <- dcast(rg, tracking_id ~ sample, value.var = "FPKM", fun.aggregate = mean)

count_dt <- merge(gene_map[, .(tracking_id, gene_short_name)], count_dt, by = "tracking_id", all.y = TRUE)
fpkm_dt  <- merge(gene_map[, .(tracking_id, gene_short_name)], fpkm_dt,  by = "tracking_id", all.y = TRUE)

write_tsv(count_dt, file.path(out_dir, "matrix", "counts_from_read_group_tracking.tsv"))
write_tsv(fpkm_dt,  file.path(out_dir, "matrix", "fpkm_from_read_group_tracking.tsv"))

coldata <- unique(rg[, .(sample, condition, replicate)])
coldata[, condition := factor(condition, levels = unique(c(control_group, treat_group, unique(condition))))]
write_tsv(coldata, file.path(out_dir, "matrix", "sample_metadata.tsv"))

n_per <- coldata[, .N, by = condition]
message("Replicates per group:\n", paste(capture.output(n_per), collapse = "\n"))

need_groups <- c(control_group, treat_group)
if (!all(need_groups %in% coldata$condition)) {
  stop("样本分组中找不到 ", paste(need_groups, collapse = " 或 "),
       "。实际分组: ", paste(unique(coldata$condition), collapse = ", "))
}

# ---------------------------------------------------------------------------
# Differential expression
# ---------------------------------------------------------------------------
can_deseq2 <- all(n_per[condition %in% need_groups, N] >= min_replicates)

deg <- NULL
deg_source <- NULL

if (can_deseq2) {
  need_pkg(c("DESeq2"), bioc = TRUE)
  suppressPackageStartupMessages(library(DESeq2))

  sample_cols <- setdiff(names(count_dt), c("tracking_id", "gene_short_name"))
  mat <- as.matrix(count_dt[, ..sample_cols])
  rownames(mat) <- count_dt$tracking_id
  mat[is.na(mat)] <- 0
  mat <- round(mat)
  storage.mode(mat) <- "integer"

  keep_samples <- coldata[condition %in% need_groups, sample]
  mat <- mat[, keep_samples, drop = FALSE]
  cd <- as.data.frame(coldata[match(keep_samples, sample)])
  rownames(cd) <- cd$sample
  cd$condition <- relevel(factor(cd$condition), ref = control_group)

  dds <- DESeqDataSetFromMatrix(countData = mat, colData = cd, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", treat_group, control_group), alpha = padj_cutoff)
  res_dt <- as.data.table(as.data.frame(res), keep.rownames = "tracking_id")
  res_dt <- merge(res_dt, gene_map[, .(tracking_id, gene_short_name)], by = "tracking_id", all.x = TRUE)
  setnames(res_dt,
           old = c("log2FoldChange", "pvalue", "padj", "gene_short_name"),
           new = c("log2FC", "p_value", "padj", "gene"),
           skip_absent = TRUE)
  deg <- res_dt
  deg_source <- "DESeq2_from_read_group_tracking"
  message("DESeq2 contrast: ", treat_group, " vs ", control_group)
} else {
  message("每组重复 < ", min_replicates, "，跳过 DESeq2，改用 Cuffdiff 差异表")
  if (!is.na(diff_path)) {
    deg <- read_cuffdiff_deg(diff_path, control_group, treat_group)
    deg_source <- "cuffdiff_gene_exp.diff"
  } else if (!is.na(excel_path)) {
    need_pkg("readxl")
    deg <- as.data.table(readxl::read_excel(excel_path))
    deg <- normalize_deg_names(deg)
    if ("status" %in% names(deg)) deg <- deg[status == "OK"]
    deg_source <- "excel_fallback"
    warning("未找到 gene_exp.diff，已回退 Excel: ", basename(excel_path))
  } else {
    stop("无法做差异分析：缺少足够重复做 DESeq2，且没有 gene_exp.diff / Excel。")
  }
}

# 无穷 log2FC（一组表达为 0）不用于火山图坐标，但仍可标显著
if ("log2FC" %in% names(deg)) {
  deg[, log2FC := as.numeric(log2FC)]
  deg[, log2FC_plot := log2FC]
  max_finite <- suppressWarnings(max(abs(deg$log2FC[is.finite(deg$log2FC)]), na.rm = TRUE))
  if (is.finite(max_finite)) {
    deg[!is.finite(log2FC_plot) & log2FC_plot > 0, log2FC_plot := max_finite]
    deg[!is.finite(log2FC_plot) & log2FC_plot < 0, log2FC_plot := -max_finite]
  }
}
if (!"padj" %in% names(deg) && "q_value" %in% names(deg)) deg[, padj := q_value]
if (!"p_value" %in% names(deg) && "pvalue" %in% names(deg)) deg[, p_value := pvalue]

deg[, direction := "NS"]
deg[!is.na(padj) & padj < padj_cutoff & !is.na(log2FC) & log2FC >= log2fc_cutoff, direction := "Up"]
deg[!is.na(padj) & padj < padj_cutoff & !is.na(log2FC) & log2FC <= -log2fc_cutoff, direction := "Down"]
deg[, deg_source := deg_source]

write_tsv(deg, file.path(out_dir, "deg", paste0("deg_", treat_group, "_vs_", control_group, ".tsv")))
sig <- deg[direction %in% c("Up", "Down")]
write_tsv(sig, file.path(out_dir, "deg", paste0("sig_deg_", treat_group, "_vs_", control_group, ".tsv")))
message("Significant genes: ", nrow(sig), " (Up ", sum(sig$direction == "Up"),
        ", Down ", sum(sig$direction == "Down"), ") via ", deg_source)

# ---------------------------------------------------------------------------
# Volcano
# ---------------------------------------------------------------------------
need_pkg(c("ggplot2", "ggrepel"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

volc <- deg[is.finite(log2FC_plot) & !is.na(padj)]
volc[, neglog10_padj := -log10(pmax(padj, .Machine$double.xmin))]
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
    label_dt[, label := vapply(gene, first_symbol, character(1))]
  } else {
    label_dt[, label := tracking_id]
  }
  label_dt <- label_dt[!is.na(label)]
}

p_volc <- ggplot(volc, aes(x = log2FC_plot, y = neglog10_padj, color = direction)) +
  geom_point(alpha = 0.7, size = 1.2) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = 2, linewidth = 0.3) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = 2, linewidth = 0.3) +
  scale_color_manual(values = c(Up = "#D62728", Down = "#1F77B4", NS = "grey70")) +
  labs(
    title = paste0(treat_group, " vs ", control_group, " volcano (", deg_source, ")"),
    x = "log2 fold change",
    y = "-log10 adjusted p-value"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.title = element_blank())

if (nrow(label_dt)) {
  p_volc <- p_volc + geom_text_repel(
    data = label_dt,
    aes(x = log2FC_plot, y = -log10(pmax(padj, .Machine$double.xmin)), label = label),
    size = 3, max.overlaps = 20, show.legend = FALSE, color = "black"
  )
}

volc_pdf <- file.path(out_dir, "plots", paste0("volcano_", treat_group, "_vs_", control_group, ".pdf"))
ggsave(volc_pdf, p_volc, width = 7, height = 6)
message("Wrote ", volc_pdf)

# ---------------------------------------------------------------------------
# GO enrichment
# ---------------------------------------------------------------------------
need_pkg(c("clusterProfiler", "enrichplot", species_orgdb), bioc = TRUE)
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
})

gene_col <- if ("gene" %in% names(sig)) sig$gene else sig$tracking_id
uni_col  <- if ("gene" %in% names(deg)) deg$gene else deg$tracking_id
sig_symbols <- unique(split_cuff_genes(gene_col))
universe_symbols <- unique(split_cuff_genes(uni_col))
if (!length(sig_symbols)) {
  warning("没有显著基因，跳过 GO 分析。可放宽 padj / |log2FC| 阈值。")
  quit(save = "no", status = 0)
}

map_ids <- function(symbols) {
  suppressMessages(bitr(
    symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = species_orgdb
  ))
}

sig_map <- map_ids(sig_symbols)
uni_map <- map_ids(universe_symbols)
message("Mapped ", nrow(sig_map), "/", length(sig_symbols), " significant symbols to Entrez")

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

for (ont in go_ontologies) {
  ego <- run_go(sig_map$ENTREZID, ont, uni_map$ENTREZID)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message("GO ", ont, ": no enriched terms")
    next
  }
  res_go <- as.data.table(as.data.frame(ego))
  write_tsv(res_go, file.path(out_dir, "go", paste0("GO_", ont, "_", treat_group, "_vs_", control_group, ".tsv")))

  pdf(file.path(out_dir, "plots", paste0("GO_", ont, "_dotplot.pdf")), width = 8, height = 6)
  print(dotplot(ego, showCategory = 15) + ggtitle(paste("GO", ont, treat_group, "vs", control_group)))
  dev.off()

  pdf(file.path(out_dir, "plots", paste0("GO_", ont, "_barplot.pdf")), width = 8, height = 6)
  print(barplot(ego, showCategory = 15) + ggtitle(paste("GO", ont, treat_group, "vs", control_group)))
  dev.off()
}

# Up / Down 分开的 GO（便于解释敲低效应）
for (dirn in c("Up", "Down")) {
  genes_dir <- unique(split_cuff_genes(sig[direction == dirn, gene]))
  if (length(genes_dir) < 5) {
    message("Skip GO for ", dirn, " (n=", length(genes_dir), ")")
    next
  }
  mapped <- map_ids(genes_dir)
  ego_bp <- run_go(mapped$ENTREZID, "BP", uni_map$ENTREZID)
  if (is.null(ego_bp) || nrow(as.data.frame(ego_bp)) == 0) next
  write_tsv(
    as.data.table(as.data.frame(ego_bp)),
    file.path(out_dir, "go", paste0("GO_BP_", dirn, "_", treat_group, "_vs_", control_group, ".tsv"))
  )
}

message("Done. Results in ", normalizePath(out_dir, mustWork = FALSE))
