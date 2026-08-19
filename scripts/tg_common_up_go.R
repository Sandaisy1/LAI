#!/usr/bin/env Rscript
# 共同上调基因 GO（独立脚本，不修改 tg_vs_ntc_deg_go.R）
#   1) 单独比较 TG_sh1 vs NTC、TG_sh5 vs NTC
#   2) 取两组都上调的基因
#   3a) 线性 FC>=1 / 1.25 / 1.5 / 2（两侧都达标）各自 GO + KEGG + Reactome/WikiPathways
#   3b) 各组上调前 50/75/100/150/200/250/300 的交集，各自 GO + KEGG + Reactome/WikiPathways
#
# 用法:
#   Rscript scripts/tg_common_up_go.R [data_dir] [out_dir]
# 默认: E:/R/TG_BRCA/TG

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else "E:/R/TG_BRCA/TG"
out_dir  <- if (length(args) >= 2) args[[2]] else file.path(data_dir, "results")

control_group <- "NTC"
treat_groups  <- c("TG_sh1", "TG_sh5")
padj_cutoff   <- 0.05
min_replicates <- 2
species_orgdb <- "org.Hs.eg.db"
go_ontologies <- c("BP", "MF", "CC")
go_fc_cutoffs <- c("up_FC1" = 1, "up_FC1.25" = 1.25, "up_FC1.5" = 1.5, "up_FC2" = 2)
top_ns <- c(50L, 75L, 100L, 150L, 200L, 250L, 300L)

root <- file.path(out_dir, "common_up")
dir.create(file.path(root, "pairwise"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "go_fc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "go_topn"), recursive = TRUE, showWarnings = FALSE)

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

write_tsv <- function(x, path) {
  fwrite(x, path, sep = "\t")
  message("Wrote ", path)
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
  if (!"padj" %in% names(dt)) dt[, padj := NA_real_]
  dt[, log2FC := as.numeric(log2FC)]
  dt
}

explode_symbols <- function(deg) {
  dt <- copy(deg)
  if (!"gene" %in% names(dt)) dt[, gene := tracking_id]
  dt[, gene := as.character(gene)]
  dt[, rid := .I]
  exploded <- dt[, .(gene = unlist(strsplit(gene, ",\\s*"))), by = rid]
  exploded[, gene := trimws(gene)]
  keep <- setdiff(names(dt), "gene")
  out <- merge(exploded, dt[, ..keep], by = "rid", all.x = TRUE)
  out <- out[!is.na(gene) & !gene %in% c("", "-", ".", "NA")]
  out[, log2FC := as.numeric(log2FC)]
  out <- out[is.finite(log2FC)]
  out[order(-log2FC)][, .SD[1], by = gene]
}

rank_up <- function(sym_dt) {
  up <- sym_dt[log2FC > 0]
  if (any(!is.na(up$padj))) up <- up[order(-log2FC, padj)] else up <- up[order(-log2FC)]
  up[, up_rank := seq_len(.N)]
  up
}

# ---------------------------------------------------------------------------
# Pairwise DEG: reuse existing tables, otherwise compute from tracking
# ---------------------------------------------------------------------------
existing_paths <- c(
  TG_sh1 = file.path(out_dir, "pairwise", "TG_sh1_vs_NTC", "deg", "deg_TG_sh1_vs_NTC.tsv"),
  TG_sh5 = file.path(out_dir, "pairwise", "TG_sh5_vs_NTC", "deg", "deg_TG_sh5_vs_NTC.tsv")
)

load_existing <- function(treat) {
  p <- existing_paths[[treat]]
  if (file.exists(p)) return(normalize_deg_names(fread(p)))
  alt <- list.files(out_dir, pattern = paste0("^deg_", treat, "_vs_NTC\\.tsv$"),
                    recursive = TRUE, full.names = TRUE)
  alt <- alt[!grepl("sig_deg_", basename(alt))]
  if (length(alt)) return(normalize_deg_names(fread(alt[[1]])))
  NULL
}

deg_list <- lapply(treat_groups, load_existing)
names(deg_list) <- treat_groups

need_compute <- any(vapply(deg_list, is.null, logical(1)))
if (need_compute) {
  message("未找到完整 pairwise DEG，从 genes.read_group_tracking 重新计算 TG_sh1/TG_sh5 vs NTC")
  tracking_path <- find_input(data_dir, "^genes\\.read_group_tracking$")
  fpkm_map_path <- find_input(data_dir, "^genes\\.fpkm_tracking$", required = FALSE)
  diff_path <- find_input(data_dir, "^gene_exp\\.diff$", required = FALSE)

  rg <- fread(tracking_path)
  req_cols <- c("tracking_id", "condition", "replicate", "raw_frags", "FPKM")
  if (length(setdiff(req_cols, names(rg)))) {
    stop("genes.read_group_tracking 缺少列: ", paste(setdiff(req_cols, names(rg)), collapse = ", "))
  }
  rg[, `:=`(
    condition = as.character(condition),
    replicate = as.character(replicate),
    sample = paste(condition, replicate, sep = "_"),
    raw_frags = as.numeric(raw_frags),
    FPKM = as.numeric(FPKM)
  )]

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
    if (!"gene_short_name" %in% names(gene_map)) gene_map[, gene_short_name := tracking_id]
  }
  gene_map[, gene_short_name := sanitize_gene_symbols(gene_short_name)]

  count_dt <- dcast(rg, tracking_id ~ sample, value.var = "raw_frags", fun.aggregate = sum)
  fpkm_dt  <- dcast(rg, tracking_id ~ sample, value.var = "FPKM", fun.aggregate = mean)
  count_dt <- merge(gene_map[, .(tracking_id, gene_short_name)], count_dt, by = "tracking_id", all.y = TRUE)
  fpkm_dt  <- merge(gene_map[, .(tracking_id, gene_short_name)], fpkm_dt, by = "tracking_id", all.y = TRUE)
  coldata <- unique(rg[, .(sample, condition, replicate)])
  n_per <- coldata[, .N, by = condition]
  n_tbl <- setNames(n_per$N, as.character(n_per$condition))

  sample_cols <- setdiff(names(count_dt), c("tracking_id", "gene_short_name"))
  count_mat <- as.matrix(count_dt[, ..sample_cols])
  rownames(count_mat) <- count_dt$tracking_id
  count_mat[is.na(count_mat)] <- 0
  count_mat <- round(count_mat)
  storage.mode(count_mat) <- "integer"
  fpkm_cols <- setdiff(names(fpkm_dt), c("tracking_id", "gene_short_name"))
  fpkm_mat <- as.matrix(fpkm_dt[, ..fpkm_cols])
  rownames(fpkm_mat) <- fpkm_dt$tracking_id
  coldata_df <- as.data.frame(coldata)

  run_one_deg <- function(treat) {
    n_ctrl <- n_tbl[[control_group]]
    n_trt <- n_tbl[[treat]]
    if (isTRUE(n_ctrl >= min_replicates && n_trt >= min_replicates)) {
      need_pkg("DESeq2", bioc = TRUE)
      suppressPackageStartupMessages(library(DESeq2))
      keep <- coldata_df$condition %in% c(control_group, treat)
      mat <- count_mat[, coldata_df$sample[keep], drop = FALSE]
      cd <- droplevels(coldata_df[keep, , drop = FALSE])
      rownames(cd) <- cd$sample
      cd$condition <- relevel(factor(cd$condition), ref = control_group)
      dds <- DESeqDataSetFromMatrix(countData = mat, colData = cd, design = ~ condition)
      dds <- dds[rowSums(counts(dds)) >= 10, ]
      dds <- DESeq(dds)
      res <- results(dds, contrast = c("condition", treat, control_group), alpha = padj_cutoff)
      res_dt <- as.data.table(as.data.frame(res), keep.rownames = "tracking_id")
      res_dt <- merge(res_dt, gene_map[, .(tracking_id, gene_short_name)], by = "tracking_id", all.x = TRUE)
      setnames(res_dt,
               old = c("log2FoldChange", "pvalue", "padj", "gene_short_name"),
               new = c("log2FC", "p_value", "padj", "gene"),
               skip_absent = TRUE)
      return(normalize_deg_names(res_dt))
    }
    if (!is.na(diff_path) && file.exists(diff_path)) {
      cuff <- normalize_deg_names(fread(diff_path))
      if ("status" %in% names(cuff)) cuff <- cuff[status == "OK"]
      if (all(c("sample_1", "sample_2") %in% names(cuff))) {
        cuff <- cuff[
          (sample_1 == control_group & sample_2 == treat) |
            (sample_1 == treat & sample_2 == control_group)
        ]
        cuff[sample_1 == treat & sample_2 == control_group, log2FC := -as.numeric(log2FC)]
      }
      if (nrow(cuff) && "log2FC" %in% names(cuff)) return(cuff)
    }
    ctrl_s <- coldata[condition == control_group, sample]
    trt_s <- coldata[condition == treat, sample]
    a <- fpkm_mat[, ctrl_s, drop = FALSE]
    b <- fpkm_mat[, trt_s, drop = FALSE]
    log2fc <- log2((rowMeans(b, na.rm = TRUE) + 0.1) / (rowMeans(a, na.rm = TRUE) + 0.1))
    normalize_deg_names(data.table(
      tracking_id = rownames(fpkm_mat),
      gene = gene_map$gene_short_name[match(rownames(fpkm_mat), gene_map$tracking_id)],
      log2FC = log2fc,
      padj = NA_real_
    ))
  }

  for (treat in treat_groups) {
    if (is.null(deg_list[[treat]])) {
      if (!treat %in% as.character(coldata$condition)) {
        stop("tracking 中没有 ", treat, "，无法与 NTC 比较")
      }
      deg_list[[treat]] <- run_one_deg(treat)
    }
  }
}

# 写出两组单独比较结果
need_pkg(c("ggplot2", "ggrepel"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

for (treat in treat_groups) {
  deg <- copy(deg_list[[treat]])
  contrast_id <- paste0(treat, "_vs_", control_group)
  dest <- file.path(root, "pairwise", contrast_id)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  deg[, fold_change := 2^log2FC]
  deg[, direction := "NS"]
  deg[!is.na(padj) & padj < padj_cutoff & log2FC > 0, direction := "Up"]
  deg[!is.na(padj) & padj < padj_cutoff & log2FC < 0, direction := "Down"]
  deg[is.na(padj) & log2FC > 0, direction := "Up"]
  deg[is.na(padj) & log2FC < 0, direction := "Down"]
  write_tsv(deg, file.path(dest, paste0("deg_", contrast_id, ".tsv")))
  write_tsv(deg[direction == "Up"], file.path(dest, paste0("up_genes_", contrast_id, ".tsv")))

  volc <- deg[is.finite(log2FC)]
  volc[, neglog10_padj := -log10(pmax(fifelse(is.na(padj), 1, padj), .Machine$double.xmin))]
  p <- ggplot(volc, aes(x = log2FC, y = neglog10_padj, color = direction)) +
    geom_point(alpha = 0.7, size = 1.1) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = 2, linewidth = 0.3) +
    scale_color_manual(values = c(Up = "#D62728", Down = "#1F77B4", NS = "grey70")) +
    labs(title = paste(treat, "vs", control_group), x = "log2 fold change", y = "-log10 adjusted p") +
    theme_bw(base_size = 12) +
    theme(legend.title = element_blank())
  ggsave(file.path(dest, paste0("volcano_", contrast_id, ".pdf")), p, width = 7, height = 6)
  message("Pairwise done: ", contrast_id)
  deg_list[[treat]] <- deg
}

# ---------------------------------------------------------------------------
# Common upregulated genes
# ---------------------------------------------------------------------------
sym_sh1 <- explode_symbols(deg_list[["TG_sh1"]])
sym_sh5 <- explode_symbols(deg_list[["TG_sh5"]])
up_sh1 <- rank_up(sym_sh1)
up_sh5 <- rank_up(sym_sh5)

sh1_m <- copy(sym_sh1)
sh5_m <- copy(sym_sh5)
setnames(sh1_m, old = c("log2FC", "padj"), new = c("log2FC_sh1", "padj_sh1"), skip_absent = TRUE)
setnames(sh5_m, old = c("log2FC", "padj"), new = c("log2FC_sh5", "padj_sh5"), skip_absent = TRUE)

common <- merge(
  sh1_m[, .(gene, log2FC_sh1, padj_sh1, fold_change_sh1 = 2^log2FC_sh1)],
  sh5_m[, .(gene, log2FC_sh5, padj_sh5, fold_change_sh5 = 2^log2FC_sh5)],
  by = "gene"
)
common[, `:=`(
  min_FC = pmin(fold_change_sh1, fold_change_sh5),
  mean_log2FC = (log2FC_sh1 + log2FC_sh5) / 2
)]
common_up <- common[log2FC_sh1 > 0 & log2FC_sh5 > 0]
if (any(!is.na(common_up$padj_sh1)) && any(!is.na(common_up$padj_sh5))) {
  common_up_sig <- common_up[padj_sh1 < padj_cutoff & padj_sh5 < padj_cutoff]
  if (nrow(common_up_sig)) common_up <- common_up_sig
}
common_up <- common_up[order(-min_FC, -mean_log2FC)]
write_tsv(common, file.path(root, "all_genes_both_contrasts.tsv"))
write_tsv(common_up, file.path(root, "common_upregulated_genes.tsv"))
message("Common upregulated genes: ", nrow(common_up))

# ---------------------------------------------------------------------------
# GO helpers
# ---------------------------------------------------------------------------
need_pkg("ggplot2")
need_pkg(c("clusterProfiler", "enrichplot", "ReactomePA", species_orgdb), bioc = TRUE)
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
})

universe_symbols <- unique(c(sym_sh1$gene, sym_sh5$gene))
map_ids <- function(symbols) {
  if (!length(symbols)) return(data.frame(SYMBOL = character(), ENTREZID = character()))
  suppressMessages(bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = species_orgdb))
}
uni_map <- map_ids(universe_symbols)

run_go <- function(entrez, ont) {
  enrichGO(
    gene = unique(entrez),
    universe = unique(uni_map$ENTREZID),
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
source(kegg_helper, local = FALSE)

go_one_set <- function(symbols, dest_dir, tag, gene_dt = NULL) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  symbols <- unique(as.character(symbols))
  if (is.null(gene_dt)) {
    gene_dt <- common_up[gene %in% symbols][order(-min_FC)]
  } else {
    gene_dt <- gene_dt[gene %in% symbols]
  }
  write_tsv(gene_dt, file.path(dest_dir, paste0("genes_", tag, ".tsv")))
  message(tag, ": ", length(symbols), " common up symbols")
  if (length(symbols) < 5) {
    message("skip GO/KEGG/pathway for ", tag, " (need >=5 genes)")
    return(invisible(NULL))
  }
  mapped <- map_ids(symbols)
  if (!nrow(mapped)) return(invisible(NULL))
  for (ont in go_ontologies) {
    ego <- run_go(mapped$ENTREZID, ont)
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
      message(tag, " GO ", ont, ": no terms")
      next
    }
    write_tsv(as.data.table(as.data.frame(ego)), file.path(dest_dir, paste0("GO_", ont, "_", tag, ".tsv")))
    pdf(file.path(dest_dir, paste0("GO_", ont, "_", tag, "_dotplot.pdf")), width = 8, height = 6)
    print(dotplot(ego, showCategory = 15) + ggtitle(paste("common up", tag, ont)))
    dev.off()
    pdf(file.path(dest_dir, paste0("GO_", ont, "_", tag, "_barplot.pdf")), width = 8, height = 6)
    print(barplot(ego, showCategory = 15) + ggtitle(paste("common up", tag, ont)))
    dev.off()
  }
  save_kegg_and_pathways(
    entrez = mapped$ENTREZID,
    universe = uni_map$ENTREZID,
    dest_dir = dest_dir,
    tag = tag,
    title_prefix = paste("common up", tag)
  )
}

# 1) FoldChange：两侧线性 FC 都 >= 阈值
fc_counts <- list()
for (tag in names(go_fc_cutoffs)) {
  fc <- as.numeric(go_fc_cutoffs[[tag]])
  keep <- common_up[min_FC >= fc, gene]
  go_one_set(keep, file.path(root, "go_fc", tag), tag)
  fc_counts[[tag]] <- data.table(method = "FC", cutoff = tag, linear_FC_min = fc, n_genes = length(unique(keep)))
}
write_tsv(rbindlist(fc_counts), file.path(root, "go_fc", "common_up_FC_gene_counts.tsv"))

# 2) 排名：各组上调 top-N 的交集（不是把两组基因并在一起）
write_tsv(up_sh1, file.path(root, "pairwise", "TG_sh1_up_ranked_by_log2FC.tsv"))
write_tsv(up_sh5, file.path(root, "pairwise", "TG_sh5_up_ranked_by_log2FC.tsv"))

topn_counts <- list()
for (n in top_ns) {
  tag <- paste0("top", n)
  n1 <- min(n, nrow(up_sh1))
  n5 <- min(n, nrow(up_sh5))
  set1 <- up_sh1[seq_len(n1), gene]
  set5 <- up_sh5[seq_len(n5), gene]
  inter <- intersect(set1, set5)
  inter_dt <- merge(
    up_sh1[gene %in% inter, .(gene, log2FC_sh1 = log2FC, rank_sh1 = up_rank, padj_sh1 = padj)],
    up_sh5[gene %in% inter, .(gene, log2FC_sh5 = log2FC, rank_sh5 = up_rank, padj_sh5 = padj)],
    by = "gene"
  )
  inter_dt[, `:=`(min_FC = pmin(2^log2FC_sh1, 2^log2FC_sh5), mean_log2FC = (log2FC_sh1 + log2FC_sh5) / 2)]
  inter_dt <- inter_dt[order(-min_FC, -mean_log2FC)]
  if (n1 < n || n5 < n) {
    warning(tag, ": 单组上调基因不足 ", n, "（sh1=", n1, ", sh5=", n5, "）")
  }
  go_one_set(inter, file.path(root, "go_topn", tag), tag, gene_dt = inter_dt)
  topn_counts[[tag]] <- data.table(
    method = "topN_intersection",
    cutoff = tag,
    n_requested = n,
    n_sh1 = n1,
    n_sh5 = n5,
    n_common = length(inter)
  )
}
write_tsv(rbindlist(topn_counts), file.path(root, "go_topn", "common_up_topn_gene_counts.tsv"))

message("Done. Common-up GO in ", normalizePath(root, mustWork = FALSE))
