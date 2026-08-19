#!/usr/bin/env Rscript
# Additional GO analysis (does NOT modify scripts/tg_vs_ntc_deg_go.R):
#   For each of the three contrasts, take the top 50 / 100 / 150 / 200 / 250 / 300
#   upregulated genes (ranked by log2FC descending) and run GO + KEGG + Reactome/WikiPathways separately.
#
# 用法（先跑完原始 DEG 脚本）:
#   Rscript scripts/tg_vs_ntc_deg_go.R "E:/R/TG_BRCA/TG" "E:/R/TG_BRCA/TG/results"
#   Rscript scripts/tg_up_topn_go.R "E:/R/TG_BRCA/TG" "E:/R/TG_BRCA/TG/results"

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else "E:/R/TG_BRCA/TG"
out_dir  <- if (length(args) >= 2) args[[2]] else file.path(data_dir, "results")

contrasts <- c(
  TG_sh1_vs_NTC = file.path(out_dir, "pairwise", "TG_sh1_vs_NTC", "deg", "deg_TG_sh1_vs_NTC.tsv"),
  TG_sh5_vs_NTC = file.path(out_dir, "pairwise", "TG_sh5_vs_NTC", "deg", "deg_TG_sh5_vs_NTC.tsv"),
  TG_sh1_sh5_avg_vs_NTC = file.path(
    out_dir, "pooled_avg", "TG_sh1_sh5_avg_vs_NTC", "deg", "deg_TG_sh1_sh5_avg_vs_NTC.tsv"
  )
)
top_ns <- c(50L, 100L, 150L, 200L, 250L, 300L)
go_ontologies <- c("BP", "MF", "CC")
species_orgdb <- "org.Hs.eg.db"
topn_root <- file.path(out_dir, "go_topn")
dir.create(topn_root, recursive = TRUE, showWarnings = FALSE)

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

find_deg_fallback <- function(out_dir, contrast_id) {
  hits <- list.files(
    out_dir, pattern = paste0("^deg_", contrast_id, "\\.tsv$"),
    recursive = TRUE, full.names = TRUE
  )
  hits <- hits[!grepl("sig_deg_", basename(hits))]
  if (!length(hits)) return(NA_character_)
  hits[[1]]
}

rank_up_genes <- function(deg) {
  deg <- copy(deg)
  deg[, log2FC := as.numeric(log2FC)]
  if (!"padj" %in% names(deg)) deg[, padj := NA_real_]
  if (!"fold_change" %in% names(deg)) deg[, fold_change := 2^log2FC]
  up <- deg[!is.na(log2FC) & is.finite(log2FC) & log2FC > 0]
  if (!nrow(up)) return(up)
  if (any(!is.na(up$padj))) {
    up <- up[order(-log2FC, padj)]
  } else {
    up <- up[order(-log2FC)]
  }
  up[, up_rank := seq_len(.N)]
  up
}

need_pkg("ggplot2")
need_pkg(c("clusterProfiler", "enrichplot", "ReactomePA", species_orgdb), bioc = TRUE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
})

map_ids <- function(symbols) {
  if (!length(symbols)) {
    return(data.frame(SYMBOL = character(), ENTREZID = character()))
  }
  suppressMessages(bitr(
    symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = species_orgdb
  ))
}

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
source(kegg_helper, local = FALSE)

summary_rows <- list()

for (contrast_id in names(contrasts)) {
  deg_path <- contrasts[[contrast_id]]
  if (!file.exists(deg_path)) {
    alt <- find_deg_fallback(out_dir, contrast_id)
    if (is.na(alt)) {
      warning("找不到 ", contrast_id, " 的 DEG 表。请先运行 scripts/tg_vs_ntc_deg_go.R。跳过。")
      next
    }
    deg_path <- alt
  }
  message("==== ", contrast_id, " from ", deg_path, " ====")
  deg <- normalize_deg_names(fread(deg_path))
  if (!"log2FC" %in% names(deg)) {
    warning(contrast_id, " 缺少 log2FC，跳过")
    next
  }

  up_all <- rank_up_genes(deg)
  uni_col <- if ("gene" %in% names(deg)) deg$gene else deg$tracking_id
  uni_map <- map_ids(unique(split_cuff_genes(uni_col)))
  contrast_dir <- file.path(topn_root, contrast_id)
  dir.create(contrast_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(up_all, file.path(contrast_dir, paste0("all_upregulated_ranked_", contrast_id, ".tsv")))

  for (n in top_ns) {
    tag <- paste0("top", n)
    tag_dir <- file.path(contrast_dir, tag)
    dir.create(tag_dir, recursive = TRUE, showWarnings = FALSE)
    n_use <- min(n, nrow(up_all))
    top <- if (n_use) up_all[seq_len(n_use)] else up_all[0]
    if (n_use < n) {
      warning(contrast_id, " ", tag, ": 上调基因只有 ", n_use, " 个，不足 ", n)
    }
    write_tsv(top, file.path(tag_dir, paste0("up_genes_", tag, "_", contrast_id, ".tsv")))

    gene_col <- if ("gene" %in% names(top)) top$gene else top$tracking_id
    symbols <- unique(split_cuff_genes(gene_col))
    summary_rows[[paste(contrast_id, tag)]] <- data.table(
      contrast = contrast_id,
      top_n_requested = n,
      n_up_rows = n_use,
      n_symbols = length(symbols)
    )
    message(contrast_id, " ", tag, ": ", n_use, " rows / ", length(symbols), " symbols")

    if (length(symbols) < 5) {
      message("skip GO/KEGG/pathway for ", contrast_id, " ", tag, " (need >=5 unique symbols)")
      next
    }
    sig_map <- map_ids(symbols)
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
}

if (length(summary_rows)) {
  write_tsv(rbindlist(summary_rows), file.path(topn_root, "up_topn_gene_counts.tsv"))
}

message("Done. Top-N GO in ", normalizePath(topn_root, mustWork = FALSE))
