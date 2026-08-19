# Shared KEGG / Reactome / WikiPathways enrichment (human).
# Load clusterProfiler, enrichplot, org.Hs.eg.db, ggplot2 first.

.ensure_pkg <- function(pkg, bioc = TRUE) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  invisible(TRUE)
}

.clean_ids <- function(x) {
  x <- unique(as.character(x))
  x[!is.na(x) & nzchar(x)]
}

.save_enrich_result <- function(res, dest_dir, file_prefix, plot_title) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(invisible(NULL))
  res <- tryCatch(
    clusterProfiler::setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
    error = function(e) res
  )
  out_tsv <- file.path(dest_dir, paste0(file_prefix, ".tsv"))
  data.table::fwrite(data.table::as.data.table(as.data.frame(res)), out_tsv, sep = "\t")
  message("Wrote ", out_tsv)
  pdf(file.path(dest_dir, paste0(file_prefix, "_dotplot.pdf")), width = 8, height = 6)
  print(enrichplot::dotplot(res, showCategory = 15) + ggplot2::ggtitle(plot_title))
  dev.off()
  pdf(file.path(dest_dir, paste0(file_prefix, "_barplot.pdf")), width = 8, height = 6)
  print(barplot(res, showCategory = 15) + ggplot2::ggtitle(plot_title))
  dev.off()
  invisible(res)
}

save_kegg_enrichment <- function(entrez, universe, dest_dir, tag, title_prefix = tag,
                                 organism = "hsa") {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  gene_ids <- .clean_ids(entrez)
  uni_ids <- .clean_ids(universe)
  if (length(gene_ids) < 5) {
    message(tag, " skip KEGG (need >=5 Entrez IDs)")
    return(invisible(NULL))
  }
  ek <- tryCatch({
    clusterProfiler::enrichKEGG(
      gene = gene_ids,
      universe = if (length(uni_ids)) uni_ids else NULL,
      organism = organism,
      keyType = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  }, error = function(e) {
    warning(tag, " KEGG failed: ", conditionMessage(e))
    NULL
  })
  if (is.null(ek) || nrow(as.data.frame(ek)) == 0) {
    message(tag, " KEGG: no terms")
    return(invisible(NULL))
  }
  .save_enrich_result(ek, dest_dir, paste0("KEGG_", tag), paste(title_prefix, "KEGG"))
}

save_reactome_enrichment <- function(entrez, universe, dest_dir, tag, title_prefix = tag) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  gene_ids <- .clean_ids(entrez)
  uni_ids <- .clean_ids(universe)
  if (length(gene_ids) < 5) {
    message(tag, " skip Reactome (need >=5 Entrez IDs)")
    return(invisible(NULL))
  }
  .ensure_pkg("ReactomePA")
  er <- tryCatch({
    ReactomePA::enrichPathway(
      gene = gene_ids,
      universe = if (length(uni_ids)) uni_ids else NULL,
      organism = "human",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
  }, error = function(e) {
    warning(tag, " Reactome failed: ", conditionMessage(e))
    NULL
  })
  if (is.null(er) || nrow(as.data.frame(er)) == 0) {
    message(tag, " Reactome: no terms")
    return(invisible(NULL))
  }
  .save_enrich_result(er, dest_dir, paste0("Reactome_", tag), paste(title_prefix, "Reactome"))
}

save_wikipathways_enrichment <- function(entrez, universe, dest_dir, tag, title_prefix = tag) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  gene_ids <- .clean_ids(entrez)
  uni_ids <- .clean_ids(universe)
  if (length(gene_ids) < 5) {
    message(tag, " skip WikiPathways (need >=5 Entrez IDs)")
    return(invisible(NULL))
  }
  ew <- tryCatch({
    clusterProfiler::enrichWP(
      gene = gene_ids,
      universe = if (length(uni_ids)) uni_ids else NULL,
      organism = "Homo sapiens",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  }, error = function(e) {
    warning(tag, " WikiPathways failed: ", conditionMessage(e))
    NULL
  })
  if (is.null(ew) || nrow(as.data.frame(ew)) == 0) {
    message(tag, " WikiPathways: no terms")
    return(invisible(NULL))
  }
  .save_enrich_result(ew, dest_dir, paste0("WikiPathways_", tag), paste(title_prefix, "WikiPathways"))
}

# KEGG + Reactome + WikiPathways on the same gene set (one extra pathway layer beyond GO)
save_kegg_and_pathways <- function(entrez, universe, dest_dir, tag, title_prefix = tag) {
  save_kegg_enrichment(entrez, universe, dest_dir, tag, title_prefix)
  save_reactome_enrichment(entrez, universe, dest_dir, tag, title_prefix)
  save_wikipathways_enrichment(entrez, universe, dest_dir, tag, title_prefix)
  invisible(NULL)
}
