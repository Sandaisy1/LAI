# Shared KEGG enrichment (human). Load clusterProfiler, enrichplot, org.Hs.eg.db, ggplot2 first.

save_kegg_enrichment <- function(entrez, universe, dest_dir, tag, title_prefix = tag,
                                 organism = "hsa") {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  gene_ids <- unique(as.character(entrez))
  uni_ids <- unique(as.character(universe))
  gene_ids <- gene_ids[!is.na(gene_ids) & nzchar(gene_ids)]
  uni_ids <- uni_ids[!is.na(uni_ids) & nzchar(uni_ids)]
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

  ek <- tryCatch(
    clusterProfiler::setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
    error = function(e) ek
  )
  out_tsv <- file.path(dest_dir, paste0("KEGG_", tag, ".tsv"))
  data.table::fwrite(data.table::as.data.table(as.data.frame(ek)), out_tsv, sep = "\t")
  message("Wrote ", out_tsv)

  pdf(file.path(dest_dir, paste0("KEGG_", tag, "_dotplot.pdf")), width = 8, height = 6)
  print(enrichplot::dotplot(ek, showCategory = 15) + ggplot2::ggtitle(paste(title_prefix, "KEGG")))
  dev.off()
  pdf(file.path(dest_dir, paste0("KEGG_", tag, "_barplot.pdf")), width = 8, height = 6)
  print(barplot(ek, showCategory = 15) + ggplot2::ggtitle(paste(title_prefix, "KEGG")))
  dev.off()
  invisible(ek)
}
