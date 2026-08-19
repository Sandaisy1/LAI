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

# ---------------------------------------------------------------------------
# GSEA on a full ranked list (not a truncated top-N / FC subset)
# ---------------------------------------------------------------------------
make_gsea_genelist <- function(dt, symbol_col = NULL, score_col = NULL) {
  dt <- data.table::as.data.table(dt)
  if (is.null(symbol_col)) {
    symbol_col <- if ("gene" %in% names(dt)) "gene" else "tracking_id"
  }
  if (is.null(score_col)) {
    if ("stat" %in% names(dt)) {
      score_col <- "stat"
    } else if (all(c("log2FC", "p_value") %in% names(dt))) {
      dt[, `__gsea_score` := sign(as.numeric(log2FC)) *
           -log10(pmax(as.numeric(p_value), .Machine$double.xmin))]
      score_col <- "__gsea_score"
    } else if ("log2FC" %in% names(dt)) {
      score_col <- "log2FC"
    } else if ("mean_log2FC" %in% names(dt)) {
      score_col <- "mean_log2FC"
    } else {
      stop("GSEA 需要 stat、log2FC 或 mean_log2FC 作为排序指标")
    }
  }
  tmp <- dt[, .(
    symbol = unlist(strsplit(as.character(get(symbol_col)), ",\\s*")),
    score = as.numeric(get(score_col))
  )]
  tmp[, symbol := trimws(symbol)]
  tmp <- tmp[!is.na(symbol) & !symbol %in% c("", "-", ".", "NA") & is.finite(score)]
  mapped <- tryCatch(
    clusterProfiler::bitr(unique(tmp$symbol), fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = "org.Hs.eg.db"),
    error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
  )
  if (!nrow(mapped)) return(numeric(0))
  tmp <- merge(tmp, mapped, by.x = "symbol", by.y = "SYMBOL")
  tmp <- tmp[order(-abs(score))]
  tmp <- tmp[!duplicated(ENTREZID)]
  gl <- tmp$score
  names(gl) <- as.character(tmp$ENTREZID)
  sort(gl, decreasing = TRUE)
}

.save_gsea_result <- function(res, dest_dir, file_prefix, plot_title) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    message(file_prefix, ": no GSEA terms")
    return(invisible(NULL))
  }
  res <- tryCatch(
    clusterProfiler::setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
    error = function(e) res
  )
  out_tsv <- file.path(dest_dir, paste0(file_prefix, ".tsv"))
  data.table::fwrite(data.table::as.data.table(as.data.frame(res)), out_tsv, sep = "\t")
  message("Wrote ", out_tsv)
  pdf(file.path(dest_dir, paste0(file_prefix, "_dotplot.pdf")), width = 8, height = 6)
  p_dot <- tryCatch(
    enrichplot::dotplot(res, showCategory = 15, split = ".sign") +
      ggplot2::facet_grid(. ~ .sign) + ggplot2::ggtitle(plot_title),
    error = function(e) enrichplot::dotplot(res, showCategory = 15) + ggplot2::ggtitle(plot_title)
  )
  print(p_dot)
  dev.off()
  tryCatch({
    pdf(file.path(dest_dir, paste0(file_prefix, "_ridgeplot.pdf")), width = 8, height = 7)
    print(enrichplot::ridgeplot(res, showCategory = 15) + ggplot2::ggtitle(plot_title))
    dev.off()
  }, error = function(e) invisible(NULL))
  n_term <- min(3L, nrow(as.data.frame(res)))
  if (n_term >= 1) {
    tryCatch({
      pdf(file.path(dest_dir, paste0(file_prefix, "_gseaplot.pdf")), width = 8, height = 8)
      print(enrichplot::gseaplot2(res, geneSetID = seq_len(n_term), pvalue_table = TRUE, title = plot_title))
      dev.off()
    }, error = function(e) invisible(NULL))
  }
  invisible(res)
}

save_gsea_analyses <- function(deg_dt, dest_dir, tag, title_prefix = tag) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  geneList <- make_gsea_genelist(deg_dt)
  if (length(geneList) < 50) {
    message(tag, " skip GSEA (ranked Entrez list too short: ", length(geneList), ")")
    return(invisible(NULL))
  }
  data.table::fwrite(
    data.table::data.table(ENTREZID = names(geneList), rank_metric = unname(geneList)),
    file.path(dest_dir, paste0("GSEA_ranked_genelist_", tag, ".tsv")),
    sep = "\t"
  )
  message(tag, " GSEA ranked genes: ", length(geneList))
  set.seed(1)

  for (ont in c("BP", "MF", "CC")) {
    res <- tryCatch({
      clusterProfiler::gseGO(
        geneList = geneList,
        ont = ont,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = FALSE,
        eps = 0
      )
    }, error = function(e) {
      warning(tag, " gseGO ", ont, " failed: ", conditionMessage(e))
      NULL
    })
    .save_gsea_result(res, dest_dir, paste0("GSEA_GO_", ont, "_", tag), paste(title_prefix, "GSEA GO", ont))
  }

  ek <- tryCatch({
    clusterProfiler::gseKEGG(
      geneList = geneList,
      organism = "hsa",
      keyType = "ncbi-geneid",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE,
      eps = 0
    )
  }, error = function(e) {
    warning(tag, " gseKEGG failed: ", conditionMessage(e))
    NULL
  })
  .save_gsea_result(ek, dest_dir, paste0("GSEA_KEGG_", tag), paste(title_prefix, "GSEA KEGG"))

  .ensure_pkg("ReactomePA")
  er <- tryCatch({
    ReactomePA::gsePathway(
      geneList = geneList,
      organism = "human",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE,
      eps = 0
    )
  }, error = function(e) {
    warning(tag, " gsePathway failed: ", conditionMessage(e))
    NULL
  })
  .save_gsea_result(er, dest_dir, paste0("GSEA_Reactome_", tag), paste(title_prefix, "GSEA Reactome"))

  ew <- tryCatch({
    clusterProfiler::gseWP(
      geneList = geneList,
      organism = "Homo sapiens",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      eps = 0
    )
  }, error = function(e) {
    warning(tag, " gseWP failed: ", conditionMessage(e))
    NULL
  })
  .save_gsea_result(ew, dest_dir, paste0("GSEA_WikiPathways_", tag), paste(title_prefix, "GSEA WikiPathways"))
  invisible(geneList)
}
