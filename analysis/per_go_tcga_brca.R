#!/usr/bin/env Rscript
# Per-GO (NOT pooled) TCGA-BRCA analysis:
#   1) each listed GO term vs clinical / survival
#   2) genes (and proteins) negatively correlated with each GO score
#
# Required files in --data-dir (prefix match; .tsv/.txt/.gz OK):
#   TCGA-BRCA.clinical
#   TCGA-BRCA.protein
#   TCGA-BRCA.star_fpkm
#   TCGA-BRCA.survival
#   gencode.v36.annotation.gtf.gene
#
# Packages:
#   install.packages(c("data.table", "survival"))
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("org.Hs.eg.db", "AnnotationDbi", "GSVA", "GO.db"))
#
# Usage:
#   Rscript analysis/per_go_tcga_brca.R
#   Rscript analysis/per_go_tcga_brca.R --data-dir /path/to/data --out-dir results

suppressPackageStartupMessages({
  .need <- c("data.table")
  .miss <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(.miss)) {
    stop("Please install: install.packages(c(", paste(sprintf('\"%s\"', .miss), collapse = ", "), "))")
  }
  library(data.table)
})

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) args[[i + 1L]] else default
}

DATA_DIR <- get_arg("--data-dir", "data")
OUT_DIR <- get_arg("--out-dir", "results")
FDR_CUTOFF <- as.numeric(get_arg("--fdr", "0.05"))
MIN_RHO <- as.numeric(get_arg("--min-abs-rho", "0"))  # extra |rho| filter; 0 = none
MIN_GENES_PER_GO <- 5L
TUMOR_SAMPLE_TYPE <- "01"

# Each GO is an independent hypothesis. Never union these gene sets.
GO_TERMS <- c(
  "GO:0023041",  # neuronal signal transduction
  "GO:1904457",  # positive regulation of neuronal action potential
  "GO:1904340",  # positive regulation of dopaminergic neuron differentiation
  "GO:2001224",  # positive regulation of neuron migration
  "GO:2001222",  # regulation of neuron migration
  "GO:0019227",  # neuronal action potential propagation
  "GO:0019228",  # neuronal action potential
  "GO:1902847",  # regulation of neuronal signal transduction
  "GO:0031102",  # neuron projection regeneration
  "GO:0097492",  # sympathetic neuron axon guidance
  "GO:0097491",  # sympathetic neuron projection guidance
  "GO:0097374",  # sensory neuron axon guidance
  "GO:0007158",  # neuron cell-cell adhesion
  "GO:1902667",  # regulation of axon guidance
  "GO:0031103",  # axon regeneration
  "GO:0007411",  # axon guidance
  "GO:0007409"   # axonogenesis
)

GO_SAFE <- function(go_id) gsub(":", "_", go_id, fixed = TRUE)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

find_data_file <- function(prefix, dir = DATA_DIR) {
  if (!dir.exists(dir)) {
    stop("Data directory not found: ", dir,
         "\nPlace TCGA-BRCA.* and gencode.* files there, or pass --data-dir")
  }
  hits <- list.files(dir, full.names = TRUE, recursive = TRUE)
  hits <- hits[grepl(prefix, basename(hits), fixed = TRUE)]
  hits <- hits[!grepl("\\.(md|R|r)$", hits)]
  if (!length(hits)) {
    stop("No file matching prefix '", prefix, "' under ", dir)
  }
  # Prefer uncompressed tsv/txt if both exist
  hits <- hits[order(nchar(basename(hits)), basename(hits))]
  hits[[1L]]
}

strip_ensembl_version <- function(x) sub("\\.[0-9]+$", "", x)

tcga_sample <- function(x) substr(gsub("\\.", "-", x), 1L, 15L)
tcga_patient <- function(x) substr(gsub("\\.", "-", x), 1L, 12L)
tcga_sample_type <- function(x) substr(gsub("\\.", "-", x), 14L, 15L)

read_xena_matrix <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE, data.table = TRUE)
  id_col <- names(dt)[[1L]]
  ids <- as.character(dt[[id_col]])
  mat <- as.matrix(dt[, -1L, with = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- ids
  colnames(mat) <- tcga_sample(colnames(mat))
  mat
}

read_xena_table <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE, data.table = TRUE)
  sample_col <- names(dt)[tolower(names(dt)) %in% c("sample", "sampleid", "id", "barcode")]
  if (!length(sample_col)) sample_col <- names(dt)[[1L]] else sample_col <- sample_col[[1L]]
  dt[, sample := tcga_sample(get(sample_col))]
  dt[]
}

pick_col <- function(dt, patterns) {
  nms <- names(dt)
  for (p in patterns) {
    hit <- nms[grepl(p, nms, ignore.case = TRUE)]
    if (length(hit)) return(hit[[1L]])
  }
  NA_character_
}

ensure_bioc <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop(
      "Missing Bioconductor packages: ", paste(miss, collapse = ", "), "\n",
      "Install with:\n",
      "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
      "  BiocManager::install(c(", paste(sprintf('\"%s\"', miss), collapse = ", "), "))"
    )
  }
}

# ---------------------------------------------------------------------------
# GO gene sets: one list element per GO ID (never a pooled union)
# ---------------------------------------------------------------------------

build_go_gene_sets <- function(go_ids, symbol_universe) {
  ensure_bioc(c("org.Hs.eg.db", "AnnotationDbi"))
  sets <- vector("list", length(go_ids))
  names(sets) <- go_ids
  meta <- data.table(go_id = go_ids, go_name = NA_character_, n_annotated = 0L, n_in_expr = 0L, status = "ok")

  for (i in seq_along(go_ids)) {
    go_id <- go_ids[[i]]
    mapped <- tryCatch(
      AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = go_id,
        columns = c("SYMBOL", "ENSEMBL"),
        keytype = "GOALL"
      ),
      error = function(e) {
        warning(go_id, " mapping failed: ", conditionMessage(e))
        NULL
      }
    )
    all_symbols <- if (is.null(mapped) || !nrow(mapped)) character() else unique(na.omit(mapped$SYMBOL))
    symbols <- intersect(all_symbols, symbol_universe)
    meta$n_annotated[i] <- length(all_symbols)
    meta$n_in_expr[i] <- length(symbols)
    go_name <- NA_character_
    if (requireNamespace("GO.db", quietly = TRUE)) {
      go_name <- tryCatch(AnnotationDbi::Term(GO.db::GOTERM[[go_id]]), error = function(e) NA_character_)
    }
    meta$go_name[i] <- if (is.null(go_name) || !length(go_name)) NA_character_ else unname(go_name[[1L]])
    if (length(symbols) < MIN_GENES_PER_GO) {
      meta$status[i] <- paste0("skip: <", MIN_GENES_PER_GO, " genes in expression matrix")
      sets[[go_id]] <- NULL
    } else {
      sets[[go_id]] <- symbols
    }
  }
  sets <- sets[!vapply(sets, is.null, logical(1))]
  list(sets = sets, meta = meta)
}

# ---------------------------------------------------------------------------
# Pathway scores: one score vector per GO
# ---------------------------------------------------------------------------

ssgsea_per_go <- function(expr_symbol, go_sets) {
  ensure_bioc("GSVA")
  expr_symbol <- as.matrix(expr_symbol)
  if (requireNamespace("GSVA", quietly = TRUE) &&
      exists("ssgseaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
    param <- GSVA::ssgseaParam(expr_symbol, go_sets)
    scores <- GSVA::gsva(param, verbose = FALSE)
  } else {
    scores <- GSVA::gsva(expr_symbol, go_sets, method = "ssgsea", verbose = FALSE)
  }
  as.matrix(scores)
}

# ---------------------------------------------------------------------------
# Clinical / survival association for a single GO score
# ---------------------------------------------------------------------------

cox_one <- function(score, time, event, go_id, endpoint) {
  ok <- is.finite(score) & is.finite(time) & is.finite(event) & time > 0
  if (sum(ok) < 20L || length(unique(event[ok])) < 2L) {
    return(data.table(
      go_id = go_id, endpoint = endpoint, n = sum(ok),
      hr = NA_real_, hr_lo = NA_real_, hr_hi = NA_real_, p_value = NA_real_
    ))
  }
  fit <- survival::coxph(survival::Surv(time[ok], event[ok]) ~ score[ok])
  s <- summary(fit)
  data.table(
    go_id = go_id,
    endpoint = endpoint,
    n = sum(ok),
    hr = unname(s$conf.int[1, 1]),
    hr_lo = unname(s$conf.int[1, 3]),
    hr_hi = unname(s$conf.int[1, 4]),
    p_value = unname(s$waldtest["pvalue"])
  )
}

associate_one_go_clinical <- function(go_id, score, clinical, survival) {
  stopifnot(is.numeric(score))
  sc <- data.table(sample = names(score), go_score = as.numeric(score))
  sc[, patient := tcga_patient(sample)]

  merge_score <- function(other) {
    if (is.null(other) || !nrow(other)) return(NULL)
    m <- merge(sc, other, by = "sample")
    if (nrow(m) >= 20L) return(m)
    if (!"patient" %in% names(other)) other <- copy(other)[, patient := tcga_patient(sample)]
    merge(sc, unique(other, by = "patient"), by = "patient")
  }

  out_surv <- list()
  if (!is.null(survival)) {
    sv <- merge_score(survival)
    if (is.null(sv) || !nrow(sv)) sv <- data.table()
    os_time <- pick_col(sv, c("^OS\\.time$", "^os_time$", "overall_survival.*time"))
    os_evt <- pick_col(sv, c("^OS$", "^os$", "overall_survival$"))
    pfi_time <- pick_col(sv, c("^PFI\\.time$", "^pfi_time$", "progression.free.*time"))
    pfi_evt <- pick_col(sv, c("^PFI$", "^pfi$"))
    if (!requireNamespace("survival", quietly = TRUE)) {
      warning("Package 'survival' not installed; skipping Cox models for ", go_id)
    } else if (nrow(sv)) {
      if (!is.na(os_time) && !is.na(os_evt)) {
        out_surv[[length(out_surv) + 1L]] <- cox_one(
          sv$go_score, as.numeric(sv[[os_time]]), as.numeric(sv[[os_evt]]), go_id, "OS"
        )
      }
      if (!is.na(pfi_time) && !is.na(pfi_evt)) {
        out_surv[[length(out_surv) + 1L]] <- cox_one(
          sv$go_score, as.numeric(sv[[pfi_time]]), as.numeric(sv[[pfi_evt]]), go_id, "PFI"
        )
      }
    }
  }

  out_clin <- list()
  if (!is.null(clinical)) {
    cl <- merge_score(clinical)
    skip <- c("sample", "go_score", "patient")
    if (is.null(cl) || !nrow(cl)) cl <- data.table()
    for (col in setdiff(names(cl), skip)) {
      x <- cl[[col]]
      if (is.list(x)) next
      num <- suppressWarnings(as.numeric(as.character(x)))
      is_num <- mean(is.finite(num)) > 0.7
      if (!is_num) {
        n_uniq <- length(unique(as.character(x)[!is.na(x) & as.character(x) != ""]))
        if (n_uniq < 2L || n_uniq > 15L) next
      }
      if (is_num) {
        ok <- is.finite(cl$go_score) & is.finite(num)
        if (sum(ok) < 20L) next
        ct <- suppressWarnings(cor.test(cl$go_score[ok], num[ok], method = "spearman", exact = FALSE))
        out_clin[[length(out_clin) + 1L]] <- data.table(
          go_id = go_id, variable = col, test = "spearman",
          n = sum(ok), statistic = unname(ct$estimate), p_value = ct$p.value
        )
      } else {
        f <- factor(as.character(x))
        f[f %in% c("", "NA", "Unknown", "unknown", "[Unknown]", "[Not Available]", "[Not Evaluated]")] <- NA
        f <- droplevels(f)
        ok <- is.finite(cl$go_score) & !is.na(f)
        if (sum(ok) < 20L || nlevels(f[ok]) < 2L) next
        if (nlevels(f[ok]) == 2L) {
          wt <- wilcox.test(cl$go_score[ok] ~ f[ok], exact = FALSE)
          out_clin[[length(out_clin) + 1L]] <- data.table(
            go_id = go_id, variable = col, test = "wilcoxon",
            n = sum(ok), statistic = unname(wt$statistic), p_value = wt$p.value
          )
        } else {
          kt <- kruskal.test(cl$go_score[ok] ~ f[ok])
          out_clin[[length(out_clin) + 1L]] <- data.table(
            go_id = go_id, variable = col, test = "kruskal",
            n = sum(ok), statistic = unname(kt$statistic), p_value = kt$p.value
          )
        }
      }
    }
  }

  list(
    survival = if (length(out_surv)) rbindlist(out_surv) else NULL,
    clinical = if (length(out_clin)) rbindlist(out_clin) else NULL
  )
}

# ---------------------------------------------------------------------------
# Genes negatively correlated with ONE GO score (FDR corrected within that GO)
# ---------------------------------------------------------------------------

fast_neg_corr_features <- function(go_id, score, mat, feature_type = "gene") {
  common <- intersect(names(score), colnames(mat))
  score <- as.numeric(score[common])
  mat <- mat[, common, drop = FALSE]
  ok_s <- is.finite(score)
  score <- score[ok_s]
  mat <- mat[, ok_s, drop = FALSE]
  # Rank-based Spearman via Pearson on ranks
  score_r <- rank(score, ties.method = "average")
  mat_r <- apply(mat, 1L, function(x) {
    x[!is.finite(x)] <- NA
    if (sum(is.finite(x)) < 20L) return(rep(NA_real_, length(x)))
    rank(x, ties.method = "average", na.last = "keep")
  })
  mat_r <- t(mat_r)
  # Pairwise complete: drop features with NA ranks
  keep <- rowSums(is.finite(mat_r)) == ncol(mat_r)
  if (any(keep)) {
    rs <- apply(mat_r[keep, , drop = FALSE], 1L, sd, na.rm = TRUE)
    keep[which(keep)] <- is.finite(rs) & rs > 0
  }
  mat_r <- mat_r[keep, , drop = FALSE]
  mat <- mat[keep, , drop = FALSE]
  n <- ncol(mat_r)
  if (n < 20L || !nrow(mat_r)) {
    return(data.table(
      go_id = character(), feature = character(), feature_type = character(),
      rho = numeric(), p_value = numeric(), fdr = numeric()
    ))
  }
  score_z <- scale(score_r)[, 1]
  feat_z <- t(scale(t(mat_r)))
  rho <- as.numeric(feat_z %*% score_z / (n - 1))
  # t statistic for Spearman (approximate)
  tstat <- rho * sqrt((n - 2) / pmax(1e-12, 1 - rho^2))
  pval <- 2 * pt(-abs(tstat), df = n - 2)
  dt <- data.table(
    go_id = go_id,
    feature = rownames(mat_r),
    feature_type = feature_type,
    rho = rho,
    p_value = pval
  )
  dt <- dt[is.finite(rho) & is.finite(p_value)]
  dt[, fdr := p.adjust(p_value, method = "BH")]
  dt <- dt[rho < 0 & fdr < FDR_CUTOFF]
  if (MIN_RHO > 0) dt <- dt[abs(rho) >= MIN_RHO]
  setorder(dt, rho)
  dt[]
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  dir.create(file.path(OUT_DIR, "per_go"), recursive = TRUE, showWarnings = FALSE)

  message("Locating input files under ", DATA_DIR)
  f_clinical <- find_data_file("TCGA-BRCA.clinical")
  f_fpkm     <- find_data_file("TCGA-BRCA.star_fpkm")
  f_surv     <- find_data_file("TCGA-BRCA.survival")
  f_gtf      <- find_data_file("gencode.v36.annotation.gtf.gene")
  f_protein  <- tryCatch(find_data_file("TCGA-BRCA.protein"), error = function(e) {
    message("Protein file not found; skipping protein correlations")
    NULL
  })
  message("  clinical: ", f_clinical)
  message("  protein:  ", if (is.null(f_protein)) "(missing)" else f_protein)
  message("  fpkm:     ", f_fpkm)
  message("  survival: ", f_surv)
  message("  gencode:  ", f_gtf)

  message("Reading expression and annotation ...")
  gtf <- fread(f_gtf, sep = "\t", header = TRUE)
  gene_id_col <- pick_col(gtf, c("gene_id", "^id$"))
  gene_name_col <- pick_col(gtf, c("gene_name", "gene_symbol", "symbol"))
  gene_type_col <- pick_col(gtf, c("gene_type", "gene_biotype", "biotype"))
  if (is.na(gene_id_col) || is.na(gene_name_col)) {
    stop("gencode file must contain gene_id and gene_name columns. Found: ",
         paste(names(gtf), collapse = ", "))
  }
  gtf[, gene_id_clean := strip_ensembl_version(get(gene_id_col))]
  gtf[, gene_symbol := as.character(get(gene_name_col))]
  if (!is.na(gene_type_col)) {
    gtf <- gtf[get(gene_type_col) %in% c("protein_coding", "protein-coding")]
  }
  id2sym <- unique(gtf[, .(gene_id_clean, gene_symbol)])
  id2sym <- id2sym[gene_symbol != "" & !is.na(gene_symbol)]

  expr <- read_xena_matrix(f_fpkm)
  rownames(expr) <- strip_ensembl_version(rownames(expr))
  tumor <- tcga_sample_type(colnames(expr)) == TUMOR_SAMPLE_TYPE
  expr <- expr[, tumor, drop = FALSE]
  message("  tumor samples in FPKM: ", ncol(expr))

  # Map to symbols; average duplicate symbols
  map <- id2sym[match(rownames(expr), gene_id_clean)]
  keep <- !is.na(map$gene_symbol)
  expr <- expr[keep, , drop = FALSE]
  sym <- map$gene_symbol[keep]
  expr_log <- log2(pmax(expr, 0) + 1)
  expr_sym <- rowsum(expr_log, group = sym, reorder = FALSE)
  n_iso <- as.vector(table(factor(sym, levels = rownames(expr_sym))))
  expr_sym <- expr_sym / n_iso
  raw_sym <- rowsum(expr, group = sym, reorder = FALSE) / n_iso
  keep_var <- rowMeans(raw_sym > 1, na.rm = TRUE) >= 0.2
  expr_sym <- expr_sym[keep_var, , drop = FALSE]
  message("  unique gene symbols after low-expression filter: ", nrow(expr_sym))

  message("Building per-GO gene sets (GOALL; no pooling) ...")
  go_built <- build_go_gene_sets(GO_TERMS, rownames(expr_sym))
  fwrite(go_built$meta, file.path(OUT_DIR, "go_gene_set_sizes.csv"))
  if (!length(go_built$sets)) stop("No GO term had >= ", MIN_GENES_PER_GO, " genes in the expression matrix.")
  skipped <- go_built$meta[status != "ok"]
  if (nrow(skipped)) {
    message("Skipped GO terms:\n", paste(sprintf("  %s (%s)", skipped$go_id, skipped$status), collapse = "\n"))
  }

  message("ssGSEA for each GO term ...")
  scores <- ssgsea_per_go(expr_sym, go_built$sets)
  score_dt <- as.data.table(t(scores), keep.rownames = "sample")
  fwrite(score_dt, file.path(OUT_DIR, "go_ssgsea_scores.tsv"), sep = "\t")

  message("Reading clinical / survival / protein ...")
  clinical <- read_xena_table(f_clinical)
  survival <- read_xena_table(f_surv)
  protein <- if (is.null(f_protein)) NULL else tryCatch(read_xena_matrix(f_protein), error = function(e) {
    warning("Could not read protein matrix: ", conditionMessage(e))
    NULL
  })
  if (!is.null(protein)) {
    prot_tumor <- tcga_sample_type(colnames(protein)) == TUMOR_SAMPLE_TYPE
    protein <- protein[, prot_tumor, drop = FALSE]
  }

  all_surv <- list()
  all_clin <- list()
  all_neg_genes <- list()
  all_neg_prot <- list()

  go_ids <- rownames(scores)
  for (go_id in go_ids) {
    message("=== ", go_id, " (", which(go_ids == go_id), "/", length(go_ids), ") ===")
    score_vec <- setNames(as.numeric(scores[go_id, ]), colnames(scores))
    assoc <- associate_one_go_clinical(go_id, score_vec, clinical, survival)
    if (!is.null(assoc$survival)) {
      fwrite(assoc$survival, file.path(OUT_DIR, "per_go", paste0(GO_SAFE(go_id), "_survival.csv")))
      all_surv[[go_id]] <- assoc$survival
    }
    if (!is.null(assoc$clinical)) {
      assoc$clinical[, fdr := p.adjust(p_value, method = "BH")]
      fwrite(assoc$clinical, file.path(OUT_DIR, "per_go", paste0(GO_SAFE(go_id), "_clinical.csv")))
      all_clin[[go_id]] <- assoc$clinical
    }

    message("  negative gene correlations ...")
    neg_g <- fast_neg_corr_features(go_id, score_vec, expr_sym, "gene")
    fwrite(neg_g, file.path(OUT_DIR, "per_go", paste0(GO_SAFE(go_id), "_neg_corr_genes.csv")))
    all_neg_genes[[go_id]] <- neg_g
    message("    n_neg_genes = ", nrow(neg_g))

    if (!is.null(protein) && ncol(protein) > 20L) {
      message("  negative protein correlations ...")
      neg_p <- fast_neg_corr_features(go_id, score_vec, protein, "protein")
      fwrite(neg_p, file.path(OUT_DIR, "per_go", paste0(GO_SAFE(go_id), "_neg_corr_proteins.csv")))
      all_neg_prot[[go_id]] <- neg_p
      message("    n_neg_proteins = ", nrow(neg_p))
    }
  }

  if (length(all_surv)) fwrite(rbindlist(all_surv), file.path(OUT_DIR, "summary_survival_by_go.csv"))
  if (length(all_clin)) fwrite(rbindlist(all_clin), file.path(OUT_DIR, "summary_clinical_by_go.csv"))
  if (length(all_neg_genes)) {
    fwrite(rbindlist(all_neg_genes), file.path(OUT_DIR, "summary_neg_corr_genes_by_go.csv"))
    ntab <- rbindlist(lapply(names(all_neg_genes), function(g) {
      data.table(go_id = g, n_neg_genes = nrow(all_neg_genes[[g]]))
    }))
    fwrite(ntab, file.path(OUT_DIR, "summary_neg_gene_counts.csv"))
  }
  if (length(all_neg_prot)) fwrite(rbindlist(all_neg_prot), file.path(OUT_DIR, "summary_neg_corr_proteins_by_go.csv"))

  sink(file.path(OUT_DIR, "sessionInfo.txt"))
  print(sessionInfo())
  sink()
  message("Done. Per-GO tables are in ", file.path(OUT_DIR, "per_go"))
}

if (!interactive()) {
  main()
}
