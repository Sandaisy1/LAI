#!/usr/bin/env Rscript
# Per-GO (NOT pooled) TCGA-BRCA analysis:
#   1) each listed GO term vs clinical / survival
#   2) genes negatively correlated with EACH GO term activity
#
# Expected inputs in data_dir (prefix match is OK):
#   TCGA-BRCA.clinical, TCGA-BRCA.protein, TCGA-BRCA.star_fpkm,
#   TCGA-BRCA.survival, gencode.v36.annotation.gtf.gene

suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(survival)
})

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else "data"
out_dir  <- if (length(args) >= 2) args[[2]] else "results"
dir.create(file.path(out_dir, "clinical"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "negcorr"), recursive = TRUE, showWarnings = FALSE)

go_ids <- c(
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

min_genes <- 5
fdr_cutoff <- 0.05
rho_cutoff <- 0          # keep strictly negative correlations

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
find_input <- function(dir, pattern) {
  hits <- list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(hits) == 0) {
    stop("Cannot find input matching '", pattern, "' under ", dir)
  }
  hits[[1]]
}

patient_id <- function(x) {
  x <- gsub("\\.", "-", as.character(x))
  substr(x, 1, 12)
}

sample_type_code <- function(x) {
  x <- gsub("\\.", "-", as.character(x))
  parts <- strsplit(x, "-", fixed = TRUE)
  vapply(parts, function(p) if (length(p) >= 4) substr(p[[4]], 1, 2) else NA_character_, character(1))
}

read_matrix_like <- function(path) {
  dt <- fread(path)
  first <- names(dt)[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- as.character(dt[[first]])
  storage.mode(mat) <- "numeric"
  colnames(mat) <- gsub("\\.", "-", colnames(mat))
  mat
}

pick_col <- function(df, candidates) {
  nms <- names(df)
  exact <- candidates[candidates %in% nms]
  if (length(exact)) return(exact[[1]])
  low <- tolower(nms)
  for (cand in candidates) {
    hit <- which(low == tolower(cand))
    if (length(hit)) return(nms[[hit[[1]]]])
  }
  NULL
}

gsva_each_go <- function(expr_mat, genesets) {
  if (isTRUE(requireNamespace("GSVA", quietly = TRUE) &&
             exists("gsvaParam", where = asNamespace("GSVA"), inherits = FALSE))) {
    param <- GSVA::gsvaParam(exprData = expr_mat, geneSets = genesets, kcdf = "Gaussian")
    return(GSVA::gsva(param, verbose = FALSE))
  }
  GSVA::gsva(expr_mat, genesets, method = "gsva", kcdf = "Gaussian", verbose = FALSE)
}

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) {
    return(list(rho = NA_real_, p = NA_real_, n = sum(ok)))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

# ---------------------------------------------------------------------------
# Load annotation and expression
# ---------------------------------------------------------------------------
message("Loading inputs from: ", data_dir)
gtf_path <- find_input(data_dir, "gencode\\.v36\\.annotation\\.gtf\\.gene")
expr_path <- find_input(data_dir, "TCGA-BRCA\\.star_fpkm")
clin_path <- find_input(data_dir, "TCGA-BRCA\\.clinical")
surv_path <- find_input(data_dir, "TCGA-BRCA\\.survival")
prot_path <- tryCatch(find_input(data_dir, "TCGA-BRCA\\.protein"), error = function(e) NULL)

gtf <- fread(gtf_path)
ensembl_col <- pick_col(gtf, c("gene_id", "ensembl_id", "Ensembl", "id", names(gtf)[1]))
symbol_col  <- pick_col(gtf, c("gene_name", "symbol", "GeneSymbol", "hgnc_symbol"))
if (is.null(ensembl_col) || is.null(symbol_col)) {
  stop("gencode annotation must contain gene_id/ensembl and gene_name/symbol columns")
}
gtf[, ensembl_id := sub("\\..*$", "", get(ensembl_col))]
gene_map <- unique(gtf[, .(ensembl_id, symbol = get(symbol_col))])

expr <- read_matrix_like(expr_path)
rownames(expr) <- sub("\\..*$", "", rownames(expr))
expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]

# Keep primary tumor samples (01) when barcodes include sample type
stype <- sample_type_code(colnames(expr))
if (any(stype == "01", na.rm = TRUE)) {
  expr <- expr[, which(stype == "01"), drop = FALSE]
}

# One column per patient
pid <- patient_id(colnames(expr))
keep <- !duplicated(pid) & nchar(pid) >= 12
expr <- expr[, keep, drop = FALSE]
colnames(expr) <- pid[keep]

# log2(FPKM + 1); skip if already log-like
if (max(expr, na.rm = TRUE) > 100) {
  expr <- log2(expr + 1)
}

# Map ENSEMBL -> SYMBOL for GSVA gene sets; keep ENSEMBL matrix for genome-wide corr
symbol_for <- gene_map$symbol[match(rownames(expr), gene_map$ensembl_id)]
expr_symbol <- expr
rownames(expr_symbol) <- ifelse(is.na(symbol_for) | symbol_for == "" | symbol_for == "-",
                                rownames(expr), symbol_for)
expr_symbol <- expr_symbol[!duplicated(rownames(expr_symbol)), , drop = FALSE]

# ---------------------------------------------------------------------------
# GO gene sets: ONE list entry per GO ID (never collapse)
# ---------------------------------------------------------------------------
message("Fetching genes for each GO term via org.Hs.eg.db GOALL")
go_genesets <- lapply(go_ids, function(go_id) {
  mapped <- tryCatch(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = go_id,
      columns = c("SYMBOL", "ENSEMBL"),
      keytype = "GOALL"
    ),
    error = function(e) NULL
  )
  if (is.null(mapped) || nrow(mapped) == 0) return(character(0))
  unique(na.omit(mapped$SYMBOL))
})
names(go_genesets) <- go_ids

go_sizes <- data.table(
  GO_ID = go_ids,
  n_annotated = vapply(go_genesets, length, integer(1)),
  n_in_expr = vapply(go_genesets, function(g) sum(g %in% rownames(expr_symbol)), integer(1))
)
fwrite(go_sizes, file.path(out_dir, "go_geneset_sizes.tsv"), sep = "\t")

usable <- go_sizes$GO_ID[go_sizes$n_in_expr >= min_genes]
skipped <- setdiff(go_ids, usable)
if (length(skipped)) {
  message("Skipping GO terms with < ", min_genes, " genes in expression matrix: ",
          paste(skipped, collapse = ", "))
}
if (!length(usable)) stop("No GO term has enough genes in the expression matrix.")

# ---------------------------------------------------------------------------
# Pathway activity: GSVA matrix with one ROW PER GO (not a pooled signature)
# ---------------------------------------------------------------------------
message("Computing GSVA scores separately for each GO term")
gsva_scores <- gsva_each_go(expr_symbol, go_genesets[usable])
# gsva_scores: GO x samples
fwrite(
  data.table(GO_ID = rownames(gsva_scores), as.data.table(gsva_scores)),
  file.path(out_dir, "gsva_scores_per_GO.tsv"),
  sep = "\t"
)

samples <- intersect(colnames(gsva_scores), colnames(expr))

# ---------------------------------------------------------------------------
# Clinical + survival
# ---------------------------------------------------------------------------
clin <- fread(clin_path)
surv <- fread(surv_path)

clin_id_col <- pick_col(clin, c("submitter_id", "bcr_patient_barcode", "patient", "sample", "ID", names(clin)[1]))
surv_id_col <- pick_col(surv, c("sample", "submitter_id", "bcr_patient_barcode", "patient", "ID", names(surv)[1]))
clin[, patient := patient_id(get(clin_id_col))]
surv[, patient := patient_id(get(surv_id_col))]
clin <- clin[!duplicated(patient)]
surv <- surv[!duplicated(patient)]

time_col  <- pick_col(surv, c("OS.time", "OS.time", "_OS", "os_time", "overall_survival_time"))
event_col <- pick_col(surv, c("OS", "OS.event", "_EVENT", "os_event", "vital_status"))
if (!is.null(time_col)) surv[, os_time := as.numeric(get(time_col))]
if (!is.null(event_col)) {
  ev <- surv[[event_col]]
  if (is.character(ev) || is.factor(ev)) {
    surv[, os_event := as.integer(tolower(as.character(ev)) %in% c("1", "dead", "deceased", "event", "true"))]
  } else {
    surv[, os_event := as.integer(ev)]
  }
}

clin_num_cols <- names(clin)[vapply(clin, is.numeric, logical(1))]
clin_num_cols <- setdiff(clin_num_cols, c("patient"))

is_cat <- vapply(clin, function(x) {
  if (is.numeric(x)) return(FALSE)
  u <- unique(as.character(x)[!is.na(x) & as.character(x) != ""])
  length(u) >= 2 && length(u) <= 12
}, logical(1))
clin_cat_cols <- setdiff(names(clin)[is_cat], c("patient", clin_id_col))

# Optional protein matrix, aligned to patients
prot <- NULL
if (!is.null(prot_path)) {
  prot <- tryCatch(read_matrix_like(prot_path), error = function(e) NULL)
  if (!is.null(prot)) {
    looks_tcga <- function(x) grepl("^TCGA", gsub("\\.", "-", as.character(x)))
    if (mean(looks_tcga(rownames(prot))) > mean(looks_tcga(colnames(prot)))) {
      prot <- t(prot)
    }
    colnames(prot) <- patient_id(colnames(prot))
    prot <- prot[, !duplicated(colnames(prot)), drop = FALSE]
  }
}

# ---------------------------------------------------------------------------
# Task 1: EACH GO vs clinical (loop by GO_ID)
# ---------------------------------------------------------------------------
message("Task 1: correlating EACH GO pathway with clinical variables")
clin_rows <- list()

for (go_id in usable) {
  score <- as.numeric(gsva_scores[go_id, samples])
  names(score) <- samples
  go_clin_dir <- file.path(out_dir, "clinical", gsub(":", "_", go_id))
  dir.create(go_clin_dir, recursive = TRUE, showWarnings = FALSE)

  # numeric clinical variables
  for (cc in clin_num_cols) {
    idx <- intersect(names(score), clin$patient)
    st <- safe_spearman(score[idx], clin[[cc]][match(idx, clin$patient)])
    clin_rows[[length(clin_rows) + 1]] <- data.table(
      GO_ID = go_id, feature_type = "clinical_numeric", feature = cc,
      rho = st$rho, pvalue = st$p, n = st$n
    )
  }

  # categorical clinical variables (stage, subtype, etc.)
  for (cc in clin_cat_cols) {
    idx <- intersect(names(score), clin$patient)
    grp <- factor(as.character(clin[[cc]][match(idx, clin$patient)]))
    y <- score[idx]
    ok <- !is.na(grp) & is.finite(y)
    if (sum(ok) < 10 || nlevels(droplevels(grp[ok])) < 2) next
    kt <- tryCatch(kruskal.test(y[ok] ~ grp[ok]), error = function(e) NULL)
    if (is.null(kt)) next
    clin_rows[[length(clin_rows) + 1]] <- data.table(
      GO_ID = go_id, feature_type = "clinical_categorical", feature = cc,
      rho = unname(kt$statistic), pvalue = kt$p.value, n = sum(ok)
    )
  }

  # survival Cox (continuous GSVA score)
  if (!is.null(surv$os_time) && !is.null(surv$os_event)) {
    sdf <- surv[patient %in% names(score), .(patient, os_time, os_event)]
    sdf[, gsva := score[patient]]
    sdf <- sdf[is.finite(os_time) & is.finite(os_event) & is.finite(gsva) & os_time > 0]
    if (nrow(sdf) >= 20 && length(unique(sdf$os_event)) == 2) {
      fit <- tryCatch(coxph(Surv(os_time, os_event) ~ gsva, data = sdf), error = function(e) NULL)
      if (!is.null(fit)) {
        s <- summary(fit)
        clin_rows[[length(clin_rows) + 1]] <- data.table(
          GO_ID = go_id, feature_type = "survival_cox", feature = "OS",
          rho = unname(s$coefficients[1, "coef"]),
          pvalue = unname(s$coefficients[1, "Pr(>|z|)"]),
          n = nrow(sdf)
        )
      }
    }
  }

  # optional: each protein vs this GO score
  if (!is.null(prot)) {
    common <- intersect(names(score), colnames(prot))
    if (length(common) >= 10) {
      prot_rho <- cor(t(prot[, common, drop = FALSE]), score[common],
                      method = "spearman", use = "pairwise.complete.obs")[, 1]
      prot_tab <- data.table(GO_ID = go_id, protein = names(prot_rho), rho = as.numeric(prot_rho))
      fwrite(prot_tab, file.path(go_clin_dir, "protein_spearman.tsv"), sep = "\t")
    }
  }
}

clin_all <- rbindlist(clin_rows, fill = TRUE)
if (nrow(clin_all)) {
  clin_all[, fdr := p.adjust(pvalue, method = "BH"), by = GO_ID]
  fwrite(clin_all, file.path(out_dir, "clinical", "each_GO_vs_clinical.tsv"), sep = "\t")
  for (go_id in unique(clin_all$GO_ID)) {
    fwrite(
      clin_all[GO_ID == go_id],
      file.path(out_dir, "clinical", paste0(gsub(":", "_", go_id), "_clinical.tsv")),
      sep = "\t"
    )
  }
}

# ---------------------------------------------------------------------------
# Task 2: genes negatively correlated with EACH GO activity
# ---------------------------------------------------------------------------
message("Task 2: finding genes negatively correlated with EACH GO pathway")
expr_use <- expr[, samples, drop = FALSE]
# drop zero-variance genes
keep_genes <- apply(expr_use, 1, function(v) sd(v, na.rm = TRUE) > 0)
expr_use <- expr_use[keep_genes, , drop = FALSE]

neg_summary <- list()

for (go_id in usable) {
  score <- as.numeric(gsva_scores[go_id, samples])
  rho <- as.numeric(cor(t(expr_use), score, method = "spearman", use = "pairwise.complete.obs"))
  n <- ncol(expr_use)
  # Spearman p from rho (large-n approximation)
  tstat <- rho * sqrt((n - 2) / pmax(1e-12, 1 - rho^2))
  pval <- 2 * pt(-abs(tstat), df = n - 2)
  fdr <- p.adjust(pval, method = "BH")

  tab <- data.table(
    GO_ID = go_id,
    ensembl_id = rownames(expr_use),
    symbol = gene_map$symbol[match(rownames(expr_use), gene_map$ensembl_id)],
    rho = rho,
    pvalue = pval,
    fdr = fdr
  )
  tab[, in_this_GO := symbol %in% go_genesets[[go_id]]]
  neg <- tab[is.finite(rho) & rho < rho_cutoff & fdr < fdr_cutoff][order(rho)]
  fwrite(tab, file.path(out_dir, "negcorr", paste0(gsub(":", "_", go_id), "_all_genes.tsv")), sep = "\t")
  fwrite(neg, file.path(out_dir, "negcorr", paste0(gsub(":", "_", go_id), "_negative_genes.tsv")), sep = "\t")

  neg_summary[[go_id]] <- data.table(
    GO_ID = go_id,
    n_tested = nrow(tab),
    n_negative_fdr = nrow(neg),
    min_rho = if (nrow(neg)) min(neg$rho) else NA_real_
  )
  message("  ", go_id, ": ", nrow(neg), " negatively correlated genes (FDR < ", fdr_cutoff, ")")
}

fwrite(rbindlist(neg_summary), file.path(out_dir, "negcorr", "negative_gene_counts_per_GO.tsv"), sep = "\t")
message("Done. Per-GO clinical tables: ", file.path(out_dir, "clinical"))
message("Done. Per-GO negative-gene tables: ", file.path(out_dir, "negcorr"))
