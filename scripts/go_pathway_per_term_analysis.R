#!/usr/bin/env Rscript
# Per-GO (NOT pooled) TCGA-BRCA analysis, without GSVA / AnnotationDbi:
#   1) each listed GO term vs clinical / survival
#   2) genes negatively correlated with EACH GO term activity
#
# Pathway activity = combined z-score (mean of gene-wise z-scores in that GO).
# GO genes: local TSV/GMT, else QuickGO, else NCBI gene2go + go-basic.obo.
#
# Expected inputs in data_dir (prefix match is OK):
#   TCGA-BRCA.clinical, TCGA-BRCA.protein, TCGA-BRCA.star_fpkm,
#   TCGA-BRCA.survival, gencode.v36.annotation.gtf.gene

suppressPackageStartupMessages({
  library(data.table)
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
cache_dir <- file.path(data_dir, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

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
  if (!length(hits)) stop("Cannot find input matching '", pattern, "' under ", dir)
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

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(list(rho = NA_real_, p = NA_real_, n = sum(ok)))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

download_if_missing <- function(url, dest) {
  if (file.exists(dest) && isTRUE(file.info(dest)$size > 1000)) return(dest)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  message("Downloading ", url)
  download.file(url, dest, mode = "wb", quiet = TRUE)
  dest
}

# Combined z-score for EACH gene set separately (GSVA replacement).
# expr_mat: genes x samples; genesets: named list, one GO per element.
zscore_each_go <- function(expr_mat, genesets) {
  z <- t(scale(t(expr_mat)))
  z[!is.finite(z)] <- NA
  scores <- matrix(NA_real_, nrow = length(genesets), ncol = ncol(expr_mat),
                   dimnames = list(names(genesets), colnames(expr_mat)))
  for (gs in names(genesets)) {
    g <- intersect(genesets[[gs]], rownames(z))
    if (!length(g)) next
    scores[gs, ] <- colMeans(z[g, , drop = FALSE], na.rm = TRUE)
  }
  scores
}

read_local_genesets <- function(path, wanted) {
  if (grepl("\\.gmt$", path, ignore.case = TRUE)) {
    lines <- readLines(path, warn = FALSE)
    out <- lapply(wanted, function(x) character(0))
    names(out) <- wanted
    for (ln in lines) {
      p <- strsplit(ln, "\t", fixed = TRUE)[[1]]
      if (length(p) < 3) next
      gid <- p[[1]]
      if (!grepl("^GO:", gid) && grepl("GO:[0-9]+", p[[2]])) {
        gid <- regmatches(p[[2]], regexpr("GO:[0-9]+", p[[2]]))
      }
      if (gid %in% wanted) out[[gid]] <- unique(p[-(1:2)])
    }
    return(out)
  }
  dt <- fread(path)
  go_col <- pick_col(dt, c("GO_ID", "go_id", "GO", "term"))
  sym_col <- pick_col(dt, c("symbol", "gene_name", "Symbol", "gene"))
  if (is.null(go_col) || is.null(sym_col)) {
    stop("Local gene-set file must have GO_ID and symbol columns, or be GMT")
  }
  raw <- split(as.character(dt[[sym_col]]), as.character(dt[[go_col]]))
  out <- lapply(wanted, function(g) unique(na.omit(as.character(raw[[g]]))))
  names(out) <- wanted
  out
}

# QuickGO: term + descendants (same idea as GOALL), no AnnotationDbi.
fetch_quickgo_one <- function(go_id) {
  q <- paste0(
    "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch",
    "?goId=", utils::URLencode(go_id, reserved = TRUE),
    "&taxonId=9606&goUsage=descendants&geneProductType=protein"
  )
  tmp <- tempfile(fileext = ".tsv")
  ok <- FALSE
  if (nzchar(Sys.which("curl"))) {
    st <- suppressWarnings(system2(
      "curl",
      c("-sL", "-H", "Accept: text/tsv", "--fail", "-o", tmp, q),
      stdout = FALSE, stderr = FALSE
    ))
    ok <- identical(st, 0L) && file.exists(tmp) && isTRUE(file.info(tmp)$size > 20)
  }
  if (!ok) {
    ok <- tryCatch({
      download.file(q, tmp, mode = "wb", quiet = TRUE, headers = c(Accept = "text/tsv"))
      file.exists(tmp) && isTRUE(file.info(tmp)$size > 20)
    }, error = function(e) FALSE)
  }
  if (!ok) return(character(0))
  dt <- tryCatch(fread(tmp, sep = "\t"), error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(character(0))
  sym_col <- pick_col(dt, c("SYMBOL", "symbol", "GENE PRODUCT ID"))
  if (is.null(sym_col)) return(character(0))
  unique(na.omit(as.character(dt[[sym_col]])))
}

# NCBI gene2go + OBO descendants (offline after first download).
parse_go_children <- function(obo_path) {
  lines <- readLines(obo_path, warn = FALSE)
  children <- new.env(parent = emptyenv())
  current_id <- NA_character_
  in_term <- FALSE
  for (ln in lines) {
    if (identical(ln, "[Term]")) {
      in_term <- TRUE
      current_id <- NA_character_
    } else if (!nzchar(ln)) {
      in_term <- FALSE
    } else if (in_term && startsWith(ln, "id: GO:")) {
      current_id <- sub("^id: ", "", ln)
    } else if (in_term && startsWith(ln, "is_a: GO:") && !is.na(current_id)) {
      parent <- sub("^is_a: (GO:[0-9]+).*", "\\1", ln)
      children[[parent]] <- c(children[[parent]], current_id)
    }
  }
  children
}

go_descendants <- function(go_id, children_env) {
  out <- go_id
  queue <- go_id
  seen <- new.env(parent = emptyenv())
  seen[[go_id]] <- TRUE
  while (length(queue)) {
    x <- queue[[1]]
    queue <- queue[-1]
    kids <- children_env[[x]]
    if (is.null(kids)) next
    for (k in kids) {
      if (is.null(seen[[k]])) {
        seen[[k]] <- TRUE
        out <- c(out, k)
        queue <- c(queue, k)
      }
    }
  }
  unique(out)
}

fetch_ncbi_genesets <- function(wanted, cache_dir) {
  obo <- download_if_missing(
    "https://purl.obolibrary.org/obo/go/go-basic.obo",
    file.path(cache_dir, "go-basic.obo")
  )
  g2g_path <- download_if_missing(
    "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz",
    file.path(cache_dir, "gene2go.gz")
  )
  info_path <- download_if_missing(
    "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
    file.path(cache_dir, "Homo_sapiens.gene_info.gz")
  )
  message("Parsing GO parent/child relations")
  children <- parse_go_children(obo)
  term_map <- lapply(wanted, go_descendants, children_env = children)
  names(term_map) <- wanted
  all_terms <- unique(unlist(term_map, use.names = FALSE))

  message("Reading NCBI gene2go (human)")
  g2g <- fread(g2g_path)
  nms <- names(g2g)
  names(g2g) <- gsub("^#", "", nms)
  tax_col <- pick_col(g2g, c("tax_id", "taxon"))
  gene_col <- pick_col(g2g, c("GeneID", "Gene_ID"))
  go_col <- pick_col(g2g, c("GO_ID", "GO"))
  qual_col <- pick_col(g2g, c("Qualifier", "qualifier"))
  g2g <- g2g[as.character(get(tax_col)) == "9606" & get(go_col) %in% all_terms]
  if (!is.null(qual_col)) {
    g2g <- g2g[!grepl("(^|[|])NOT([|]|$)", as.character(get(qual_col)))]
  }

  info <- fread(info_path)
  names(info) <- gsub("^#", "", names(info))
  id_col <- pick_col(info, c("GeneID"))
  sym_col <- pick_col(info, c("Symbol", "symbol"))
  id2sym <- setNames(as.character(info[[sym_col]]), as.character(info[[id_col]]))

  out <- lapply(wanted, function(go_id) {
    ids <- unique(as.character(g2g[[gene_col]][g2g[[go_col]] %in% term_map[[go_id]]]))
    unique(na.omit(id2sym[ids]))
  })
  names(out) <- wanted
  out
}

empty_genesets <- function(wanted) {
  out <- lapply(wanted, function(x) character(0))
  names(out) <- wanted
  out
}

load_go_genesets <- function(wanted, data_dir, cache_dir) {
  local_hits <- list.files(
    data_dir,
    pattern = "go_genesets\\.(tsv|txt|csv|gmt)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(local_hits)) {
    message("Using local gene sets: ", local_hits[[1]])
    return(read_local_genesets(local_hits[[1]], wanted))
  }

  message("Fetching each GO gene set from QuickGO (no AnnotationDbi)")
  out <- empty_genesets(wanted)
  n_ok <- 0L
  for (go_id in wanted) {
    syms <- tryCatch(fetch_quickgo_one(go_id), error = function(e) character(0))
    out[[go_id]] <- setdiff(unique(syms), c("", "-", "NA"))
    if (length(out[[go_id]])) n_ok <- n_ok + 1L
    message("  ", go_id, ": ", length(out[[go_id]]), " genes")
  }
  if (n_ok > 0L) return(out)

  message("QuickGO unavailable; falling back to NCBI gene2go + go-basic.obo")
  fetch_ncbi_genesets(wanted, cache_dir)
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

stype <- sample_type_code(colnames(expr))
if (any(stype == "01", na.rm = TRUE)) {
  expr <- expr[, which(stype == "01"), drop = FALSE]
}

pid <- patient_id(colnames(expr))
keep <- !duplicated(pid) & nchar(pid) >= 12
expr <- expr[, keep, drop = FALSE]
colnames(expr) <- pid[keep]

if (max(expr, na.rm = TRUE) > 100) {
  expr <- log2(expr + 1)
}

symbol_for <- gene_map$symbol[match(rownames(expr), gene_map$ensembl_id)]
expr_symbol <- expr
rownames(expr_symbol) <- ifelse(is.na(symbol_for) | symbol_for == "" | symbol_for == "-",
                                rownames(expr), symbol_for)
expr_symbol <- expr_symbol[!duplicated(rownames(expr_symbol)), , drop = FALSE]

# ---------------------------------------------------------------------------
# GO gene sets: ONE list entry per GO ID (never collapse)
# ---------------------------------------------------------------------------
go_genesets <- load_go_genesets(go_ids, data_dir, cache_dir)
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
if (!length(usable)) {
  stop("No GO term has enough genes. Put data/go_genesets.tsv (columns GO_ID, symbol) and retry.")
}

# ---------------------------------------------------------------------------
# Pathway activity: one ROW PER GO (z-score mean, not a pooled signature)
# ---------------------------------------------------------------------------
message("Computing combined z-scores separately for each GO term")
pathway_scores <- zscore_each_go(expr_symbol, go_genesets[usable])
fwrite(
  data.table(GO_ID = rownames(pathway_scores), as.data.table(pathway_scores)),
  file.path(out_dir, "pathway_zscores_per_GO.tsv"),
  sep = "\t"
)
samples <- intersect(colnames(pathway_scores), colnames(expr))

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

clin_num_cols <- setdiff(names(clin)[vapply(clin, is.numeric, logical(1))], "patient")
is_cat <- vapply(clin, function(x) {
  if (is.numeric(x)) return(FALSE)
  u <- unique(as.character(x)[!is.na(x) & as.character(x) != ""])
  length(u) >= 2 && length(u) <= 12
}, logical(1))
clin_cat_cols <- setdiff(names(clin)[is_cat], c("patient", clin_id_col))

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
  score <- as.numeric(pathway_scores[go_id, samples])
  names(score) <- samples
  go_clin_dir <- file.path(out_dir, "clinical", gsub(":", "_", go_id))
  dir.create(go_clin_dir, recursive = TRUE, showWarnings = FALSE)

  for (cc in clin_num_cols) {
    idx <- intersect(names(score), clin$patient)
    st <- safe_spearman(score[idx], clin[[cc]][match(idx, clin$patient)])
    clin_rows[[length(clin_rows) + 1]] <- data.table(
      GO_ID = go_id, feature_type = "clinical_numeric", feature = cc,
      rho = st$rho, pvalue = st$p, n = st$n
    )
  }

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

  if (!is.null(surv$os_time) && !is.null(surv$os_event)) {
    sdf <- surv[patient %in% names(score), .(patient, os_time, os_event)]
    sdf[, pathway_score := score[patient]]
    sdf <- sdf[is.finite(os_time) & is.finite(os_event) & is.finite(pathway_score) & os_time > 0]
    if (nrow(sdf) >= 20 && length(unique(sdf$os_event)) == 2) {
      fit <- tryCatch(coxph(Surv(os_time, os_event) ~ pathway_score, data = sdf), error = function(e) NULL)
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
keep_genes <- apply(expr_use, 1, function(v) sd(v, na.rm = TRUE) > 0)
expr_use <- expr_use[keep_genes, , drop = FALSE]

neg_summary <- list()
for (go_id in usable) {
  score <- as.numeric(pathway_scores[go_id, samples])
  rho <- as.numeric(cor(t(expr_use), score, method = "spearman", use = "pairwise.complete.obs"))
  n <- ncol(expr_use)
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
