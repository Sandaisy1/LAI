#!/usr/bin/env Rscript
# =============================================================================
# 乳腺癌多器官转移单细胞：肿瘤类高表达基因 → 神经浸润
#
# 主队列：Klughammer et al. Nat Med 2024  scRNAseq.h5ad（SCP2702）
# 神经补齐：Gonzalez et al. Cell 2022  GSE186344
# 这是 scRNA，不是空间邻域。独立于 TG_RNAseq_*.R 与 BRCA_met_ST_*.R。
#
#   setwd("E:/R/BRCA_met_scRNA_nerve")
#   source("BRCA_met_scRNA_nerve_infiltration.R")
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")

need_cran <- c("dplyr", "tidyr", "tibble", "stringr", "ggplot2", "matrixStats")
need_bioc <- c("limma", "org.Hs.eg.db")
opt_cran <- c("writexl", "pheatmap", "ggrepel", "Matrix")
opt_bioc <- c("DESeq2", "clusterProfiler", "enrichplot", "fgsea", "msigdbr",
              "zellkonverter", "SingleCellExperiment", "SummarizedExperiment")

install_if_missing <- function(pkgs, bioc = FALSE, required = FALSE) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) == 0) return(invisible(TRUE))
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      tryCatch(install.packages("BiocManager", repos = "https://cloud.r-project.org"),
               error = function(e) NULL)
    }
    tryCatch(BiocManager::install(miss, update = FALSE, ask = FALSE),
             error = function(e) message("Bioc install failed: ", e$message))
  } else {
    tryCatch(install.packages(miss, repos = "https://cloud.r-project.org"),
             error = function(e) message("CRAN install failed: ", e$message))
  }
  still <- miss[!vapply(miss, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still) > 0 && required) stop("缺少 R 包: ", paste(still, collapse = ", "))
  if (length(still) > 0) message("可选包未装: ", paste(still, collapse = ", "))
  invisible(TRUE)
}
install_if_missing(need_cran, FALSE, TRUE)
install_if_missing(opt_cran, FALSE, FALSE)
install_if_missing(need_bioc, TRUE, FALSE)
install_if_missing(opt_bioc, TRUE, FALSE)
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
for (p in c(need_cran, opt_cran, need_bioc, opt_bioc)) {
  if (has_pkg(p)) suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# -----------------------------------------------------------------------------
# 1. 路径
# -----------------------------------------------------------------------------
looks_like_scrna <- function(d) {
  any(file.exists(file.path(d, c("scRNAseq.h5ad", "scRNA.h5ad")))) ||
    file.exists(file.path(d, "Klughammer_SCP2702", "scRNAseq.h5ad")) ||
    dir.exists(file.path(d, "GSE186344")) ||
    length(list.files(d, pattern = "matrix\\.mtx(\\.gz)?$", recursive = TRUE)) > 0 ||
    length(list.files(d, pattern = "\\.h5ad$", recursive = TRUE)) > 0
}

resolve_dir <- function() {
  env_dir <- Sys.getenv("BRCA_MET_SCRNA_DIR", unset = "")
  cands <- unique(c(env_dir, "E:/R/BRCA_met_scRNA_nerve", "E:\\R\\BRCA_met_scRNA_nerve", getwd()))
  cands <- cands[nzchar(cands)]
  for (d in cands) if (dir.exists(d) && looks_like_scrna(d)) {
    return(normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  if (dir.exists("E:/R/BRCA_met_scRNA_nerve")) {
    return(normalizePath("E:/R/BRCA_met_scRNA_nerve", winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

met_dir <- resolve_dir()
result_dir <- file.path(met_dir, "results")
log_dir <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
log_msg("scRNA data directory: ", met_dir)

p_cutoff <- 0.01
fc_cutoffs_local <- c("FC_1.25" = 1.25, "FC_1.5" = 1.5, "FC_2" = 2)
max_cells <- suppressWarnings(as.integer(Sys.getenv("BRCA_MET_SCRNA_MAX_CELLS", "150000")))
if (is.na(max_cells) || max_cells < 1000) max_cells <- 150000

# -----------------------------------------------------------------------------
# 2. 签名
# -----------------------------------------------------------------------------
tumor_epi_genes <- c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1", "SCGB2A2", "ESR1", "ERBB2")
immune_genes <- c("PTPRC", "CD3D", "CD3E", "CD8A", "CD4", "MS4A1", "CD79A", "CD68", "LYZ", "NKG7")
schwann_genes <- c("SOX10", "MPZ", "PMP22", "S100B", "PLP1", "NGFR", "NCAM1", "MBP", "L1CAM")
glia_genes <- c("GFAP", "AQP4", "ALDH1L1", "OLIG2", "MOG", "MBP", "PLP1", "CSPG4", "TMEM119", "P2RY12")
neuron_genes <- c("RBFOX3", "MAP2", "SYP", "SNAP25", "SLC17A7", "GAD1", "SYT1")

tumor_ligands <- data.frame(
  gene = c("NGF", "BDNF", "NTF3", "NTF4", "GDNF", "ARTN", "NRTN", "PSPN",
           "NTN1", "NTN4", "SLIT2", "SLIT3", "SEMA3A", "SEMA3C", "SEMA3F",
           "EFNA1", "EFNA5", "EFNB2", "NRG1", "CXCL12", "VEGFA", "MDK", "PTN",
           "TGFB1", "WNT5A", "L1CAM", "NCAM1", "CCL2"),
  receptor = c("NTRK1/NGFR", "NTRK2", "NTRK3", "NTRK2", "GFRA1/RET", "GFRA3/RET",
               "GFRA2/RET", "GFRA4/RET", "DCC/UNC5", "UNC5/ITGB4", "ROBO1/2", "ROBO",
               "NRP1/PLXNA", "NRP2/PLXNA", "NRP2", "EPHA", "EPHA", "EPHB",
               "ERBB2/3/4", "CXCR4", "VEGFR/NRP1", "PTPRZ1/ALK", "PTPRZ1",
               "TGFBR", "FZD/ROR2", "L1CAM/CNTN", "NCAM1", "CCR2"),
  stringsAsFactors = FALSE
)
receptor_genes <- unique(unlist(strsplit(gsub("/", ",", tumor_ligands$receptor), "[, ]+")))
receptor_genes <- receptor_genes[nzchar(receptor_genes) & !grepl("^(UNC5|ROBO|EPHA|EPHB|TGFBR|VEGFR)$", receptor_genes)]
receptor_genes <- unique(c(receptor_genes, "NTRK1", "NTRK2", "NTRK3", "NGFR", "GFRA1", "GFRA2",
                           "GFRA3", "RET", "DCC", "UNC5A", "ROBO1", "ROBO2", "NRP1", "NRP2",
                           "CXCR4", "ERBB2", "ERBB3", "PTPRZ1", "CCR2"))
nerve_term_pat <- paste("axon", "neur", "nerve", "schwann", "synap", "neurotroph",
                        "semaphorin", "ephrin", "glia", "oligodend", "astrocyte", sep = "|")

this_script_dir <- function() {
  ofile <- NULL
  if (sys.nframe() > 0) {
    for (i in sys.nframe():1) {
      e <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
      if (!is.null(e)) { ofile <- e; break }
    }
  }
  if (is.null(ofile)) getwd() else dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE))
}
load_pipeline_functions_only <- function() {
  cands <- c(file.path(this_script_dir(), "TG_RNAseq_pipeline.R"), "TG_RNAseq_pipeline.R",
             "E:/R/TG_BRCA/TG/TG_RNAseq_pipeline.R")
  pipe <- cands[file.exists(cands)][1]
  if (is.na(pipe) || !nzchar(pipe)) return(FALSE)
  lines <- readLines(pipe, warn = FALSE)
  main_at <- grep("^# 10\\. 主流程", lines)[1]
  if (is.na(main_at) || main_at < 2) return(FALSE)
  eval(parse(text = lines[seq_len(main_at - 1)]), envir = .GlobalEnv)
  TRUE
}
have_pipeline <- FALSE
if (!identical(Sys.getenv("BRCA_MET_SCRNA_SKIP_PIPELINE"), "1")) {
  have_pipeline <- tryCatch(load_pipeline_functions_only(), error = function(e) {
    log_msg("未载入 TG 作图函数: ", e$message); FALSE
  })
}

# -----------------------------------------------------------------------------
# 3. 元数据
# -----------------------------------------------------------------------------
norm_organ <- function(x) {
  x <- tolower(paste(x, collapse = " "))
  if (grepl("brain|cns|cerebr|brm|bcbm", x)) return("brain")
  if (grepl("lung|pulmon", x)) return("lung")
  if (grepl("bone|osseo|marrow|femur|spine|vertebra", x)) return("bone")
  if (grepl("liver|hepat", x)) return("liver")
  if (grepl("axilla|lymph", x)) return("axilla")
  if (grepl("chest|wall", x)) return("chest_wall")
  if (grepl("breast|mammar|primary", x)) return("breast")
  if (grepl("skin|derm", x)) return("skin")
  if (grepl("neck|cervic", x)) return("neck")
  "other"
}
pick_col <- function(df, keys) {
  nms <- names(df)
  low <- tolower(gsub("[^a-z0-9]", "", nms))
  for (k in keys) {
    kk <- tolower(gsub("[^a-z0-9]", "", k))
    hit <- which(low == kk | grepl(kk, low))
    if (length(hit) > 0) return(nms[hit[1]])
  }
  NA_character_
}
is_tumor_label <- function(x) {
  grepl("malig|tumor|tumour|cancer|epithel|brca|carcinoma|mtc", x, ignore.case = TRUE) &
    !grepl("immune|tcell|t cell|macro|fibro|endoth|bcell|plasma|nk |astro|oligo|microglia|neuron",
           x, ignore.case = TRUE)
}
is_immune_label <- function(x) {
  grepl("immune|tcell|t cell|cd8|cd4|nk|macro|mono|dc|dendrit|bcell|b cell|plasma|myeloid|lymph|treg",
        x, ignore.case = TRUE)
}
is_neural_label <- function(x) {
  grepl("schwann|nerve|neuron|astro|oligodend|glia|microglia|opc", x, ignore.case = TRUE)
}
present_genes <- function(genes, universe) intersect(unique(genes), universe)
log1p_cpm <- function(cnts) {
  lib <- pmax(rowSums(cnts), 1)
  log2(sweep(cnts, 1, lib, "/") * 1e4 + 1)
}
module_score <- function(logmat, genes) {
  g <- present_genes(genes, colnames(logmat))
  if (length(g) < 2) return(rep(NA_real_, nrow(logmat)))
  z <- scale(logmat[, g, drop = FALSE])
  z[is.na(z)] <- 0
  rowMeans(z)
}

# -----------------------------------------------------------------------------
# 4. 读入
# -----------------------------------------------------------------------------
read_tsv_flexible <- function(path) utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)

sce_to_list <- function(sce, source_file = "") {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) return(NULL)
  assays <- SummarizedExperiment::assays(sce)
  mat <- as.matrix(assays[[if ("counts" %in% names(assays)) "counts" else 1]])
  cd <- as.data.frame(SummarizedExperiment::colData(sce))
  rd <- as.data.frame(SummarizedExperiment::rowData(sce))
  gene_names <- rownames(rd)
  for (cc in c("gene_name", "Gene", "symbol", "feature_name")) {
    if (cc %in% names(rd)) { gene_names <- as.character(rd[[cc]]); break }
  }
  if (ncol(mat) == nrow(cd)) {
    mat <- t(mat)
    colnames(mat) <- gene_names
    obs <- cd
  } else if (nrow(mat) == nrow(cd)) {
    colnames(mat) <- if (length(gene_names) == ncol(mat)) gene_names else colnames(mat)
    obs <- cd
  } else {
    log_msg("h5ad 维度与 colData 对不上"); return(NULL)
  }
  obs$source_file <- source_file
  storage.mode(mat) <- "double"
  mat[is.na(mat)] <- 0
  list(cnts = mat, obs = obs)
}

read_h5ad <- function(path) {
  log_msg("Reading h5ad: ", path)
  if (has_pkg("zellkonverter") && has_pkg("SingleCellExperiment")) {
    sce <- tryCatch(zellkonverter::readH5AD(path, reader = "R", use_hdf5 = FALSE),
                    error = function(e) {
                      log_msg("zellkonverter R reader failed: ", e$message)
                      tryCatch(zellkonverter::readH5AD(path), error = function(e2) NULL)
                    })
    if (!is.null(sce)) return(sce_to_list(sce, basename(path)))
  }
  if (has_pkg("reticulate")) {
    obj <- tryCatch({
      ad <- reticulate::import("anndata", convert = FALSE)
      ad$read_h5ad(path)
    }, error = function(e) { log_msg("anndata failed: ", e$message); NULL })
    if (!is.null(obj)) {
      X <- tryCatch(as.matrix(reticulate::py_to_r(obj$X)), error = function(e) NULL)
      obs <- tryCatch(reticulate::py_to_r(obj$obs), error = function(e) NULL)
      var <- tryCatch(reticulate::py_to_r(obj$var), error = function(e) NULL)
      if (!is.null(X) && !is.null(obs) && nrow(X) == nrow(obs)) {
        genes <- if (!is.null(var) && "gene_name" %in% names(var)) var$gene_name else colnames(X)
        colnames(X) <- as.character(genes)
        obs$source_file <- basename(path)
        storage.mode(X) <- "double"
        X[is.na(X)] <- 0
        return(list(cnts = X, obs = as.data.frame(obs)))
      }
    }
  }
  log_msg("无法读取 ", path)
  NULL
}

gzfile_if <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
}

read_10x_dir <- function(d) {
  mtx <- list.files(d, pattern = "matrix\\.mtx(\\.gz)?$", full.names = TRUE)[1]
  bc <- list.files(d, pattern = "barcodes\\.(tsv|txt)(\\.gz)?$", full.names = TRUE)[1]
  ft <- list.files(d, pattern = "(features|genes)\\.(tsv|txt)(\\.gz)?$", full.names = TRUE)[1]
  if (is.na(mtx) || is.na(bc) || is.na(ft)) return(NULL)
  log_msg("Reading 10x: ", d)
  if (!has_pkg("Matrix")) {
    log_msg("需要 Matrix 包读 mtx"); return(NULL)
  }
  con_m <- if (grepl("\\.gz$", mtx, ignore.case = TRUE)) gzcon(file(mtx, "rb")) else file(mtx, "rb")
  on.exit(try(close(con_m), silent = TRUE), add = TRUE)
  m <- Matrix::readMM(con_m)
  barcodes <- readLines(gzfile_if(bc))
  feat <- utils::read.delim(gzfile_if(ft), header = FALSE, stringsAsFactors = FALSE)
  genes <- if (ncol(feat) >= 2) feat[[2]] else feat[[1]]
  if (nrow(m) == length(barcodes) && ncol(m) == length(genes)) m <- Matrix::t(m)
  if (ncol(m) != length(barcodes)) {
    log_msg("10x 维度异常: ", d); return(NULL)
  }
  mat <- as.matrix(Matrix::t(m))
  colnames(mat) <- make.unique(as.character(genes))
  rownames(mat) <- make.unique(as.character(barcodes))
  obs <- data.frame(barcode = rownames(mat), source_file = basename(d), stringsAsFactors = FALSE)
  storage.mode(mat) <- "double"
  list(cnts = mat, obs = obs)
}

find_h5ad <- function(root) unique(list.files(root, pattern = "\\.h5ad$", recursive = TRUE, full.names = TRUE))
find_10x_dirs <- function(root) {
  mt <- list.files(root, pattern = "matrix\\.mtx(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  unique(dirname(mt))
}

standardize_obs <- function(obs, fallback_patient = NA_character_, fallback_organ = "other") {
  obs <- as.data.frame(obs, stringsAsFactors = FALSE)
  pcol <- pick_col(obs, c("patient", "donor_id", "donor", "sample_id", "biosample",
                          "replicate", "condition", "sample", "orig.ident"))
  ocol <- pick_col(obs, c("biopsy_site", "site", "organ", "anatomic_site", "tissue",
                          "metastasis_site", "location", "type.of.tumor"))
  ccol <- pick_col(obs, c("cell_type", "celltype", "labels_cl_unif2_broad", "labels_unif",
                          "compartment", "annotation", "annot", "cell_ontology_class"))
  obs$patient <- if (!is.na(pcol)) as.character(obs[[pcol]]) else as.character(fallback_patient)
  if (all(is.na(obs$patient) | obs$patient == "")) obs$patient <- as.character(fallback_patient)
  org_raw <- if (!is.na(ocol)) as.character(obs[[ocol]]) else rep(fallback_organ, nrow(obs))
  if (length(org_raw) == 1) org_raw <- rep(org_raw, nrow(obs))
  miss <- is.na(org_raw) | org_raw == ""
  if (any(miss)) {
    extra <- paste(obs$patient, obs$source_file)
    org_raw[miss] <- extra[miss]
  }
  obs$organ <- vapply(as.character(org_raw), norm_organ, character(1))
  if (!is.na(fallback_organ) && fallback_organ == "brain") {
    obs$organ[obs$organ %in% c("other", "breast")] <- "brain"
  }
  obs$cell_label <- if (!is.na(ccol)) as.character(obs[[ccol]]) else ""
  obs
}

downsample_cells <- function(cnts, obs) {
  if (nrow(cnts) <= max_cells) return(list(cnts = cnts, obs = obs))
  log_msg("细胞数 ", nrow(cnts), " > ", max_cells, "，保留全部肿瘤/神经，其余下采样")
  # 粗分：后面 classify 会再标
  keep_force <- rep(FALSE, nrow(obs))
  if ("cell_label" %in% names(obs)) {
    keep_force <- is_tumor_label(obs$cell_label) | is_neural_label(obs$cell_label)
  }
  rest <- which(!keep_force)
  n_keep_rest <- max(1000, max_cells - sum(keep_force))
  if (length(rest) > n_keep_rest) {
    set.seed(1)
    rest <- sample(rest, n_keep_rest)
  }
  idx <- sort(unique(c(which(keep_force), rest)))
  list(cnts = cnts[idx, , drop = FALSE], obs = obs[idx, , drop = FALSE])
}

# -----------------------------------------------------------------------------
# 5. 细胞分类（无空间）
# -----------------------------------------------------------------------------
classify_cells <- function(cnts, obs, assay_name = "scRNA") {
  keep <- rowSums(cnts) >= 200
  if (sum(keep) < 50) {
    log_msg(assay_name, " 有效细胞太少"); return(NULL)
  }
  cnts <- cnts[keep, , drop = FALSE]
  obs <- obs[keep, , drop = FALSE]
  logmat <- log1p_cpm(cnts)
  obs$score_tumor <- module_score(logmat, tumor_epi_genes)
  obs$score_immune <- module_score(logmat, immune_genes)
  obs$score_schwann <- module_score(logmat, schwann_genes)
  obs$score_glia <- module_score(logmat, glia_genes)
  obs$score_neuron <- module_score(logmat, neuron_genes)
  lab <- obs$cell_label
  obs$is_tumor <- is_tumor_label(lab)
  obs$is_immune <- is_immune_label(lab)
  obs$is_neural <- is_neural_label(lab)
  no_lab <- is.na(lab) | lab == ""
  if (any(no_lab)) {
    obs$is_tumor[no_lab] <- !is.na(obs$score_tumor[no_lab]) &
      obs$score_tumor[no_lab] > pmax(obs$score_immune[no_lab], obs$score_glia[no_lab], 0.4, na.rm = TRUE)
    obs$is_immune[no_lab] <- !is.na(obs$score_immune[no_lab]) &
      obs$score_immune[no_lab] > pmax(obs$score_tumor[no_lab], 0.4, na.rm = TRUE)
    obs$is_neural[no_lab] <-
      (!is.na(obs$score_glia[no_lab]) & obs$score_glia[no_lab] > 0.8) |
      (!is.na(obs$score_neuron[no_lab]) & obs$score_neuron[no_lab] > 0.8) |
      (!is.na(obs$score_schwann[no_lab]) & obs$score_schwann[no_lab] > 0.8)
  }
  obs$is_tumor[is.na(obs$is_tumor)] <- FALSE
  obs$is_neural[is.na(obs$is_neural)] <- FALSE
  obs$is_immune[is.na(obs$is_immune)] <- FALSE
  obs$is_tumor[obs$is_neural | obs$is_immune] <- FALSE
  list(cnts = cnts, logmat = logmat, obs = obs, assay = assay_name)
}

# -----------------------------------------------------------------------------
# 6. DE / 作图
# -----------------------------------------------------------------------------
empty_de <- function() {
  data.frame(gene = character(), log2FC = numeric(), AveExpr = numeric(),
             pvalue = numeric(), padj = numeric(), padj_BH = numeric(),
             stringsAsFactors = FALSE)
}
limma_two_group <- function(logmat, group, patient = NULL, label = "") {
  if (is.null(logmat) || nrow(logmat) < 4 || ncol(logmat) < 20) return(empty_de())
  group <- factor(group, levels = c("far", "near"))
  if (nlevels(droplevels(group)) < 2) return(empty_de())
  mat <- t(logmat)
  design <- if (!is.null(patient) && length(unique(patient)) > 1) {
    stats::model.matrix(~ patient + group)
  } else stats::model.matrix(~ group)
  fit <- tryCatch({
    if (!has_pkg("limma")) return(NULL)
    limma::eBayes(limma::lmFit(mat, design), trend = TRUE, robust = TRUE)
  }, error = function(e) { log_msg("limma failed ", label, ": ", e$message); NULL })
  if (is.null(fit)) return(empty_de())
  coefn <- grep("groupnear", colnames(design), value = TRUE)
  if (length(coefn) == 0) coefn <- colnames(design)[ncol(design)]
  tt <- limma::topTable(fit, coef = coefn, number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), log2FC = tt$logFC, AveExpr = tt$AveExpr,
             pvalue = tt$P.Value, padj = tt$P.Value, padj_BH = tt$adj.P.Val,
             stringsAsFactors = FALSE)
}
mean_fc_two_group <- function(logmat, near_id, far_id) {
  if (length(near_id) < 3 || length(far_id) < 3) return(empty_de())
  log2FC <- colMeans(logmat[near_id, , drop = FALSE]) - colMeans(logmat[far_id, , drop = FALSE])
  data.frame(gene = colnames(logmat), log2FC = as.numeric(log2FC),
             AveExpr = colMeans(logmat[c(near_id, far_id), , drop = FALSE]),
             pvalue = NA_real_, padj = NA_real_, padj_BH = NA_real_,
             stringsAsFactors = FALSE)
}
deseq2_pseudobulk <- function(pb_counts, sample_info, label = "") {
  if (is.null(pb_counts) || ncol(pb_counts) < 4) return(empty_de())
  si <- sample_info
  si$group <- factor(si$group, levels = c("far", "near"))
  si$patient <- factor(si$patient)
  if (!has_pkg("DESeq2") || length(unique(si$patient)) < 2) {
    logm <- log2(sweep(t(pb_counts), 1, pmax(colSums(pb_counts), 1), "/") * 1e6 + 1)
    return(limma_two_group(logm, si$group, si$patient, label))
  }
  dds <- tryCatch({
    mode(pb_counts) <- "integer"
    dds <- DESeq2::DESeqDataSetFromMatrix(pb_counts, si, design = ~ patient + group)
    dds <- dds[rowSums(DESeq2::counts(dds)) >= 20, ]
    DESeq2::DESeq(dds, quiet = TRUE)
  }, error = function(e) { log_msg("DESeq2 failed ", label, ": ", e$message); NULL })
  if (is.null(dds)) {
    logm <- log2(sweep(t(pb_counts), 1, pmax(colSums(pb_counts), 1), "/") * 1e6 + 1)
    return(limma_two_group(logm, si$group, si$patient, paste(label, "limma_fallback")))
  }
  res <- DESeq2::results(dds, contrast = c("group", "near", "far"))
  data.frame(gene = rownames(res), log2FC = as.numeric(res$log2FoldChange),
             AveExpr = as.numeric(res$baseMean), pvalue = as.numeric(res$pvalue),
             padj = as.numeric(res$pvalue), padj_BH = as.numeric(res$padj),
             stringsAsFactors = FALSE)
}

basic_volcano <- function(de, title, outfile, fc_line = 1.25) {
  if (nrow(de) == 0) return(invisible(NULL))
  df <- de
  df$y <- if (any(!is.na(df$pvalue))) -log10(pmax(df$pvalue, 1e-300)) else abs(df$log2FC)
  df$col <- "ns"
  df$col[!is.na(df$pvalue) & df$pvalue < p_cutoff & df$log2FC >= log2(fc_line)] <- "up"
  if (all(is.na(df$pvalue))) df$col[df$log2FC >= log2(fc_line)] <- "up"
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  if (has_pkg("ggplot2")) {
    p <- ggplot2::ggplot(df, ggplot2::aes(log2FC, y, color = col)) +
      ggplot2::geom_point(alpha = 0.5, size = 0.8) +
      ggplot2::scale_color_manual(values = c(ns = "grey70", up = "#D62828")) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(title = title, y = "-log10(p)", x = "log2FC tumor_neural_context / other")
    ggplot2::ggsave(paste0(outfile, ".pdf"), p, width = 7, height = 6)
    ggplot2::ggsave(paste0(outfile, ".png"), p, width = 7, height = 6, dpi = 150)
  }
}

extract_focus_nerve <- function(outdir, tag) {
  csvs <- list.files(outdir, pattern = "ORA_.*\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(csvs) == 0) return(invisible(NULL))
  fdir <- file.path(outdir, "Focused_nerve_invasion")
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  chunks <- list()
  for (f in csvs) {
    df <- tryCatch(utils::read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    desc <- if ("Description" %in% names(df)) df$Description else df[[1]]
    hit <- grepl(nerve_term_pat, desc, ignore.case = TRUE)
    if (!any(hit)) next
    sub <- df[hit, , drop = FALSE]
    pcol <- intersect(c("pvalue", "p.adjust", "pval", "padj"), names(df))[1]
    if (length(pcol) == 1 && !is.na(pcol)) {
      ord <- order(df[[pcol]])
      sub$genome_wide_rank <- match(which(hit), ord)
    } else sub$genome_wide_rank <- seq_len(nrow(sub))
    sub$source_file <- basename(f)
    chunks[[length(chunks) + 1]] <- sub
  }
  if (length(chunks) == 0) {
    writeLines("no nerve-related ORA terms", file.path(fdir, paste0(tag, "_FOCUS_nerve_invasion_EMPTY.txt")))
    return(invisible(NULL))
  }
  hit <- do.call(rbind, lapply(chunks, function(d) d[, intersect(names(d), Reduce(intersect, lapply(chunks, names))), drop = FALSE]))
  utils::write.csv(hit, file.path(fdir, paste0(tag, "_FOCUS_nerve_invasion.csv")), row.names = FALSE)
}

emit_comparison <- function(comp_name, de, heat_mat, sample_info) {
  de <- de[!is.na(de$log2FC), , drop = FALSE]
  de$padj <- de$pvalue
  base <- file.path(result_dir, comp_name)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(de, file.path(base, "DE_full.csv"), row.names = FALSE)
  have_p <- any(!is.na(de$pvalue))
  log_msg(comp_name, " genes=", nrow(de), " have_p=", have_p)
  writeLines(
    c("scRNA 比较：上调 FC >= 1.25 / 1.5 / 2（先 p < 0.01）。",
      "不做 FC=1，不做 TopRank。不是空间近/远 spot。",
      "专项在 Focused_nerve_invasion/，不改全库 p。"),
    file.path(base, "00_READ_ME_先看这里.txt")
  )
  invisible(lapply(file.path(base, "FoldChange", names(fc_cutoffs_local)),
                   dir.create, recursive = TRUE, showWarnings = FALSE))
  gsea_cache <- list()
  use_pipe <- have_pipeline && exists("emit_subset_analysis", mode = "function")
  assign("result_dir", result_dir, envir = .GlobalEnv)
  assign("padj_cutoff", p_cutoff, envir = .GlobalEnv)
  for (nm in names(fc_cutoffs_local)) {
    fc <- unname(fc_cutoffs_local[[nm]])
    keep <- !is.na(de$log2FC) & de$log2FC > 0 & (2^de$log2FC >= fc)
    if (have_p) keep <- keep & !is.na(de$pvalue) & de$pvalue < p_cutoff
    sub <- de[keep, , drop = FALSE]
    if (nrow(sub) > 0) sub <- sub[order(sub$log2FC, decreasing = TRUE), , drop = FALSE]
    od <- file.path(base, "FoldChange", nm)
    if (use_pipe) {
      tryCatch(
        emit_subset_analysis(comp_name, sub, nm, paste0(comp_name, " | up FC >= ", fc),
                             od, de, heat_mat, sample_info, gsea_cache, fc_line = fc),
        error = function(e) log_msg("ERROR subset ", comp_name, " ", nm, ": ", e$message)
      )
    } else {
      utils::write.csv(sub, file.path(od, paste0(nm, "_DE_selected_genes.csv")), row.names = FALSE)
      basic_volcano(de, paste(comp_name, nm), file.path(od, paste0(nm, "_volcano")), fc)
    }
    tryCatch(extract_focus_nerve(od, nm), error = function(e) log_msg("focus: ", e$message))
  }
  if (use_pipe && exists("build_gsea_cache", mode = "function")) {
    tryCatch({
      gsea_cache <- build_gsea_cache(de)
      full_gsea_dir <- file.path(base, "00_GSEA_all_genes_NOT_FC_or_topN")
      dir.create(full_gsea_dir, recursive = TRUE, showWarnings = FALSE)
      writeLines("全基因 GSEA，不是 FC 分层。", file.path(full_gsea_dir, "00_README.txt"))
      for (nm in c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Hallmark")) {
        if (exists("plot_gsea_object", mode = "function")) {
          plot_gsea_object(gsea_cache[[nm]], file.path(full_gsea_dir, paste0("allGenes_GSEA_", nm)),
                           paste("GSEA", nm, "|", comp_name, "| ALL genes"))
        }
      }
    }, error = function(e) log_msg("GSEA failed: ", e$message))
  }
  invisible(de)
}

# -----------------------------------------------------------------------------
# 7. 载入各数据集
# -----------------------------------------------------------------------------
objects <- list()
add_obj <- function(obj, name) {
  if (is.null(obj)) return(invisible(NULL))
  objects[[name]] <<- obj
  log_msg(name, " cells=", nrow(obj$obs),
          " tumor=", sum(obj$obs$is_tumor),
          " immune=", sum(obj$obs$is_immune),
          " neural=", sum(obj$obs$is_neural),
          " organs=", paste(sort(unique(obj$obs$organ)), collapse = ","))
}

ingest_blob <- function(blob, name, fallback_organ = "other") {
  if (is.null(blob)) return(invisible(NULL))
  blob$obs <- standardize_obs(blob$obs, fallback_patient = name, fallback_organ = fallback_organ)
  ds <- downsample_cells(blob$cnts, blob$obs)
  obj <- classify_cells(ds$cnts, ds$obs, name)
  add_obj(obj, name)
}

# Klughammer scRNA（跳过 slide_seq / merfish）
klugh_roots <- c(file.path(met_dir, "Klughammer_SCP2702"), met_dir)
for (root in klugh_roots[dir.exists(klugh_roots)]) {
  h5s <- find_h5ad(root)
  spatial <- grepl("slide_seq|merfish|exseq|codex|visium|xenium", basename(h5s), ignore.case = TRUE)
  pref <- grepl("scrna|snrna|single", basename(h5s), ignore.case = TRUE) |
    grepl("^scRNAseq\\.h5ad$", basename(h5s))
  h5s <- unique(c(h5s[pref & !spatial], file.path(root, "scRNAseq.h5ad")))
  h5s <- h5s[file.exists(h5s)]
  if (length(h5s) == 0) h5s <- find_h5ad(root)[!spatial]

  for (f in unique(h5s)) {
    ingest_blob(tryCatch(read_h5ad(f), error = function(e) { log_msg(e$message); NULL }),
                paste0("Klughammer_", tools::file_path_sans_ext(basename(f))))
  }
}

# GSE186344：全部当脑转移
gse_dir <- file.path(met_dir, "GSE186344")
if (dir.exists(gse_dir)) {
  for (f in find_h5ad(gse_dir)) {
    ingest_blob(tryCatch(read_h5ad(f), error = function(e) NULL),
                paste0("GSE186344_", tools::file_path_sans_ext(basename(f))), "brain")
  }
  for (d in find_10x_dirs(gse_dir)) {
    ingest_blob(tryCatch(read_10x_dir(d), error = function(e) NULL),
                paste0("GSE186344_", basename(d)), "brain")
  }
}

if (length(objects) == 0) {
  log_msg("没有读到 scRNA 矩阵。请看 BRCA_met_scRNA_nerve_DOWNLOAD.txt")
}

# -----------------------------------------------------------------------------
# 8. 按样本标记「是否含神经细胞」，肿瘤细胞伪 bulk
# -----------------------------------------------------------------------------
qc_rows <- list()
pb_list <- list()
si_list <- list()
lr_rows <- list()

for (nm in names(objects)) {
  obj <- objects[[nm]]
  sp <- obj$obs
  sp$sample_key <- paste(sp$patient, sp$organ, sep = "__")
  neural_n <- tapply(sp$is_neural, sp$sample_key, sum)
  sp$sample_has_neural <- neural_n[sp$sample_key] >= 10
  obj$obs <- sp
  objects[[nm]] <- obj

  qc_rows[[nm]] <- do.call(rbind, lapply(split(sp, sp$sample_key), function(d) {
    data.frame(assay = nm, sample_key = d$sample_key[1], patient = d$patient[1],
               organ = d$organ[1], n = nrow(d), n_tumor = sum(d$is_tumor),
               n_immune = sum(d$is_immune), n_neural = sum(d$is_neural),
               sample_has_neural = d$sample_has_neural[1], stringsAsFactors = FALSE)
  }))

  tum <- which(sp$is_tumor)
  if (length(tum) < 20) {
    log_msg(nm, " 肿瘤细胞不足，跳过 DE")
    next
  }
  # 伪 bulk：patient x organ x (near=含神经样本的肿瘤 / far=不含)
  keys <- paste(sp$patient[tum], sp$organ[tum], ifelse(sp$sample_has_neural[tum], "near", "far"), sep = "__")
  uk <- unique(keys)
  if (length(uk) >= 2 && length(unique(ifelse(sp$sample_has_neural[tum], "near", "far"))) == 2) {
    pb <- sapply(uk, function(k) colSums(obj$cnts[tum[keys == k], , drop = FALSE]), simplify = "matrix")
    if (is.null(dim(pb))) pb <- matrix(pb, ncol = 1, dimnames = list(names(pb), uk))
    si <- do.call(rbind, strsplit(uk, "__", fixed = TRUE))
    si <- data.frame(sample = uk, patient = si[, 1], organ = si[, 2], group = si[, 3],
                     assay = nm, stringsAsFactors = FALSE)
    rownames(si) <- uk
    pb_list[[nm]] <- pb
    si_list[[nm]] <- si
  }

  # 配体–受体：肿瘤配体均值 vs 神经受体均值
  lig <- present_genes(tumor_ligands$gene, colnames(obj$logmat))
  rec <- present_genes(receptor_genes, colnames(obj$logmat))
  neu <- which(sp$is_neural)
  if (length(lig) > 0 && length(tum) >= 10) {
    lr_rows[[nm]] <- data.frame(
      assay = nm,
      gene = lig,
      tumor_mean = colMeans(obj$logmat[tum, lig, drop = FALSE]),
      neural_n = length(neu),
      receptor_mean_on_neural = if (length(neu) >= 5 && length(rec) > 0) {
        mean(colMeans(obj$logmat[neu, rec, drop = FALSE]))
      } else NA_real_,
      frac_tumor_in_neural_samples = mean(sp$sample_has_neural[tum]),
      stringsAsFactors = FALSE
    )
  }
}

if (length(qc_rows) > 0) {
  utils::write.csv(do.call(rbind, qc_rows), file.path(log_dir, "cell_class_counts.csv"), row.names = FALSE)
}

organ_de <- list()
# 合并伪 bulk（基因取交集），按器官做 DE
if (length(pb_list) > 0) {
  genes <- Reduce(intersect, lapply(pb_list, rownames))
  pb_all <- do.call(cbind, lapply(pb_list, function(m) m[genes, , drop = FALSE]))
  si_all <- do.call(rbind, si_list)
  si_all <- si_all[colnames(pb_all), , drop = FALSE]
  for (org in sort(unique(si_all$organ))) {
    si_o <- si_all[si_all$organ == org, , drop = FALSE]
    if (nrow(si_o) < 2 || length(unique(si_o$group)) < 2) {
      log_msg("器官 ", org, " 没有同时存在「含神经 / 不含神经」的肿瘤伪 bulk，跳过")
      next
    }
    pb_o <- pb_all[, si_o$sample, drop = FALSE]
    if (length(unique(si_o$patient)) < 2) {
      log_msg("器官 ", org, " 病人 < 2，不伪造 p，改用均值 FC")
      logm <- log2(sweep(t(pb_o), 1, pmax(colSums(pb_o), 1), "/") * 1e6 + 1)
      near <- which(si_o$group == "near")
      far <- which(si_o$group == "far")
      de_o <- mean_fc_two_group(logm, near, far)
    } else {
      de_o <- deseq2_pseudobulk(round(pmax(pb_o, 0)), si_o, paste("organ", org))
    }
    heat <- log2(sweep(pb_o, 2, pmax(colSums(pb_o), 1), "/") * 1e6 + 1)
    emit_comparison(paste0("organ_", org, "_tumor_in_neural_samples_vs_other"), de_o, heat, si_o)
    organ_de[[org]] <- de_o
  }
}

target_organs <- intersect(c("lung", "brain", "bone"), names(organ_de))
if (length(target_organs) >= 2) {
  get_up <- function(de, fc) {
    keep <- !is.na(de$log2FC) & de$log2FC > 0 & (2^de$log2FC >= fc)
    if (any(!is.na(de$pvalue))) keep <- keep & !is.na(de$pvalue) & de$pvalue < p_cutoff
    de$gene[keep]
  }
  for (nm in names(fc_cutoffs_local)) {
    sets <- lapply(organ_de[target_organs], get_up, fc = unname(fc_cutoffs_local[[nm]]))
    inter <- Reduce(intersect, sets)
    log_msg("common up lung/brain/bone ", nm, " n=", length(inter))
    utils::write.csv(data.frame(gene = inter),
                     file.path(result_dir, paste0("common_up_lung_brain_bone_", nm, ".csv")),
                     row.names = FALSE)
  }
} else {
  log_msg("肺/脑/骨中「含神经样本的肿瘤 vs 其他」够用的器官 < 2。脑请补 GSE186344。")
}

# -----------------------------------------------------------------------------
# 9. 候选分子：肿瘤配体 + 神经受体 + 含神经样本上调
# -----------------------------------------------------------------------------
lig_all <- if (length(lr_rows) > 0) do.call(rbind, lr_rows) else NULL
if (!is.null(lig_all)) {
  agg <- do.call(rbind, lapply(split(lig_all, lig_all$gene), function(d) {
    data.frame(gene = d$gene[1], n_assays = nrow(d),
               mean_tumor_expr = mean(d$tumor_mean, na.rm = TRUE),
               mean_frac_neural_samples = mean(d$frac_tumor_in_neural_samples, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
} else {
  agg <- data.frame(gene = tumor_ligands$gene, n_assays = 0,
                    mean_tumor_expr = NA, mean_frac_neural_samples = NA)
}
cand <- merge(tumor_ligands, agg, by = "gene", all.x = TRUE)
for (org in names(organ_de)) {
  d <- organ_de[[org]]
  cand[[paste0(org, "_neuralctx_log2FC")]] <- d$log2FC[match(cand$gene, d$gene)]
  cand[[paste0(org, "_neuralctx_pvalue")]] <- d$pvalue[match(cand$gene, d$gene)]
}
cand$up_in_neural_context <- FALSE
fc_org <- grep("_neuralctx_log2FC$", names(cand), value = TRUE)
if (length(fc_org) > 0) {
  up_any <- rep(FALSE, nrow(cand))
  for (fcn in fc_org) {
    org <- sub("_neuralctx_log2FC$", "", fcn)
    pvn <- paste0(org, "_neuralctx_pvalue")
    fc <- cand[[fcn]]
    pv <- if (pvn %in% names(cand)) cand[[pvn]] else rep(NA_real_, nrow(cand))
    up_any <- up_any | (!is.na(fc) & fc > 0 & (is.na(pv) | pv < p_cutoff))
  }
  cand$up_in_neural_context <- up_any
}
cand$priority <- 1 * (!is.na(cand$mean_tumor_expr) & cand$mean_tumor_expr > 0) +
  2 * (!is.na(cand$mean_frac_neural_samples) & cand$mean_frac_neural_samples > 0.5) +
  3 * cand$up_in_neural_context
cand <- cand[order(-cand$priority, -cand$mean_tumor_expr), ]
utils::write.csv(cand, file.path(result_dir, "01_CANDIDATE_MOLECULES_tumor_to_nerve.csv"), row.names = FALSE)
if (has_pkg("writexl")) {
  tryCatch(writexl::write_xlsx(cand, file.path(result_dir, "01_CANDIDATE_MOLECULES_tumor_to_nerve.xlsx")),
           error = function(e) NULL)
}

protocol <- c(
  "============================================================",
  "乳腺癌转移单细胞：肿瘤高表达基因与神经浸润",
  "数据目录: ", met_dir,
  "============================================================",
  "",
  "这是 scRNA，不能证明空间上的神经浸润。",
  "主队列 Klughammer scRNAseq.h5ad：多病人多器官，肿瘤+免疫；脑/肺各约 1 例。",
  "神经细胞请用 GSE186344。不要用小鼠 PDX，也不要跑空间邻域脚本。",
  "",
  "证据链：恶性上皮高表达配体 + 神经/胶质有受体 + 含神经样本里该配体更高。",
  "先 p < 0.01，再上调 FC >= 1.25 / 1.5 / 2。病人和器官不合并。",
  "器官间肿瘤差异在文件名带 organ_adaptation_NOT_infiltration 时，不要写成 PNI。",
  "",
  "先看 results/01_CANDIDATE_MOLECULES_tumor_to_nerve.csv",
  "以及 organ_*_tumor_in_neural_samples_vs_other/",
  "Focused_nerve_invasion/ 保留原始 p 与 genome_wide_rank。"
)
writeLines(protocol, file.path(result_dir, "00_PROTOCOL_神经浸润分析方案.txt"))
log_msg("Done. ", file.path(result_dir, "01_CANDIDATE_MOLECULES_tumor_to_nerve.csv"))
