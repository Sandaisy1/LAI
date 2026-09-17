#!/usr/bin/env Rscript
# =============================================================================
# 乳腺癌多器官转移空间转录组：肿瘤类高表达基因 → 神经浸润
#
# 主队列：Klughammer et al. Nat Med 2024（SCP2702）
# 神经补齐：GSE325935 人源脑转移 ST；Wang 2024 病理 Nerve 标注
# 独立于 TG_RNAseq_*.R，不修改原 Cuffdiff 流程。
#
# 用法：
#   setwd("E:/R/BRCA_met_ST_nerve")
#   source("BRCA_met_ST_nerve_infiltration.R")
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")

# -----------------------------------------------------------------------------
# 0. 轻量依赖（作图/富集可缺；读入失败才停）
# -----------------------------------------------------------------------------
need_cran <- c("dplyr", "tidyr", "tibble", "stringr", "ggplot2", "matrixStats")
need_bioc <- c("limma", "org.Hs.eg.db")
opt_cran <- c("writexl", "pheatmap", "ggrepel", "hdf5r")
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
  if (length(still) > 0) message("可选包未装，相关步骤将跳过: ", paste(still, collapse = ", "))
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
looks_like_klughammer <- function(d) {
  any(file.exists(file.path(d, c("slide_seq.h5ad", "scRNAseq.h5ad", "merfish.h5ad")))) ||
    any(file.exists(file.path(d, "Klughammer_SCP2702",
                              c("slide_seq.h5ad", "scRNAseq.h5ad", "merfish.h5ad")))) ||
    length(list.files(d, pattern = "counts\\.tsv$", recursive = TRUE)) > 0
}
looks_like_wang <- function(d) {
  dir.exists(file.path(d, "Robjects", "annotsBySpot")) ||
    dir.exists(file.path(d, "annotsBySpot")) ||
    file.exists(file.path(d, "Robjects.tar"))
}

resolve_met_dir <- function() {
  env_dir <- Sys.getenv("BRCA_MET_ST_DIR", unset = "")
  cands <- unique(c(
    env_dir,
    "E:/R/BRCA_met_ST_nerve", "E:\\R\\BRCA_met_ST_nerve",
    file.path(getwd(), "Klughammer_SCP2702"),
    getwd()
  ))
  cands <- cands[nzchar(cands)]
  for (d in cands) {
    if (dir.exists(d) && (looks_like_klughammer(d) || looks_like_wang(d) ||
                          dir.exists(file.path(d, "GSE325935")) ||
                          dir.exists(file.path(d, "Wang_TNBC")))) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  if (dir.exists("E:/R/BRCA_met_ST_nerve")) {
    return(normalizePath("E:/R/BRCA_met_ST_nerve", winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

met_dir <- resolve_met_dir()
result_dir <- file.path(met_dir, "results")
log_dir <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
log_msg("Data directory: ", met_dir)

p_cutoff <- 0.01
fc_cutoffs_local <- c("FC_1.25" = 1.25, "FC_1.5" = 1.5, "FC_2" = 2)

# -----------------------------------------------------------------------------
# 2. 签名与配体（肿瘤→神经）
# -----------------------------------------------------------------------------
tumor_epi_genes <- c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1", "SCGB2A2", "ESR1", "ERBB2")
immune_genes <- c("PTPRC", "CD3D", "CD3E", "CD8A", "CD4", "MS4A1", "CD79A", "CD68", "LYZ", "NKG7")
schwann_genes <- c("SOX10", "MPZ", "PMP22", "S100B", "PLP1", "NGFR", "NCAM1", "MBP", "L1CAM")
glia_genes <- c("GFAP", "AQP4", "ALDH1L1", "OLIG2", "MOG", "MBP", "PLP1", "CSPG4")
neuron_genes <- c("RBFOX3", "MAP2", "SYP", "SNAP25", "SLC17A7", "GAD1")

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
  rationale = c(
    "NGF–TrkA/NGFR，PNI / 感觉神经", "BDNF–TrkB", "NT-3–TrkC", "NT-4–TrkB",
    "GDNF 家族，促轴突与 PNI", "ARTN–GFRA3", "neurturin", "persephin",
    "netrin-1 轴突导向", "netrin-4", "SLIT–ROBO", "SLIT3",
    "SEMA3A", "SEMA3C 导向与侵袭", "SEMA3F", "ephrin-A1", "ephrin-A5", "ephrin-B2",
    "NRG1–ERBB，Schwann 互作", "CXCL12–CXCR4", "VEGF 血管/神经", "midkine",
    "pleiotrophin", "TGF-β / EMT", "WNT5A", "L1CAM 黏附", "NCAM1", "CCL2 趋化"
  ),
  stringsAsFactors = FALSE
)
nerve_term_pat <- paste(
  "axon", "neur", "nerve", "schwann", "synap", "neurotroph", "semaphorin",
  "ephrin", "slit", "netrin", "glia", "oligodend", "astrocyte",
  sep = "|"
)

# -----------------------------------------------------------------------------
# 3. 可选载入 TG 流程的作图函数（不跑主流程）
# -----------------------------------------------------------------------------
this_script_dir <- function() {
  ofile <- NULL
  if (sys.nframe() > 0) {
    for (i in sys.nframe():1) {
      e <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
      if (!is.null(e)) {
        ofile <- e
        break
      }
    }
  }
  if (is.null(ofile)) getwd() else dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE))
}

load_pipeline_functions_only <- function() {
  cands <- c(
    file.path(this_script_dir(), "TG_RNAseq_pipeline.R"),
    "TG_RNAseq_pipeline.R",
    file.path(getwd(), "TG_RNAseq_pipeline.R"),
    "E:/R/TG_BRCA/TG/TG_RNAseq_pipeline.R"
  )
  pipe <- cands[file.exists(cands)][1]
  if (is.na(pipe) || !nzchar(pipe)) return(FALSE)
  lines <- readLines(pipe, warn = FALSE)
  main_at <- grep("^# 10\\. 主流程", lines)[1]
  if (is.na(main_at) || main_at < 2) return(FALSE)
  eval(parse(text = lines[seq_len(main_at - 1)]), envir = .GlobalEnv)
  TRUE
}

have_pipeline <- FALSE
if (!identical(Sys.getenv("BRCA_MET_ST_SKIP_PIPELINE"), "1")) {
  have_pipeline <- tryCatch(load_pipeline_functions_only(), error = function(e) {
    log_msg("未载入 TG 作图函数: ", e$message)
    FALSE
  })
}

# -----------------------------------------------------------------------------
# 4. 元数据清洗
# -----------------------------------------------------------------------------
norm_organ <- function(x) {
  x <- tolower(paste(x, collapse = " "))
  if (grepl("brain|cns|cerebr", x)) return("brain")
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
  grepl("malig|tumor|tumour|cancer|epithel|brca|carcinoma|epithelial", x, ignore.case = TRUE) &
    !grepl("immune|tcell|t cell|macro|fibro|endoth|bcell|plasma|nk ", x, ignore.case = TRUE)
}
is_immune_label <- function(x) {
  grepl("immune|tcell|t cell|cd8|cd4|nk|macro|mono|dc|dendrit|bcell|b cell|plasma|myeloid|lymph",
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

min_dist_to <- function(xy, idx_from, idx_to) {
  if (length(idx_to) == 0) return(rep(Inf, length(idx_from)))
  src <- as.matrix(xy[idx_from, , drop = FALSE])
  tgt <- as.matrix(xy[idx_to, , drop = FALSE])
  apply(src, 1, function(p) min(sqrt((p[1] - tgt[, 1])^2 + (p[2] - tgt[, 2])^2)))
}

# -----------------------------------------------------------------------------
# 5. 读入：h5ad / counts.tsv+annot.tsv / 10x mtx
# -----------------------------------------------------------------------------
read_tsv_flexible <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

sce_to_list <- function(sce, source_file = "") {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) return(NULL)
  assays <- SummarizedExperiment::assays(sce)
  raw <- assays[[if ("counts" %in% names(assays)) "counts" else 1]]
  mat <- as.matrix(raw)
  cd <- as.data.frame(SummarizedExperiment::colData(sce))
  rd <- as.data.frame(SummarizedExperiment::rowData(sce))
  gene_names <- rownames(rd)
  for (cc in c("gene_name", "Gene", "symbol", "feature_name")) {
    if (cc %in% names(rd)) {
      gene_names <- as.character(rd[[cc]])
      break
    }
  }
  # SCE 惯例是 genes x cells；本脚本用 spots/cells x genes
  if (ncol(mat) == nrow(cd)) {
    mat <- t(mat)
    colnames(mat) <- gene_names
    obs <- cd
  } else if (nrow(mat) == nrow(cd)) {
    colnames(mat) <- if (length(gene_names) == ncol(mat)) gene_names else colnames(mat)
    obs <- cd
  } else {
    log_msg("h5ad 维度与 colData 对不上")
    return(NULL)
  }
  if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    red <- SingleCellExperiment::reducedDims(sce)
    nm <- intersect(c("spatial", "X_spatial"), names(red))
    if (length(nm) > 0) {
      xy <- as.matrix(red[[nm[1]]])
      obs$x <- xy[, 1]
      obs$y <- xy[, 2]
    }
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
    }, error = function(e) {
      log_msg("anndata/python failed: ", e$message)
      NULL
    })
    if (!is.null(obj)) {
      X <- tryCatch(as.matrix(reticulate::py_to_r(obj$X)), error = function(e) NULL)
      obs <- tryCatch(reticulate::py_to_r(obj$obs), error = function(e) NULL)
      var <- tryCatch(reticulate::py_to_r(obj$var), error = function(e) NULL)
      if (!is.null(X) && !is.null(obs)) {
        if (nrow(X) == nrow(obs)) {
          genes <- if (!is.null(var) && "gene_name" %in% names(var)) var$gene_name else colnames(X)
          colnames(X) <- as.character(genes)
          obsm <- tryCatch(reticulate::py_to_r(obj$obsm), error = function(e) NULL)
          if (!is.null(obsm)) {
            for (nm in c("spatial", "X_spatial")) {
              if (nm %in% names(obsm)) {
                xy <- as.matrix(obsm[[nm]])
                obs$x <- xy[, 1]
                obs$y <- xy[, 2]
              }
            }
          }
          obs$source_file <- basename(path)
          storage.mode(X) <- "double"
          X[is.na(X)] <- 0
          return(list(cnts = X, obs = as.data.frame(obs)))
        }
      }
    }
  }
  log_msg("无法读取 ", path, "。请安装 zellkonverter 或 Python anndata，或改用 counts.tsv+annot.tsv。")
  NULL
}

read_counts_annot_pair <- function(count_f, annot_f) {
  log_msg("Reading TSV pair: ", count_f)
  cnt <- read_tsv_flexible(count_f)
  ann <- read_tsv_flexible(annot_f)
  idc <- names(cnt)[1]
  rownames(cnt) <- make.unique(as.character(cnt[[idc]]))
  mat <- as.matrix(cnt[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  mat[is.na(mat)] <- 0
  ida <- names(ann)[1]
  rownames(ann) <- make.unique(as.character(ann[[ida]]))
  common <- intersect(rownames(mat), rownames(ann))
  if (length(common) < 10) {
    # 可能是 genes x spots
    if (ncol(mat) == nrow(ann)) {
      mat <- t(mat)
      rownames(mat) <- rownames(ann)
      common <- rownames(ann)
    }
  }
  common <- intersect(rownames(mat), rownames(ann))
  if (length(common) < 10) {
    log_msg("counts/annot 对不上: ", count_f)
    return(NULL)
  }
  list(cnts = mat[common, , drop = FALSE],
       obs = ann[common, , drop = FALSE])
}

find_tsv_pairs <- function(root) {
  fs <- list.files(root, pattern = "^counts\\.tsv$", recursive = TRUE, full.names = TRUE)
  out <- list()
  for (f in fs) {
    ann <- file.path(dirname(f), "annot.tsv")
    if (!file.exists(ann)) ann <- file.path(dirname(f), "annotation.tsv")
    if (file.exists(ann)) out[[length(out) + 1]] <- list(counts = f, annot = ann)
  }
  out
}

find_h5ad <- function(root) {
  unique(list.files(root, pattern = "\\.h5ad$", recursive = TRUE, full.names = TRUE))
}

standardize_obs <- function(obs, fallback_patient = NA_character_, fallback_organ = "other") {
  obs <- as.data.frame(obs, stringsAsFactors = FALSE)
  pcol <- pick_col(obs, c("patient", "donor_id", "donor", "sample_id", "biosample",
                          "HTAN_Biospecimen_ID", "replicate", "condition", "sample"))
  ocol <- pick_col(obs, c("biopsy_site", "site", "organ", "anatomic_site", "tissue",
                          "metastasis_site", "location"))
  ccol <- pick_col(obs, c("cell_type", "celltype", "labels_cl_unif2_broad",
                          "labels_unif", "labels_cl_unif", "compartment", "compartments",
                          "annotation", "annot"))
  xcol <- pick_col(obs, c("x", "array_col", "X", "spatial_1", "pxl_col_in_fullres"))
  ycol <- pick_col(obs, c("y", "array_row", "Y", "spatial_2", "pxl_row_in_fullres"))
  scol <- pick_col(obs, c("slide", "slide_id", "puck", "section", "source_file"))

  obs$patient <- if (!is.na(pcol)) as.character(obs[[pcol]]) else fallback_patient
  if (all(is.na(obs$patient) | obs$patient == "")) {
    obs$patient <- fallback_patient
  }
  org_raw <- if (!is.na(ocol)) as.character(obs[[ocol]]) else rep(fallback_organ, nrow(obs))
  if (length(org_raw) == 1) org_raw <- rep(org_raw, nrow(obs))
  miss <- is.na(org_raw) | org_raw == "" | org_raw == "other"
  if (any(miss)) {
    extra <- paste(obs$patient, if (!is.na(scol)) obs[[scol]] else "")
    org_raw[miss] <- extra[miss]
  }
  obs$organ <- vapply(as.character(org_raw), norm_organ, character(1))
  obs$cell_label <- if (!is.na(ccol)) as.character(obs[[ccol]]) else ""
  obs$x <- if (!is.na(xcol)) suppressWarnings(as.numeric(obs[[xcol]])) else NA_real_
  obs$y <- if (!is.na(ycol)) suppressWarnings(as.numeric(obs[[ycol]])) else NA_real_
  obs$slide_id <- if (!is.na(scol)) as.character(obs[[scol]]) else "slide1"
  obs$slide_id[is.na(obs$slide_id) | obs$slide_id == ""] <- "slide1"
  obs
}

# -----------------------------------------------------------------------------
# 6. 分类：肿瘤 / 免疫 / 神经 + 近/远神经
# -----------------------------------------------------------------------------
classify_spatial <- function(cnts, obs, assay_name = "spatial") {
  keep <- rowSums(cnts) >= 30
  if (sum(keep) < 40) {
    log_msg(assay_name, " 有效观测太少，跳过")
    return(NULL)
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
  # 无标签时用签名（肿瘤必须高于免疫/神经，避免把巨噬细胞当肿瘤）
  no_lab <- is.na(lab) | lab == ""
  if (any(no_lab)) {
    obs$is_tumor[no_lab] <- !is.na(obs$score_tumor[no_lab]) &
      obs$score_tumor[no_lab] > pmax(obs$score_immune[no_lab], obs$score_schwann[no_lab],
                                     obs$score_glia[no_lab], 0.4, na.rm = TRUE)
    obs$is_immune[no_lab] <- !is.na(obs$score_immune[no_lab]) &
      obs$score_immune[no_lab] > pmax(obs$score_tumor[no_lab], 0.4, na.rm = TRUE)
    obs$is_neural[no_lab] <- (!is.na(obs$score_schwann[no_lab]) & obs$score_schwann[no_lab] > 0.8) |
      (!is.na(obs$score_glia[no_lab]) & obs$score_glia[no_lab] > 0.8)
  }
  # 切片内 Schwann/胶质 z>1 也算神经
  obs$is_schwann_high <- FALSE
  for (sl in unique(obs$slide_id)) {
    w <- which(obs$slide_id == sl)
    for (sc in c("score_schwann", "score_glia")) {
      v <- obs[[sc]][w]
      if (all(is.na(v))) next
      sdv <- stats::sd(v, na.rm = TRUE)
      if (is.na(sdv) || sdv == 0) sdv <- 1
      obs$is_neural[w] <- obs$is_neural[w] | (v > mean(v, na.rm = TRUE) + sdv)
    }
  }
  obs$is_tumor[is.na(obs$is_tumor)] <- FALSE
  obs$is_neural[is.na(obs$is_neural)] <- FALSE
  obs$is_immune[is.na(obs$is_immune)] <- FALSE
  # 不要把神经/免疫标成肿瘤
  obs$is_tumor[obs$is_neural | obs$is_immune] <- FALSE

  obs$dist_neural <- Inf
  has_xy <- sum(is.finite(obs$x) & is.finite(obs$y)) > 20
  if (has_xy) {
    xy <- as.matrix(obs[, c("x", "y")])
    for (sl in unique(obs$slide_id)) {
      w <- which(obs$slide_id == sl)
      neu <- w[obs$is_neural[w]]
      if (length(neu) > 0) obs$dist_neural[w] <- min_dist_to(xy, w, neu)
    }
  }
  obs$near_nerve <- FALSE
  obs$far_nerve <- FALSE
  for (sl in unique(obs$slide_id)) {
    w <- which(obs$slide_id == sl & obs$is_tumor)
    d <- obs$dist_neural[w]
    ok <- is.finite(d)
    if (sum(ok) < 8) next
    q <- stats::quantile(d[ok], probs = c(0.25, 0.75), na.rm = TRUE)
    obs$near_nerve[w][ok] <- d[ok] <= q[[1]]
    obs$far_nerve[w][ok] <- d[ok] >= q[[2]]
  }
  list(cnts = cnts, logmat = logmat, obs = obs, has_xy = has_xy, assay = assay_name)
}

spot_map_plot <- function(obs, title, outfile) {
  if (!has_pkg("ggplot2") || !sum(is.finite(obs$x) & is.finite(obs$y))) return(invisible(NULL))
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  df <- obs
  df$cls <- "other"
  df$cls[df$is_immune] <- "immune"
  df$cls[df$is_tumor] <- "tumor_far"
  df$cls[df$near_nerve] <- "tumor_near"
  df$cls[df$is_neural] <- "neural"
  pal <- c(other = "grey85", immune = "#54A24B", tumor_far = "#4C78A8",
           tumor_near = "#D62828", neural = "#2A9D8F")
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, color = cls)) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::facet_wrap(~paste(patient, organ, slide_id, sep = " | ")) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(title = title, color = NULL)
  ggplot2::ggsave(paste0(outfile, ".pdf"), p, width = 10, height = 7)
  ggplot2::ggsave(paste0(outfile, ".png"), p, width = 10, height = 7, dpi = 140)
}

# -----------------------------------------------------------------------------
# 7. 差异表达
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
  } else {
    stats::model.matrix(~ group)
  }
  fit <- tryCatch({
    if (!has_pkg("limma")) return(NULL)
    limma::eBayes(limma::lmFit(mat, design), trend = TRUE, robust = TRUE)
  }, error = function(e) {
    log_msg("limma failed ", label, ": ", e$message)
    NULL
  })
  if (is.null(fit)) return(empty_de())
  coefn <- grep("groupnear", colnames(design), value = TRUE)
  if (length(coefn) == 0) coefn <- colnames(design)[ncol(design)]
  tt <- limma::topTable(fit, coef = coefn, number = Inf, sort.by = "none")
  data.frame(
    gene = rownames(tt), log2FC = tt$logFC, AveExpr = tt$AveExpr,
    pvalue = tt$P.Value, padj = tt$P.Value, padj_BH = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
}

mean_fc_two_group <- function(logmat, near_id, far_id) {
  if (length(near_id) < 3 || length(far_id) < 3) return(empty_de())
  log2FC <- colMeans(logmat[near_id, , drop = FALSE]) - colMeans(logmat[far_id, , drop = FALSE])
  data.frame(
    gene = colnames(logmat), log2FC = as.numeric(log2FC),
    AveExpr = colMeans(logmat[c(near_id, far_id), , drop = FALSE]),
    pvalue = NA_real_, padj = NA_real_, padj_BH = NA_real_,
    stringsAsFactors = FALSE
  )
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
  }, error = function(e) {
    log_msg("DESeq2 failed ", label, ": ", e$message)
    NULL
  })
  if (is.null(dds)) {
    logm <- log2(sweep(t(pb_counts), 1, pmax(colSums(pb_counts), 1), "/") * 1e6 + 1)
    return(limma_two_group(logm, si$group, si$patient, paste(label, "limma_fallback")))
  }
  res <- DESeq2::results(dds, contrast = c("group", "near", "far"))
  data.frame(
    gene = rownames(res),
    log2FC = as.numeric(res$log2FoldChange),
    AveExpr = as.numeric(res$baseMean),
    pvalue = as.numeric(res$pvalue),
    padj = as.numeric(res$pvalue),
    padj_BH = as.numeric(res$padj),
    stringsAsFactors = FALSE
  )
}

pseudobulk_near_far <- function(cnts, obs) {
  keep <- obs$near_nerve | obs$far_nerve
  if (sum(keep) < 20) return(NULL)
  df <- obs[keep, , drop = FALSE]
  mat <- cnts[keep, , drop = FALSE]
  df$key <- paste(df$patient, df$organ, ifelse(df$near_nerve, "near", "far"), sep = "__")
  keys <- unique(df$key)
  pb <- sapply(keys, function(k) colSums(mat[df$key == k, , drop = FALSE]), simplify = "matrix")
  if (is.null(dim(pb))) {
    pb <- matrix(pb, ncol = 1, dimnames = list(names(pb), keys))
  }
  si <- do.call(rbind, strsplit(keys, "__", fixed = TRUE))
  si <- data.frame(sample = keys, patient = si[, 1], organ = si[, 2],
                   group = si[, 3], stringsAsFactors = FALSE)
  rownames(si) <- keys
  list(pb = pb, si = si)
}

# -----------------------------------------------------------------------------
# 8. 作图：FC 分层 + 神经专项（不改全库 p）
# -----------------------------------------------------------------------------
basic_volcano <- function(de, title, outfile, fc_line = 1.25) {
  if (nrow(de) == 0) return(invisible(NULL))
  df <- de
  df$y <- if (any(!is.na(df$pvalue))) -log10(pmax(df$pvalue, 1e-300)) else abs(df$log2FC)
  df$col <- "ns"
  df$col[!is.na(df$pvalue) & df$pvalue < p_cutoff & df$log2FC >= log2(fc_line)] <- "up"
  if (all(is.na(df$pvalue))) {
    df$col[df$log2FC >= log2(fc_line)] <- "up"
  }
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  if (has_pkg("ggplot2")) {
    p <- ggplot2::ggplot(df, ggplot2::aes(log2FC, y, color = col)) +
      ggplot2::geom_point(alpha = 0.5, size = 0.8) +
      ggplot2::scale_color_manual(values = c(ns = "grey70", up = "#D62828")) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(title = title, y = "-log10(p)", x = "log2FC tumor_near / tumor_far")
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
    desc <- if ("Description" %in% names(df)) df$Description else if ("ID" %in% names(df)) df$ID else df[[1]]
    hit <- grepl(nerve_term_pat, desc, ignore.case = TRUE)
    if (!any(hit)) next
    sub <- df[hit, , drop = FALSE]
    pcol <- intersect(c("pvalue", "p.adjust", "pval", "padj"), names(df))[1]
    if (length(pcol) == 1 && !is.na(pcol)) {
      ord <- order(df[[pcol]])
      sub$genome_wide_rank <- match(which(hit), ord)
    } else {
      sub$genome_wide_rank <- seq_len(nrow(sub))
    }
    sub$source_file <- basename(f)
    chunks[[length(chunks) + 1]] <- sub
  }
  if (length(chunks) == 0) {
    writeLines("no nerve-related ORA terms; original genome-wide p unchanged",
               file.path(fdir, paste0(tag, "_FOCUS_nerve_invasion_EMPTY.txt")))
    return(invisible(NULL))
  }
  hit <- do.call(rbind, lapply(chunks, function(d) {
    d[, intersect(names(d), unique(unlist(lapply(chunks, names)))), drop = FALSE]
  }))
  utils::write.csv(hit, file.path(fdir, paste0(tag, "_FOCUS_nerve_invasion.csv")), row.names = FALSE)
}

emit_comparison <- function(comp_name, de, heat_mat, sample_info) {
  de <- de[!is.na(de$log2FC), , drop = FALSE]
  de$padj <- de$pvalue
  base <- file.path(result_dir, comp_name)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(de, file.path(base, "DE_full.csv"), row.names = FALSE)
  have_p <- any(!is.na(de$pvalue))
  log_msg(comp_name, " genes=", nrow(de), " have_p=", have_p,
          " up p<", p_cutoff, " n=",
          sum(!is.na(de$pvalue) & de$pvalue < p_cutoff & de$log2FC > 0))
  writeLines(
    c("本比较只做上调 FC >= 1.25 / 1.5 / 2（先 p < 0.01）。",
      "不做 FC=1，也不做 TopRank。",
      "分层图在 FoldChange/FC_1.25  FC_1.5  FC_2。",
      "全基因 GSEA 在 00_GSEA_all_genes_NOT_FC_or_topN。",
      "神经浸润专项在 Focused_nerve_invasion/（不改全库 p）。"),
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
        emit_subset_analysis(
          comp_name, sub, nm, paste0(comp_name, " | up FC >= ", fc),
          od, de, heat_mat, sample_info, gsea_cache, fc_line = fc
        ),
        error = function(e) log_msg("ERROR subset ", comp_name, " ", nm, ": ", e$message)
      )
    } else {
      utils::write.csv(sub, file.path(od, paste0(nm, "_DE_selected_genes.csv")), row.names = FALSE)
      basic_volcano(de, paste(comp_name, nm), file.path(od, paste0(nm, "_volcano")), fc)
    }
    tryCatch(extract_focus_nerve(od, nm), error = function(e) log_msg("focus nerve: ", e$message))
  }

  if (use_pipe && exists("build_gsea_cache", mode = "function")) {
    tryCatch({
      gsea_cache <- build_gsea_cache(de)
      full_gsea_dir <- file.path(base, "00_GSEA_all_genes_NOT_FC_or_topN")
      dir.create(full_gsea_dir, recursive = TRUE, showWarnings = FALSE)
      writeLines("全基因 GSEA，不是 FC 分层结果。", file.path(full_gsea_dir, "00_README.txt"))
      for (nm in c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Hallmark")) {
        if (exists("plot_gsea_object", mode = "function")) {
          plot_gsea_object(gsea_cache[[nm]], file.path(full_gsea_dir, paste0("allGenes_GSEA_", nm)),
                           paste("GSEA", nm, "|", comp_name, "| ALL genes, NOT FC subset"))
        }
      }
    }, error = function(e) log_msg("full-list GSEA failed: ", e$message))
  }
  invisible(de)
}

# -----------------------------------------------------------------------------
# 9. 对一份空间对象做近/远神经分析（按器官、按病人）
# -----------------------------------------------------------------------------
de_store <- list()
qc_rows <- list()
ligand_rows <- list()
pb_all <- NULL
si_all <- NULL

run_one_spatial <- function(obj, prefix) {
  if (is.null(obj)) return(invisible(NULL))
  sp <- obj$obs
  qc_key <- paste(prefix, paste(unique(sp$patient), collapse = "_"))
  qc_rows[[qc_key]] <<- data.frame(
    assay = prefix,
    n_obs = nrow(sp),
    n_patients = length(unique(sp$patient)),
    organs = paste(sort(unique(sp$organ)), collapse = ","),
    n_tumor = sum(sp$is_tumor),
    n_immune = sum(sp$is_immune),
    n_neural = sum(sp$is_neural),
    n_tumor_near = sum(sp$near_nerve),
    n_tumor_far = sum(sp$far_nerve),
    has_xy = obj$has_xy,
    stringsAsFactors = FALSE
  )
  tryCatch(
    spot_map_plot(sp, paste(prefix, "tumor / immune / neural"),
                  file.path(result_dir, "00_QC_maps", paste0(prefix, "_map"))),
    error = function(e) log_msg("map failed: ", e$message)
  )
  if (!obj$has_xy) {
    log_msg(prefix, " 没有空间坐标，不能做神经浸润邻域。只写 QC。")
    return(invisible(NULL))
  }
  if (sum(sp$is_neural) < 5) {
    log_msg(prefix, " 神经/胶质/Schwann 观测 < 5，跳过邻域 DEG（器官=",
            paste(unique(sp$organ), collapse = ","), "）")
    return(invisible(NULL))
  }

  # 病人 1-vs-1（不合并）
  for (pid in unique(sp$patient)) {
    for (org in unique(sp$organ[sp$patient == pid])) {
      w <- which(sp$patient == pid & sp$organ == org)
      near <- w[sp$near_nerve[w]]
      far <- w[sp$far_nerve[w]]
      if (length(near) < 5 || length(far) < 8) next
      grp <- ifelse(sp$near_nerve[w], "near", ifelse(sp$far_nerve[w], "far", NA))
      keep <- !is.na(grp)
      de_i <- limma_two_group(obj$logmat[w[keep], , drop = FALSE], grp[keep],
                              label = paste(prefix, pid, org))
      if (nrow(de_i) == 0 || all(is.na(de_i$pvalue))) {
        de_i <- mean_fc_two_group(obj$logmat, near, far)
        log_msg(prefix, " ", pid, " ", org, " 单切片无法估计 p，按 FC 分层")
      }
      nm <- paste0(prefix, "_", gsub("[^A-Za-z0-9]+", "_", pid), "_", org, "_tumor_near_nerve_vs_far")
      heat <- t(obj$logmat[c(near, far), , drop = FALSE])
      si <- data.frame(
        sample = colnames(heat),
        group = c(rep("near", length(near)), rep("far", length(far))),
        stringsAsFactors = FALSE
      )
      emit_comparison(nm, de_i, heat, si)
      de_store[[nm]] <<- de_i
    }
  }

  pb <- pseudobulk_near_far(obj$cnts, sp)
  if (!is.null(pb)) {
    if (is.null(pb_all)) {
      pb_all <<- pb$pb
      si_all <<- pb$si
    } else {
      gn <- intersect(rownames(pb_all), rownames(pb$pb))
      if (length(gn) > 50) {
        pb_all <<- cbind(pb_all[gn, , drop = FALSE], pb$pb[gn, , drop = FALSE])
        si_all <<- rbind(si_all, pb$si)
      }
    }
  }

  lig <- present_genes(tumor_ligands$gene, colnames(obj$logmat))
  tum <- which(sp$is_tumor)
  if (length(lig) > 0 && length(tum) >= 15) {
    cors <- sapply(lig, function(g) {
      stats::cor(obj$logmat[tum, g], -pmin(sp$dist_neural[tum], stats::median(sp$dist_neural[is.finite(sp$dist_neural)]) * 4),
                 method = "spearman", use = "complete.obs")
    })
    ligand_rows[[paste(prefix, paste(unique(sp$patient), collapse = ","))]] <<- data.frame(
      assay = prefix,
      organ = paste(sort(unique(sp$organ)), collapse = ","),
      gene = lig,
      mean_in_tumor = colMeans(obj$logmat[tum, lig, drop = FALSE]),
      spearman_vs_near_nerve = as.numeric(cors),
      stringsAsFactors = FALSE
    )
  }
}

# -----------------------------------------------------------------------------
# 10. 载入 Klughammer / BCBM / Wang
# -----------------------------------------------------------------------------
klugh_root <- c(
  file.path(met_dir, "Klughammer_SCP2702"),
  met_dir,
  file.path(met_dir, "SCP2702")
)
klugh_root <- klugh_root[dir.exists(klugh_root)][1]
if (is.na(klugh_root)) klugh_root <- met_dir

h5s <- if (dir.exists(klugh_root)) find_h5ad(klugh_root) else character()
# 空间优先，scRNA 放后
ord <- order(!grepl("slide_seq|merfish|exseq|visium|xenium|codex", basename(h5s), ignore.case = TRUE))
h5s <- h5s[ord]
# 跳过超大 scRNA 直到空间跑完
spatial_h5 <- h5s[grepl("slide_seq|merfish|exseq|visium|xenium", basename(h5s), ignore.case = TRUE)]
scrna_h5 <- setdiff(h5s, spatial_h5)

pairs <- if (dir.exists(klugh_root)) find_tsv_pairs(klugh_root) else list()
loaded_any <- FALSE

for (pr in pairs) {
  blob <- tryCatch(read_counts_annot_pair(pr$counts, pr$annot), error = function(e) {
    log_msg("TSV pair failed: ", e$message)
    NULL
  })
  if (is.null(blob)) next
  loaded_any <- TRUE
  parent <- basename(dirname(dirname(pr$counts)))
  blob$obs$source_file <- basename(dirname(pr$counts))
  blob$obs <- standardize_obs(blob$obs, fallback_patient = parent, fallback_organ = parent)
  obj <- classify_spatial(blob$cnts, blob$obs, assay_name = paste0("tsv_", parent))
  run_one_spatial(obj, paste0("Klughammer_", gsub("[^A-Za-z0-9]+", "_", parent)))
}

if (length(pairs) == 0) {
  for (f in spatial_h5) {
    blob <- tryCatch(read_h5ad(f), error = function(e) {
      log_msg("h5ad failed ", f, ": ", e$message)
      NULL
    })
    if (is.null(blob)) next
    loaded_any <- TRUE
    blob$obs <- standardize_obs(blob$obs, fallback_patient = tools::file_path_sans_ext(basename(f)),
                                fallback_organ = basename(f))
    obj <- classify_spatial(blob$cnts, blob$obs, assay_name = basename(f))
    run_one_spatial(obj, paste0("Klughammer_", tools::file_path_sans_ext(basename(f))))
  }
}

# GSE325935 / 任意 Visium 目录
extra_st <- c(file.path(met_dir, "GSE325935"), file.path(met_dir, "BCBM_ST"))
for (d in extra_st[dir.exists(extra_st)]) {
  log_msg("Scanning extra ST: ", d)
  for (f in find_h5ad(d)) {
    blob <- tryCatch(read_h5ad(f), error = function(e) NULL)
    if (is.null(blob)) next
    loaded_any <- TRUE
    blob$obs <- standardize_obs(blob$obs, fallback_organ = "brain")
    blob$obs$organ[blob$obs$organ == "other"] <- "brain"
    obj <- classify_spatial(blob$cnts, blob$obs, assay_name = basename(f))
    run_one_spatial(obj, paste0("BCBM_", tools::file_path_sans_ext(basename(f))))
  }
  for (pr in find_tsv_pairs(d)) {
    blob <- tryCatch(read_counts_annot_pair(pr$counts, pr$annot), error = function(e) NULL)
    if (is.null(blob)) next
    loaded_any <- TRUE
    blob$obs <- standardize_obs(blob$obs, fallback_organ = "brain")
    blob$obs$organ[blob$obs$organ == "other"] <- "brain"
    obj <- classify_spatial(blob$cnts, blob$obs, "BCBM_tsv")
    run_one_spatial(obj, "BCBM_tsv")
  }
}

# Wang 病理 Nerve：若本机已有 E:/R/Nerve 或 Wang_TNBC，只提示用另一脚本，避免重复 92 例循环
wang_dirs <- c(file.path(met_dir, "Wang_TNBC"), "E:/R/Nerve", "E:\\R\\Nerve")
wang_hit <- wang_dirs[vapply(wang_dirs, function(d) dir.exists(d) && looks_like_wang(d), logical(1))]
if (length(wang_hit) > 0) {
  log_msg("检测到 Wang TNBC ST: ", wang_hit[1],
          " —— 病理 Nerve 金标准请另跑 Wang_ST_nerve_infiltration.R，不要和转移队列混成一张 DEG 表。")
}

if (!loaded_any) {
  log_msg("没有读到 Klughammer / BCBM 空间矩阵。请看 BRCA_met_ST_nerve_DOWNLOAD.txt")
}

# -----------------------------------------------------------------------------
# 11. 按器官伪 bulk（病人作协变量；肺/脑/骨分开，再取共同上调）
# -----------------------------------------------------------------------------
organ_de <- list()
if (!is.null(pb_all) && ncol(pb_all) >= 4) {
  for (org in sort(unique(si_all$organ))) {
    si_o <- si_all[si_all$organ == org, , drop = FALSE]
    if (nrow(si_o) < 4 || length(unique(si_o$group)) < 2) {
      log_msg("器官 ", org, " 伪 bulk 不足，跳过")
      next
    }
    pb_o <- pb_all[, si_o$sample, drop = FALSE]
    de_o <- deseq2_pseudobulk(round(pmax(pb_o, 0)), si_o, paste("organ", org))
    heat <- log2(sweep(pb_o, 2, pmax(colSums(pb_o), 1), "/") * 1e6 + 1)
    emit_comparison(paste0("organ_", org, "_tumor_near_nerve_vs_far"), de_o, heat, si_o)
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
    log_msg("common up lung/brain/bone ", nm, " n=", length(inter),
            " organs=", paste(target_organs, collapse = ","))
    utils::write.csv(
      data.frame(gene = inter, stringsAsFactors = FALSE),
      file.path(result_dir, paste0("common_up_lung_brain_bone_", nm, ".csv")),
      row.names = FALSE
    )
  }
  genes_all <- Reduce(intersect, lapply(organ_de[target_organs], function(d) {
    keep <- d$log2FC > 0
    if (any(!is.na(d$pvalue))) keep <- keep & !is.na(d$pvalue) & d$pvalue < p_cutoff
    d$gene[keep]
  }))
  if (length(genes_all) > 0) {
    tab <- data.frame(gene = genes_all, stringsAsFactors = FALSE)
    for (org in target_organs) {
      d <- organ_de[[org]]
      tab[[paste0("log2FC_", org)]] <- d$log2FC[match(tab$gene, d$gene)]
      tab[[paste0("pvalue_", org)]] <- d$pvalue[match(tab$gene, d$gene)]
    }
    fc_cols <- grep("^log2FC_", names(tab), value = TRUE)
    tab$log2FC <- rowMeans(as.matrix(tab[, fc_cols, drop = FALSE]), na.rm = TRUE)
    pv_cols <- grep("^pvalue_", names(tab), value = TRUE)
    tab$pvalue <- apply(as.matrix(tab[, pv_cols, drop = FALSE]), 1, function(x) {
      if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    })
    tab$AveExpr <- NA_real_
    tab$padj <- tab$pvalue
    tab$padj_BH <- NA_real_
    heat <- matrix(tab$log2FC, ncol = 1, dimnames = list(tab$gene, "mean_log2FC"))
    si <- data.frame(sample = "mean_log2FC", group = "near", stringsAsFactors = FALSE)
    emit_comparison("common_up_lung_brain_bone", tab, heat, si)
  }
} else {
  log_msg("肺/脑/骨中空间神经邻域足够的器官 < 2，跳过三器官共同上调。请补 GSE325935（脑）或 GSE298286（骨）。")
}

# -----------------------------------------------------------------------------
# 12. 候选分子总表
# -----------------------------------------------------------------------------
if (length(qc_rows) > 0) {
  utils::write.csv(do.call(rbind, qc_rows), file.path(log_dir, "spot_class_counts.csv"),
                   row.names = FALSE)
}

lig_all <- if (length(ligand_rows) > 0) do.call(rbind, ligand_rows) else NULL
if (!is.null(lig_all)) {
  agg <- do.call(rbind, lapply(split(lig_all, lig_all$gene), function(d) {
    data.frame(
      gene = d$gene[1],
      n_assays = nrow(d),
      mean_tumor_expr = mean(d$mean_in_tumor, na.rm = TRUE),
      mean_spearman_near_nerve = mean(d$spearman_vs_near_nerve, na.rm = TRUE),
      n_pos_cor = sum(d$spearman_vs_near_nerve > 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
} else {
  agg <- data.frame(gene = tumor_ligands$gene, n_assays = 0,
                    mean_tumor_expr = NA, mean_spearman_near_nerve = NA, n_pos_cor = 0)
}

cand <- merge(tumor_ligands, agg, by = "gene", all.x = TRUE)
# 用各器官伪 bulk 注释候选
for (org in names(organ_de)) {
  d <- organ_de[[org]]
  cand[[paste0(org, "_near_log2FC")]] <- d$log2FC[match(cand$gene, d$gene)]
  cand[[paste0(org, "_near_pvalue")]] <- d$pvalue[match(cand$gene, d$gene)]
}
cand$up_spatial_near <- FALSE
fc_org <- grep("_near_log2FC$", names(cand), value = TRUE)
if (length(fc_org) > 0) {
  up_any <- rep(FALSE, nrow(cand))
  for (fcn in fc_org) {
    org <- sub("_near_log2FC$", "", fcn)
    pvn <- paste0(org, "_near_pvalue")
    fc <- cand[[fcn]]
    pv <- if (pvn %in% names(cand)) cand[[pvn]] else rep(NA_real_, nrow(cand))
    up_any <- up_any | (!is.na(fc) & fc > 0 & !is.na(pv) & pv < p_cutoff)
  }
  cand$up_spatial_near <- up_any
}
cand$priority <-
  1 * (!is.na(cand$mean_spearman_near_nerve) & cand$mean_spearman_near_nerve > 0.1) +
  3 * isTRUE(cand$up_spatial_near)
cand <- cand[order(-cand$priority, -cand$mean_spearman_near_nerve), ]
utils::write.csv(cand, file.path(result_dir, "01_CANDIDATE_MOLECULES_tumor_to_nerve.csv"),
                 row.names = FALSE)
if (has_pkg("writexl")) {
  tryCatch(writexl::write_xlsx(cand, file.path(result_dir, "01_CANDIDATE_MOLECULES_tumor_to_nerve.xlsx")),
           error = function(e) NULL)
}

# -----------------------------------------------------------------------------
# 13. 方案说明
# -----------------------------------------------------------------------------
protocol <- c(
  "============================================================",
  "乳腺癌转移空间转录组：肿瘤高表达基因与神经浸润",
  "数据目录: ", met_dir,
  "============================================================",
  "",
  "公开库没有一份 ST 同时有：多病人 + 肺/脑/骨 + 肿瘤细胞 + 免疫 + 神经。",
  "主队列 Klughammer Nat Med 2024（SCP2702）：60 病人，多器官，恶性+免疫，15 例空间。",
  "脑/肺空间各约 1 例，神经细胞不够 → 补 GSE325935（人脑转移 ST）。",
  "病理 Nerve 金标准：Wang 2024 TNBC（原发灶，不要冒充转移）。",
  "下载清单：BRCA_met_ST_nerve_DOWNLOAD.txt",
  "",
  "【定义】",
  "  肿瘤类 = 恶性上皮（标签或 EPCAM/KRT 签名），不是免疫/成纤维",
  "  神经   = Schwann/Nerve 或中枢胶质（GFAP/MBP/PLP1 等）；神经元常极少，缺则写日志",
  "  浸润   = 肿瘤–神经空间邻域（近 25% 距离 vs 远 25%），不是 bulk 共表达",
  "",
  "【统计】先 p < 0.01，再上调 FC >= 1.25 / 1.5 / 2；不做 FC=1、不做 topN。",
  "  病人、器官不合并后再做 DEG。单切片不伪造 p。",
  "  关联 ≠ 因果：候选须肿瘤侧高表达 + 空间邻接 + 神经侧受体。",
  "",
  "【先看】",
  "  results/01_CANDIDATE_MOLECULES_tumor_to_nerve.csv",
  "  results/organ_brain_tumor_near_nerve_vs_far/",
  "  results/organ_lung_tumor_near_nerve_vs_far/",
  "  results/organ_bone_tumor_near_nerve_vs_far/",
  "  results/common_up_lung_brain_bone/   （器官都够才有）",
  "  Focused_nerve_invasion/  保留原始 p 与 genome_wide_rank",
  "",
  "ORA 在 GO/ Pathway/ KEGG/，文件名 ORA_ 开头。",
  "GSEA 在 GSEA/；全基因 GSEA 只在 00_GSEA_all_genes_NOT_FC_or_topN。"
)
writeLines(protocol, file.path(result_dir, "00_PROTOCOL_神经浸润分析方案.txt"))
log_msg("Done. Open: ", file.path(result_dir, "00_PROTOCOL_神经浸润分析方案.txt"))
log_msg("and: ", file.path(result_dir, "01_CANDIDATE_MOLECULES_tumor_to_nerve.csv"))
