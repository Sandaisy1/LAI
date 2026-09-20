#!/usr/bin/env Rscript
# =============================================================================
# GSE110590_Human_breast.R
# Siegel 2018 UNC RAP：配对原发 vs 肺/骨转移 RNA-seq（GSE110590）
#
# 数据目录（脚本与结果默认同目录）：
#   E:/R/Human breast cancer/GSE110590
#
# 四个问题（只分析「原发相对配对转移更低」= 转移/原发 FC > 1）：
#   1) 配对原发 vs 肺转移、配对原发 vs 骨转移
#   2) 器官特异：只肺或只骨（配对 DE 取差集；本队列无「只转骨」患者）
#   3) 轴突导向 / 施旺细胞 / 神经营养 三套签名给原发打神经浸润分
#   4) 三种神经浸润分别与肺转移的关系
#
# 阈值：p < 0.05，FC > 1 与 FC > 1.25。不做 top50–300。
# 配对：同一患者 1 个原发对 1 个该器官转移（多灶则对该器官取均值）。
# 热图列顺序 Primary_i, Met_i。
#
# 本数据集：RSEM 上分位数标准化 log2（约 17k 基因 × 83 样品）。
#   有原发的肺对约 10 例；骨对约 4 例；A8 无原发列，不能配对。
# 不要改 TG_RNAseq_*.R。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")
options(clusterProfiler.download.method = "auto")

# -----------------------------------------------------------------------------
# 0. 包
# -----------------------------------------------------------------------------
cran_required <- c(
  "dplyr", "tidyr", "tibble", "stringr", "ggplot2", "ggrepel",
  "pheatmap", "RColorBrewer", "matrixStats", "writexl"
)
cran_optional <- c("ggvenn", "GSVA")
bioc_required <- c(
  "limma", "clusterProfiler", "org.Hs.eg.db",
  "enrichplot", "AnnotationDbi", "fgsea"
)
bioc_optional <- c("ReactomePA", "msigdbr", "pathview")

install_if_missing <- function(pkgs, bioc = FALSE, required = TRUE) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) == 0) return(invisible(TRUE))
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    tryCatch(
      BiocManager::install(miss, update = FALSE, ask = FALSE),
      error = function(e) message("Bioconductor install failed: ", e$message)
    )
  } else {
    tryCatch(
      install.packages(miss, repos = "https://cloud.r-project.org"),
      error = function(e) message("CRAN install failed: ", e$message)
    )
  }
  still <- miss[!vapply(miss, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still) > 0 && required) {
    stop("缺少必需 R 包: ", paste(still, collapse = ", "))
  }
  if (length(still) > 0) {
    message("可选包未安装，相关分析将跳过: ", paste(still, collapse = ", "))
  }
  invisible(TRUE)
}

install_if_missing(cran_required, bioc = FALSE, required = TRUE)
install_if_missing(cran_optional, bioc = FALSE, required = FALSE)
install_if_missing(bioc_required, bioc = TRUE, required = TRUE)
install_if_missing(bioc_optional, bioc = TRUE, required = FALSE)

safe_library <- function(pkgs) {
  for (p in pkgs) {
    if (requireNamespace(p, quietly = TRUE)) {
      suppressPackageStartupMessages(library(p, character.only = TRUE))
    }
  }
}
safe_library(c(cran_required, cran_optional, bioc_required, bioc_optional))
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

# -----------------------------------------------------------------------------
# 1. 路径
# -----------------------------------------------------------------------------
resolve_project_dir <- function() {
  env_dir <- Sys.getenv("GSE110590_DIR", unset = "")
  script_dir <- tryCatch({
    of <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(of)) dirname(normalizePath(sub("^--file=", "", of[1]))) else NA_character_
  }, error = function(e) NA_character_)
  candidates <- c(
    env_dir,
    "E:/R/Human breast cancer/GSE110590",
    "E:\\R\\Human breast cancer\\GSE110590",
    script_dir,
    file.path(getwd(), "GSE110590"),
    getwd()
  )
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  hit_file <- function(d) {
    if (!dir.exists(d)) return(FALSE)
    ff <- list.files(d, full.names = FALSE, ignore.case = TRUE)
    any(grepl("GSE110590|RAP_A16|log2\\.sne|series_matrix|GPL11154|GPL16791", ff, ignore.case = TRUE))
  }
  for (d in candidates) {
    if (hit_file(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

project_dir <- resolve_project_dir()
result_dir  <- file.path(project_dir, "results_GSE110590_Human_breast")
log_dir     <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("GSE110590_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_msg("Working directory: ", project_dir)
log_msg("Results: ", result_dir)

p_cutoff  <- 0.05
fc_cutoffs <- c("FC_1" = 1, "FC_1.25" = 1.25)

# -----------------------------------------------------------------------------
# 2. 读 GEO 文件
# -----------------------------------------------------------------------------
find_first <- function(patterns) {
  ff <- list.files(project_dir, full.names = TRUE, recursive = TRUE)
  ff <- ff[!grepl("results_GSE110590", ff, ignore.case = TRUE)]
  base <- basename(ff)
  for (pat in patterns) {
    hit <- ff[grepl(pat, base, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

find_all <- function(patterns) {
  ff <- list.files(project_dir, full.names = TRUE, recursive = TRUE)
  ff <- ff[!grepl("results_GSE110590", ff, ignore.case = TRUE)]
  base <- basename(ff)
  hit <- character(0)
  for (pat in patterns) {
    hit <- c(hit, ff[grepl(pat, base, ignore.case = TRUE)])
  }
  unique(hit)
}

strip_quotes <- function(x) gsub('^\"|\"$', "", x)

parse_series_matrix_pheno <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  titles <- gsms <- source <- NULL
  chars <- list()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "!series_matrix_table_begin")) break
    if (startsWith(line, "!Sample_title\t")) {
      titles <- strip_quotes(strsplit(line, "\t", fixed = TRUE)[[1]][-1])
    } else if (startsWith(line, "!Sample_geo_accession\t")) {
      gsms <- strip_quotes(strsplit(line, "\t", fixed = TRUE)[[1]][-1])
    } else if (startsWith(line, "!Sample_source_name_ch1\t")) {
      source <- strip_quotes(strsplit(line, "\t", fixed = TRUE)[[1]][-1])
    } else if (startsWith(line, "!Sample_characteristics_ch1\t")) {
      vals <- strip_quotes(strsplit(line, "\t", fixed = TRUE)[[1]][-1])
      field <- sub(":.*$", "", vals[1])
      chars[[field]] <- sub("^[^:]+:\\s*", "", vals)
    }
  }
  if (is.null(titles)) stop("series matrix 没有 !Sample_title: ", path)
  n <- length(titles)
  pick <- function(nm, default = NA_character_) {
    if (nm %in% names(chars)) chars[[nm]] else rep(default, n)
  }
  data.frame(
    title = titles,
    gsm = if (is.null(gsms)) rep(NA_character_, n) else gsms,
    patient_geo = pick("patient id"),
    tumor_location = pick("tumor location"),
    source = if (is.null(source)) rep(NA_character_, n) else source,
    stringsAsFactors = FALSE
  )
}

is_primary_name <- function(nm) {
  u <- toupper(nm)
  if (grepl("MET", u) && !grepl("PRIMT|PRIMARY|PTCORE", u)) return(FALSE)
  grepl("PRIMT|PRIMARY|PTcore|(^|-)PT(-|$)|PT-FFPE", nm, ignore.case = TRUE)
}

infer_kind_tissue <- function(nm, loc = NA_character_, src = NA_character_) {
  u <- toupper(nm)
  locu <- toupper(trimws(ifelse(is.na(loc), "", loc)))
  srcu <- toupper(trimws(ifelse(is.na(src), "", src)))
  if (identical(srcu, "BREAST CANCER")) srcu <- ""
  if (is_primary_name(nm)) {
    return(c(kind = "primary", tissue = "Breast"))
  }
  name_map <- list(
    Lung = "LUNG|(^|-)(RLL|LLL|LUL|RUL)(-|$)",
    Bone = "RIB|SPIN|SKULL|\\bBONE\\b",
    Brain = "BRAIN|CELEB|DURA|OCCI",
    Liver = "LIV",
    LN = "LN-|AX-LN|LYMPH|SUBCARINAL|MEDIASTN",
    Adrenal = "ADREN",
    Pleura = "PLEURA",
    Pancreas = "PANCREAS",
    Kidney = "KIDNY|KIDNEY",
    Ovary = "OVA",
    Skin = "SKIN",
    SoftTissue = "SOFTTISSUE",
    ChestWall = "CHESTWALL|CHESTMET|CHEST-?MET"
  )
  for (org in names(name_map)) {
    if (grepl(name_map[[org]], u)) {
      return(c(kind = "met", tissue = org))
    }
  }
  loc_map <- c(
    PRIMARY = "Breast", LUNG = "Lung", RIB = "Bone", SPINAL = "Bone",
    SKULL = "Bone", BONE = "Bone", BRAIN = "Brain", LIVER = "Liver",
    ADRENAL = "Adrenal", PLEURA = "Pleura", OVARY = "Ovary", SKIN = "Skin",
    PANCREAS = "Pancreas", SOFTTISSUE = "SoftTissue",
    `AUXILLARY LYMPH NODE` = "LN", `LYMPH NODE` = "LN",
    `SUBCARINAL LYMPH NODE` = "LN"
  )
  hit <- loc_map[locu]
  if (!is.na(hit) && nzchar(hit)) {
    kind <- if (identical(unname(hit), "Breast")) "primary" else "met"
    return(c(kind = kind, tissue = unname(hit)))
  }
  hit2 <- loc_map[srcu]
  if (!is.na(hit2) && nzchar(hit2)) {
    kind <- if (identical(unname(hit2), "Breast")) "primary" else "met"
    return(c(kind = kind, tissue = unname(hit2)))
  }
  if (grepl("MET", u)) return(c(kind = "met", tissue = "Other"))
  c(kind = "other", tissue = "Other")
}

match_expr_to_title <- function(cn, titles) {
  if (cn %in% titles) return(cn)
  # TSV 有 A8-OVA-MET-RNA-RNA，GEO title 是 A8-OVA-MET-RNA
  stripped <- sub("(-RNA)+$", "", cn)
  tstrip <- sub("(-RNA)+$", "", titles)
  hit <- titles[tstrip == stripped]
  if (length(hit) == 1) return(hit[1])
  hit <- titles[toupper(titles) == toupper(cn)]
  if (length(hit) == 1) return(hit[1])
  NA_character_
}

read_rap_matrix <- function(path) {
  log_msg("Reading RAP log2 matrix: ", path)
  df <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    utils::read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  raw_id <- as.character(df[[1]])
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  parts <- strsplit(raw_id, "|", fixed = TRUE)
  n <- length(parts)
  symbol <- character(n)
  entrez <- character(n)
  for (i in seq_len(n)) {
    p <- parts[[i]]
    sym <- if (length(p) >= 1) p[1] else ""
    ent <- if (length(p) >= 2) p[2] else ""
    entrez[i] <- ent
    if (sym %in% c("?", "", "NA") || is.na(sym)) {
      symbol[i] <- if (nzchar(ent)) paste0("ENTREZID_", ent) else paste0("row", i)
    } else {
      symbol[i] <- sym
    }
  }
  if (any(duplicated(symbol))) {
    means <- rowMeans(mat, na.rm = TRUE)
    ord <- order(means, decreasing = TRUE)
    mat <- mat[ord, , drop = FALSE]
    symbol <- symbol[ord]
    entrez <- entrez[ord]
    keep <- !duplicated(symbol)
    mat <- mat[keep, , drop = FALSE]
    symbol <- symbol[keep]
    entrez <- entrez[keep]
  }
  rownames(mat) <- symbol
  attr(mat, "entrez") <- setNames(entrez, symbol)
  mat
}

expr_path <- find_first(c(
  "GSE110590_RAP_A16_log2\\.sne\\.tsv(\\.gz)?$",
  "RAP_A16_log2",
  "log2\\.sne\\.tsv"
))
pheno_paths <- find_all(c(
  "GSE110590-GPL11154_series_matrix",
  "GSE110590-GPL16791_series_matrix",
  "GSE110590_GPL11154",
  "GSE110590_GPL16791",
  "series_matrix"
))
if (is.na(expr_path) || !file.exists(expr_path)) {
  stop("找不到 GSE110590_RAP_A16_log2.sne.tsv.gz，请放在: ", project_dir)
}

raw <- read_rap_matrix(expr_path)
gene_entrez <- attr(raw, "entrez")
log_msg("RAP matrix: ", nrow(raw), " genes x ", ncol(raw), " samples from ", basename(expr_path))

pheno_geo <- NULL
if (length(pheno_paths) > 0) {
  pheno_geo <- dplyr::bind_rows(lapply(pheno_paths, parse_series_matrix_pheno))
  pheno_geo <- pheno_geo[!duplicated(pheno_geo$title), ]
  log_msg("Series matrix samples: ", nrow(pheno_geo), " from ",
          paste(basename(pheno_paths), collapse = ", "))
} else {
  log_msg("WARNING: 未找到 series matrix，仅用列名推断器官")
}

cn <- colnames(raw)
title_hit <- if (is.null(pheno_geo)) rep(NA_character_, length(cn)) else {
  vapply(cn, match_expr_to_title, character(1), titles = pheno_geo$title)
}
geo_row <- if (is.null(pheno_geo)) {
  list(patient_geo = rep(NA_character_, length(cn)),
       tumor_location = rep(NA_character_, length(cn)),
       source = rep(NA_character_, length(cn)),
       gsm = rep(NA_character_, length(cn)),
       title = cn)
} else {
  idx <- match(title_hit, pheno_geo$title)
  list(
    patient_geo = pheno_geo$patient_geo[idx],
    tumor_location = pheno_geo$tumor_location[idx],
    source = pheno_geo$source[idx],
    gsm = pheno_geo$gsm[idx],
    title = ifelse(is.na(title_hit), cn, title_hit)
  )
}

patient_from_name <- sub("^((A[0-9]+)).*$", "\\1", cn)
patient <- ifelse(!is.na(geo_row$patient_geo) & nzchar(geo_row$patient_geo),
                  geo_row$patient_geo, patient_from_name)

kt <- t(mapply(infer_kind_tissue, cn, geo_row$tumor_location, geo_row$source,
               USE.NAMES = FALSE))
pheno <- data.frame(
  sample = cn,
  gsm = geo_row$gsm,
  title = geo_row$title,
  patient = patient,
  kind = as.character(kt[, "kind"]),
  tissue = as.character(kt[, "tissue"]),
  tumor_location = geo_row$tumor_location,
  source = geo_row$source,
  stringsAsFactors = FALSE
)
rownames(pheno) <- pheno$sample
log_msg("kind: ", paste(names(table(pheno$kind)), table(pheno$kind), sep = "=", collapse = "; "))
log_msg("tissue: ", paste(names(table(pheno$tissue)), table(pheno$tissue), sep = "=", collapse = "; "))

storage.mode(raw) <- "double"
raw[!is.finite(raw)] <- NA
mx <- stats::quantile(as.numeric(raw), 0.99, na.rm = TRUE)
if (is.finite(mx) && mx > 40) {
  log_msg("Values look like raw counts (q99=", signif(mx, 3), "); log2(x+1)")
  logmat <- log2(pmax(raw, 0) + 1)
} else {
  log_msg("Values look log-scale (q99=", signif(mx, 3), "); use as-is (RSEM UQN log2)")
  logmat <- raw
}
keep_gene <- rowSums(is.finite(logmat)) >= 3 &
  matrixStats::rowSds(logmat, na.rm = TRUE) > 1e-6
log_msg("Low-variance filter: keep ", sum(keep_gene), " / ", nrow(logmat))
logmat <- logmat[keep_gene, , drop = FALSE]
if (!is.null(gene_entrez)) gene_entrez <- gene_entrez[rownames(logmat)]
pheno_use <- pheno
heat_mat <- t(scale(t(logmat)))
heat_mat[!is.finite(heat_mat)] <- 0
utils::write.csv(pheno_use, file.path(result_dir, "00_logs", "sample_annotation.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. 一一对应：同一患者 1 原发 + 该器官转移（多灶取均值）
# -----------------------------------------------------------------------------
pick_primary <- function(rows) {
  if (nrow(rows) == 1) return(rows$sample[1])
  rows$sample[1]
}

build_pairs <- function(organ, organ_labels) {
  mets <- pheno_use[pheno_use$kind == "met" & pheno_use$tissue %in% organ_labels, , drop = FALSE]
  prim <- pheno_use[pheno_use$kind == "primary", , drop = FALSE]
  out <- list()
  for (pid in intersect(unique(mets$patient), unique(prim$patient))) {
    p_rows <- prim[prim$patient == pid, , drop = FALSE]
    m_rows <- mets[mets$patient == pid, , drop = FALSE]
    p_id <- pick_primary(p_rows)
    m_ids <- m_rows$sample
    out[[length(out) + 1]] <- data.frame(
      pair_id = pid,
      patient = pid,
      organ = organ,
      primary = p_id,
      metastasis = if (length(m_ids) == 1) m_ids else paste(m_ids, collapse = ";"),
      n_met_samples = length(m_ids),
      primary_title = p_rows$title[match(p_id, p_rows$sample)],
      met_title = paste(m_rows$title, collapse = ";"),
      met_tissue = paste(unique(m_rows$tissue), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) {
    return(data.frame(
      pair_id = character(), patient = character(), organ = character(),
      primary = character(), metastasis = character(), n_met_samples = integer(),
      primary_title = character(), met_title = character(),
      met_tissue = character(), stringsAsFactors = FALSE
    ))
  }
  dplyr::bind_rows(out)
}

pairs_lung <- build_pairs("Lung", "Lung")
pairs_bone <- build_pairs("Bone", c("Bone", "Rib", "Spine", "Spinal", "Skull"))
pair_dir <- file.path(result_dir, "00_sample_pairing")
dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(pairs_lung, file.path(pair_dir, "pairs_primary_vs_lung_1to1.csv"), row.names = FALSE)
utils::write.csv(pairs_bone, file.path(pair_dir, "pairs_primary_vs_bone_1to1.csv"), row.names = FALSE)
n_prim <- sum(pheno_use$kind == "primary")
log_msg("Primaries: ", n_prim, " patients ", paste(sort(unique(pheno_use$patient[pheno_use$kind == "primary"])), collapse = ", "))
log_msg("Paired primary-lung: ", nrow(pairs_lung), " patients")
log_msg("Paired primary-bone: ", nrow(pairs_bone), " patients")
if (nrow(pairs_lung) > 0) log_msg("Lung pairs: ", paste(pairs_lung$patient, collapse = ", "))
if (nrow(pairs_bone) > 0) log_msg("Bone pairs: ", paste(pairs_bone$patient, collapse = ", "))
no_pt <- setdiff(unique(pheno_use$patient[pheno_use$kind == "met"]),
                 unique(pheno_use$patient[pheno_use$kind == "primary"]))
if (length(no_pt) > 0) {
  log_msg("Patients with met but no primary in matrix (cannot 1-to-1): ", paste(no_pt, collapse = ", "))
  writeLines(
    c("These patients have metastases in GSE110590 but no primary column, so they are excluded from paired DE:",
      paste(no_pt, collapse = ", "),
      "A8 is the typical example (lung/spine/liver/ovary mets, no PRIM/PT column)."),
    file.path(pair_dir, "UNPAIRED_MET_ONLY_PATIENTS.txt")
  )
}

# 每个配对患者：原发向量 vs 该器官转移（多灶则均值）
pair_vectors <- function(pairs) {
  prim_mat <- matrix(NA_real_, nrow = nrow(logmat), ncol = nrow(pairs),
                     dimnames = list(rownames(logmat), pairs$patient))
  met_mat <- prim_mat
  for (i in seq_len(nrow(pairs))) {
    p <- pairs$primary[i]
    m <- strsplit(pairs$metastasis[i], ";", fixed = TRUE)[[1]]
    m <- intersect(m, colnames(logmat))
    prim_mat[, i] <- logmat[, p]
    met_mat[, i] <- if (length(m) == 1) logmat[, m] else rowMeans(logmat[, m, drop = FALSE])
  }
  list(primary = prim_mat, met = met_mat, patients = pairs$patient)
}

# -----------------------------------------------------------------------------
# 5. 配对 limma：转移 vs 原发；FC = 转移/原发（>1 表示原发更低）
# -----------------------------------------------------------------------------
paired_limma <- function(pairs, label) {
  if (nrow(pairs) < 2) {
    log_msg(label, ": fewer than 2 pairs, skip limma")
    return(NULL)
  }
  if (nrow(pairs) < 3) {
    log_msg("WARNING ", label, ": only ", nrow(pairs),
            " pairs. p-values have very few residual df; interpret cautiously.")
  }
  vec <- pair_vectors(pairs)
  expr <- cbind(vec$primary, vec$met)
  colnames(expr) <- c(paste0(vec$patients, "_Primary"), paste0(vec$patients, "_Met"))
  patient <- factor(c(vec$patients, vec$patients))
  tissue  <- factor(c(rep("Primary", ncol(vec$primary)), rep("Met", ncol(vec$met))),
                    levels = c("Primary", "Met"))
  design <- stats::model.matrix(~ patient + tissue)
  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit, coef = "tissueMet", number = Inf, sort.by = "none")
  fc_pairs <- 2^(vec$met - vec$primary)
  colnames(fc_pairs) <- paste0("FC_", vec$patients)
  de <- data.frame(
    gene = rownames(tt),
    log2FC = tt$logFC,
    FC = 2^tt$logFC,
    AveExpr = tt$AveExpr,
    t = tt$t,
    pvalue = tt$P.Value,
    padj = tt$adj.P.Val,
    n_pairs = nrow(pairs),
    n_pairs_FC_gt_1 = rowSums(fc_pairs > 1, na.rm = TRUE),
    n_pairs_FC_gt_1.25 = rowSums(fc_pairs > 1.25, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  de <- cbind(de, as.data.frame(fc_pairs))
  de <- de[order(de$pvalue, -de$log2FC), ]
  rownames(de) <- NULL
  attr(de, "pair_expr") <- vec
  attr(de, "pairs") <- pairs
  attr(de, "label") <- label
  de
}

select_up_in_met <- function(de, fc_min, p_min = p_cutoff) {
  if (is.null(de) || !is.data.frame(de) || nrow(de) == 0) {
    return(data.frame())
  }
  de[!is.na(de$pvalue) & de$pvalue < p_min & de$FC > fc_min, , drop = FALSE]
}

# -----------------------------------------------------------------------------
# 6. 绘图 / ORA / GSEA
# -----------------------------------------------------------------------------
save_gg <- function(plot, path_stub, width = 8, height = 6) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".pdf"), plot, width = width, height = height),
           error = function(e) log_msg("pdf ggsave failed: ", e$message))
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".png"), plot, width = width, height = height, dpi = 300),
           error = function(e) log_msg("png ggsave failed: ", e$message))
}

note_empty <- function(stub, msg) writeLines(msg, paste0(stub, "_EMPTY.txt"))

try_save_plot <- function(fun, stub, width = 9, height = 7) {
  p <- tryCatch(fun(), error = function(e) {
    log_msg("Plot failed (", basename(stub), "): ", e$message)
    NULL
  })
  if (is.null(p)) return(invisible(FALSE))
  save_gg(p, stub, width = width, height = height)
}

plot_volcano <- function(de, highlight, title, outfile, fc_line = 1) {
  df <- de
  df$y <- -log10(pmax(df$pvalue, 1e-300))
  df$set <- ifelse(df$gene %in% highlight, "selected", "other")
  df$label <- ifelse(df$gene %in% utils::head(highlight, 15), df$gene, NA)
  lfc_line <- log2(fc_line)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = y, color = set)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.4) +
    ggplot2::scale_color_manual(values = c(other = "grey70", selected = "#D62828")) +
    ggplot2::geom_vline(xintercept = c(-lfc_line, lfc_line), linetype = 2, color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = 2, color = "grey40") +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = title,
      x = "log2FC (matched metastasis / primary); >0 = lower in primary",
      y = "-log10(p value)", color = NULL
    )
  save_gg(p, outfile)
}

plot_de_bar <- function(sub, title, outfile) {
  if (nrow(sub) == 0) return(invisible(NULL))
  df <- sub[order(sub$log2FC, decreasing = TRUE), , drop = FALSE]
  if (nrow(df) > 60) df <- rbind(utils::head(df, 30), utils::tail(df, 30))
  df$gene <- factor(df$gene, levels = rev(unique(df$gene)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = gene, y = log2FC)) +
    ggplot2::geom_col(fill = "#D62828", width = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title, x = NULL, y = "log2FC (met / primary)")
  save_gg(p, outfile, width = 8, height = max(5, min(16, 0.22 * nrow(df) + 2)))
}

group_heatmap <- function(mat, samples, groups, genes, title, outfile) {
  genes <- intersect(genes, rownames(mat))
  samples <- intersect(samples, colnames(mat))
  if (length(genes) > 200) genes <- genes[seq_len(200)]
  if (length(genes) < 2 || length(samples) < 2) {
    log_msg("Group heatmap skipped: ", title)
    return(invisible(NULL))
  }
  sub <- mat[genes, samples, drop = FALSE]
  ann <- data.frame(Group = groups[match(samples, names(groups))], row.names = samples)
  default_pal <- c(
    low = "#4C78A8", high = "#D62828",
    Primary = "#4C78A8", Metastasis = "#D62828",
    lung_only = "#D62828", pleura_only = "#4C78A8", bone_only = "#2A9D8F",
    lung_and_bone = "#E9C46A"
  )
  ug <- unique(as.character(ann$Group))
  pal_g <- default_pal[names(default_pal) %in% ug]
  missing <- setdiff(ug, names(pal_g))
  if (length(missing) > 0) {
    extra_cols <- RColorBrewer::brewer.pal(max(3, length(missing)), "Set2")[seq_along(missing)]
    pal_g <- c(pal_g, stats::setNames(extra_cols, missing))
  }
  pal <- list(Group = pal_g)
  draw <- function() {
    pheatmap::pheatmap(
      sub, scale = "row", annotation_col = ann, annotation_colors = pal,
      cluster_cols = TRUE, cluster_rows = TRUE,
      show_rownames = nrow(sub) <= 80, fontsize_row = 6, main = title,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100)
    )
  }
  grDevices::pdf(paste0(outfile, ".pdf"), width = 8, height = max(6, min(18, 0.18 * nrow(sub) + 3)))
  on.exit({
    while (grDevices::dev.cur() > 1) grDevices::dev.off()
  }, add = TRUE)
  draw()
  grDevices::dev.off()
  grDevices::png(paste0(outfile, ".png"), width = 2400, height = max(1800, 40 * nrow(sub) + 400), res = 300)
  draw()
  grDevices::dev.off()
}

paired_heatmap <- function(de_obj, genes, title, outfile) {
  vec <- attr(de_obj, "pair_expr")
  genes <- intersect(genes, rownames(vec$primary))
  if (length(genes) > 200) genes <- genes[seq_len(200)]
  if (length(genes) < 2) {
    log_msg("Heatmap skipped (<2 genes): ", title)
    return(invisible(NULL))
  }
  mat <- matrix(NA_real_, nrow = length(genes), ncol = 2 * length(vec$patients),
                dimnames = list(genes, NULL))
  cn <- character(0)
  ann <- data.frame(Patient = character(0), Tissue = character(0), stringsAsFactors = FALSE)
  col_i <- 1
  for (j in seq_along(vec$patients)) {
    mat[, col_i] <- vec$primary[genes, j]
    mat[, col_i + 1] <- vec$met[genes, j]
    cn <- c(cn, paste0(vec$patients[j], "_P"), paste0(vec$patients[j], "_M"))
    ann <- rbind(ann, data.frame(
      Patient = vec$patients[j], Tissue = "Primary", stringsAsFactors = FALSE
    ))
    ann <- rbind(ann, data.frame(
      Patient = vec$patients[j], Tissue = "Metastasis", stringsAsFactors = FALSE
    ))
    col_i <- col_i + 2
  }
  colnames(mat) <- cn
  rownames(ann) <- cn
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  pal <- list(
    Tissue = c(Primary = "#4C78A8", Metastasis = "#D62828")
  )
  draw <- function() {
    pheatmap::pheatmap(
      z, scale = "none", annotation_col = ann, annotation_colors = pal,
      cluster_cols = FALSE, cluster_rows = TRUE,
      show_rownames = nrow(z) <= 80, fontsize_row = 6, main = title,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100)
    )
  }
  grDevices::pdf(paste0(outfile, ".pdf"), width = max(8, 0.45 * ncol(z) + 4),
                 height = max(6, min(18, 0.18 * nrow(z) + 3)))
  on.exit({
    while (grDevices::dev.cur() > 1) grDevices::dev.off()
  }, add = TRUE)
  draw()
  grDevices::dev.off()
  grDevices::png(paste0(outfile, ".png"), width = max(1800, 80 * ncol(z)),
                 height = max(1800, 40 * nrow(z) + 400), res = 300)
  draw()
  grDevices::dev.off()
}

pick_official_symbol <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.[0-9]+$", "", x)
  if (length(x) != 1 || is.na(x) || x %in% c("", "-", ".", "NA")) return(NA_character_)
  parts <- unlist(strsplit(x, "[,;|/]+"))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) return(NA_character_)
  parts[1]
}

map_to_entrez <- function(symbols) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  if (length(symbols) == 0) return(data.frame(gene = character(), entrez = character()))
  ent <- if (exists("gene_entrez", inherits = TRUE) && !is.null(gene_entrez)) {
    unname(gene_entrez[symbols])
  } else {
    rep(NA_character_, length(symbols))
  }
  out <- data.frame(gene = symbols, entrez = ent, stringsAsFactors = FALSE)
  need <- is.na(out$entrez) | !nzchar(out$entrez)
  if (any(need)) {
    query <- unique(vapply(out$gene[need], pick_official_symbol, character(1), USE.NAMES = FALSE))
    query <- query[!is.na(query) & nzchar(query) & !startsWith(query, "ENTREZID_")]
    if (length(query) > 0) {
      m <- tryCatch(
        clusterProfiler::bitr(query, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
        error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
      )
      if (nrow(m) > 0) {
        m <- m[!duplicated(m[[1]]), ]
        out$entrez[need] <- m[[2]][match(out$gene[need], m[[1]])]
      }
    }
    still <- is.na(out$entrez) | !nzchar(out$entrez)
    from_name <- sub("^ENTREZID_", "", out$gene[still])
    ok <- grepl("^[0-9]+$", from_name)
    out$entrez[still][ok] <- from_name[ok]
  }
  out <- out[!is.na(out$entrez) & nzchar(out$entrez), , drop = FALSE]
  out[!duplicated(out$gene), , drop = FALSE]
}

ranked_entrez <- function(de) {
  mp <- map_to_entrez(de$gene)
  de2 <- merge(de, mp, by = "gene")
  de2 <- de2[!is.na(de2$entrez) & !is.na(de2$log2FC), ]
  de2 <- de2[order(abs(de2$log2FC), decreasing = TRUE), ]
  de2 <- de2[!duplicated(de2$entrez), ]
  stats <- de2$log2FC
  names(stats) <- de2$entrez
  sort(stats, decreasing = TRUE)
}

enrich_or_relax <- function(strict_fun, relax_fun, label) {
  obj <- tryCatch(strict_fun(), error = function(e) {
    log_msg(label, " strict failed: ", e$message)
    NULL
  })
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    attr(obj, "relaxed") <- FALSE
    return(obj)
  }
  obj2 <- tryCatch(relax_fun(), error = function(e) {
    log_msg(label, " relaxed failed: ", e$message)
    NULL
  })
  if (!is.null(obj2)) attr(obj2, "relaxed") <- TRUE
  obj2
}

plot_ora_object <- function(x, stub, title, fold_change = NULL) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) {
    note_empty(stub, "no enrichment terms")
    return(invisible(NULL))
  }
  df <- as.data.frame(x)
  utils::write.csv(df, paste0(stub, ".csv"), row.names = FALSE)
  nshow <- min(15, nrow(df))
  try_save_plot(function() enrichplot::dotplot(x, showCategory = nshow) + ggplot2::ggtitle(title),
                paste0(stub, "_dotplot"), 9, 7)
  try_save_plot(function() enrichplot::barplot(x, showCategory = nshow) + ggplot2::ggtitle(title),
                paste0(stub, "_barplot"), 9, 7)
}

run_ora_plots <- function(genes, de_sub, outdir, label, tag) {
  go_dir <- file.path(outdir, "GO")
  pw_dir <- file.path(outdir, "Pathway")
  kg_dir <- file.path(outdir, "KEGG")
  dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(kg_dir, recursive = TRUE, showWarnings = FALSE)
  pref <- paste0(tag, "_")
  mp <- map_to_entrez(genes)
  entrez <- unique(mp$entrez)
  fc_sym <- setNames(de_sub$log2FC, de_sub$gene)
  if (length(entrez) < 3) {
    writeLines(paste("mapped_entrez", length(entrez)), file.path(outdir, paste0(pref, "ORA_skipped.txt")))
    note_empty(file.path(go_dir, paste0(pref, "ORA_GO")), "too few mapped genes")
    note_empty(file.path(pw_dir, paste0(pref, "ORA_Pathway")), "too few mapped genes")
    note_empty(file.path(kg_dir, paste0(pref, "ORA_KEGG")), "too few mapped genes")
    return(invisible(NULL))
  }
  for (ont in c("BP", "MF", "CC")) {
    ego <- enrich_or_relax(
      function() clusterProfiler::enrichGO(
        gene = entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
        pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
      ),
      function() clusterProfiler::enrichGO(
        gene = entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
        pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
      ),
      paste("enrichGO", ont)
    )
    plot_ora_object(ego, file.path(go_dir, paste0(pref, "ORA_GO_", ont)),
                    paste(label, "| ORA GO", ont), fold_change = fc_sym)
  }
  ek <- enrich_or_relax(
    function() clusterProfiler::enrichKEGG(
      gene = entrez, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.2
    ),
    function() clusterProfiler::enrichKEGG(
      gene = entrez, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1
    ),
    "enrichKEGG"
  )
  if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) {
    ek <- tryCatch(clusterProfiler::setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
                   error = function(e) ek)
  }
  plot_ora_object(ek, file.path(kg_dir, paste0(pref, "ORA_KEGG")),
                  paste(label, "| ORA KEGG"), fold_change = fc_sym)
  if (has_pkg("ReactomePA")) {
    er <- enrich_or_relax(
      function() ReactomePA::enrichPathway(
        gene = entrez, organism = "human", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
      ),
      function() ReactomePA::enrichPathway(
        gene = entrez, organism = "human", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
      ),
      "enrichPathway"
    )
    plot_ora_object(er, file.path(pw_dir, paste0(pref, "ORA_Reactome_pathway")),
                    paste(label, "| ORA Reactome"), fold_change = fc_sym)
  } else {
    note_empty(file.path(pw_dir, paste0(pref, "ORA_Reactome_pathway")), "ReactomePA not installed")
  }
  writeLines(
    c("GO/Pathway/KEGG 是 ORA，不是 GSEA。",
      "GSEA 在同级 GSEA/ 目录，文件名以 GSEA_ 开头。"),
    file.path(outdir, paste0(pref, "00_ORA_is_not_GSEA.txt"))
  )
}

run_gsea_plots <- function(full_de, sub, outdir, tag, label) {
  gsea_dir <- file.path(outdir, "GSEA")
  dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
  pref <- paste0(tag, "_")
  stats <- ranked_entrez(full_de)
  if (length(stats) < 8) {
    note_empty(file.path(gsea_dir, paste0(pref, "GSEA")), "too few ranked genes")
    return(invisible(NULL))
  }
  kegg <- enrich_or_relax(
    function() clusterProfiler::gseKEGG(
      geneList = stats, organism = "hsa", minGSSize = 5, maxGSSize = 500,
      pvalueCutoff = 0.05, verbose = FALSE, eps = 0
    ),
    function() clusterProfiler::gseKEGG(
      geneList = stats, organism = "hsa", minGSSize = 3, maxGSSize = 500,
      pvalueCutoff = 1, verbose = FALSE, eps = 0
    ),
    paste("gseKEGG", tag)
  )
  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
    kegg <- tryCatch(clusterProfiler::setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
                     error = function(e) kegg)
  }
  plot_ora_object(kegg, file.path(gsea_dir, paste0(pref, "GSEA_KEGG")),
                  paste(label, "| GSEA KEGG"))
  if (has_pkg("msigdbr") && has_pkg("fgsea")) {
    hm <- tryCatch({
      m <- tryCatch(
        msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
        error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H")
      )
      gs <- split(as.character(m$entrez_gene), m$gs_name)
      fgsea::fgsea(pathways = gs, stats = stats, minSize = 5, maxSize = 500)
    }, error = function(e) {
      log_msg("fgsea hallmark failed: ", e$message)
      NULL
    })
    if (!is.null(hm) && nrow(hm) > 0) {
      hm <- hm[order(hm$pval), ]
      utils::write.csv(as.data.frame(hm), file.path(gsea_dir, paste0(pref, "GSEA_Hallmark.csv")),
                       row.names = FALSE)
    }
  }
}

emit_subset <- function(comp_name, de_full, sub, tag, title, outdir, fc_line) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(paste("comparison:", comp_name),
      paste("subset:", tag),
      paste("title:", title),
      paste("n_genes:", if (is.null(sub)) 0 else nrow(sub)),
      "FC = matched metastasis / matched primary (same patient).",
      "Selected genes are LOW in the primary relative to the paired metastasis."),
    file.path(outdir, paste0("00_", tag, "_THIS_FOLDER.txt"))
  )
  if (is.null(sub) || nrow(sub) == 0) {
    writeLines("no genes", file.path(outdir, paste0(tag, "_EMPTY.txt")))
    log_msg(comp_name, " ", tag, ": EMPTY")
    return(invisible(NULL))
  }
  utils::write.csv(sub, file.path(outdir, paste0(tag, "_DE_selected_genes.csv")), row.names = FALSE)
  tryCatch(writexl::write_xlsx(sub, file.path(outdir, paste0(tag, "_DE_selected_genes.xlsx"))),
           error = function(e) log_msg("xlsx write failed: ", e$message))
  log_msg(comp_name, " ", tag, ": n = ", nrow(sub))
  tryCatch(plot_de_bar(sub, paste(title, "| genes low in primary"), file.path(outdir, paste0(tag, "_DE_log2FC_barplot"))),
           error = function(e) log_msg("barplot failed: ", e$message))
  tryCatch(plot_volcano(de_full, sub$gene, title, file.path(outdir, paste0(tag, "_volcano")), fc_line = fc_line),
           error = function(e) log_msg("volcano failed: ", e$message))
  hm_mode <- attr(de_full, "heatmap_mode")
  if (identical(hm_mode, "group")) {
    tryCatch(group_heatmap(
      attr(de_full, "heat_mat"), attr(de_full, "heat_samples"),
      attr(de_full, "heat_groups"), sub$gene, title,
      file.path(outdir, paste0(tag, "_heatmap_high_vs_low"))
    ), error = function(e) {
      while (grDevices::dev.cur() > 1) grDevices::dev.off()
      log_msg("group heatmap failed: ", e$message)
    })
  } else {
    tryCatch(paired_heatmap(de_full, sub$gene, title, file.path(outdir, paste0(tag, "_heatmap_paired_1to1"))),
             error = function(e) {
               while (grDevices::dev.cur() > 1) grDevices::dev.off()
               log_msg("heatmap failed: ", e$message)
             })
  }
  tryCatch(run_ora_plots(sub$gene, sub, outdir, title, tag),
           error = function(e) log_msg("ORA failed: ", e$message))
  tryCatch(run_gsea_plots(de_full, sub, outdir, tag, title),
           error = function(e) log_msg("GSEA failed: ", e$message))
}

analyze_paired <- function(comp_name, pairs, folder) {
  base <- file.path(result_dir, folder)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  de <- paired_limma(pairs, comp_name)
  if (is.null(de)) {
    writeLines("not enough pairs", file.path(base, "SKIPPED_not_enough_pairs.txt"))
    return(NULL)
  }
  utils::write.csv(de, file.path(base, paste0(comp_name, "_full_genome_paired_DE.csv")), row.names = FALSE)
  tryCatch(writexl::write_xlsx(de, file.path(base, paste0(comp_name, "_full_genome_paired_DE.xlsx"))),
           error = function(e) log_msg(e$message))
  gsea_all <- file.path(base, "00_GSEA_all_genes_NOT_FC_or_topN")
  dir.create(gsea_all, recursive = TRUE, showWarnings = FALSE)
  tryCatch(run_gsea_plots(de, de, gsea_all, "ALLgenes", paste(comp_name, "all genes")),
           error = function(e) log_msg("all-gene GSEA failed: ", e$message))
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    sub <- select_up_in_met(de, fc_min = fc)
    outdir <- file.path(base, "FoldChange", nm)
    title <- paste0(comp_name, " | p<", p_cutoff, " & FC>", fc, " (low in primary)")
    emit_subset(comp_name, de, sub, nm, title, outdir, fc_line = fc)
  }
  de
}

# -----------------------------------------------------------------------------
# Q1 配对肺 / 骨
# -----------------------------------------------------------------------------
de_lung <- analyze_paired("paired_primary_vs_lung", pairs_lung, "01_lung_paired_met_vs_primary")
de_bone <- analyze_paired("paired_primary_vs_bone", pairs_bone, "02_bone_paired_met_vs_primary")

# -----------------------------------------------------------------------------
# Q2 器官特异：肺对过阈值但骨对不过，或反过来
# -----------------------------------------------------------------------------
organ_dir <- file.path(result_dir, "03_organ_specific")
dir.create(organ_dir, recursive = TRUE, showWarnings = FALSE)

write_setdiff <- function(a, b, name_a, name_b, out_stub, de_ref, folder_tag) {
  genes_a <- if (is.null(a) || nrow(a) == 0) character(0) else a$gene
  genes_b <- if (is.null(b) || nrow(b) == 0) character(0) else b$gene
  only_a <- setdiff(genes_a, genes_b)
  empty <- if (is.null(a)) data.frame() else a[0, ]
  sub <- if (is.null(a) || length(only_a) == 0) empty else a[a$gene %in% only_a, , drop = FALSE]
  title <- paste0(name_a, " specific vs ", name_b)
  emit_subset(title, de_ref, sub, folder_tag, title, out_stub, fc_line = 1)
  invisible(sub)
}

if (!is.null(de_lung) && !is.null(de_bone)) {
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    lung_sub <- select_up_in_met(de_lung, fc)
    bone_sub <- select_up_in_met(de_bone, fc)
    write_setdiff(
      lung_sub, bone_sub, "Lung", "Bone",
      file.path(organ_dir, nm, "lung_specific_not_bone"),
      de_lung, paste0(nm, "_lung_specific")
    )
    write_setdiff(
      bone_sub, lung_sub, "Bone", "Lung",
      file.path(organ_dir, nm, "bone_specific_not_lung"),
      de_bone, paste0(nm, "_bone_specific")
    )
    both <- intersect(lung_sub$gene, bone_sub$gene)
    overlap <- if (length(both) == 0) lung_sub[0, ] else lung_sub[lung_sub$gene %in% both, ]
    emit_subset(
      "lung_and_bone_shared", de_lung, overlap, paste0(nm, "_shared"),
      paste0("Shared lung & bone | p<", p_cutoff, " FC>", fc),
      file.path(organ_dir, nm, "shared_lung_and_bone"), fc
    )
    if (has_pkg("ggvenn")) {
      vdf <- list(Lung = lung_sub$gene, Bone = bone_sub$gene)
      p <- ggvenn::ggvenn(vdf, fill_color = c("#4C78A8", "#D62828")) +
        ggplot2::ggtitle(paste("Organ-specific DE", nm, "p<", p_cutoff))
      save_gg(p, file.path(organ_dir, nm, paste0(nm, "_venn_lung_vs_bone")))
    }
  }
} else {
  writeLines(
    "Need both lung and bone paired DE to call organ-specific genes.",
    file.path(organ_dir, "SKIPPED.txt")
  )
}

lung_only_pid <- setdiff(pairs_lung$patient, pairs_bone$patient)
bone_only_pid <- setdiff(pairs_bone$patient, pairs_lung$patient)
dual_pid <- intersect(pairs_lung$patient, pairs_bone$patient)
tropism <- data.frame(
  patient = c(lung_only_pid, bone_only_pid, dual_pid),
  tropism = c(rep("lung_only", length(lung_only_pid)),
              rep("bone_only", length(bone_only_pid)),
              rep("lung_and_bone", length(dual_pid))),
  stringsAsFactors = FALSE
)
utils::write.csv(tropism, file.path(organ_dir, "patients_lung_only_vs_bone_only.csv"), row.names = FALSE)
log_msg("Lung-only paired patients: ", paste(lung_only_pid, collapse = ", "))
log_msg("Bone-only paired patients: ", paste(bone_only_pid, collapse = ", "))
log_msg("Both lung and bone paired: ", paste(dual_pid, collapse = ", "))

run_tropism_de <- function(pid_a, pid_b, lab_a, lab_b, pairs_a, pairs_b, out_stub, title_low_in_a) {
  if (length(pid_a) < 2 || length(pid_b) < 2) return(invisible(NULL))
  prim_ids <- c(pairs_a$primary[match(pid_a, pairs_a$patient)],
                pairs_b$primary[match(pid_b, pairs_b$patient)])
  if (any(is.na(prim_ids))) {
    log_msg("Tropism DE skipped: missing primary IDs")
    return(invisible(NULL))
  }
  grp <- factor(c(rep(lab_a, length(pid_a)), rep(lab_b, length(pid_b))),
                levels = c(lab_b, lab_a))
  design <- stats::model.matrix(~ grp)
  fit <- limma::lmFit(logmat[, prim_ids, drop = FALSE], design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  coef <- grep(paste0("grp", lab_a), colnames(design), value = TRUE, fixed = TRUE)
  if (length(coef) == 0) coef <- colnames(design)[ncol(design)]
  tt <- limma::topTable(fit, coef = coef, number = Inf)
  trop_de <- data.frame(
    gene = rownames(tt), log2FC = tt$logFC, FC = 2^tt$logFC,
    pvalue = tt$P.Value, padj = tt$adj.P.Val, AveExpr = tt$AveExpr,
    stringsAsFactors = FALSE
  )
  trop_de <- trop_de[order(trop_de$pvalue), ]
  utils::write.csv(trop_de, paste0(out_stub, "_full_DE.csv"), row.names = FALSE)
  attr(trop_de, "heatmap_mode") <- "group"
  attr(trop_de, "heat_mat") <- logmat
  attr(trop_de, "heat_samples") <- prim_ids
  attr(trop_de, "heat_groups") <- stats::setNames(as.character(grp), prim_ids)
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    sub <- trop_de[!is.na(trop_de$pvalue) & trop_de$pvalue < p_cutoff & trop_de$FC < (1 / fc), ]
    if (nrow(sub) > 0) sub$FC_low_in_group_a <- 1 / sub$FC
    emit_subset(
      title_low_in_a, trop_de, sub, nm,
      paste0(title_low_in_a, " | p<", p_cutoff, " FC>", fc),
      file.path(out_stub, nm), fc
    )
  }
  invisible(trop_de)
}

if (length(lung_only_pid) >= 2 && length(bone_only_pid) >= 2) {
  run_tropism_de(
    lung_only_pid, bone_only_pid, "lung_only", "bone_only",
    pairs_lung, pairs_bone,
    file.path(organ_dir, "primary_lung_only_vs_bone_only"),
    "Primary lower in lung-only vs bone-only"
  )
} else {
  writeLines(
    c("No bone-only patients with a matched primary (all bone-paired patients also have lung).",
      paste("lung_only:", paste(lung_only_pid, collapse = ", ")),
      paste("bone_only:", paste(bone_only_pid, collapse = ", ")),
      paste("lung_and_bone:", paste(dual_pid, collapse = ", ")),
      "Falling back to lung-only vs dual lung+bone tropism if both groups have n>=2."),
    file.path(organ_dir, "NOTE_no_bone_only_patients.txt")
  )
  if (length(lung_only_pid) >= 2 && length(dual_pid) >= 2) {
    run_tropism_de(
      lung_only_pid, dual_pid, "lung_only", "lung_and_bone",
      pairs_lung, pairs_lung,
      file.path(organ_dir, "primary_lung_only_vs_dual_lung_bone"),
      "Primary lower in lung-only vs dual lung+bone"
    )
  }
}

# -----------------------------------------------------------------------------
# Q3 神经浸润签名（原发）
# GEO 无 PNI 病理；用轴突导向 / 施旺 / 神经营养 三套基因在原发打分
# -----------------------------------------------------------------------------
pni_curated <- list(
  axon_guidance = c(
    "SEMA3A", "SEMA3B", "SEMA3C", "SEMA3D", "SEMA3E", "SEMA3F", "SEMA3G",
    "SEMA4A", "SEMA4B", "SEMA4C", "SEMA4D", "SEMA4F", "SEMA4G",
    "SEMA5A", "SEMA5B", "SEMA6A", "SEMA6B", "SEMA6C", "SEMA6D", "SEMA7A",
    "PLXNA1", "PLXNA2", "PLXNA3", "PLXNA4", "PLXNB1", "PLXNB2", "PLXNB3",
    "PLXNC1", "PLXND1", "NRP1", "NRP2",
    "EPHA1", "EPHA2", "EPHA3", "EPHA4", "EPHA5", "EPHA6", "EPHA7", "EPHA8",
    "EPHB1", "EPHB2", "EPHB3", "EPHB4", "EPHB6",
    "EFNA1", "EFNA2", "EFNA3", "EFNA4", "EFNA5", "EFNB1", "EFNB2", "EFNB3",
    "NTN1", "NTN4", "DCC", "UNC5A", "UNC5B", "UNC5C", "UNC5D",
    "SLIT1", "SLIT2", "SLIT3", "ROBO1", "ROBO2", "ROBO3", "ROBO4",
    "NTNG1", "NTNG2", "RHOA", "RAC1", "CDC42", "PAK1", "GSK3B", "MAPK1"
  ),
  schwann = c(
    "SOX10", "S100B", "MPZ", "MBP", "PMP22", "MAG", "PRX", "EGR2", "POU3F1",
    "ERBB2", "ERBB3", "NRG1", "NCAM1", "L1CAM", "NGFR", "GFAP", "GAP43",
    "MAL", "GJB1", "PMP2", "DRP2", "PLP1", "MPZL1", "CDH19", "FOXD3",
    "SOX2", "OCT6", "POU3F2", "ITGA4", "ITGB8", "LAMA2", "LAMB2"
  ),
  neurotrophic = c(
    "NGF", "BDNF", "NTF3", "NTF4", "NTRK1", "NTRK2", "NTRK3", "NGFR",
    "GDNF", "NRTN", "ARTN", "PSPN", "GFRA1", "GFRA2", "GFRA3", "GFRA4",
    "RET", "CNTF", "CNTFR", "LIF", "LIFR", "OSM", "OSMR", "IGF1", "IGF1R",
    "VEGFA", "SORT1"
  )
)

expand_msig <- function(key_words) {
  if (!has_pkg("msigdbr")) return(character(0))
  m <- tryCatch({
    tryCatch(
      msigdbr::msigdbr(species = "Homo sapiens", collection = "C5"),
      error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "C5")
    )
  }, error = function(e) NULL)
  if (is.null(m) || !("gs_name" %in% names(m))) return(character(0))
  hit <- m[grepl(key_words, m$gs_name, ignore.case = TRUE), ]
  unique(as.character(hit$gene_symbol))
}

pni_sets <- list(
  axon_guidance = unique(c(pni_curated$axon_guidance, expand_msig("AXON_GUIDANCE"))),
  schwann = unique(c(pni_curated$schwann, expand_msig("SCHWANN"))),
  neurotrophic = unique(c(pni_curated$neurotrophic, expand_msig("NEUROTROPH")))
)

zscore_sig <- function(mat, genes) {
  genes <- intersect(unique(genes), rownames(mat))
  if (length(genes) < 3) {
    return(list(score = setNames(rep(NA_real_, ncol(mat)), colnames(mat)), genes = genes))
  }
  z <- t(scale(t(mat[genes, , drop = FALSE])))
  z[!is.finite(z)] <- NA
  list(score = colMeans(z, na.rm = TRUE), genes = genes)
}

prim_samples <- pheno_use$sample[pheno_use$kind == "primary"]
pni_dir <- file.path(result_dir, "04_PNI_primary")
dir.create(pni_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  c("GEO has no pathology PNI field. Scores are z-means of axon guidance / Schwann / neurotrophic genes on RAP primaries.",
    "GSE110590 is full-transcriptome RSEM (not a NanoString panel)."),
  file.path(pni_dir, "NOTE_signature_based_PNI.txt")
)

score_df <- data.frame(
  sample = prim_samples,
  patient = pheno_use$patient[match(prim_samples, pheno_use$sample)],
  has_paired_lung = pheno_use$patient[match(prim_samples, pheno_use$sample)] %in% pairs_lung$patient,
  has_paired_bone = pheno_use$patient[match(prim_samples, pheno_use$sample)] %in% pairs_bone$patient,
  stringsAsFactors = FALSE
)
used_genes <- list()
for (nm in names(pni_sets)) {
  sc <- zscore_sig(logmat[, prim_samples, drop = FALSE], pni_sets[[nm]])
  score_df[[paste0("score_", nm)]] <- sc$score[prim_samples]
  used_genes[[nm]] <- sc$genes
  writeLines(sc$genes, file.path(pni_dir, paste0("signature_genes_", nm, ".txt")))
  log_msg("PNI signature ", nm, ": ", length(sc$genes), " genes in matrix")
}
utils::write.csv(score_df, file.path(pni_dir, "primary_PNI_signature_scores.csv"), row.names = FALSE)

unpaired_limma <- function(expr, group, coef_name) {
  group <- factor(group)
  design <- stats::model.matrix(~ group)
  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  coef <- grep(coef_name, colnames(design), value = TRUE)
  if (length(coef) == 0) coef <- colnames(design)[ncol(design)]
  tt <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
  data.frame(
    gene = rownames(tt), log2FC = tt$logFC, FC = 2^tt$logFC,
    AveExpr = tt$AveExpr, t = tt$t, pvalue = tt$P.Value, padj = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
}

# 高神经浸润 vs 低：基因在高分组更低 => log2FC(high/low) < 0，FC_low/high > 阈值
pni_de_list <- list()
for (nm in names(pni_sets)) {
  base <- file.path(pni_dir, nm)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  genes_present <- intersect(pni_sets[[nm]], rownames(logmat))
  if (length(genes_present) < 3) {
    writeLines(paste("too few signature genes in matrix:", paste(genes_present, collapse = ", ")),
               file.path(base, "SKIPPED_too_few_signature_genes.txt"))
    log_msg("PNI ", nm, " skipped: only ", length(genes_present), " genes")
    next
  }
  sc <- score_df[[paste0("score_", nm)]]
  ok <- prim_samples[!is.na(sc)]
  sc2 <- sc[match(ok, prim_samples)]
  hi <- ok[sc2 >= stats::median(sc2, na.rm = TRUE)]
  lo <- ok[sc2 < stats::median(sc2, na.rm = TRUE)]
  log_msg("PNI ", nm, " high n=", length(hi), " low n=", length(lo))
  expr <- logmat[, c(lo, hi), drop = FALSE]
  grp <- factor(c(rep("low", length(lo)), rep("high", length(hi))), levels = c("low", "high"))
  de <- unpaired_limma(expr, grp, "high")
  de <- de[order(de$pvalue), ]
  # 高 PNI 原发里更低：FC_high/low < 1
  de$FC_low_over_high <- 1 / de$FC
  pni_de_list[[nm]] <- de
  utils::write.csv(de, file.path(base, paste0(nm, "_high_vs_low_full_DE.csv")), row.names = FALSE)
  grp_ids <- c(lo, hi)
  grp_lab <- setNames(c(rep("low", length(lo)), rep("high", length(hi))), grp_ids)
  for (fcnm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[fcnm]])
    sub <- de[!is.na(de$pvalue) & de$pvalue < p_cutoff & de$FC_low_over_high > fc, ]
    if (nrow(sub) > 0) {
      sub$log2FC <- log2(sub$FC_low_over_high)
      sub$FC <- sub$FC_low_over_high
    }
    full_for_volcano <- de
    full_for_volcano$log2FC <- log2(pmax(full_for_volcano$FC_low_over_high, 1e-8))
    full_for_volcano$FC <- full_for_volcano$FC_low_over_high
    attr(full_for_volcano, "heatmap_mode") <- "group"
    attr(full_for_volcano, "heat_mat") <- logmat
    attr(full_for_volcano, "heat_samples") <- grp_ids
    attr(full_for_volcano, "heat_groups") <- grp_lab
    emit_subset(
      paste0("PNI_", nm), full_for_volcano, sub, fcnm,
      paste0("PNI ", nm, " | genes lower in high-score primaries | p<", p_cutoff, " FC>", fc),
      file.path(base, "FoldChange", fcnm), fc
    )
  }
}

# -----------------------------------------------------------------------------
# Q4 三种神经浸润 vs 肺转移（分别做）
# -----------------------------------------------------------------------------
q4_dir <- file.path(result_dir, "05_PNI_vs_lung_met")
dir.create(q4_dir, recursive = TRUE, showWarnings = FALSE)

score_long <- tidyr::pivot_longer(
  score_df,
  cols = dplyr::starts_with("score_"),
  names_to = "signature",
  values_to = "score"
)
score_long$signature <- sub("^score_", "", score_long$signature)
score_long$lung_status <- ifelse(score_long$has_paired_lung,
                                 "primary_with_paired_lung",
                                 "primary_without_paired_lung")

p <- ggplot2::ggplot(score_long, ggplot2::aes(x = lung_status, y = score, fill = signature)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  ggplot2::geom_jitter(width = 0.12, size = 1.4, alpha = 0.8) +
  ggplot2::facet_wrap(~ signature, scales = "free_y") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)) +
  ggplot2::labs(
    title = "Primary PNI signature scores vs paired lung metastasis",
    x = NULL, y = "z-score (primary tumor)"
  )
save_gg(p, file.path(q4_dir, "PNI_scores_by_lung_pair_status"), width = 11, height = 5)

wilcox_rows <- list()
for (nm in names(pni_sets)) {
  sub <- score_df[!is.na(score_df[[paste0("score_", nm)]]), ]
  a <- sub[[paste0("score_", nm)]][sub$has_paired_lung]
  b <- sub[[paste0("score_", nm)]][!sub$has_paired_lung]
  wt <- tryCatch(stats::wilcox.test(a, b, exact = FALSE), error = function(e) NULL)
  wilcox_rows[[nm]] <- data.frame(
    signature = nm,
    n_primary_with_paired_lung = length(a),
    n_primary_without_paired_lung = length(b),
    median_with_lung = stats::median(a, na.rm = TRUE),
    median_without_lung = stats::median(b, na.rm = TRUE),
    wilcox_p = if (is.null(wt)) NA_real_ else unname(wt$p.value),
    stringsAsFactors = FALSE
  )
}
wilcox_tab <- dplyr::bind_rows(wilcox_rows)
utils::write.csv(wilcox_tab, file.path(q4_dir, "PNI_score_wilcox_paired_lung_vs_not.csv"), row.names = FALSE)

# 配对肺患者：原发 vs 肺转移 的三种签名分数（同一患者一一对应）
if (nrow(pairs_lung) >= 2) {
  vec <- pair_vectors(pairs_lung)
  paired_score <- list()
  for (nm in names(pni_sets)) {
    sc_p <- zscore_sig(vec$primary, pni_sets[[nm]])$score
    sc_m <- zscore_sig(vec$met, pni_sets[[nm]])$score
    df <- data.frame(
      patient = names(sc_p),
      signature = nm,
      score_primary = unname(sc_p),
      score_lung_met = unname(sc_m),
      delta_met_minus_primary = unname(sc_m - sc_p),
      stringsAsFactors = FALSE
    )
    paired_score[[nm]] <- df
    wt <- tryCatch(stats::wilcox.test(df$score_lung_met, df$score_primary, paired = TRUE, exact = FALSE),
                   error = function(e) NULL)
    log_msg("Paired PNI ", nm, " lung vs primary Wilcoxon p=",
            if (is.null(wt)) "NA" else signif(wt$p.value, 3))
    pdf <- tidyr::pivot_longer(df, cols = c("score_primary", "score_lung_met"),
                               names_to = "tissue", values_to = "score")
    pdf$tissue <- ifelse(pdf$tissue == "score_primary", "Primary", "Lung met")
    pp <- ggplot2::ggplot(pdf, ggplot2::aes(x = tissue, y = score, group = patient)) +
      ggplot2::geom_line(color = "grey60") +
      ggplot2::geom_point(ggplot2::aes(color = tissue), size = 2.5) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(
        title = paste0("1-to-1 PNI ", nm, " : primary vs paired lung"),
        x = NULL, y = "signature z-score"
      )
    save_gg(pp, file.path(q4_dir, nm, paste0(nm, "_paired_primary_vs_lung_score")))
  }
  utils::write.csv(dplyr::bind_rows(paired_score),
                   file.path(q4_dir, "paired_PNI_scores_primary_vs_lung.csv"),
                   row.names = FALSE)
}

# 基因交集：原发相对肺转移更低 ∩ 高神经浸润原发里更低
if (!is.null(de_lung)) {
  for (nm in names(pni_de_list)) {
    for (fcnm in names(fc_cutoffs)) {
      fc <- unname(fc_cutoffs[[fcnm]])
      lung_sub <- select_up_in_met(de_lung, fc)
      pni_de <- pni_de_list[[nm]]
      pni_sub <- pni_de[!is.na(pni_de$pvalue) & pni_de$pvalue < p_cutoff & pni_de$FC_low_over_high > fc, ]
      both <- intersect(lung_sub$gene, pni_sub$gene)
      ov <- if (length(both) == 0) lung_sub[0, ] else lung_sub[lung_sub$gene %in% both, ]
      emit_subset(
        paste0("PNI_", nm, "_AND_lung"),
        de_lung, ov, fcnm,
        paste0(nm, " PNI-low genes also low in primary vs paired lung | ", fcnm),
        file.path(q4_dir, nm, "overlap_with_lung_paired", fcnm), fc
      )
    }
  }
}

# -----------------------------------------------------------------------------
# 总表
# -----------------------------------------------------------------------------
summarize_n <- function(de, fc) {
  if (is.null(de) || nrow(de) == 0) return(0L)
  nrow(select_up_in_met(de, fc))
}
summary_tab <- data.frame(
  question = c(
    "Q1_lung_paired", "Q1_bone_paired",
    "Q2_lung_specific_FC1.25", "Q2_bone_specific_FC1.25",
    "Q3_axon_guidance", "Q3_schwann", "Q3_neurotrophic"
  ),
  n_pairs_or_primaries = c(
    nrow(pairs_lung), nrow(pairs_bone),
    nrow(pairs_lung), nrow(pairs_bone),
    sum(!is.na(score_df$score_axon_guidance)),
    sum(!is.na(score_df$score_schwann)),
    sum(!is.na(score_df$score_neurotrophic))
  ),
  n_genes_p0.05_FC_gt_1 = c(
    summarize_n(de_lung, 1), summarize_n(de_bone, 1),
    NA, NA, NA, NA, NA
  ),
  n_genes_p0.05_FC_gt_1.25 = c(
    summarize_n(de_lung, 1.25), summarize_n(de_bone, 1.25),
    NA, NA, NA, NA, NA
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(summary_tab, file.path(result_dir, "00_run_summary.csv"), row.names = FALSE)
writeLines(
  c(
    "GSE110590_Human_breast.R finished.",
    paste("Pairs lung:", nrow(pairs_lung), paste(pairs_lung$patient, collapse = ", ")),
    paste("Pairs bone:", nrow(pairs_bone), paste(pairs_bone$patient, collapse = ", ")),
    "Heatmaps keep patient order: Primary_i next to Met_i (cluster_cols = FALSE).",
    "No Top50-300. Thresholds: p < 0.05 and FC > 1 / 1.25 only.",
    "PNI is signature-based (no pathology PNI in GEO series matrix).",
    "A8 has mets but no primary in this matrix and is excluded from paired DE.",
    "No bone-only patients: all bone-paired patients also have a lung pair."
  ),
  file.path(result_dir, "00_README_results.txt")
)
log_msg("Done. Results in ", result_dir)
