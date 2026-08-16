#!/usr/bin/env Rscript
# =============================================================================
# TCGA BRCA Nature 2012（cBioPortal）——每个 GO 通路独立分析
#
# 1) 各 GO 通路得分 vs PAM50 分型、远处转移(M)、淋巴结转移(N) ：气泡图
# 2) 各 GO 通路得分 vs 全基因组表达：Spearman 负相关基因
#
# 默认数据目录：E:/R/TCGA_BRCA
# 在 RStudio 中：确认 DATA_DIR 后直接 Source 本文件
# 命令行：
#   Rscript scripts/analyze_go_subtype_metastasis.R
#   Rscript scripts/analyze_go_subtype_metastasis.R --data-dir "E:/R/TCGA_BRCA" --out-dir "E:/R/TCGA_BRCA/results"
#   Rscript scripts/analyze_go_subtype_metastasis.R --demo
#
# 不需要 AnnotationDbi / GSVA。可选 CRAN 包：ggplot2、data.table。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
invisible(Sys.setlocale("LC_COLLATE", "C"))

# -----------------------------------------------------------------------------
# 可改参数（RStudio 用户主要改这里）
# -----------------------------------------------------------------------------

DATA_DIR <- "E:/R/TCGA_BRCA"
OUT_DIR  <- file.path(DATA_DIR, "results")
CACHE_DIR <- file.path(DATA_DIR, "cache", "go_genes")
FDR_CUTOFF <- 0.05
MIN_GENES <- 2L
MIN_GROUP_N <- 5L
RUN_DEMO <- FALSE

# 用户原文 GO_2001224s 视为 GO:2001224。用小写 c()，不要写成 C()。
GO_IDS <- c(
  "GO:0007158",
  "GO:0007409",
  "GO:0007411",
  "GO:0008089",
  "GO:0019227",
  "GO:0019228",
  "GO:0021872",
  "GO:0021879",
  "GO:0023041",
  "GO:0031102",
  "GO:0031103",
  "GO:0033564",
  "GO:0036514",
  "GO:0036518",
  "GO:0048168",
  "GO:0048169",
  "GO:0060071",
  "GO:0070286",
  "GO:0071542",
  "GO:0097374",
  "GO:0097491",
  "GO:0097492",
  "GO:0098930",
  "GO:0140058",
  "GO:1902667",
  "GO:1902847",
  "GO:1904340",
  "GO:1904457",
  "GO:1904936",
  "GO:2001222",
  "GO:2001224"
)

GO_NAMES <- c(
  "GO:0007158" = "neuron cell-cell adhesion",
  "GO:0007409" = "axonogenesis",
  "GO:0007411" = "axon guidance",
  "GO:0008089" = "anterograde axonal transport",
  "GO:0019227" = "neuronal action potential propagation",
  "GO:0019228" = "neuronal action potential",
  "GO:0021872" = "forebrain generation of neurons",
  "GO:0021879" = "forebrain neuron differentiation",
  "GO:0023041" = "neuronal signal transduction",
  "GO:0031102" = "neuron projection regeneration",
  "GO:0031103" = "axon regeneration",
  "GO:0033564" = "anterior/posterior axon guidance",
  "GO:0036514" = "dopaminergic neuron axon guidance",
  "GO:0036518" = "chemorepulsion of dopaminergic neuron axon",
  "GO:0048168" = "regulation of neuronal synaptic plasticity",
  "GO:0048169" = "regulation of long-term neuronal synaptic plasticity",
  "GO:0060071" = "Wnt PCP pathway",
  "GO:0070286" = "axonemal dynein complex assembly",
  "GO:0071542" = "dopaminergic neuron differentiation",
  "GO:0097374" = "sensory neuron axon guidance",
  "GO:0097491" = "sympathetic neuron projection guidance",
  "GO:0097492" = "sympathetic neuron axon guidance",
  "GO:0098930" = "axonal transport",
  "GO:0140058" = "neuron projection arborization",
  "GO:1902667" = "regulation of axon guidance",
  "GO:1902847" = "regulation of neuronal signal transduction",
  "GO:1904340" = "positive regulation of dopaminergic neuron differentiation",
  "GO:1904457" = "positive regulation of neuronal action potential",
  "GO:1904936" = "interneuron migration",
  "GO:2001222" = "regulation of neuron migration",
  "GO:2001224" = "positive regulation of neuron migration"
)

# -----------------------------------------------------------------------------
# 工具
# -----------------------------------------------------------------------------

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

try_install <- function(pkgs) {
  for (p in pkgs) {
    if (!has_pkg(p)) {
      log_msg("安装 CRAN 包: ", p)
      tryCatch(
        install.packages(p, repos = "https://cloud.r-project.org"),
        error = function(e) log_msg("安装失败 ", p, ": ", conditionMessage(e))
      )
    }
  }
}

parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  defaults <- list(
    data_dir = DATA_DIR,
    out_dir = OUT_DIR,
    cache_dir = CACHE_DIR,
    fdr = FDR_CUTOFF,
    min_genes = MIN_GENES,
    demo = RUN_DEMO
  )
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (key == "--demo") {
      defaults$demo <- TRUE
      i <- i + 1L
    } else if (i < length(argv) && startsWith(key, "--")) {
      val <- argv[[i + 1L]]
      name <- gsub("-", "_", sub("^--", "", key))
      if (!name %in% names(defaults)) stop("未知参数: ", key)
      if (is.numeric(defaults[[name]])) val <- as.numeric(val)
      defaults[[name]] <- val
      i <- i + 2L
    } else {
      stop("无法解析参数: ", key)
    }
  }
  defaults
}

norm_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("\\.", "-", x)
  x[x %in% c("", "NA", "NAN", "NULL", "NONE", "[NOT AVAILABLE]", "[NOT APPLICABLE]")] <- NA_character_
  x
}

patient_id <- function(x) substr(norm_id(x), 1, 12)

go_label <- function(gid) {
  nm <- unname(GO_NAMES[gid])
  if (is.na(nm) || !nzchar(nm)) gid else paste0(gid, "\n", nm)
}

go_file_stem <- function(gid) gsub(":", "_", gid)

na_blank <- function(x) {
  x <- as.character(x)
  x[trimws(x) == "" | toupper(trimws(x)) %in% c("NA", "NAN", "UNKNOWN", "N/A")] <- NA_character_
  x
}

# -----------------------------------------------------------------------------
# 寻找并读取 cBioPortal 文件
# -----------------------------------------------------------------------------

list_data_files <- function(data_dir) {
  if (!dir.exists(data_dir)) return(character())
  list.files(data_dir, recursive = TRUE, full.names = TRUE)
}

pick_file <- function(files, patterns, exclude = NULL) {
  base <- basename(files)
  keep <- files
  if (!is.null(exclude)) {
    keep <- keep[!grepl(exclude, basename(keep), ignore.case = TRUE)]
  }
  keep <- keep[!grepl("^meta_", basename(keep), ignore.case = TRUE)]
  for (p in patterns) {
    hit <- keep[grepl(p, basename(keep), ignore.case = TRUE)]
    if (length(hit)) return(hit[[1]])
  }
  NA_character_
}

read_cbioportal_table <- function(path) {
  log_msg("读取: ", path)
  if (has_pkg("data.table")) {
    df <- data.table::fread(
      path, sep = "\t", header = TRUE, data.table = FALSE,
      check.names = FALSE, quote = "", comment.char = "#"
    )
  } else {
    df <- utils::read.delim(
      path, sep = "\t", header = TRUE, check.names = FALSE,
      stringsAsFactors = FALSE, comment.char = "#", quote = ""
    )
  }
  names(df) <- gsub("^\\s+|\\s+$", "", names(df))
  df
}

read_expr_matrix <- function(path) {
  df <- read_cbioportal_table(path)
  nms <- names(df)
  drop <- which(toupper(nms) %in% c(
    "HUGO_SYMBOL", "ENTREZ_GENE_ID", "ENTREZGENEID", "GENE_SYMBOL",
    "GENE", "ID", "LOCUS_ID"
  ))
  gene_col <- if ("Hugo_Symbol" %in% nms) {
    "Hugo_Symbol"
  } else if (length(drop)) {
    nms[[drop[[1]]]]
  } else {
    nms[[1]]
  }
  sample_cols <- setdiff(seq_along(nms), match(gene_col, nms))
  if (length(drop)) sample_cols <- setdiff(sample_cols, drop)
  genes <- as.character(df[[gene_col]])
  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "double"
  colnames(mat) <- norm_id(colnames(mat))
  keep <- !is.na(genes) & nzchar(genes) & genes != "NA"
  mat <- mat[keep, , drop = FALSE]
  genes <- genes[keep]
  if (any(duplicated(genes))) {
    log_msg("合并重复基因符号: ", sum(duplicated(genes)))
    split_idx <- split(seq_along(genes), genes)
    mat <- do.call(rbind, lapply(split_idx, function(ii) {
      if (length(ii) == 1L) mat[ii, , drop = FALSE] else {
        matrix(colMeans(mat[ii, , drop = FALSE], na.rm = TRUE), nrow = 1,
               dimnames = list(NULL, colnames(mat)))
      }
    }))
    rownames(mat) <- names(split_idx)
  } else {
    rownames(mat) <- genes
  }
  finite_n <- rowSums(is.finite(mat))
  mat <- mat[finite_n >= 5L, , drop = FALSE]
  mx <- suppressWarnings(max(mat, na.rm = TRUE))
  if (is.finite(mx) && mx > 100) {
    log_msg("表达值偏大，做 log2(x - min + 1) 变换")
    minv <- min(mat, na.rm = TRUE)
    mat <- log2(mat - minv + 1)
  }
  mat
}

read_case_list <- function(path) {
  if (is.na(path) || !file.exists(path)) return(character())
  lines <- readLines(path, warn = FALSE)
  hit <- grep("^case_list_ids:", lines, ignore.case = TRUE)
  if (!length(hit)) return(character())
  ids <- sub("^case_list_ids:\\s*", "", lines[[hit[[1]]]], ignore.case = TRUE)
  unique(norm_id(unlist(strsplit(ids, "[,[:space:]]+"))))
}

subtype_from_case_lists <- function(data_dir) {
  mapping <- c(
    basal = "Basal-like",
    her2 = "HER2-enriched",
    luma = "Luminal A",
    lumb = "Luminal B",
    claudin = "Claudin-low",
    luminal = NA_character_
  )
  files <- list_data_files(data_dir)
  out <- list()
  for (key in names(mapping)) {
    if (is.na(mapping[[key]])) next
    path <- pick_file(files, paste0("^cases_", key, "(\\.txt)?$"))
    ids <- read_case_list(path)
    if (length(ids)) out[[mapping[[key]]]] <- ids
  }
  out
}

normalize_pam50 <- function(x) {
  x <- na_blank(x)
  y <- tolower(gsub("[^[:alnum:]]+", "", x, perl = TRUE))
  out <- rep(NA_character_, length(x))
  out[y %in% c("luminala", "luma", "la")] <- "Luminal A"
  out[y %in% c("luminalb", "lumb", "lb")] <- "Luminal B"
  out[y %in% c("basallike", "basal")] <- "Basal-like"
  out[y %in% c("her2enriched", "her2", "her2plus")] <- "HER2-enriched"
  out[y %in% c("normallike", "normal")] <- "Normal-like"
  out[y %in% c("claudinlow", "claudin")] <- "Claudin-low"
  leftover <- is.na(out) & !is.na(x)
  out[leftover] <- x[leftover]
  out
}

normalize_tnm <- function(x, kind = c("N", "M", "T")) {
  kind <- match.arg(kind)
  x <- toupper(na_blank(x))
  x <- gsub("\\s+", "", x)
  if (kind == "N") {
    out <- ifelse(grepl("N3", x), "N3",
           ifelse(grepl("N2", x), "N2",
           ifelse(grepl("N1", x), "N1",
           ifelse(grepl("N0", x), "N0", NA_character_))))
    pos <- grepl("POSITIVE|^POS$|^P$", x)
    neg <- grepl("NEGATIVE|^NEG$|^N$", x)
    out[is.na(out) & pos] <- "N+"
    out[is.na(out) & neg] <- "N0"
    out
  } else if (kind == "M") {
    out <- ifelse(grepl("M1", x) | grepl("POSITIVE|^YES|^Y$|^POS$", x), "M1",
           ifelse(grepl("M0", x) | grepl("NEGATIVE|^NO$|^N$|^NEG$", x), "M0",
                  NA_character_))
    out
  } else {
    out <- ifelse(grepl("T4", x), "T4",
           ifelse(grepl("T3", x), "T3",
           ifelse(grepl("T2", x), "T2",
           ifelse(grepl("T1", x), "T1", NA_character_))))
    out
  }
}

pick_col <- function(df, candidates) {
  nms <- toupper(gsub("[^[:alnum:]]+", "_", names(df), perl = TRUE))
  cand <- toupper(gsub("[^[:alnum:]]+", "_", candidates, perl = TRUE))
  for (c in cand) {
    i <- match(c, nms)
    if (!is.na(i)) return(names(df)[i])
  }
  for (c in cand) {
    hit <- which(grepl(c, nms))
    if (length(hit)) return(names(df)[hit[[1]]])
  }
  NA_character_
}

read_subtypes_file <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  df <- tryCatch(read_cbioportal_table(path), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  id_col <- pick_col(df, c("SAMPLE_ID", "PATIENT_ID", "SAMPLE", "ID"))
  sub_col <- pick_col(df, c("PAM50_SUBTYPE", "SUBTYPE", "PAM50", "CLAUDIN_SUBTYPE"))
  if (!is.na(id_col) && !is.na(sub_col)) {
    return(data.frame(
      SAMPLE_ID = norm_id(df[[id_col]]),
      SUBTYPE = normalize_pam50(df[[sub_col]]),
      stringsAsFactors = FALSE
    ))
  }
  NULL
}

build_clinical <- function(data_dir, sample_ids) {
  files <- list_data_files(data_dir)
  patient_path <- pick_file(files, c("^data_clinical_patient"))
  sample_path  <- pick_file(files, c("^data_clinical_sample"))
  subtype_path <- pick_file(files, c("^data_subtypes"))

  clin <- data.frame(SAMPLE_ID = sample_ids, stringsAsFactors = FALSE)
  clin$PATIENT_ID <- patient_id(clin$SAMPLE_ID)

  if (!is.na(sample_path)) {
    samp <- read_cbioportal_table(sample_path)
    sid <- pick_col(samp, c("SAMPLE_ID", "SAMPLE"))
    if (is.na(sid)) sid <- names(samp)[[1]]
    samp$._sid <- norm_id(samp[[sid]])
    samp$._pid <- patient_id(samp$._sid)
    idx <- match(clin$SAMPLE_ID, samp$._sid)
    miss <- is.na(idx)
    if (any(miss)) idx[miss] <- match(clin$PATIENT_ID[miss], samp$._pid)
    extra <- setdiff(names(samp), c("._sid", "._pid", sid))
    for (nm in extra) clin[[nm]] <- samp[[nm]][idx]
  }

  if (!is.na(patient_path)) {
    pat <- read_cbioportal_table(patient_path)
    pid <- pick_col(pat, c("PATIENT_ID", "PATIENT"))
    if (is.na(pid)) pid <- names(pat)[[1]]
    pat$._pid <- patient_id(pat[[pid]])
    idx <- match(clin$PATIENT_ID, pat$._pid)
    extra <- setdiff(names(pat), c("._pid", pid))
    for (nm in extra) {
      if (!nm %in% names(clin)) clin[[nm]] <- pat[[nm]][idx]
    }
  }

  sub_tbl <- read_subtypes_file(subtype_path)
  pam_col <- pick_col(clin, c("PAM50_SUBTYPE", "SUBTYPE", "PAM50"))
  clin$SUBTYPE <- NA_character_
  if (!is.na(pam_col)) clin$SUBTYPE <- normalize_pam50(clin[[pam_col]])
  if (!is.null(sub_tbl)) {
    idx <- match(clin$SAMPLE_ID, sub_tbl$SAMPLE_ID)
    miss <- is.na(idx)
    if (any(miss)) idx[miss] <- match(clin$PATIENT_ID[miss], patient_id(sub_tbl$SAMPLE_ID))
    fill <- is.na(clin$SUBTYPE) & !is.na(idx)
    clin$SUBTYPE[fill] <- sub_tbl$SUBTYPE[idx[fill]]
  }
  case_map <- subtype_from_case_lists(data_dir)
  for (lab in names(case_map)) {
    ids <- case_map[[lab]]
    hit <- clin$SAMPLE_ID %in% ids | clin$PATIENT_ID %in% patient_id(ids)
    clin$SUBTYPE[hit & is.na(clin$SUBTYPE)] <- lab
  }

  m_col <- pick_col(clin, c("METASTASIS", "METASTASIS_CODED", "M_STAGE", "AJCC_METASTASIS_PATHOLOGIC_PM"))
  n_col <- pick_col(clin, c("NODES", "NODE_CODED", "N_STAGE", "AJCC_NODES_PATHOLOGIC_PN"))
  t_col <- pick_col(clin, c("TUMOR_STAGE", "TUMOR_T1_CODED", "T_STAGE"))
  clin$M_STAGE <- if (!is.na(m_col)) normalize_tnm(clin[[m_col]], "M") else NA_character_
  clin$N_STAGE <- if (!is.na(n_col)) normalize_tnm(clin[[n_col]], "N") else NA_character_
  clin$T_STAGE <- if (!is.na(t_col)) normalize_tnm(clin[[t_col]], "T") else NA_character_
  clin$N_PLUS <- ifelse(is.na(clin$N_STAGE), NA_character_,
                        ifelse(clin$N_STAGE == "N0", "N0", "N+"))
  clin$MET_GROUP <- NA_character_
  none <- !is.na(clin$M_STAGE) & !is.na(clin$N_STAGE) & clin$M_STAGE == "M0" & clin$N_STAGE == "N0"
  ln   <- !is.na(clin$M_STAGE) & !is.na(clin$N_STAGE) & clin$M_STAGE == "M0" & clin$N_STAGE != "N0"
  dist <- !is.na(clin$M_STAGE) & clin$M_STAGE == "M1"
  clin$MET_GROUP[none] <- "No met (M0 N0)"
  clin$MET_GROUP[ln]   <- "Lymph node only (M0 N+)"
  clin$MET_GROUP[dist] <- "Distant (M1)"
  clin
}

# -----------------------------------------------------------------------------
# QuickGO 基因集（每个 GO 单独下载，含 descendants）
# -----------------------------------------------------------------------------

download_url <- function(url, dest, extra_headers = NULL) {
  hdr <- c(`User-Agent` = "LAI-tcga-brca-go")
  if (length(extra_headers)) hdr <- c(hdr, extra_headers)
  ok <- tryCatch({
    utils::download.file(
      url, destfile = dest, quiet = TRUE, mode = "wb",
      method = "libcurl", headers = hdr
    )
    file.exists(dest) && file.info(dest)$size > 0
  }, error = function(e) FALSE)
  if (ok) return(TRUE)
  if (nzchar(Sys.which("curl"))) {
    args <- c("-sS", "-L", "--fail", "-o", dest, "-A", "LAI-tcga-brca-go")
    if ("Accept" %in% names(hdr)) args <- c(args, "-H", paste0("Accept: ", hdr[["Accept"]]))
    args <- c(args, url)
    st <- suppressWarnings(system2("curl", args, stdout = FALSE, stderr = FALSE))
    return(isTRUE(st == 0L) && file.exists(dest) && file.info(dest)$size > 0)
  }
  FALSE
}

fetch_one_go_from_quickgo <- function(go_id) {
  url <- paste0(
    "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?",
    "goId=", utils::URLencode(go_id, reserved = TRUE),
    "&taxonId=9606&goUsage=descendants"
  )
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  if (!download_url(url, tmp, extra_headers = c(Accept = "text/tsv"))) return(NULL)
  first <- tryCatch(readLines(tmp, n = 1L, warn = FALSE), error = function(e) "")
  if (!length(first) || grepl("^\\s*\\{", first)) return(NULL)
  df <- utils::read.delim(tmp, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(df)) return(character())
  nms <- toupper(gsub("[^[:alnum:]]+", "_", names(df), perl = TRUE))
  names(df) <- nms
  sym_col <- intersect(c("SYMBOL", "GENE_PRODUCT_SYMBOL"), names(df))
  if (!length(sym_col)) return(character())
  symbol <- unique(as.character(df[[sym_col[[1]]]]))
  symbol[!is.na(symbol) & nzchar(symbol) & symbol != "-"]
}

get_genes_for_one_go <- function(go_id, cache_dir) {
  cache_file <- file.path(cache_dir, paste0(go_file_stem(go_id), "_genes.tsv"))
  if (file.exists(cache_file) && file.info(cache_file)$size > 0) {
    log_msg(go_id, " 使用缓存: ", cache_file)
    df <- utils::read.delim(cache_file, stringsAsFactors = FALSE)
    return(unique(na.omit(as.character(df$symbol))))
  }
  log_msg(go_id, " 从 QuickGO 下载基因集（含 descendants）")
  symbols <- tryCatch(fetch_one_go_from_quickgo(go_id), error = function(e) {
    log_msg(go_id, " QuickGO 失败: ", conditionMessage(e))
    NULL
  })
  if (is.null(symbols)) {
    stop("无法获取 ", go_id, " 的基因集。请检查网络，或手动写入: ", cache_file)
  }
  ensure_dir(cache_dir)
  utils::write.table(
    data.frame(symbol = symbols, stringsAsFactors = FALSE),
    cache_file, sep = "\t", quote = FALSE, row.names = FALSE
  )
  symbols
}

# -----------------------------------------------------------------------------
# 通路得分：通路内基因跨样本 z-score 后取均值
# -----------------------------------------------------------------------------

mean_z_score <- function(expr, genes) {
  g <- intersect(genes, rownames(expr))
  mat <- expr[g, , drop = FALSE]
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  stats::setNames(colMeans(z, na.rm = TRUE), colnames(expr))
}

# -----------------------------------------------------------------------------
# 关联检验
# -----------------------------------------------------------------------------

rank_biserial <- function(a, b) {
  # 正值：a 的得分更高
  if (length(a) < 2L || length(b) < 2L) return(NA_real_)
  wt <- suppressWarnings(stats::wilcox.test(a, b, exact = FALSE, correct = FALSE))
  n1 <- length(a); n2 <- length(b)
  u <- as.numeric(wt$statistic)
  (2 * u) / (n1 * n2) - 1
}

test_group_vs_rest <- function(score, group, min_n = MIN_GROUP_N) {
  group <- as.character(group)
  ok <- !is.na(score) & !is.na(group) & nzchar(group)
  score <- score[ok]; group <- group[ok]
  levs <- names(sort(table(group), decreasing = TRUE))
  rows <- list()
  if (length(levs) >= 3L && min(table(group)) >= 2L) {
    kt <- suppressWarnings(stats::kruskal.test(score ~ factor(group)))
    rows[[length(rows) + 1L]] <- data.frame(
      comparison = "overall Kruskal-Wallis",
      group = "all",
      n_group = length(score),
      n_ref = NA_integer_,
      mean_group = NA_real_,
      mean_ref = NA_real_,
      delta = NA_real_,
      effect = NA_real_,
      p_value = kt$p.value,
      stringsAsFactors = FALSE
    )
  }
  for (lv in levs) {
    in_g <- group == lv
    if (sum(in_g) < min_n || sum(!in_g) < min_n) next
    a <- score[in_g]; b <- score[!in_g]
    wt <- suppressWarnings(stats::wilcox.test(a, b, exact = FALSE))
    rows[[length(rows) + 1L]] <- data.frame(
      comparison = paste0(lv, " vs rest"),
      group = lv,
      n_group = length(a),
      n_ref = length(b),
      mean_group = mean(a, na.rm = TRUE),
      mean_ref = mean(b, na.rm = TRUE),
      delta = mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE),
      effect = rank_biserial(a, b),
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

test_binary <- function(score, label, pos, neg, min_n = MIN_GROUP_N) {
  a <- score[label == pos]
  b <- score[label == neg]
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < min_n || length(b) < min_n) return(NULL)
  wt <- suppressWarnings(stats::wilcox.test(a, b, exact = FALSE))
  data.frame(
    comparison = paste0(pos, " vs ", neg),
    group = pos,
    n_group = length(a),
    n_ref = length(b),
    mean_group = mean(a),
    mean_ref = mean(b),
    delta = mean(a) - mean(b),
    effect = rank_biserial(a, b),
    p_value = wt$p.value,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# 气泡图
# -----------------------------------------------------------------------------

plot_bubble <- function(df, title, outfile, y_order = NULL) {
  df <- df[!is.na(df$p_value) & df$comparison != "overall Kruskal-Wallis", , drop = FALSE]
  if (!nrow(df)) {
    log_msg("无可画气泡图的行: ", title)
    return(invisible(NULL))
  }
  df$neglog10p <- pmin(-log10(pmax(df$p_value, 1e-20)), 12)
  df$go_axis <- paste0(df$go_id, "  ", df$go_name)
  if (!is.null(y_order)) {
    df$group <- factor(df$group, levels = intersect(y_order, unique(df$group)))
    df <- df[!is.na(df$group), , drop = FALSE]
  }
  df$sig <- ifelse(!is.na(df$fdr) & df$fdr < FDR_CUTOFF, "FDR<0.05", "n.s.")
  ensure_dir(dirname(outfile))
  pdf_file <- outfile
  png_file <- sub("\\.pdf$", ".png", outfile)

  draw_base <- function() {
    gos <- unique(df$go_axis)
    grps <- if (is.factor(df$group)) levels(df$group) else unique(as.character(df$group))
    xs <- match(df$group, grps)
    ys <- match(df$go_axis, gos)
    col <- ifelse(df$delta >= 0, "#B2182B", "#2166AC")
    cex <- 0.4 + 2.2 * (df$neglog10p / max(df$neglog10p, na.rm = TRUE))
    op <- par(mar = c(8, 18, 3, 2), xpd = NA)
    on.exit(par(op), add = TRUE)
    plot(xs, ys, pch = 21, bg = col, col = ifelse(df$sig == "FDR<0.05", "black", "grey60"),
         cex = cex, xlim = c(0.5, length(grps) + 0.5),
         ylim = c(0.5, length(gos) + 0.5), xaxt = "n", yaxt = "n",
         xlab = "", ylab = "", main = title)
    axis(1, at = seq_along(grps), labels = grps, las = 2, cex.axis = 0.8)
    axis(2, at = seq_along(gos), labels = gos, las = 2, cex.axis = 0.55)
  }

  if (has_pkg("ggplot2")) {
    suppressPackageStartupMessages(library(ggplot2, quietly = TRUE))
    p <- ggplot(df, aes(x = .data$group, y = .data$go_axis, size = .data$neglog10p, color = .data$delta)) +
      geom_point(alpha = 0.9) +
      scale_color_gradient2(
        low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
        name = "delta\n(group - rest)"
      ) +
      scale_size_continuous(name = "-log10(p)", range = c(1.5, 9)) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 12, face = "bold")
      ) +
      labs(x = NULL, y = NULL, title = title,
           subtitle = "点大小 = -log10(p)；颜色 = 该组通路得分均值 − 其余样本")
    ggsave(pdf_file, p, width = 10, height = max(7, 0.22 * length(unique(df$go_id)) + 3))
    ggsave(png_file, p, width = 10, height = max(7, 0.22 * length(unique(df$go_id)) + 3), dpi = 150)
  } else {
    grDevices::pdf(pdf_file, width = 11, height = max(8, 0.22 * length(unique(df$go_id)) + 3))
    draw_base()
    grDevices::dev.off()
    grDevices::png(png_file, width = 1400, height = max(900, 28 * length(unique(df$go_id))), res = 120)
    draw_base()
    grDevices::dev.off()
  }
  log_msg("写出气泡图: ", pdf_file)
}

# -----------------------------------------------------------------------------
# 负相关基因
# -----------------------------------------------------------------------------

spearman_negative <- function(score, expr, fdr = FDR_CUTOFF) {
  common <- intersect(names(score), colnames(expr))
  sc <- score[common]
  mat <- expr[, common, drop = FALSE]
  keep <- is.finite(sc)
  sc <- sc[keep]
  mat <- mat[, keep, drop = FALSE]
  n <- ncol(mat)
  if (n < 10L) return(NULL)
  rho <- suppressWarnings(stats::cor(t(mat), sc, method = "spearman", use = "pairwise.complete.obs"))[, 1]
  nn <- rowSums(is.finite(mat))
  tstat <- rho * sqrt((nn - 2) / pmax(1e-12, 1 - rho^2))
  pval <- 2 * stats::pt(-abs(tstat), df = pmax(nn - 2, 1))
  pval[!is.finite(pval)] <- NA_real_
  fdr_v <- stats::p.adjust(pval, method = "BH")
  out <- data.frame(
    gene = names(rho),
    rho = as.numeric(rho),
    n = as.integer(nn),
    p_value = as.numeric(pval),
    fdr = as.numeric(fdr_v),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$rho) & out$rho < 0 & is.finite(out$fdr) & out$fdr < fdr, , drop = FALSE]
  out[order(out$rho, out$fdr), , drop = FALSE]
}

# -----------------------------------------------------------------------------
# Demo：不依赖本地 TCGA 文件
# -----------------------------------------------------------------------------

write_demo_data <- function(root) {
  ensure_dir(root)
  set.seed(1)
  n <- 80L
  samples <- sprintf("TCGA-A2-%04d-01", seq_len(n))
  patients <- substr(samples, 1, 12)
  subtype <- rep(c("Luminal A", "Luminal B", "Basal-like", "HER2-enriched"), length.out = n)
  mstage <- c(rep("M0", 70), rep("M1", 10))
  nstage <- c(rep("N0", 30), rep("N1", 25), rep("N2", 15), rep("N3", 10))
  samp <- data.frame(
    SAMPLE_ID = samples,
    PATIENT_ID = patients,
    PAM50_SUBTYPE = subtype,
    METASTASIS_CODED = ifelse(mstage == "M1", "Positive", "Negative"),
    NODES = nstage,
    NODE_CODED = ifelse(nstage == "N0", "Negative", "Positive"),
    TUMOR_STAGE = rep(c("T1", "T2", "T3"), length.out = n),
    stringsAsFactors = FALSE
  )
  pat <- data.frame(
    PATIENT_ID = patients,
    SEX = "Female",
    METASTASIS = mstage,
    OS_MONTHS = round(runif(n, 5, 120), 1),
    OS_STATUS = sample(c("0:LIVING", "1:DECEASED"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  genes <- c(paste0("GOG", sprintf("%03d", 1:15)), paste0("NEG", sprintf("%03d", 1:8)),
             paste0("BG", sprintf("%03d", 1:40)))
  latent <- rnorm(n)
  latent[subtype == "Basal-like"] <- latent[subtype == "Basal-like"] + 1.2
  latent[mstage == "M1"] <- latent[mstage == "M1"] + 0.8
  mat <- matrix(rnorm(length(genes) * n), nrow = length(genes), dimnames = list(genes, samples))
  mat[1:15, ] <- mat[1:15, ] + matrix(rep(latent, each = 15), nrow = 15)
  mat[16:23, ] <- mat[16:23, ] - matrix(rep(latent, each = 8), nrow = 8)
  expr <- data.frame(Hugo_Symbol = rownames(mat), Entrez_Gene_Id = seq_len(nrow(mat)), mat,
                     check.names = FALSE, stringsAsFactors = FALSE)
  utils::write.table(samp, file.path(root, "data_clinical_sample.txt"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(pat, file.path(root, "data_clinical_patient.txt"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(expr, file.path(root, "data_mrna_agilent_microarray.txt"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  cl_dir <- file.path(root, "case_lists")
  ensure_dir(cl_dir)
  write_cl <- function(name, ids) {
    writeLines(c(
      "cancer_study_identifier: brca_tcga_pub",
      paste0("case_list_ids: ", paste(ids, collapse = " "))
    ), file.path(cl_dir, name))
  }
  write_cl("cases_luma.txt", samples[subtype == "Luminal A"])
  write_cl("cases_lumb.txt", samples[subtype == "Luminal B"])
  write_cl("cases_basal.txt", samples[subtype == "Basal-like"])
  write_cl("cases_her2.txt", samples[subtype == "HER2-enriched"])
  cache <- file.path(root, "cache", "go_genes")
  ensure_dir(cache)
  demo_gos <- c("GO:0007411", "GO:0007158", "GO:2001224")
  for (gid in demo_gos) {
    utils::write.table(
      data.frame(symbol = paste0("GOG", sprintf("%03d", 1:15))),
      file.path(cache, paste0(go_file_stem(gid), "_genes.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  list(root = root, go_ids = demo_gos)
}

# -----------------------------------------------------------------------------
# 主流程
# -----------------------------------------------------------------------------

run_analysis <- function(data_dir, out_dir, cache_dir, go_ids, fdr, min_genes, demo = FALSE) {
  if (!isTRUE(demo)) try_install(c("ggplot2"))
  ensure_dir(out_dir)
  ensure_dir(file.path(out_dir, "neg_genes"))
  ensure_dir(cache_dir)

  files <- list_data_files(data_dir)
  expr_path <- pick_file(
    files,
    c("^data_mrna_agilent_microarray(\\.txt)?$"),
    exclude = "zscore"
  )
  if (is.na(expr_path)) {
    expr_path <- pick_file(files, c("data_mrna_agilent_microarray"))
  }
  if (is.na(expr_path)) {
    stop(
      "在 ", data_dir, " 找不到表达文件 data_mrna_agilent_microarray。\n",
      "请确认 cBioPortal 数据包已解压到该目录。"
    )
  }

  expr <- read_expr_matrix(expr_path)
  log_msg("表达矩阵: ", nrow(expr), " 基因 x ", ncol(expr), " 样本")
  clin <- build_clinical(data_dir, colnames(expr))
  rownames(clin) <- clin$SAMPLE_ID

  log_msg("PAM50 分型计数:")
  print(table(clin$SUBTYPE, useNA = "ifany"))
  log_msg("远处转移 M 分期计数（本队列几乎无解剖部位信息，仅 M0/M1）:")
  print(table(clin$M_STAGE, useNA = "ifany"))
  log_msg("淋巴结 N 分期计数（作为区域转移程度）:")
  print(table(clin$N_STAGE, useNA = "ifany"))
  log_msg("转移组合（无转移 / 仅淋巴结 / 远处）:")
  print(table(clin$MET_GROUP, useNA = "ifany"))

  note <- paste0(
    "TCGA BRCA Nature 2012（cBioPortal brca_tcga_pub）临床文件没有骨/肺/肝/脑等",
    "解剖转移部位。转移分析使用：METASTASIS (M0/M1)、NODES (N0-N3)、",
    "以及 No met / Lymph node only / Distant 三分类。"
  )
  writeLines(note, file.path(out_dir, "NOTE_metastasis_sites.txt"))

  score_mat <- NULL
  assoc_sub <- list()
  assoc_met <- list()
  skipped <- list()
  sizes <- list()

  for (gid in go_ids) {
    log_msg("==== ", gid, " ", GO_NAMES[gid], " ====")
    symbols <- get_genes_for_one_go(gid, cache_dir)
    in_expr <- intersect(symbols, rownames(expr))
    sizes[[gid]] <- data.frame(
      go_id = gid, go_name = unname(GO_NAMES[gid]),
      n_quickgo = length(symbols), n_in_expr = length(in_expr),
      stringsAsFactors = FALSE
    )
    if (length(in_expr) < min_genes) {
      log_msg(gid, " 表达矩阵中仅 ", length(in_expr), " 个基因，跳过（不并入其他通路）")
      skipped[[gid]] <- data.frame(
        go_id = gid, go_name = unname(GO_NAMES[gid]),
        n_in_expr = length(in_expr), reason = "fewer than min_genes in expression matrix",
        stringsAsFactors = FALSE
      )
    } else {
      score <- mean_z_score(expr, in_expr)
      if (is.null(score_mat)) {
        score_mat <- matrix(score, nrow = 1, dimnames = list(gid, names(score)))
      } else {
        score_mat <- rbind(score_mat, score[colnames(score_mat)])
        rownames(score_mat)[nrow(score_mat)] <- gid
      }

      sub_res <- test_group_vs_rest(score, clin[names(score), "SUBTYPE"])
      if (!is.null(sub_res)) {
        sub_res$go_id <- gid
        sub_res$go_name <- unname(GO_NAMES[gid])
        sub_res$phenotype <- "PAM50 subtype"
        assoc_sub[[gid]] <- sub_res
      }

      met_rows <- list(
        test_group_vs_rest(score, clin[names(score), "MET_GROUP"]),
        test_group_vs_rest(score, clin[names(score), "N_STAGE"]),
        test_binary(score, clin[names(score), "M_STAGE"], "M1", "M0"),
        test_binary(score, clin[names(score), "N_PLUS"], "N+", "N0")
      )
      met_rows <- met_rows[!vapply(met_rows, is.null, logical(1))]
      if (length(met_rows)) {
        met_res <- do.call(rbind, met_rows)
        met_res$go_id <- gid
        met_res$go_name <- unname(GO_NAMES[gid])
        met_res$phenotype <- "metastasis"
        assoc_met[[gid]] <- met_res
      }

      neg <- spearman_negative(score, expr, fdr = fdr)
      if (is.null(neg) || !nrow(neg)) {
        log_msg(gid, " 无 FDR<", fdr, " 的负相关基因")
        neg <- data.frame(
          gene = character(), rho = numeric(), n = integer(),
          p_value = numeric(), fdr = numeric(), stringsAsFactors = FALSE
        )
      } else {
        log_msg(gid, " 负相关基因: ", nrow(neg))
      }
      neg$go_id <- gid
      neg$go_name <- unname(GO_NAMES[gid])
      neg$in_this_go <- neg$gene %in% in_expr
      utils::write.table(
        neg, file.path(out_dir, "neg_genes", paste0(go_file_stem(gid), "_negative_genes.tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
  }

  size_df <- do.call(rbind, sizes)
  utils::write.table(size_df, file.path(out_dir, "go_gene_set_sizes.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  if (length(skipped)) {
    utils::write.table(do.call(rbind, skipped), file.path(out_dir, "skipped_go.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (!is.null(score_mat)) {
    score_out <- data.frame(go_id = rownames(score_mat), score_mat, check.names = FALSE)
    utils::write.table(score_out, file.path(out_dir, "go_sample_scores.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }

  add_fdr <- function(df) {
    if (is.null(df) || !nrow(df)) return(df)
    pair <- df$comparison != "overall Kruskal-Wallis"
    df$fdr <- NA_real_
    if (any(pair)) df$fdr[pair] <- stats::p.adjust(df$p_value[pair], method = "BH")
    df
  }

  if (length(assoc_sub)) {
    sub_df <- add_fdr(do.call(rbind, assoc_sub))
    rownames(sub_df) <- NULL
    utils::write.table(sub_df, file.path(out_dir, "association_subtype.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
    plot_bubble(
      sub_df,
      "GO pathway scores vs PAM50 subtype (one-vs-rest)",
      file.path(out_dir, "bubble_subtype.pdf"),
      y_order = c("Luminal A", "Luminal B", "HER2-enriched", "Basal-like",
                  "Normal-like", "Claudin-low")
    )
  }

  if (length(assoc_met)) {
    met_df <- add_fdr(do.call(rbind, assoc_met))
    rownames(met_df) <- NULL
    utils::write.table(met_df, file.path(out_dir, "association_metastasis.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
    bubble_met <- met_df[met_df$group %in% c(
      "M1", "N+", "N1", "N2", "N3",
      "Distant (M1)", "Lymph node only (M0 N+)", "No met (M0 N0)"
    ), , drop = FALSE]
    plot_bubble(
      bubble_met,
      "GO pathway scores vs metastasis (M / N / combined)",
      file.path(out_dir, "bubble_metastasis.pdf"),
      y_order = c("M1", "N+", "N1", "N2", "N3",
                  "Distant (M1)", "Lymph node only (M0 N+)", "No met (M0 N0)")
    )
  }

  neg_files <- list.files(file.path(out_dir, "neg_genes"), pattern = "_negative_genes\\.tsv$",
                          full.names = TRUE)
  if (length(neg_files)) {
    all_neg <- do.call(rbind, lapply(neg_files, function(f) {
      df <- utils::read.delim(f, stringsAsFactors = FALSE)
      if (!nrow(df)) return(NULL)
      df
    }))
    if (!is.null(all_neg) && nrow(all_neg)) {
      utils::write.table(all_neg, file.path(out_dir, "neg_genes_all.tsv"),
                         sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }

  log_msg("完成。结果目录: ", out_dir)
  invisible(out_dir)
}

main <- function() {
  args <- parse_args()
  go_ids <- GO_IDS
  data_dir <- args$data_dir
  out_dir <- args$out_dir
  cache_dir <- args$cache_dir
  if (isTRUE(args$demo)) {
    demo_root <- file.path(tempdir(), "tcga_brca_demo")
    demo <- write_demo_data(demo_root)
    data_dir <- demo$root
    out_dir <- file.path(demo_root, "results")
    cache_dir <- file.path(demo_root, "cache", "go_genes")
    go_ids <- demo$go_ids
    log_msg("DEMO 模式，数据写在: ", data_dir)
  }
  run_analysis(
    data_dir = data_dir,
    out_dir = out_dir,
    cache_dir = cache_dir,
    go_ids = go_ids,
    fdr = args$fdr,
    min_genes = as.integer(args$min_genes),
    demo = isTRUE(args$demo)
  )
}

if (sys.nframe() == 0L || identical(Sys.getenv("R_SCRIPT_RUN"), "1")) {
  main()
} else if (!interactive() && length(commandArgs(trailingOnly = TRUE))) {
  main()
} else if (interactive()) {
  main()
}
