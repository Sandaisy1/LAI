#!/usr/bin/env Rscript
# =============================================================================
# TCGA-BRCA：每个 GO 通路单独 vs 临床，并寻找与该通路负相关的基因
#
# 硬性约束：GO_TERMS 中每一个 ID 独立取基因、独立打分、独立做临床相关与
# 负相关基因筛选。禁止把所有 GO 基因并成一个集合再分析一次。
#
# 用法：
#   Rscript R/run_go_tcga_brca.R --data-dir /path/to/data --out-dir ./results
#   Rscript R/run_go_tcga_brca.R --demo --out-dir ./results_demo
#   Rscript R/run_go_tcga_brca.R --install-deps   # 仅安装依赖
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANGUAGE = "en")

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  if (file.exists("R/00_config.R")) return(normalizePath("R"))
  normalizePath(getwd())
}

SCRIPT_DIR <- get_script_dir()
ROOT_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."))
source(file.path(SCRIPT_DIR, "00_config.R"), local = FALSE)

# -----------------------------------------------------------------------------
# 命令行
# -----------------------------------------------------------------------------
parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  cfg <- list(
    data_dir = file.path(ROOT_DIR, "data"),
    out_dir = file.path(ROOT_DIR, "results"),
    gmt = NULL,
    demo = FALSE,
    install_deps = FALSE,
    method = "zmean" # zmean 或 ssgsea
  )
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("--demo", "-demo")) {
      cfg$demo <- TRUE
    } else if (a %in% c("--install-deps", "--install_deps")) {
      cfg$install_deps <- TRUE
    } else if (a %in% c("--data-dir", "--data_dir") && i < length(args)) {
      i <- i + 1L; cfg$data_dir <- args[[i]]
    } else if (a %in% c("--out-dir", "--out_dir") && i < length(args)) {
      i <- i + 1L; cfg$out_dir <- args[[i]]
    } else if (a == "--gmt" && i < length(args)) {
      i <- i + 1L; cfg$gmt <- args[[i]]
    } else if (a == "--method" && i < length(args)) {
      i <- i + 1L; cfg$method <- args[[i]]
    } else if (a %in% c("-h", "--help")) {
      cat(paste0(
        "Usage: Rscript R/run_go_tcga_brca.R [options]\n",
        "  --data-dir DIR     含 TCGA-BRCA.* 与 gencode 注释的目录\n",
        "  --out-dir DIR      输出目录（每个 GO 一个子文件夹）\n",
        "  --gmt FILE         可选 GMT，当未安装 org.Hs.eg.db 时按 GO ID 取基因\n",
        "  --method zmean|ssgsea  通路活性方法（ssgsea 需要 GSVA）\n",
        "  --demo             不读真实数据，用模拟数据验证“每个 GO 独立”\n",
        "  --install-deps     安装 CRAN/Bioconductor 依赖后退出\n"
      ))
      quit(save = "no", status = 0)
    } else {
      stop("未知参数: ", a, "  （使用 --help）")
    }
    i <- i + 1L
  }
  cfg
}

# -----------------------------------------------------------------------------
# 依赖
# -----------------------------------------------------------------------------
CRAN_PKGS <- c("data.table", "survival")
BIOC_PKGS <- c("AnnotationDbi", "org.Hs.eg.db", "GO.db")
OPTIONAL_PKGS <- c("GSVA", "ggplot2")

install_deps <- function() {
  repos <- "https://cloud.r-project.org"
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = repos)
  }
  cran_need <- CRAN_PKGS[!vapply(CRAN_PKGS, requireNamespace, logical(1), quietly = TRUE)]
  if (length(cran_need)) install.packages(cran_need, repos = repos)
  bioc_need <- BIOC_PKGS[!vapply(BIOC_PKGS, requireNamespace, logical(1), quietly = TRUE)]
  if (length(bioc_need)) BiocManager::install(bioc_need, update = FALSE, ask = FALSE)
  opt_need <- OPTIONAL_PKGS[!vapply(OPTIONAL_PKGS, requireNamespace, logical(1), quietly = TRUE)]
  if (length(opt_need)) {
    message("可选包未装（不影响主流程）: ", paste(opt_need, collapse = ", "))
    try(BiocManager::install(opt_need, update = FALSE, ask = FALSE), silent = TRUE)
  }
  message("依赖安装步骤完成。")
}

need_pkg <- function(pkg, hard = TRUE) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok && hard) {
    stop("缺少 R 包 '", pkg, "'。请先运行: Rscript R/run_go_tcga_brca.R --install-deps")
  }
  ok
}

# -----------------------------------------------------------------------------
# IO 与 ID 处理
# -----------------------------------------------------------------------------
harmonize_barcode <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("\\.", "-", x)
  x <- gsub("^X(?=TCGA-)", "", x, perl = TRUE)
  x
}

tcga_patient <- function(x) substr(harmonize_barcode(x), 1L, 12L)
tcga_sample15 <- function(x) substr(harmonize_barcode(x), 1L, 15L)
tcga_sample_type <- function(x) substr(harmonize_barcode(x), 14L, 15L)
strip_ensembl_version <- function(x) sub("\\.[0-9]+$", "", as.character(x))

find_data_file <- function(data_dir, key) {
  if (!dir.exists(data_dir)) {
    stop("数据目录不存在: ", data_dir)
  }
  files <- list.files(data_dir, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  files <- files[!grepl("/\\.", files)]
  hit <- files[grepl(key, basename(files), ignore.case = TRUE)]
  hit <- hit[!grepl("\\.(r|R|md|pdf|png)$", hit)]
  if (!length(hit)) {
    stop("在 ", data_dir, " 中找不到文件名包含 '", key, "' 的数据。")
  }
  if (length(hit) > 1L) {
    hit <- hit[order(nchar(basename(hit)))]
    message("匹配到多个文件，使用: ", hit[[1]])
  }
  hit[[1]]
}

fread_auto <- function(path) {
  if (need_pkg("data.table", hard = FALSE)) {
    return(data.table::fread(path, data.table = FALSE, check.names = FALSE))
  }
  con_fun <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile else file
  con <- con_fun(path, open = "rt")
  on.exit(close(con), add = TRUE)
  utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE)
}

is_tcga_id <- function(x) {
  grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", harmonize_barcode(x))
}

matrix_from_first_col <- function(df) {
  if (ncol(df) < 2L) stop("矩阵文件列数不足")
  ids <- as.character(df[[1]])
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- ids
  colnames(mat) <- harmonize_barcode(colnames(mat))
  mat
}

# Xena 有时是 feature × sample（列名为 TCGA barcode）
orient_feature_sample <- function(df, prefer_sample_cols = TRUE) {
  cn <- colnames(df)
  rn_guess <- as.character(df[[1]])
  n_col_tcga <- sum(is_tcga_id(cn))
  n_row_tcga <- sum(is_tcga_id(rn_guess))
  if (prefer_sample_cols && n_col_tcga >= 5L && n_col_tcga > n_row_tcga) {
    return(matrix_from_first_col(df))
  }
  if (n_row_tcga >= 5L && n_row_tcga > n_col_tcga) {
    rn <- harmonize_barcode(rn_guess)
    mat <- as.matrix(df[, -1, drop = FALSE])
    storage.mode(mat) <- "numeric"
    rownames(mat) <- rn
    return(t(mat))
  }
  matrix_from_first_col(df)
}

read_clinical_table <- function(path) {
  df <- fread_auto(path)
  cn <- colnames(df)
  n_col_tcga <- sum(is_tcga_id(cn))
  if (n_col_tcga >= 5L) {
    feat <- as.character(df[[1]])
    mat <- as.matrix(df[, -1, drop = FALSE])
    out <- as.data.frame(t(mat), stringsAsFactors = FALSE)
    colnames(out) <- feat
    out$sampleID <- rownames(out)
    rownames(out) <- NULL
    return(out)
  }
  id_col <- grep("sample|submitter|barcode|patient", cn, ignore.case = TRUE, value = TRUE)[1]
  first_is_tcga <- is_tcga_id(as.character(df[[1]])[seq_len(min(20L, nrow(df)))])
  if ("sampleID" %in% colnames(df)) {
    df$sampleID <- as.character(df$sampleID)
  } else if (mean(first_is_tcga) >= 0.5) {
    df$sampleID <- as.character(df[[1]])
  } else if (length(id_col) && !is.na(id_col)) {
    df$sampleID <- as.character(df[[id_col]])
  } else {
    df$sampleID <- as.character(df[[1]])
  }
  df
}

read_gencode <- function(path) {
  df <- fread_auto(path)
  cn <- tolower(colnames(df))
  colnames(df) <- cn
  # 也兼容原始 GTF
  if (ncol(df) >= 9 && (identical(cn[1], "v1") || grepl("gtf$", path) || any(cn %in% c("gene_id", "geneid")))) {
    if ("gene_id" %in% cn && "gene_name" %in% cn) {
      out <- data.frame(
        gene_id = strip_ensembl_version(df$gene_id),
        gene_name = as.character(df$gene_name),
        gene_type = if ("gene_type" %in% cn) df$gene_type else if ("gene_biotype" %in% cn) df$gene_biotype else NA_character_,
        stringsAsFactors = FALSE
      )
      return(unique(out))
    }
  }
  id_col <- which(cn %in% c("gene_id", "geneid", "ensembl_id", "ensembl", "id"))[1]
  name_col <- which(cn %in% c("gene_name", "symbol", "gene_symbol", "hgnc_symbol"))[1]
  type_col <- which(cn %in% c("gene_type", "gene_biotype", "biotype", "type"))[1]
  if (is.na(id_col)) id_col <- 1L
  if (is.na(name_col)) name_col <- min(2L, ncol(df))
  data.frame(
    gene_id = strip_ensembl_version(df[[id_col]]),
    gene_name = as.character(df[[name_col]]),
    gene_type = if (!is.na(type_col)) as.character(df[[type_col]]) else NA_character_,
    stringsAsFactors = FALSE
  )
}

collapse_duplicate_genes <- function(mat) {
  id <- strip_ensembl_version(rownames(mat))
  if (!anyDuplicated(id)) {
    rownames(mat) <- id
    return(mat)
  }
  message("合并重复 Ensembl ID（去版本号后），保留均值最高的一行")
  means <- rowMeans(mat, na.rm = TRUE)
  keep <- !duplicated(id[order(-means)])
  ord <- order(-means)
  mat <- mat[ord, , drop = FALSE]
  id <- id[ord]
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- id[keep]
  mat
}

filter_primary_tumor <- function(mat, sample_type_prefix) {
  if (is.null(sample_type_prefix) || !nzchar(sample_type_prefix)) return(mat)
  st <- tcga_sample_type(colnames(mat))
  keep <- st == sample_type_prefix
  if (!any(keep) && all(st %in% c("", "NA"))) {
    message("无法从 barcode 解析样本类型，保留全部样本")
    return(mat)
  }
  if (!any(keep)) {
    warning("没有样本类型为 ", sample_type_prefix, " 的列，保留全部样本")
    return(mat)
  }
  mat[, keep, drop = FALSE]
}

maybe_log2_fpkm <- function(mat) {
  mx <- suppressWarnings(max(mat, na.rm = TRUE))
  if (is.finite(mx) && mx > 100) {
    message("FPKM 最大值 ", signif(mx, 3), "，进行 log2(x+1) 转换")
    mat <- log2(mat + 1)
  }
  mat
}

# -----------------------------------------------------------------------------
# GO -> 基因（每个 GO ID 单独查询，绝不做 union）
# -----------------------------------------------------------------------------
read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  out <- list()
  for (ln in lines) {
    sp <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(sp) < 3L) next
    id <- sp[[1]]
    genes <- unique(sp[-c(1, 2)])
    genes <- genes[nzchar(genes)]
    out[[id]] <- genes
  }
  out
}

go_term_name <- function(go_id, fallback) {
  if (need_pkg("GO.db", hard = FALSE) && need_pkg("AnnotationDbi", hard = FALSE)) {
    got <- try(AnnotationDbi::Term(GO.db::GOTERM[[go_id]]), silent = TRUE)
    if (!inherits(got, "try-error") && length(got) && !is.na(got)) return(unname(got))
  }
  fallback
}

genes_for_one_go <- function(go_id, gmt = NULL, keytype = ANALYSIS_PARAMS$go_keytype) {
  # 返回 data.frame(ensembl, symbol)，仅该 GO（及 GOALL 子注释），不含其他所列 GO
  if (!is.null(gmt) && go_id %in% names(gmt)) {
    g <- unique(gmt[[go_id]])
    return(data.frame(ensembl = NA_character_, symbol = g, stringsAsFactors = FALSE))
  }
  if (!need_pkg("org.Hs.eg.db", hard = FALSE) || !need_pkg("AnnotationDbi", hard = FALSE)) {
    stop("无法映射 ", go_id, "：请安装 org.Hs.eg.db 或提供包含该 GO ID 的 --gmt 文件")
  }
  db <- org.Hs.eg.db::org.Hs.eg.db
  kt <- keytype
  keys_have <- try(AnnotationDbi::keys(db, keytype = kt), silent = TRUE)
  if (inherits(keys_have, "try-error") || !(go_id %in% keys_have)) {
    kt <- "GO"
  }
  df <- try(
    AnnotationDbi::select(db, keys = go_id, columns = c("ENSEMBL", "SYMBOL"), keytype = kt),
    silent = TRUE
  )
  if (inherits(df, "try-error") || is.null(df) || !nrow(df)) {
    return(data.frame(ensembl = character(), symbol = character(), stringsAsFactors = FALSE))
  }
  data.frame(
    ensembl = strip_ensembl_version(df$ENSEMBL),
    symbol = as.character(df$SYMBOL),
    stringsAsFactors = FALSE
  )
}

map_go_genes_to_expr <- function(go_genes_df, expr, gencode) {
  ens <- unique(go_genes_df$ensembl[!is.na(go_genes_df$ensembl) & nzchar(go_genes_df$ensembl)])
  sym <- unique(go_genes_df$symbol[!is.na(go_genes_df$symbol) & nzchar(go_genes_df$symbol)])
  hit <- character()
  if (length(ens)) hit <- c(hit, intersect(ens, rownames(expr)))
  if (!is.null(gencode) && length(sym)) {
    g2 <- gencode[gencode$gene_name %in% sym, , drop = FALSE]
    hit <- c(hit, intersect(g2$gene_id, rownames(expr)))
  }
  # 表达矩阵可能是 gene symbol
  if (length(sym)) hit <- c(hit, intersect(sym, rownames(expr)))
  unique(hit)
}

# -----------------------------------------------------------------------------
# 通路活性：只使用传入的该 GO 基因
# -----------------------------------------------------------------------------
pathway_score_zmean <- function(expr, genes) {
  genes <- intersect(genes, rownames(expr))
  sub <- expr[genes, , drop = FALSE]
  z <- t(scale(t(sub)))
  z[!is.finite(z)] <- 0
  colMeans(z, na.rm = TRUE)
}

pathway_score_ssgsea <- function(expr, genes, go_id) {
  if (!need_pkg("GSVA", hard = FALSE)) {
    message("未安装 GSVA，", go_id, " 改用 z-score 均值")
    return(pathway_score_zmean(expr, genes))
  }
  glist <- list()
  glist[[go_id]] <- intersect(genes, rownames(expr))
  mat <- as.matrix(expr)
  sc <- try({
    if (exists("ssgseaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
      GSVA::gsva(GSVA::ssgseaParam(mat, glist), verbose = FALSE)
    } else {
      GSVA::gsva(mat, glist, method = "ssgsea", verbose = FALSE)
    }
  }, silent = TRUE)
  if (inherits(sc, "try-error")) {
    message("ssGSEA 失败，", go_id, " 改用 z-score 均值: ", sc)
    return(pathway_score_zmean(expr, genes))
  }
  sc_vec <- as.numeric(sc[1, ])
  names(sc_vec) <- colnames(sc)
  sc_vec
}

score_one_go <- function(expr, genes, go_id, method = "zmean") {
  genes <- intersect(genes, rownames(expr))
  if (length(genes) < ANALYSIS_PARAMS$min_go_genes) {
    return(NULL)
  }
  if (identical(method, "ssgsea")) {
    sc <- pathway_score_ssgsea(expr, genes, go_id)
  } else {
    sc <- pathway_score_zmean(expr, genes)
  }
  sc
}

# -----------------------------------------------------------------------------
# 样本对齐
# -----------------------------------------------------------------------------
align_named_to_table <- function(score, tbl, id_candidates = c("sampleID", "sample", "bcr_patient_barcode")) {
  if (is.null(tbl) || !nrow(tbl)) return(NULL)
  tbl <- tbl
  if (!"sampleID" %in% colnames(tbl)) {
    for (cc in id_candidates) {
      if (cc %in% colnames(tbl)) {
        tbl$sampleID <- as.character(tbl[[cc]])
        break
      }
    }
  }
  if (!"sampleID" %in% colnames(tbl)) tbl$sampleID <- as.character(tbl[[1]])
  tbl$sampleID <- harmonize_barcode(tbl$sampleID)
  tbl$patient_id <- tcga_patient(tbl$sampleID)
  tbl$sample15 <- tcga_sample15(tbl$sampleID)
  sc <- data.frame(
    sampleID = names(score),
    pathway_score = as.numeric(score),
    patient_id = tcga_patient(names(score)),
    sample15 = tcga_sample15(names(score)),
    stringsAsFactors = FALSE
  )
  m <- merge(sc, tbl, by = "sample15", suffixes = c("", "_clin"))
  if (nrow(m) >= 10L) return(m)
  m2 <- merge(sc, tbl, by = "patient_id", suffixes = c("", "_clin"))
  m2
}

# -----------------------------------------------------------------------------
# 临床关联：针对“这一条”通路得分
# -----------------------------------------------------------------------------
is_id_like <- function(x, name) {
  grepl("id$|barcode|uuid|submitter|bcr_patient|sample", name, ignore.case = TRUE) ||
    (is.character(x) && mean(is_tcga_id(x), na.rm = TRUE) > 0.5)
}

summarize_clinical_one <- function(merged, go_id, go_name) {
  skip <- c("sampleID", "sampleID_clin", "patient_id", "patient_id_clin",
            "sample15", "sample15_clin", "pathway_score")
  vars <- setdiff(colnames(merged), skip)
  vars <- vars[!grepl("^(OS|DSS|PFI|DFI)($|[._])", vars, ignore.case = TRUE)]
  rows <- list()
  for (v in vars) {
    x <- merged[[v]]
    if (is_id_like(x, v)) next
    if (is.list(x)) next
    x_num <- suppressWarnings(as.numeric(as.character(x)))
    numeric_frac <- mean(is.finite(x_num) | is.na(x))
    n_uniq <- length(unique(x[!is.na(x) & x != "" & x != "NA"]))
    if (n_uniq <= 1L) next
    if (mean(is.na(x) | x == "" | x == "NA") > 0.6) next
    if (n_uniq > 30L && !(is.numeric(x) || (is.finite(numeric_frac) && numeric_frac > 0.8))) next

    use_numeric <- is.numeric(x) || (n_uniq > 8L && mean(is.finite(x_num)) > 0.8)
    if (use_numeric) {
      ok <- is.finite(x_num) & is.finite(merged$pathway_score)
      if (sum(ok) < 10L) next
      ct <- suppressWarnings(stats::cor.test(merged$pathway_score[ok], x_num[ok], method = "spearman", exact = FALSE))
      rows[[length(rows) + 1L]] <- data.frame(
        go_id = go_id, go_name = go_name, variable = v, test = "spearman",
        n = sum(ok), statistic = unname(ct$estimate), p_value = ct$p.value,
        extra = paste0("rho=", signif(unname(ct$estimate), 4)),
        stringsAsFactors = FALSE
      )
    } else {
      f <- factor(as.character(x))
      f <- droplevels(f[!(is.na(f) | f %in% c("", "NA", "Unknown", "unknown", "[Not Available]"))])
      df <- merged[as.character(x) %in% levels(f), , drop = FALSE]
      f <- factor(as.character(df[[v]]))
      if (nlevels(f) < 2L || nrow(df) < 10L) next
      if (nlevels(f) == 2L) {
        wt <- stats::wilcox.test(pathway_score ~ f, data = transform(df, f = f))
        rows[[length(rows) + 1L]] <- data.frame(
          go_id = go_id, go_name = go_name, variable = v, test = "wilcoxon",
          n = nrow(df), statistic = unname(wt$statistic), p_value = wt$p.value,
          extra = paste(levels(f), collapse = " vs "),
          stringsAsFactors = FALSE
        )
      } else {
        kt <- stats::kruskal.test(pathway_score ~ f, data = transform(df, f = f))
        rows[[length(rows) + 1L]] <- data.frame(
          go_id = go_id, go_name = go_name, variable = v, test = "kruskal",
          n = nrow(df), statistic = unname(kt$statistic), p_value = kt$p.value,
          extra = paste0("k=", nlevels(f)),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) {
    return(data.frame(
      go_id = go_id, go_name = go_name, variable = NA_character_, test = NA_character_,
      n = NA_integer_, statistic = NA_real_, p_value = NA_real_, fdr = NA_real_, extra = NA_character_,
      stringsAsFactors = FALSE
    )[0, ])
  }
  out <- do.call(rbind, rows)
  out$fdr <- p.adjust(out$p_value, method = "BH")
  out[order(out$p_value), ]
}

cox_one_endpoint <- function(merged, time_col, event_col, go_id, go_name, endpoint) {
  if (!all(c(time_col, event_col) %in% colnames(merged))) return(NULL)
  if (!need_pkg("survival", hard = FALSE)) return(NULL)
  time <- suppressWarnings(as.numeric(merged[[time_col]]))
  event <- suppressWarnings(as.numeric(merged[[event_col]]))
  ok <- is.finite(time) & time > 0 & event %in% c(0, 1) & is.finite(merged$pathway_score)
  if (sum(ok) < 20L || sum(event[ok] == 1) < 5L) return(NULL)
  d <- data.frame(time = time[ok], event = event[ok], score = merged$pathway_score[ok])
  fit <- try(survival::coxph(survival::Surv(time, event) ~ score, data = d), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  s <- summary(fit)
  hr <- unname(s$coefficients[1, "exp(coef)"])
  pcox <- unname(s$coefficients[1, "Pr(>|z|)"])
  ci <- s$conf.int[1, c("lower .95", "upper .95")]
  d$group <- factor(ifelse(d$score >= stats::quantile(d$score, ANALYSIS_PARAMS$km_quantile, na.rm = TRUE),
                           "High", "Low"), levels = c("Low", "High"))
  sdif <- try(survival::survdiff(survival::Surv(time, event) ~ group, data = d), silent = TRUE)
  km_p <- NA_real_
  if (!inherits(sdif, "try-error")) {
    km_p <- 1 - pchisq(sdif$chisq, length(sdif$n) - 1)
  }
  data.frame(
    go_id = go_id, go_name = go_name, endpoint = endpoint,
    n = nrow(d), n_events = sum(d$event == 1),
    cox_hr = hr, cox_ci_low = unname(ci[1]), cox_ci_high = unname(ci[2]),
    cox_p = pcox, km_logrank_p = km_p,
    stringsAsFactors = FALSE
  )
}

survival_one <- function(merged, go_id, go_name, out_pdf = NULL) {
  ends <- list(
    OS  = c("OS.time", "OS"),
    DSS = c("DSS.time", "DSS"),
    PFI = c("PFI.time", "PFI"),
    DFI = c("DFI.time", "DFI")
  )
  # 宽松列名
  alt <- list(
    OS = list(time = c("OS.time", "OS_time", "overall_survival_time", "days_to_death"),
              event = c("OS", "OS_status", "vital_status"))
  )
  rows <- lapply(names(ends), function(ep) {
    cox_one_endpoint(merged, ends[[ep]][1], ends[[ep]][2], go_id, go_name, ep)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    return(data.frame(
      go_id = character(), go_name = character(), endpoint = character(),
      n = integer(), n_events = integer(), cox_hr = numeric(),
      cox_ci_low = numeric(), cox_ci_high = numeric(), cox_p = numeric(),
      km_logrank_p = numeric(), stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  if (!is.null(out_pdf) && "OS.time" %in% colnames(merged) && "OS" %in% colnames(merged) &&
      need_pkg("survival", hard = FALSE)) {
    time <- suppressWarnings(as.numeric(merged$OS.time))
    event <- suppressWarnings(as.numeric(merged$OS))
    ok <- is.finite(time) & time > 0 & event %in% c(0, 1) & is.finite(merged$pathway_score)
    if (sum(ok) >= 20L) {
      d <- data.frame(time = time[ok], event = event[ok], score = merged$pathway_score[ok])
      d$group <- factor(ifelse(d$score >= median(d$score), "High", "Low"), levels = c("Low", "High"))
      fit <- survival::survfit(survival::Surv(time, event) ~ group, data = d)
      grDevices::pdf(out_pdf, width = 6, height = 5)
      plot(fit, col = c("#3B75AF", "#C44E52"), xlab = "Time", ylab = "OS probability",
           main = paste0(go_id, "\n", go_name), lwd = 2)
      legend("topright", legend = c("Low score", "High score"), col = c("#3B75AF", "#C44E52"), lty = 1, lwd = 2)
      grDevices::dev.off()
    }
  }
  out
}

# -----------------------------------------------------------------------------
# 负相关基因：全基因组 vs 这一条通路得分（不是 vs 所有 GO 的综合得分）
# -----------------------------------------------------------------------------
spearman_vs_score <- function(mat, score) {
  common <- intersect(colnames(mat), names(score))
  mat <- mat[, common, drop = FALSE]
  score <- score[common]
  n <- length(score)
  score_rank <- rank(score, ties.method = "average")
  # 行=特征, 列=样本 → 对每行 rank 后与 score_rank 做 Pearson
  feat_rank <- t(apply(mat, 1, function(v) {
    v[!is.finite(v)] <- NA
    if (sum(is.finite(v)) < 8L) return(rep(NA_real_, length(v)))
    rank(v, ties.method = "average", na.last = "keep")
  }))
  # 处理 NA：对每个特征用成对完整观测近似（此处要求缺失很少）
  cors <- as.numeric(stats::cor(t(feat_rank), score_rank, use = "pairwise.complete.obs"))
  names(cors) <- rownames(mat)
  den <- 1 - cors^2
  den[den < 1e-12] <- 1e-12
  tstat <- cors * sqrt((n - 2) / den)
  p <- 2 * stats::pt(-abs(tstat), df = pmax(n - 2, 1))
  p[!is.finite(cors)] <- NA
  data.frame(
    feature = names(cors),
    spearman_rho = cors,
    p_value = p,
    stringsAsFactors = FALSE
  )
}

negative_features_one <- function(mat, score, go_id, go_name, member_ids, gencode = NULL, feature_kind = "gene") {
  if (is.null(mat) || is.null(score) || nrow(mat) == 0L) {
    return(NULL)
  }
  tab <- spearman_vs_score(mat, score)
  tab$fdr <- p.adjust(tab$p_value, method = "BH")
  tab$go_id <- go_id
  tab$go_name <- go_name
  tab$in_query_go <- tab$feature %in% member_ids
  tab$feature_kind <- feature_kind
  if (!is.null(gencode) && identical(feature_kind, "gene")) {
    map <- gencode[match(tab$feature, gencode$gene_id), , drop = FALSE]
    tab$gene_symbol <- map$gene_name
    tab$gene_type <- map$gene_type
  } else if (identical(feature_kind, "gene")) {
    tab$gene_symbol <- tab$feature
    tab$gene_type <- NA_character_
  }
  keep <- is.finite(tab$spearman_rho) &
    tab$spearman_rho < 0 &
    tab$spearman_rho <= ANALYSIS_PARAMS$min_spearman_rho_neg &
    is.finite(tab$fdr) &
    tab$fdr <= ANALYSIS_PARAMS$fdr_cutoff
  out <- tab[keep, , drop = FALSE]
  out[order(out$spearman_rho), ]
}

# -----------------------------------------------------------------------------
# 单个 GO 的完整分析（唯一入口；主循环只调用本函数）
# -----------------------------------------------------------------------------
analyze_one_go <- function(go_id, go_name, expr, clinical, survival, protein,
                           gencode, gmt, method, out_root) {
  message("========== 独立分析 ", go_id, " | ", go_name, " ==========")
  go_dir <- file.path(out_root, gsub(":", "_", go_id, fixed = TRUE))
  dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)

  gdf <- genes_for_one_go(go_id, gmt = gmt, keytype = ANALYSIS_PARAMS$go_keytype)
  member <- map_go_genes_to_expr(gdf, expr, gencode)
  utils::write.table(
    data.frame(go_id = go_id, go_name = go_name, mapped_gene = member, stringsAsFactors = FALSE),
    file = file.path(go_dir, "go_member_genes.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  if (length(member) < ANALYSIS_PARAMS$min_go_genes) {
    warning(go_id, " 映射到表达矩阵的基因数 = ", length(member), "，跳过")
    return(list(go_id = go_id, skipped = TRUE, n_member = length(member)))
  }

  score <- score_one_go(expr, member, go_id, method = method)
  if (is.null(score)) {
    return(list(go_id = go_id, skipped = TRUE, n_member = length(member)))
  }
  score_df <- data.frame(
    sampleID = names(score),
    patient_id = tcga_patient(names(score)),
    go_id = go_id,
    go_name = go_name,
    pathway_score = as.numeric(score),
    n_member_genes = length(member),
    stringsAsFactors = FALSE
  )
  utils::write.table(score_df, file = file.path(go_dir, "pathway_score.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)

  clin_merged <- align_named_to_table(score, clinical)
  if (!is.null(survival)) {
    surv_merged <- align_named_to_table(score, survival)
    if (!is.null(clin_merged) && !is.null(surv_merged)) {
      keep_surv <- setdiff(colnames(surv_merged), colnames(clin_merged))
      clin_merged <- merge(clin_merged, surv_merged[, c("sample15", keep_surv), drop = FALSE],
                           by = "sample15", all.x = TRUE)
    } else if (is.null(clin_merged)) {
      clin_merged <- surv_merged
    }
  }

  clin_res <- NULL
  surv_res <- NULL
  if (!is.null(clin_merged) && nrow(clin_merged) >= 10L) {
    clin_res <- summarize_clinical_one(clin_merged, go_id, go_name)
    utils::write.table(clin_res, file = file.path(go_dir, "clinical_association.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
    surv_res <- survival_one(clin_merged, go_id, go_name,
                             out_pdf = file.path(go_dir, "os_km.pdf"))
    if (nrow(surv_res)) {
      utils::write.table(surv_res, file = file.path(go_dir, "survival_cox.tsv"),
                         sep = "\t", quote = FALSE, row.names = FALSE)
    }
  } else {
    warning(go_id, " 未能与临床/生存表对齐到足够样本")
  }

  # 低表达过滤仅用于全基因组相关，不影响该 GO 成员打分
  frac <- rowMeans(is.finite(expr) & expr > 0)
  expr_f <- expr[frac >= ANALYSIS_PARAMS$min_expr_frac, , drop = FALSE]
  neg_genes <- negative_features_one(
    expr_f, score, go_id, go_name, member_ids = member,
    gencode = gencode, feature_kind = "gene"
  )
  neg_path <- file.path(go_dir, "negative_correlated_genes.tsv")
  if (!is.null(neg_genes) && nrow(neg_genes)) {
    utils::write.table(neg_genes, file = neg_path, sep = "\t", quote = FALSE, row.names = FALSE)
  } else if (!is.null(neg_genes)) {
    utils::write.table(neg_genes, file = neg_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  neg_prot <- NULL
  if (!is.null(protein)) {
    neg_prot <- negative_features_one(
      protein, score, go_id, go_name, member_ids = character(),
      gencode = NULL, feature_kind = "protein"
    )
    if (!is.null(neg_prot)) {
      utils::write.table(neg_prot, file = file.path(go_dir, "negative_correlated_proteins.tsv"),
                         sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }

  list(
    go_id = go_id, go_name = go_name, skipped = FALSE,
    n_member = length(member), n_samples = length(score),
    n_neg_genes = if (is.null(neg_genes)) 0L else nrow(neg_genes),
    n_neg_proteins = if (is.null(neg_prot)) 0L else nrow(neg_prot),
    clinical = clin_res, survival = surv_res, score = score
  )
}

# -----------------------------------------------------------------------------
# 真实数据加载
# -----------------------------------------------------------------------------
load_real_data <- function(data_dir) {
  expr_path <- find_data_file(data_dir, DATA_FILE_KEYS$expr)
  clin_path <- find_data_file(data_dir, DATA_FILE_KEYS$clinical)
  surv_path <- find_data_file(data_dir, DATA_FILE_KEYS$survival)
  genc_path <- find_data_file(data_dir, DATA_FILE_KEYS$gencode)
  prot_path <- try(find_data_file(data_dir, DATA_FILE_KEYS$protein), silent = TRUE)

  message("表达: ", expr_path)
  message("临床: ", clin_path)
  message("生存: ", surv_path)
  message("注释: ", genc_path)

  expr <- collapse_duplicate_genes(orient_feature_sample(fread_auto(expr_path)))
  expr <- filter_primary_tumor(expr, ANALYSIS_PARAMS$sample_type_prefix)
  expr <- maybe_log2_fpkm(expr)
  expr[!is.finite(expr)] <- NA

  gencode <- read_gencode(genc_path)
  gencode <- unique(gencode)
  if (isTRUE(ANALYSIS_PARAMS$protein_coding_only) && "gene_type" %in% colnames(gencode)) {
    pc <- gencode$gene_id[grepl("protein_coding", gencode$gene_type, ignore.case = TRUE)]
    if (length(pc) >= 1000L) {
      keep <- intersect(rownames(expr), pc)
      if (length(keep) >= 1000L) expr <- expr[keep, , drop = FALSE]
    }
  }

  clinical <- read_clinical_table(clin_path)
  survival <- read_clinical_table(surv_path)

  protein <- NULL
  if (!inherits(prot_path, "try-error")) {
    message("蛋白: ", prot_path)
    protein <- orient_feature_sample(fread_auto(prot_path))
    colnames(protein) <- harmonize_barcode(colnames(protein))
  }

  list(expr = expr, clinical = clinical, survival = survival, protein = protein, gencode = gencode)
}

# -----------------------------------------------------------------------------
# Demo：三个彼此不同的 GO，确认不会被合并
# -----------------------------------------------------------------------------
make_demo <- function() {
  set.seed(ANALYSIS_PARAMS$seed)
  n_s <- 40L
  samples <- sprintf("TCGA-A2-%04d-01A", seq_len(n_s))
  genes <- c(
    paste0("ENSG", sprintf("%011d", 1:15)),
    paste0("ENSG", sprintf("%011d", 101:130))
  )
  expr <- matrix(rnorm(length(genes) * n_s), nrow = length(genes), dimnames = list(genes, samples))
  # GO A 成员在部分样本上抬高，并与 age 正相关；背景基因与 A 负相关
  a_genes <- paste0("ENSG", sprintf("%011d", 1:5))
  b_genes <- paste0("ENSG", sprintf("%011d", 6:10))
  c_genes <- paste0("ENSG", sprintf("%011d", 11:15))
  age <- 40 + 30 * seq_len(n_s) / n_s
  expr[a_genes, ] <- expr[a_genes, ] + matrix(rep(scale(age), each = 5), nrow = 5)
  expr[b_genes, ] <- expr[b_genes, ] - matrix(rep(scale(age), each = 5), nrow = 5)
  anti <- paste0("ENSG", sprintf("%011d", 101:110))
  expr[anti, ] <- expr[anti, ] - matrix(rep(colMeans(expr[a_genes, , drop = FALSE]), each = 10), nrow = 10)

  gencode <- data.frame(
    gene_id = genes,
    gene_name = genes,
    gene_type = "protein_coding",
    stringsAsFactors = FALSE
  )
  clinical <- data.frame(
    sampleID = samples,
    age = age,
    stage = rep(c("I", "II", "III"), length.out = n_s),
    stringsAsFactors = FALSE
  )
  survival <- data.frame(
    sampleID = samples,
    OS = as.integer(rbinom(n_s, 1, 0.4)),
    OS.time = round(runif(n_s, 100, 3000)),
    stringsAsFactors = FALSE
  )
  gmt <- list()
  gmt[["GO:0007409"]] <- a_genes
  gmt[["GO:0007411"]] <- b_genes
  gmt[["GO:0031103"]] <- c_genes
  demo_terms <- c(
    "GO:0007409" = "axonogenesis",
    "GO:0007411" = "axon guidance",
    "GO:0031103" = "axon regeneration"
  )
  list(expr = expr, clinical = clinical, survival = survival, protein = NULL,
       gencode = gencode, gmt = gmt, go_terms = demo_terms)
}

assert_demo_independence <- function(results, demo) {
  # 三个 GO 必须都产生各自得分；A 与 B 的得分不应被做成同一个向量
  ids <- vapply(results, `[[`, character(1), "go_id")
  stopifnot(length(unique(ids)) == 3L)
  sc_a <- results[[which(ids == "GO:0007409")]]$score
  sc_b <- results[[which(ids == "GO:0007411")]]$score
  rho <- suppressWarnings(stats::cor(sc_a, sc_b, method = "spearman"))
  if (!is.na(rho) && rho > 0.95) {
    stop("Demo 失败：两个 GO 的通路得分几乎相同，说明发生了基因集合并")
  }
  message("Demo 自检通过：3 个 GO 独立打分（A vs B Spearman rho = ", signif(rho, 3), "）")
}

# -----------------------------------------------------------------------------
# 主流程
# -----------------------------------------------------------------------------
main <- function() {
  cfg <- parse_args()
  if (cfg$install_deps) {
    install_deps()
    return(invisible(TRUE))
  }
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(ANALYSIS_PARAMS$seed)

  gmt <- NULL
  go_terms <- GO_TERMS
  if (cfg$demo) {
    message("*** DEMO 模式：模拟数据，仅用于验证每个 GO 独立分析 ***")
    dat <- make_demo()
    gmt <- dat$gmt
    go_terms <- dat$go_terms
  } else {
    if (!need_pkg("data.table", hard = FALSE)) {
      message("建议安装 data.table 以加快 FPKM 读取")
    }
    dat <- load_real_data(cfg$data_dir)
    if (!is.null(cfg$gmt)) gmt <- read_gmt(cfg$gmt)
  }

  if (!cfg$demo && is.null(gmt)) {
    need_pkg("org.Hs.eg.db", hard = TRUE)
    need_pkg("AnnotationDbi", hard = TRUE)
  }

  results <- vector("list", length(go_terms))
  names(results) <- names(go_terms)

  # 关键：逐个 GO 调用 analyze_one_go，循环体内不累积基因集
  for (i in seq_along(go_terms)) {
    go_id <- names(go_terms)[[i]]
    go_name <- go_term_name(go_id, unname(go_terms[[i]]))
    results[[i]] <- analyze_one_go(
      go_id = go_id,
      go_name = go_name,
      expr = dat$expr,
      clinical = dat$clinical,
      survival = dat$survival,
      protein = dat$protein,
      gencode = dat$gencode,
      gmt = gmt,
      method = cfg$method,
      out_root = cfg$out_dir
    )
  }

  clin_all <- do.call(rbind, lapply(results, function(r) {
    if (isTRUE(r$skipped) || is.null(r$clinical)) return(NULL)
    r$clinical
  }))
  surv_all <- do.call(rbind, lapply(results, function(r) {
    if (isTRUE(r$skipped) || is.null(r$survival) || !nrow(r$survival)) return(NULL)
    r$survival
  }))
  overview <- do.call(rbind, lapply(results, function(r) {
    data.frame(
      go_id = r$go_id,
      go_name = if (!is.null(r$go_name)) r$go_name else GO_TERMS[[r$go_id]],
      skipped = isTRUE(r$skipped),
      n_member_genes = if (!is.null(r$n_member)) r$n_member else NA_integer_,
      n_samples = if (!is.null(r$n_samples)) r$n_samples else NA_integer_,
      n_neg_genes = if (!is.null(r$n_neg_genes)) r$n_neg_genes else NA_integer_,
      n_neg_proteins = if (!is.null(r$n_neg_proteins)) r$n_neg_proteins else NA_integer_,
      stringsAsFactors = FALSE
    )
  }))
  utils::write.table(overview, file = file.path(cfg$out_dir, "summary_by_go.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  if (!is.null(clin_all) && nrow(clin_all)) {
    utils::write.table(clin_all, file = file.path(cfg$out_dir, "summary_clinical_long.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (!is.null(surv_all) && nrow(surv_all)) {
    utils::write.table(surv_all, file = file.path(cfg$out_dir, "summary_survival_long.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }

  sink(file.path(cfg$out_dir, "sessionInfo.txt"))
  print(utils::sessionInfo())
  sink()

  if (cfg$demo) assert_demo_independence(results, dat)
  message("完成。每个 GO 的结果在: ", normalizePath(cfg$out_dir))
  invisible(results)
}

if (identical(environment(), globalenv()) && !interactive()) {
  main()
}
