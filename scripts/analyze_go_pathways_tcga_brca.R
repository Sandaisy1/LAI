#!/usr/bin/env Rscript
# =============================================================================
# 每个 GO 通路独立分析（禁止合并为“所有通路”）
#
# 1) 每个 GO 与临床 / 生存的关联，各自输出一份结果
# 2) 每个 GO 通路得分与表达矩阵中每个基因的 Spearman 负相关，各自输出一份结果
#
# 用法：
#   Rscript scripts/analyze_go_pathways_tcga_brca.R --data-dir data --out-dir results
#   Rscript scripts/analyze_go_pathways_tcga_brca.R --demo
#
# 依赖（真实数据模式，只需 CRAN，不需要 AnnotationDbi / GSVA）：
#   install.packages(c("data.table", "survival"))
# GO 基因集：QuickGO 接口（按每个 GO 单独下载，结果缓存到 data/cache/go_genes/）
# 通路得分：通路内基因的均值 z-score（替代 GSVA）
# =============================================================================

suppressPackageStartupMessages({
  options(stringsAsFactors = FALSE, warn = 1)
})

GO_IDS <- c(
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

GO_NAMES <- c(
  "GO:0023041" = "neuronal signal transduction",
  "GO:1904457" = "positive regulation of neuronal action potential",
  "GO:1904340" = "positive regulation of dopaminergic neuron differentiation",
  "GO:2001224" = "positive regulation of neuron migration",
  "GO:2001222" = "regulation of neuron migration",
  "GO:0019227" = "neuronal action potential propagation",
  "GO:0019228" = "neuronal action potential",
  "GO:1902847" = "regulation of neuronal signal transduction",
  "GO:0031102" = "neuron projection regeneration",
  "GO:0097492" = "sympathetic neuron axon guidance",
  "GO:0097491" = "sympathetic neuron projection guidance",
  "GO:0097374" = "sensory neuron axon guidance",
  "GO:0007158" = "neuron cell-cell adhesion",
  "GO:1902667" = "regulation of axon guidance",
  "GO:0031103" = "axon regeneration",
  "GO:0007411" = "axon guidance",
  "GO:0007409" = "axonogenesis"
)

DATA_STEMS <- c(
  clinical = "TCGA-BRCA.clinical",
  protein  = "TCGA-BRCA.protein",
  expr     = "TCGA-BRCA.star_fpkm",
  survival = "TCGA-BRCA.survival",
  gencode  = "gencode.v36.annotation.gtf.gene"
)

# -----------------------------------------------------------------------------
# 工具函数
# -----------------------------------------------------------------------------

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
}

parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  defaults <- list(
    data_dir = "data",
    out_dir = "results",
    fdr = 0.05,
    min_rho = 0,
    min_genes = 5,
    primary_tumor_only = TRUE,
    demo = FALSE,
    install_deps = FALSE
  )
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (key %in% c("--demo", "--install-deps", "--all-samples")) {
      if (key == "--demo") defaults$demo <- TRUE
      if (key == "--install-deps") defaults$install_deps <- TRUE
      if (key == "--all-samples") defaults$primary_tumor_only <- FALSE
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

ensure_dir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

strip_ensembl_version <- function(x) sub("\\.\\d+$", "", x)

patient_barcode <- function(x) substr(gsub("\\.", "-", x), 1, 12)

sample_type_code <- function(x) substr(gsub("\\.", "-", x), 14, 15)

is_tcga_id <- function(x) grepl("^TCGA-", gsub("\\.", "-", x))

safe_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

install_deps_if_needed <- function() {
  cran <- c("data.table", "survival")
  for (p in cran) {
    if (!safe_pkg(p)) install.packages(p, repos = "https://cloud.r-project.org")
  }
}

find_data_file <- function(stem, data_dir) {
  if (!dir.exists(data_dir)) return(NA_character_)
  all_files <- list.files(data_dir, recursive = TRUE, full.names = TRUE)
  all_files <- all_files[!grepl("\\.gitkeep$", all_files)]
  base <- basename(all_files)
  hit <- all_files[startsWith(base, stem) | grepl(stem, base, fixed = TRUE)]
  if (!length(hit)) return(NA_character_)
  hit[[1]]
}

fread_auto <- function(path) {
  if (safe_pkg("data.table")) {
    data.table::fread(path, data.table = FALSE, check.names = FALSE)
  } else {
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

pick_id_column <- function(df) {
  nms <- names(df)
  prefer <- c(
    "sampleID", "sample", "Sample", "SAMPLE",
    "bcr_patient_barcode", "submitter_id", "_PATIENT", "patient",
    "barcode", "ID", "id"
  )
  for (p in prefer) {
    if (p %in% nms) return(p)
  }
  scores <- vapply(nms, function(nm) {
    x <- as.character(df[[nm]])
    mean(is_tcga_id(x), na.rm = TRUE)
  }, numeric(1))
  if (max(scores, na.rm = TRUE) >= 0.5) return(names(scores)[which.max(scores)])
  nms[[1]]
}

# 读成 genes x samples 数值矩阵
read_expr_matrix <- function(path) {
  log_msg("读取表达: ", path)
  df <- fread_auto(path)
  id_col <- names(df)[[1]]
  row_ids <- as.character(df[[id_col]])
  col_ids <- names(df)[-1]
  row_tcga <- mean(is_tcga_id(row_ids), na.rm = TRUE)
  col_tcga <- mean(is_tcga_id(col_ids), na.rm = TRUE)
  num_df <- df[, -1, drop = FALSE]
  num_df[] <- lapply(num_df, function(x) suppressWarnings(as.numeric(x)))
  if (col_tcga >= 0.5 || (col_tcga >= row_tcga)) {
    mat <- as.matrix(num_df)
    rownames(mat) <- row_ids
    colnames(mat) <- col_ids
  } else {
    mat <- t(as.matrix(num_df))
    colnames(mat) <- row_ids
    rownames(mat) <- col_ids
  }
  storage.mode(mat) <- "double"
  mat
}

read_table_with_ids <- function(path) {
  log_msg("读取表: ", path)
  df <- fread_auto(path)
  id_col <- pick_id_column(df)
  ids <- as.character(df[[id_col]])
  df[[id_col]] <- NULL
  rownames(df) <- make.unique(ids)
  df
}

read_gencode <- function(path) {
  log_msg("读取基因注释: ", path)
  first <- readLines(path, n = 20)
  first <- first[!startsWith(first, "#")]
  if (!length(first)) stop("gencode 文件为空: ", path)
  n_tab <- length(strsplit(first[[1]], "\t", fixed = TRUE)[[1]])
  looks_gtf <- n_tab >= 9 && grepl("gene_id", first[[1]])
  if (looks_gtf) {
    df <- fread_auto(path)
    if (ncol(df) >= 9) {
      attr_col <- df[[9]]
      gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", attr_col)
      gene_name <- ifelse(
        grepl("gene_name", attr_col),
        sub('.*gene_name "([^"]+)".*', "\\1", attr_col),
        NA_character_
      )
      keep <- df[[3]] == "gene" | !("gene" %in% df[[3]])
      out <- data.frame(
        gene_id = strip_ensembl_version(gene_id[keep]),
        gene_name = gene_name[keep],
        stringsAsFactors = FALSE
      )
      return(unique(out))
    }
  }
  df <- fread_auto(path)
  nms <- tolower(names(df))
  id_i <- which(nms %in% c("gene_id", "ensembl", "ensembl_id", "geneid"))[1]
  name_i <- which(nms %in% c("gene_name", "symbol", "gene_symbol", "hgnc_symbol"))[1]
  if (is.na(id_i)) id_i <- 1L
  if (is.na(name_i)) name_i <- min(2L, ncol(df))
  data.frame(
    gene_id = strip_ensembl_version(as.character(df[[id_i]])),
    gene_name = as.character(df[[name_i]]),
    stringsAsFactors = FALSE
  )
}

keep_primary_tumor <- function(mat) {
  codes <- sample_type_code(colnames(mat))
  keep <- is.na(codes) | codes == "" | codes == "01"
  mat <- mat[, keep, drop = FALSE]
  pid <- patient_barcode(colnames(mat))
  mat[, !duplicated(pid), drop = FALSE]
}

align_samples <- function(score, table_df) {
  sc_names <- names(score)
  tbl_names <- rownames(table_df)
  map_sc <- patient_barcode(sc_names)
  map_tbl <- patient_barcode(tbl_names)
  common_p <- intersect(map_sc, map_tbl)
  if (!length(common_p)) return(list(score = numeric(), table = table_df[0, , drop = FALSE]))
  sc_idx <- match(common_p, map_sc)
  tbl_idx <- match(common_p, map_tbl)
  list(score = score[sc_idx], table = table_df[tbl_idx, , drop = FALSE])
}

# -----------------------------------------------------------------------------
# GO 基因集：每个 GO 单独取（QuickGO，不依赖 AnnotationDbi）
# -----------------------------------------------------------------------------

go_term_name <- function(go_id) {
  unname(GO_NAMES[go_id])
}

go_cache_path <- function(go_id, cache_dir) {
  file.path(cache_dir, paste0(gsub(":", "_", go_id), "_genes.tsv"))
}

read_cached_go_genes <- function(path) {
  if (!file.exists(path) || !file.info(path)$size) {
    return(NULL)
  }
  df <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  list(
    ensembl = unique(strip_ensembl_version(na.omit(as.character(df$ensembl)))),
    symbol = unique(na.omit(as.character(df$symbol)))
  )
}

write_cached_go_genes <- function(path, genes) {
  ensure_dir(dirname(path))
  sym <- unique(as.character(genes$symbol))
  ens <- unique(as.character(genes$ensembl))
  sym <- sym[!is.na(sym) & nzchar(sym)]
  ens <- ens[!is.na(ens) & nzchar(ens)]
  n <- max(length(sym), length(ens), 0L)
  df <- data.frame(
    symbol = if (n) c(sym, rep(NA_character_, n - length(sym))) else character(),
    ensembl = if (n) c(ens, rep(NA_character_, n - length(ens))) else character(),
    stringsAsFactors = FALSE
  )
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

extract_human_ensembl <- function(...) {
  txt <- paste(unlist(list(...)), collapse = " ")
  hits <- unlist(regmatches(txt, gregexpr("ENSG[0-9]+", txt)))
  unique(strip_ensembl_version(hits))
}

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
  if (!download_url(url, tmp, extra_headers = c(Accept = "text/tsv"))) {
    return(NULL)
  }
  first <- tryCatch(readLines(tmp, n = 1L, warn = FALSE), error = function(e) "")
  if (!length(first) || grepl("^\\s*\\{", first)) return(NULL)
  df <- utils::read.delim(tmp, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(df)) {
    return(list(ensembl = character(), symbol = character()))
  }
  nms <- toupper(gsub("[^A-Za-z0-9]+", "_", names(df)))
  names(df) <- nms
  sym_col <- intersect(c("SYMBOL", "GENE_PRODUCT_SYMBOL"), names(df))
  id_col <- intersect(c("GENE_PRODUCT_ID", "GENEPRODUCTID"), names(df))
  with_col <- intersect(c("WITH_FROM", "WITHFROM"), names(df))
  symbol <- if (length(sym_col)) as.character(df[[sym_col[[1]]]]) else character()
  symbol <- unique(symbol[!is.na(symbol) & nzchar(symbol) & symbol != "-"])
  ensembl <- extract_human_ensembl(
    if (length(id_col)) df[[id_col[[1]]]] else NULL,
    if (length(with_col)) df[[with_col[[1]]]] else NULL
  )
  list(ensembl = ensembl, symbol = symbol)
}

get_genes_for_one_go <- function(go_id, cache_dir) {
  cache_file <- go_cache_path(go_id, cache_dir)
  if (file.exists(cache_file)) {
    cached <- read_cached_go_genes(cache_file)
    if (!is.null(cached)) {
      log_msg(go_id, " 使用缓存基因集: ", cache_file)
      return(cached)
    }
  }
  log_msg(go_id, " 从 QuickGO 下载基因集（含 descendant terms）")
  fetched <- tryCatch(fetch_one_go_from_quickgo(go_id), error = function(e) {
    log_msg(go_id, " QuickGO 失败: ", conditionMessage(e))
    NULL
  })
  if (is.null(fetched)) {
    stop(
      "无法获取 ", go_id, " 的基因集（不使用 AnnotationDbi，且 QuickGO 下载失败）。\n",
      "请检查网络，或手动把 symbol/ensembl 两列写入: ", cache_file
    )
  }
  write_cached_go_genes(cache_file, fetched)
  fetched
}

map_symbols_via_gencode <- function(genes, gencode_df) {
  if (is.null(gencode_df) || !nrow(gencode_df)) return(genes)
  ens_from_sym <- gencode_df$gene_id[match(genes$symbol, gencode_df$gene_name)]
  genes$ensembl <- unique(c(genes$ensembl, strip_ensembl_version(na.omit(ens_from_sym))))
  genes
}

match_genes_to_expr <- function(gene_ids, expr_ids) {
  expr_ids <- as.character(expr_ids)
  raw <- unique(as.character(gene_ids))
  hit <- intersect(raw, expr_ids)
  if (length(hit)) return(hit)
  intersect(strip_ensembl_version(raw), strip_ensembl_version(expr_ids))
}

# -----------------------------------------------------------------------------
# 通路得分：一次只对一个 GO 的基因集打分
# -----------------------------------------------------------------------------

mean_z_score <- function(expr, genes) {
  g <- intersect(genes, rownames(expr))
  mat <- expr[g, , drop = FALSE]
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  stats::setNames(colMeans(z, na.rm = TRUE), colnames(expr))
}

one_go_pathway_score <- function(expr, genes, go_id) {
  g <- intersect(genes, rownames(expr))
  if (length(g) < 2L) return(NULL)
  # 替代 GSVA：每个基因跨样本 z-score 后取均值，得到该 GO 的样本得分
  mean_z_score(expr, g)
}

# -----------------------------------------------------------------------------
# 临床 / 生存：针对“这一条”通路得分
# -----------------------------------------------------------------------------

is_id_like_column <- function(nm, x) {
  nm_l <- tolower(nm)
  if (grepl("barcode|submitter|uuid|patient_id|sample_id|^sample$|^id$", nm_l)) return(TRUE)
  if (is.character(x) || is.factor(x)) {
    u <- unique(as.character(x[!is.na(x) & x != "" & x != "NA"]))
    if (length(u) > 20 && mean(is_tcga_id(u)) > 0.3) return(TRUE)
  }
  FALSE
}

clinical_association_one_go <- function(score, clinical) {
  al <- align_samples(score, clinical)
  score <- al$score
  clinical <- al$table
  if (length(score) < 10L) {
    return(data.frame(
      variable = character(), n = integer(), n_groups = integer(),
      test = character(), statistic = numeric(), p_value = numeric(),
      extra = character(), stringsAsFactors = FALSE
    ))
  }
  rows <- list()
  for (nm in names(clinical)) {
    x <- clinical[[nm]]
    if (is_id_like_column(nm, x)) next
    if (is.list(x) && !is.data.frame(x)) next
    x_num <- suppressWarnings(as.numeric(as.character(x)))
    numeric_ok <- mean(!is.na(x_num)) >= 0.5 && length(unique(x_num[!is.na(x_num)])) >= 5
    if (numeric_ok) {
      keep <- is.finite(score) & is.finite(x_num)
      if (sum(keep) < 10L) next
      ct <- suppressWarnings(stats::cor.test(score[keep], x_num[keep], method = "spearman", exact = FALSE))
      rows[[length(rows) + 1L]] <- data.frame(
        variable = nm, n = sum(keep), n_groups = NA_integer_,
        test = "spearman", statistic = unname(ct$estimate), p_value = ct$p.value,
        extra = sprintf("rho=%.4f", unname(ct$estimate)),
        stringsAsFactors = FALSE
      )
    } else {
      x_chr <- as.character(x)
      x_chr[x_chr %in% c("", "NA", "Unknown", "unknown", "[Not Available]", "not reported")] <- NA
      keep <- is.finite(score) & !is.na(x_chr)
      if (sum(keep) < 10L) next
      grp <- factor(x_chr[keep])
      tab <- table(grp)
      grp <- droplevels(grp[grp %in% names(tab)[tab >= 5]])
      sc <- score[keep][x_chr[keep] %in% levels(grp)]
      if (nlevels(grp) < 2L || length(sc) < 10L) next
      if (nlevels(grp) == 2L) {
        wt <- stats::wilcox.test(sc ~ grp)
        rows[[length(rows) + 1L]] <- data.frame(
          variable = nm, n = length(sc), n_groups = 2L,
          test = "wilcoxon", statistic = unname(wt$statistic), p_value = wt$p.value,
          extra = paste(names(table(grp)), table(grp), sep = "=", collapse = "; "),
          stringsAsFactors = FALSE
        )
      } else if (nlevels(grp) <= 12L) {
        kt <- stats::kruskal.test(sc ~ grp)
        rows[[length(rows) + 1L]] <- data.frame(
          variable = nm, n = length(sc), n_groups = nlevels(grp),
          test = "kruskal", statistic = unname(kt$statistic), p_value = kt$p.value,
          extra = paste(names(table(grp)), table(grp), sep = "=", collapse = "; "),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) {
    return(data.frame(
      variable = character(), n = integer(), n_groups = integer(),
      test = character(), statistic = numeric(), p_value = numeric(),
      extra = character(), stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  out$fdr <- stats::p.adjust(out$p_value, method = "BH")
  out[order(out$p_value), ]
}

survival_association_one_go <- function(score, survival) {
  empty <- data.frame(
    endpoint = character(), n = integer(), n_events = integer(),
    cox_hr = numeric(), cox_p = numeric(),
    logrank_p = numeric(), median_cutoff = numeric(),
    stringsAsFactors = FALSE
  )
  if (!safe_pkg("survival")) {
    log_msg("未安装 survival，跳过生存分析")
    return(empty)
  }
  al <- align_samples(score, survival)
  score <- al$score
  surv <- al$table
  if (length(score) < 20L) return(empty)
  nms <- names(surv)
  nms_l <- tolower(nms)
  endpoints <- c("OS", "DSS", "PFI", "DFI", "PFS", "DFS")
  rows <- list()
  for (ep in endpoints) {
    ev_i <- which(nms_l == tolower(ep) | nms_l == paste0(tolower(ep), "_event"))
    tm_i <- which(nms_l %in% c(paste0(tolower(ep), ".time"), paste0(tolower(ep), "_time"), paste0(tolower(ep), ".days")))
    if (!length(ev_i) || !length(tm_i)) next
    event <- suppressWarnings(as.numeric(as.character(surv[[ev_i[[1]]]])))
    time <- suppressWarnings(as.numeric(as.character(surv[[tm_i[[1]]]])))
    keep <- is.finite(score) & is.finite(time) & is.finite(event) & time > 0
    if (sum(keep) < 20L || sum(event[keep] == 1) < 5L) next
    sc <- score[keep]
    sdf <- data.frame(time = time[keep], event = event[keep], score = sc)
    fit <- tryCatch(survival::coxph(survival::Surv(time, event) ~ score, data = sdf), error = function(e) NULL)
    if (is.null(fit)) next
    sfit <- summary(fit)
    hr <- unname(sfit$coefficients[1, "exp(coef)"])
    cox_p <- unname(sfit$coefficients[1, "Pr(>|z|)"])
    cut <- stats::median(sc, na.rm = TRUE)
    sdf$group <- factor(ifelse(sc >= cut, "high", "low"), levels = c("low", "high"))
    lr <- tryCatch(survival::survdiff(survival::Surv(time, event) ~ group, data = sdf), error = function(e) NULL)
    logrank_p <- if (is.null(lr)) NA_real_ else stats::pchisq(lr$chisq, df = 1, lower.tail = FALSE)
    rows[[length(rows) + 1L]] <- data.frame(
      endpoint = ep, n = nrow(sdf), n_events = sum(sdf$event == 1),
      cox_hr = hr, cox_p = cox_p, logrank_p = logrank_p, median_cutoff = cut,
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(empty)
  do.call(rbind, rows)
}

# -----------------------------------------------------------------------------
# 负相关基因：通路得分 vs 表达矩阵中每一个基因（不含本通路基因）
# -----------------------------------------------------------------------------

spearman_vs_matrix <- function(score, mat) {
  common <- intersect(names(score), colnames(mat))
  score <- score[common]
  mat <- mat[, common, drop = FALSE]
  n <- length(score)
  # 每个基因（行）在样本间取秩，再与通路得分的秩做 Pearson = Spearman
  score_r <- rank(score, na.last = "keep")
  mat_r <- t(apply(mat, 1, function(v) rank(v, na.last = "keep")))
  if (nrow(mat) == 1L) mat_r <- matrix(mat_r, nrow = 1L, dimnames = dimnames(mat))
  sc <- as.numeric(scale(score_r))
  mr <- t(scale(t(mat_r)))
  mr[!is.finite(mr)] <- 0
  sc[!is.finite(sc)] <- 0
  rho <- as.numeric(mr %*% sc / (n - 1))
  names(rho) <- rownames(mat)
  tstat <- rho * sqrt((n - 2) / pmax(1e-16, 1 - rho^2))
  pval <- 2 * stats::pt(-abs(tstat), df = n - 2)
  data.frame(
    gene = rownames(mat),
    rho = rho,
    p_value = pval,
    stringsAsFactors = FALSE
  )
}

negative_genes_one_go <- function(score, expr, pathway_genes, fdr_cutoff, min_rho, symbol_map = NULL) {
  tab <- spearman_vs_matrix(score, expr)
  in_pw <- tab$gene %in% pathway_genes |
    strip_ensembl_version(tab$gene) %in% strip_ensembl_version(pathway_genes)
  tab$in_pathway <- in_pw
  tab <- tab[!in_pw, , drop = FALSE]
  tab$fdr <- stats::p.adjust(tab$p_value, method = "BH")
  if (!is.null(symbol_map)) {
    key <- strip_ensembl_version(tab$gene)
    tab$symbol <- symbol_map[key]
    if (all(is.na(tab$symbol))) tab$symbol <- symbol_map[tab$gene]
  }
  keep <- is.finite(tab$rho) & tab$rho < 0 & tab$fdr <= fdr_cutoff & tab$rho <= min_rho
  out <- tab[keep, , drop = FALSE]
  out[order(out$rho), ]
}

negative_proteins_one_go <- function(score, protein, fdr_cutoff, min_rho) {
  if (is.null(protein) || !ncol(protein) || !nrow(protein)) return(NULL)
  al_names <- intersect(names(score), colnames(protein))
  if (length(al_names) < 10L) {
    # 蛋白矩阵可能是 samples x proteins
    return(NULL)
  }
  tab <- spearman_vs_matrix(score, protein)
  tab$fdr <- stats::p.adjust(tab$p_value, method = "BH")
  keep <- is.finite(tab$rho) & tab$rho < 0 & tab$fdr <= fdr_cutoff & tab$rho <= min_rho
  out <- tab[keep, , drop = FALSE]
  names(out)[names(out) == "gene"] <- "protein"
  out[order(out$rho), ]
}

write_tsv <- function(df, path) {
  ensure_dir(dirname(path))
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

go_file_stub <- function(go_id) gsub(":", "_", go_id)

annotate_go <- function(df, gid, gname) {
  n <- nrow(df)
  df$go_id <- rep(gid, n)
  df$go_name <- rep(gname, n)
  df
}

# -----------------------------------------------------------------------------
# Demo：不读 TCGA，用合成数据验证“每个 GO 独立”
# -----------------------------------------------------------------------------

run_demo <- function(out_dir, fdr, min_rho, min_genes) {
  log_msg("DEMO 模式：合成数据，只演示每个 GO 独立循环，结果不可用于发表")
  set.seed(1)
  n_sample <- 40L
  patients <- sprintf("TCGA-A1-%04d", seq_len(n_sample))
  samples <- paste0(patients, "-01A")
  genes <- paste0("GENE", sprintf("%03d", seq_len(60)))
  expr <- matrix(rnorm(length(genes) * n_sample), nrow = length(genes), dimnames = list(genes, samples))
  go_ids <- GO_IDS[1:3]
  genesets <- list(genes[1:8], genes[9:16], genes[17:24])
  names(genesets) <- go_ids
  # 让 GO1 得分与 GENE050 强负相关，且不影响 GO2
  expr[genes[1:8], ] <- expr[genes[1:8], ] + rep(seq(-2, 2, length.out = n_sample), each = 8)
  expr["GENE050", ] <- -colMeans(expr[genes[1:8], , drop = FALSE]) + rnorm(n_sample, sd = 0.05)
  clinical <- data.frame(
    age = 40 + seq_len(n_sample) + rnorm(n_sample),
    stage = rep(c("I", "II", "III"), length.out = n_sample),
    row.names = samples
  )
  survival <- data.frame(
    OS = rbinom(n_sample, 1, 0.4),
    OS.time = sample(200:2000, n_sample, replace = TRUE),
    row.names = samples
  )
  expr <- log2(pmax(expr - min(expr) + 0.1, 0.1))
  demo_out <- file.path(out_dir, "demo")
  run_each_go(
    go_ids = go_ids,
    genesets_on_expr = genesets,
    expr = expr,
    clinical = clinical,
    survival = survival,
    protein = NULL,
    symbol_map = stats::setNames(genes, genes),
    out_dir = demo_out,
    fdr = fdr,
    min_rho = min_rho,
    min_genes = min_genes
  )
  neg1 <- utils::read.delim(
    file.path(demo_out, "neg_genes", paste0(go_file_stub(go_ids[[1]]), "_negative_genes.tsv")),
    check.names = FALSE
  )
  if (!"GENE050" %in% neg1$gene) {
    stop("DEMO 失败：GO1 的负相关基因中应包含 GENE050（每个通路独立计算）")
  }
  log_msg("DEMO 检查通过：", go_ids[[1]], " 独立识别到负相关基因 GENE050")
}

# -----------------------------------------------------------------------------
# 主循环：一个 GO 一次，写一个 GO 的文件
# -----------------------------------------------------------------------------

run_each_go <- function(go_ids, genesets_on_expr, expr, clinical, survival, protein,
                        symbol_map, out_dir, fdr, min_rho, min_genes) {
  dir_clin <- file.path(out_dir, "clinical")
  dir_surv <- file.path(out_dir, "survival")
  dir_neg <- file.path(out_dir, "neg_genes")
  dir_prot <- file.path(out_dir, "neg_proteins")
  dir_score <- file.path(out_dir, "scores")
  lapply(c(dir_clin, dir_surv, dir_neg, dir_prot, dir_score), ensure_dir)

  skip_log <- list()
  clin_all <- list()
  surv_all <- list()
  neg_counts <- list()

  for (gid in go_ids) {
    gname <- go_term_name(gid)
    genes_i <- genesets_on_expr[[gid]]
    n_i <- length(genes_i)
    log_msg("==== ", gid, " | ", gname, " | 匹配到 ", n_i, " 个表达基因 ====")
    stub <- go_file_stub(gid)

    if (n_i < min_genes) {
      msg <- sprintf("基因数 %s < %s，跳过（不并入其他通路）", n_i, min_genes)
      log_msg(msg)
      skip_log[[gid]] <- data.frame(go_id = gid, go_name = gname, n_genes = n_i, reason = msg)
      next
    }

    score_i <- one_go_pathway_score(expr, genes_i, gid)
    if (is.null(score_i)) {
      msg <- "无法计算通路得分，跳过"
      log_msg(msg)
      skip_log[[gid]] <- data.frame(go_id = gid, go_name = gname, n_genes = n_i, reason = msg)
      next
    }
    if (is.null(names(score_i))) names(score_i) <- colnames(expr)

    score_df <- data.frame(
      sample = names(score_i),
      patient = patient_barcode(names(score_i)),
      score = as.numeric(score_i),
      score_method = "mean_z",
      go_id = gid,
      go_name = gname,
      stringsAsFactors = FALSE
    )
    write_tsv(score_df, file.path(dir_score, paste0(stub, "_score.tsv")))

    clin_i <- annotate_go(clinical_association_one_go(score_i, clinical), gid, gname)
    write_tsv(clin_i, file.path(dir_clin, paste0(stub, "_clinical.tsv")))
    clin_all[[gid]] <- clin_i

    surv_i <- annotate_go(survival_association_one_go(score_i, survival), gid, gname)
    write_tsv(surv_i, file.path(dir_surv, paste0(stub, "_survival.tsv")))
    surv_all[[gid]] <- surv_i

    neg_i <- annotate_go(
      negative_genes_one_go(score_i, expr, genes_i, fdr, min_rho, symbol_map),
      gid, gname
    )
    write_tsv(neg_i, file.path(dir_neg, paste0(stub, "_negative_genes.tsv")))
    neg_counts[[gid]] <- data.frame(
      go_id = gid, go_name = gname, n_genes_in_pathway = n_i,
      n_negative_genes = nrow(neg_i),
      stringsAsFactors = FALSE
    )
    log_msg(gid, " 负相关基因: ", nrow(neg_i))

    if (!is.null(protein)) {
        prot_i <- negative_proteins_one_go(score_i, protein, fdr, min_rho)
      if (!is.null(prot_i)) {
        prot_i <- annotate_go(prot_i, gid, gname)
        write_tsv(prot_i, file.path(dir_prot, paste0(stub, "_negative_proteins.tsv")))
      }
    }
  }

  if (length(clin_all)) {
    write_tsv(do.call(rbind, clin_all), file.path(out_dir, "summary_clinical_all_go.tsv"))
  }
  if (length(surv_all)) {
    write_tsv(do.call(rbind, surv_all), file.path(out_dir, "summary_survival_all_go.tsv"))
  }
  if (length(neg_counts)) {
    write_tsv(do.call(rbind, neg_counts), file.path(out_dir, "summary_neg_gene_counts.tsv"))
  }
  if (length(skip_log)) {
    write_tsv(do.call(rbind, skip_log), file.path(out_dir, "skipped_go.tsv"))
  }
  log_msg("完成。每个 GO 的结果在: ", normalizePath(out_dir, mustWork = FALSE))
}

prepare_symbol_map <- function(gencode_df, expr_ids) {
  map <- stats::setNames(rep(NA_character_, length(expr_ids)), strip_ensembl_version(expr_ids))
  if (!is.null(gencode_df) && nrow(gencode_df)) {
    map[strip_ensembl_version(gencode_df$gene_id)] <- gencode_df$gene_name
    if ("gene_name" %in% names(gencode_df)) {
      map[gencode_df$gene_name] <- gencode_df$gene_name
    }
  }
  map
}

harmonize_expr_rownames <- function(expr, gencode_df) {
  rn <- rownames(expr)
  rownames(expr) <- strip_ensembl_version(rn)
  if (!is.null(gencode_df) && mean(grepl("^ENSG", rownames(expr))) < 0.5) {
    return(expr)
  }
  expr
}

maybe_log_fpkm <- function(expr) {
  mx <- max(expr, na.rm = TRUE)
  if (is.finite(mx) && mx > 100) {
    log_msg("FPKM 看起来未取对数，使用 log2(FPKM + 1)")
    expr <- log2(expr + 1)
  } else {
    log_msg("表达值范围较小，假定已是 log 空间，不再 log2")
  }
  expr
}

protein_as_gene_matrix <- function(protein_df) {
  nms <- names(protein_df)
  rn <- rownames(protein_df)
  row_tcga <- mean(is_tcga_id(rn), na.rm = TRUE)
  col_tcga <- mean(is_tcga_id(nms), na.rm = TRUE)
  if (!is.finite(row_tcga)) row_tcga <- 0
  if (!is.finite(col_tcga)) col_tcga <- 0
  num <- protein_df
  num[] <- lapply(num, function(x) suppressWarnings(as.numeric(as.character(x))))
  if (col_tcga >= 0.5 || col_tcga >= row_tcga) {
    # 蛋白 x 样本
    mat <- as.matrix(num)
    rownames(mat) <- rn
  } else {
    # 样本 x 蛋白
    mat <- t(as.matrix(num))
    colnames(mat) <- rn
  }
  storage.mode(mat) <- "double"
  mat
}

run_real <- function(opt) {
  paths <- lapply(DATA_STEMS, find_data_file, data_dir = opt$data_dir)
  required <- c("clinical", "expr", "survival", "gencode")
  missing <- required[is.na(unlist(paths[required]))]
  if (length(missing)) {
    stop(
      "在 ", opt$data_dir, " 中找不到: ",
      paste(DATA_STEMS[missing], collapse = ", "),
      "\n请把数据文件放到该目录，或改用 --demo 先跑通流程。"
    )
  }
  expr <- read_expr_matrix(paths$expr)
  if (opt$primary_tumor_only) {
    before <- ncol(expr)
    expr <- keep_primary_tumor(expr)
    log_msg("保留原发瘤且每患者一个样本: ", before, " -> ", ncol(expr))
  }
  gencode_df <- read_gencode(paths$gencode)
  expr <- harmonize_expr_rownames(expr, gencode_df)
  expr <- maybe_log_fpkm(expr)
  clinical <- read_table_with_ids(paths$clinical)
  survival <- read_table_with_ids(paths$survival)
  protein <- NULL
  if (!is.na(paths$protein)) {
    protein <- protein_as_gene_matrix(read_table_with_ids(paths$protein))
  } else {
    log_msg("未找到蛋白文件，跳过蛋白负相关")
  }

  log_msg("按每个 GO 分别从 QuickGO 取基因集（本 term + descendants；不使用 AnnotationDbi）")
  cache_dir <- file.path(opt$data_dir, "cache", "go_genes")
  ensure_dir(cache_dir)
  genesets_on_expr <- list()
  gene_inventory <- list()
  expr_ids <- rownames(expr)
  symbol_map <- prepare_symbol_map(gencode_df, expr_ids)

  for (gid in GO_IDS) {
    raw <- map_symbols_via_gencode(get_genes_for_one_go(gid, cache_dir), gencode_df)
    hit_ens <- match_genes_to_expr(raw$ensembl, expr_ids)
    hit_sym <- match_genes_to_expr(raw$symbol, expr_ids)
    hit <- unique(c(hit_ens, hit_sym))
    genesets_on_expr[[gid]] <- hit
    gene_inventory[[gid]] <- data.frame(
      go_id = gid,
      go_name = go_term_name(gid),
      n_annotated_ensembl = length(raw$ensembl),
      n_annotated_symbol = length(raw$symbol),
      n_in_expr = length(hit),
      stringsAsFactors = FALSE
    )
    log_msg(gid, " 注释基因 ", length(raw$ensembl), "/", length(raw$symbol),
            "，表达矩阵命中 ", length(hit))
  }
  ensure_dir(opt$out_dir)
  write_tsv(do.call(rbind, gene_inventory), file.path(opt$out_dir, "go_gene_inventory.tsv"))

  run_each_go(
    go_ids = GO_IDS,
    genesets_on_expr = genesets_on_expr,
    expr = expr,
    clinical = clinical,
    survival = survival,
    protein = protein,
    symbol_map = symbol_map,
    out_dir = opt$out_dir,
    fdr = opt$fdr,
    min_rho = opt$min_rho,
    min_genes = opt$min_genes
  )
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  opt <- parse_args(argv)
  if (opt$install_deps) install_deps_if_needed()
  if (opt$demo) {
    run_demo(opt$out_dir, opt$fdr, opt$min_rho, opt$min_genes)
  } else {
    run_real(opt)
  }
  invisible(NULL)
}

if (sys.nframe() == 0L) {
  main()
}
