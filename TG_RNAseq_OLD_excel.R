#!/usr/bin/env Rscript
# =============================================================================
# OLD盘 Excel 表：NTC vs TG_sh1 差异分析
#
# 脚本位置（仓库根目录，和 TG_RNAseq_pipeline.R 放一起）：
#   TG_RNAseq_OLD_excel.R
#
# 用法（在 R / RStudio 里）：
#   setwd("E:/R/TG_BRCA/TG")
#   source("TG_RNAseq_OLD_excel.R")
#
# Excel 放在：
#   E:/R/TG_BRCA/TG/OLD盘/*.xlsx  （或 .xls）
# 表头需含截图中的列：
#   test_id, gene_id, gene, locus, sample_1, sample_2, status,
#   value_1, value_2, log2(fold_change), test_stat, p_value, q_value, significant
# sample_1=NTC，sample_2=TG_sh1。只保留 status==OK。
# 三套思路：上调 FC；上调 topN；p<0.01 与 p<0.05（不区分上下调，表和图都做标记）
#
# 结果：
#   E:/R/TG_BRCA/TG/OLD盘/results_TG_sh1_vs_NTC/
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")
options(clusterProfiler.download.method = "auto")

# -----------------------------------------------------------------------------
# 0. 依赖
# -----------------------------------------------------------------------------
cran_required <- c(
  "readxl", "writexl", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer", "stringr"
)
bioc_required <- c(
  "clusterProfiler", "org.Hs.eg.db", "enrichplot", "AnnotationDbi",
  "fgsea", "msigdbr", "GSVA"
)
bioc_optional <- c(
  "ReactomePA", "pathview",
  "GenomicRanges", "GenomicFeatures", "IRanges", "S4Vectors", "GenomeInfoDb",
  "TxDb.Hsapiens.UCSC.hg38.knownGene", "TxDb.Hsapiens.UCSC.hg19.knownGene"
)

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
  if (length(still) > 0) message("可选包未安装，相关分析将跳过: ", paste(still, collapse = ", "))
  invisible(TRUE)
}

install_if_missing(cran_required, bioc = FALSE, required = TRUE)
install_if_missing(bioc_required, bioc = TRUE, required = TRUE)
install_if_missing(bioc_optional, bioc = TRUE, required = FALSE)

safe_library <- function(pkgs) {
  for (p in pkgs) {
    if (requireNamespace(p, quietly = TRUE)) {
      suppressPackageStartupMessages(library(p, character.only = TRUE))
    }
  }
}
safe_library(c(cran_required, bioc_required, bioc_optional))
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

# -----------------------------------------------------------------------------
# 1. 路径
# -----------------------------------------------------------------------------
resolve_old_dir <- function() {
  env_dir <- Sys.getenv("TG_RNASEQ_OLD_DIR", unset = "")
  candidates <- unique(c(
    env_dir,
    "E:/R/TG_BRCA/TG/OLD盘",
    "E:\\R\\TG_BRCA\\TG\\OLD盘",
    file.path("E:/R/TG_BRCA/TG", "OLD"),
    file.path(getwd(), "OLD盘"),
    file.path(getwd(), "OLD"),
    getwd()
  ))
  candidates <- candidates[nzchar(candidates)]
  has_excel <- function(d) {
    dir.exists(d) && length(list.files(d, pattern = "\\.(xlsx|xls)$", ignore.case = TRUE)) > 0
  }
  for (d in candidates) {
    if (has_excel(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  for (d in candidates) {
    if (dir.exists(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

old_dir <- resolve_old_dir()
result_dir <- file.path(old_dir, "results_TG_sh1_vs_NTC")
log_dir <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("OLD_excel_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_msg("OLD盘目录: ", old_dir)
log_msg("结果目录: ", result_dir)

fc_cutoffs <- c("FC_1" = 1, "FC_1.25" = 1.25, "FC_1.5" = 1.5, "FC_2" = 2)
top_ns <- c(50, 75, 100, 150, 200, 250, 300, 4000)
p_cutoffs <- c("p_lt_0.01" = 0.01, "p_lt_0.05" = 0.05)

# -----------------------------------------------------------------------------
# 2. 读 Excel（截图表头）
# -----------------------------------------------------------------------------
norm_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

find_col <- function(df, ...) {
  keys <- unique(norm_key(unlist(list(...), use.names = FALSE)))
  n <- norm_key(names(df))
  for (k in keys) {
    hit <- which(n == k)
    if (length(hit) > 0) return(names(df)[hit[1]])
  }
  for (k in keys) {
    hit <- which(grepl(k, n, fixed = TRUE))
    if (length(hit) > 0) return(names(df)[hit[1]])
  }
  NA_character_
}

pick_official_symbol <- function(x) {
  x <- trimws(as.character(x))
  if (length(x) != 1 || is.na(x) || x %in% c("", "-", ".", "NA")) return(NA_character_)
  parts <- unlist(strsplit(x, "[,;|/]+"))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts) & !parts %in% c("-", ".", "NA")]
  if (length(parts) == 0) return(NA_character_)
  is_fusion <- vapply(parts, function(t) {
    bits <- strsplit(t, "-", fixed = TRUE)[[1]]
    length(bits) == 2 && (bits[1] %in% parts || bits[2] %in% parts)
  }, logical(1))
  if (any(!is_fusion)) parts <- parts[!is_fusion]
  score <- vapply(parts, function(s) {
    if (grepl("^(XLOC|TCONS|CUFF)_", s, ignore.case = TRUE)) return(0)
    if (grepl("^MIR[0-9]", s, ignore.case = TRUE)) return(1)
    if (grepl("^LOC[0-9]+$", s, ignore.case = TRUE)) return(2)
    3
  }, numeric(1))
  parts[which.max(score)]
}

find_excel_file <- function(dir) {
  files <- list.files(dir, pattern = "\\.(xlsx|xls)$", full.names = TRUE, ignore.case = TRUE)
  files <- files[!grepl("^~\\$", basename(files))]
  files <- files[!grepl("^results_", basename(files), ignore.case = TRUE)]
  if (length(files) == 0) return(NA_character_)
  score <- vapply(files, function(f) {
    b <- tolower(basename(f))
    s <- 0
    if (grepl("gene_exp|diff|cuff|ntc|tg_sh|express|表达|差异", b)) s <- s + 10
    if (grepl("sh1", b)) s <- s + 3
    s
  }, numeric(1))
  files[order(-score, nchar(basename(files)))][1]
}

sheet_has_de_cols <- function(df) {
  n <- norm_key(names(df))
  has_fc <- any(grepl("log2", n) & grepl("fold", n)) || "log2_fold_change" %in% n
  has_p <- "p_value" %in% n || "pvalue" %in% n
  has_gene <- any(n %in% c("gene", "gene_id", "gene_short_name", "symbol"))
  has_fc && has_p && has_gene
}

read_old_excel <- function(path) {
  sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) "Sheet1")
  picked <- NULL
  for (sh in sheets) {
    df <- tryCatch(
      as.data.frame(readxl::read_excel(path, sheet = sh), stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(df) || ncol(df) < 5) next
    if (sheet_has_de_cols(df)) {
      picked <- df
      log_msg("使用工作表: ", sh)
      break
    }
    if (is.null(picked)) picked <- df
  }
  if (is.null(picked)) stop("无法读取 Excel: ", path)
  picked
}

excel_path <- find_excel_file(old_dir)
if (is.na(excel_path)) {
  stop(
    "在目录里找不到 .xlsx/.xls：", old_dir,
    "\n请把截图那份 Excel 放到 E:/R/TG_BRCA/TG/OLD盘/ 后再运行 TG_RNAseq_OLD_excel.R"
  )
}
log_msg("Excel 文件: ", excel_path)

raw <- read_old_excel(excel_path)
log_msg("原始行列: ", nrow(raw), " x ", ncol(raw), " ; 列名: ", paste(names(raw), collapse = ", "))

col_gene <- find_col(raw, "gene", "gene_short_name", "symbol", "gene_name")
col_gene_id <- find_col(raw, "gene_id", "test_id", "tracking_id")
col_status <- find_col(raw, "status")
col_v1 <- find_col(raw, "value_1", "value1", "fpkm_1")
col_v2 <- find_col(raw, "value_2", "value2", "fpkm_2")
col_lfc <- find_col(raw, "log2_fold_change", "log2(fold_change)", "log2foldchange", "log2fc")
col_p <- find_col(raw, "p_value", "pvalue", "pval")
col_q <- find_col(raw, "q_value", "qvalue", "padj", "fdr")
col_s1 <- find_col(raw, "sample_1", "sample1")
col_s2 <- find_col(raw, "sample_2", "sample2")
col_stat <- find_col(raw, "test_stat")
col_sig <- find_col(raw, "significant")
col_locus <- find_col(raw, "locus")
col_test <- find_col(raw, "test_id")

need <- c(col_gene = col_gene, col_lfc = col_lfc, col_p = col_p)
if (anyNA(need)) {
  stop("Excel 缺少必要列。已有: ", paste(names(raw), collapse = ", "),
       "；缺失: ", paste(names(need)[is.na(need)], collapse = ", "))
}

de <- data.frame(
  test_id = if (!is.na(col_test)) as.character(raw[[col_test]]) else as.character(seq_len(nrow(raw))),
  gene_id = if (!is.na(col_gene_id)) as.character(raw[[col_gene_id]]) else NA_character_,
  gene_raw = as.character(raw[[col_gene]]),
  locus = if (!is.na(col_locus)) as.character(raw[[col_locus]]) else NA_character_,
  sample_1 = if (!is.na(col_s1)) as.character(raw[[col_s1]]) else "NTC",
  sample_2 = if (!is.na(col_s2)) as.character(raw[[col_s2]]) else "TG_sh1",
  status = if (!is.na(col_status)) as.character(raw[[col_status]]) else "OK",
  value_1 = if (!is.na(col_v1)) suppressWarnings(as.numeric(raw[[col_v1]])) else NA_real_,
  value_2 = if (!is.na(col_v2)) suppressWarnings(as.numeric(raw[[col_v2]])) else NA_real_,
  log2FC = suppressWarnings(as.numeric(raw[[col_lfc]])),
  test_stat = if (!is.na(col_stat)) suppressWarnings(as.numeric(raw[[col_stat]])) else NA_real_,
  pvalue = suppressWarnings(as.numeric(raw[[col_p]])),
  q_value = if (!is.na(col_q)) suppressWarnings(as.numeric(raw[[col_q]])) else NA_real_,
  significant = if (!is.na(col_sig)) as.character(raw[[col_sig]]) else NA_character_,
  stringsAsFactors = FALSE
)

n_before <- nrow(de)
de <- de[!is.na(de$status) & toupper(de$status) == "OK", , drop = FALSE]
log_msg("status==OK: ", nrow(de), " / ", n_before, " (丢掉 NOTEST 等 ", n_before - nrow(de), " 行)")

de$gene <- vapply(de$gene_raw, pick_official_symbol, character(1), USE.NAMES = FALSE)
missing_sym <- is.na(de$gene) | !nzchar(de$gene)
de$gene[missing_sym] <- ifelse(
  !is.na(de$test_id[missing_sym]) & nzchar(de$test_id[missing_sym]),
  de$test_id[missing_sym],
  de$gene_id[missing_sym]
)
de$FC <- 2^de$log2FC
de$AveExpr <- rowMeans(cbind(de$value_1, de$value_2), na.rm = TRUE)

finite <- is.finite(de$log2FC)
if (!all(finite)) {
  cap <- suppressWarnings(max(abs(de$log2FC[finite]), na.rm = TRUE))
  if (!is.finite(cap) || cap < 10) cap <- 10
  de$log2FC[!finite & de$log2FC > 0] <- cap + 1
  de$log2FC[!finite & de$log2FC < 0] <- -(cap + 1)
  de$FC <- 2^de$log2FC
  log_msg("非有限 log2FC 已截断到 ±", cap + 1, " ，条数=", sum(!finite))
}

de <- de[!is.na(de$gene) & nzchar(de$gene) & !is.na(de$log2FC), , drop = FALSE]
de <- de[order(de$pvalue, -abs(de$log2FC), na.last = TRUE), , drop = FALSE]
de <- de[!duplicated(de$gene), , drop = FALSE]
log_msg("去重后基因数: ", nrow(de))

s1 <- unique(de$sample_1)
s1 <- s1[!is.na(s1) & nzchar(s1)][1]
s2 <- unique(de$sample_2)
s2 <- s2[!is.na(s2) & nzchar(s2)][1]
if (is.na(s1) || !nzchar(s1)) s1 <- "NTC"
if (is.na(s2) || !nzchar(s2)) s2 <- "TG_sh1"
comp_name <- paste0(s2, "_vs_", s1)
log_msg("比较: ", s2, " vs ", s1, "  (目录名 ", comp_name, ")")

n_xloc_gene <- sum(grepl("^(XLOC|TCONS|CUFF)_", de$gene, ignore.case = TRUE))
log_msg(
  "gene 列为 XLOC/TCONS/CUFF 的行: ", n_xloc_gene,
  " / ", nrow(de),
  "。这些是 Cufflinks 组装位点 ID，后面的数字是不同基因组区间，不是同一基因的不同转录本。"
)

# -----------------------------------------------------------------------------
# 3. ID 映射与绘图（全表只映射一次；bitr 未映射警告不再每个子集刷屏）
# -----------------------------------------------------------------------------
bitr_quiet <- function(keys, fromType) {
  keys <- unique(as.character(keys))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  if (length(keys) == 0) return(NULL)
  suppressWarnings(suppressMessages(
    tryCatch(
      clusterProfiler::bitr(keys, fromType = fromType, toType = "ENTREZID", OrgDb = org.Hs.eg.db),
      error = function(e) NULL
    )
  ))
}

append_map <- function(acc, src, query_col) {
  if (is.null(src) || nrow(src) == 0) return(acc)
  df <- data.frame(
    gene = as.character(src[[query_col]]),
    entrez = as.character(src$ENTREZID),
    via = query_col,
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$gene) & nzchar(df$gene) & !is.na(df$entrez) & nzchar(df$entrez), , drop = FALSE]
  if (nrow(df) == 0) return(acc)
  rbind(acc, df)
}

build_gene_entrez_map <- function(symbols, extra_alias = NULL) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  acc <- data.frame(gene = character(), entrez = character(), via = character(), stringsAsFactors = FALSE)

  acc <- append_map(acc, bitr_quiet(symbols, "SYMBOL"), "SYMBOL")
  mapped <- unique(acc$gene)
  rest <- setdiff(symbols, mapped)

  if (length(rest) > 0) {
    acc <- append_map(acc, bitr_quiet(rest, "ALIAS"), "ALIAS")
    mapped <- unique(acc$gene)
    rest <- setdiff(symbols, mapped)
  }

  ens <- rest[grepl("^ENS[GT]", rest, ignore.case = TRUE)]
  if (length(ens) > 0) {
    acc <- append_map(acc, bitr_quiet(ens, "ENSEMBL"), "ENSEMBL")
    mapped <- unique(acc$gene)
    rest <- setdiff(symbols, mapped)
  }

  if (length(rest) > 0) {
    all_sym <- tryCatch(AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL"), error = function(e) character())
    if (length(all_sym) > 0) {
      hit <- match(toupper(rest), toupper(all_sym))
      ok <- !is.na(hit)
      if (any(ok)) {
        canon <- all_sym[hit[ok]]
        m_case <- bitr_quiet(canon, "SYMBOL")
        if (!is.null(m_case) && nrow(m_case) > 0) {
          m_case$QUERY <- rest[ok][match(m_case$SYMBOL, canon)]
          acc <- append_map(acc, m_case, "QUERY")
        }
        mapped <- unique(acc$gene)
        rest <- setdiff(symbols, mapped)
      }
    }
  }

  extra_alias <- unique(as.character(extra_alias))
  extra_alias <- extra_alias[!is.na(extra_alias) & nzchar(extra_alias)]
  if (length(extra_alias) > 0) {
    acc <- append_map(acc, bitr_quiet(extra_alias, "ALIAS"), "ALIAS")
  }

  acc <- acc[!duplicated(acc$gene), , drop = FALSE]
  unmapped <- setdiff(symbols, acc$gene)
  list(map = acc[, c("gene", "entrez")], unmapped = unmapped, via = acc)
}

idmap <- build_gene_entrez_map(de$gene, extra_alias = NULL)
gene_entrez_map <- idmap$map
still <- setdiff(unique(de$gene), gene_entrez_map$gene)
if (length(still) > 0) {
  raws <- as.character(de$gene_raw[match(still, de$gene)])
  parts <- strsplit(raws, "[,;|/]+")
  long <- data.frame(
    gene = rep(still, lengths(parts)),
    part = trimws(unlist(parts)),
    stringsAsFactors = FALSE
  )
  long <- long[!is.na(long$part) & nzchar(long$part), , drop = FALSE]
  m_sym <- bitr_quiet(unique(long$part), "SYMBOL")
  m_alias <- bitr_quiet(unique(long$part), "ALIAS")
  conv <- data.frame(part = character(), ENTREZID = character(), stringsAsFactors = FALSE)
  if (!is.null(m_sym) && nrow(m_sym) > 0) {
    conv <- rbind(conv, data.frame(part = as.character(m_sym$SYMBOL), ENTREZID = as.character(m_sym$ENTREZID)))
  }
  if (!is.null(m_alias) && nrow(m_alias) > 0) {
    conv <- rbind(conv, data.frame(part = as.character(m_alias$ALIAS), ENTREZID = as.character(m_alias$ENTREZID)))
  }
  if (nrow(conv) > 0) {
    long$entrez <- conv$ENTREZID[match(long$part, conv$part)]
    long <- long[!is.na(long$entrez), , drop = FALSE]
    long <- long[!duplicated(long$gene), , drop = FALSE]
    gene_entrez_map <- rbind(
      gene_entrez_map,
      data.frame(gene = long$gene, entrez = as.character(long$entrez), stringsAsFactors = FALSE)
    )
    gene_entrez_map <- gene_entrez_map[!duplicated(gene_entrez_map$gene), , drop = FALSE]
  }
}

is_cuff_id <- function(x) {
  grepl("^(XLOC|TCONS|CUFF)_", as.character(x), ignore.case = TRUE)
}

parse_locus_df <- function(locus) {
  locus <- as.character(locus)
  m <- stringr::str_match(locus, "(?i)^(chr)?([^:]+):([0-9]+)-([0-9]+)")
  chr_raw <- m[, 3]
  chr <- ifelse(is.na(chr_raw), NA_character_, paste0("chr", sub("(?i)^chr", "", chr_raw, perl = TRUE)))
  data.frame(
    chr = chr,
    start = suppressWarnings(as.integer(m[, 4])),
    end = suppressWarnings(as.integer(m[, 5])),
    stringsAsFactors = FALSE
  )
}

load_txdb_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NULL)
  getExportedValue(pkg, pkg)
}

overlap_entrez <- function(chr, start, end, txdb) {
  n <- length(chr)
  ans <- rep(NA_character_, n)
  ok <- !is.na(chr) & !is.na(start) & !is.na(end) & end >= start
  if (!any(ok) || is.null(txdb)) return(ans)
  gr <- GenomicRanges::GRanges(
    seqnames = chr[ok],
    ranges = IRanges::IRanges(pmin(start[ok], end[ok]), pmax(start[ok], end[ok]))
  )
  g <- GenomicFeatures::genes(txdb)
  if (has_pkg("GenomeInfoDb")) {
    gr <- tryCatch(GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse"), error = function(e) gr)
    g <- tryCatch(GenomeInfoDb::keepStandardChromosomes(g, pruning.mode = "coarse"), error = function(e) g)
    suppressWarnings(try(GenomeInfoDb::seqlevelsStyle(gr) <- GenomeInfoDb::seqlevelsStyle(g)[1], silent = TRUE))
  }
  ol <- GenomicRanges::findOverlaps(gr, g, ignore.strand = TRUE)
  if (length(ol) == 0) return(ans)
  qh <- S4Vectors::queryHits(ol)
  sh <- S4Vectors::subjectHits(ol)
  ov <- GenomicRanges::pintersect(gr[qh], g[sh], ignore.strand = TRUE)
  w <- GenomicRanges::width(ov)
  entrez <- as.character(GenomicRanges::mcols(g)$gene_id)
  if (length(entrez) != length(g) || all(is.na(entrez) | !nzchar(entrez))) {
    entrez <- as.character(names(g))
  }
  ord <- order(qh, -w)
  keep <- !duplicated(qh[ord])
  q_ok <- which(ok)
  ans[q_ok[qh[ord][keep]]] <- entrez[sh[ord][keep]]
  ans
}

entrez_to_symbol <- function(entrez) {
  entrez <- as.character(entrez)
  out <- rep(NA_character_, length(entrez))
  ok <- !is.na(entrez) & nzchar(entrez)
  if (!any(ok)) return(out)
  sym <- tryCatch(
    AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = unique(entrez[ok]), column = "SYMBOL",
      keytype = "ENTREZID", multiVals = "first"
    ),
    error = function(e) NULL
  )
  if (is.null(sym)) return(out)
  out[ok] <- unname(sym[entrez[ok]])
  out
}

de$xloc_rescued <- FALSE
de$map_via <- ifelse(de$gene %in% gene_entrez_map$gene, "SYMBOL/ALIAS", NA_character_)
need_idx <- which(
  (is_cuff_id(de$gene) | is_cuff_id(de$gene_raw) | is_cuff_id(de$test_id)) &
    !(de$gene %in% gene_entrez_map$gene)
)
if (length(need_idx) == 0) {
  need_idx <- which(!(de$gene %in% gene_entrez_map$gene) & !is.na(de$locus) & nzchar(de$locus))
}

can_overlap <- has_pkg("GenomicRanges") && has_pkg("GenomicFeatures") && has_pkg("IRanges")
if (length(need_idx) > 0 && !can_overlap) {
  log_msg("XLOC 有 ", length(need_idx), " 个未映射，但 GenomicRanges/TxDb 未安装，无法用 locus 补注释")
}
if (length(need_idx) > 0 && can_overlap) {
  loc <- parse_locus_df(de$locus[need_idx])
  tx38 <- load_txdb_pkg("TxDb.Hsapiens.UCSC.hg38.knownGene")
  tx19 <- load_txdb_pkg("TxDb.Hsapiens.UCSC.hg19.knownGene")
  rate38 <- if (is.null(tx38)) -1 else mean(!is.na(overlap_entrez(loc$chr, loc$start, loc$end, tx38)))
  rate19 <- if (is.null(tx19)) -1 else mean(!is.na(overlap_entrez(loc$chr, loc$start, loc$end, tx19)))
  if (rate38 < 0 && rate19 < 0) {
    log_msg("未安装 TxDb hg38/hg19，XLOC 无法用坐标补基因名")
  } else {
    use_hg19 <- rate19 > rate38
    txdb <- if (use_hg19) tx19 else tx38
    build <- if (use_hg19) "hg19" else "hg38"
    log_msg(sprintf(
      "XLOC 用 locus 重叠已知基因：hg38 overlap=%.1f%%, hg19 overlap=%.1f%%，选用 %s",
      100 * max(rate38, 0), 100 * max(rate19, 0), build
    ))
    hit_entrez <- overlap_entrez(loc$chr, loc$start, loc$end, txdb)
    hit_sym <- entrez_to_symbol(hit_entrez)
    ok <- !is.na(hit_entrez) & nzchar(hit_entrez)
    n_rescue <- sum(ok)
    log_msg("XLOC/未映射位点经 locus 补到 Entrez: ", n_rescue, " / ", length(need_idx))
    if (n_rescue > 0) {
      idx <- need_idx[ok]
      de$xloc_rescued[idx] <- TRUE
      de$map_via[idx] <- paste0("locus_", build)
      new_sym <- hit_sym[ok]
      new_ent <- as.character(hit_entrez[ok])
      fallback <- is.na(new_sym) | !nzchar(new_sym)
      new_sym[fallback] <- paste0("ENTREZ_", new_ent[fallback])
      de$gene[idx] <- new_sym
      add <- data.frame(gene = new_sym, entrez = new_ent, stringsAsFactors = FALSE)
      gene_entrez_map <- rbind(gene_entrez_map, add)
      gene_entrez_map <- gene_entrez_map[!duplicated(gene_entrez_map$gene), , drop = FALSE]
      loc_report <- data.frame(
        test_id = de$test_id[idx],
        gene_raw = de$gene_raw[idx],
        locus = de$locus[idx],
        genome = build,
        entrez = new_ent,
        gene_symbol = new_sym,
        stringsAsFactors = FALSE
      )
      utils::write.csv(loc_report, file.path(log_dir, "XLOC_rescued_by_locus.csv"), row.names = FALSE)
    }
  }
}

de <- de[order(de$pvalue, -abs(de$log2FC), na.last = TRUE), , drop = FALSE]
de <- de[!duplicated(de$gene), , drop = FALSE]
gene_entrez_map <- gene_entrez_map[gene_entrez_map$gene %in% de$gene, , drop = FALSE]
log_msg("locus 补注释并按官方符号去重后基因数: ", nrow(de))

heat_mat <- cbind(de$value_1, de$value_2)
colnames(heat_mat) <- c(s1, s2)
rownames(heat_mat) <- de$gene
heat_mat[is.na(heat_mat)] <- 0
if (max(heat_mat, na.rm = TRUE) > 50) {
  heat_log <- log2(heat_mat + 1)
} else {
  heat_log <- heat_mat
}
sample_info <- data.frame(
  sample = c(s1, s2),
  group = c("NTC", "TG_sh1"),
  stringsAsFactors = FALSE
)

n_input <- length(unique(de$gene))
n_map <- sum(de$gene %in% gene_entrez_map$gene)
n_fail <- n_input - n_map
n_xloc_left <- sum(is_cuff_id(de$gene))
pct_fail <- if (n_input > 0) 100 * n_fail / n_input else 0
idmap$unmapped <- setdiff(unique(de$gene), gene_entrez_map$gene)
log_msg(sprintf(
  "Entrez 映射: %d / %d 成功 (%.1f%%)；未映射 %d (%.1f%%)，其中仍为 XLOC 的 %d 个（基因组上对不上已知基因的新组装位点）。",
  n_map, n_input, 100 - pct_fail, n_fail, pct_fail, n_xloc_left
))
utils::write.csv(
  merge(
    de[, c("gene", "gene_raw", "test_id", "locus", "xloc_rescued", "map_via"), drop = FALSE],
    gene_entrez_map, by = "gene", all.x = TRUE
  ),
  file.path(log_dir, "ID_mapping_all_genes.csv"),
  row.names = FALSE
)
if (n_fail > 0) {
  um <- de[de$gene %in% idmap$unmapped, c("gene", "gene_raw", "test_id", "locus"), drop = FALSE]
  utils::write.csv(um, file.path(log_dir, "ID_unmapped_genes.csv"), row.names = FALSE)
}
writeLines(
  c(
    "为什么 XLOC_000001、XLOC_000002 不能直接转成 Entrez：",
    "  XLOC_ 是 Cufflinks 给每个组装出来的基因组区间编的流水号。",
    "  后面的数字不同 = 染色体上不同的一段（locus 不同），不是同一个基因的别名。",
    "  NCBI org.Hs.eg.db 里没有 XLOC 这种 ID，所以 bitr(SYMBOL→ENTREZID) 一定会失败。",
    "",
    "脚本做法：gene 列已是官方符号的，用 SYMBOL/ALIAS 转 Entrez；",
    "仍是 XLOC 的，用 Excel 的 locus（如 chr1:11873-14409）去和 hg38/hg19 已知基因重叠。",
    "重叠上的会改成官方符号并进入 GO/KEGG；对不上的是新位点，留在差异表/火山图/热图，不进通路。",
    "",
    sprintf("mapped=%d  unmapped=%d  still_XLOC=%d  fail_rate=%.2f%%", n_map, n_fail, n_xloc_left, pct_fail),
    "全表映射: ID_mapping_all_genes.csv",
    "未映射名单: ID_unmapped_genes.csv",
    "坐标救回的 XLOC: XLOC_rescued_by_locus.csv（若有）"
  ),
  file.path(log_dir, "ID_mapping_README.txt")
)

map_to_entrez <- function(symbols) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  if (length(symbols) == 0) return(data.frame(gene = character(), entrez = character()))
  hit <- gene_entrez_map[gene_entrez_map$gene %in% symbols, , drop = FALSE]
  if (nrow(hit) == 0) return(data.frame(gene = character(), entrez = character()))
  hit[!duplicated(hit$gene), , drop = FALSE]
}

ranked_entrez <- function(df) {
  mp <- map_to_entrez(df$gene)
  de2 <- merge(df, mp, by = "gene")
  de2 <- de2[!is.na(de2$entrez) & !is.na(de2$log2FC), ]
  de2 <- de2[order(abs(de2$log2FC), decreasing = TRUE), ]
  de2 <- de2[!duplicated(de2$entrez), ]
  stats <- de2$log2FC
  names(stats) <- de2$entrez
  sort(stats, decreasing = TRUE)
}

save_gg <- function(plot, path_stub, width = 8, height = 6) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".pdf"), plot, width = width, height = height),
           error = function(e) log_msg("pdf ggsave failed: ", e$message))
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".png"), plot, width = width, height = height, dpi = 300),
           error = function(e) log_msg("png ggsave failed: ", e$message))
}

try_save_plot <- function(fun, stub, width = 9, height = 7) {
  p <- tryCatch(fun(), error = function(e) {
    log_msg("Plot failed (", basename(stub), "): ", e$message)
    NULL
  })
  if (is.null(p)) return(invisible(FALSE))
  tryCatch({
    save_gg(p, stub, width = width, height = height)
    TRUE
  }, error = function(e) {
    log_msg("ggsave failed (", basename(stub), "): ", e$message)
    FALSE
  })
}

note_empty <- function(stub, msg) {
  dir.create(dirname(stub), recursive = TRUE, showWarnings = FALSE)
  writeLines(msg, paste0(stub, "_EMPTY.txt"))
}

plot_volcano <- function(full_de, highlight, title, outfile, fc_line = 1, p_line = 0.05) {
  df <- full_de
  df$y <- -log10(pmax(df$pvalue, 1e-300))
  df$set <- "other"
  up <- highlight[highlight %in% df$gene[df$log2FC >= 0]]
  down <- highlight[highlight %in% df$gene[df$log2FC < 0]]
  df$set[df$gene %in% up] <- "up_in_subset"
  df$set[df$gene %in% down] <- "down_in_subset"
  lab_n <- min(15, length(highlight))
  df$label <- ifelse(df$gene %in% utils::head(highlight, lab_n), df$gene, NA)
  lfc_line <- log2(fc_line)
  pal <- c(other = "grey70", up_in_subset = "#D62828", down_in_subset = "#1D4ED8")
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = y, color = set)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.3) +
    ggplot2::scale_color_manual(
      values = pal, drop = FALSE,
      labels = c(other = "other", up_in_subset = "up in subset", down_in_subset = "down in subset")
    ) +
    ggplot2::geom_vline(xintercept = c(-lfc_line, lfc_line), linetype = 2, color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(p_line), linetype = 2, color = "grey40") +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = title,
      subtitle = paste0(
        "highlight up n=", length(up), ", down n=", length(down),
        " | dashed p=", p_line, " | FC line=", fc_line
      ),
      x = "log2 Fold Change (Excel)", y = "-log10(p_value)", color = NULL
    )
  if (isTRUE(abs(p_line - 0.05) < 1e-12)) {
    p <- p + ggplot2::geom_hline(yintercept = -log10(0.01), linetype = 3, color = "grey55")
  }
  save_gg(p, outfile)
}

plot_heatmap <- function(mat, si, genes, title, outfile, lfc = NULL) {
  genes <- intersect(genes, rownames(mat))
  if (length(genes) > 400) {
    log_msg("热图只画该子集 |log2FC| 最高的 400 / ", length(genes), " 个基因: ", title)
    if (is.null(lfc)) {
      lfc <- setNames(rep(0, length(genes)), genes)
    }
    ord <- order(abs(as.numeric(lfc[genes])), decreasing = TRUE, na.last = TRUE)
    genes <- genes[ord][seq_len(400)]
  }
  if (length(genes) < 2) {
    log_msg("Heatmap skipped (<2 genes): ", title)
    note_empty(outfile, "fewer than 2 genes")
    return(invisible(NULL))
  }
  sub <- mat[genes, , drop = FALSE]
  ann <- data.frame(Group = si$group, row.names = si$sample)
  ann <- ann[colnames(sub), , drop = FALSE]
  pal <- c(NTC = "#4C78A8", TG_sh1 = "#F58518")
  draw_hm <- function() {
    args <- list(
      mat = sub, scale = "row", annotation_col = ann,
      annotation_colors = list(Group = pal[names(pal) %in% unique(ann$Group)]),
      show_rownames = nrow(sub) <= 80, fontsize_row = 6, main = title,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100),
      cluster_cols = ncol(sub) >= 2
    )
    tryCatch(
      do.call(pheatmap::pheatmap, c(args, list(clustering_distance_rows = "correlation"))),
      error = function(e) do.call(pheatmap::pheatmap, c(args, list(clustering_distance_rows = "euclidean")))
    )
  }
  grDevices::pdf(paste0(outfile, ".pdf"), width = 8, height = max(6, min(18, 0.18 * nrow(sub) + 3)))
  on.exit({
    while (grDevices::dev.cur() > 1) grDevices::dev.off()
  }, add = TRUE)
  draw_hm()
  grDevices::dev.off()
  grDevices::png(paste0(outfile, ".png"), width = 2400, height = max(1800, 40 * nrow(sub) + 400), res = 300)
  draw_hm()
  grDevices::dev.off()
}

plot_de_bar <- function(sub, title, outfile) {
  if (nrow(sub) == 0) return(invisible(NULL))
  df <- sub[order(sub$log2FC, decreasing = TRUE), , drop = FALSE]
  if (nrow(df) > 60) df <- rbind(utils::head(df, 30), utils::tail(df, 30))
  df$gene <- factor(df$gene, levels = rev(unique(df$gene)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = gene, y = log2FC, fill = log2FC >= 0)) +
    ggplot2::geom_col(width = 0.8, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c("FALSE" = "#1D4ED8", "TRUE" = "#D62828")) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title, x = NULL, y = "log2 Fold Change")
  save_gg(p, outfile, width = 8, height = max(5, min(16, 0.22 * nrow(df) + 2)))
}

# -----------------------------------------------------------------------------
# 4. ORA / GSEA / GSVA
# -----------------------------------------------------------------------------
enrich_or_relax <- function(fun_strict, fun_relax, label) {
  obj <- tryCatch(fun_strict(), error = function(e) {
    log_msg(label, " strict failed: ", e$message)
    NULL
  })
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    attr(obj, "relaxed") <- FALSE
    return(obj)
  }
  obj <- tryCatch(fun_relax(), error = function(e) {
    log_msg(label, " relaxed failed: ", e$message)
    NULL
  })
  if (!is.null(obj)) attr(obj, "relaxed") <- TRUE
  obj
}

title_maybe_relaxed <- function(obj, base) {
  if (isTRUE(attr(obj, "relaxed"))) paste0(base, " (relaxed cutoff)") else base
}

msig_tbl <- function(collection, category = NULL, subcategory = NULL) {
  if (!is.null(category)) {
    tryCatch(
      msigdbr::msigdbr(species = "Homo sapiens", collection = category, subcollection = subcategory),
      error = function(e1) {
        tryCatch(
          msigdbr::msigdbr(species = "Homo sapiens", category = category, subcategory = subcategory),
          error = function(e2) NULL
        )
      }
    )
  } else {
    tryCatch(
      msigdbr::msigdbr(species = "Homo sapiens", collection = collection),
      error = function(e1) {
        tryCatch(
          msigdbr::msigdbr(species = "Homo sapiens", category = collection),
          error = function(e2) NULL
        )
      }
    )
  }
}

msig_hallmark_map <- function() {
  msig <- msig_tbl("H")
  if (is.null(msig) || nrow(msig) == 0) stop("msigdbr Hallmark 为空")
  gs <- if ("gs_name" %in% names(msig)) msig$gs_name else msig[[grep("gs_name|gs_id", names(msig))[1]]]
  id <- if ("ncbi_gene" %in% names(msig)) msig$ncbi_gene else if ("entrez_gene" %in% names(msig)) msig$entrez_gene else msig$entrez_gene
  data.frame(term = gs, entrez = as.character(id), stringsAsFactors = FALSE)
}

msig_symbol_sets <- function() {
  out <- list()
  hm <- msig_tbl("H")
  if (!is.null(hm) && nrow(hm) > 0) {
    sym <- if ("gene_symbol" %in% names(hm)) hm$gene_symbol else hm$symbol
    out$Hallmark <- split(as.character(sym), hm$gs_name)
  }
  kegg <- msig_tbl(NULL, category = "C2", subcategory = "CP:KEGG")
  if (is.null(kegg) || nrow(kegg) == 0) {
    kegg <- msig_tbl(NULL, category = "C2", subcategory = "KEGG")
  }
  if (!is.null(kegg) && nrow(kegg) > 0) {
    sym <- if ("gene_symbol" %in% names(kegg)) kegg$gene_symbol else kegg$symbol
    out$KEGG <- split(as.character(sym), kegg$gs_name)
  }
  out
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
  x2 <- tryCatch(enrichplot::pairwise_termsim(x), error = function(e) NULL)
  if (!is.null(x2) && nrow(df) >= 2) {
    try_save_plot(function() enrichplot::emapplot(x2, showCategory = min(20, nrow(df))) + ggplot2::ggtitle(title),
                  paste0(stub, "_emapplot"), 10, 8)
  }
  try_save_plot(function() enrichplot::cnetplot(
    x, showCategory = min(8, nshow), foldChange = fold_change, circular = FALSE
  ) + ggplot2::ggtitle(title), paste0(stub, "_cnetplot"), 10, 8)
}

plot_gsea_object <- function(x, stub, title) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) {
    note_empty(stub, "no GSEA terms")
    return(invisible(NULL))
  }
  df <- as.data.frame(x)
  utils::write.csv(df, paste0(stub, ".csv"), row.names = FALSE)
  nshow <- min(15, nrow(df))
  try_save_plot(function() {
    p <- enrichplot::dotplot(x, showCategory = nshow, split = ".sign")
    tryCatch(p + ggplot2::facet_grid(. ~ .sign) + ggplot2::ggtitle(title),
             error = function(e) p + ggplot2::ggtitle(title))
  }, paste0(stub, "_dotplot"), 10, 7)
  try_save_plot(function() enrichplot::ridgeplot(x, showCategory = nshow) + ggplot2::ggtitle(title),
                paste0(stub, "_ridgeplot"), 10, 8)
  ncurve <- min(5, nrow(df))
  try_save_plot(function() enrichplot::gseaplot2(x, geneSetID = seq_len(ncurve), pvalue_table = TRUE, title = title),
                paste0(stub, "_gseaplot"), 10, 8)
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
    log_msg("ORA skipped, mapped genes < 3: ", outdir)
    note_empty(file.path(go_dir, paste0(pref, "ORA_GO")), "too few mapped genes")
    note_empty(file.path(pw_dir, paste0(pref, "ORA_Pathway")), "too few mapped genes")
    note_empty(file.path(kg_dir, paste0(pref, "ORA_KEGG")), "too few mapped genes")
    return(invisible(NULL))
  }
  universe <- unique(map_to_entrez(de$gene)$entrez)

  for (ont in c("BP", "MF", "CC")) {
    ego <- enrich_or_relax(
      function() clusterProfiler::enrichGO(
        gene = entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
        universe = universe, pAdjustMethod = "BH", pvalueCutoff = 0.05,
        qvalueCutoff = 0.2, readable = TRUE
      ),
      function() clusterProfiler::enrichGO(
        gene = entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
        universe = universe, pAdjustMethod = "BH", pvalueCutoff = 1,
        qvalueCutoff = 1, readable = TRUE
      ),
      paste("enrichGO", ont)
    )
    plot_ora_object(ego, file.path(go_dir, paste0(pref, "ORA_GO_", ont)),
                    title_maybe_relaxed(ego, paste(label, "| ORA GO", ont)), fold_change = fc_sym)
  }

  ek <- enrich_or_relax(
    function() clusterProfiler::enrichKEGG(
      gene = entrez, organism = "hsa", universe = universe,
      pvalueCutoff = 0.05, qvalueCutoff = 0.2
    ),
    function() clusterProfiler::enrichKEGG(
      gene = entrez, organism = "hsa", universe = universe,
      pvalueCutoff = 1, qvalueCutoff = 1
    ),
    "enrichKEGG"
  )
  if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) {
    ek <- tryCatch(clusterProfiler::setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
                   error = function(e) ek)
  }
  plot_ora_object(ek, file.path(kg_dir, paste0(pref, "ORA_KEGG")),
                  title_maybe_relaxed(ek, paste(label, "| ORA KEGG")), fold_change = fc_sym)

  if (has_pkg("ReactomePA")) {
    er <- enrich_or_relax(
      function() ReactomePA::enrichPathway(
        gene = entrez, organism = "human", universe = universe,
        pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
      ),
      function() ReactomePA::enrichPathway(
        gene = entrez, organism = "human", universe = universe,
        pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
      ),
      "enrichPathway"
    )
    plot_ora_object(er, file.path(pw_dir, paste0(pref, "ORA_Reactome_pathway")),
                    title_maybe_relaxed(er, paste(label, "| ORA Reactome")), fold_change = fc_sym)
  } else {
    note_empty(file.path(pw_dir, paste0(pref, "ORA_Reactome_pathway")), "ReactomePA not installed")
  }

  hm <- enrich_or_relax(
    function() {
      term2gene <- msig_hallmark_map()
      clusterProfiler::enricher(entrez, TERM2GENE = term2gene, universe = universe,
                                pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    },
    function() {
      term2gene <- msig_hallmark_map()
      clusterProfiler::enricher(entrez, TERM2GENE = term2gene, universe = universe,
                                pvalueCutoff = 1, qvalueCutoff = 1)
    },
    "Hallmark"
  )
  if (!is.null(hm) && nrow(as.data.frame(hm)) > 0) {
    hm <- tryCatch(clusterProfiler::setReadable(hm, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
                   error = function(e) hm)
  }
  plot_ora_object(hm, file.path(pw_dir, paste0(pref, "ORA_MSigDB_Hallmark_pathway")),
                  title_maybe_relaxed(hm, paste(label, "| ORA Hallmark")), fold_change = fc_sym)
  writeLines(
    c("This GO/Pathway/KEGG folder is ORA, NOT GSEA.",
      "GSEA files are in GSEA/ and start with GSEA_.",
      "GSVA files are in GSVA/ and start with GSVA_."),
    file.path(outdir, paste0(pref, "00_ORA_is_not_GSEA.txt"))
  )
}

run_gsea_plots <- function(sub, outdir, tag, label) {
  gsea_dir <- file.path(outdir, "GSEA")
  dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
  pref <- paste0(tag, "_")
  stats <- ranked_entrez(sub)
  if (length(stats) < 8) {
    note_empty(file.path(gsea_dir, paste0(pref, "GSEA")), "too few ranked genes")
    return(invisible(NULL))
  }
  term2gene <- tryCatch(msig_hallmark_map(), error = function(e) NULL)
  if (!is.null(term2gene)) {
    hm <- enrich_or_relax(
      function() clusterProfiler::GSEA(
        geneList = stats, TERM2GENE = term2gene, minGSSize = 5,
        maxGSSize = 500, pvalueCutoff = 0.05, eps = 0, verbose = FALSE
      ),
      function() clusterProfiler::GSEA(
        geneList = stats, TERM2GENE = term2gene, minGSSize = 3,
        maxGSSize = 500, pvalueCutoff = 1, eps = 0, verbose = FALSE
      ),
      paste("Hallmark GSEA", tag)
    )
    plot_gsea_object(hm, file.path(gsea_dir, paste0(pref, "GSEA_Hallmark")),
                     paste(label, "| GSEA Hallmark"))
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
    paste("KEGG GSEA", tag)
  )
  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
    kegg <- tryCatch(clusterProfiler::setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
                     error = function(e) kegg)
  }
  plot_gsea_object(kegg, file.path(gsea_dir, paste0(pref, "GSEA_KEGG")),
                   paste(label, "| GSEA KEGG"))
  for (ont in c("BP", "MF", "CC")) {
    go <- enrich_or_relax(
      function() clusterProfiler::gseGO(
        geneList = stats, OrgDb = org.Hs.eg.db, ont = ont, keyType = "ENTREZID",
        minGSSize = 10, pvalueCutoff = 0.05, verbose = FALSE, eps = 0
      ),
      function() clusterProfiler::gseGO(
        geneList = stats, OrgDb = org.Hs.eg.db, ont = ont, keyType = "ENTREZID",
        minGSSize = 5, pvalueCutoff = 1, verbose = FALSE, eps = 0
      ),
      paste("gseGO", ont, tag)
    )
    plot_gsea_object(go, file.path(gsea_dir, paste0(pref, "GSEA_GO_", ont)),
                     paste(label, "| GSEA GO", ont))
  }
}

run_gsva_one <- function(expr, gsets) {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"
  gsets <- lapply(gsets, function(x) intersect(unique(as.character(x)), rownames(expr)))
  gsets <- gsets[vapply(gsets, length, integer(1)) >= 3]
  if (length(gsets) == 0) return(NULL)
  if ("gsvaParam" %in% getNamespaceExports("GSVA")) {
    param <- GSVA::gsvaParam(exprData = expr, geneSets = gsets, kcdf = "Gaussian")
    GSVA::gsva(param, verbose = FALSE)
  } else {
    GSVA::gsva(expr, gsets, method = "gsva", kcdf = "Gaussian", verbose = FALSE)
  }
}

run_gsva_plots <- function(expr_mat, genes_keep, outdir, tag, title) {
  gsva_dir <- file.path(outdir, "GSVA")
  dir.create(gsva_dir, recursive = TRUE, showWarnings = FALSE)
  pref <- paste0(tag, "_GSVA_")
  genes_keep <- intersect(genes_keep, rownames(expr_mat))
  if (length(genes_keep) < 10) {
    note_empty(file.path(gsva_dir, paste0(pref, "all")), "too few genes for GSVA")
    return(invisible(NULL))
  }
  sub <- expr_mat[genes_keep, , drop = FALSE]
  sets <- tryCatch(msig_symbol_sets(), error = function(e) {
    log_msg("msig symbol sets failed: ", e$message)
    NULL
  })
  if (is.null(sets) || length(sets) == 0) {
    note_empty(file.path(gsva_dir, paste0(pref, "all")), "no gene sets")
    return(invisible(NULL))
  }
  for (nm in names(sets)) {
    scores <- tryCatch(run_gsva_one(sub, sets[[nm]]), error = function(e) {
      log_msg("GSVA ", nm, " failed: ", e$message)
      NULL
    })
    if (is.null(scores) || nrow(scores) == 0) {
      note_empty(file.path(gsva_dir, paste0(pref, nm)), "GSVA returned empty")
      next
    }
    scores <- as.matrix(scores)
    utils::write.csv(
      data.frame(pathway = rownames(scores), scores, check.names = FALSE),
      file.path(gsva_dir, paste0(pref, nm, "_scores.csv")),
      row.names = FALSE
    )
    delta <- scores[, 2] - scores[, 1]
    bar_df <- data.frame(
      pathway = rownames(scores),
      delta = as.numeric(delta),
      stringsAsFactors = FALSE
    )
    bar_df <- bar_df[order(abs(bar_df$delta), decreasing = TRUE), ]
    bar_df <- utils::head(bar_df, 30)
    bar_df$pathway <- factor(bar_df$pathway, levels = rev(bar_df$pathway))
    p <- ggplot2::ggplot(bar_df, ggplot2::aes(x = pathway, y = delta, fill = delta >= 0)) +
      ggplot2::geom_col(show.legend = FALSE) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c("FALSE" = "#1D4ED8", "TRUE" = "#D62828")) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::labs(
        title = paste(title, "| GSVA", nm, "|", colnames(scores)[2], "-", colnames(scores)[1]),
        x = NULL, y = "GSVA score difference"
      )
    save_gg(p, file.path(gsva_dir, paste0(pref, nm, "_delta_barplot")), 10, max(6, 0.22 * nrow(bar_df) + 2))
    nshow <- min(40, nrow(scores))
    keep_idx <- order(abs(delta), decreasing = TRUE)[seq_len(nshow)]
    hm <- scores[keep_idx, , drop = FALSE]
    grDevices::pdf(file.path(gsva_dir, paste0(pref, nm, "_heatmap.pdf")),
                   width = 8, height = max(6, min(16, 0.22 * nrow(hm) + 3)))
    pheatmap::pheatmap(
      hm, cluster_cols = FALSE, scale = "none",
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100),
      main = paste(title, "| GSVA", nm), fontsize_row = 7
    )
    grDevices::dev.off()
    grDevices::png(file.path(gsva_dir, paste0(pref, nm, "_heatmap.png")),
                   width = 2400, height = max(1800, 40 * nrow(hm) + 400), res = 300)
    pheatmap::pheatmap(
      hm, cluster_cols = FALSE, scale = "none",
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100),
      main = paste(title, "| GSVA", nm), fontsize_row = 7
    )
    grDevices::dev.off()
  }
}

# -----------------------------------------------------------------------------
# 5. 子集分析
# -----------------------------------------------------------------------------
select_by_fc <- function(df, fc) {
  df[!is.na(df$log2FC) & df$log2FC >= log2(fc), , drop = FALSE]
}

select_by_topn <- function(df, n) {
  x <- df[!is.na(df$log2FC) & df$log2FC > 0, , drop = FALSE]
  x <- x[order(x$log2FC, decreasing = TRUE), , drop = FALSE]
  if (nrow(x) < n) {
    log_msg("上调基因不足 ", n, "，实际取 ", nrow(x))
  }
  utils::head(x, n)
}

select_by_p <- function(df, pcut) {
  df[!is.na(df$pvalue) & df$pvalue < pcut, , drop = FALSE]
}

emit_subset <- function(sub, tag, title, outdir, fc_line = 1, p_line = 0.05) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  n_up <- sum(!is.na(sub$log2FC) & sub$log2FC >= 0)
  n_down <- sum(!is.na(sub$log2FC) & sub$log2FC < 0)
  sub$subset_tag <- tag
  sub$mark_this_subset <- TRUE
  writeLines(
    c(paste("comparison:", comp_name),
      paste("subset:", tag),
      paste("title:", title),
      paste("n_genes:", nrow(sub)),
      paste("n_up:", n_up),
      paste("n_down:", n_down),
      paste("excel:", excel_path),
      "思路3（Pvalue）不区分上调下调：上下调都进入本档分析，表里用 direction / mark_p_lt_* 标记。"),
    file.path(outdir, paste0("00_", tag, "_THIS_FOLDER.txt"))
  )
  utils::write.csv(sub, file.path(outdir, paste0(tag, "_DE_selected_genes.csv")), row.names = FALSE)
  tryCatch(writexl::write_xlsx(sub, file.path(outdir, paste0(tag, "_DE_selected_genes.xlsx"))),
           error = function(e) log_msg("xlsx write failed: ", e$message))
  log_msg(comp_name, " ", tag, ": n = ", nrow(sub))
  if (nrow(sub) == 0) {
    writeLines("no genes", file.path(outdir, paste0(tag, "_EMPTY.txt")))
    return(invisible(NULL))
  }
  tryCatch(plot_de_bar(sub, paste(title, "| DE genes"), file.path(outdir, paste0(tag, "_DE_log2FC_barplot"))),
           error = function(e) log_msg("DE barplot failed: ", e$message))
  tryCatch(plot_volcano(de, sub$gene, title, file.path(outdir, paste0(tag, "_volcano")), fc_line = fc_line, p_line = p_line),
           error = function(e) log_msg("volcano failed: ", e$message))
  tryCatch(plot_heatmap(
    heat_log, sample_info, sub$gene, title,
    file.path(outdir, paste0(tag, "_heatmap")),
    lfc = setNames(sub$log2FC, sub$gene)
  ),
           error = function(e) {
             while (grDevices::dev.cur() > 1) grDevices::dev.off()
             log_msg("heatmap failed: ", e$message)
           })
  tryCatch(run_ora_plots(sub$gene, sub, outdir, title, tag),
           error = function(e) log_msg("ORA failed: ", e$message))
  tryCatch(run_gsea_plots(sub, outdir, tag, title),
           error = function(e) log_msg("GSEA failed: ", e$message))
  tryCatch(run_gsva_plots(heat_log, sub$gene, outdir, tag, title),
           error = function(e) {
             while (grDevices::dev.cur() > 1) grDevices::dev.off()
             log_msg("GSVA failed: ", e$message)
           })
}

# -----------------------------------------------------------------------------
# 6. 主流程：4 个 FC + 8 个 topN + p<0.01 / p<0.05（思路3，不区分上下调）
# -----------------------------------------------------------------------------
de$mark_p_lt_0.01 <- !is.na(de$pvalue) & de$pvalue < 0.01
de$mark_p_lt_0.05 <- !is.na(de$pvalue) & de$pvalue < 0.05
de$direction <- ifelse(
  is.na(de$log2FC), NA_character_,
  ifelse(de$log2FC > 0, "up", ifelse(de$log2FC < 0, "down", "zero"))
)
log_msg(
  "思路3标记: p<0.01 n=", sum(de$mark_p_lt_0.01),
  " (up=", sum(de$mark_p_lt_0.01 & de$direction == "up"),
  ", down=", sum(de$mark_p_lt_0.01 & de$direction == "down"), "); ",
  "p<0.05 n=", sum(de$mark_p_lt_0.05),
  " (up=", sum(de$mark_p_lt_0.05 & de$direction == "up"),
  ", down=", sum(de$mark_p_lt_0.05 & de$direction == "down"), ")"
)

base <- file.path(result_dir, comp_name)
dir.create(base, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(de, file.path(base, "DE_full_from_excel.csv"), row.names = FALSE)
tryCatch(writexl::write_xlsx(de, file.path(base, "DE_full_from_excel.xlsx")),
         error = function(e) log_msg("DE_full xlsx failed: ", e$message))

writeLines(
  c("OLD盘 Excel：NTC vs TG_sh1",
    paste("输入文件:", excel_path),
    "只保留 status==OK，直接用表中 log2(fold_change) 和 p_value。",
    "",
    "三套思路都会做，不是只做倍数。脚本按顺序跑，文件夹会陆续出现：",
    "  1) TopRank/top50 ... top300, top4000   （上调排名，先跑，较快）",
    "  2) FoldChange/FC_2  FC_1.5  FC_1.25  FC_1  （上调倍数；FC_1 基因最多、最慢）",
    "  3) Pvalue/p_lt_0.01  与  Pvalue/p_lt_0.05  （p<0.01、p<0.05，上下调都保留并做标记）",
    "",
    "如果现在只能看到 FoldChange，多半是还在跑 FC 各档的 GO/KEGG/GSEA/GSVA，",
    "请看 00_PROGRESS.txt，不要中途关掉 R。",
    "",
    "每个非空子文件夹里应有：差异表、火山图、热图、GO/、Pathway/、KEGG/、GSEA/、GSVA/。",
    "00_GSEA_all_genes_NOT_FC_or_topN 是全表基因 GSEA，不是分层结果。",
    "bitr 未映射警告已抑制：全表只映射一次，见 00_logs/ID_mapping_README.txt。"),
  file.path(base, "00_READ_ME_先看这里.txt")
)

fc_dirs <- file.path(base, "FoldChange", names(fc_cutoffs))
top_dirs <- file.path(base, "TopRank", paste0("top", top_ns))
p_dirs <- file.path(base, "Pvalue", names(p_cutoffs))
invisible(lapply(c(fc_dirs, top_dirs, p_dirs), dir.create, recursive = TRUE, showWarnings = FALSE))
placeholder <- function(d, label) {
  f <- file.path(d, "00_WAITING_脚本还没跑到这一档.txt")
  if (!file.exists(file.path(d, paste0("00_", label, "_THIS_FOLDER.txt")))) {
    writeLines(
      c("这个文件夹一开始就会建好，避免只看到 FoldChange。",
        "脚本按 TopRank → FoldChange → Pvalue 的顺序往里面写图。",
        "请看上一级的 00_PROGRESS.txt。"),
      f
    )
  }
}
for (nm in names(fc_cutoffs)) placeholder(file.path(base, "FoldChange", nm), nm)
for (n in top_ns) placeholder(file.path(base, "TopRank", paste0("top", n)), paste0("top", n))
for (nm in names(p_cutoffs)) placeholder(file.path(base, "Pvalue", nm), nm)

progress_file <- file.path(base, "00_PROGRESS.txt")
write_progress <- function(...) {
  writeLines(
    c(paste("更新时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste(..., collapse = ""),
      "",
      "计划顺序：TopRank（8 档）→ FoldChange（4 档，从 FC_2 到 FC_1）→ Pvalue/p_lt_0.01 与 p_lt_0.05 → 全基因 GSEA/GSVA",
      "只看到 FoldChange 不等于没做排名和 p 值，通常是倍数分析还没跑完。"),
    progress_file
  )
}

jobs <- list()
for (n in top_ns) {
  tag <- paste0("top", n)
  jobs[[length(jobs) + 1]] <- list(
    kind = "TopRank", tag = tag,
    dir = file.path(base, "TopRank", tag),
    fc_line = 1, p_line = 0.05,
    title = paste(comp_name, "| upregulated", tag),
    sub = select_by_topn(de, n)
  )
}
for (nm in rev(names(fc_cutoffs))) {
  fc <- unname(fc_cutoffs[[nm]])
  jobs[[length(jobs) + 1]] <- list(
    kind = "FoldChange", tag = nm,
    dir = file.path(base, "FoldChange", nm),
    fc_line = fc, p_line = 0.05,
    title = paste(comp_name, "|", nm, "| up FC >=", fc),
    sub = select_by_fc(de, fc)
  )
}
for (nm in names(p_cutoffs)) {
  pc <- unname(p_cutoffs[[nm]])
  jobs[[length(jobs) + 1]] <- list(
    kind = "Pvalue", tag = nm,
    dir = file.path(base, "Pvalue", nm),
    fc_line = 1, p_line = pc,
    title = paste(comp_name, "|", nm, "| p_value <", pc, "(up and down, marked)"),
    sub = select_by_p(de, pc)
  )
}

n_jobs <- length(jobs)
write_progress("已建好 TopRank / FoldChange / Pvalue 三个文件夹，开始跑第 1 / ", n_jobs, " 档")
for (i in seq_len(n_jobs)) {
  job <- jobs[[i]]
  write_progress(
    "正在做 ", i, " / ", n_jobs, " : ", job$kind, " / ", job$tag,
    " (n=", nrow(job$sub), ")"
  )
  wait_f <- file.path(job$dir, "00_WAITING_脚本还没跑到这一档.txt")
  if (file.exists(wait_f)) unlink(wait_f)
  tryCatch(
    emit_subset(job$sub, job$tag, job$title, job$dir, fc_line = job$fc_line, p_line = job$p_line),
    error = function(e) {
      while (grDevices::dev.cur() > 1) grDevices::dev.off()
      log_msg("ERROR subset ", job$kind, " ", job$tag, ": ", e$message)
      writeLines(paste("failed:", e$message), file.path(job$dir, paste0(job$tag, "_FAILED.txt")))
    }
  )
}

full_gsea_dir <- file.path(base, "00_GSEA_all_genes_NOT_FC_or_topN")
dir.create(full_gsea_dir, recursive = TRUE, showWarnings = FALSE)
writeLines("全表基因 GSEA / GSVA，不是 FC/topN 分层结果。",
           file.path(full_gsea_dir, "00_README.txt"))
write_progress("分层已写完，正在做全表基因 GSEA/GSVA")
tryCatch(run_gsea_plots(de, full_gsea_dir, "allGenes", paste(comp_name, "| all genes")),
         error = function(e) log_msg("full GSEA failed: ", e$message))
tryCatch(run_gsva_plots(heat_log, de$gene, full_gsea_dir, "allGenes", paste(comp_name, "| all genes")),
         error = function(e) log_msg("full GSVA failed: ", e$message))

write_progress("全部完成。请打开 TopRank、FoldChange、Pvalue 三个文件夹。")
log_msg("完成。请打开: ", base)
log_msg("先看 00_READ_ME_先看这里.txt 和 00_PROGRESS.txt")
