#!/usr/bin/env Rscript
#
# TG knockdown RNA-seq pipeline (NTC, TG_sh1, TG_sh5)
#
# Three comparison tracks:
#   1) TG_sh1 vs NTC and TG_sh5 vs NTC
#   2) mean(TG_sh1, TG_sh5) vs NTC
#   3) genes upregulated in BOTH sh1 vs NTC and sh5 vs NTC
#
# Two filter strategies (upregulated only), applied to every track:
#   A) linear FC >= 1, 1.25, 1.5, 2
#   B) top 50 / 75 / 100 / 150 / 200 / 250 / 300 by log2FC
#
# For each filter: DEG table, volcano (or sh1-vs-sh5 scatter for common-up),
# heatmap, GO ORA, Reactome ORA, KEGG ORA. GSEA is run once per comparison
# on the full ranked list.
#
# Usage (Windows R console — do not rely on Documents as the working directory):
#   setwd("E:/R/TG_BRCA/TG")
#   source("R/analyze_tg_rnaseq.R")
# or:
#   Rscript R/analyze_tg_rnaseq.R "E:/R/TG_BRCA/TG" "E:/R/TG_BRCA/TG/results"
#
# If 00_config.R is missing (e.g. the script was pasted into the console),
# built-in defaults are used and data are read from E:/R/TG_BRCA/TG when present.

# ---------------------------------------------------------------------------
# Built-in config (overridden if 00_config.R is found)
# ---------------------------------------------------------------------------
CTRL_GROUP <- "NTC"
TREAT_GROUPS <- c("TG_sh1", "TG_sh5")
ALL_GROUPS <- c(CTRL_GROUP, TREAT_GROUPS)
FC_CUTOFFS <- c(1, 1.25, 1.5, 2)
TOP_N <- c(50, 75, 100, 150, 200, 250, 300)
PSEUDOCOUNT <- 0.1
MIN_FPKM <- 1
HEATMAP_MAX_GENES <- 100
ORA_MIN_GENES <- 10
GSEA_MIN_GENES <- 15
GSEA_PVALUE_CUTOFF <- 1
SPECIES_ORGDB <- "org.Hs.eg.db"
KEGG_ORG <- "hsa"
P_ADJUST <- "BH"
ENRICH_P_CUTOFF <- 0.05
ENRICH_Q_CUTOFF <- 0.20
EXCEL_PATTERN <- "shTG"
DIFF_FILENAME <- "gene_exp.diff"
FPKM_FILENAME <- "genes.fpkm_tracking"

default_data_dir <- function() {
  candidates <- c(
    Sys.getenv("TG_RNASEQ_DIR", unset = ""),
    "E:/R/TG_BRCA/TG",
    "E:\\R\\TG_BRCA\\TG",
    file.path(".", "data")
  )
  candidates <- candidates[nzchar(candidates)]
  for (p in candidates) {
    if (dir.exists(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  normalizePath("data", winslash = "/", mustWork = FALSE)
}

default_out_dir <- function(data_dir) {
  file.path(data_dir, "results")
}

script_dirs_guess <- function() {
  dirs <- character()
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_full, value = TRUE)
  if (length(file_arg)) {
    fp <- sub("^--file=", "", file_arg[[1]])
    if (nzchar(fp) && file.exists(fp)) {
      dirs <- c(dirs, dirname(normalizePath(fp, winslash = "/", mustWork = FALSE)))
    }
  }
  if (sys.nframe() > 0) {
    for (i in sys.nframe():1) {
      ev <- sys.frame(i)
      if (exists("ofile", envir = ev, inherits = FALSE)) {
        of <- get("ofile", envir = ev)
        if (is.character(of) && length(of) && nzchar(of[[1]]) && file.exists(of[[1]])) {
          dirs <- c(dirs, dirname(normalizePath(of[[1]], winslash = "/", mustWork = FALSE)))
        }
      }
    }
  }
  wd <- tryCatch(normalizePath(getwd(), winslash = "/", mustWork = FALSE), error = function(e) getwd())
  unique(c(
    dirs,
    file.path(dirs, ".."),
    wd,
    file.path(wd, "R"),
    "E:/R/TG_BRCA/TG",
    "E:/R/TG_BRCA/TG/R"
  ))
}

config_file <- NA_character_
for (d in script_dirs_guess()) {
  if (!is.character(d) || !nzchar(d)) next
  cand <- file.path(d, "00_config.R")
  if (file.exists(cand)) {
    config_file <- cand
    break
  }
}
if (!is.na(config_file)) {
  message("Loading config: ", config_file)
  sys.source(config_file, envir = environment())
} else {
  message("00_config.R not found; using built-in defaults.")
  message("Data directory default: E:/R/TG_BRCA/TG (if that folder exists).")
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (is.atomic(a) && all(is.na(a)))) b else a

need_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

optional_library <- function(pkg) {
  if (need_pkg(pkg)) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    TRUE
  } else {
    message("Optional package not installed: ", pkg)
    FALSE
  }
}

# Cuffdiff gene_short_name is often empty (XLOC_*), or several symbols joined by commas.
is_xloc_id <- function(x) {
  grepl("^(XLOC_|CUFF\\.|TCONS_)", as.character(x), ignore.case = TRUE)
}

split_gene_tokens <- function(x) {
  x <- as.character(x[[1]])
  x <- gsub("///", ",", x)
  toks <- unlist(strsplit(x, "[,;|/]+"), use.names = FALSE)
  toks <- trimws(toks)
  toks <- toks[!is.na(toks) & nzchar(toks) & toks != "-"]
  toks <- toks[!tolower(toks) %in% c("na", "nan", "none", "n/a")]
  unique(toks)
}

primary_gene_symbol <- function(x) {
  vapply(as.character(x), function(one) {
    if (is.na(one) || !nzchar(one)) return(NA_character_)
    toks <- split_gene_tokens(one)
    if (!length(toks)) return(NA_character_)
    ann <- toks[!is_xloc_id(toks)]
    if (length(ann)) ann[[1]] else toks[[1]]
  }, character(1), USE.NAMES = FALSE)
}

flatten_symbols <- function(symbols) {
  unique(unlist(lapply(as.character(symbols), function(s) {
    toks <- split_gene_tokens(s)
    toks[!is_xloc_id(toks)]
  }), use.names = FALSE))
}

collapse_expr_rows <- function(mat, genes) {
  raw <- as.character(genes)
  primary <- primary_gene_symbol(raw)
  keep <- !is.na(primary) & nzchar(primary)
  mat <- mat[keep, , drop = FALSE]
  raw <- raw[keep]
  primary <- primary[keep]
  if (!length(primary)) stop("No gene labels after cleaning Cuffdiff names")
  ord <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
  mat <- mat[ord, , drop = FALSE]
  raw <- raw[ord]
  primary <- primary[ord]
  uniq <- !duplicated(primary)
  mat <- mat[uniq, , drop = FALSE]
  raw <- raw[uniq]
  primary <- primary[uniq]
  rownames(mat) <- primary
  attr(mat, "gene_raw") <- setNames(raw, primary)
  mat
}

# ---------------------------------------------------------------------------
# Input discovery and readers
# ---------------------------------------------------------------------------

find_first <- function(dir, pattern, fixed_name = NULL) {
  if (!is.null(fixed_name)) {
    direct <- file.path(dir, fixed_name)
    if (file.exists(direct)) return(direct)
  }
  hits <- list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (!length(hits)) return(NA_character_)
  hits[[1]]
}

normalize_colnames <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[()]", "", x)
  x <- gsub("[\\. ]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

read_table_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!need_pkg("readxl")) stop("Install readxl to read Excel: ", path)
    sheets <- readxl::excel_sheets(path)
    out <- lapply(sheets, function(sh) {
      df <- as.data.frame(readxl::read_excel(path, sheet = sh), stringsAsFactors = FALSE)
      attr(df, "sheet") <- sh
      df
    })
    names(out) <- sheets
    return(out)
  }
  sep <- if (ext %in% c("csv")) "," else "\t"
  df <- utils::read.delim(path, sep = sep, header = TRUE, check.names = FALSE,
                          stringsAsFactors = FALSE, na.strings = c("", "NA", "nan", "NaN"))
  list(main = df)
}

is_cuffdiff_diff <- function(df) {
  nms <- normalize_colnames(names(df))
  all(c("sample_1", "sample_2", "value_1", "value_2") %in% nms) ||
    all(c("log2_fold_change", "sample_1", "sample_2") %in% nms)
}

tidy_diff_df <- function(df) {
  nms <- names(df)
  nms_n <- normalize_colnames(nms)
  names(df) <- nms_n
  rename_map <- c(
    "log2fold_change" = "log2_fold_change",
    "log2_foldchange" = "log2_fold_change",
    "pvalue" = "p_value",
    "qvalue" = "q_value",
    "gene_short_name" = "gene"
  )
  for (old in names(rename_map)) {
    if (old %in% names(df) && !rename_map[[old]] %in% names(df)) {
      names(df)[names(df) == old] <- rename_map[[old]]
    }
  }
  if (!"gene" %in% names(df)) {
    if ("gene_short_name" %in% names(df)) df$gene <- df$gene_short_name
    else if ("gene_id" %in% names(df)) df$gene <- df$gene_id
    else df$gene <- df[[1]]
  }
  raw <- ifelse(is.na(df$gene) | trimws(as.character(df$gene)) == "",
                as.character(df$gene_id %||% df$test_id %||% seq_len(nrow(df))),
                trimws(as.character(df$gene)))
  df$gene_raw <- raw
  df$gene <- primary_gene_symbol(raw)
  num_cols <- intersect(c("value_1", "value_2", "log2_fold_change", "test_stat",
                          "p_value", "q_value"), names(df))
  for (cc in num_cols) {
    df[[cc]] <- suppressWarnings(as.numeric(gsub("inf", "Inf", as.character(df[[cc]]), ignore.case = TRUE)))
  }
  df
}

split_diff_by_comparison <- function(df) {
  df <- tidy_diff_df(df)
  if (!all(c("sample_1", "sample_2") %in% names(df))) {
    return(list(all = df))
  }
  key <- paste(df$sample_1, df$sample_2, sep = "_vs_")
  split(df, key)
}

guess_group_name <- function(x) {
  x <- trimws(as.character(x))
  x_u <- toupper(gsub("[^A-Za-z0-9]", "", x))
  if (grepl("NTC", x_u)) return("NTC")
  if (grepl("SH1|SHRNA1|TGSH1", x_u) || grepl("SH1", x_u)) return("TG_sh1")
  if (grepl("SH5|SHRNA5|TGSH5", x_u) || grepl("SH5", x_u)) return("TG_sh5")
  x
}

read_fpkm_tracking <- function(path) {
  df <- utils::read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
                          stringsAsFactors = FALSE)
  nms <- names(df)
  gene_col <- intersect(c("gene_short_name", "gene", "gene_id", "tracking_id"), nms)
  if (!length(gene_col)) stop("Cannot find gene column in ", path)
  fpkm_cols <- grep("_FPKM$", nms, value = TRUE)
  if (!length(fpkm_cols)) {
    fpkm_cols <- grep("FPKM", nms, value = TRUE, ignore.case = TRUE)
  }
  if (!length(fpkm_cols)) stop("Cannot find FPKM columns in ", path)
  genes <- as.character(df[[gene_col[[1]]]])
  if ("gene_id" %in% nms) {
    gid <- as.character(df$gene_id)
    genes <- ifelse(is.na(genes) | trimws(genes) == "" | genes == "-", gid, genes)
  }
  mat <- as.matrix(df[, fpkm_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  colnames(mat) <- vapply(sub("_FPKM$", "", colnames(mat)), guess_group_name, character(1))
  collapse_expr_rows(mat, genes)
}

fpkm_from_diff <- function(diff_df) {
  diff_df <- tidy_diff_df(diff_df)
  if (all(c("sample_1", "sample_2") %in% names(diff_df))) {
    g1 <- guess_group_name(diff_df$sample_1[[1]])
    g2 <- guess_group_name(diff_df$sample_2[[1]])
  } else {
    g1 <- CTRL_GROUP
    g2 <- TREAT_GROUPS[[1]]
  }
  mat <- cbind(as.numeric(diff_df$value_1), as.numeric(diff_df$value_2))
  colnames(mat) <- c(g1, g2)
  genes <- diff_df$gene_raw %||% diff_df$gene
  keep <- !is.na(genes) & genes != ""
  mat <- mat[keep, , drop = FALSE]
  genes <- genes[keep]
  d <- diff_df[keep, , drop = FALSE]
  mat <- collapse_expr_rows(mat, genes)
  d$gene <- primary_gene_symbol(d$gene_raw %||% d$gene)
  d <- d[match(rownames(mat), d$gene), , drop = FALSE]
  list(mat = mat, diff = d)
}

merge_fpkm <- function(pieces) {
  genes <- Reduce(union, lapply(pieces, rownames))
  groups <- unique(unlist(lapply(pieces, colnames)))
  out <- matrix(NA_real_, nrow = length(genes), ncol = length(groups),
                dimnames = list(genes, groups))
  raw_map <- character()
  for (m in pieces) {
    out[rownames(m), colnames(m)] <- m
    gr <- attr(m, "gene_raw")
    if (!is.null(gr)) {
      raw_map[names(gr)] <- unname(gr)
    }
  }
  if (length(raw_map)) attr(out, "gene_raw") <- raw_map[rownames(out)]
  out
}

load_inputs <- function(data_dir) {
  find_excel <- function(dir) {
    hits <- list.files(dir, pattern = EXCEL_PATTERN, full.names = TRUE, ignore.case = TRUE)
    hits <- hits[grepl("\\.(xlsx|xls|csv)$", hits, ignore.case = TRUE)]
    if (!length(hits)) return(NA_character_)
    xlsx <- hits[grepl("\\.xlsx$", hits, ignore.case = TRUE)]
    if (length(xlsx)) return(xlsx[[1]])
    hits[[1]]
  }
  excel_path <- find_excel(data_dir)
  diff_path <- find_first(data_dir, paste0("^", DIFF_FILENAME, "$"), DIFF_FILENAME)
  fpkm_path <- find_first(data_dir, paste0("^", FPKM_FILENAME, "$"), FPKM_FILENAME)

  fpkm_list <- list()
  diff_map <- list()

  ingest_diff <- function(df, src) {
    parts <- split_diff_by_comparison(df)
    for (nm in names(parts)) {
      d <- parts[[nm]]
      parsed <- fpkm_from_diff(d)
      fpkm_list[[length(fpkm_list) + 1]] <<- parsed$mat
      g1 <- colnames(parsed$mat)[1]
      g2 <- colnames(parsed$mat)[2]
      key <- paste(g2, "vs", g1, sep = "_")
      d2 <- tidy_diff_df(d)
      d2$source <- src
      diff_map[[key]] <<- d2
    }
  }

  if (!is.na(excel_path)) {
    message("Reading Excel/CSV: ", excel_path)
    sheets <- read_table_auto(excel_path)
    for (nm in names(sheets)) {
      df <- sheets[[nm]]
      if (is_cuffdiff_diff(df)) ingest_diff(df, paste0("excel:", nm))
    }
  }
  if (!is.na(diff_path)) {
    message("Reading Cuffdiff: ", diff_path)
    df <- read_table_auto(diff_path)[[1]]
    if (is_cuffdiff_diff(df)) ingest_diff(df, "gene_exp.diff")
  }
  if (!is.na(fpkm_path)) {
    message("Reading FPKM tracking: ", fpkm_path)
    fpkm_list[[length(fpkm_list) + 1]] <- read_fpkm_tracking(fpkm_path)
  }

  if (!length(fpkm_list)) {
    stop("No expression matrix found in ", data_dir,
         ". Put shTG Excel, gene_exp.diff, or genes.fpkm_tracking there.")
  }

  fpkm <- merge_fpkm(fpkm_list)
  missing_groups <- setdiff(ALL_GROUPS, colnames(fpkm))
  if (length(missing_groups)) {
    stop("Expression matrix is missing groups: ", paste(missing_groups, collapse = ", "),
         ". Current columns: ", paste(colnames(fpkm), collapse = ", "),
         ". Add genes.fpkm_tracking so TG_sh1 and TG_sh5 are both present.")
  }
  rawmap <- attr(fpkm, "gene_raw")
  fpkm <- fpkm[, ALL_GROUPS, drop = FALSE]
  fpkm[is.na(fpkm)] <- 0
  if (!is.null(rawmap)) attr(fpkm, "gene_raw") <- rawmap[rownames(fpkm)]
  list(fpkm = fpkm, diff_map = diff_map, excel_path = excel_path,
       diff_path = diff_path, fpkm_path = fpkm_path)
}

# ---------------------------------------------------------------------------
# DEG tables
# ---------------------------------------------------------------------------

match_diff_stats <- function(genes, treat, ctrl, diff_map) {
  keys <- c(
    paste(treat, "vs", ctrl, sep = "_"),
    paste(ctrl, "vs", treat, sep = "_")
  )
  hit <- intersect(keys, names(diff_map))
  out <- data.frame(
    gene = genes,
    p_value = NA_real_,
    q_value = NA_real_,
    cuffdiff_log2fc = NA_real_,
    status = NA_character_,
    stringsAsFactors = FALSE
  )
  if (!length(hit)) return(out)
  d <- diff_map[[hit[[1]]]]
  idx <- match(genes, d$gene)
  flip <- grepl(paste0("^", ctrl, "_vs_"), hit[[1]])
  out$p_value <- d$p_value[idx]
  out$q_value <- d$q_value[idx]
  out$cuffdiff_log2fc <- d$log2_fold_change[idx]
  if (flip) out$cuffdiff_log2fc <- -out$cuffdiff_log2fc
  if ("status" %in% names(d)) out$status <- d$status[idx]
  out
}

make_deg <- function(fpkm, treat_vec, ctrl_name, label, diff_map = list(),
                     treat_name = NULL) {
  ctrl <- fpkm[, ctrl_name]
  if (is.matrix(treat_vec) || is.data.frame(treat_vec)) {
    treat <- rowMeans(as.matrix(treat_vec), na.rm = TRUE)
  } else {
    treat <- as.numeric(treat_vec)
  }
  names(treat) <- rownames(fpkm)
  log2fc <- log2((treat + PSEUDOCOUNT) / (ctrl + PSEUDOCOUNT))
  fc <- 2^log2fc
  stats <- match_diff_stats(rownames(fpkm),
                            treat = treat_name %||% label,
                            ctrl = ctrl_name, diff_map = diff_map)
  rawmap <- attr(fpkm, "gene_raw")
  gene_raw <- if (!is.null(rawmap)) as.character(rawmap[rownames(fpkm)]) else rownames(fpkm)
  gene_raw[is.na(gene_raw) | !nzchar(gene_raw)] <- rownames(fpkm)
  deg <- data.frame(
    gene = rownames(fpkm),
    gene_raw = gene_raw,
    ctrl_fpkm = as.numeric(ctrl),
    treat_fpkm = as.numeric(treat),
    log2FC = as.numeric(log2fc),
    FC = as.numeric(fc),
    p_value = stats$p_value,
    q_value = stats$q_value,
    status = stats$status,
    annotated = !is_xloc_id(rownames(fpkm)),
    stringsAsFactors = FALSE
  )
  if (label == "meanSH_vs_NTC") {
    s1 <- match_diff_stats(deg$gene, "TG_sh1", "NTC", diff_map)
    s5 <- match_diff_stats(deg$gene, "TG_sh5", "NTC", diff_map)
    p1 <- pmax(s1$p_value, .Machine$double.xmin)
    p5 <- pmax(s5$p_value, .Machine$double.xmin)
    ok <- is.finite(s1$p_value) & is.finite(s5$p_value)
    fisher <- rep(NA_real_, nrow(deg))
    fisher[ok] <- stats::pchisq(-2 * (log(p1[ok]) + log(p5[ok])), df = 4, lower.tail = FALSE)
    deg$p_value <- fisher
    deg$q_value <- ifelse(is.na(deg$p_value), NA_real_, p.adjust(deg$p_value, method = P_ADJUST))
  }
  deg$neglog10p <- ifelse(is.finite(deg$p_value) & deg$p_value > 0,
                          -log10(pmax(deg$p_value, .Machine$double.xmin)), NA_real_)
  deg[order(deg$log2FC, decreasing = TRUE), ]
}

up_by_fc <- function(deg, fc_cutoff) {
  deg[!is.na(deg$FC) & deg$FC >= fc_cutoff & deg$log2FC >= 0, , drop = FALSE]
}

up_by_topn <- function(deg, n) {
  d <- deg[!is.na(deg$log2FC) & deg$log2FC >= 0, , drop = FALSE]
  d <- d[order(d$log2FC, decreasing = TRUE), , drop = FALSE]
  head(d, n)
}

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

safe_ggsave <- function(path, plot, width = 7, height = 6) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (need_pkg("ggplot2")) {
    ggplot2::ggsave(path, plot, width = width, height = height, dpi = 300)
  }
}

plot_volcano <- function(deg, highlight_genes, title, outfile,
                         fc_cutoff = 1.25) {
  if (!need_pkg("ggplot2")) return(invisible(NULL))
  df <- deg
  log2_cut <- log2(fc_cutoff)
  df$highlight <- df$gene %in% highlight_genes
  has_p <- mean(is.finite(df$neglog10p)) > 0.3
  if (!has_p) {
    df$neglog10p <- abs(df$log2FC)
    ylab <- "|log2FC| (no p-value; y shows effect size)"
  } else {
    ylab <- "-log10(p)"
  }
  df$set <- "other"
  df$set[df$log2FC >= log2_cut] <- "up"
  df$set[df$log2FC <= -log2_cut] <- "down"
  df$set[df$highlight] <- "selected_up"
  df$set <- factor(df$set, levels = c("selected_up", "up", "down", "other"))
  pal <- c(selected_up = "#B2182B", up = "#F4A582", down = "#2166AC", other = "grey75")
  lab <- df[df$highlight & !is_xloc_id(df$gene) & !grepl(",", df$gene), , drop = FALSE]
  lab <- head(lab[order(lab$log2FC, decreasing = TRUE), ], 20)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = neglog10p, color = set)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.4, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = pal, drop = FALSE) +
    ggplot2::geom_vline(xintercept = c(-log2_cut, log2_cut), linetype = 2, colour = "grey40") +
    ggplot2::labs(title = title, x = "log2 fold change", y = ylab, color = NULL) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
  if (nrow(lab) && need_pkg("ggrepel")) {
    p <- p + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = gene),
                                      size = 3, max.overlaps = 30, colour = "black")
  }
  safe_ggsave(outfile, p)
}

plot_common_scatter <- function(sh1, sh5, highlight_genes, title, outfile) {
  if (!need_pkg("ggplot2")) return(invisible(NULL))
  genes <- intersect(sh1$gene, sh5$gene)
  df <- data.frame(
    gene = genes,
    log2FC_sh1 = sh1$log2FC[match(genes, sh1$gene)],
    log2FC_sh5 = sh5$log2FC[match(genes, sh5$gene)],
    stringsAsFactors = FALSE
  )
  df$highlight <- df$gene %in% highlight_genes
  lab <- head(df[df$highlight & !is_xloc_id(df$gene) & !grepl(",", df$gene), ], 20)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC_sh1, y = log2FC_sh5, color = highlight)) +
    ggplot2::geom_point(alpha = 0.65, size = 1.4, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#B2182B"),
                                labels = c("FALSE" = "other", "TRUE" = "common up")) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    ggplot2::labs(title = title, x = "log2FC TG_sh1 vs NTC", y = "log2FC TG_sh5 vs NTC",
                  color = NULL) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
  if (nrow(lab) && need_pkg("ggrepel")) {
    p <- p + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = gene),
                                      size = 3, max.overlaps = 30, colour = "black")
  }
  safe_ggsave(outfile, p)
}

plot_heatmap <- function(fpkm, genes, outfile, title) {
  genes <- intersect(as.character(genes), rownames(fpkm))
  genes <- genes[!is_xloc_id(genes) & !grepl("[,;|/]", genes)]
  if (exists("MIN_FPKM") && is.finite(MIN_FPKM) && length(genes)) {
    mx <- apply(fpkm[genes, ALL_GROUPS, drop = FALSE], 1, max, na.rm = TRUE)
    genes <- genes[is.finite(mx) & mx >= MIN_FPKM]
  }
  if (length(genes) < 2) {
    writeLines("Fewer than 2 annotated genes passing min FPKM; heatmap skipped.",
               con = sub("\\.pdf$", ".txt", outfile))
    return(invisible(NULL))
  }
  if (length(genes) > HEATMAP_MAX_GENES) {
    fc <- log2((rowMeans(fpkm[genes, TREAT_GROUPS, drop = FALSE]) + PSEUDOCOUNT) /
                 (fpkm[genes, CTRL_GROUP] + PSEUDOCOUNT))
    genes <- names(sort(fc, decreasing = TRUE))[seq_len(HEATMAP_MAX_GENES)]
  }
  title <- sprintf("%s | heatmap %d annotated genes", title, length(genes))
  mat <- log2(fpkm[genes, ALL_GROUPS, drop = FALSE] + 1)
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  if (need_pkg("pheatmap")) {
    grDevices::pdf(outfile, width = 6, height = max(5, length(genes) * 0.12 + 2))
    pheatmap::pheatmap(mat, scale = "row", cluster_cols = FALSE,
                       main = title, border_color = NA,
                       fontsize_row = ifelse(length(genes) > 60, 4, 7))
    grDevices::dev.off()
  } else {
    utils::write.csv(mat, sub("\\.pdf$", ".csv", outfile))
  }
}

# ---------------------------------------------------------------------------
# ID mapping and enrichment
# ---------------------------------------------------------------------------

symbol_to_entrez <- function(symbols) {
  symbols <- flatten_symbols(symbols)
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  if (!length(symbols) || !need_pkg("clusterProfiler") || !need_pkg("org.Hs.eg.db")) {
    return(data.frame(SYMBOL = character(), ENTREZID = character(), stringsAsFactors = FALSE))
  }
  mapped <- tryCatch(
    clusterProfiler::bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db::org.Hs.eg.db),
    error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
  )
  mapped[!duplicated(mapped$SYMBOL), , drop = FALSE]
}

expand_deg_for_rank <- function(deg) {
  src <- if ("gene_raw" %in% names(deg)) deg$gene_raw else deg$gene
  src[is.na(src) | !nzchar(src)] <- deg$gene[is.na(src) | !nzchar(src)]
  pieces <- vector("list", nrow(deg))
  for (i in seq_len(nrow(deg))) {
    toks <- split_gene_tokens(src[i])
    toks <- toks[!is_xloc_id(toks)]
    if (!length(toks)) next
    pieces[[i]] <- data.frame(gene = toks, log2FC = deg$log2FC[i], stringsAsFactors = FALSE)
  }
  x <- do.call(rbind, pieces)
  if (is.null(x) || !nrow(x)) {
    return(data.frame(gene = character(), log2FC = numeric(), stringsAsFactors = FALSE))
  }
  stats::aggregate(log2FC ~ gene, data = x, FUN = mean)
}

write_empty_enrich <- function(path, reason) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(data.frame(reason = reason), path, row.names = FALSE)
}

save_plot_file <- function(path, plot, width = 8, height = 6) {
  safe_ggsave(path, plot, width = width, height = height)
  png <- sub("\\.pdf$", ".png", path)
  if (need_pkg("ggplot2") && grepl("\\.pdf$", path)) {
    try(ggplot2::ggsave(png, plot, width = width, height = height, dpi = 150), silent = TRUE)
  }
}

save_enrich_result <- function(res, prefix, err = NULL) {
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  if (!is.null(err)) {
    write_empty_enrich(paste0(prefix, ".csv"), paste("error:", err))
    return(invisible(NULL))
  }
  empty <- is.null(res) ||
    (is.data.frame(res) && !nrow(res)) ||
    (inherits(res, "enrichResult") && !nrow(as.data.frame(res))) ||
    (inherits(res, "gseaResult") && !nrow(as.data.frame(res)))
  if (empty) {
    write_empty_enrich(paste0(prefix, ".csv"), "no_terms_returned")
    return(invisible(NULL))
  }
  df <- as.data.frame(res)
  utils::write.csv(df, paste0(prefix, "_all.csv"), row.names = FALSE)
  if ("p.adjust" %in% names(df)) {
    sig <- df[!is.na(df$p.adjust) & df$p.adjust < ENRICH_Q_CUTOFF, , drop = FALSE]
    utils::write.csv(sig, paste0(prefix, "_significant.csv"), row.names = FALSE)
  }
  if (need_pkg("enrichplot") && need_pkg("ggplot2") && nrow(df)) {
    n <- min(20, nrow(df))
    tryCatch({
      save_plot_file(paste0(prefix, "_dotplot.pdf"),
                     enrichplot::dotplot(res, showCategory = n) + ggplot2::ggtitle(basename(prefix)),
                     width = 9, height = 7)
    }, error = function(e) writeLines(conditionMessage(e), paste0(prefix, "_dotplot_error.txt")))
    tryCatch({
      save_plot_file(paste0(prefix, "_barplot.pdf"),
                     enrichplot::barplot(res, showCategory = n) + ggplot2::ggtitle(basename(prefix)),
                     width = 9, height = 7)
    }, error = function(e) writeLines(conditionMessage(e), paste0(prefix, "_barplot_error.txt")))
  }
}

run_ora <- function(genes, universe_symbols, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (length(genes) < ORA_MIN_GENES) {
    write_empty_enrich(file.path(out_dir, "ORA_skipped.csv"),
                       paste("only", length(genes), "genes; need", ORA_MIN_GENES))
    return(invisible(NULL))
  }
  if (!need_pkg("clusterProfiler") || !need_pkg("org.Hs.eg.db")) {
    write_empty_enrich(file.path(out_dir, "ORA_skipped.csv"),
                       "clusterProfiler/org.Hs.eg.db not installed")
    return(invisible(NULL))
  }
  map <- symbol_to_entrez(genes)
  univ <- symbol_to_entrez(universe_symbols)
  if (nrow(map) < ORA_MIN_GENES) {
    write_empty_enrich(file.path(out_dir, "ORA_skipped.csv"),
                       paste("only", nrow(map), "mapped Entrez IDs"))
    return(invisible(NULL))
  }
  utils::write.csv(map, file.path(out_dir, "gene_id_map.csv"), row.names = FALSE)
  entrez <- unique(map$ENTREZID)
  universe <- unique(univ$ENTREZID)
  for (ont in c("BP", "MF", "CC")) {
    res <- tryCatch(
      clusterProfiler::enrichGO(
        gene = entrez, OrgDb = org.Hs.eg.db::org.Hs.eg.db, ont = ont,
        keyType = "ENTREZID", pAdjustMethod = P_ADJUST,
        pvalueCutoff = ENRICH_P_CUTOFF, qvalueCutoff = ENRICH_Q_CUTOFF,
        universe = if (length(universe) > length(entrez)) universe else NULL,
        readable = TRUE
      ),
      error = function(e) NULL
    )
    save_enrich_result(res, file.path(out_dir, paste0("GO_", ont)))
  }
  kegg <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene = entrez, organism = KEGG_ORG, pvalueCutoff = ENRICH_P_CUTOFF,
      qvalueCutoff = ENRICH_Q_CUTOFF,
      universe = if (length(universe) > length(entrez)) universe else NULL
    ),
    error = function(e) e
  )
  if (inherits(kegg, "error")) {
    save_enrich_result(NULL, file.path(out_dir, "KEGG"), err = conditionMessage(kegg))
  } else {
    if (!is.null(kegg) && need_pkg("clusterProfiler")) {
      kegg <- tryCatch(clusterProfiler::setReadable(kegg, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                                    keyType = "ENTREZID"), error = function(e) kegg)
    }
    save_enrich_result(kegg, file.path(out_dir, "KEGG"))
  }
  if (need_pkg("ReactomePA")) {
    rea <- tryCatch(
      ReactomePA::enrichPathway(
        gene = entrez, organism = "human", pvalueCutoff = ENRICH_P_CUTOFF,
        qvalueCutoff = ENRICH_Q_CUTOFF, readable = TRUE,
        universe = if (length(universe) > length(entrez)) universe else NULL
      ),
      error = function(e) NULL
    )
    save_enrich_result(rea, file.path(out_dir, "Reactome"))
  } else {
    write_empty_enrich(file.path(out_dir, "Reactome.csv"), "ReactomePA not installed")
  }
}

run_gsea <- function(deg, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ranked <- expand_deg_for_rank(deg)
  diag <- data.frame(
    n_input_rows = nrow(deg),
    n_xloc_input = sum(is_xloc_id(deg$gene)),
    n_symbols_after_split = nrow(ranked),
    n_mapped = NA_integer_,
    n_ranked = NA_integer_,
    stringsAsFactors = FALSE
  )
  if (!need_pkg("clusterProfiler") || !need_pkg("org.Hs.eg.db")) {
    diag$note <- "clusterProfiler/org.Hs.eg.db not installed"
    utils::write.csv(diag, file.path(out_dir, "GSEA_diagnostics.csv"), row.names = FALSE)
    write_empty_enrich(file.path(out_dir, "GSEA_skipped.csv"),
                       "clusterProfiler/org.Hs.eg.db not installed")
    return(invisible(NULL))
  }
  map <- symbol_to_entrez(ranked$gene)
  diag$n_mapped <- nrow(map)
  if (nrow(map) < GSEA_MIN_GENES) {
    diag$note <- "too few mapped Entrez IDs (XLOC/comma names dropped)"
    utils::write.csv(diag, file.path(out_dir, "GSEA_diagnostics.csv"), row.names = FALSE)
    write_empty_enrich(file.path(out_dir, "GSEA_skipped.csv"),
                       paste("only", nrow(map), "mapped Entrez IDs"))
    return(invisible(NULL))
  }
  utils::write.csv(map, file.path(out_dir, "gsea_gene_id_map.csv"), row.names = FALSE)
  stats <- ranked$log2FC[match(map$SYMBOL, ranked$gene)]
  names(stats) <- map$ENTREZID
  stats <- stats[is.finite(stats)]
  stats <- tapply(stats, names(stats), mean)
  stats <- sort(stats, decreasing = TRUE)
  diag$n_ranked <- length(stats)
  utils::write.csv(diag, file.path(out_dir, "GSEA_diagnostics.csv"), row.names = FALSE)
  if (length(stats) < GSEA_MIN_GENES) {
    write_empty_enrich(file.path(out_dir, "GSEA_skipped.csv"), "too few finite ranks")
    return(invisible(NULL))
  }

  save_gsea <- function(res_or_err, prefix) {
    if (inherits(res_or_err, "error") || inherits(res_or_err, "simpleError")) {
      save_enrich_result(NULL, prefix, err = conditionMessage(res_or_err))
      return(invisible(NULL))
    }
    save_enrich_result(res_or_err, prefix)
    if (is.null(res_or_err)) return(invisible(NULL))
    df <- tryCatch(as.data.frame(res_or_err), error = function(e) data.frame())
    if (!nrow(df) || !need_pkg("enrichplot") || !need_pkg("ggplot2")) return(invisible(NULL))
    nshow <- min(8, nrow(df))
    for (i in seq_len(nshow)) {
      tryCatch({
        p <- enrichplot::gseaplot2(res_or_err, geneSetID = i, title = df$Description[i])
        save_plot_file(sprintf("%s_gseaplot_%02d.pdf", prefix, i), p, width = 8, height = 6)
      }, error = function(e) {
        writeLines(conditionMessage(e), sprintf("%s_gseaplot_%02d_error.txt", prefix, i))
      })
    }
  }

  gsea_cut <- if (exists("GSEA_PVALUE_CUTOFF")) GSEA_PVALUE_CUTOFF else 1
  for (ont in c("BP", "MF", "CC")) {
    res <- tryCatch(
      clusterProfiler::gseGO(
        geneList = stats, OrgDb = org.Hs.eg.db::org.Hs.eg.db, ont = ont,
        keyType = "ENTREZID", minGSSize = 10, maxGSSize = 500,
        pvalueCutoff = gsea_cut, verbose = FALSE, eps = 0
      ),
      error = function(e) e
    )
    save_gsea(res, file.path(out_dir, paste0("GSEA_GO_", ont)))
  }
  kegg <- tryCatch(
    clusterProfiler::gseKEGG(
      geneList = stats, organism = KEGG_ORG, minGSSize = 10, maxGSSize = 500,
      pvalueCutoff = gsea_cut, verbose = FALSE, eps = 0
    ),
    error = function(e) e
  )
  save_gsea(kegg, file.path(out_dir, "GSEA_KEGG"))
}

# ---------------------------------------------------------------------------
# One filter / one comparison
# ---------------------------------------------------------------------------

write_deg_table <- function(deg, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(deg, path, row.names = FALSE)
}

analyze_gene_set <- function(cmp_label, filter_label, deg_all, selected, fpkm,
                             out_dir, volcano_mode = "standard",
                             sh1 = NULL, sh5 = NULL, fc_cutoff = 1.25) {
  message("  ", cmp_label, " / ", filter_label, " : ", nrow(selected), " up genes")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_deg_table(selected, file.path(out_dir, "up_genes.csv"))
  writeLines(selected$gene, file.path(out_dir, "up_genes.txt"))
  title <- paste0(cmp_label, " | ", filter_label, " | n=", nrow(selected))
  if (identical(volcano_mode, "common_scatter")) {
    plot_common_scatter(sh1, sh5, selected$gene, title,
                        file.path(out_dir, "scatter_sh1_vs_sh5.pdf"))
  } else {
    plot_volcano(deg_all, selected$gene, title,
                 file.path(out_dir, "volcano.pdf"), fc_cutoff = fc_cutoff)
  }
  plot_heatmap(fpkm, selected$gene, file.path(out_dir, "heatmap.pdf"), title)
  ora_genes <- if ("gene_raw" %in% names(selected)) selected$gene_raw else selected$gene
  ora_univ <- if ("gene_raw" %in% names(deg_all)) deg_all$gene_raw else deg_all$gene
  run_ora(ora_genes, ora_univ, file.path(out_dir, "enrichment"))
}

run_two_strategies <- function(cmp_label, deg_all, fpkm, out_root,
                               volcano_mode = "standard", sh1 = NULL, sh5 = NULL,
                               deg_sh1 = NULL, deg_sh5 = NULL) {
  cmp_dir <- file.path(out_root, cmp_label)
  dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)
  write_deg_table(deg_all, file.path(cmp_dir, "all_genes_ranked.csv"))

  if (!identical(volcano_mode, "common_scatter")) {
    run_gsea(deg_all, file.path(cmp_dir, "GSEA_full_ranked_list"))
  } else {
    # common-up track: rank by mean log2FC of the two knockdowns
    common_rank <- deg_all
    run_gsea(common_rank, file.path(cmp_dir, "GSEA_full_ranked_list"))
  }

  for (fc in FC_CUTOFFS) {
    lab <- sprintf("FC_ge_%s", format(fc, nsmall = 2, trim = TRUE))
    if (identical(volcano_mode, "common_scatter")) {
      g1 <- up_by_fc(deg_sh1, fc)$gene
      g5 <- up_by_fc(deg_sh5, fc)$gene
      sel_genes <- intersect(g1, g5)
      selected <- deg_all[deg_all$gene %in% sel_genes, , drop = FALSE]
    } else {
      selected <- up_by_fc(deg_all, fc)
    }
    analyze_gene_set(cmp_label, lab, deg_all, selected, fpkm,
                     file.path(cmp_dir, lab), volcano_mode, sh1, sh5, fc_cutoff = fc)
  }
  for (n in TOP_N) {
    lab <- sprintf("top_%d", n)
    if (identical(volcano_mode, "common_scatter")) {
      g1 <- up_by_topn(deg_sh1, n)$gene
      g5 <- up_by_topn(deg_sh5, n)$gene
      sel_genes <- intersect(g1, g5)
      selected <- deg_all[match(sel_genes, deg_all$gene), , drop = FALSE]
      selected <- selected[!is.na(selected$gene), , drop = FALSE]
    } else {
      selected <- up_by_topn(deg_all, n)
    }
    analyze_gene_set(cmp_label, lab, deg_all, selected, fpkm,
                     file.path(cmp_dir, lab), volcano_mode, sh1, sh5,
                     fc_cutoff = 1)
  }
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function(data_dir = NULL, out_dir = NULL) {
  if (is.null(data_dir) || !nzchar(data_dir)) data_dir <- default_data_dir()
  if (is.null(out_dir) || !nzchar(out_dir)) out_dir <- default_out_dir(data_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("Data dir: ", data_dir)
  message("Out dir:  ", out_dir)

  optional_library("ggplot2")
  optional_library("ggrepel")
  optional_library("pheatmap")
  optional_library("clusterProfiler")
  optional_library("org.Hs.eg.db")
  optional_library("enrichplot")
  optional_library("ReactomePA")

  inp <- load_inputs(data_dir)
  fpkm <- inp$fpkm
  diff_map <- inp$diff_map
  utils::write.csv(as.data.frame(fpkm), file.path(out_dir, "expression_fpkm_matrix.csv"))
  message("FPKM groups: ", paste(colnames(fpkm), collapse = ", "),
          " ; genes: ", nrow(fpkm))

  deg_sh1 <- make_deg(fpkm, fpkm[, "TG_sh1"], "NTC", "sh1_vs_NTC", diff_map,
                      treat_name = "TG_sh1")
  deg_sh5 <- make_deg(fpkm, fpkm[, "TG_sh5"], "NTC", "sh5_vs_NTC", diff_map,
                      treat_name = "TG_sh5")
  deg_mean <- make_deg(fpkm, fpkm[, TREAT_GROUPS, drop = FALSE], "NTC",
                       "meanSH_vs_NTC", diff_map, treat_name = "meanSH")

  genes <- intersect(deg_sh1$gene, deg_sh5$gene)
  deg_common <- data.frame(
    gene = genes,
    gene_raw = deg_sh1$gene_raw[match(genes, deg_sh1$gene)],
    ctrl_fpkm = deg_sh1$ctrl_fpkm[match(genes, deg_sh1$gene)],
    treat_fpkm = rowMeans(cbind(
      deg_sh1$treat_fpkm[match(genes, deg_sh1$gene)],
      deg_sh5$treat_fpkm[match(genes, deg_sh5$gene)]
    )),
    log2FC_sh1 = deg_sh1$log2FC[match(genes, deg_sh1$gene)],
    log2FC_sh5 = deg_sh5$log2FC[match(genes, deg_sh5$gene)],
    stringsAsFactors = FALSE
  )
  deg_common$log2FC <- rowMeans(cbind(deg_common$log2FC_sh1, deg_common$log2FC_sh5))
  deg_common$FC <- 2^deg_common$log2FC
  deg_common$p_value <- NA_real_
  deg_common$q_value <- NA_real_
  deg_common$neglog10p <- abs(deg_common$log2FC)
  deg_common$annotated <- !is_xloc_id(deg_common$gene)
  deg_common <- deg_common[order(deg_common$log2FC, decreasing = TRUE), ]

  run_two_strategies("sh1_vs_NTC", deg_sh1, fpkm, out_dir)
  run_two_strategies("sh5_vs_NTC", deg_sh5, fpkm, out_dir)
  run_two_strategies("meanSH_vs_NTC", deg_mean, fpkm, out_dir)
  run_two_strategies("common_up_sh1_sh5", deg_common, fpkm, out_dir,
                     volcano_mode = "common_scatter",
                     sh1 = deg_sh1, sh5 = deg_sh5,
                     deg_sh1 = deg_sh1, deg_sh5 = deg_sh5)

  common_n <- function(fc) {
    length(intersect(up_by_fc(deg_sh1, fc)$gene, up_by_fc(deg_sh5, fc)$gene))
  }
  summary <- data.frame(
    comparison = c("sh1_vs_NTC", "sh5_vs_NTC", "meanSH_vs_NTC", "common_up_sh1_sh5"),
    n_genes = c(nrow(deg_sh1), nrow(deg_sh5), nrow(deg_mean), nrow(deg_common)),
    n_up_FC1 = c(nrow(up_by_fc(deg_sh1, 1)), nrow(up_by_fc(deg_sh5, 1)),
                 nrow(up_by_fc(deg_mean, 1)), common_n(1)),
    n_up_FC1.25 = c(nrow(up_by_fc(deg_sh1, 1.25)), nrow(up_by_fc(deg_sh5, 1.25)),
                    nrow(up_by_fc(deg_mean, 1.25)), common_n(1.25)),
    n_up_FC1.5 = c(nrow(up_by_fc(deg_sh1, 1.5)), nrow(up_by_fc(deg_sh5, 1.5)),
                   nrow(up_by_fc(deg_mean, 1.5)), common_n(1.5)),
    n_up_FC2 = c(nrow(up_by_fc(deg_sh1, 2)), nrow(up_by_fc(deg_sh5, 2)),
                 nrow(up_by_fc(deg_mean, 2)), common_n(2)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(summary, file.path(out_dir, "summary_up_counts.csv"), row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  message("Done. Results written to ", out_dir)
}

# Run in both Rscript and interactive paste/source(). Set TG_RNASEQ_NO_AUTO <- TRUE
# before sourcing if you only want to load functions.
if (!exists("TG_RNASEQ_NO_AUTO") || isFALSE(get("TG_RNASEQ_NO_AUTO"))) {
  args <- commandArgs(trailingOnly = TRUE)
  data_dir <- if (length(args) >= 1) args[[1]] else NULL
  out_dir <- if (length(args) >= 2) args[[2]] else NULL
  main(data_dir, out_dir)
}
