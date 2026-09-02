# =============================================================================
# 去除各类免疫球蛋白后的丰度排名气泡图（独立脚本，不必先 source 原脚本）
#
#   setwd("E:/天府/实验管理/课题/赵章寻/血清蛋白质组学")
#   source("serum_proteomics_ranked_bubble_no_Ig.R", encoding = "UTF-8")
#
# 结果：results/serum_proteomics_bubble_no_Ig/
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

cran_required <- c("dplyr", "tidyr", "tibble", "stringr", "ggplot2", "readr")
miss <- cran_required[!vapply(cran_required, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) install.packages(miss, repos = "https://cloud.r-project.org")
for (p in cran_required) suppressPackageStartupMessages(library(p, character.only = TRUE))

default_proteomics_dirs <- c(
  "E:/天府/实验管理/课题/赵章寻/血清蛋白质组学",
  "E:\\天府\\实验管理\\课题\\赵章寻\\血清蛋白质组学"
)
top_ns <- c(20L, 30L, 50L)
min_detect_frac <- 0.5
meta_cols <- c(
  "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
  "First.Protein.Description", "N.Proteins", "Proteotypic",
  "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
  "Precursor.Id", "Precursor.Lib.Index", "Q.Value"
)

resolve_proteomics_dir <- function() {
  env_dir <- Sys.getenv("SERUM_PROTEOMICS_DIR", unset = "")
  script_dir <- tryCatch({
    ofile <- sys.frame(1)$ofile
    if (!is.null(ofile) && nzchar(ofile)) dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)) else NA_character_
  }, error = function(e) NA_character_)
  candidates <- unique(c(env_dir, default_proteomics_dirs, script_dir, getwd()))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  for (d in candidates) {
    if (!dir.exists(d)) next
    if (file.exists(file.path(d, "report.pg_matrix")) ||
        file.exists(file.path(d, "report.pr_matrix")) ||
        file.exists(file.path(d, "report.pg_matrix.tsv"))) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

project_dir <- resolve_proteomics_dir()
result_dir <- file.path(project_dir, "results", "serum_proteomics_bubble_no_Ig")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(result_dir, paste0("log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

pick_official_symbol <- function(x) {
  x <- trimws(as.character(x))
  if (length(x) != 1 || is.na(x) || x %in% c("", "-", ".", "NA")) return(NA_character_)
  parts <- trimws(unlist(strsplit(x, "[,;|/]+")))
  parts <- parts[nzchar(parts) & !parts %in% c("-", ".", "NA")]
  if (length(parts) == 0) return(NA_character_)
  score <- vapply(parts, function(s) {
    if (grepl("^(XLOC|TCONS|CUFF)_", s, ignore.case = TRUE)) return(0)
    if (grepl("^LOC[0-9]+$", s, ignore.case = TRUE)) return(1)
    2
  }, numeric(1))
  parts[which.max(score)]
}

is_immunoglobulin_symbol <- function(sym) {
  s <- toupper(trimws(as.character(sym)))
  if (!nzchar(s) || s %in% c("NA", "-", ".")) return(FALSE)
  if (s %in% c("JCHAIN", "IGJ")) return(TRUE)
  grepl("^(IGH|IGK|IGL)", s)
}

is_immunoglobulin_text <- function(...) {
  blob <- tolower(paste(..., collapse = " "))
  if (!nzchar(trimws(blob))) return(FALSE)
  if (grepl("immunoglobulin superfamily", blob) &&
      !grepl("\\bigh[agmdvejk]|\\bigk[cvlj]|\\bigl[cvlj]", blob)) {
    return(FALSE)
  }
  grepl("immunoglobulin\\s+(heavy|kappa|lambda|alpha|gamma|mu|delta|epsilon)", blob) ||
    grepl("\\big\\s*(heavy|kappa|lambda|gamma|alpha|mu)\\b", blob) ||
    grepl("\\big\\s+(gamma|alpha|mu|kappa|lambda)\\b", blob)
}

immunoglobulin_mask <- function(meta) {
  n <- nrow(meta)
  if (n == 0) return(logical(0))
  gene_raw <- if ("Genes" %in% names(meta)) as.character(meta$Genes) else rep("", n)
  names_ <- if ("Protein.Names" %in% names(meta)) as.character(meta$Protein.Names) else rep("", n)
  desc <- if ("First.Protein.Description" %in% names(meta)) as.character(meta$First.Protein.Description) else rep("", n)
  ids <- if ("Protein.Ids" %in% names(meta)) as.character(meta$Protein.Ids) else rep("", n)
  pg <- if ("Protein.Group" %in% names(meta)) as.character(meta$Protein.Group) else rep("", n)
  vapply(seq_len(n), function(i) {
    parts <- trimws(unlist(strsplit(paste(c(gene_raw[i], names_[i], pg[i]), collapse = ","), "[,;|/ ]+")))
    parts <- parts[nzchar(parts)]
    any(vapply(parts, is_immunoglobulin_symbol, logical(1))) ||
      is_immunoglobulin_text(gene_raw[i], names_[i], desc[i], ids[i], pg[i])
  }, logical(1))
}

find_matrix_file <- function(dir, stems) {
  for (s in stems) {
    for (name in c(s, paste0(s, ".tsv"), paste0(s, ".txt"))) {
      p <- file.path(dir, name)
      if (file.exists(p)) return(p)
    }
  }
  NA_character_
}

is_intensity_col <- function(nm, df) {
  if (nm %in% meta_cols) return(FALSE)
  if (grepl("^(Q\\.Value|PG\\.Q|GG\\.Q|Protein\\.Q|Lib\\.|Ms1\\.|Normalisation)", nm)) return(FALSE)
  x <- df[[nm]]
  if (is.numeric(x) || is.integer(x)) return(TRUE)
  if (is.character(x) || is.factor(x)) {
    num <- suppressWarnings(as.numeric(as.character(x)))
    return(mean(is.finite(num)) >= 0.3)
  }
  FALSE
}

read_diann_matrix <- function(path) {
  log_msg("读取矩阵: ", path)
  df <- as.data.frame(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE), stringsAsFactors = FALSE)
  int_cols <- names(df)[vapply(names(df), is_intensity_col, logical(1), df = df)]
  if (length(int_cols) == 0) stop("没有找到样品强度列: ", path)
  for (col in int_cols) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  attr(df, "intensity_cols") <- int_cols
  df
}

aggregate_pr_to_pg <- function(pr) {
  int_cols <- attr(pr, "intensity_cols")
  key <- if ("Protein.Group" %in% names(pr)) "Protein.Group" else names(pr)[1]
  log_msg("pr_matrix 按 ", key, " 取肽段强度中位数，聚合到蛋白")
  num <- pr[, c(key, int_cols), drop = FALSE]
  agg <- num |>
    dplyr::group_by(.data[[key]]) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(int_cols), ~ median(.x, na.rm = TRUE)), .groups = "drop")
  extra_keys <- intersect(c("Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description"), names(pr))
  if (length(extra_keys) > 0) {
    meta <- pr[, c(key, extra_keys), drop = FALSE] |>
      dplyr::group_by(.data[[key]]) |>
      dplyr::summarise(dplyr::across(dplyr::everything(), ~ {
        u <- unique(as.character(.x))
        u <- u[!is.na(u) & nzchar(u)]
        if (length(u) == 0) NA_character_ else u[1]
      }), .groups = "drop")
    agg <- dplyr::left_join(agg, meta, by = key)
  }
  agg <- as.data.frame(agg)
  attr(agg, "intensity_cols") <- int_cols
  agg
}

load_protein_matrix <- function(dir) {
  pg_path <- find_matrix_file(dir, c("report.pg_matrix"))
  pr_path <- find_matrix_file(dir, c("report.pr_matrix"))
  pg <- NULL
  if (!is.na(pg_path)) {
    pg <- read_diann_matrix(pg_path)
    log_msg("pg_matrix: ", nrow(pg), " 行, ", length(attr(pg, "intensity_cols")), " 个样品列")
    if (nrow(pg) < 5) {
      log_msg("pg_matrix 蛋白行过少，回退 pr_matrix")
      pg <- NULL
    }
  }
  if (!is.null(pg)) return(pg)
  if (is.na(pr_path)) stop("找不到 report.pg_matrix 或 report.pr_matrix，目录: ", dir)
  aggregate_pr_to_pg(read_diann_matrix(pr_path))
}

normalize_sample_token <- function(x) {
  x <- basename(gsub("\\\\", "/", as.character(x)))
  sub("\\.(raw|wiff|mzML|d)$", "", x, ignore.case = TRUE)
}

preprocess <- function(mat, int_cols) {
  raw <- as.matrix(mat[, int_cols, drop = FALSE])
  storage.mode(raw) <- "double"
  raw[!is.finite(raw) | raw < 0] <- NA_real_
  keep <- rowMeans(is.finite(raw) & raw > 0) >= min_detect_frac
  log_msg("过滤前 ", nrow(raw), " 蛋白，过滤后 ", sum(keep), "；样品: ", paste(int_cols, collapse = ", "))
  mat <- mat[keep, , drop = FALSE]
  raw <- raw[keep, , drop = FALSE]
  log2x <- log2(raw + 1)
  med <- apply(log2x, 2, median, na.rm = TRUE)
  log2x <- sweep(log2x, 2, med - median(med, na.rm = TRUE), "-")
  list(meta = mat, log2 = log2x, sample_cols = int_cols)
}

rank_proteins <- function(meta, log2x, sample_cols) {
  mean_ab <- rowMeans(log2x, na.rm = TRUE)
  gene_raw <- if ("Genes" %in% names(meta)) meta$Genes else NA_character_
  pg <- if ("Protein.Group" %in% names(meta)) meta$Protein.Group else seq_len(nrow(meta))
  gene <- vapply(as.character(gene_raw), pick_official_symbol, character(1), USE.NAMES = FALSE)
  gene[is.na(gene) | !nzchar(gene)] <- as.character(pg[is.na(gene) | !nzchar(gene)])
  gene[duplicated(gene)] <- paste0(gene[duplicated(gene)], "_", which(duplicated(gene)))
  out <- data.frame(
    Protein.Group = as.character(pg),
    gene = gene,
    description = if ("First.Protein.Description" %in% names(meta)) as.character(meta$First.Protein.Description) else NA_character_,
    mean_abundance = as.numeric(mean_ab),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(sample_cols)) {
    out[[paste0("sample_", normalize_sample_token(sample_cols[[i]]))]] <- as.numeric(log2x[, i])
  }
  out$abundance_rank <- rank(-out$mean_abundance, ties.method = "first")
  out[order(out$abundance_rank), ]
}

plot_rank_abundance_bubble <- function(df, n, outfile, label_n = 12L) {
  sub <- df[df$abundance_rank <= n, ]
  sub <- sub[order(sub$abundance_rank), ]
  lab <- sub[sub$abundance_rank <= min(as.integer(label_n), n), ]
  p <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$abundance_rank, y = .data$mean_abundance)) +
    ggplot2::geom_line(color = "#9ECAE1", size = 0.45) +
    ggplot2::geom_point(size = 2.8, color = "#2C7FB8", alpha = 0.92) +
    ggplot2::labs(
      title = paste0("去除免疫球蛋白后 血清蛋白丰度排名气泡图 (top", n, ")"),
      subtitle = "去掉 IGH/IGK/IGL/JCHAIN 等后，两样品平均丰度排名；气泡大小一致",
      x = "丰度排名",
      y = "蛋白丰度值"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold", hjust = 0))
  if (nrow(lab) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = lab, ggplot2::aes(label = .data$gene),
        size = 3, fontface = "italic", max.overlaps = 40, min.segment.length = 0
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = lab, ggplot2::aes(label = .data$gene),
        vjust = -0.9, size = 3, fontface = "italic", check_overlap = TRUE
      )
    }
  }
  ggplot2::ggsave(outfile, p, width = 8.2, height = 5.4, dpi = 150)
}

log_msg("工作目录: ", project_dir)
mat <- load_protein_matrix(project_dir)
int_cols <- attr(mat, "intensity_cols")
prep <- preprocess(mat, int_cols)
ig <- immunoglobulin_mask(prep$meta)
dropped <- prep$meta[ig, , drop = FALSE]
drop_path <- file.path(result_dir, "removed_immunoglobulins.csv")
if (nrow(dropped) > 0) {
  keep_cols <- intersect(
    c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description"),
    names(dropped)
  )
  readr::write_csv(dropped[, keep_cols, drop = FALSE], drop_path)
}
log_msg("去除免疫球蛋白 ", sum(ig), " / ", length(ig), " -> ", drop_path)
prep$meta <- prep$meta[!ig, , drop = FALSE]
prep$log2 <- prep$log2[!ig, , drop = FALSE]
ranked <- rank_proteins(prep$meta, prep$log2, prep$sample_cols)
log_msg("去 Ig 后按两样品平均丰度排名，蛋白数 ", nrow(ranked))
readr::write_csv(ranked, file.path(result_dir, "protein_abundance_ranking.csv"))
ns <- unique(c(top_ns[top_ns <= nrow(ranked)], nrow(ranked)))
for (n_use in ns) {
  if (n_use < 1) next
  tag <- if (n_use == nrow(ranked)) "all" else paste0("top", n_use)
  readr::write_csv(ranked[ranked$abundance_rank <= n_use, ], file.path(result_dir, paste0(tag, "_ranked_proteins.csv")))
  plot_rank_abundance_bubble(ranked, n_use, file.path(result_dir, paste0(tag, "_abundance_rank_bubble.pdf")))
  log_msg("完成 ", tag)
}
log_msg("全部完成 -> ", result_dir)
