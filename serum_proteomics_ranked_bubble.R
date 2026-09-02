#!/usr/bin/env Rscript
# =============================================================================
# 两组病人血清蛋白质组：按蛋白丰度做排名气泡图
# 输入：DIA-NN report.pg_matrix（首选）或 report.pr_matrix
# 分组：sample_annotation.csv（列 sample,group，恰好两组）
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

cran_required <- c("dplyr", "tidyr", "tibble", "stringr", "ggplot2", "readr")

install_if_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) == 0) return(invisible(TRUE))
  install.packages(miss, repos = "https://cloud.r-project.org")
  still <- miss[!vapply(miss, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still) > 0) stop("缺少必需 R 包: ", paste(still, collapse = ", "))
}

install_if_missing(cran_required)
for (p in cran_required) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# -----------------------------------------------------------------------------
# 路径
# -----------------------------------------------------------------------------
resolve_proteomics_dir <- function() {
  env_dir <- Sys.getenv("SERUM_PROTEOMICS_DIR", unset = "")
  candidates <- unique(c(
    env_dir,
    "E:/R/TG_BRCA/TG",
    "E:\\R\\TG_BRCA\\TG",
    file.path(getwd(), "serum_proteomics"),
    getwd()
  ))
  candidates <- candidates[nzchar(candidates)]
  for (d in candidates) {
    if (!dir.exists(d)) next
    if (file.exists(file.path(d, "report.pg_matrix")) ||
        file.exists(file.path(d, "report.pr_matrix")) ||
        file.exists(file.path(d, "report.pg_matrix.tsv")) ||
        file.exists(file.path(d, "sample_annotation.csv"))) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

project_dir <- resolve_proteomics_dir()
result_dir <- file.path(project_dir, "results", "serum_proteomics_bubble")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(result_dir, paste0("log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

top_ns <- c(20L, 30L, 50L)
min_detect_frac <- 0.5

meta_cols <- c(
  "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
  "First.Protein.Description", "N.Proteins", "Proteotypic",
  "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
  "Precursor.Id", "Precursor.Lib.Index", "Q.Value"
)

# -----------------------------------------------------------------------------
# 基因名清洗（复合名取一个官方符号）
# -----------------------------------------------------------------------------
pick_official_symbol <- function(x) {
  x <- trimws(as.character(x))
  if (length(x) != 1 || is.na(x) || x %in% c("", "-", ".", "NA")) return(NA_character_)
  parts <- unlist(strsplit(x, "[,;|/]+"))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts) & !parts %in% c("-", ".", "NA")]
  if (length(parts) == 0) return(NA_character_)
  score <- vapply(parts, function(s) {
    if (grepl("^(XLOC|TCONS|CUFF)_", s, ignore.case = TRUE)) return(0)
    if (grepl("^LOC[0-9]+$", s, ignore.case = TRUE)) return(1)
    2
  }, numeric(1))
  parts[which.max(score)]
}

# -----------------------------------------------------------------------------
# 读 DIA-NN 矩阵
# -----------------------------------------------------------------------------
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
  if (grepl("^(Q\\.Value|PG\\.Q|GG\\.Q|Protein\\.Q|Lib\\.|Ms1\\.|Normalisation)", nm)) {
    return(FALSE)
  }
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
  df <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  int_cols <- names(df)[vapply(names(df), is_intensity_col, logical(1), df = df)]
  if (length(int_cols) == 0) stop("没有找到样品强度列: ", path)
  for (col in int_cols) {
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  }
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
    n_prot <- nrow(pg)
    n_int <- length(attr(pg, "intensity_cols"))
    log_msg("pg_matrix: ", n_prot, " 行, ", n_int, " 个样品列")
    if (n_prot < 5) {
      log_msg("pg_matrix 蛋白行过少，回退 pr_matrix")
      pg <- NULL
    }
  }
  if (!is.null(pg)) return(pg)
  if (is.na(pr_path)) stop("找不到 report.pg_matrix 或 report.pr_matrix，目录: ", dir)
  pr <- read_diann_matrix(pr_path)
  aggregate_pr_to_pg(pr)
}

# -----------------------------------------------------------------------------
# 样品注释（恰好两组）
# -----------------------------------------------------------------------------
normalize_sample_token <- function(x) {
  x <- as.character(x)
  x <- sub("\\\\", "/", x)
  x <- basename(x)
  x <- sub("\\.(raw|wiff|mzML|d)$", "", x, ignore.case = TRUE)
  x
}

load_annotation <- function(dir, intensity_cols) {
  anno_path <- file.path(dir, "sample_annotation.csv")
  tokens <- vapply(intensity_cols, normalize_sample_token, character(1), USE.NAMES = FALSE)
  if (!file.exists(anno_path)) {
    stop(
      "缺少 sample_annotation.csv（需要列 sample,group，恰好两组病人）。\n",
      "当前强度列: ", paste(tokens, collapse = ", ")
    )
  }
  anno <- readr::read_csv(anno_path, show_col_types = FALSE)
  anno <- as.data.frame(anno)
  names(anno) <- tolower(names(anno))
  if (!all(c("sample", "group") %in% names(anno))) {
    stop("sample_annotation.csv 必须包含 sample, group 两列")
  }
  anno$sample_token <- vapply(anno$sample, normalize_sample_token, character(1))
  anno$group <- as.character(anno$group)
  groups <- unique(anno$group)
  if (length(groups) != 2) {
    stop("必须恰好两组病人，当前: ", paste(groups, collapse = ", "))
  }
  map <- match(tokens, anno$sample_token)
  if (any(is.na(map))) {
    map2 <- match(tokens, anno$sample)
    map[is.na(map)] <- map2[is.na(map)]
  }
  if (any(is.na(map))) {
    miss <- tokens[is.na(map)]
    stop("注释表对不上这些样品列: ", paste(miss, collapse = ", "))
  }
  sample_info <- data.frame(
    column = intensity_cols,
    sample = tokens,
    group = anno$group[map],
    stringsAsFactors = FALSE
  )
  sample_info$group <- factor(sample_info$group, levels = groups)
  sample_info
}

# -----------------------------------------------------------------------------
# 过滤 + log2 + 样品中位数中心化
# -----------------------------------------------------------------------------
preprocess <- function(mat, sample_info) {
  int_cols <- sample_info$column
  raw <- as.matrix(mat[, int_cols, drop = FALSE])
  storage.mode(raw) <- "double"
  raw[!is.finite(raw) | raw < 0] <- NA_real_

  keep <- rep(FALSE, nrow(raw))
  for (g in levels(sample_info$group)) {
    cols <- sample_info$column[sample_info$group == g]
    frac <- rowMeans(is.finite(raw[, cols, drop = FALSE]) & raw[, cols, drop = FALSE] > 0, na.rm = FALSE)
    keep <- keep | (frac >= min_detect_frac)
  }
  log_msg("过滤前 ", nrow(raw), " 蛋白，过滤后 ", sum(keep))
  mat <- mat[keep, , drop = FALSE]
  raw <- raw[keep, , drop = FALSE]

  log2x <- log2(raw + 1)
  med <- apply(log2x, 2, median, na.rm = TRUE)
  log2x <- sweep(log2x, 2, med - median(med, na.rm = TRUE), "-")
  list(meta = mat, log2 = log2x)
}

# -----------------------------------------------------------------------------
# 两组统计：丰度排名 + FC；有重复才算 p
# -----------------------------------------------------------------------------
rank_proteins <- function(meta, log2x, sample_info) {
  g <- levels(sample_info$group)
  g1 <- g[1]
  g2 <- g[2]
  c1 <- sample_info$column[sample_info$group == g1]
  c2 <- sample_info$column[sample_info$group == g2]
  m1 <- rowMeans(log2x[, c1, drop = FALSE], na.rm = TRUE)
  m2 <- rowMeans(log2x[, c2, drop = FALSE], na.rm = TRUE)
  mean_ab <- rowMeans(log2x, na.rm = TRUE)
  log2fc <- m2 - m1
  n1 <- length(c1)
  n2 <- length(c2)
  pval <- rep(NA_real_, length(mean_ab))
  if (n1 >= 2 && n2 >= 2) {
    for (i in seq_along(mean_ab)) {
      a <- log2x[i, c1]
      b <- log2x[i, c2]
      a <- a[is.finite(a)]
      b <- b[is.finite(b)]
      if (length(a) >= 2 && length(b) >= 2 && (sd(a) > 0 || sd(b) > 0)) {
        pval[i] <- tryCatch(t.test(b, a, var.equal = FALSE)$p.value, error = function(e) NA_real_)
      }
    }
  } else {
    log_msg("每组样品不足 2，不伪造 p 值，只按丰度/FC 排名")
  }
  gene_raw <- if ("Genes" %in% names(meta)) meta$Genes else NA_character_
  pg <- if ("Protein.Group" %in% names(meta)) meta$Protein.Group else rownames(meta)
  gene <- vapply(as.character(gene_raw), pick_official_symbol, character(1), USE.NAMES = FALSE)
  gene[is.na(gene) | !nzchar(gene)] <- as.character(pg[is.na(gene) | !nzchar(gene)])
  gene[duplicated(gene)] <- paste0(gene[duplicated(gene)], "_", which(duplicated(gene)))

  out <- data.frame(
    Protein.Group = as.character(pg),
    gene = gene,
    description = if ("First.Protein.Description" %in% names(meta)) as.character(meta$First.Protein.Description) else NA_character_,
    mean_log2 = as.numeric(mean_ab),
    mean_group1 = as.numeric(m1),
    mean_group2 = as.numeric(m2),
    log2FC = as.numeric(log2fc),
    FoldChange = 2^as.numeric(log2fc),
    pvalue = pval,
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "mean_group1"] <- paste0("mean_", g1)
  names(out)[names(out) == "mean_group2"] <- paste0("mean_", g2)
  out$abundance_rank <- rank(-out$mean_log2, ties.method = "first")
  out$fc_rank <- rank(-abs(out$log2FC), ties.method = "first")
  out$neglog10p <- ifelse(is.finite(out$pvalue) & out$pvalue > 0, -log10(out$pvalue), NA_real_)
  out <- out[order(out$abundance_rank), ]
  attr(out, "group1") <- g1
  attr(out, "group2") <- g2
  out
}

# -----------------------------------------------------------------------------
# 作图
# -----------------------------------------------------------------------------
theme_bubble <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      axis.text.y = ggplot2::element_text(size = 9, face = "italic")
    )
}

plot_ranked_fc_bubble <- function(df, n, g1, g2, outfile) {
  sub <- df[df$abundance_rank <= n, ]
  sub <- sub[order(sub$abundance_rank, decreasing = TRUE), ]
  sub$gene <- factor(sub$gene, levels = sub$gene)
  title <- paste0("血清蛋白丰度排名气泡图 (top", n, ")")
  subtitle <- paste0("x = log2FC(", g2, " / ", g1, ")；点大小 = 平均 log2 丰度")
  p <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$log2FC, y = .data$gene, size = .data$mean_log2, color = .data$log2FC)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(alpha = 0.88) +
    ggplot2::scale_size_continuous(name = "mean log2\nabundance", range = c(2.5, 11)) +
    ggplot2::scale_color_gradient2(
      name = "log2FC",
      low = "#2166AC", mid = "grey75", high = "#B2182B", midpoint = 0
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, x = paste0("log2 Fold Change (", g2, " / ", g1, ")"), y = paste0("Protein (abundance rank 1–", n, ")")) +
    theme_bubble()
  ggplot2::ggsave(outfile, p, width = 8.5, height = max(4.5, 0.28 * n + 2.2), dpi = 150, limitsize = FALSE)
}

plot_two_group_bubble <- function(df, n, g1, g2, outfile) {
  sub <- df[df$abundance_rank <= n, ]
  sub <- sub[order(sub$abundance_rank, decreasing = TRUE), ]
  long <- tidyr::pivot_longer(
    sub,
    cols = c(paste0("mean_", g1), paste0("mean_", g2)),
    names_to = "group",
    values_to = "group_mean"
  )
  long$group <- sub("^mean_", "", long$group)
  long$group <- factor(long$group, levels = c(g1, g2))
  long$gene <- factor(long$gene, levels = unique(sub$gene))
  title <- paste0("两组病人血清丰度气泡图 (top", n, ")")
  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$group, y = .data$gene, size = .data$group_mean, color = .data$group)) +
    ggplot2::geom_point(alpha = 0.88) +
    ggplot2::scale_size_continuous(name = "group mean\nlog2 abundance", range = c(2.5, 11)) +
    ggplot2::scale_color_manual(values = setNames(c("#4C78A8", "#F58518"), c(g1, g2)), name = "Group") +
    ggplot2::labs(title = title, subtitle = "同一蛋白两个点 = 两组平均丰度；按全样品平均丰度降序", x = NULL, y = paste0("Protein (abundance rank 1–", n, ")")) +
    theme_bubble()
  ggplot2::ggsave(outfile, p, width = 7.2, height = max(4.5, 0.28 * n + 2.2), dpi = 150, limitsize = FALSE)
}

# -----------------------------------------------------------------------------
# 主流程
# -----------------------------------------------------------------------------
log_msg("工作目录: ", project_dir)
mat <- load_protein_matrix(project_dir)
sample_info <- load_annotation(project_dir, attr(mat, "intensity_cols"))
log_msg("组别: ", paste(levels(sample_info$group), collapse = " vs "))
prep <- preprocess(mat, sample_info)
ranked <- rank_proteins(prep$meta, prep$log2, sample_info)
g1 <- attr(ranked, "group1")
g2 <- attr(ranked, "group2")

rank_path <- file.path(result_dir, "protein_abundance_ranking.csv")
readr::write_csv(ranked, rank_path)
log_msg("写出排名表: ", rank_path)

for (n in top_ns) {
  n_use <- min(n, nrow(ranked))
  if (n_use < 1) next
  tag <- paste0("top", n_use)
  sub_path <- file.path(result_dir, paste0(tag, "_ranked_proteins.csv"))
  readr::write_csv(ranked[ranked$abundance_rank <= n_use, ], sub_path)
  plot_ranked_fc_bubble(
    ranked, n_use, g1, g2,
    file.path(result_dir, paste0(tag, "_abundance_rank_bubble.pdf"))
  )
  plot_two_group_bubble(
    ranked, n_use, g1, g2,
    file.path(result_dir, paste0(tag, "_two_group_abundance_bubble.pdf"))
  )
  log_msg("完成 ", tag)
}

log_msg("全部完成 -> ", result_dir)
