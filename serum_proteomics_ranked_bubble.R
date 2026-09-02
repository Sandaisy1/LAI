#!/usr/bin/env Rscript
# =============================================================================
# 两组样品取平均丰度，按排名画等大气泡图（x=丰度排名，y=蛋白丰度值）
# 输入：DIA-NN report.pg_matrix（首选）或 report.pr_matrix
# 数据目录：E:/天府/实验管理/课题/赵章寻/血清蛋白质组学
# 样品：GP_WJZ_11、GP_WJZ_18 取平均后再排名
#
# 运行（当前工作目录必须能看到这个文件，否则会报“无法打开链接”）：
#   source("run_serum_proteomics_bubble.R", encoding = "UTF-8")
# 或：
#   setwd("E:/天府/实验管理/课题/赵章寻/血清蛋白质组学")  # 先把本文件复制过去
#   source("serum_proteomics_ranked_bubble.R", encoding = "UTF-8")
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
default_proteomics_dirs <- c(
  "E:/天府/实验管理/课题/赵章寻/血清蛋白质组学",
  "E:\\天府\\实验管理\\课题\\赵章寻\\血清蛋白质组学"
)

resolve_proteomics_dir <- function() {
  env_dir <- Sys.getenv("SERUM_PROTEOMICS_DIR", unset = "")
  script_dir <- tryCatch({
    ofile <- sys.frame(1)$ofile
    if (!is.null(ofile) && nzchar(ofile)) {
      dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE))
    } else {
      NA_character_
    }
  }, error = function(e) NA_character_)
  candidates <- unique(c(
    env_dir,
    default_proteomics_dirs,
    script_dir,
    file.path(getwd(), "serum_proteomics"),
    getwd()
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
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

top_ns <- c(20L, 30L, 50L)
min_detect_frac <- 0.5

log_file <- NULL
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  if (!is.null(log_file)) cat(msg, "\n", file = log_file, append = TRUE)
}

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

is_immunoglobulin_symbol <- function(sym) {
  s <- toupper(trimws(as.character(sym)))
  if (!nzchar(s) || s %in% c("NA", "-", ".")) return(FALSE)
  if (s %in% c("JCHAIN", "IGJ")) return(TRUE)
  grepl("^(IGH|IGK|IGL)", s)
}

is_trypsin_symbol <- function(sym) {
  s <- toupper(trimws(as.character(sym)))
  s <- gsub("-", "_", s)
  s <- sub("^(CON__|SP\\||TR\\|)", "", s)
  if (grepl("\\|", s)) s <- sub(".*\\|", "", s)
  if (s %in% c(
    "PRSS1", "PRSS2", "PRSS3", "TRY1", "TRY2", "TRY3", "TRYP", "TRYP_PIG",
    "TRY1_BOVIN", "TRY2_BOVIN", "TRY3_BOVIN",
    "P00761", "P00760", "P00763", "P07477", "P07478", "P35030"
  )) return(TRUE)
  grepl("^(PRSS[123]|TRY[123]|TRYP)(_|$)", s)
}

is_immunoglobulin_text <- function(...) {
  blob <- tolower(paste(..., collapse = " "))
  if (!nzchar(trimws(blob))) return(FALSE)
  if (grepl("铁蛋白重链|ferritin heavy", blob)) return(FALSE)
  if (grepl("immunoglobulin superfamily", blob) &&
      !grepl("\\bigh[agmdvejk]|\\bigk[cvlj]|\\bigl[cvlj]", blob)) {
    return(FALSE)
  }
  grepl("immunoglobulin\\s+(heavy|kappa|lambda|alpha|gamma|mu|delta|epsilon)", blob) ||
    grepl("immunoglobulin heavy chain", blob) ||
    grepl("ig heavy chain", blob) ||
    grepl("免疫球蛋白重链", blob) ||
    grepl("免疫球蛋白", blob) ||
    grepl("\\big\\s*(heavy|kappa|lambda|gamma|alpha|mu)\\b", blob) ||
    grepl("\\big\\s+(gamma|alpha|mu|kappa|lambda)\\b", blob) ||
    grepl("heavy chain c region", blob)
}

is_trypsin_text <- function(...) {
  blob <- tolower(paste(..., collapse = " "))
  if (!nzchar(trimws(blob))) return(FALSE)
  if (grepl("antitrypsin|anti-trypsin|trypsin inhibitor|抗胰蛋白酶", blob)) return(FALSE)
  if (grepl("tryptophan", blob)) return(FALSE)
  grepl("\\btrypsin(ogen)?\\b|胰蛋白酶", blob)
}

exclusion_reasons <- function(meta) {
  n <- nrow(meta)
  if (n == 0) return(character(0))
  gene_raw <- if ("Genes" %in% names(meta)) as.character(meta$Genes) else rep("", n)
  names_ <- if ("Protein.Names" %in% names(meta)) as.character(meta$Protein.Names) else rep("", n)
  desc <- if ("First.Protein.Description" %in% names(meta)) as.character(meta$First.Protein.Description) else rep("", n)
  ids <- if ("Protein.Ids" %in% names(meta)) as.character(meta$Protein.Ids) else rep("", n)
  pg <- if ("Protein.Group" %in% names(meta)) as.character(meta$Protein.Group) else rep("", n)
  vapply(seq_len(n), function(i) {
    parts <- trimws(unlist(strsplit(paste(c(gene_raw[i], names_[i], desc[i], ids[i], pg[i]), collapse = ","), "[,;|/ ]+")))
    parts <- parts[nzchar(parts)]
    tags <- character(0)
    if (any(vapply(parts, is_immunoglobulin_symbol, logical(1))) ||
        is_immunoglobulin_text(gene_raw[i], names_[i], desc[i], ids[i], pg[i])) {
      tags <- c(tags, "immunoglobulin")
    }
    if (any(vapply(parts, is_trypsin_symbol, logical(1))) ||
        is_trypsin_text(gene_raw[i], names_[i], desc[i], ids[i], pg[i])) {
      tags <- c(tags, "trypsin")
    }
    paste(tags, collapse = ";")
  }, character(1))
}

immunoglobulin_mask <- function(meta) {
  grepl("immunoglobulin", exclusion_reasons(meta))
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

write_annotation_template <- function(path, tokens) {
  if (length(unique(tokens)) == 2) {
    anno <- data.frame(sample = unique(tokens), group = unique(tokens), stringsAsFactors = FALSE)
  } else {
    anno <- data.frame(sample = tokens, group = "", stringsAsFactors = FALSE)
  }
  readr::write_excel_csv(anno, path)
  log_msg("已写出 ", path)
}

load_annotation <- function(dir, intensity_cols) {
  anno_path <- file.path(dir, "sample_annotation.csv")
  tokens <- vapply(intensity_cols, normalize_sample_token, character(1), USE.NAMES = FALSE)
  if (!file.exists(anno_path)) {
    write_annotation_template(anno_path, tokens)
    if (length(unique(tokens)) != 2) {
      stop(
        "缺少 sample_annotation.csv。已按样品列写出模板，请把 group 填成恰好两组后重新 source。\n",
        "文件: ", anno_path, "\n",
        "当前强度列: ", paste(tokens, collapse = ", ")
      )
    }
    log_msg("矩阵恰好两列样品，按 1-vs-1 分组（", paste(unique(tokens), collapse = " vs "), "），不伪造 p 值")
  }
  anno <- readr::read_csv(anno_path, show_col_types = FALSE)
  anno <- as.data.frame(anno)
  names(anno) <- tolower(names(anno))
  if (!all(c("sample", "group") %in% names(anno))) {
    stop("sample_annotation.csv 必须包含 sample, group 两列")
  }
  anno$sample_token <- vapply(anno$sample, normalize_sample_token, character(1))
  anno$group <- trimws(as.character(anno$group))
  groups <- unique(anno$group[!is.na(anno$group) & nzchar(anno$group) & !anno$group %in% c("NA", "NaN")])
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
# 过滤 + log2 + 样品中位数中心化；两样品后面再取平均
# -----------------------------------------------------------------------------
preprocess <- function(mat, int_cols) {
  raw <- as.matrix(mat[, int_cols, drop = FALSE])
  storage.mode(raw) <- "double"
  raw[!is.finite(raw) | raw < 0] <- NA_real_
  frac <- rowMeans(is.finite(raw) & raw > 0)
  keep <- frac >= min_detect_frac
  log_msg("过滤前 ", nrow(raw), " 蛋白，过滤后 ", sum(keep),
          "；样品: ", paste(int_cols, collapse = ", "))
  mat <- mat[keep, , drop = FALSE]
  raw <- raw[keep, , drop = FALSE]
  log2x <- log2(raw + 1)
  med <- apply(log2x, 2, median, na.rm = TRUE)
  log2x <- sweep(log2x, 2, med - median(med, na.rm = TRUE), "-")
  list(meta = mat, log2 = log2x, sample_cols = int_cols)
}

# -----------------------------------------------------------------------------
# 两样品平均丰度，再按该值降序排名
# -----------------------------------------------------------------------------
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
    description = if ("First.Protein.Description" %in% names(meta)) {
      as.character(meta$First.Protein.Description)
    } else {
      NA_character_
    },
    mean_abundance = as.numeric(mean_ab),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(sample_cols)) {
    out[[paste0("sample_", normalize_sample_token(sample_cols[[i]]))]] <- as.numeric(log2x[, i])
  }
  out$abundance_rank <- rank(-out$mean_abundance, ties.method = "first")
  out[order(out$abundance_rank), ]
}

# -----------------------------------------------------------------------------
# 作图：x=丰度排名，y=蛋白丰度值，气泡大小一致
# -----------------------------------------------------------------------------
theme_bubble <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0)
    )
}

plot_rank_abundance_bubble <- function(df, n, outfile, label_n = 12L,
                                       title_prefix = "血清蛋白丰度排名气泡图",
                                       subtitle = "两样品平均丰度后排名；气泡大小一致") {
  sub <- df[df$abundance_rank <= n, ]
  sub <- sub[order(sub$abundance_rank), ]
  lab <- sub[sub$abundance_rank <= min(as.integer(label_n), n), ]
  title <- paste0(title_prefix, " (top", n, ")")
  p <- ggplot2::ggplot(sub, ggplot2::aes(x = .data$abundance_rank, y = .data$mean_abundance)) +
    ggplot2::geom_line(color = "#9ECAE1", size = 0.45) +
    ggplot2::geom_point(size = 2.8, color = "#2C7FB8", alpha = 0.92) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "丰度排名",
      y = "蛋白丰度值"
    ) +
    theme_bubble()
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

# -----------------------------------------------------------------------------
# 主流程（去 Ig 脚本可设 SERUM_PROTEOMICS_SKIP_MAIN 后自行调用）
# -----------------------------------------------------------------------------
run_serum_abundance_bubble <- function(
  result_subdir = "serum_proteomics_bubble",
  drop_immunoglobulin = FALSE,
  title_prefix = "血清蛋白丰度排名气泡图",
  subtitle = "两样品平均丰度后排名；气泡大小一致"
) {
  project_dir <- resolve_proteomics_dir()
  result_dir <- file.path(project_dir, "results", result_subdir)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <<- file.path(result_dir, paste0("log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  log_msg("工作目录: ", project_dir)
  mat <- load_protein_matrix(project_dir)
  int_cols <- attr(mat, "intensity_cols")
  prep <- preprocess(mat, int_cols)
  if (isTRUE(drop_immunoglobulin)) {
    reasons <- exclusion_reasons(prep$meta)
    drop <- nzchar(reasons)
    dropped <- prep$meta[drop, , drop = FALSE]
    drop_path <- file.path(result_dir, "removed_Ig_trypsin.csv")
    if (nrow(dropped) > 0) {
      keep_cols <- intersect(c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description"), names(dropped))
      out <- dropped[, keep_cols, drop = FALSE]
      out$reason <- reasons[drop]
      readr::write_csv(out, drop_path)
    }
    log_msg("去除免疫球蛋白/重链与胰蛋白酶 ", sum(drop), " / ", length(drop), " -> ", drop_path)
    prep$meta <- prep$meta[!drop, , drop = FALSE]
    prep$log2 <- prep$log2[!drop, , drop = FALSE]
  }
  ranked <- rank_proteins(prep$meta, prep$log2, prep$sample_cols)
  log_msg("按两样品平均丰度排名，蛋白数 ", nrow(ranked))
  rank_path <- file.path(result_dir, "protein_abundance_ranking.csv")
  readr::write_csv(ranked, rank_path)
  log_msg("写出排名表: ", rank_path)
  ns <- unique(c(top_ns[top_ns <= nrow(ranked)], nrow(ranked)))
  for (n_use in ns) {
    if (n_use < 1) next
    tag <- if (n_use == nrow(ranked)) "all" else paste0("top", n_use)
    readr::write_csv(ranked[ranked$abundance_rank <= n_use, ], file.path(result_dir, paste0(tag, "_ranked_proteins.csv")))
    plot_rank_abundance_bubble(
      ranked, n_use,
      file.path(result_dir, paste0(tag, "_abundance_rank_bubble.pdf")),
      title_prefix = title_prefix,
      subtitle = subtitle
    )
    log_msg("完成 ", tag)
  }
  log_msg("全部完成 -> ", result_dir)
  invisible(result_dir)
}

if (!exists("SERUM_PROTEOMICS_SKIP_MAIN") || !isTRUE(SERUM_PROTEOMICS_SKIP_MAIN)) {
  run_serum_abundance_bubble()
}