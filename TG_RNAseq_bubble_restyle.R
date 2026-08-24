#!/usr/bin/env Rscript
# =============================================================================
# 按 TG_bubble_plot_style.R 重画气泡图（不重跑 enrichGO）
#
#   setwd("E:/R/TG_BRCA/TG")
#   # 先改 TG_bubble_plot_style.R 里的气泡大小 / 坐标字体
#   source("TG_RNAseq_bubble_restyle.R")
#
# 也可改某张图旁边的 *_plotdata.csv（删行或改顺序）后再运行本脚本。
# 被 TG_RNAseq_pipeline.R 以 skip_run 方式加载时，只提供作图函数。
# =============================================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("需要 ggplot2。请先安装：install.packages(\"ggplot2\")")
}

bubble_style_defaults <- function() {
  list(
    bubble_size_min = 6,
    bubble_size_max = 18,
    axis_text_y_size = 10,
    axis_text_x_size = 11,
    axis_title_size = 12,
    title_size = 12,
    legend_text_size = 10,
    legend_title_size = 11,
    base_size = 12,
    plot_width = 9,
    plot_height = NA_real_,
    point_stroke = 0.5,
    only_this_csv = ""
  )
}

find_bubble_style_file <- function() {
  env_f <- Sys.getenv("TG_BUBBLE_STYLE", unset = "")
  cwd <- tryCatch(normalizePath(getwd(), winslash = "/", mustWork = FALSE), error = function(e) getwd())
  candidates <- c(
    env_f,
    file.path(cwd, "TG_bubble_plot_style.R"),
    "E:/R/TG_BRCA/TG/TG_bubble_plot_style.R",
    file.path(cwd, "results", "00_bubble_restyle", "TG_bubble_plot_style.R")
  )
  if (exists("project_dir", inherits = TRUE)) {
    pd <- get("project_dir", inherits = TRUE)
    candidates <- c(
      file.path(pd, "TG_bubble_plot_style.R"),
      file.path(pd, "results", "00_bubble_restyle", "TG_bubble_plot_style.R"),
      candidates
    )
  }
  candidates <- unique(candidates[nzchar(candidates)])
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

load_bubble_style <- function(style_file = NULL) {
  st <- bubble_style_defaults()
  if (is.null(style_file) || !nzchar(style_file)) style_file <- find_bubble_style_file()
  if (!is.na(style_file) && file.exists(style_file)) {
    env <- list2env(st, parent = emptyenv())
    sys.source(style_file, envir = env)
    for (nm in names(st)) {
      if (exists(nm, envir = env, inherits = FALSE)) st[[nm]] <- get(nm, envir = env)
    }
    message("气泡图样式: ", style_file)
  } else {
    message("未找到 TG_bubble_plot_style.R，使用默认气泡大小和字体")
  }
  st
}

parse_gene_ratio_vec <- function(x) {
  vapply(as.character(x), function(s) {
    p <- strsplit(s, "/", fixed = TRUE)[[1]]
    if (length(p) < 2) return(suppressWarnings(as.numeric(s)))
    as.numeric(p[1]) / as.numeric(p[2])
  }, numeric(1), USE.NAMES = FALSE)
}

prepare_ora_bubble_df <- function(ora, n_show, padj_cutoff = 0.2) {
  df <- ora
  if (is.null(df) || nrow(df) == 0) return(NULL)
  if (!"GeneRatio_num" %in% names(df)) {
    df$GeneRatio_num <- if (is.numeric(df$GeneRatio)) df$GeneRatio else parse_gene_ratio_vec(df$GeneRatio)
  }
  df <- df[is.finite(df$GeneRatio_num), , drop = FALSE]
  df <- df[is.finite(df$p.adjust) & df$p.adjust < padj_cutoff, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  if (!"pvalue" %in% names(df)) df$pvalue <- NA_real_
  df <- df[order(-df$GeneRatio_num, df$p.adjust, df$pvalue), , drop = FALSE]
  df <- utils::head(df, n_show)
  gid <- if ("go_id" %in% names(df)) df$go_id else df$ID
  df$go_id <- as.character(gid)
  df$label <- paste0(df$Description, " (", df$go_id, ")")
  df
}

export_ora_bubble_plotdata <- function(df, title, outfile) {
  out <- df
  keep <- intersect(
    c("label", "Description", "go_id", "ID", "ONTOLOGY", "GeneRatio", "GeneRatio_num",
      "pvalue", "p.adjust", "qvalue", "Count", "geneID", "BgRatio"),
    names(out)
  )
  out <- out[, keep, drop = FALSE]
  out$plot_title <- title
  out$n_terms <- nrow(out)
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(out, paste0(outfile, "_plotdata.csv"), row.names = FALSE)
  if (requireNamespace("writexl", quietly = TRUE)) {
    tryCatch(
      writexl::write_xlsx(out, paste0(outfile, "_plotdata.xlsx")),
      error = function(e) NULL
    )
  }
  invisible(out)
}

draw_ora_bubble_gg <- function(df, title, style = NULL) {
  if (is.null(style)) style <- load_bubble_style()
  plot_df <- df
  if (!"label" %in% names(plot_df) || any(!nzchar(as.character(plot_df$label)))) {
    gid <- if ("go_id" %in% names(plot_df)) plot_df$go_id else plot_df$ID
    plot_df$label <- paste0(plot_df$Description, " (", gid, ")")
  }
  if (!"GeneRatio_num" %in% names(plot_df) || any(!is.finite(plot_df$GeneRatio_num))) {
    plot_df$GeneRatio_num <- if (is.numeric(plot_df$GeneRatio)) {
      plot_df$GeneRatio
    } else {
      parse_gene_ratio_vec(plot_df$GeneRatio)
    }
  }
  plot_df$label <- factor(plot_df$label, levels = rev(unique(as.character(plot_df$label))))
  plot_df$p_adjust <- plot_df$p.adjust
  plot_df$p_adjust[!is.finite(plot_df$p_adjust)] <- 1
  plot_df$p_adjust <- pmax(plot_df$p_adjust, 1e-300)
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = GeneRatio_num, y = label, size = Count, fill = p_adjust)
  ) +
    ggplot2::geom_point(shape = 21, color = "grey30", stroke = style$point_stroke) +
    ggplot2::scale_size_continuous(range = c(style$bubble_size_min, style$bubble_size_max)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.16))) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.55)) +
    ggplot2::theme_bw(base_size = style$base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = style$title_size),
      axis.text.y = ggplot2::element_text(size = style$axis_text_y_size),
      axis.text.x = ggplot2::element_text(size = style$axis_text_x_size),
      axis.title.x = ggplot2::element_text(size = style$axis_title_size),
      legend.text = ggplot2::element_text(size = style$legend_text_size),
      legend.title = ggplot2::element_text(size = style$legend_title_size),
      plot.margin = ggplot2::margin(6, 10, 6, 6),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = title,
      x = "GeneRatio",
      y = NULL,
      fill = "p.adjust",
      size = "Count"
    )
  rng <- range(plot_df$p_adjust[plot_df$p_adjust > 0], na.rm = TRUE)
  diverging <- c("blue", "white", "red")
  if (is.finite(rng[1]) && rng[1] > 0 && rng[2] / rng[1] >= 10) {
    p <- p + ggplot2::scale_fill_gradientn(colours = diverging, trans = "log10")
  } else {
    p <- p + ggplot2::scale_fill_gradientn(colours = diverging)
  }
  p
}

bubble_auto_height <- function(n_rows) {
  max(5.2, min(10.5, 0.32 * n_rows + 2.4))
}

save_ora_bubble_gg <- function(plot, outfile, style, n_rows) {
  width <- style$plot_width
  height <- style$plot_height
  if (length(height) != 1 || !is.finite(height)) height <- bubble_auto_height(n_rows)
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  note <- function(...) {
    if (exists("log_msg", mode = "function")) log_msg(...) else message(paste0(...))
  }
  tryCatch(
    ggplot2::ggsave(paste0(outfile, ".pdf"), plot, width = width, height = height),
    error = function(e) note("pdf ggsave failed: ", e$message)
  )
  tryCatch(
    ggplot2::ggsave(paste0(outfile, ".png"), plot, width = width, height = height, dpi = 300),
    error = function(e) note("png ggsave failed: ", e$message)
  )
}

read_ora_bubble_plotdata <- function(csv_path) {
  df <- utils::read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(NULL)
  if (!"p.adjust" %in% names(df) && "p_adjust" %in% names(df)) df$p.adjust <- df$p_adjust
  if (!"Count" %in% names(df) && "count" %in% names(df)) df$Count <- df$count
  df
}

plotdata_outfile_stub <- function(csv_path) {
  sub("_plotdata\\.csv$", "", csv_path, ignore.case = TRUE)
}

find_bubble_plotdata_files <- function(root, only_this_csv = "") {
  if (nzchar(only_this_csv)) {
    f <- only_this_csv
    if (!file.exists(f)) stop("找不到 only_this_csv: ", f)
    return(normalizePath(f, winslash = "/", mustWork = TRUE))
  }
  if (!dir.exists(root)) return(character())
  list.files(root, pattern = "_plotdata\\.csv$", full.names = TRUE, recursive = TRUE)
}

write_bubble_restyle_kit <- function(result_root, style_src = NA_character_) {
  kit <- file.path(result_root, "00_bubble_restyle")
  dir.create(kit, recursive = TRUE, showWarnings = FALSE)
  if (is.na(style_src) || !file.exists(style_src)) style_src <- find_bubble_style_file()
  if (!is.na(style_src) && file.exists(style_src)) {
    file.copy(style_src, file.path(kit, "TG_bubble_plot_style.R"), overwrite = TRUE)
  }
  writeLines(
    c(
      "改气泡大小、坐标字体（不用重跑 GO）：",
      "1. 编辑项目目录或本文件夹里的 TG_bubble_plot_style.R",
      "     bubble_size_min / bubble_size_max  = 气泡大小",
      "     axis_text_y_size                   = 左侧通路名字体",
      "     axis_text_x_size                   = 底部刻度字体",
      "2. setwd(\"E:/R/TG_BRCA/TG\")",
      "3. source(\"TG_RNAseq_bubble_restyle.R\")",
      "",
      "每张气泡图旁边的 *_plotdata.csv / .xlsx 是作图数据。",
      "可在 Excel 里删行或改顺序，保存回 CSV 后再运行 restyle。",
      "改 p.adjust 阈值或重做富集，需要重新运行 pipeline。"
    ),
    file.path(kit, "00_READ_ME.txt")
  )
  kit
}

restyle_ora_bubbles <- function(result_root = NULL, style_file = NULL) {
  if (is.null(result_root) || !nzchar(result_root)) {
    if (exists("result_dir", inherits = TRUE)) {
      result_root <- get("result_dir", inherits = TRUE)
    } else {
      result_root <- file.path(getwd(), "results")
    }
  }
  style <- load_bubble_style(style_file)
  files <- find_bubble_plotdata_files(result_root, style$only_this_csv)
  if (length(files) == 0) {
    stop("在 ", result_root, " 下没有找到 *_plotdata.csv。请先运行 TG_RNAseq_pipeline.R")
  }
  n_ok <- 0L
  for (f in files) {
    df <- tryCatch(read_ora_bubble_plotdata(f), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    title <- if ("plot_title" %in% names(df) && nzchar(df$plot_title[1])) {
      df$plot_title[1]
    } else {
      basename(plotdata_outfile_stub(f))
    }
    stub <- plotdata_outfile_stub(f)
    p <- draw_ora_bubble_gg(df, title, style)
    save_ora_bubble_gg(p, stub, style, nrow(df))
    n_ok <- n_ok + 1L
    message("redrawn: ", stub)
  }
  message("Done. Redrew ", n_ok, " bubble plot(s).")
  invisible(n_ok)
}

if (!isTRUE(getOption("tg.bubble.skip_run", FALSE))) {
  restyle_ora_bubbles()
}
