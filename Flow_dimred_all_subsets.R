#!/usr/bin/env Rscript
# =============================================================================
# 三 panel 全亚群总览（不修改 P1/P2/P3 各自的降维）
#
# P1/P2/P3 抗体不同，禁止拼成一张表达矩阵再 UMAP/tSNE。
# 本脚本在三个 panel 已经圈完亚群之后，只汇总各管内频率，做一个
# 「所有亚群」的 T vs T6 比较。百分数是该 panel 管子里的比例，
# 不是同一群细胞，也不能加总成全血组成。
#
# 细亚群以对应 panel 为准：
#   P1 = T/NK，P2 = B，P3 = 髓系
# P1 的 B/Myeloid、P3 的 T/B/NK 只是 dump 通道，图里标成 dump。
#
# 用法（三个 panel 跑完之后）：
#   setwd("E:/R/flow J-LJY WJZ ZZX")
#   source("Flow_dimred_pipeline.R")          # 各 panel 降维 + 自动调用本脚本
#   source("Flow_dimred_all_subsets.R")       # 只重出总览，不重跑 UMAP
#
# 结果：results_flow/all_subsets/
# =============================================================================

load_flow_pipeline_functions <- function() {
  if (exists("save_gg", mode = "function") && exists("lineage_frequencies", mode = "function")) {
    return(invisible(TRUE))
  }
  pipe <- "Flow_dimred_pipeline.R"
  cands <- c(pipe, file.path(getwd(), pipe))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    cands <- c(file.path(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))), pipe), cands)
  }
  ofile <- NULL
  if (sys.nframe() > 0) {
    for (i in sys.nframe():1) {
      ofile <- sys.frame(i)$ofile
      if (!is.null(ofile) && nzchar(ofile)) break
    }
  }
  if (!is.null(ofile) && nzchar(ofile)) {
    cands <- c(file.path(dirname(normalizePath(ofile, mustWork = FALSE)), pipe), cands)
  }
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit) || !nzchar(hit)) stop("找不到 Flow_dimred_pipeline.R")
  old <- Sys.getenv("FLOW_FUNCTIONS_ONLY", unset = NA)
  Sys.setenv(FLOW_FUNCTIONS_ONLY = "1")
  source(hit, local = FALSE)
  if (is.na(old)) Sys.unsetenv("FLOW_FUNCTIONS_ONLY") else Sys.setenv(FLOW_FUNCTIONS_ONLY = old)
  invisible(TRUE)
}

# dump = 该 panel 不是用来细拆这一系的；细亚群看对应 panel
subset_role <- function(panel_id, lineage) {
  lin <- as.character(lineage)
  if (identical(panel_id, "P1") && lin %in% c("B", "Myeloid")) return("dump")
  if (identical(panel_id, "P3") && lin %in% c("T", "B", "NK", "NKT", "CD4", "CD8")) return("dump")
  "focus"
}

read_panel_embeddings <- function(result_dir, panel_id) {
  p <- file.path(result_dir, panel_id, paste0(panel_id, "_cell_embeddings.csv"))
  if (!file.exists(p)) return(NULL)
  df <- read_embed_csv(p)
  need <- c("sample", "group", "lineage")
  if (!all(need %in% names(df))) return(NULL)
  df[, need, drop = FALSE]
}

collect_all_subset_frequencies <- function(cells_by_panel) {
  rows <- lapply(names(cells_by_panel), function(pn) {
    cells <- cells_by_panel[[pn]]
    if (is.null(cells) || !is.data.frame(cells) || nrow(cells) < 1) return(NULL)
    tab <- lineage_frequencies(cells)
    tab$panel <- pn
    tab$celltype <- celltype_label(tab$lineage, pn)
    tab$role <- vapply(as.character(tab$lineage), function(x) subset_role(pn, x), character(1))
    tab$subset_label <- paste0(pn, " · ", tab$celltype)
    tab
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

all_subset_stats <- function(freq) {
  pieces <- lapply(split(freq, freq$panel), function(d) {
    s <- compare_group_freq(d, "subset_label")
    meta <- unique(d[, c("subset_label", "panel", "lineage", "celltype", "role")])
    merge(s, meta, by = "subset_label", all.x = TRUE)
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out[order(out$panel, out$p_value, na.last = TRUE), ]
}

plot_all_freq_facet <- function(freq, title) {
  meta <- unique(freq[, c("panel", "celltype", "subset_label")])
  meta <- meta[order(meta$panel, meta$celltype), ]
  freq$subset_label <- factor(freq$subset_label, levels = unique(meta$subset_label))
  ggplot2::ggplot(freq, ggplot2::aes(x = celltype, y = percent, fill = group)) +
    ggplot2::stat_summary(fun = mean, geom = "col",
                          position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.8), size = 1.4, alpha = 0.85) +
    ggplot2::facet_wrap(~panel, scales = "free_x", nrow = 1) +
    ggplot2::scale_fill_manual(values = pal_group) +
    theme_dr() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8)) +
    ggplot2::labs(title = title, x = NULL, y = "% of cells in that panel")
}

plot_all_stacked <- function(freq, title) {
  mean_df <- aggregate(percent ~ group + panel + celltype, data = freq, FUN = mean)
  names(mean_df)[names(mean_df) == "percent"] <- "mean_percent"
  levs <- unique(mean_df$celltype)
  ggplot2::ggplot(mean_df, ggplot2::aes(x = group, y = mean_percent, fill = celltype)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::facet_wrap(~panel, nrow = 1) +
    ggplot2::scale_fill_manual(values = celltype_colors(levs)) +
    theme_dr() +
    ggplot2::labs(title = title, x = NULL, y = "Mean % of cells in that panel")
}

plot_all_heatmap <- function(stats_df, title) {
  long <- rbind(
    data.frame(subset_label = stats_df$subset_label, panel = stats_df$panel,
               group = "T", percent = stats_df$mean_T, stringsAsFactors = FALSE),
    data.frame(subset_label = stats_df$subset_label, panel = stats_df$panel,
               group = "T6", percent = stats_df$mean_T6, stringsAsFactors = FALSE)
  )
  ord <- stats_df$subset_label[order(stats_df$panel, -pmax(stats_df$mean_T, stats_df$mean_T6))]
  long$subset_label <- factor(long$subset_label, levels = rev(unique(ord)))
  long$group <- factor(long$group, levels = c("T", "T6"))
  ggplot2::ggplot(long, ggplot2::aes(x = group, y = subset_label, fill = percent)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#08306B") +
    theme_dr() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8)) +
    ggplot2::labs(title = title, x = NULL, y = NULL, fill = "Mean %")
}

plot_all_lollipop <- function(stats_df, title) {
  d <- stats_df[order(stats_df$delta_T6_minus_T), ]
  d$subset_label <- factor(d$subset_label, levels = d$subset_label)
  ggplot2::ggplot(d, ggplot2::aes(x = delta_T6_minus_T, y = subset_label, color = panel)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = delta_T6_minus_T,
                                       y = subset_label, yend = subset_label),
                          colour = "grey75") +
    ggplot2::geom_point(size = 2.4) +
    theme_dr() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8)) +
    ggplot2::labs(title = title, x = "T6 - T (percentage points)", y = NULL)
}

export_all_subsets_analysis <- function(result_dir) {
  if (missing(result_dir) || !nzchar(result_dir)) {
    if (exists("result_dir", envir = .GlobalEnv, inherits = FALSE)) {
      result_dir <- get("result_dir", envir = .GlobalEnv)
    } else {
      result_dir <- file.path(getwd(), "results_flow")
    }
  }
  panels <- c("P1", "P2", "P3")
  cells_by_panel <- lapply(panels, function(pn) read_panel_embeddings(result_dir, pn))
  names(cells_by_panel) <- panels
  n_ok <- sum(vapply(cells_by_panel, function(x) !is.null(x) && nrow(x) > 0, logical(1)))
  if (n_ok < 1) {
    log_msg("all-subsets: no panel embeddings under ", result_dir, " ; run Flow_dimred_pipeline.R first")
    return(invisible(NULL))
  }
  freq <- collect_all_subset_frequencies(cells_by_panel)
  if (is.null(freq) || !nrow(freq)) {
    log_msg("all-subsets: empty frequency table")
    return(invisible(NULL))
  }
  stats <- all_subset_stats(freq)
  out_dir <- file.path(result_dir, "all_subsets")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  note <- c(
    "P1/P2/P3 cannot be concatenated for a joint UMAP: the antibody panels differ.",
    "Numbers below are % of cells inside each panel tube, not a whole-blood composition.",
    "Do not add P1+P2+P3 percentages together.",
    "Focus subsets: P1 T/NK, P2 B, P3 myeloid.",
    "Dump channels (coarse only): P1 B and Myeloid; P3 T, B, NK.",
    paste("Panels found:", paste(unique(freq$panel), collapse = ", "))
  )
  writeLines(note, file.path(out_dir, "ALL_SUBSETS_NOTE.txt"))
  utils::write.csv(freq, file.path(out_dir, "all_subsets_frequency_by_sample.csv"), row.names = FALSE)
  utils::write.csv(stats, file.path(out_dir, "all_subsets_T6_vs_T_stats.csv"), row.names = FALSE)

  n_lab <- length(unique(freq$subset_label))
  w_facet <- max(12, 1.1 * n_lab / 2)
  save_gg(plot_all_freq_facet(freq, "All subsets (within-panel %)  T vs T6"),
          file.path(out_dir, "all_subsets_frequency_T6_vs_T"),
          width = w_facet, height = 5.8)
  save_gg(plot_all_stacked(freq, "Within-panel composition (sums to 100% per panel)"),
          file.path(out_dir, "all_subsets_composition_stacked"),
          width = 11, height = 5.6)

  focus <- freq[freq$role == "focus", , drop = FALSE]
  stats_f <- stats[stats$role == "focus", , drop = FALSE]
  if (nrow(focus) > 0) {
    save_gg(plot_all_freq_facet(focus, "Focus subsets only  (P1 T/NK, P2 B, P3 myeloid)"),
            file.path(out_dir, "all_subsets_focus_frequency_T6_vs_T"),
            width = max(11, 1.15 * length(unique(focus$subset_label)) / 2), height = 5.8)
  }
  if (nrow(stats_f) > 0) {
    save_gg(plot_all_heatmap(stats_f, "Focus subset mean %  (within panel)"),
            file.path(out_dir, "all_subsets_focus_mean_heatmap"),
            width = 6.5, height = max(6, 0.28 * nrow(stats_f) + 2))
    save_gg(plot_all_lollipop(stats_f, "Focus subsets  T6 minus T"),
            file.path(out_dir, "all_subsets_focus_delta_lollipop"),
            width = 8.5, height = max(6, 0.28 * nrow(stats_f) + 2))
  }

  log_msg("All-subset summary: ", out_dir,
          " (", n_ok, " panels, ", n_lab, " subset labels)")
  invisible(list(frequency = freq, stats = stats, out_dir = out_dir))
}

load_flow_pipeline_functions()

if (!identical(toupper(Sys.getenv("FLOW_ALL_SUBSETS_FROM_PIPELINE", "0")), "1") &&
    !identical(toupper(Sys.getenv("FLOW_FUNCTIONS_ONLY", "0")), "1")) {
  export_all_subsets_analysis(if (exists("result_dir")) result_dir else file.path(getwd(), "results_flow"))
}
