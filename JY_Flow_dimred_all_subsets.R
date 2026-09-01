#!/usr/bin/env Rscript
# =============================================================================
# JY：P1/P2/P3 亚群频率总览（JY-NNK / JY-EVNK）
# 只汇总本方案 results_flow/，不读 E:/R/fuction of cell 或 Internation。
# 入口仍是 source("JY_Flow_dimred_pipeline.R")
# =============================================================================

jy_keep_cand <- function(p) {
  s <- gsub("\\", "/", as.character(p))
  if (grepl("cell-ljy", s, ignore.case = TRUE)) return(TRUE)
  if (grepl("Internation cell immune", s, ignore.case = TRUE)) return(FALSE)
  if (grepl("(^|/)fuction of cell(/|$)|(^|/)function of cell(/|$)", s, ignore.case = TRUE)) return(FALSE)
  TRUE
}

load_jy_engine_functions <- function() {
  if (isTRUE(get0("jy_engine_loaded", ifnotfound = FALSE))) {
    return(invisible(TRUE))
  }
  pipe <- "JY_flow_engine.R"
  cands <- c(file.path(getwd(), pipe), pipe)
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
  cands <- unique(cands)
  cands <- cands[vapply(cands, jy_keep_cand, logical(1))]
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit) || !nzchar(hit)) {
    stop("找不到 JY_flow_engine.R。JY 方案不使用 Flow_dimred_pipeline.R。")
  }
  source(hit, local = FALSE)
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
  keep <- unique(c(need, "bio_sample", "tech_rep"))
  df[, intersect(keep, names(df)), drop = FALSE]
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

# 总览图用技术重复平均后再去掉极端生物学重复的值；原始管内频率另存 CSV
all_subset_freq_for_plots <- function(freq) {
  if (is.null(freq) || !nrow(freq)) return(freq)
  bio <- aggregate_freq_by_bio(freq, "subset_label")
  bio <- maybe_trim_bio(bio, "subset_label")
  dropped <- attr(bio, "dropped")
  meta_cols <- intersect(c("subset_label", "panel", "lineage", "celltype", "role"), names(freq))
  meta <- unique(freq[, meta_cols, drop = FALSE])
  out <- merge(bio, meta, by = "subset_label", all.x = TRUE)
  attr(out, "dropped") <- dropped
  out
}

plot_all_freq_facet <- function(freq, title) {
  meta <- unique(freq[, c("panel", "celltype", "subset_label")])
  meta <- meta[order(meta$panel, meta$celltype), ]
  freq$subset_label <- factor(freq$subset_label, levels = unique(meta$subset_label))
  ylab <- if (exists("flow_should_trim_bio", mode = "function") && flow_should_trim_bio()) {
    "% of cells in that panel (bio-rep; max/min dropped)"
  } else {
    "% of cells in that panel"
  }
  ggplot2::ggplot(freq, ggplot2::aes(x = celltype, y = percent, fill = group)) +
    ggplot2::stat_summary(fun = mean, geom = "col",
                          position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.8), size = 1.4, alpha = 0.85) +
    ggplot2::facet_wrap(~panel, scales = "free_x", nrow = 1) +
    ggplot2::scale_fill_manual(values = pal_group) +
    theme_dr() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8)) +
    ggplot2::labs(title = title, x = NULL, y = ylab)
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
  g1 <- flow_ctrl_group
  g2 <- flow_trt_group
  m1 <- stats_df[[paste0("mean_", g1)]]
  m2 <- stats_df[[paste0("mean_", g2)]]
  long <- rbind(
    data.frame(subset_label = stats_df$subset_label, panel = stats_df$panel,
               group = g1, percent = m1, stringsAsFactors = FALSE),
    data.frame(subset_label = stats_df$subset_label, panel = stats_df$panel,
               group = g2, percent = m2, stringsAsFactors = FALSE)
  )
  ord <- stats_df$subset_label[order(stats_df$panel, -pmax(m1, m2))]
  long$subset_label <- factor(long$subset_label, levels = rev(unique(ord)))
  long$group <- factor(long$group, levels = flow_group_levels)
  ggplot2::ggplot(long, ggplot2::aes(x = group, y = subset_label, fill = percent)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#08306B") +
    theme_dr() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8)) +
    ggplot2::labs(title = title, x = NULL, y = NULL, fill = "Mean %")
}

plot_all_lollipop <- function(stats_df, title) {
  delta_col <- paste0("delta_", flow_trt_group, "_minus_", flow_ctrl_group)
  d <- stats_df[order(stats_df[[delta_col]]), ]
  d$subset_label <- factor(d$subset_label, levels = d$subset_label)
  d$delta <- d[[delta_col]]
  ggplot2::ggplot(d, ggplot2::aes(x = delta, y = subset_label, color = panel)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = delta,
                                       y = subset_label, yend = subset_label),
                          colour = "grey75") +
    ggplot2::geom_point(size = 2.4) +
    theme_dr() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8)) +
    ggplot2::labs(title = title, x = paste0(flow_trt_group, " - ", flow_ctrl_group, " (percentage points)"), y = NULL)
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
    log_msg("all-subsets: no panel embeddings under ", result_dir, " ; run JY_Flow_dimred_pipeline.R first")
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
    "Stats: average tech reps, then drop 1 extreme bio-rep (max or min) per group.",
    paste("Panels found:", paste(unique(freq$panel), collapse = ", "))
  )
  writeLines(note, file.path(out_dir, "ALL_SUBSETS_NOTE.txt"))
  utils::write.csv(freq, file.path(out_dir, "all_subsets_frequency_by_sample.csv"), row.names = FALSE)
  utils::write.csv(stats, file.path(out_dir, "all_subsets_JY_NNK_vs_JY_EVNK_stats.csv"), row.names = FALSE)
  freq_plot <- all_subset_freq_for_plots(freq)
  if (!is.null(freq_plot) && nrow(freq_plot)) {
    utils::write.csv(freq_plot, file.path(out_dir, "all_subsets_frequency_by_bio_trimmed.csv"), row.names = FALSE)
    dropped <- attr(freq_plot, "dropped")
    if (!is.null(dropped) && nrow(dropped)) {
      utils::write.csv(dropped, file.path(out_dir, "all_subsets_dropped_bio_extremes.csv"), row.names = FALSE)
    }
  }

  n_lab <- length(unique(freq$subset_label))
  w_facet <- max(12, 1.1 * n_lab / 2)
  save_gg(plot_all_freq_facet(freq_plot, "All subsets (within-panel %)  JY-NNK / JY-EVNK"),
          file.path(out_dir, "all_subsets_frequency_JY_NNK_vs_JY_EVNK"),
          width = w_facet, height = 5.8)
  save_gg(plot_all_stacked(freq_plot, "Within-panel composition (mean of remaining bio-reps)"),
          file.path(out_dir, "all_subsets_composition_stacked"),
          width = 11, height = 5.6)

  focus <- freq_plot[freq_plot$role == "focus", , drop = FALSE]
  stats_f <- stats[stats$role == "focus", , drop = FALSE]
  if (nrow(focus) > 0) {
    save_gg(plot_all_freq_facet(focus, "Focus subsets only  (P1 T/NK, P2 B, P3 myeloid)"),
            file.path(out_dir, "all_subsets_focus_frequency_JY_NNK_vs_JY_EVNK"),
            width = max(11, 1.15 * length(unique(focus$subset_label)) / 2), height = 5.8)
  }
  if (nrow(stats_f) > 0) {
    save_gg(plot_all_heatmap(stats_f, "Focus subset mean %  (within panel)"),
            file.path(out_dir, "all_subsets_focus_mean_heatmap"),
            width = 6.5, height = max(6, 0.28 * nrow(stats_f) + 2))
    save_gg(plot_all_lollipop(stats_f, "Focus subsets  JY-NNK minus JY-EVNK"),
            file.path(out_dir, "all_subsets_focus_delta_lollipop"),
            width = 8.5, height = max(6, 0.28 * nrow(stats_f) + 2))
  }

  log_msg("All-subset summary: ", out_dir,
          " (", n_ok, " panels, ", n_lab, " subset labels)")
  invisible(list(frequency = freq, stats = stats, out_dir = out_dir))
}

load_jy_engine_functions()

if (!identical(toupper(Sys.getenv("FLOW_ALL_SUBSETS_FROM_PIPELINE", "0")), "1") &&
    !identical(toupper(Sys.getenv("FLOW_FUNCTIONS_ONLY", "0")), "1")) {
  export_all_subsets_analysis(if (exists("result_dir")) result_dir else file.path(getwd(), "results_flow"))
}
