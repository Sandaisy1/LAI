#!/usr/bin/env Rscript
# =============================================================================
# Internation cell immune：His+ 靶细胞 + P1/P3 免疫亚群
#
# 独立于原来的 Flow_dimred_pipeline.R（E:/R/fuction of cell）。
# 圈门和分析思路相同，染色按 ICI_flow_panel_map.json（图 1–3）。
# 没有 Panel 2。主目的是找 His-FITC 靶细胞，再看它们落在哪些免疫亚群里。
#
# 组别：ZZX-EV（EV-1 / EV-2 / EV-3）vs ZZX-H（H-1 / H-2 / H-3），比较 H vs EV，n=3。
#
# 用法：
#   setwd("E:/R/Internation cell immune")
#   source("ICI_Flow_dimred_pipeline.R")
# 无 FCS 时可 Sys.setenv(FLOW_DEMO = "1")
#
# 结果：E:/R/Internation cell immune/results_flow/
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANGUAGE = "en")

ici_get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)))
  }
  ofile <- NULL
  n <- sys.nframe()
  if (n > 0) {
    for (i in n:1) {
      ofile <- sys.frame(i)$ofile
      if (!is.null(ofile) && nzchar(ofile)) break
    }
  }
  if (!is.null(ofile) && nzchar(ofile)) {
    return(dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

ici_script_dir <- ici_get_script_dir()
ici_primary_data_dir <- "E:/R/Internation cell immune"

load_original_flow_functions <- function() {
  if (exists("hierarchical_gate", mode = "function") && exists("analyze_one_panel", mode = "function")) {
    return(invisible(TRUE))
  }
  pipe <- "Flow_dimred_pipeline.R"
  cands <- c(
    file.path(ici_script_dir, pipe),
    file.path(getwd(), pipe),
    pipe
  )
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit) || !nzchar(hit)) stop("找不到 Flow_dimred_pipeline.R；ICI 脚本要和原流程放在同一仓库")
  old <- Sys.getenv("FLOW_FUNCTIONS_ONLY", unset = NA)
  Sys.setenv(FLOW_FUNCTIONS_ONLY = "1")
  source(hit, local = FALSE)
  if (is.na(old)) Sys.unsetenv("FLOW_FUNCTIONS_ONLY") else Sys.setenv(FLOW_FUNCTIONS_ONLY = old)
  invisible(TRUE)
}

load_original_flow_functions()
# ICI 没有技术重复、n=3，不要套用免疫亚群分析里「去掉一个极端生物学重复」的规则
flow_trim_bio_extremes <- FALSE

find_ici_panel_map <- function() {
  names <- c("ICI_flow_panel_map.json", "ICI_flow_panel_map.json.txt")
  dirs <- unique(c(
    getwd(),
    ici_script_dir,
    ici_primary_data_dir,
    "E:\\R\\Internation cell immune"
  ))
  dirs <- dirs[nzchar(dirs)]
  for (d in dirs) {
    if (!dir.exists(d)) next
    for (nm in names) {
      p <- file.path(d, nm)
      if (file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
  }
  stop("找不到 ICI_flow_panel_map.json（应与 ICI_Flow_dimred_pipeline.R 同目录，或放在 E:/R/Internation cell immune）")
}

resolve_ici_dir <- function() {
  env_dir <- Sys.getenv("ICI_FLOW_DIR", unset = "")
  preferred <- c(
    env_dir,
    ici_primary_data_dir,
    "E:\\R\\Internation cell immune",
    getwd()
  )
  preferred <- unique(preferred[nzchar(preferred)])
  for (d in preferred) {
    if (dir.exists(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  ici_script_dir
}

# 覆盖原流程的数据目录与 panel 表；不改原文件
panel_map <- read_panel_map_file(find_ici_panel_map())
log_msg("ICI panel map: ", find_ici_panel_map())
if (!is.null(panel_map$qc$asinh_cofactor)) asinh_cofactor <- as.numeric(panel_map$qc$asinh_cofactor)
flow_ctrl_group <- "EV"
flow_trt_group <- "H"
flow_group_levels <- c("EV", "H")
flow_cohort <- "ZZX"
pal_group <- setNames(c("#1A1A1A", "#E31A1C"), flow_group_levels)

flow_primary_data_dir <- ici_primary_data_dir
project_dir <- resolve_ici_dir()
result_dir <- file.path(project_dir, "results_flow")
log_dir <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# 样品写成 EV-1 / H-2，匹配 ZZX-EV / ZZX-H 的三次重复（无技术重复时）
parse_fcs_filename_flow <- parse_fcs_filename

parse_ici_filename_fallback <- function(path) {
  b <- basename(path)
  m <- regexec(
    "^(?:ZZX[_-]?)?(EV|H)[-_ ]?([123])(?:[-_ ]([12]))?[-_ ]+(?:PANEL[-_ ]?)?P?0?([13])[-_ ].*(unmixed|raw)\\.fcs$",
    b,
    ignore.case = TRUE
  )
  hit <- regmatches(b, m)[[1]]
  if (length(hit) < 5L) return(NULL)
  grp <- toupper(hit[2])
  bio <- hit[3]
  tech <- if (nzchar(hit[4])) hit[4] else NA_character_
  panel_n <- hit[5]
  kind <- if (length(hit) >= 6L) tolower(hit[6]) else "unmixed"
  bio_sample <- paste0(grp, "-", bio)
  list(
    file = b,
    path = path,
    cohort = flow_cohort,
    group = grp,
    replicate = bio,
    tech_rep = tech,
    bio_sample = bio_sample,
    sample = if (!is.na(tech) && nzchar(tech)) paste0(bio_sample, "-", tech) else bio_sample,
    panel = paste0("P", panel_n),
    kind = kind
  )
}

parse_fcs_filename <- function(path) {
  x <- parse_fcs_filename_flow(path)
  if (is.null(x)) x <- parse_ici_filename_fallback(path)
  if (is.null(x)) return(NULL)
  if (is.na(x$panel) || !x$panel %in% c("P1", "P3")) {
    log_msg("skip Panel ", x$panel, " (Internation sheets are P1/P3 only): ", basename(path))
    return(NULL)
  }
  x$bio_sample <- paste0(x$group, "-", x$replicate)
  if (is.na(x$tech_rep) || !nzchar(as.character(x$tech_rep))) {
    x$sample <- x$bio_sample
  } else {
    x$sample <- paste0(x$bio_sample, "-", x$tech_rep)
  }
  x
}

# 找靶细胞：QC 只去双联体和死细胞，不要先切 CD45+ / 小淋巴门把 His+ CD45- 丢掉
qc_filter_matrix_flow <- qc_filter_matrix
qc_filter_matrix <- function(exprs, names, map, panel_id = NA) {
  keep <- rep(TRUE, nrow(exprs))
  fsc <- grep("^FSC-A$|^FSC_A$|^FSC\\.A$", names, ignore.case = TRUE)
  ssc <- grep("^SSC-A$|^SSC_A$|^SSC\\.A$", names, ignore.case = TRUE)
  if (length(fsc) && length(ssc)) {
    fs <- exprs[, fsc[1]]
    ss <- exprs[, ssc[1]]
    keep <- keep & fs > quantile(fs, 0.01, na.rm = TRUE) & ss > quantile(ss, 0.01, na.rm = TRUE)
    keep <- keep & ss <= quantile(ss, 0.97, na.rm = TRUE)
    keep <- keep & fs >= quantile(fs, 0.02, na.rm = TRUE) & fs <= quantile(fs, 0.99, na.rm = TRUE)
  }
  fsch <- grep("^FSC-H$|^FSC_H$", names, ignore.case = TRUE)
  if (length(fsc) && length(fsch)) {
    ratio <- exprs[, fsc[1]] / pmax(exprs[, fsch[1]], 1)
    keep <- keep & ratio > 0.5 & ratio < 2
  }
  fscw <- grep("^FSC-W$|^FSC_W$", names, ignore.case = TRUE)
  if (length(fscw) && any(keep)) {
    fw <- exprs[, fscw[1]]
    keep <- keep & fw <= quantile(fw[keep], 0.95, na.rm = TRUE)
  }
  ld_idx <- map$channel_index[map$marker == "L/D"]
  if (length(ld_idx) && !is.na(ld_idx[1])) {
    ld <- exprs[, ld_idx[1]]
    keep <- keep & ld <= quantile(ld, 0.95, na.rm = TRUE)
  }
  keep[is.na(keep)] <- FALSE
  keep
}

# His+ CD45- = 靶细胞；His+ CD45+ 仍走原来的免疫亚群，只打 his_pos 标记
ici_apply_his <- function(mat, h) {
  his <- colv(mat, "His")
  cd45 <- colv(mat, "CD45")
  n <- length(h$subset)
  if (length(his) != n || !any(is.finite(his))) return(h)
  cut_his <- axis_pos_cut(his)
  his_pos <- is.finite(his) & is.finite(cut_his) & his >= cut_his
  if (any(is.finite(cd45))) {
    cut45 <- axis_pos_cut(cd45)
    cd45_pos <- is.finite(cd45) & is.finite(cut45) & cd45 >= cut45
  } else {
    cd45_pos <- rep(TRUE, n)
  }
  tgt <- his_pos & !cd45_pos
  if (any(tgt)) {
    h$major[tgt] <- "Target"
    h$subset[tgt] <- "Target"
  }
  h
}

hierarchical_gate_flow <- hierarchical_gate
hierarchical_gate <- function(mat, panel_id) {
  mat <- as.matrix(mat)
  n <- nrow(mat)
  peek <- ici_apply_his(mat, list(major = rep("other", n), subset = rep("other", n)))
  tgt <- peek$major == "Target"
  tgt[is.na(tgt)] <- FALSE
  out <- peek
  rest <- !tgt
  if (any(rest)) {
    h <- hierarchical_gate_flow(mat[rest, , drop = FALSE], panel_id)
    out$major[rest] <- h$major
    out$subset[rest] <- h$subset
  }
  out
}

pal_celltype <- c(pal_celltype,
  "His+ target" = "#00ACC1",
  "Target" = "#00ACC1"
)

celltype_label_flow <- celltype_label
celltype_label <- function(lineage, panel_id) {
  lab <- as.character(lineage)
  out <- ifelse(lab %in% c("Target", "His_target"), "His+ target", celltype_label_flow(lineage, panel_id))
  unname(out)
}

parent_mask_flow <- parent_mask
parent_mask <- function(cells, parent) {
  if (identical(as.character(parent), "Target") || identical(as.character(parent), "His")) {
    n <- nrow(cells)
    lin <- if ("lineage" %in% names(cells)) as.character(cells$lineage) else rep("", n)
    lin[is.na(lin)] <- ""
    cl <- if ("cluster_lineage" %in% names(cells)) as.character(cells$cluster_lineage) else rep("", n)
    cl[is.na(cl)] <- ""
    return(lin %in% c("Target", "His_target") | cl == "Target")
  }
  parent_mask_flow(cells, parent)
}

subset_plot_specs_flow <- subset_plot_specs
subset_plot_specs <- function(panel_id) {
  mk <- function(lineage, x, y, parent, ylab, use_major = FALSE,
                 gate = "box", x_hi = NA, y_hi = NA) {
    list(lineage = lineage, x = x, y = y, parent = parent, ylab = ylab,
         use_major = use_major, gate = gate, x_hi = x_hi, y_hi = y_hi)
  }
  his <- list(
    mk("Target", "His", "CD45", "all", "His+ CD45- target (%)", gate = "quad", x_hi = TRUE, y_hi = FALSE)
  )
  orig <- tryCatch(subset_plot_specs_flow(panel_id), error = function(e) list())
  c(his, orig)
}

demo_means_p1 <- function() {
  mk <- dr_marker_names("P1")
  base <- setNames(rep(0.2, length(mk)), mk)
  pop <- function(...) {
    v <- base
    alt <- list(...)
    for (nm in names(alt)) v[[nm]] <- alt[[nm]]
    v
  }
  list(
    CD4_naive = pop(CD3 = 3.2, CD4 = 3.0, CD8 = 0.1, CD62L = 3.2, CD44 = 0.3, CD45 = 3.0, His = 0.3),
    CD4_TCM = pop(CD3 = 3.2, CD4 = 3.0, CD62L = 2.8, CD44 = 2.8, CD45 = 3.0, His = 0.3),
    CD4_TEM = pop(CD3 = 3.2, CD4 = 3.0, CD62L = 0.3, CD44 = 3.1, CD45 = 3.0, His = 0.3),
    Treg = pop(CD3 = 3.1, CD4 = 3.0, CD25 = 3.2, CD69 = 0.4, CD45 = 3.0, His = 0.3),
    CD4_act = pop(CD3 = 3.2, CD4 = 3.0, CD69 = 3.1, CD62L = 3.2, CD44 = 0.3, CD45 = 3.0, His = 0.3),
    CD8_naive = pop(CD3 = 3.3, CD8 = 3.1, CD4 = 0.1, CD62L = 3.1, CD44 = 0.3, CD45 = 3.0, His = 0.3),
    CD8_TEM = pop(CD3 = 3.3, CD8 = 3.1, CD44 = 3.0, CD62L = 0.3, CD45 = 3.0, His = 0.3),
    NK = pop(CD3 = 0.1, NKp46 = 3.0, `NK1.1` = 3.1, CD11B = 2.4, CD45 = 3.0, His = 0.3, NKG2D = 2.6),
    NK_act = pop(CD3 = 0.1, NKp46 = 3.0, CD69 = 3.1, CD45 = 3.0, His = 0.3),
    NKT_CD4 = pop(CD3 = 3.0, NKp46 = 2.6, CD4 = 3.0, CD45 = 3.0, His = 0.3),
    Myeloid = pop(CD11B = 3.2, CD3 = 0.1, CD45 = 3.0, His = 0.3),
    Target = pop(His = 3.2, CD45 = 0.3, CD3 = 0.1, CD11B = 0.2, NKp46 = 0.2)
  )
}

demo_means_p3 <- function() {
  mk <- dr_marker_names("P3")
  base <- setNames(rep(0.2, length(mk)), mk)
  pop <- function(...) {
    v <- base
    alt <- list(...)
    for (nm in names(alt)) v[[nm]] <- alt[[nm]]
    v
  }
  list(
    Neutrophil = pop(CD11B = 3.2, LY6G = 3.3, LY6C = 1.5, CD3 = 0.1, CD19 = 0.1, CD45 = 3.0, His = 0.3),
    Mono_Ly6Chi = pop(CD11B = 3.1, LY6C = 3.2, LY6G = 0.2, `F4/80` = 0.5, CD45 = 3.0, His = 0.3),
    Macrophage = pop(CD11B = 3.0, `F4/80` = 3.2, `I-A/I-E` = 1.8, CD45 = 3.0, His = 0.3),
    M2_like = pop(CD11B = 3.0, `F4/80` = 3.0, CD206 = 3.1, CD45 = 3.0, His = 0.3),
    DC = pop(CD11C = 3.2, `I-A/I-E` = 3.3, CD11B = 0.4, CD45 = 3.0, His = 0.3, CD103 = 0.2),
    cDC1 = pop(CD11C = 3.1, CD103 = 3.0, `I-A/I-E` = 3.0, CD11B = 0.4, CD45 = 3.0, His = 0.3),
    Eosinophil = pop(`Siglec-F` = 3.2, CCR3 = 2.8, CD11B = 2.6, LY6G = 0.2, CD45 = 3.0, His = 0.3),
    Mast = pop(FceRI = 3.1, CD11B = 2.8, CD45 = 3.0, His = 0.3),
    Target = pop(His = 3.2, CD45 = 0.3, CD3 = 0.1, CD11B = 0.2, CD19 = 0.1)
  )
}

demo_props <- function(panel_id, group) {
  ctrl <- identical(as.character(group), flow_ctrl_group)
  if (panel_id == "P1") {
    if (ctrl) {
      return(c(CD4_naive = 0.16, CD4_TCM = 0.08, CD4_TEM = 0.08, Treg = 0.05, CD4_act = 0.06,
               CD8_naive = 0.10, CD8_TEM = 0.08, NK = 0.10, NK_act = 0.04, NKT_CD4 = 0.04,
               Myeloid = 0.06, Target = 0.15))
    }
    return(c(CD4_naive = 0.10, CD4_TCM = 0.07, CD4_TEM = 0.08, Treg = 0.05, CD4_act = 0.08,
             CD8_naive = 0.07, CD8_TEM = 0.08, NK = 0.09, NK_act = 0.06, NKT_CD4 = 0.04,
             Myeloid = 0.05, Target = 0.23))
  }
  if (ctrl) {
    return(c(Neutrophil = 0.14, Mono_Ly6Chi = 0.10, Macrophage = 0.12, M2_like = 0.08,
             DC = 0.08, cDC1 = 0.08, Eosinophil = 0.08, Mast = 0.07, Target = 0.25))
  }
  c(Neutrophil = 0.10, Mono_Ly6Chi = 0.08, Macrophage = 0.10, M2_like = 0.10,
    DC = 0.09, cDC1 = 0.08, Eosinophil = 0.08, Mast = 0.07, Target = 0.30)
}

ici_his_mask <- function(cells) {
  lin <- if ("lineage" %in% names(cells)) as.character(cells$lineage) else rep("", nrow(cells))
  lin[is.na(lin)] <- ""
  tgt <- lin == "Target"
  his <- if ("His" %in% names(cells)) as.numeric(cells$His) else rep(NA_real_, nrow(cells))
  if (any(is.finite(his)) && !any(tgt)) {
    cut_his <- axis_pos_cut(his)
    tgt <- is.finite(his) & is.finite(cut_his) & his >= cut_his
    cd45 <- if ("CD45" %in% names(cells)) as.numeric(cells$CD45) else rep(NA_real_, nrow(cells))
    if (any(is.finite(cd45))) {
      cut45 <- axis_pos_cut(cd45)
      tgt <- tgt & !(is.finite(cd45) & is.finite(cut45) & cd45 >= cut45)
    }
  }
  tgt[is.na(tgt)] <- FALSE
  tgt
}

plot_ici_his_split <- function(cells, x, y, xlab, ylab, title, pal) {
  plot_df <- cells
  plot_df$group <- factor(plot_df$group, levels = flow_group_levels)
  levs <- intersect(names(pal), unique(as.character(plot_df$celltype)))
  plot_df$celltype <- factor(plot_df$celltype, levels = levs)
  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = celltype)) +
    ggplot2::geom_point(size = 0.72, alpha = 0.92, stroke = 0) +
    ggplot2::facet_wrap(~group, ncol = 2, scales = "fixed") +
    ggplot2::scale_color_manual(values = pal[levs], drop = FALSE) +
    ggplot2::labs(title = title, x = xlab, y = ylab, color = NULL) +
    ggplot2::coord_fixed(ratio = 1, clip = "off") +
    theme_split_dr() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(1.15, "lines"),
      legend.key = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 18, 10, 10)
    )
}

export_ici_his_stats <- function(cells, panel_id, out_dir) {
  if (!nrow(cells) || !("His" %in% names(cells))) return(invisible(NULL))
  tgt <- ici_his_mask(cells)
  his <- as.numeric(cells$His)
  cut_his <- axis_pos_cut(his)
  his_pos <- is.finite(his) & is.finite(cut_his) & his >= cut_his
  his_pos[is.na(his_pos)] <- FALSE
  smp <- unique(as.character(cells$sample))
  rows <- list()
  for (s in smp) {
    ii <- which(as.character(cells$sample) == s)
    if (!length(ii)) next
    grp <- as.character(cells$group[ii[1]])
    bio <- if ("bio_sample" %in% names(cells)) as.character(cells$bio_sample[ii[1]]) else s
    n <- length(ii)
    rows[[length(rows) + 1]] <- data.frame(
      sample = s, bio_sample = bio, group = grp, panel = panel_id,
      metric = "His_pos_of_live", n_cells = n,
      percent = 100 * mean(his_pos[ii]),
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1]] <- data.frame(
      sample = s, bio_sample = bio, group = grp, panel = panel_id,
      metric = "His_target_CD45neg", n_cells = n,
      percent = 100 * mean(tgt[ii]),
      stringsAsFactors = FALSE
    )
  }
  tab <- do.call(rbind, rows)
  his_dir <- file.path(out_dir, "target_His")
  dir.create(his_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(tab, file.path(his_dir, paste0(panel_id, "_His_target_by_sample.csv")), row.names = FALSE)
  stats <- compare_group_freq(tab, "metric")
  utils::write.csv(stats, file.path(his_dir, paste0(panel_id, "_His_target_H_vs_EV_stats.csv")), row.names = FALSE)
  lin_tab <- lineage_frequencies(cells)
  lin_tab$his_pos_pct <- NA_real_
  for (i in seq_len(nrow(lin_tab))) {
    hit <- as.character(cells$sample) == lin_tab$sample[i] & as.character(cells$lineage) == lin_tab$lineage[i]
    if (any(hit)) lin_tab$his_pos_pct[i] <- 100 * mean(his_pos[hit])
  }
  utils::write.csv(lin_tab, file.path(his_dir, paste0(panel_id, "_His_within_subset_by_sample.csv")), row.names = FALSE)
  log_msg(panel_id, " His+ target tables: ", his_dir)
  pal_his <- c("His+ target" = "#00ACC1", "His+ CD45+" = "#F9A825", "His-" = "#B0B0B0")
  if ("His" %in% names(cells) && "tSNE1" %in% names(cells)) {
    plot_df <- cells
    plot_df$celltype <- ifelse(tgt, "His+ target", ifelse(his_pos, "His+ CD45+", "His-"))
    n_keys <- length(unique(as.character(plot_df$celltype)))
    save_split_dr(
      plot_ici_his_split(plot_df, "tSNE1", "tSNE2", "tSNE-1", "tSNE-2",
                         paste(panel_id, "  EV | H  His+ target"), pal_his),
      file.path(his_dir, paste0(panel_id, "_H_vs_EV_tSNE_His_target")),
      n_keys
    )
    if ("UMAP1" %in% names(cells)) {
      save_split_dr(
        plot_ici_his_split(plot_df, "UMAP1", "UMAP2", "UMAP-1", "UMAP-2",
                           paste(panel_id, "  EV | H  His+ target"), pal_his),
        file.path(his_dir, paste0(panel_id, "_H_vs_EV_UMAP_His_target")),
        n_keys
      )
    }
  }
  invisible(tab)
}

export_dimred_plots_flow <- export_dimred_plots
export_dimred_plots <- function(cells, med, annot, freq_df, panel_id, out_dir,
                                umap_is_pca = FALSE, tsne_is_pca = FALSE) {
  export_dimred_plots_flow(cells, med, annot, freq_df, panel_id, out_dir, umap_is_pca, tsne_is_pca)
  tryCatch(
    export_ici_his_stats(cells, panel_id, out_dir),
    error = function(e) log_msg(panel_id, " His target export failed: ", e$message)
  )
}

# -----------------------------------------------------------------------------
# 主流程（P1 / P3；不要跑原来的 P2）
# -----------------------------------------------------------------------------
if (identical(toupper(Sys.getenv("FLOW_FUNCTIONS_ONLY", "0")), "1") ||
    identical(toupper(Sys.getenv("ICI_FUNCTIONS_ONLY", "0")), "1")) {
  log_msg("ICI FLOW_FUNCTIONS_ONLY=1: skip analysis")
} else {

log_msg("ICI flow dir: ", project_dir)
log_msg("ICI results: ", result_dir)
log_msg("ICI purpose: His+ target cells; immune gates reuse original P1/P3 logic with this sheet")
file_tab <- list_unmixed_files(project_dir)
use_demo <- demo_flag
if (nrow(file_tab) == 0) {
  log_msg("No *_unmixed.fcs in ", project_dir)
  log_msg("Put ZZX-EV / ZZX-H unmixed files (EV-1_P1, H-2_P3, …) in ", ici_primary_data_dir,
          ", or set FLOW_DEMO=1")
  use_demo <- TRUE
  log_msg("Auto-fallback to DEMO so the script can still export figure templates")
} else {
  log_msg("Found unmixed files:\n", paste(file_tab$file, collapse = "\n"))
  for (pn in c("P1", "P3")) {
    n_pn <- sum(file_tab$panel == pn)
    log_msg(pn, " files parsed: ", n_pn)
  }
  if (any(file_tab$group == "H" & grepl("(^|[_-])EV", file_tab$file, ignore.case = TRUE))) {
    stop("Filename parser classified an EV/ZZX-EV file as H; refusing to continue")
  }
  ensure_flowcore()
}

if (use_demo) {
  writeLines("DEMO / synthetic cells. Do not use as biological results.",
             file.path(log_dir, "DEMO_WARNING.txt"))
}

panels <- c("P1", "P3")
summaries <- list()
for (pn in panels) {
  log_msg("ICI ", pn, ": His+ CD45- = target; remaining CD45+ use original lineage gates")
  log_msg(pn, " dimred: *_major_split by major class; *_lineage_split all fine subsets on the same embedding; dimred_by_major/ per-class; missing stain channels skip that item only")
  summaries[[pn]] <- tryCatch(
    analyze_one_panel(pn, file_tab, use_demo),
    error = function(e) {
      log_msg("Panel ", pn, " failed: ", e$message)
      NULL
    }
  )
}

sum_path <- file.path(result_dir, "H_vs_EV_lineage_stats_all_panels.csv")
lin_rows <- lapply(names(summaries), function(pn) {
  x <- summaries[[pn]]
  if (is.null(x) || is.null(x$stats_lineage)) return(NULL)
  cbind(panel = pn, x$stats_lineage)
})
lin_rows <- Filter(Negate(is.null), lin_rows)
if (length(lin_rows) > 0) {
  utils::write.csv(do.call(rbind, lin_rows), sum_path, row.names = FALSE)
  log_msg("Summary table: ", sum_path)
}

extra <- file.path(ici_script_dir, "Flow_dimred_all_subsets.R")
if (file.exists(extra)) {
  tryCatch({
    Sys.setenv(FLOW_ALL_SUBSETS_FROM_PIPELINE = "1")
    sys.source(extra, envir = .GlobalEnv)
    export_all_subsets_analysis(result_dir)
  }, error = function(e) log_msg("all-subsets summary failed: ", e$message))
  Sys.unsetenv("FLOW_ALL_SUBSETS_FROM_PIPELINE")
}

traj <- file.path(ici_script_dir, "Flow_dimred_trajectory.R")
if (file.exists(traj)) {
  tryCatch({
    Sys.setenv(FLOW_TRAJECTORY_FROM_PIPELINE = "1")
    sys.source(traj, envir = .GlobalEnv)
    export_all_panel_trajectories(result_dir)
  }, error = function(e) log_msg("trajectory summary failed: ", e$message))
  Sys.unsetenv("FLOW_TRAJECTORY_FROM_PIPELINE")
}

log_msg("ICI done. Copy ICI_Flow_dimred_pipeline.R + ICI_flow_panel_map.json into ",
        ici_primary_data_dir, " and source there.")
}
