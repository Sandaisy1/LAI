#!/usr/bin/env Rscript
# =============================================================================
# Internation cell immune：His+ 靶细胞 + P1/P2/P3 免疫亚群
#
# 独立于原来的 Flow_dimred_pipeline.R（E:/R/fuction of cell）。
# P1/P2/P3 圈门三步（与免疫降维同一套物理门，分析母群换成 His-FITC+）：
#   1）去除粘连体 → 圈定淋巴细胞（P1 紧淋巴 / P2 宽单核 / P3 不加淋巴门）→ 死活排除
#   2）圈定 His-FITC 阳性细胞群（分析母群，不要求 CD45+）
#   3）在 His+ 上按免疫细胞降维方案做各类亚群（染色以 ICI_flow_panel_map.json 为准）
# P2 是 B 细胞 7 色（CD19/CD27/IgG/IgM + His），没有 IgD/BLIMP/CD40。
#
# 组别：ZZX-EV（EV-1 / EV-2 / EV-3）vs ZZX-H（H-1 / H-2 / H-3），比较 H vs EV，n=3。
#
# 用法：把本方案全部 ICI_* 文件放到 Internation 目录（不要去 fuction of cell 找原流程）
#   setwd("E:/R/Internation cell immune")
#   source("ICI_Flow_dimred_pipeline.R")
# 总结果已出、只补功能状态图：
#   source("ICI_Flow_dimred_functional_state.R")
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

# 只加载本目录的 ICI_flow_engine.R，禁止去 fuction of cell / Flow_dimred_pipeline.R
load_ici_engine <- function() {
  pipe <- "ICI_flow_engine.R"
  cands <- unique(c(
    file.path(ici_script_dir, pipe),
    file.path(getwd(), pipe)
  ))
  bad <- grepl("fuction of cell|function of cell", cands, ignore.case = TRUE)
  cands <- cands[!bad]
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit) || !nzchar(hit)) {
    stop(
      "找不到 ICI_flow_engine.R。靶细胞方案与免疫亚群方案完全独立，",
      "不会去 E:/R/fuction of cell 找 Flow_dimred_pipeline.R。\n",
      "请把这些文件一起放到 E:/R/Internation cell immune：\n",
      "  ICI_Flow_dimred_pipeline.R\n",
      "  ICI_flow_engine.R\n",
      "  ICI_flow_panel_map.json\n",
      "  ICI_Flow_dimred_all_subsets.R\n",
      "  ICI_Flow_dimred_trajectory.R\n",
      "  ICI_Flow_dimred_functional_state.R"
    )
  }
  source(hit, local = FALSE)
  if (exists("log_msg", mode = "function")) log_msg("ICI engine (this scheme only): ", hit)
  invisible(TRUE)
}

load_ici_engine()
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

# 覆盖本方案的数据目录与 ICI 染色表（不读、不改 fuction of cell 的文件）
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
    "^(?:ZZX[_-]?)?(EV|H)[-_ ]?([123])(?:[-_ ]([12]))?[-_ ]+(?:PANEL[-_ ]?)?P?0?([123])[-_ ].*(unmixed|raw)\\.fcs$",
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
  if (is.na(x$panel) || !x$panel %in% c("P1", "P2", "P3")) {
    log_msg("skip Panel ", x$panel, " (Internation sheets are P1/P2/P3): ", basename(path))
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

# QC 用引擎：单细胞 → 淋巴散射（P1 紧 / P2 宽 / P3 无）→ 活细胞，不切 CD45+。
# 分析母群 His-FITC+ 在 load_panel_cells 里圈。下面只保留免疫亚群方案带 CD45 的 QC，供对照测试。
qc_filter_matrix_flow <- function(exprs, names, map, panel_id = NA) {
  keep <- qc_filter_matrix(exprs, names, map, panel_id)
  cd45_idx <- map$channel_index[map$marker == "CD45"]
  if (length(cd45_idx) && !is.na(cd45_idx[1])) {
    cd <- asinh(pmax(as.numeric(exprs[, cd45_idx[1]]), 0) / 150)
    keep <- keep & cd >= qc_cd45_cut(cd)
  }
  keep[is.na(keep)] <- FALSE
  keep
}

# 图上的细胞已经全是 His+，不要再单独标一个 His+ target 亚群。
# 未圈中的 His+ 与免疫降维一样叫 other / dump，不要改成 Target。

# ICI P2 没有 IgD / BLIMP：CD19 母门后用 CD27 × IgM/IgG 对应原 P2 的 naive / unswitched / switched
gate_p2_major_flow <- gate_p2_major
gate_p2_major <- function(mat) {
  mat <- as.matrix(mat)
  n <- nrow(mat)
  out <- rep("Naive_B", n)
  if (n < 40) return(out)
  b <- rep(TRUE, n)
  if ("CD19" %in% colnames(mat)) {
    b_hi <- gate_k2_high(mat, seq_len(n), "CD19", 0.12)
    if (any(b_hi) && mean(b_hi) > 0.08 && mean(b_hi) < 0.97) {
      lo_med <- median(colv(mat, "CD19")[!b_hi], na.rm = TRUE)
      hi_med <- median(colv(mat, "CD19")[b_hi], na.rm = TRUE)
      if (is.finite(lo_med) && is.finite(hi_med) && hi_med > lo_med + 1.0 && lo_med < 0.9) {
        b <- b_hi
        out[!b] <- "other"
      }
    }
  }
  ib <- which(b)
  if (length(ib) >= 40 && "CD27" %in% colnames(mat)) {
    cd27_pos <- p2_pos_mask(colv(mat, "CD27")[ib])
    cd27_pos[!is.finite(colv(mat, "CD27")[ib])] <- FALSE
    igg <- colv(mat, "IgG")[ib]
    igm <- colv(mat, "IgM")[ib]
    igg_hi <- if (any(is.finite(igg))) p2_pos_mask(igg) else rep(FALSE, length(ib))
    igm_hi <- if (any(is.finite(igm))) p2_pos_mask(igm) else rep(FALSE, length(ib))
    igg_hi[!is.finite(igg)] <- FALSE
    igm_hi[!is.finite(igm)] <- FALSE
    lab <- rep("Naive_B", length(ib))
    lab[!cd27_pos] <- "Naive_B"
    lab[cd27_pos & igg_hi] <- "Switched_B"
    lab[cd27_pos & !igg_hi] <- "Unswitched_B"
    out[ib] <- lab
  }
  p2_assign_cd19neg_plasma(mat, out)
}

celltype_label_flow <- celltype_label
celltype_label <- function(lineage, panel_id) {
  lab <- as.character(lineage)
  out <- ifelse(lab %in% c("Target", "His_target", "His+ target"), "other",
                celltype_label_flow(lineage, panel_id))
  unname(out)
}

subset_plot_specs_flow <- subset_plot_specs
subset_plot_specs <- function(panel_id) {
  if (identical(panel_id, "P2")) {
    mk <- function(lineage, x, y, parent, ylab, use_major = FALSE,
                   gate = "box", x_hi = NA, y_hi = NA) {
      list(lineage = lineage, x = x, y = y, parent = parent, ylab = ylab,
           use_major = use_major, gate = gate, x_hi = x_hi, y_hi = y_hi)
    }
    return(list(
      mk("Naive_B", "IgM", "CD27", "CD19", "Naive B in CD19+ (%)", gate = "half_y", y_hi = FALSE),
      mk("Unswitched_B", "IgM", "CD27", "CD19", "Unswitched memory B in CD19+ (%)", gate = "quad", x_hi = TRUE, y_hi = TRUE),
      mk("Switched_B", "IgG", "CD27", "CD19", "Switched memory B in CD19+ (%)", gate = "hi_hi"),
      mk("MZ_B", "IgM", "CD27", "Naive_Unswitched", "MZ B in naive/unswitched (%)", gate = "hi_x", x_hi = TRUE)
    ))
  }
  subset_plot_specs_flow(panel_id)
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
    CD4_naive = pop(CD3 = 3.2, CD4 = 3.0, CD8 = 0.1, CD62L = 3.2, CD44 = 0.3, CD45 = 3.0, His = 3.2),
    CD4_TCM = pop(CD3 = 3.2, CD4 = 3.0, CD62L = 2.8, CD44 = 2.8, CD45 = 3.0, His = 3.2),
    CD4_TEM = pop(CD3 = 3.2, CD4 = 3.0, CD62L = 0.3, CD44 = 3.1, CD45 = 3.0, His = 3.2),
    Treg = pop(CD3 = 3.1, CD4 = 3.0, CD25 = 3.2, CD69 = 0.4, CD45 = 3.0, His = 3.2),
    CD4_act = pop(CD3 = 3.2, CD4 = 3.0, CD69 = 3.1, CD62L = 3.2, CD44 = 0.3, CD45 = 3.0, His = 3.2),
    CD8_naive = pop(CD3 = 3.3, CD8 = 3.1, CD4 = 0.1, CD62L = 3.1, CD44 = 0.3, CD45 = 3.0, His = 3.2),
    CD8_TEM = pop(CD3 = 3.3, CD8 = 3.1, CD44 = 3.0, CD62L = 0.3, CD45 = 3.0, His = 3.2),
    NK = pop(CD3 = 0.1, NKp46 = 3.0, `NK1.1` = 3.1, CD11B = 2.4, CD45 = 3.0, His = 3.2, NKG2D = 2.6),
    NK_act = pop(CD3 = 0.1, NKp46 = 3.0, CD69 = 3.1, CD45 = 3.0, His = 3.2),
    NKT_CD4 = pop(CD3 = 3.0, NKp46 = 2.6, CD4 = 3.0, CD45 = 3.0, His = 3.2),
    Myeloid = pop(CD11B = 3.2, CD3 = 0.1, CD45 = 3.0, His = 3.2),
    Other = pop(His = 3.2, CD45 = 0.3, CD3 = 0.1, CD11B = 0.2, NKp46 = 0.2)
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
    Neutrophil = pop(CD11B = 3.2, LY6G = 3.3, LY6C = 1.5, CD3 = 0.1, CD19 = 0.1, CD45 = 3.0, His = 3.2),
    Mono_Ly6Chi = pop(CD11B = 3.1, LY6C = 3.2, LY6G = 0.2, `F4/80` = 0.5, CD45 = 3.0, His = 3.2),
    Macrophage = pop(CD11B = 3.0, `F4/80` = 3.2, `I-A/I-E` = 1.8, CD45 = 3.0, His = 3.2),
    M2_like = pop(CD11B = 3.0, `F4/80` = 3.0, CD206 = 3.1, CD45 = 3.0, His = 3.2),
    DC = pop(CD11C = 3.2, `I-A/I-E` = 3.3, CD11B = 0.4, CD45 = 3.0, His = 3.2, CD103 = 0.2),
    cDC1 = pop(CD11C = 3.1, CD103 = 3.0, `I-A/I-E` = 3.0, CD11B = 0.4, CD45 = 3.0, His = 3.2),
    Eosinophil = pop(`Siglec-F` = 3.2, CCR3 = 2.8, CD11B = 2.6, LY6G = 0.2, CD45 = 3.0, His = 3.2),
    Mast = pop(FceRI = 3.1, CD11B = 2.8, CD45 = 3.0, His = 3.2),
    Other = pop(His = 3.2, CD45 = 0.3, CD3 = 0.1, CD11B = 0.2, CD19 = 0.1)
  )
}

demo_means_p2 <- function() {
  mk <- dr_marker_names("P2")
  base <- setNames(rep(0.2, length(mk)), mk)
  pop <- function(...) {
    v <- base
    alt <- list(...)
    for (nm in names(alt)) v[[nm]] <- alt[[nm]]
    v
  }
  list(
    Naive_B = pop(CD19 = 3.2, CD27 = 0.3, IgM = 2.4, IgG = 0.2, CD45 = 3.0, His = 3.2),
    Unswitched_B = pop(CD19 = 3.1, CD27 = 3.0, IgM = 2.8, IgG = 0.2, CD45 = 3.0, His = 3.2),
    MZ_B = pop(CD19 = 3.1, CD27 = 0.3, IgM = 3.3, IgG = 0.2, CD45 = 3.0, His = 3.2),
    Switched_B = pop(CD19 = 3.0, CD27 = 2.8, IgG = 3.1, IgM = 0.3, CD45 = 3.0, His = 3.2),
    Other = pop(His = 3.2, CD45 = 0.3, CD19 = 0.2, CD27 = 0.2, IgM = 0.2, IgG = 0.2)
  )
}

demo_props <- function(panel_id, group) {
  ctrl <- identical(as.character(group), flow_ctrl_group)
  if (panel_id == "P1") {
    if (ctrl) {
      return(c(CD4_naive = 0.16, CD4_TCM = 0.08, CD4_TEM = 0.08, Treg = 0.05, CD4_act = 0.06,
               CD8_naive = 0.10, CD8_TEM = 0.08, NK = 0.10, NK_act = 0.04, NKT_CD4 = 0.04,
               Myeloid = 0.06, Other = 0.15))
    }
    return(c(CD4_naive = 0.10, CD4_TCM = 0.07, CD4_TEM = 0.08, Treg = 0.05, CD4_act = 0.08,
             CD8_naive = 0.07, CD8_TEM = 0.08, NK = 0.09, NK_act = 0.06, NKT_CD4 = 0.04,
             Myeloid = 0.05, Other = 0.23))
  }
  if (panel_id == "P2") {
    if (ctrl) {
      return(c(Naive_B = 0.28, Unswitched_B = 0.18, MZ_B = 0.10, Switched_B = 0.22, Other = 0.22))
    }
    return(c(Naive_B = 0.20, Unswitched_B = 0.16, MZ_B = 0.08, Switched_B = 0.24, Other = 0.32))
  }
  if (ctrl) {
    return(c(Neutrophil = 0.14, Mono_Ly6Chi = 0.10, Macrophage = 0.12, M2_like = 0.08,
             DC = 0.08, cDC1 = 0.08, Eosinophil = 0.08, Mast = 0.07, Other = 0.25))
  }
  c(Neutrophil = 0.10, Mono_Ly6Chi = 0.08, Macrophage = 0.10, M2_like = 0.10,
    DC = 0.09, cDC1 = 0.08, Eosinophil = 0.08, Mast = 0.07, Other = 0.30)
}

export_ici_his_stats <- function(cells, panel_id, out_dir) {
  if (!nrow(cells)) return(invisible(NULL))
  his_dir <- file.path(out_dir, "target_His")
  dir.create(his_dir, recursive = TRUE, showWarnings = FALSE)
  lin_tab <- lineage_frequencies(cells)
  lin_tab$his_pos_pct <- 100
  utils::write.csv(lin_tab, file.path(his_dir, paste0(panel_id, "_subset_pct_of_His_parent_by_sample.csv")),
                   row.names = FALSE)
  log_msg(panel_id, " His+ parent subset tables: ", his_dir, " (all plotted cells are His+; no His+ target subset)")
  invisible(lin_tab)
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
# 主流程（P1 / P2 / P3）
# -----------------------------------------------------------------------------
if (identical(toupper(Sys.getenv("FLOW_FUNCTIONS_ONLY", "0")), "1") ||
    identical(toupper(Sys.getenv("ICI_FUNCTIONS_ONLY", "0")), "1")) {
  log_msg("ICI FLOW_FUNCTIONS_ONLY=1: skip analysis")
} else {

log_msg("ICI flow dir: ", project_dir)
log_msg("ICI results: ", result_dir)
log_msg("ICI purpose: QC singlets/lymph/live, then His-FITC+ parent; P1/P2/P3 subsets follow immune dimred gates")
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
  for (pn in c("P1", "P2", "P3")) {
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

panels <- c("P1", "P2", "P3")
summaries <- list()
for (pn in panels) {
  log_msg("ICI ", pn, ": His-FITC+ is the analysis parent; subsets follow immune dimred P1/P2/P3 gates")
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

find_ici_extra <- function(nm) {
  cands <- unique(c(
    file.path(ici_script_dir, nm),
    file.path(getwd(), nm)
  ))
  cands <- cands[!grepl("fuction of cell|function of cell", cands, ignore.case = TRUE)]
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit) || !nzchar(hit)) NA_character_ else hit
}
extra <- find_ici_extra("ICI_Flow_dimred_all_subsets.R")
if (!is.na(extra) && nzchar(extra)) {
  tryCatch({
    Sys.setenv(FLOW_ALL_SUBSETS_FROM_PIPELINE = "1")
    sys.source(extra, envir = .GlobalEnv)
    export_all_subsets_analysis(result_dir)
  }, error = function(e) log_msg("all-subsets summary failed: ", e$message))
  Sys.unsetenv("FLOW_ALL_SUBSETS_FROM_PIPELINE")
}

traj <- find_ici_extra("ICI_Flow_dimred_trajectory.R")
if (!is.na(traj) && nzchar(traj)) {
  tryCatch({
    Sys.setenv(FLOW_TRAJECTORY_FROM_PIPELINE = "1")
    sys.source(traj, envir = .GlobalEnv)
    export_all_panel_trajectories(result_dir)
  }, error = function(e) log_msg("trajectory summary failed: ", e$message))
  Sys.unsetenv("FLOW_TRAJECTORY_FROM_PIPELINE")
}

log_msg("ICI done. Keep ICI_* files in ", ici_primary_data_dir,
        "; do not source Flow_dimred_pipeline.R from fuction of cell.")
}
