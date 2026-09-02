#!/usr/bin/env Rscript
# =============================================================================
# T-T6：免疫细胞亚群降维（比较 T6 vs T；P1 T/NK，P2 B，P3 髓系）
#
# 独立于 Flow_* / ICI_* / JY_* / JZ_*。
# 圈门与原免疫亚群方案一致（单细胞 → P1 紧淋巴 / P2 宽单核 / P3 不加淋巴门 → 活 → CD45+）。
# 三个生物学重复全部保留（n=3），不去极端值。染色见 T_T6_flow_panel_map.json。
# P1 Perforin 是 FITC，P2 IgD 是 FITC，P3 iNOS 是 AF488：都不是 His，不要套 ICI 的 His+ 门。
#
# 组别：T（T-1 / T-2 / T-3）vs T6（T6-1 / T6-2 / T6-3）
# 文件：T-1_P1.fcs / T-1_P1_unmixed.fcs / T6-2_P3_unmixed.fcs
#
# 用法：把本方案全部 T_T6_* 文件放到数据目录
#   setwd("E:/R/flow J")
#   source("T_T6_Flow_dimred_pipeline.R")
# 总结果已经出来、只补各亚群活化 / 抑制 / 耗竭（不重跑 UMAP）：
#   source("T_T6_Flow_dimred_functional_state.R")
# 无 FCS 时可 Sys.setenv(FLOW_DEMO = "1")
#
# 结果：E:/R/flow J/results_flow/
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANGUAGE = "en")

tt6_get_script_dir <- function() {
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

tt6_script_dir <- tt6_get_script_dir()
tt6_primary_data_dir <- "E:/R/flow J"

tt6_keep_cand <- function(p) {
  s <- gsub("\\", "/", as.character(p), fixed = TRUE)
  if (grepl("flow J|flowJ|flow_J|/flow-j", s, ignore.case = TRUE)) return(TRUE)
  if (grepl("cell-ljy|cell-wjz", s, ignore.case = TRUE)) return(FALSE)
  if (grepl("Internation cell immune", s, ignore.case = TRUE)) return(FALSE)
  if (grepl("fuction of cell|function of cell", s, ignore.case = TRUE)) return(FALSE)
  TRUE
}

# 只加载本目录的 T_T6_flow_engine.R，禁止去 fuction of cell / Internation 找原流程
load_tt6_engine <- function() {
  pipe <- "T_T6_flow_engine.R"
  cands <- unique(c(
    file.path(tt6_script_dir, pipe),
    file.path(getwd(), pipe)
  ))
  cands <- cands[vapply(cands, tt6_keep_cand, logical(1))]
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit) || !nzchar(hit)) {
    stop(
      "找不到 T_T6_flow_engine.R。T-T6 方案与免疫亚群 / 靶细胞 / JY / JZ 方案完全独立，",
      "不会去 E:/R/fuction of cell 找 Flow_dimred_pipeline.R，也不会用 ICI_* / JY_* / JZ_*。\n",
      "请把这些文件一起放到 E:/R/flow J：\n",
      "  T_T6_Flow_dimred_pipeline.R\n",
      "  T_T6_flow_engine.R\n",
      "  T_T6_flow_panel_map.json\n",
      "  T_T6_Flow_dimred_all_subsets.R\n",
      "  T_T6_Flow_dimred_trajectory.R\n",
      "  T_T6_Flow_dimred_functional_state.R"
    )
  }
  source(hit, local = FALSE)
  if (exists("log_msg", mode = "function")) log_msg("T-T6 engine (this scheme only): ", hit)
  invisible(TRUE)
}

load_tt6_engine()
# T-T6：三个生物学重复全部保留，统计 n=3，不去极端值
flow_trim_bio_extremes <- FALSE

find_tt6_panel_map <- function() {
  names <- c("T_T6_flow_panel_map.json", "T_T6_flow_panel_map.json.txt")
  dirs <- unique(c(
    getwd(),
    tt6_script_dir,
    tt6_primary_data_dir,
    "E:\\R\\flow J"
  ))
  dirs <- dirs[nzchar(dirs)]
  dirs <- dirs[vapply(dirs, tt6_keep_cand, logical(1))]
  for (d in dirs) {
    if (!dir.exists(d)) next
    for (nm in names) {
      p <- file.path(d, nm)
      if (file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
  }
  stop("找不到 T_T6_flow_panel_map.json（应与 T_T6_Flow_dimred_pipeline.R 同目录，或放在 E:/R/flow J）")
}

resolve_tt6_dir <- function() {
  env_dir <- Sys.getenv("T_T6_FLOW_DIR", unset = "")
  if (!nzchar(env_dir)) env_dir <- Sys.getenv("TT6_FLOW_DIR", unset = "")
  preferred <- c(
    env_dir,
    tt6_primary_data_dir,
    "E:\\R\\flow J",
    "E:/R/flow J",
    getwd()
  )
  preferred <- unique(preferred[nzchar(preferred)])
  preferred <- preferred[vapply(preferred, tt6_keep_cand, logical(1))]
  for (d in preferred) {
    if (dir.exists(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  tt6_script_dir
}

# 覆盖本方案的数据目录与染色表（不读、不改 fuction of cell 或 Internation 的文件）
panel_map <- read_panel_map_file(find_tt6_panel_map())
log_msg("T-T6 panel map: ", find_tt6_panel_map())
if (!is.null(panel_map$qc$asinh_cofactor)) asinh_cofactor <- as.numeric(panel_map$qc$asinh_cofactor)
flow_ctrl_group <- "T"
flow_trt_group <- "T6"
flow_group_levels <- c("T", "T6")
if (!is.null(panel_map$groups) && length(panel_map$groups) >= 2) {
  flow_ctrl_group <- as.character(panel_map$groups[[1]])
  flow_trt_group <- as.character(panel_map$groups[[2]])
  flow_group_levels <- c(flow_ctrl_group, flow_trt_group)
}
flow_cohort <- "T-T6"
pal_group <- setNames(c("#1A1A1A", "#E31A1C"), flow_group_levels)
pal_group_shape <- setNames(c(16, 15), flow_group_levels)

flow_primary_data_dir <- tt6_primary_data_dir
project_dir <- resolve_tt6_dir()
result_dir <- file.path(project_dir, "results_flow")
log_dir <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

find_tt6_extra <- function(nm) {
  cands <- unique(c(
    file.path(tt6_script_dir, nm),
    file.path(getwd(), nm)
  ))
  cands <- cands[vapply(cands, tt6_keep_cand, logical(1))]
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit) || !nzchar(hit)) NA_character_ else hit
}

# -----------------------------------------------------------------------------
# 主流程（P1 / P2 / P3）
# -----------------------------------------------------------------------------
if (identical(toupper(Sys.getenv("FLOW_FUNCTIONS_ONLY", "0")), "1") ||
    identical(toupper(Sys.getenv("T_T6_FUNCTIONS_ONLY", "0")), "1") ||
    identical(toupper(Sys.getenv("TT6_FUNCTIONS_ONLY", "0")), "1")) {
  log_msg("T-T6 FLOW_FUNCTIONS_ONLY=1: skip analysis")
} else {

log_msg("T-T6 flow dir: ", project_dir)
log_msg("T-T6 results: ", result_dir)
log_msg("T-T6 purpose: immune-subset dimred, comparison T6 vs T (same gates as original P1/P2/P3)")
log_msg("T-T6 stats: keep all 3 bio-reps per group (n=3); do not drop an extreme")
file_tab <- list_unmixed_files(project_dir)
use_demo <- demo_flag
if (nrow(file_tab) == 0) {
  if (demo_flag) {
    log_msg("No unmixed FCS; FLOW_DEMO=1 -> synthetic data for plot export")
    use_demo <- TRUE
  } else {
    log_msg("No *_unmixed.fcs in ", project_dir)
    log_msg("Put T / T6 files (T-1_P1, T6-1_P1_unmixed.fcs, …) in ", tt6_primary_data_dir,
            ", or set FLOW_DEMO=1")
    use_demo <- TRUE
    log_msg("Auto-fallback to DEMO so the script can still export figure templates")
  }
} else {
  log_msg("Found unmixed files:\n", paste(file_tab$file, collapse = "\n"))
  for (pn in c("P1", "P2", "P3")) {
    n_pn <- sum(file_tab$panel == pn)
    log_msg(pn, " files parsed: ", n_pn)
    if (n_pn == 0) {
      log_msg(pn, " missing. Put files named like T-1_", pn, ".fcs or T6-1_", pn, "_unmixed.fcs")
    }
  }
  if (any(file_tab$group == "T" & grepl("T6", file_tab$file, ignore.case = TRUE))) {
    stop("Filename parser classified a T6 file as T; refusing to continue")
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
  summaries[[pn]] <- tryCatch(
    analyze_one_panel(pn, file_tab, use_demo),
    error = function(e) {
      log_msg("Panel ", pn, " failed: ", e$message)
      NULL
    }
  )
}

sum_path <- file.path(result_dir, "T6_vs_T_lineage_stats_all_panels.csv")
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

extra <- find_tt6_extra("T_T6_Flow_dimred_all_subsets.R")
if (!is.na(extra) && nzchar(extra)) {
  tryCatch({
    Sys.setenv(FLOW_ALL_SUBSETS_FROM_PIPELINE = "1")
    sys.source(extra, envir = .GlobalEnv)
    export_all_subsets_analysis(result_dir)
  }, error = function(e) log_msg("all-subsets summary failed: ", e$message))
  Sys.unsetenv("FLOW_ALL_SUBSETS_FROM_PIPELINE")
}

traj <- find_tt6_extra("T_T6_Flow_dimred_trajectory.R")
if (!is.na(traj) && nzchar(traj)) {
  tryCatch({
    Sys.setenv(FLOW_TRAJECTORY_FROM_PIPELINE = "1")
    sys.source(traj, envir = .GlobalEnv)
    export_all_panel_trajectories(result_dir)
  }, error = function(e) log_msg("trajectory summary failed: ", e$message))
  Sys.unsetenv("FLOW_TRAJECTORY_FROM_PIPELINE")
}

log_msg("T-T6 done. Keep T_T6_* files in ", tt6_primary_data_dir,
        "; do not source Flow_*, ICI_*, JY_*, or JZ_*.")
invisible(summaries)

}
