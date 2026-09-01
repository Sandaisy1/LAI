#!/usr/bin/env Rscript
# =============================================================================
# JY：免疫细胞亚群降维（比较 JY-NNK / JY-EVNK；P1 T/NK，P2 B，P3 髓系）
#
# 独立于 E:/R/fuction of cell 的 Flow_dimred_pipeline.R，也独立于 ICI_*。
# 圈门、联合降维、技术重复平均后再去掉 1 个极端生物学重复（n=2）与原免疫亚群方案一致。
# 染色仍是原 P1/P2/P3（JY_flow_panel_map.json），不是 His 靶细胞表。
#
# 组别：JY-EVNK（EVNK-1/2/3 × 技术重复 EVNK1-1/EVNK1-2）
#       JY-NNK（NNK-1/2/3 × 技术重复 NNK1-1/NNK1-2）
#
# 用法：把本方案全部 JY_* 文件放到数据目录（不要去 fuction of cell 找原流程）
#   setwd("E:/R/fuction of cell-ljy")
#   source("JY_Flow_dimred_pipeline.R")
# 总结果已经出来、只补 NKT/B 活化与耗竭：
#   source("JY_Flow_dimred_functional_state.R")
# 无 FCS 时可 Sys.setenv(FLOW_DEMO = "1")
#
# 结果：E:/R/fuction of cell-ljy/results_flow/
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANGUAGE = "en")

jy_get_script_dir <- function() {
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

jy_script_dir <- jy_get_script_dir()
jy_primary_data_dir <- "E:/R/fuction of cell-ljy"

jy_keep_cand <- function(p) {
  s <- gsub("\\", "/", as.character(p), fixed = TRUE)
  if (grepl("cell-ljy", s, ignore.case = TRUE)) return(TRUE)
  if (grepl("cell-wjz", s, ignore.case = TRUE)) return(FALSE)
  if (grepl("Internation cell immune", s, ignore.case = TRUE)) return(FALSE)
  if (grepl("fuction of cell|function of cell", s, ignore.case = TRUE)) return(FALSE)
  TRUE
}

# 只加载本目录的 JY_flow_engine.R，禁止去 fuction of cell / Internation 找原流程
load_jy_engine <- function() {
  pipe <- "JY_flow_engine.R"
  cands <- unique(c(
    file.path(jy_script_dir, pipe),
    file.path(getwd(), pipe)
  ))
  cands <- cands[vapply(cands, jy_keep_cand, logical(1))]
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit) || !nzchar(hit)) {
    stop(
      "找不到 JY_flow_engine.R。JY 方案与免疫亚群 / 靶细胞方案完全独立，",
      "不会去 E:/R/fuction of cell 找 Flow_dimred_pipeline.R，也不会用 ICI_*。\n",
      "请把这些文件一起放到 E:/R/fuction of cell-ljy：\n",
      "  JY_Flow_dimred_pipeline.R\n",
      "  JY_flow_engine.R\n",
      "  JY_flow_panel_map.json\n",
      "  JY_Flow_dimred_all_subsets.R\n",
      "  JY_Flow_dimred_trajectory.R\n",
      "  JY_Flow_dimred_functional_state.R"
    )
  }
  source(hit, local = FALSE)
  if (exists("log_msg", mode = "function")) log_msg("JY engine (this scheme only): ", hit)
  invisible(TRUE)
}

load_jy_engine()
# 与原免疫亚群方案相同：技术重复平均后去掉 1 个极端生物学重复，统计 n=2
flow_trim_bio_extremes <- TRUE

find_jy_panel_map <- function() {
  names <- c("JY_flow_panel_map.json", "JY_flow_panel_map.json.txt")
  dirs <- unique(c(
    getwd(),
    jy_script_dir,
    jy_primary_data_dir,
    "E:\\R\\fuction of cell-ljy"
  ))
  dirs <- dirs[nzchar(dirs)]
  dirs <- dirs[vapply(dirs, jy_keep_cand, logical(1))]
  for (d in dirs) {
    if (!dir.exists(d)) next
    for (nm in names) {
      p <- file.path(d, nm)
      if (file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
  }
  stop("找不到 JY_flow_panel_map.json（应与 JY_Flow_dimred_pipeline.R 同目录，或放在 E:/R/fuction of cell-ljy）")
}

resolve_jy_dir <- function() {
  env_dir <- Sys.getenv("JY_FLOW_DIR", unset = "")
  preferred <- c(
    env_dir,
    jy_primary_data_dir,
    "E:\\R\\fuction of cell-ljy",
    "E:/R/function of cell-ljy",
    getwd()
  )
  preferred <- unique(preferred[nzchar(preferred)])
  preferred <- preferred[vapply(preferred, jy_keep_cand, logical(1))]
  for (d in preferred) {
    if (dir.exists(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  jy_script_dir
}

# 覆盖本方案的数据目录与染色表（不读、不改 fuction of cell 或 Internation 的文件）
panel_map <- read_panel_map_file(find_jy_panel_map())
log_msg("JY panel map: ", find_jy_panel_map())
if (!is.null(panel_map$qc$asinh_cofactor)) asinh_cofactor <- as.numeric(panel_map$qc$asinh_cofactor)
flow_ctrl_group <- "JY-EVNK"
flow_trt_group <- "JY-NNK"
flow_group_levels <- c("JY-EVNK", "JY-NNK")
if (!is.null(panel_map$groups) && length(panel_map$groups) >= 2) {
  flow_ctrl_group <- as.character(panel_map$groups[[1]])
  flow_trt_group <- as.character(panel_map$groups[[2]])
  flow_group_levels <- c(flow_ctrl_group, flow_trt_group)
}
flow_cohort <- "JY"
pal_group <- setNames(c("#1A1A1A", "#E31A1C"), flow_group_levels)
pal_group_shape <- setNames(c(16, 15), flow_group_levels)

flow_primary_data_dir <- jy_primary_data_dir
project_dir <- resolve_jy_dir()
result_dir <- file.path(project_dir, "results_flow")
log_dir <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

find_jy_extra <- function(nm) {
  cands <- unique(c(
    file.path(jy_script_dir, nm),
    file.path(getwd(), nm)
  ))
  cands <- cands[vapply(cands, jy_keep_cand, logical(1))]
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit) || !nzchar(hit)) NA_character_ else hit
}

# -----------------------------------------------------------------------------
# 主流程（P1 / P2 / P3）
# -----------------------------------------------------------------------------
if (identical(toupper(Sys.getenv("FLOW_FUNCTIONS_ONLY", "0")), "1") ||
    identical(toupper(Sys.getenv("JY_FUNCTIONS_ONLY", "0")), "1")) {
  log_msg("JY FLOW_FUNCTIONS_ONLY=1: skip analysis")
} else {

log_msg("JY flow dir: ", project_dir)
log_msg("JY results: ", result_dir)
log_msg("JY purpose: immune-subset dimred, comparison JY-NNK / JY-EVNK (same gates as original P1/P2/P3)")
log_msg("JY stats: average tech reps, then drop 1 extreme bio-rep (max or min) per group; n=2")
file_tab <- list_unmixed_files(project_dir)
use_demo <- demo_flag
if (nrow(file_tab) == 0) {
  if (demo_flag) {
    log_msg("No unmixed FCS; FLOW_DEMO=1 -> synthetic data for plot export")
    use_demo <- TRUE
  } else {
    log_msg("No *_unmixed.fcs in ", project_dir)
    log_msg("Put JY-EVNK / JY-NNK unmixed files (EVNK1-1, NNK1-2, …) in ", jy_primary_data_dir,
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
      log_msg(pn, " missing. Put files named like JY-EVNK1-1_", pn, "_unmixed.fcs or EVNK1-1_", pn, "_unmixed.fcs")
    }
  }
  if (any(file_tab$group == "JY-NNK" & grepl("EVNK", file_tab$file, ignore.case = TRUE))) {
    stop("Filename parser classified an EVNK file as JY-NNK; refusing to continue")
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

sum_path <- file.path(result_dir, "JY_NNK_vs_JY_EVNK_lineage_stats_all_panels.csv")
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

extra <- find_jy_extra("JY_Flow_dimred_all_subsets.R")
if (!is.na(extra) && nzchar(extra)) {
  tryCatch({
    Sys.setenv(FLOW_ALL_SUBSETS_FROM_PIPELINE = "1")
    sys.source(extra, envir = .GlobalEnv)
    export_all_subsets_analysis(result_dir)
  }, error = function(e) log_msg("all-subsets summary failed: ", e$message))
  Sys.unsetenv("FLOW_ALL_SUBSETS_FROM_PIPELINE")
}

traj <- find_jy_extra("JY_Flow_dimred_trajectory.R")
if (!is.na(traj) && nzchar(traj)) {
  tryCatch({
    Sys.setenv(FLOW_TRAJECTORY_FROM_PIPELINE = "1")
    sys.source(traj, envir = .GlobalEnv)
    export_all_panel_trajectories(result_dir)
  }, error = function(e) log_msg("trajectory summary failed: ", e$message))
  Sys.unsetenv("FLOW_TRAJECTORY_FROM_PIPELINE")
}

log_msg("JY done. Keep JY_* files in ", jy_primary_data_dir,
        "; do not source Flow_dimred_pipeline.R or ICI_*.")
invisible(summaries)

}
