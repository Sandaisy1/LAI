#!/usr/bin/env Rscript
# =============================================================================
# JY：只重出各亚群的活化、抑制或耗竭（JY-NNK / JY-EVNK）
# 读本方案 results_flow/P1|P2|P3/*_cell_embeddings.csv，不读 FCS，不重跑 UMAP。
# 不要 source Flow_dimred_pipeline.R 或 ICI_* / JZ_*。
#
#   setwd("E:/R/fuction of cell-ljy")
#   source("JY_Flow_dimred_functional_state.R")
#
# 若刚才刚跑过主流程，请先关掉 R 再开，或至少覆盖最新 JY_flow_engine.R。
# =============================================================================

jy_keep_cand <- function(p) {
  s <- gsub("\\", "/", as.character(p), fixed = TRUE)
  if (grepl("cell-ljy", s, ignore.case = TRUE)) return(TRUE)
  if (grepl("cell-wjz", s, ignore.case = TRUE)) return(FALSE)
  if (grepl("Internation cell immune", s, ignore.case = TRUE)) return(FALSE)
  if (grepl("fuction of cell|function of cell", s, ignore.case = TRUE)) return(FALSE)
  TRUE
}

load_jy_engine_functions <- function() {
  # 同一会话里刚跑过旧主流程时 jy_engine_loaded 已是 TRUE，但函数库是旧的。
  # 缺 export_functional_state_from_results 时必须再 source 本目录的 JY_flow_engine.R。
  if (exists("export_functional_state_from_results", mode = "function") &&
      isTRUE(get0("jy_engine_loaded", ifnotfound = FALSE))) {
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
  if (!exists("export_functional_state_from_results", mode = "function")) {
    stop(
      "当前 JY_flow_engine.R 太旧，没有 export_functional_state_from_results。\n",
      "请把最新 JY_* 整套覆盖到 E:/R/fuction of cell-ljy（尤其是 JY_flow_engine.R），",
      "关掉 R 再开，然后：\n",
      "  setwd(\"E:/R/fuction of cell-ljy\")\n",
      "  source(\"JY_Flow_dimred_functional_state.R\")"
    )
  }
  invisible(TRUE)
}

load_jy_engine_functions()

if (!identical(toupper(Sys.getenv("FLOW_FUNCTIONAL_STATE_FROM_PIPELINE", "0")), "1") &&
    !identical(toupper(Sys.getenv("FLOW_FUNCTIONS_ONLY", "0")), "1") &&
    !identical(toupper(Sys.getenv("JY_FUNCTIONS_ONLY", "0")), "1")) {
  export_functional_state_from_results(
    if (exists("result_dir")) result_dir else file.path(getwd(), "results_flow")
  )
}
