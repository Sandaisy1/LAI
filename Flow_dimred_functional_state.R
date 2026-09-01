#!/usr/bin/env Rscript
# =============================================================================
# 免疫细胞亚群：只重出 NKT / B 亚群的活化、效应、耗竭
#
# 不读 FCS，不重跑 UMAP。读取已经写好的
#   results_flow/P1/P1_cell_embeddings.csv
#   results_flow/P2/P2_cell_embeddings.csv
#
# 用法（总结果已经出来之后）：
#   setwd("E:/R/fuction of cell")
#   source("Flow_dimred_functional_state.R")
#
# JY / JZ 请用各自目录里的 JY_Flow_dimred_functional_state.R / JZ_Flow_dimred_functional_state.R
# 结果：results_flow/P1/functional_state/ 与 P2/functional_state/
# =============================================================================

load_flow_pipeline_functions <- function() {
  if (exists("export_functional_state_from_results", mode = "function")) {
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
  if (!exists("export_functional_state_from_results", mode = "function")) {
    stop(
      "当前 Flow_dimred_pipeline.R 太旧，没有 export_functional_state_from_results。\n",
      "请覆盖最新 Flow_dimred_pipeline.R 与 Flow_dimred_functional_state.R 后新开 R 会话。"
    )
  }
  invisible(TRUE)
}

load_flow_pipeline_functions()

if (!identical(toupper(Sys.getenv("FLOW_FUNCTIONAL_STATE_FROM_PIPELINE", "0")), "1") &&
    !identical(toupper(Sys.getenv("FLOW_FUNCTIONS_ONLY", "0")), "1")) {
  export_functional_state_from_results(
    if (exists("result_dir")) result_dir else file.path(getwd(), "results_flow")
  )
}
