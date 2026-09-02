#!/usr/bin/env Rscript
# =============================================================================
# 免疫细胞亚群：只重出各亚群的活化、抑制或耗竭（图1 布局）
#
# 总降维已经跑完时单独运行。不读 FCS，不重跑 UMAP。读取已经写好的
#   results_flow/P1/P1_cell_embeddings.csv
#   results_flow/P2/P2_cell_embeddings.csv
#   results_flow/P3/P3_cell_embeddings.csv
#
# 用法：
#   setwd("E:/R/fuction of cell")
#   source("Flow_dimred_functional_state.R")
#
# JY / JZ 请用各自目录里的 JY_Flow_dimred_functional_state.R / JZ_Flow_dimred_functional_state.R
# 结果：results_flow/P1|P2|P3/functional_state/
# =============================================================================

flow_keep_functional_state_cand <- function(p) {
  s <- gsub("\\", "/", as.character(p), fixed = TRUE)
  if (grepl("cell-ljy", s, ignore.case = TRUE)) return(FALSE)
  if (grepl("cell-wjz", s, ignore.case = TRUE)) return(FALSE)
  if (grepl("Internation cell immune", s, ignore.case = TRUE)) return(FALSE)
  TRUE
}

load_flow_pipeline_functions <- function() {
  # 总结果已出、单独 source 本文件时：始终从本目录重新加载，覆盖会话里旧的 NKT/B-only 函数。
  pipe <- "Flow_dimred_pipeline.R"
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
  cands <- cands[vapply(cands, flow_keep_functional_state_cand, logical(1))]
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit) || !nzchar(hit)) stop("找不到 Flow_dimred_pipeline.R")
  old <- Sys.getenv("FLOW_FUNCTIONS_ONLY", unset = NA)
  Sys.setenv(FLOW_FUNCTIONS_ONLY = "1")
  source(hit, local = FALSE)
  if (is.na(old)) Sys.unsetenv("FLOW_FUNCTIONS_ONLY") else Sys.setenv(FLOW_FUNCTIONS_ONLY = old)
  if (!exists("export_functional_state_from_results", mode = "function") ||
      !exists("func_p3_state_parents", mode = "function")) {
    stop(
      "当前 Flow_dimred_pipeline.R 太旧，不能出全部亚群的功能状态图。\n",
      "请覆盖最新 Flow_dimred_pipeline.R 与 Flow_dimred_functional_state.R 后新开 R 会话，然后：\n",
      "  setwd(\"E:/R/fuction of cell\")\n",
      "  source(\"Flow_dimred_functional_state.R\")"
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
