#!/usr/bin/env Rscript
# =============================================================================
# T-T6：只重出各亚群的活化、抑制或耗竭（T6 vs T）
#
# 总降维已经跑完时单独运行。读本方案 results_flow/P1|P2|P3/*_cell_embeddings.csv，
# 不读 FCS，不重跑 UMAP。不要 source Flow_* / ICI_* / JY_* / JZ_*。
#
#   setwd("E:/R/flow J")
#   source("T_T6_Flow_dimred_functional_state.R")
# =============================================================================

tt6_keep_cand <- function(p) {
  s <- gsub("\\", "/", as.character(p), fixed = TRUE)
  if (grepl("flow J|flowJ|flow_J|/flow-j", s, ignore.case = TRUE)) return(TRUE)
  if (grepl("cell-ljy|cell-wjz", s, ignore.case = TRUE)) return(FALSE)
  if (grepl("Internation cell immune", s, ignore.case = TRUE)) return(FALSE)
  if (grepl("fuction of cell|function of cell", s, ignore.case = TRUE)) return(FALSE)
  TRUE
}

load_tt6_engine_functions <- function() {
  # 总结果已出、单独 source 本文件时：始终从本目录重新加载 T_T6_flow_engine.R，
  # 覆盖同一会话里刚跑过的旧主流程（tt6_engine_loaded 已是 TRUE 也不跳过）。
  pipe <- "T_T6_flow_engine.R"
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
  cands <- cands[vapply(cands, tt6_keep_cand, logical(1))]
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit) || !nzchar(hit)) {
    stop("找不到 T_T6_flow_engine.R。T-T6 方案不使用 Flow_* / ICI_* / JY_* / JZ_*。")
  }
  source(hit, local = FALSE)
  if (!exists("export_functional_state_from_results", mode = "function") ||
      !exists("func_p3_state_parents", mode = "function")) {
    stop(
      "当前 T_T6_flow_engine.R 太旧，不能出全部亚群的功能状态图。\n",
      "请把最新 T_T6_* 整套覆盖到 E:/R/flow J（尤其是 T_T6_flow_engine.R），",
      "然后：\n",
      "  setwd(\"E:/R/flow J\")\n",
      "  source(\"T_T6_Flow_dimred_functional_state.R\")"
    )
  }
  invisible(TRUE)
}

load_tt6_engine_functions()

if (!identical(toupper(Sys.getenv("FLOW_FUNCTIONAL_STATE_FROM_PIPELINE", "0")), "1") &&
    !identical(toupper(Sys.getenv("FLOW_FUNCTIONS_ONLY", "0")), "1") &&
    !identical(toupper(Sys.getenv("T_T6_FUNCTIONS_ONLY", "0")), "1")) {
  export_functional_state_from_results(
    if (exists("result_dir")) result_dir else file.path(getwd(), "results_flow")
  )
}
