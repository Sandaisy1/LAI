#!/usr/bin/env Rscript
# =============================================================================
# JY：只重出 NKT / B 亚群的活化、效应、耗竭（JY-NNK / JY-EVNK）
# 读本方案 results_flow/P1|P2/*_cell_embeddings.csv，不读 FCS，不重跑 UMAP。
# 不要 source Flow_dimred_pipeline.R 或 ICI_* / JZ_*。
#
#   setwd("E:/R/fuction of cell-ljy")
#   source("JY_Flow_dimred_functional_state.R")
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

load_jy_engine_functions()

if (!identical(toupper(Sys.getenv("FLOW_FUNCTIONAL_STATE_FROM_PIPELINE", "0")), "1") &&
    !identical(toupper(Sys.getenv("FLOW_FUNCTIONS_ONLY", "0")), "1") &&
    !identical(toupper(Sys.getenv("JY_FUNCTIONS_ONLY", "0")), "1")) {
  export_functional_state_from_results(
    if (exists("result_dir")) result_dir else file.path(getwd(), "results_flow")
  )
}
