#!/usr/bin/env Rscript
# =============================================================================
# Internation cell immune：只重出 His-FITC+ 各亚群的功能状态（图1 布局）
#
# 总降维已经跑完时单独运行。读本方案 results_flow/P1|P2|P3/*_cell_embeddings.csv，
# 不读 FCS，不重跑 UMAP。必须在 His+ 母群上做。不要 source Flow_* / JY_* / JZ_*。
#
#   setwd("E:/R/Internation cell immune")
#   source("ICI_Flow_dimred_functional_state.R")
#
# 结果：results_flow/P1|P2|P3/functional_state/
# =============================================================================

ici_keep_functional_state_cand <- function(p) {
  s <- gsub("\\", "/", as.character(p), fixed = TRUE)
  if (grepl("cell-ljy", s, ignore.case = TRUE)) return(FALSE)
  if (grepl("cell-wjz", s, ignore.case = TRUE)) return(FALSE)
  if (grepl("fuction of cell|function of cell", s, ignore.case = TRUE) &&
      !grepl("Internation cell immune", s, ignore.case = TRUE)) {
    return(FALSE)
  }
  TRUE
}

load_ici_engine_for_functional_state <- function() {
  pipe <- "ICI_flow_engine.R"
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
  cands <- cands[vapply(cands, ici_keep_functional_state_cand, logical(1))]
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit) || !nzchar(hit)) {
    stop("找不到 ICI_flow_engine.R。靶细胞方案不使用 Flow_dimred_pipeline.R。")
  }
  old_flow <- Sys.getenv("FLOW_FUNCTIONS_ONLY", unset = NA)
  old_ici <- Sys.getenv("ICI_FUNCTIONS_ONLY", unset = NA)
  Sys.setenv(FLOW_FUNCTIONS_ONLY = "1", ICI_FUNCTIONS_ONLY = "1")
  source(hit, local = FALSE)
  if (is.na(old_flow)) Sys.unsetenv("FLOW_FUNCTIONS_ONLY") else Sys.setenv(FLOW_FUNCTIONS_ONLY = old_flow)
  if (is.na(old_ici)) Sys.unsetenv("ICI_FUNCTIONS_ONLY") else Sys.setenv(ICI_FUNCTIONS_ONLY = old_ici)
  if (!exists("export_functional_state_from_results", mode = "function") ||
      !exists("func_p3_state_parents", mode = "function") ||
      !exists("ici_keep_his_parent_cells", mode = "function")) {
    stop(
      "当前 ICI_flow_engine.R 太旧，不能出 His+ 亚群功能状态图。\n",
      "请把最新 ICI_* 整套覆盖到 E:/R/Internation cell immune（尤其是 ICI_flow_engine.R），",
      "然后：\n",
      "  setwd(\"E:/R/Internation cell immune\")\n",
      "  source(\"ICI_Flow_dimred_functional_state.R\")"
    )
  }
  invisible(TRUE)
}

load_ici_engine_for_functional_state()
flow_trim_bio_extremes <- FALSE

if (!identical(toupper(Sys.getenv("FLOW_FUNCTIONAL_STATE_FROM_PIPELINE", "0")), "1") &&
    !identical(toupper(Sys.getenv("FLOW_FUNCTIONS_ONLY", "0")), "1") &&
    !identical(toupper(Sys.getenv("ICI_FUNCTIONS_ONLY", "0")), "1")) {
  export_functional_state_from_results(
    if (exists("result_dir")) result_dir else file.path(getwd(), "results_flow")
  )
}
