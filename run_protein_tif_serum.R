#!/usr/bin/env Rscript
# 一键运行 TIF / 血清蛋白质组 T vs N
# 用法：
#   Rscript run_protein_tif_serum.R
#   Rscript run_protein_tif_serum.R "E:/R/Protein TIF serum"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) {
  Sys.setenv(PROTEIN_TIF_SERUM_DIR = args[[1]])
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  script_dir <- getwd()
}
pipeline <- file.path(script_dir, "Protein_TIF_Serum_pipeline.R")
if (!file.exists(pipeline)) {
  stop("找不到 Protein_TIF_Serum_pipeline.R，请把本脚本和流程脚本放在同一目录")
}
source(pipeline, encoding = "UTF-8")
