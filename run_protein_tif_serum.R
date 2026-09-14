#!/usr/bin/env Rscript
# 一键运行 TIF / 血清蛋白质组 T vs N
#
# 在 R / RStudio 控制台（不要输入 Rscript）：
#   setwd("E:/R/Protein TIF serum")
#   source("run_protein_tif_serum.R")
# 或：
#   source("Protein_TIF_Serum_pipeline.R")
#
# 在 Windows 命令提示符 / PowerShell（不是 R 控制台）：
#   Rscript run_protein_tif_serum.R "E:/R/Protein TIF serum"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) {
  if (grepl("^Rscript$", args[[1]], ignore.case = TRUE)) {
    stop("Rscript 是系统命令，不能在 R 控制台里运行。请改用：source(\"run_protein_tif_serum.R\")", call. = FALSE)
  }
  Sys.setenv(PROTEIN_TIF_SERUM_DIR = args[[1]])
}

this_file <- NULL
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) > 0) {
  this_file <- sub("^--file=", "", file_arg[1])
} else {
  for (i in sys.nframe():1) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile) && nzchar(ofile)) {
      this_file <- ofile
      break
    }
  }
}

if (!is.null(this_file) && file.exists(this_file)) {
  script_dir <- dirname(normalizePath(this_file, winslash = "/", mustWork = FALSE))
} else {
  script_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

if (!nzchar(Sys.getenv("PROTEIN_TIF_SERUM_DIR", unset = ""))) {
  default_data <- "E:/R/Protein TIF serum"
  if (dir.exists(default_data) &&
      file.exists(file.path(default_data, "TIF_report.pg_matrix")) &&
      file.exists(file.path(default_data, "Serum_report.pg_matrix"))) {
    Sys.setenv(PROTEIN_TIF_SERUM_DIR = default_data)
  }
}

pipeline <- file.path(script_dir, "Protein_TIF_Serum_pipeline.R")
if (!file.exists(pipeline)) {
  pipeline <- file.path(getwd(), "Protein_TIF_Serum_pipeline.R")
}
if (!file.exists(pipeline)) {
  stop(
    "找不到 Protein_TIF_Serum_pipeline.R。\n",
    "请把 Protein_TIF_Serum_pipeline.R 和 run_protein_tif_serum.R 拷到数据目录后，在 R 里运行：\n",
    "  setwd(\"E:/R/Protein TIF serum\")\n",
    "  source(\"Protein_TIF_Serum_pipeline.R\")",
    call. = FALSE
  )
}
source(pipeline, encoding = "UTF-8")
