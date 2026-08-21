#!/usr/bin/env Rscript
# 分析代码在 OLD盘 文件夹里：
#   OLD盘/TG_RNAseq_OLD_excel.R
# 本文件只负责找到并运行它。
#
# 用法：
#   setwd("E:/R/TG_BRCA/TG")
#   source("TG_RNAseq_OLD_excel.R")
# 或打开：
#   E:/R/TG_BRCA/TG/OLD盘/TG_RNAseq_OLD_excel.R
#   然后 source("TG_RNAseq_OLD_excel.R")

candidates <- c(
  file.path("OLD盘", "TG_RNAseq_OLD_excel.R"),
  file.path(getwd(), "OLD盘", "TG_RNAseq_OLD_excel.R"),
  "E:/R/TG_BRCA/TG/OLD盘/TG_RNAseq_OLD_excel.R",
  file.path(getwd(), "TG_RNAseq_OLD_excel.R")
)
src <- candidates[file.exists(candidates)][1]
if (is.na(src) || !nzchar(src)) {
  stop(
    "找不到分析脚本。请打开这个文件：\n",
    "  OLD盘/TG_RNAseq_OLD_excel.R\n",
    "或 GitHub：LAI 仓库里的 OLD盘/TG_RNAseq_OLD_excel.R"
  )
}
message("分析代码位置: ", normalizePath(src, winslash = "/", mustWork = FALSE))
source(src, encoding = "UTF-8")
