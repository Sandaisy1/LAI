# 从 OLD盘 目录启动分析。主脚本在上一级：TG_RNAseq_OLD_excel.R
this_dir <- tryCatch(normalizePath("."), error = function(e) getwd())
Sys.setenv(TG_RNASEQ_OLD_DIR = this_dir)
candidates <- c(
  file.path(dirname(this_dir), "TG_RNAseq_OLD_excel.R"),
  file.path(this_dir, "..", "TG_RNAseq_OLD_excel.R"),
  "TG_RNAseq_OLD_excel.R",
  file.path("E:/R/TG_BRCA/TG", "TG_RNAseq_OLD_excel.R")
)
src <- candidates[file.exists(candidates)][1]
if (is.na(src) || !nzchar(src)) {
  stop("找不到 TG_RNAseq_OLD_excel.R。请先把该脚本放到 E:/R/TG_BRCA/TG/ 再运行。")
}
message("source: ", normalizePath(src, winslash = "/", mustWork = FALSE))
source(src, encoding = "UTF-8")
