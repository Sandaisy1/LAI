# 分析代码就在本文件夹：TG_RNAseq_OLD_excel.R
this_dir <- tryCatch(normalizePath("."), error = function(e) getwd())
Sys.setenv(TG_RNASEQ_OLD_DIR = this_dir)
src <- file.path(this_dir, "TG_RNAseq_OLD_excel.R")
if (!file.exists(src)) {
  src <- file.path(dirname(this_dir), "OLD盘", "TG_RNAseq_OLD_excel.R")
}
if (!file.exists(src)) {
  stop("找不到 TG_RNAseq_OLD_excel.R。它应该和本文件在同一个 OLD盘 文件夹里。")
}
message("分析代码位置: ", normalizePath(src, winslash = "/", mustWork = FALSE))
source(src, encoding = "UTF-8")
