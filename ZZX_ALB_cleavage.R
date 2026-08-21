# 用 R 启动 ALB 剪切分析（分析本身是 Python，不要把 .py 粘贴进 R Console）
#
# 在 RStudio / Cursor 的 R Console 里运行：
#   source("ZZX_ALB_cleavage.R")
#
# 不要：
#   source("ZZX_ALB_cleavage.py")
#   也不要把 from __future__ import annotations 贴进 R

data_dir <- "C:/Users/Lenovo/Desktop/ZZX"
script <- file.path(getwd(), "ZZX_ALB_cleavage.py")
if (!file.exists(script)) {
  script <- "ZZX_ALB_cleavage.py"
}
if (!file.exists(script)) {
  stop("找不到 ZZX_ALB_cleavage.py，请先 setwd() 到脚本所在目录")
}

candidates <- Sys.which(c("python", "python3", "py"))
py <- unname(candidates[nzchar(candidates)][1])
if (is.na(py) || !nzchar(py)) {
  stop(
    "没有找到 Python。请打开 Terminal（不是 R Console）运行：\n",
    "  python ZZX_ALB_cleavage.py --data-dir ", data_dir
  )
}

message("Running: ", py, " ", script, " --data-dir ", data_dir)
status <- system2(py, c(script, "--data-dir", data_dir))
if (!is.na(status) && status != 0) {
  stop("Python 脚本退出码 ", status)
}
