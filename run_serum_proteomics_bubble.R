# =============================================================================
# 血清蛋白质组排名气泡图：启动器
# 在 R 控制台整段粘贴，或先把本文件放到数据目录再 source("run_serum_proteomics_bubble.R")
#
# 不要只写 source("serum_proteomics_ranked_bubble.R")——除非当前工作目录里真有这个文件。
# =============================================================================

data_dir <- "E:/天府/实验管理/课题/赵章寻/血清蛋白质组学"
script_name <- "serum_proteomics_ranked_bubble.R"

this_dir <- tryCatch({
  ofile <- sys.frame(1)$ofile
  if (!is.null(ofile) && nzchar(ofile)) dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)) else NA_character_
}, error = function(e) NA_character_)

candidates <- unique(c(
  file.path(data_dir, script_name),
  if (!is.na(this_dir)) file.path(this_dir, script_name) else NA_character_,
  file.path(getwd(), script_name)
))
candidates <- candidates[nzchar(candidates) & !is.na(candidates)]
hit <- candidates[file.exists(candidates)][1]

if (length(hit) == 0 || is.na(hit)) {
  if (!dir.exists(data_dir)) {
    stop(
      "找不到数据目录: ", data_dir, "\n",
      "当前工作目录是: ", getwd(), "\n",
      "请确认 U 盘/网盘路径，或 Sys.setenv(SERUM_PROTEOMICS_DIR = \"你的目录\")"
    )
  }
  dest <- file.path(data_dir, script_name)
  urls <- c(
    "https://raw.githubusercontent.com/Sandaisy1/LAI/cursor/serum-proteomics-bubble-e55e/serum_proteomics_ranked_bubble.R",
    "https://raw.githubusercontent.com/Sandaisy1/LAI/main/serum_proteomics_ranked_bubble.R"
  )
  for (u in urls) {
    message("本地没有脚本，尝试下载: ", u)
    ok <- tryCatch({
      download.file(u, dest, mode = "wb", quiet = FALSE)
      TRUE
    }, error = function(e) {
      message("下载失败: ", e$message)
      FALSE
    })
    if (isTRUE(ok) && file.exists(dest) && isTRUE(file.info(dest)$size > 200)) {
      hit <- dest
      break
    }
  }
}

if (length(hit) == 0 || is.na(hit) || !file.exists(hit)) {
  stop(
    "无法打开文件 '", script_name, "': No such file or directory\n",
    "当前工作目录: ", getwd(), "\n",
    "请把 ", script_name, " 复制到数据目录后执行：\n",
    "  setwd(\"", data_dir, "\")\n",
    "  source(\"", script_name, "\", encoding = \"UTF-8\")\n",
    "或在仓库根目录（有该 .R 文件的地方）运行本启动器。"
  )
}

if (dir.exists(data_dir)) {
  setwd(data_dir)
}
message("工作目录: ", getwd())
message("脚本: ", hit)
source(hit, encoding = "UTF-8")
