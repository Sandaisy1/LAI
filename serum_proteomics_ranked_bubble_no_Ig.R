# =============================================================================
# 去除各类免疫球蛋白后的丰度排名气泡图（另写，不改原脚本出图）
#
# 把本文件和 serum_proteomics_ranked_bubble.R 放到同一目录后：
#   setwd("E:/天府/实验管理/课题/赵章寻/血清蛋白质组学")
#   source("serum_proteomics_ranked_bubble_no_Ig.R", encoding = "UTF-8")
# =============================================================================

SERUM_PROTEOMICS_SKIP_MAIN <- TRUE

find_base_script <- function() {
  this_dir <- tryCatch({
    ofile <- sys.frame(1)$ofile
    if (!is.null(ofile) && nzchar(ofile)) dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)) else NA_character_
  }, error = function(e) NA_character_)
  data_dir <- "E:/天府/实验管理/课题/赵章寻/血清蛋白质组学"
  candidates <- unique(c(
    if (!is.na(this_dir)) file.path(this_dir, "serum_proteomics_ranked_bubble.R") else NA_character_,
    file.path(getwd(), "serum_proteomics_ranked_bubble.R"),
    file.path(data_dir, "serum_proteomics_ranked_bubble.R")
  ))
  candidates <- candidates[!is.na(candidates) & file.exists(candidates)]
  if (length(candidates) == 0) {
    stop(
      "找不到 serum_proteomics_ranked_bubble.R。\n",
      "请与本文件放在同一目录后再 source(\"serum_proteomics_ranked_bubble_no_Ig.R\", encoding = \"UTF-8\")"
    )
  }
  candidates[[1]]
}

base_script <- find_base_script()
message("载入: ", base_script)
source(base_script, encoding = "UTF-8")

run_serum_abundance_bubble(
  result_subdir = "serum_proteomics_bubble_no_Ig",
  drop_immunoglobulin = TRUE,
  title_prefix = "去除免疫球蛋白后 血清蛋白丰度排名气泡图",
  subtitle = "去掉 IGH/IGK/IGL/JCHAIN 等后，两样品平均丰度排名；气泡大小一致"
)
