# =============================================================================
# 气泡图外观：只改这里的数字，然后运行
#   setwd("E:/R/TG_BRCA/TG")
#   source("TG_RNAseq_bubble_restyle.R")
# 不需要重跑差异分析或 enrichGO。
#
# 每张图旁边还有 *_plotdata.csv：可删行、改顺序后再 restyle。
# =============================================================================

# --- 气泡大小（ggplot size，不是像素）。越大点越大 ---
bubble_size_min <- 6
bubble_size_max <- 18

# --- 坐标轴字体 ---
axis_text_y_size <- 10   # 左侧通路名
axis_text_x_size <- 11   # 底部 GeneRatio 刻度
axis_title_size  <- 12   # “GeneRatio” 轴标题

# --- 图标题、图例 ---
title_size        <- 12
legend_text_size  <- 10
legend_title_size <- 11
base_size         <- 12  # 其余文字的基准字号

# --- 图幅（英寸）---
plot_width  <- 9
plot_height <- NA        # 填数字则固定高度，例如 7；NA = 按条目数自动
point_stroke <- 0.5      # 气泡描边粗细

# --- 只重画一张图时填 CSV 完整路径；留空则重画 results/ 下全部气泡图 ---
only_this_csv <- ""
