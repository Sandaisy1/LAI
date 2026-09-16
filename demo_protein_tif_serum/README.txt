演示用 DIA-NN pg_matrix（不是真实实验数据）。
真实数据请放到 E:/R/Protein TIF serum 后，在 R 控制台运行：

  setwd("E:/R/Protein TIF serum")
  source("Protein_TIF_Serum_TVsN_standalone.R")   # 独立脚本：火山图 + 上调 GO/KEGG + 排名图
  source("Protein_TIF_Serum_pipeline.R")          # 另一份脚本：同样 1–2 条，第 3 条画热图

不要在 R 控制台输入 Rscript。Windows 命令提示符才用：
  Rscript run_protein_tif_serum.R "E:/R/Protein TIF serum"
