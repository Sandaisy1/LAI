# OLD盘 Excel 分析

R 脚本不在这个文件夹里，而在**上一级仓库根目录**：

- `../TG_RNAseq_OLD_excel.R`

请把本分析用的 Excel（截图那份，含 `gene` / `value_1` / `value_2` / `log2(fold_change)` / `p_value`）放到这个 `OLD盘` 目录。

然后：

```r
setwd("E:/R/TG_BRCA/TG")          # 或仓库根目录
source("TG_RNAseq_OLD_excel.R")
```

结果写在：`results_TG_sh1_vs_NTC/`
