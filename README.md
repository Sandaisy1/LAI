# LAI — TCGA BRCA Nature 2012 GO 分析

cBioPortal 数据目录（本机）：`E:/R/TCGA_BRCA`

## 运行

把 `scripts/analyze_go_subtype_metastasis.R` 拷到数据目录，或直接指定路径：

```r
# RStudio：打开该脚本后 Source
# 或：
source("scripts/analyze_go_subtype_metastasis.R")
```

```bash
Rscript scripts/analyze_go_subtype_metastasis.R --data-dir "E:/R/TCGA_BRCA" --out-dir "E:/R/TCGA_BRCA/results"
Rscript scripts/analyze_go_subtype_metastasis.R --demo
```

结果在 `E:/R/TCGA_BRCA/results/`：`bubble_subtype.pdf`、`bubble_metastasis.pdf`、`neg_genes/`。
