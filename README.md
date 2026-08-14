# LAI

TCGA-BRCA 中 **每个** 神经元相关 GO 通路单独做临床关联，并寻找与该通路得分负相关的基因。

## 数据

把下列文件放到 `data/`（扩展名 `.tsv` / `.txt` / `.csv` / `.gz` 均可）：

- `TCGA-BRCA.clinical`
- `TCGA-BRCA.protein`（可选）
- `TCGA-BRCA.star_fpkm`
- `TCGA-BRCA.survival`
- `gencode.v36.annotation.gtf.gene`

## 运行

不需要 `AnnotationDbi` / `GSVA`。GO 基因集来自 QuickGO；通路得分为该通路基因的均值 z-score。可选 CRAN 包：`data.table`、`survival`。

```bash
# 可选：安装 CRAN 包
Rscript scripts/analyze_go_pathways_tcga_brca.R --install-deps --data-dir data --out-dir results

# 正式分析（17 个 GO 各自输出，不会合并）
Rscript scripts/analyze_go_pathways_tcga_brca.R --data-dir data --out-dir results

# 无数据时用合成数据检查流程
Rscript scripts/analyze_go_pathways_tcga_brca.R --demo
```

基因集缓存到 `data/cache/go_genes/`，断网时可复用。每个 GO 会单独写入 `results/clinical/`、`results/survival/`、`results/neg_genes/`。

