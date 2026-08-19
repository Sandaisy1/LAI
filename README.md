# LAI

细胞 RNA-seq（TG 敲低）：组别 **NTC**、**TG_sh1**、**TG_sh5**。Cursor 规则写在 `.cursor/rules/`，分析脚本在 `R/`。

## 分析约定（规则已固化）

必须做完三组比较，且每组都跑两套**只上调**筛选：

1. `TG_sh1 vs NTC` 与 `TG_sh5 vs NTC` 分别比较
2. `mean(TG_sh1, TG_sh5) vs NTC`
3. 两组成对结果的**共同上调基因（交集）**

筛选：

- 线性 FoldChange ≥ 1、1.25、1.5、2（FC≥1 即全部上调；“2 倍”不是 `|log2FC|≥2`）
- 按 log2FC 降序取上调前 50、75、100、150、200、250、300

每一套输出：差异表、火山图（共同上调为 sh1/sh5 散点图）、热图（列含三组）、GO、Reactome 通路、KEGG；GSEA 在每个比较上用全基因 log2FC 排序做一次。

Excel 里如果只有 NTC vs TG_sh1，脚本会再读同目录的 `genes.fpkm_tracking` / `gene_exp.diff` 补齐 TG_sh5。

## 输入

把 Cuffdiff 结果放到 `E:\R\TG_BRCA\TG`（或仓库 `data/`）：

- `shTG(20220723).xlsx`（或 `.xls`）：Cuffdiff 风格表达/差异表
- `gene_exp.diff`（推荐）
- `genes.fpkm_tracking`（三组 FPKM，缺 sh5 时必需）

## 运行

```bash
Rscript scripts/install_r_packages.R
Rscript R/analyze_tg_rnaseq.R "E:/R/TG_BRCA/TG" "E:/R/TG_BRCA/TG/results"
```

未传路径时依次尝试 `TG_RNASEQ_DIR` 环境变量、`E:/R/TG_BRCA/TG`、`./data`。

依赖：`ggplot2`、`ggrepel`、`pheatmap`、`readxl`、`clusterProfiler`、`org.Hs.eg.db`、`enrichplot`、`ReactomePA`。

## 输出

`results/<比较>/<筛选>/`

- 比较：`sh1_vs_NTC`、`sh5_vs_NTC`、`meanSH_vs_NTC`、`common_up_sh1_sh5`
- 筛选：`FC_ge_1.00` … `FC_ge_2.00`，`top_50` … `top_300`
- 每个比较目录下还有 `all_genes_ranked.csv` 与 `GSEA_full_ranked_list/`
