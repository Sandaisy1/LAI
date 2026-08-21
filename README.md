# TG BRCA 细胞 RNA-seq 分析（肿瘤转移）

针对 `NTC_rep0`、`NTC_rep1`、`TG_sh1`、`TG_sh5` 四个样品。两个 NTC **不在 1-vs-1 比较里合并**。分析背景限定肿瘤转移。

## 数据位置

默认读取 `E:/R/TG_BRCA/TG`。输入为 Cuffdiff 文件（不再使用 Excel）：

- `genes.read_group_tracking`（首选）
- `genes.count_tracking`
- `genes.fpkm_tracking`

自定义转移基因列表：`metastasis_custom_genes.txt`（可改，不影响全库 p 值）。

## 运行

```r
setwd("E:/R/TG_BRCA/TG")
source("TG_RNAseq_pipeline.R")                 # 比较 1–4
source("TG_RNAseq_TGsh_mean_vs_NTC_reps.R")    # 比较 5–6，不改上面的脚本
```

## 六组比较（各自单独出图）

1. 四个 1-vs-1：`TG_sh1 vs NTC_rep0`、`TG_sh5 vs NTC_rep0`、`TG_sh1 vs NTC_rep1`、`TG_sh5 vs NTC_rep1`（无 p 值时写在 `noPvalue/`）
2. `(TG_sh1 + TG_sh5)/2` vs NTC 组均值（两侧均有 2 个样品时用 limma 估计 p）
3. 共同上调：相对 `NTC_rep0` 的 sh1 与 sh5 上调基因交集
4. 共同上调：相对 `NTC_rep1` 的 sh1 与 sh5 上调基因交集
5. `mean(TG_sh1, TG_sh5)` vs `NTC_rep0`（不要混入 NTC_rep1）
6. `mean(TG_sh1, TG_sh5)` vs `NTC_rep1`（不要混入 NTC_rep0）

## 五套分层（每组都跑）

有 p 值时用 **p < 0.05** 和 **p < 0.01**（`pvalue`，不是 padj）。无法估计 p 时不要伪造，FC/topN 只出一套，目录 `noPvalue/`。

1. `p < 0.05` 且上调 FC ≥ 1 / 1.25 / 1.5 / 2
2. `p < 0.05` 后上调 top 50 / 75 / 100 / 150 / 200 / 250 / 300
3. `p < 0.01` 且上调 FC ≥ 1 / 1.25 / 1.5 / 2
4. `p < 0.01` 后同样 7 个 topN
5. 全部上调 **和** 全部下调：`AllDE/all_up/`、`AllDE/all_down/`、`AllDE/all_up_and_down/`

每个非空子集输出差异表、火山图、热图、GO、通路、KEGG、GSVA、GSEA。

## 结果目录

```
results/
  TG_sh1_vs_NTC_rep0/
    noPvalue/FoldChange/FC_1.5/
    noPvalue/TopRank/top100/
    AllDE/all_up/
    AllDE/all_down/
    AllDE/all_up_and_down/
    00_GSEA_all_genes_NOT_FC_or_topN/
    00_GSVA_all_genes_NOT_FC_or_topN/
    Focused_metastasis/
  TGsh_mean_vs_NTC/
    p0.05/  p0.01/   # 2-vs-2 有 limma p 时
    AllDE/
  TGsh_mean_vs_NTC_rep0/
  TGsh_mean_vs_NTC_rep1/
```

同一档里有三套富集，不要只看 GSEA 文件夹：

- `GO/`、`Pathway/`、`KEGG/`：过表达分析（ORA），文件名以 `ORA_` 开头
- `GSEA/`：GSEA，文件名以 `GSEA_` 开头
- `GSVA/`：GSVA，文件名以 `GSVA_` 开头

分层图示例：

```
results/TGsh_mean_vs_NTC/p0.01/FoldChange/FC_1.5/p0.01_FC_1.5_volcano.pdf
results/TG_sh1_vs_NTC_rep0/noPvalue/TopRank/top100/noPvalue_top100_heatmap.pdf
```

重新运行前建议先删掉旧的 `results/`。

## 热图上的基因名（XLOC / 逗号）

这些来自 Cuffdiff 的 `gene_short_name`，**不是 FC 算错**：

- `XLOC_003812`：Cufflinks 组装出来的位点 ID，没有官方基因符号时会保留
- `SAA2,SAA2-SAA4,SAA4`：重叠基因座被写成一条复合名；脚本会拆成官方符号（这里取 `SAA2`）

读入时会清洗复合名并优先用官方符号。没有符号的 novel locus 仍会显示为 `XLOC_`，这些行会留在差异分析里，但 GO/KEGG 通常映射不上。

## 肿瘤转移通路怎么看

**不能、也不该**去改全库 GO/KEGG 的 p 值，把转移相关通路人为抬到第一。全库排名由统计量和多重检验决定。

专项检验：

1. 每个比较目录下的 `Focused_metastasis/`：只用 EMT、迁移、侵袭、ECM、血管生成和自定义转移基因集做 GSEA。
2. 全库 GO/KEGG/GSEA 表旁边的 `*_FOCUS_metastasis.csv`：把匹配到的条目抽出来，**保留原始 p 值和 `genome_wide_rank`**。
3. 每个 FC / topN / AllDE 子文件夹里也有 `Focused_metastasis/`，是针对该基因子集的专项 ORA。

看专项结果请先看比较目录下的 `Focused_metastasis/`（全基因 GSEA）。若专项分析仍不显著，说明这些通路在本数据里没有协同变化。
