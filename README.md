# TG BRCA 细胞 RNA-seq 分析

针对 `NTC_rep0`、`NTC_rep1`、`TG_sh1`、`TG_sh5` 四个样品的 RNA-seq 分析，**下调方向**。两个 NTC **不在 1-vs-1 比较里合并**（仅设计 2 用 NTC 组均值）。

## 数据位置

默认读取 `E:/R/TG_BRCA/TG`。输入为 Cuffdiff 文件：

- `genes.read_group_tracking`（首选）
- `genes.count_tracking`
- `genes.fpkm_tracking`

## 运行

```r
setwd("E:/R/TG_BRCA/TG")
source("TG_RNAseq_pipeline.R")
```

六组比较会先一起过滤低表达、再标准化，然后各自出图。

## 六种比较（各自单独出图）

1. 四个 1-vs-1：`TG_sh1 vs NTC_rep0`、`TG_sh5 vs NTC_rep0`、`TG_sh1 vs NTC_rep1`、`TG_sh5 vs NTC_rep1`（标准化后直接算 FC，无 P 值）
2. `(TG_sh1 + TG_sh5)/2` vs NTC 组均值（`NTC_rep0` 与 `NTC_rep1`）
3. 共同下调：相对 `NTC_rep0` 的 sh1 与 sh5 下调基因交集
4. 共同下调：相对 `NTC_rep1` 的 sh1 与 sh5 下调基因交集
5. `(TG_sh1 + TG_sh5)/2` vs `NTC_rep0`
6. `(TG_sh1 + TG_sh5)/2` vs `NTC_rep1`

每个比较**先**按两种子集策略出差异表、火山图、热图、单独 GO/通路/KEGG；**全部子集图完成后再**跑 GSEA。全基因 GSEA 放在最后，报错不会跳过前面的图：

1. 下调 FoldChange：treat/control ≤ 1 / 0.8 / 0.6 / 0.5
2. 下调排名：top 50 / 75 / 100 / 150 / 200 / 250 / 300（log2FC 最负）

## 结果目录

```
results/
  TG_sh1_vs_NTC_rep0/
  TG_sh5_vs_NTC_rep0/
  TG_sh1_vs_NTC_rep1/
  TG_sh5_vs_NTC_rep1/
  TGsh_mean_vs_NTC/
  common_down_vs_NTC_rep0/
  common_down_vs_NTC_rep1/
  TGsh_mean_vs_NTC_rep0/
  TGsh_mean_vs_NTC_rep1/
```

同一档里有两套富集，不要只看 GSEA 文件夹：

- `GO/`、`Pathway/`、`KEGG/`：单独的过表达分析（ORA），文件名以 `ORA_` 开头
- `GSEA/`：GSEA 分析，文件名以 `GSEA_` 开头

**注意：** 比较目录里的 `99_GSEA_all_genes_NOT_FC_or_topN` 是全部基因的 GSEA，**最后才跑**，**不是** FC/topN 分层图。分层图在：

```
results/TG_sh1_vs_NTC_rep0/FoldChange/FC_0.5/FC_0.5_volcano.pdf
results/TG_sh1_vs_NTC_rep0/TopRank/top100/top100_heatmap.pdf
results/TG_sh1_vs_NTC_rep0/FoldChange/FC_0.5/GO/FC_0.5_ORA_GO_BP_dotplot.pdf
```

文件名和图标题都会带上 `FC_0.5` 或 `top100`。重新运行前建议先删掉旧的 `results/`。

## 热图上的基因名（XLOC / 逗号）

这些来自 Cuffdiff 的 `gene_short_name`，**不是 FC 算错**：

- `XLOC_003812`：Cufflinks 组装出来的位点 ID，没有官方基因符号时会保留
- `SAA2,SAA2-SAA4,SAA4`：重叠基因座被写成一条复合名；脚本会拆成官方符号（这里取 `SAA2`）

读入时会清洗复合名并优先用官方符号。没有符号的 novel locus 仍会显示为 `XLOC_`，这些行会留在差异分析里，但 GO/KEGG 通常映射不上。

## 细胞骨架运动 / 线粒体通路怎么看

**不能、也不该**去改全库 GO/KEGG 的 p 值，把这两类通路人为抬到第一。全库排名由统计量和多重检验决定。

可以做的是专项检验（脚本已加）：

1. 每个比较目录下的 `Focused_cytoskeleton_mito/`：只用细胞骨架运动、细胞迁移、线粒体相关基因集做 GSEA。
2. 全库 GO/KEGG/GSEA 表旁边的 `*_FOCUS_cytoskeleton_mito.csv`：把匹配到的条目抽出来，**保留原始 p 值和 `genome_wide_rank`**。
3. 每个 FC / topN 子文件夹里也有 `Focused_cytoskeleton_mito/`，是针对该基因子集的专项 ORA。

## 只跑设计 5 / 6

`TG_RNAseq_TGsh_mean_vs_NTC_reps.R` 可以单独跑 `mean(TG_sh1, TG_sh5)` vs 单个 NTC。主流程已经包含这两组；若刚刚跑过主脚本，该文件会跳过以免重复。

```r
setwd("E:/R/TG_BRCA/TG")
source("TG_RNAseq_TGsh_mean_vs_NTC_reps.R")
```
