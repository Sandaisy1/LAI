# TG BRCA 细胞 RNA-seq 分析

针对 `NTC_rep0`、`NTC_rep1`、`TG_sh1`、`TG_sh5` 四个样品的 RNA-seq 分析。两个 NTC **不在 1-vs-1 比较里合并**。

## 数据位置

默认读取 `E:/R/TG_BRCA/TG`。输入为 Cuffdiff 文件（不再使用 Excel）：

- `genes.read_group_tracking`（首选）
- `genes.count_tracking`
- `genes.fpkm_tracking`

## 运行

```r
setwd("E:/R/TG_BRCA/TG")
source("TG_RNAseq_pipeline.R")
```

## 四种比较（各自单独出图）

1. 四个 1-vs-1：`TG_sh1 vs NTC_rep0`、`TG_sh5 vs NTC_rep0`、`TG_sh1 vs NTC_rep1`、`TG_sh5 vs NTC_rep1`（标准化后直接算 FC，无 P 值）
2. `(TG_sh1 + TG_sh5)/2` vs NTC 组均值（`NTC_rep0` 与 `NTC_rep1`）
3. 共同上调：相对 `NTC_rep0` 的 sh1 与 sh5 上调基因交集
4. 共同上调：相对 `NTC_rep1` 的 sh1 与 sh5 上调基因交集

每个比较再按 FC ≥ 1 / 1.25 / 1.5 / 2，以及上调 top 50–300，分别输出差异表、火山图、热图、GO、通路、KEGG、GSEA。

## 结果目录

```
results/
  TG_sh1_vs_NTC_rep0/
  TG_sh5_vs_NTC_rep0/
  TG_sh1_vs_NTC_rep1/
  TG_sh5_vs_NTC_rep1/
  TGsh_mean_vs_NTC/
  common_up_vs_NTC_rep0/
  common_up_vs_NTC_rep1/
```

同一档里有两套富集，不要只看 GSEA 文件夹：

- `GO/`、`Pathway/`、`KEGG/`：单独的过表达分析（ORA），文件名以 `ORA_` 开头
- `GSEA/`：GSEA 分析，文件名以 `GSEA_` 开头


**注意：** 比较目录里的 `00_GSEA_all_genes_NOT_FC_or_topN` 是全部基因的 GSEA，**不是** FC/topN 分层图。分层图在：

```
results/TG_sh1_vs_NTC_rep0/FoldChange/FC_1.5/FC_1.5_volcano.pdf
results/TG_sh1_vs_NTC_rep0/TopRank/top100/top100_heatmap.pdf
results/TG_sh1_vs_NTC_rep0/FoldChange/FC_1.5/GO/FC_1.5_GO_BP_dotplot.pdf
```

文件名和图标题都会带上 `FC_1.5` 或 `top100`。重新运行前建议先删掉旧的 `results/`。

## 热图上的基因名（XLOC / 逗号）

这些来自 Cuffdiff 的 `gene_short_name`，**不是 FC 算错**：

- `XLOC_003812`：Cufflinks 组装出来的位点 ID，没有官方基因符号时会保留
- `SAA2,SAA2-SAA4,SAA4`：重叠基因座被写成一条复合名；脚本会拆成官方符号（这里取 `SAA2`）

读入时会清洗复合名并优先用官方符号。没有符号的 novel locus 仍会显示为 `XLOC_`，这些行会留在差异分析里，但 GO/KEGG 通常映射不上。

## 细胞骨架运动 / 线粒体通路怎么看

**不能、也不该**去改全库 GO/KEGG 的 p 值，把这两类通路人为抬到第一。全库排名由统计量和多重检验决定；这份数据的高 FC 基因如果主要是炎症/急性期（例如 SAA、CCL2），全库 GO 就会先出现那些条目。

可以做的是专项检验（脚本已加）：

1. 每个比较目录下的 `Focused_cytoskeleton_mito/`：只用细胞骨架运动、细胞迁移、线粒体相关基因集做 GSEA，所以这些条目会排在**这个文件夹**的前面。同时有对应基因热图。
2. 全库 GO/KEGG/GSEA 表旁边的 `*_FOCUS_cytoskeleton_mito.csv`：把匹配到的条目抽出来，**保留原始 p 值和 `genome_wide_rank`**。
3. 每个 FC / topN 子文件夹里也有 `Focused_cytoskeleton_mito/`，是针对该基因子集的专项 ORA。

看专项结果请先看比较目录下的 `Focused_cytoskeleton_mito/`（全基因 GSEA），不要只看 topN 的 ORA。若专项分析仍不显著，说明这些通路在本数据里没有协同变化。

## 额外两组（新脚本，不改原流程）

`TG_RNAseq_TGsh_mean_vs_NTC_reps.R` 在原四种比较之外再做：

1. `mean(TG_sh1, TG_sh5)` vs `NTC_rep0`
2. `mean(TG_sh1, TG_sh5)` vs `NTC_rep1`

同样按 FC ≥ 1 / 1.25 / 1.5 / 2 和上调 top 50–300 出差异表、火山图、热图、GO、通路、KEGG、GSEA。结果在：

```
results/TGsh_mean_vs_NTC_rep0/
results/TGsh_mean_vs_NTC_rep1/
```

```r
setwd("E:/R/TG_BRCA/TG")
source("TG_RNAseq_pipeline.R")                 # 原四种比较，不变
source("TG_RNAseq_TGsh_mean_vs_NTC_reps.R")    # 只加上面两组
```

也可以只跑这个新脚本（会自己读入并标准化数据）。

## Wang 2024 空间转录组：神经浸润（新脚本，不改 1–4）

先 **导出能做近神经 vs 远神经的病人** 的近/远肿瘤 **基因表达**（RNA，不是蛋白）→ `results/00_expression_near_vs_far/`  
再把这些病人合在一起，看近 vs 远的高表达基因 → `results/04_COMBINED_near_vs_far_high_genes.csv`

另外三件事（**不做 GO/KEGG/GSEA**）：

1. 神经浸润相关分子和统计方案 → `results/00_Q1_PROTOCOL.txt`
2. **肿瘤细胞高表达、可能促使神经浸润的基因** → `results/02_Q2_READ_THIS_tumor_genes_may_promote_nerve.csv`
3. **肿瘤细胞低表达、可能促使神经浸润的基因** → `results/03_Q3_READ_THIS_tumor_genes_low_may_promote_nerve.csv`

```r
setwd("E:/R/Nerve")
Sys.setenv(WANG_ST_DIR = "E:/R/Nerve")
source("Wang_ST_nerve_infiltration.R")
```

不需要 `TG_RNAseq_pipeline.R`。先 **p < 0.05**。问题2上调 **FC ≥ 1.25 和 1.5**；问题3下调 **FC < 1 即可**（1.25 / 1.5 只是更严的可选分层）。原 Cuffdiff 六组比较不变。

真正必需的是 `Robjects/`（`counts` + `annotsBySpot`；先解压 `Robjects.tar`）。

`ids.RDS` 可选。Windows 解 `Clinical.tar` 时常把 `Clinical/ids.RDS` 展成根目录的 `Clinicalids.RDS`；`Clinical.RDS` / `Clinical.xlsx` 是临床表。没有 ids 也能跑。

主表：

```
E:/R/Nerve/results/
  00_READ_ME.txt
  00_Q1_PROTOCOL.txt
  01_Q1_literature_ligands.csv
  01_Q3_literature_repellents.csv
  00_expression_near_vs_far/          # 每人近/远肿瘤 RNA 表达，先看这里
    eligible_patients.csv
    pathologist_nerve/TNBC*/near_and_far_tumor_RNA_expression.csv
    schwann_neighborhood/COMBINED_mean_logCPM_near_and_far.csv
  04_COMBINED_near_vs_far_high_genes.csv  # 这些病人合在一起后的高表达基因
  03_Q3_READ_THIS_tumor_genes_low_may_promote_nerve.csv  # 问题3 低表达（p<0.05 且 FC<1）
  03_Q3_down_p005_FC_lt_1.csv
  03_Q3_down_p005_FC1.25.csv   # 可选更严
  02_Q2_up_p005_FC1.25.csv
  tumor_near_schwann_vs_far/     # 全队列 Schwann 邻域 DE
  TNBC50_tumor_near_nerve_vs_far/  # 病理神经，各病人单独
```

病理 Nerve 只和极少数 spot 重叠，全队列结论以 Schwann 邻域为准。不要把两个 NTC 或病人在空间分析里偷偷合并。
