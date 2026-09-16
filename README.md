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

# TIF / 血清蛋白质组（T vs N）

数据在 `E:/R/Protein TIF serum`，输入是 DIA-NN 蛋白矩阵（不要用 `*.pr_matrix` 做蛋白水平差异）。Windows 资源管理器常把后缀藏起来，文件类型若是「TSV 文件」，实际名字是 `TIF_report.pg_matrix.tsv` / `Serum_report.pg_matrix.tsv`，脚本两种都认：

- `TIF_report.pg_matrix`：组织间质液，组别 `N`、`T`、`T6`
- `Serum_report.pg_matrix`：血清，`N`=`N1,N3,N7`，`T`=`T1,T3,T5`，`T6`=`T6-1,T6-2,T6-3`

只做三件事，**T6 不进入 T vs N**，TIF 与血清分开标准化：

1. TIF `T vs N`：差异蛋白、火山图、上调 GO、上调 KEGG
2. 血清 `T vs N`：同上
3. 先找 TIF T vs N **上调**蛋白，再去掉也在血清 T vs N **上调名单**里的蛋白

有两份**互不 source** 的脚本，拷到数据目录后在 R 控制台运行（不要输入 `Rscript`）：

**热图版**（结果 `results_protein/`）：

```r
setwd("E:/R/Protein TIF serum")
source("Protein_TIF_Serum_pipeline.R")
```

**独立排名图版**（结果 `results_protein_standalone/`，第 3 条画排名图）：

```r
setwd("E:/R/Protein TIF serum")
source("Protein_TIF_Serum_TVsN_standalone.R")
```

独立脚本不读取、不修改 `Protein_TIF_Serum_pipeline.R`。

若在 **Windows 命令提示符** 或 PowerShell 里运行，才用：

```bat
Rscript run_protein_tif_serum.R "E:/R/Protein TIF serum"
```

没有真实矩阵时，仓库里的 `demo_protein_tif_serum/` 可先跑通流程（演示数据，不是实验结果）。

结果目录：

```
results_protein/                  # Protein_TIF_Serum_pipeline.R
  TIF_T_vs_N/
  Serum_T_vs_N/
  TIF_specific_vs_Serum/          # 热图
results_protein_standalone/       # Protein_TIF_Serum_TVsN_standalone.R
  TIF_T_vs_N/
  Serum_T_vs_N/
  TIF_specific_vs_Serum/          # 排名图
```

第 3 组：`TIF 上调名单 − 血清上调名单`。两边都上调的蛋白在 `TIF_up_AND_Serum_up_excluded.csv`。

- 热图版：`TIF_specific_vs_Serum/heatmap_TIF_up_absent_from_Serum_up` 只用 TIF 的 T、N
- 独立排名图版：`results_protein_standalone/TIF_specific_vs_Serum/rank_TIF_up_absent_from_Serum_up`（lollipop + rank vs log2FC）

列名识别失败时，可在数据目录放 `sample_map.csv`（列：`file,assay,group,replicate`），`file` 匹配原始列名即可。样本名会先匹配 `T6` 再匹配 `T`，避免把 `T6` 当成 `T`。
