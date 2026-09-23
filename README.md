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

## Wang 2024 空间转录组：远神经肿瘤 vs 近神经肿瘤

数据在 `E:/R/Nerve`。比较是 **远神经肿瘤组织 vs 近神经肿瘤组织**（上调 = 远神经肿瘤更高）。FC **只做 >1、1.25、1.5**，没有 FC=2、没有 topN、没有下调。

请 source `Wang_ST_nerve_infiltration.R`，不要 source `Nerve_RNA_primary_lung_tropism.R`（那是另一套 GPL96 分析）。

```
E:/R/Nerve/results/
  00_请先看这里_单个病人和上调基因.txt
  00_eligible_single_patients/                 # 符合条件的单个病人
    符合条件病人_打开这个.csv
    INDEX_eligible_patients.csv
    schwann_neighborhood/TNBC32/near_and_far_tumor_RNA_expression.csv
    schwann_neighborhood/TNBC32/DE_far_tumor_vs_near_tumor/up_FC_gt_1/
    schwann_neighborhood/TNBC32/DE_far_tumor_vs_near_tumor/up_FC_1.25/
    schwann_neighborhood/TNBC32/DE_far_tumor_vs_near_tumor/up_FC_1.5/
    pathologist_nerve/TNBC50/                  # 病理神经，病人很少
  01_combined_far_tumor_vs_near_tumor/         # 这些人合在一起
    upregulated_far_tumor_vs_near_tumor_FC_gt_1.csv
    upregulated_far_tumor_vs_near_tumor_FC_1.25.csv
    upregulated_far_tumor_vs_near_tumor_FC_1.5.csv
```

旧文件夹 `tumor_near_schwann_vs_far/`、`TNBC50_tumor_near_nerve_vs_far/`、`FoldChange/FC_2/`、`TopRank/` 是上一版近 vs 远，重新跑会删掉。

```r
setwd("E:/R/Nerve")
Sys.setenv(WANG_ST_DIR = "E:/R/Nerve")
source("Wang_ST_nerve_infiltration.R")
```

不需要 `TG_RNAseq_pipeline.R`。先 **p < 0.05**，再上调 **FC > 1、1.25、1.5**。原 Cuffdiff 六组比较不变。

真正必需的是 `Robjects/`（`counts` + `annotsBySpot`；先解压 `Robjects.tar`）。

`ids.RDS` 可选。Windows 解 `Clinical.tar` 时常把 `Clinical/ids.RDS` 展成根目录的 `Clinicalids.RDS`；`Clinical.RDS` / `Clinical.xlsx` 是临床表。没有 ids 也能跑。

病理 Nerve 只和极少数 spot 重叠，全队列结论以 Schwann 邻域为准。

### 指定 32 个病人的综合分析（单独脚本）

全队列跑完后，另跑 `Wang_ST_far_vs_near_selected_patients.R`。只综合病人 2、3、4、6、11、13、14、15、20、22、27、28、30、31、37、39、51、52、53、56、61、62、67、69、74、79、81、85、86、90、93、94。比较仍是 **远神经肿瘤 vs 近神经肿瘤**，上调只做 **FC>1、1.25**（没有 1.5）。不覆盖上面的 `00_` / `01_` 结果。

```r
setwd("E:/R/Nerve")
Sys.setenv(WANG_ST_DIR = "E:/R/Nerve")
source("Wang_ST_far_vs_near_selected_patients.R")
```

```
E:/R/Nerve/results/02_selected32_far_tumor_vs_near_tumor/
  00_请先看这里.txt
  INDEX_requested_vs_used.csv
  upregulated_far_tumor_vs_near_tumor_FC_gt_1.csv
  upregulated_far_tumor_vs_near_tumor_FC_1.25.csv
  03_distance_to_Schwann/
    INDEX_FC_1.25_each_gene_each_patient.csv
    FC_1.25_each_gene/GENE/TNBC病人_GENE_high_vs_low_vs_Schwann_distance.pdf
```

近/远仍是主脚本那套：肿瘤 spot 到施旺 **≤2 spot 为近，>4 spot 为远**。1 spot ≈ 0.1 mm。

距离图只画 **`upregulated_far_tumor_vs_near_tumor_FC_1.25.csv` 里那几个上调基因**（你这份大约 10 个，例如 PCYOX1L、ABHD3、CEP152）。**每个基因、每个病人各一张**高/低表达 vs 施旺距离图，不把病人拼成一张总图。

不要看 `03_distance_to_Schwann/genes/`。那里的 ARTN、NGF、CXCL12 是上一版误画的文献配体，不是这批上调基因；更新脚本后重跑会删掉这个文件夹。
