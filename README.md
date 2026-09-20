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

## GSE165393 小鼠原位 vs 肺/骨转移

数据目录：`E:/R/Mouse Breast`，放入 GEO 补充文件 `GSE165393_AllTissues_TPM.csv.gz`。把 `GSE165393_Mouse_breast.R` 拷到该目录后：

```r
setwd("E:/R/Mouse Breast")
source("GSE165393_Mouse_breast.R")
```

脚本回答三问（低表达促进转移 = 转移灶相对原位 MFP 下调）：

1. 肺和骨髓都下调的共享基因 → `results_GSE165393/01_shared_lung_and_bone_down/` 与 `07_summary/Q1_*`
2. 只肺或只骨下调 → `02_lung_specific_down/`、`03_bone_specific_down/` 与 `07_summary/Q2_*`
3. 神经浸润拆成三套 marker（施旺细胞、神经营养因子、轴突导向）分别打分、找负相关基因，并比较三套分数与肺/骨转移的关系 → `06_neural_invasion/` 与 `07_summary/Q3_*`

先看 `results_GSE165393/07_summary/`。入选为 limma **p < 0.05**，下调 **FC < 1** 即可，并再出 **FC < 1/1.25**。BM 是骨髓来源细胞系，不是皮质骨灶。

## GSE273439 小鼠原位 vs 配对肺（Visium）

数据目录：`E:/R/Mouse Breast/GSE273439`，放入并解压 `GSE273439_RAW.tar`（四个样品的 MTX / 坐标 / PNG）。把 `GSE273439_Mouse_breast.R` 拷到该目录后：

```r
setwd("E:/R/Mouse Breast/GSE273439")
source("GSE273439_Mouse_breast.R")
```

**一一对应：** Mouse627 原位只对 Mouse627 肺，Mouse628 只对 Mouse628 肺，两只鼠不混成一组。入选为 **p < 0.05**，下调倍数 **原位/肺 ≥ 1 和 ≥ 1.25**。只做这两档 FC，**没有 Top50–300**。

1. 配对肺下调 → `results_GSE273439/01_lung_down_Mouse627/`、`01_lung_down_Mouse628/`、`01_lung_down_both_mice/` 与 `08_summary/Q1_*`
2. 肺转移灶特异下调（相对未受累肺实质）→ `02_lung_foci_specific_down/`。**本套数据没有骨**，`03_bone_not_in_GSE273439/00_NOTE.txt` 说明原因
3. 原位 TM 上三套神经分数（施旺 / 神经营养 / 轴突导向）及负相关基因 → `06_neural_invasion/`
4. 三种神经分数分别 vs 配对肺 → `07_neural_vs_lung/` 与 `08_summary/Q4_*`

先看 `results_GSE273439/08_summary/`。spot 不是独立生物学重复，主结论用「两鼠都过阈值」。跳过富集可设环境变量 `GSE273439_SKIP_ENRICHMENT=1`。

## GSE146012 小鼠原位 vs 配对肺（4T1 GFP+ bulk）

数据目录：`E:/R/Mouse Breast/GSE146012`，放入 `GSE146012_RNA_seq_logtransformed_count.txt.gz`（浏览器解压成 `.txt` 也可以）。把 `GSE146012_Mouse_breast.R` 拷到该目录后：

```r
setwd("E:/R/Mouse Breast/GSE146012")
source("GSE146012_Mouse_breast.R")
```

矩阵已经是 DESeq2 log，**不要再跑 DESeq2**。**一一对应：** Tumor-1 只对 MFP-Met-1（及 TVI-Met-1），Tumor-2 只对 Met-2；不要把两只原发和三只肺混成一组，也不要把 MFP 与 TVI 合并。Met-3 没有 Tumor-3，只相对原发均值作补充、不算真配对。入选为 **p < 0.05**（两对配对 limma；1-vs-1 无法估计 p 则只按 FC），下调倍数 **原位/肺 ≥ 1 和 ≥ 1.25**。只做这两档 FC，**没有 Top50–300**。

1. 配对肺下调 → `results_GSE146012/01_lung_down_pair1_Tumor1_vs_MFPMet1/`、`01_lung_down_pair2_Tumor2_vs_MFPMet2/`、`01_lung_down_both_pairs_MFP/`。TVI 路径单独在 `01b_*`。**本套数据没有骨**，`03_bone_not_in_GSE146012/00_NOTE.txt` 说明原因
2. 路径特异（MFP 下调而 TVI 未同时下调；不能做骨器官特异）→ `02_MFP_lung_specific_vs_TVI/`
3. 原位 Tumor 上三套神经分数（施旺 / 神经营养 / 轴突导向）及负相关基因 → `06_neural_invasion/`
4. 三种神经分数分别 vs 配对 MFP 肺 → `07_neural_vs_lung/` 与 `08_summary/Q4_*`

先看 `results_GSE146012/08_summary/`。跳过富集可设 `GSE146012_SKIP_ENRICHMENT=1`。

## GSE54773 小鼠原位衍生系 vs 配对肺/脑（芯片）

数据目录：`E:/R/Mouse Breast/GSE54773`，放入 `GSE54773_series_matrix.txt.gz`（建议同时放 `GPL6246.annot.gz`）。把 `GSE54773_Mouse_breast.R` 拷到该目录后：

```r
setwd("E:/R/Mouse Breast/GSE54773")
source("GSE54773_Mouse_breast.R")
```

矩阵已经是 RMA log2，**不要再跑 DESeq2**。**一一对应：** T2-1 只对 LM2-1 / BM2-1，T2-2 对 2 号，T2-3 对 3 号。GEO 写明 9 只独立小鼠，脚本仍按编号配对。入选为 **p < 0.05**（三对配对 limma；1-vs-1 无法估计 p 则只按 FC），下调倍数 **原位/转移 ≥ 1 和 ≥ 1.25**。只做这两档 FC，**没有 Top50–300**。

1. 配对肺下调 → `results_GSE54773/01_lung_down_pair1_T2_1_vs_Lung_1/` … `01_lung_down_all_pairs/`。脑转移在 `01b_brain_down_*`。**本套数据没有骨**，`03_bone_not_in_GSE54773/00_NOTE.txt` 说明原因
2. 器官特异（肺下调而脑未同时下调，或反过来）→ `02_lung_specific_vs_brain/`、`02b_brain_specific_vs_lung/`
3. 原位 T2 上三套神经分数（施旺 / 神经营养 / 轴突导向）及负相关基因 → `06_neural_invasion/`
4. 三种神经分数分别 vs 配对肺 → `07_neural_vs_lung/` 与 `08_summary/Q4_*`

先看 `results_GSE54773/08_summary/`。跳过富集可设 `GSE54773_SKIP_ENRICHMENT=1`。

## 小鼠转移组学对照文献

本仓库样品是人 BRCA 细胞 TG 敲低 RNA-seq，不是小鼠组织。若要用「原位 vs 肺/骨/肝/脑」公开数据做签名 overlap，见：

- `literature/mouse_breast_cancer_metastasis_datasets.md`（文献说明与选用建议）
- `literature/mouse_breast_cancer_metastasis_datasets.csv`（登录号、器官、数据类型表）

没有一份公开数据同时覆盖原位 + 四器官配对 bulk RNA-seq；多器官比较需组合 GSE165393、GSE146012、GSE37975、GSE238214、GSE54773，或用 E-MTAB-16621 的四器官 scRNA-seq（无原位瘤）。
