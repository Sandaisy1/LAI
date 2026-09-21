# TG BRCA 细胞 RNA-seq 分析

针对 `NTC_rep0`、`NTC_rep1`、`TG_sh1`、`TG_sh5` 四个样品的 RNA-seq 分析。两个 NTC **不在 1-vs-1 比较里合并**。

人源乳腺癌转移的**临床组织** RNA/蛋白文献（必须同时有原发组织和转移组织，优先肺/骨/肝/脑）见 `docs/human_brca_metastasis_clinical_omics.md`。

AURORA US `GSE209998`、Sinn `GSE145752`、UNC RAP `GSE110590`、ConvertHER `GSE92977` 的配对原发–转移分析脚本分别为 `GSE209998_Human_breast.R`、`GSE145752_Human_breast.R`、`GSE110590_Human_breast.R`、`GSE92977_Human_breast.R`（数据默认 `E:/R/Human breast cancer/<GSE>`）。

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

## GSE209998 人源原发 vs 配对肺/骨（`GSE209998_Human_breast.R`）

把 `GSE209998_AUR_129_raw_counts.txt.gz` 和 `GSE209998_series_matrix.txt.gz` 放到 `E:/R/Human breast cancer/GSE209998`，然后：

```r
setwd("E:/R/Human breast cancer/GSE209998")
source("GSE209998_Human_breast.R")
```

只做 p < 0.05 且 FC > 1 / 1.25，按患者一一对应（热图列顺序为 原发_i、转移_i）。不做 top50–300。结果在同目录 `results_GSE209998_Human_breast/`。

## GSE145752 人源原发 vs 配对肺/胸膜（`GSE145752_Human_breast.R`）

Sinn 2020 NanoString 269 基因面板。把 `GSE145752_series_matrix.txt.gz`（建议同时留 `GSE145752_RAW.tar`）放到 `E:/R/Human breast cancer/GSE145752`，把脚本也拷到该目录，然后：

```r
setwd("E:/R/Human breast cancer/GSE145752")
source("GSE145752_Human_breast.R")
```

57 对一一对应：肺转移 8 对（patient 14, 20, 31, 38, 43, 44, 55, 88），胸膜 49 对，**无骨转移**（骨分析会写成 SKIPPED）。器官特异改为肺 vs 胸膜。只做 p < 0.05 且 FC > 1 / 1.25，不做 top50–300。结果在同目录 `results_GSE145752_Human_breast/`。

## GSE110590 人源原发 vs 配对肺/骨（`GSE110590_Human_breast.R`）

Siegel *JCI* 2018 UNC RAP。把 `GSE110590_RAP_A16_log2.sne.tsv.gz` 和两个 series matrix（`GSE110590-GPL11154_series_matrix.txt.gz`、`GSE110590-GPL16791_series_matrix.txt.gz`）放到 `E:/R/Human breast cancer/GSE110590`，把脚本也拷到该目录，然后：

```r
setwd("E:/R/Human breast cancer/GSE110590")
source("GSE110590_Human_breast.R")
```

矩阵已是 RSEM 上分位数标准化 log2。按患者一一对应：肺约 10 对（A1, A2, A4, A7, A11, A12, A15, A20, A26, A28），骨约 4 对（A1, A7, A11, A12）。A8 有转移无原发，不进配对。无「只转骨」患者，器官特异以配对 DE 差集为主，原发倾向比较用只转肺 vs 肺+骨。只做 p < 0.05 且 FC > 1 / 1.25，不做 top50–300。结果在同目录 `results_GSE110590_Human_breast/`。

## GSE92977 人源原发 vs 配对肺/骨（`GSE92977_Human_breast.R`）

Cejalvo ConvertHER 2017 NanoString。把 `GSE92977_series_matrix.txt.gz` 和 `GSE92977_raw_data.txt.gz` 放到 `E:/R/Human breast cancer/GSE92977`，把脚本也拷到该目录，然后：

```r
setwd("E:/R/Human breast cancer/GSE92977")
source("GSE92977_Human_breast.R")
```

123 对一一对应（Patient N 原发对 Patient N 转移）。肺 7 对（1, 11, 33, 49, 53, 80, 114），骨 16 对（6, 12, 17, 18, 21, 22, 23, 24, 27, 31, 35, 36, 40, 60, 93, 96），无人同时有肺和骨。管家基因标准化后的 ~105 基因面板。只做 p < 0.05 且 FC > 1 / 1.25，不做 top50–300。结果在同目录 `results_GSE92977_Human_breast/`。

