# TG BRCA 细胞 RNA-seq：按列出的 GO 通路分析表达

四个样品：`NTC_rep0`、`NTC_rep1`、`TG_sh1`、`TG_sh5`。两个 NTC **不在 1-vs-1 比较里合并**。

本流程对列出的 GO 做通路表达（热图/分数），气泡图的统计量来自全基因组 enrichGO 再抽取，不在自选通路上重算 p。

## 数据

默认读取 `E:/R/TG_BRCA/TG`。输入为 Cuffdiff：

- `genes.read_group_tracking`（首选）
- `genes.count_tracking` / `genes.fpkm_tracking`
- `metastasis_custom_genes.txt`（每行：通路名称 + `GO:#######`，不要只写基因符号）

## 运行

```r
setwd("E:/R/TG_BRCA/TG")
source("TG_RNAseq_pipeline.R")                 # 比较 1–4
source("TG_RNAseq_TGsh_mean_vs_NTC_reps.R")    # 比较 5–6
```

先过滤低表达，再标准化（有 count 用 DESeq2 size factor；仅 FPKM 则 log2 后分位数标准化），然后做六组比较。不要用原始值算 FC。

## 六组比较（各自单独出图）

1. 四个 1-vs-1：`TG_sh1 vs NTC_rep0`、`TG_sh5 vs NTC_rep0`、`TG_sh1 vs NTC_rep1`、`TG_sh5 vs NTC_rep1`（无 p，不伪造）
2. `mean(TG_sh1, TG_sh5)` vs `mean(NTC_rep0, NTC_rep1)`（2-vs-2，limma 估 p）
3. 相对 `NTC_rep0` 的共同上调（sh1 与 sh5 交集）
4. 相对 `NTC_rep1` 的共同上调
5. `mean(TG_sh1, TG_sh5)` vs `NTC_rep0`（不要混入 NTC_rep1）
6. `mean(TG_sh1, TG_sh5)` vs `NTC_rep1`（不要混入 NTC_rep0）

## 五套分层（每组都跑）

用 **p 值**（不是 padj）：

1. `p < 0.05` 且上调 FC ≥ 1 / 1.25 / 1.5 / 2
2. `p < 0.05` 后上调 top 50 / 75 / 100 / 150 / 200 / 250 / 300
3. `p < 0.01` 且上调 FC ≥ 1 / 1.25 / 1.5 / 2
4. `p < 0.01` 后上调 top 50–300
5. `AllDE`：全部上调 + 下调（能估 p 时先滤 `p < 0.05`）。同时单独出 `UpDE/`（只上调）和 `DownDE/`（只下调）

无法估 p 时仍按 FC/排名分层，写 `NO_PVALUE.txt`，不伪造 p。

每个非空子集：差异基因表、火山图、热图。`GO/` 里有全库 BP/CC/MF 表和气泡图（`*_dotplot_top15.pdf` 与 `top20.pdf`）。气泡图**先做全基因组 `enrichGO`**，再抽出 `metastasis_custom_genes.txt` 中通路的 GeneRatio、p.adjust、Count（不在自选通路上重新校正 p）。最终气泡图只保留 **p.adjust < 0.05**，再按 **GeneRatio 从大到小**取前 15 与前 20（同一套排序，20 比 15 只在下面多 5 条）。全库完整表仍写出，不改 p.adjust。

全库 GO 表在各子集的 `GO/`（`*_ORA_GO_BP.csv` 等）；抽出的通路在 `CustomGO/`。

## 通路表达（不改全库 p 值）

`results/00_PathwayExpression/`：

- 每条列出 GO 的基因热图
- 通路分数热图（ssGSEA；没有 GSVA 包则用基因 mean z-score）
- 各比较目录下还有该比较中每条 GO 的 mean log2FC 柱状图

```
results/TG_sh1_vs_NTC_rep0/p0.05/FoldChange/FC_1/CustomGO/p0.05_FC_1_ORA_CustomGO_dotplot.pdf
results/00_PathwayExpression/pathway_score_heatmap.pdf
```
