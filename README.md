# TG BRCA 细胞 RNA-seq：按列出的 GO 通路分析表达

四个样品：`NTC_rep0`、`NTC_rep1`、`TG_sh1`、`TG_sh5`。两个 NTC **不在 1-vs-1 比较里合并**。

本流程对列出的 GO 做通路表达（热图/分数），气泡图的统计量来自全基因组 enrichGO 再抽取，不在自选通路上重算 p。

## 数据

默认读取 `E:/R/TG_BRCA/TG`。输入为 Cuffdiff：

- `genes.read_group_tracking` 的 **FPKM** 列（首选；Cuffdiff 已按长度和深度标准化）
- 否则 `genes.fpkm_tracking`
- 不用 `raw_frags` / `genes.count_tracking` 再做 DESeq2
- `metastasis_custom_genes.txt`（每行：通路名称 + `GO:#######`，不要只写基因符号）

## 运行

```r
setwd("E:/R/TG_BRCA/TG")
source("TG_RNAseq_pipeline.R")                 # 比较 1–4
source("TG_RNAseq_TGsh_mean_vs_NTC_reps.R")    # 比较 5–6
```

先过滤低表达，再用 Cuffdiff **FPKM** 做 `log2(x+1)` 后做六组比较。不再读 fragment count，也不做 DESeq2 size factor。不要用未过滤的值算 FC。

## 六组比较（各自单独出图）

1. 四个 1-vs-1：`TG_sh1 vs NTC_rep0`、`TG_sh5 vs NTC_rep0`、`TG_sh1 vs NTC_rep1`、`TG_sh5 vs NTC_rep1`（无 p，不伪造）
2. `mean(TG_sh1, TG_sh5)` vs `mean(NTC_rep0, NTC_rep1)`（2-vs-2，limma 估 p）
3. 相对 `NTC_rep0` 的共同上调（sh1 与 sh5 交集）
4. 相对 `NTC_rep1` 的共同上调
5. `mean(TG_sh1, TG_sh5)` vs `NTC_rep0`（不要混入 NTC_rep1）
6. `mean(TG_sh1, TG_sh5)` vs `NTC_rep1`（不要混入 NTC_rep0）

## 两套分层（每组都跑）

用 **p 值**（不是 padj），一律 **p < 0.05**：

1. 上调 FC ≥ 1 / 1.25 / 1.5 / 2
2. 上调 top 50 / 75 / 100 / 150 / 200 / 250 / 300

不再跑 p<0.01，也不再单独出 AllDE/UpDE/DownDE。无法估 p 时仍按 FC/排名分层，写 `NO_PVALUE.txt`，不伪造 p。

每个非空子集：差异基因表、火山图、热图。`GO/` 里有全库 BP/CC/MF 表、气泡图，以及 BP/CC/MF 合在一张上的柱状图（x 轴 **-lgP (p.adjust)**，红=BP、蓝=CC、绿=MF；每个 ontology 取前 10/15/20）。气泡图**先做全基因组 `enrichGO`**（`minGSSize=1`，很细的通路也会测；`maxGSSize=500`），再抽出 `metastasis_custom_genes.txt` 中通路的 GeneRatio、p.adjust、Count（不在自选通路上重新校正 p）。最终气泡图只保留 **p.adjust < 0.2**，再按 **GeneRatio 从大到小**取前 15 与前 20。柱状图同样只保留 p.adjust < 0.2，按显著性排序。全库完整表仍写出，不改 p.adjust。

全库 GO 表在各子集的 `GO/`（`*_ORA_GO_BP.csv` 等）；抽出的通路在 `CustomGO/`。每张气泡图会再导出 `*_plotdata.csv` / `.xlsx`。气泡大小和坐标字体在 `TG_RNAseq_pipeline.R` 开头改（`bubble_size_min` / `bubble_size_max`、`axis_text_y_size` / `axis_text_x_size`）。只重画已有图、不重跑 enrichGO：

```r
options(tg.rnaseq.restyle_only = TRUE)
source("TG_RNAseq_pipeline.R")
```

## 通路表达（不改全库 p 值）

`results/00_PathwayExpression/`：

- 每条列出 GO 的基因热图
- 通路分数热图（ssGSEA；没有 GSVA 包则用基因 mean z-score）
- 各比较目录 `CustomGO/` 还有总的 listed-GO mean log2FC 柱状图：上调+下调全部，以及只含上调通路的 top 10 / 15 / 20
- 各比较目录 `PathwayScore/`：ssGSEA 与 mean z 气泡图。x 轴是该比较的通路分数差（处理 − 对照），点大小是通路基因数，颜色是分数差。出全部通路，以及分数升高的 top 10 / 15 / 20。1-vs-1 不伪造 p。共同上调比较用 `mean(TG_sh1, TG_sh5)` vs 对应 NTC（ssGSEA 不能用基因交集）

```
results/TG_sh1_vs_NTC_rep0/CustomGO/TG_sh1_vs_NTC_rep0_pathway_mean_log2FC.pdf
results/TG_sh1_vs_NTC_rep0/CustomGO/TG_sh1_vs_NTC_rep0_pathway_mean_log2FC_up_top10.pdf
results/TG_sh1_vs_NTC_rep0/CustomGO/TG_sh1_vs_NTC_rep0_pathway_mean_log2FC_up_top15.pdf
results/TG_sh1_vs_NTC_rep0/CustomGO/TG_sh1_vs_NTC_rep0_pathway_mean_log2FC_up_top20.pdf
results/TG_sh1_vs_NTC_rep0/p0.05/FoldChange/FC_1/GO/p0.05_FC_1_ORA_GO_BP_CC_MF_barplot_top15.pdf
results/TG_sh1_vs_NTC_rep0/p0.05/FoldChange/FC_1/CustomGO/p0.05_FC_1_ORA_CustomGO_BP_CC_MF_barplot_top15.pdf
results/TG_sh1_vs_NTC_rep0/PathwayScore/TG_sh1_vs_NTC_rep0_ssgsea_dotplot_up_top15.pdf
results/TG_sh1_vs_NTC_rep0/PathwayScore/TG_sh1_vs_NTC_rep0_mean_z_dotplot_up_top15.pdf
results/00_PathwayExpression/pathway_score_heatmap_ssgsea.pdf
```
