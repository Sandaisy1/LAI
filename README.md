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

## 血清蛋白质组：两组病人丰度排名气泡图

输入是 DIA-NN 的 `report.pg_matrix`（蛋白组定量）。若该文件蛋白行过少，脚本会用 `report.pr_matrix` 按 `Protein.Group` 取肽段强度中位数，聚合到蛋白。不要用 Cuffdiff RNA-seq 流程处理这些矩阵。

同一目录还需要 `sample_annotation.csv`。若没有且矩阵恰好两列样品（当前数据为 `GP_WJZ_11`、`GP_WJZ_18`），脚本会自动写出该表并按 1-vs-1 继续（不算 p 值）。也可手动复制 `serum_proteomics/sample_annotation.GP_WJZ.csv` 为 `sample_annotation.csv`。默认读取 `E:\天府\实验管理\课题\赵章寻\血清蛋白质组学`（可用 `SERUM_PROTEOMICS_DIR` 覆盖；该盘不存在时回退仓库内 `serum_proteomics/`）。矩阵、注释和结果都在这个目录：

```
sample_annotation.csv   # 列：sample,group；必须恰好两组病人
report.pg_matrix        # DIA-NN 蛋白组矩阵（TSV）
report.pr_matrix        # 可选；pg_matrix 不可用时回退
```

```r
# 错误示范：当前目录没有这个文件时会报“无法打开链接”
# source("serum_proteomics_ranked_bubble.R")

# 推荐：整段粘贴，或 source 启动器（会 setwd 到数据目录）
source("run_serum_proteomics_bubble.R", encoding = "UTF-8")
```

若启动器也不在当前目录，在 R 里先切到**仓库根目录**（能看到 `serum_proteomics_ranked_bubble.R` 的地方），或把该 `.R` 复制到数据目录后再：

```r
setwd("E:/天府/实验管理/课题/赵章寻/血清蛋白质组学")
source("serum_proteomics_ranked_bubble.R", encoding = "UTF-8")
```

会去 `E:/天府/实验管理/课题/赵章寻/血清蛋白质组学` 读矩阵，结果写在同目录 `results/`。

```bash
python3 serum_proteomics_ranked_bubble.py
```

按 **两样品平均丰度** 降序排名。横轴是丰度排名，纵轴是蛋白丰度值，气泡大小一致。默认 top20 / top30 / top50 以及全部蛋白：

```
results/serum_proteomics_bubble/
  protein_abundance_ranking.csv
  top20_abundance_rank_bubble.pdf|.png
  all_abundance_rank_bubble.pdf|.png
```

这是蛋白丰度排名气泡图（两样品平均，等大气泡），不是 GO/KEGG 气泡图，也不是两组 FC 图。
