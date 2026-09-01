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

# 流式降维（免疫细胞亚群，H vs EV，P1 / P2 / P3）

这是 **`E:/R/fuction of cell`** 的免疫细胞亚群降维。找 His+ 靶细胞是另一套实验，走 `ICI_Flow_dimred_pipeline.R`（`E:/R/Internation cell immune`）。**JY-NNK / JY-EVNK** 走 `JY_Flow_dimred_pipeline.R`（`E:/R/fuction of cell-ljy`）。**JZ-AB / JZ-EVB** 走 `JZ_Flow_dimred_pipeline.R`（`E:/R/fuction of cell-wjz`）。几套不要混用。

小鼠流式只用 `*_unmixed.fcs`。分析结果写在同一目录的 `results_flow/`。三个 panel **分开**做 UMAP/tSNE，比较 **ZZX_EV**（EV1/2/3，每管两个技术重复 EV1-1/EV1-2）与 **ZZX_H**（H1/2/3，H1-1/H1-2）。统计：技术重复先平均；三个生物学重复里去掉离中位数更远的最大值或最小值，用剩下 n=2 做均值/SD/检验。tSNE 仍用全部细胞。

```r
setwd("E:/R/fuction of cell")
# 把 Flow_dimred_pipeline.R、Flow_dimred_all_subsets.R、Flow_dimred_trajectory.R 与 flow_panel_map.json 放在该目录（或仓库根）
source("Flow_dimred_pipeline.R")
```

文件名：`ZZX_EV1-1_P1_unmixed.fcs`、`EV1-2_P3_unmixed.fcs`、`H-1-1_P2_unmixed.fcs` 都可以（`ZZX_` 可有可无）。旧名 `EV1_P1_unmixed.fcs`（无技术重复号）仍能解析。文件夹名 `ZZX_EV` / `ZZX_H` 也可用来认组别。

若 Windows 把 json 显示成“文本文档”，真实文件名多半是 `flow_panel_map.json.txt`。关掉“隐藏已知文件类型的扩展名”后改名为 `flow_panel_map.json`；新版脚本两种名字都能读。

每个 panel 的图在 `E:/R/fuction of cell/results_flow/P1/`、`P2/`、`P3/`（PDF + PNG）：

- `*_H_vs_EV_tSNE_major_split` / `*_UMAP_major_split`：**图1**，左 EV、右 H，共用联合 tSNE，点按 **免疫大类**（CD4 T / CD8 T / NK / NKT / B / Myeloid）着色；图注在右侧，画布加宽以免裁切。联合降维会加重 CD4/CD8/CD3 等谱系标志，避免 naive CD4 与 naive CD8 因共享 CD62L/CD44 而揉在一起。圈门是 CD4+CD8− / CD4−CD8+ 象限，不是「谁的荧光更高」。
- `*_H_vs_EV_tSNE_lineage_split` / `*_UMAP_lineage_split`：同一套 embedding 上的**全体细胞、每一个细亚群**（不是只画六个大类）。
- `dimred_by_major/`：每个大类单独再降维，点按该大类的**细亚群**着色（同样 EV | H）。类内用高对比定性色，不要 CD8 全绿 / CD4 全红棕 / NK 全紫 / NKT 全黄橙。细胞少时点会放大、画布缩小，避免大空白上颜色糊成一团。
- `*_UMAP_by_group` / `*_tSNE_by_group`：EV vs H
- `*_UMAP_by_cluster` / `*_UMAP_by_lineage`
- `subset_stats/`：每个亚群一张图，上为 EV（黑）vs H（红）柱状图；下为 FlowJo 风格 2D 图。门是铺到坐标轴的完整象限/半平面（CD62L/CD44 标四个象限），EV 与 H **各自切阈值**。不是只框 10–90% 命中细胞的小矩形。
- `functional_state/`：NKT 以及 Naive/Unswitched/Switched/Atypical/Activated B **各自单独**看活化或耗竭。图按多列排（上柱状图、下对照在上处理在下的等高线）。P1 NKT：CD69 / IFN-γ / TNF-α / GZMB，以及 PD-L1 / LAG-3 / TIM-3。P2 各 B 亚群：CD86 / CD80 / CD40（该 panel 没有耗竭通道，不伪造）。JY/JZ 同样出在各自 `results_flow/`，组别色仍是对照黑、处理红。总降维已经跑完时不必重跑 UMAP，在对应数据目录 `source("Flow_dimred_functional_state.R")`（JY/JZ 用各自的 `JY_Flow_dimred_functional_state.R` / `JZ_Flow_dimred_functional_state.R`），它只读 `P1_cell_embeddings.csv` / `P2_cell_embeddings.csv`。
- `gating/<样品>/`：**每个 FCS 单独**的完整圈门图 + `*_per_sample_gate_cuts.csv`。圈门按每个样品单独做完同一套逻辑。P1 分层：naive / **T_CM** / **T_SCM**（CD27+ CD95+）/ **T_EM early·late** / **SLEC** / **MPEC** / **T_EFF** / **exhausted**（PD-L1、LAG-3、TIM-3），以及 CD69 活化。NK 在 CD11b **之前**用 NKp46 圈出，再按 CD69 / 杀伤 / 耗竭 / **CD27 vs CD11b**（immature、DP、mature）拆；NKT 再分成 CD4 NKT 与 DN NKT。NKG2D 不当亚群，与 IFN-g/TNF-a/GZMB 一起出 MFI 表。P2：宽单核门 → CD45+ → CD19+ → IgD vs CD27 的 Naive / Unswitched / Switched；再分 MZ、Plasmablast、Plasma；CD40/CD80/CD86 出 MFI 表，不当第 1 层亚群。P3：CD3/CD19/NK1.1 大类 → 三阴髓系按 CD11B 分；嗜酸要 Siglec-F+CCR3+；肥大细胞 FceRI+CD200R3+；F4/80 high 巨噬 vs Ly6C hi 单核；CD11B− 上 CD11C+MHC-II+ 再分 cDC1/cDC2。
- `markers/`：各通道在 UMAP 上的着色
- `*_cluster_marker_heatmap`、`*_annotation_heatmap`（按注释亚群的标志物中位数，另有 `*_annotation_heatmap_major`）、`*_cluster_frequency_H_vs_EV`、`*_H_vs_EV_dimred_overview`

三个 panel 的抗体不同，**不能**拼成一张矩阵做联合 UMAP。`Flow_dimred_all_subsets.R` 在各自圈完亚群后，把频率汇总到 `results_flow/all_subsets/`（H vs EV 柱状图、堆叠组成、热图、H−EV 差值）。百分数是该 panel 管子里的比例，P1 的 B/髓系、P3 的 T/B/NK 只是 dump。主流程跑完会自动出这份总览；也可单独 `source("Flow_dimred_all_subsets.R")`。

每个 panel 的轨迹和降维同一套结构：`Flow_dimred_trajectory.R` → `results_flow/P1/trajectory/` 等。先画 **全体大类** 树（`P1_major_trajectory`），再画各大类的亚群树（P1 的 CD4/CD8/NK/NKT，P2 的 B，P3 的髓系）。坐标是该类内的 Component 1/2，点按细亚群（或大类）着色，黑线是分支骨架。根节点按惯例是 naive T、NK immature、Naive B 或 Ly6C hi 单核。也可单独 `source("Flow_dimred_trajectory.R")`。日志里 `skip trajectory (... subsets= 1)` 或 `skip major-class trajectory` 表示节点不够画树，不是分析失败。图注在右侧，不要被裁掉。

无 FCS 时可 `Sys.setenv(FLOW_DEMO = "1")` 导出演示图，不能当正式结果。

## Internation cell immune（另一套实验，找 His+ 靶细胞）

数据在 `E:/R/Internation cell immune`。染色是 **P1（T/NK）+ P2（B）+ P3（髓系）**，都带 **His-FITC**。主目的是找靶细胞。圈门思路与免疫亚群方案相同，但染色按 `ICI_flow_panel_map.json`，**不要改** `Flow_dimred_pipeline.R` / `flow_panel_map.json`。

ICI **P2** 是 7 色：L/D、CD45、CD19、CD27、IgG-RB613、IgM-RB780、His-FITC（没有 IgD / BLIMP / CD40）。圈门仍是 CD19+ 为 B 母门，再用 CD27×IgM/IgG 对应原 P2 的 naive / unswitched / switched；缺的项跳过，不跳过整个 panel。

组别：**ZZX-EV**（EV-1、EV-2、EV-3）vs **ZZX-H**（H-1、H-2、H-3），比较 H vs EV，n=3。

两套方案完全独立。把下面文件**整套**拷到 `E:/R/Internation cell immune`（不要去 `E:/R/fuction of cell` 找原流程，也不要 `source` `Flow_dimred_pipeline.R`）：

- `ICI_Flow_dimred_pipeline.R`（入口）
- `ICI_flow_engine.R`（本方案自己的分析函数库）
- `ICI_flow_panel_map.json`
- `ICI_Flow_dimred_all_subsets.R`
- `ICI_Flow_dimred_trajectory.R`

免疫亚群那套 `Flow_*` 继续只放在 `E:/R/fuction of cell`，两边互不干扰。

```r
setwd("E:/R/Internation cell immune")
source("ICI_Flow_dimred_pipeline.R")
```

QC 只去双联体和死细胞，**不要求 CD45+**，也 **不用 P1 紧淋巴门**，以免丢掉 His+ CD45− 靶细胞。活细胞里 His+ CD45− 标成 **His+ target**；His+ CD45+ 仍按原来的 T/NK/B/髓系亚群命名，并另出各亚群内 His+ 比例。

降维和轨迹与免疫亚群方案同一套：`*_major_split` 按大类着色（含 His+ target），`*_lineage_split` 是全体细亚群，`dimred_by_major/` 再按大类重降维（高对比色）；轨迹先画全体大类树再画各类亚群树。每个 panel 出免疫细胞注释热图 `*_annotation_heatmap`（含 His+ target）。ICI P1 没有 CD19，缺通道按阴性处理，不要再出现 `NAs are not allowed in subscripted assignments`；缺的染色通道该项可以不分析，但不能跳过该 panel 其余图。结果在同目录 `results_flow/`，靶细胞表和图在 `target_His/`。

## JY（第三套实验，免疫亚群降维，比较 JY-NNK / JY-EVNK）

数据在 `E:/R/fuction of cell-ljy`。分析与上面的免疫亚群降维一致（同一套 P1/P2/P3 染色、圈门、图件），只换目录和组别。

组别：**JY-EVNK**（生物学重复 EVNK-1、EVNK-2、EVNK-3；每个生物学重复两个技术重复，如 EVNK1-1、EVNK1-2）与 **JY-NNK**（NNK-1、NNK-2、NNK-3；技术重复如 NNK1-1、NNK1-2）。比较是 **JY-NNK / JY-EVNK**（NNK 相对 EVNK）。统计同样是技术重复先平均，再去掉 1 个极端生物学重复，n=2。tSNE/UMAP 仍用全部管子。

四套方案完全独立。把下面文件**整套**拷到 `E:/R/fuction of cell-ljy`（不要去 `E:/R/fuction of cell` 找原流程，也不要 `source` `Flow_dimred_pipeline.R`、`ICI_*` 或 `JZ_*`）：

- `JY_Flow_dimred_pipeline.R`（入口）
- `JY_flow_engine.R`（本方案自己的分析函数库）
- `JY_flow_panel_map.json`
- `JY_Flow_dimred_all_subsets.R`
- `JY_Flow_dimred_trajectory.R`
- `JY_Flow_dimred_functional_state.R`（总结果已出时，只补 NKT/B 活化与耗竭）

```r
setwd("E:/R/fuction of cell-ljy")
source("JY_Flow_dimred_pipeline.R")
```

文件名：`EVNK1-1_P1_unmixed.fcs`、`JY-NNK-2-2_P3_unmixed.fcs`、`NNK-3-1_P2_unmixed.fcs` 都可以（`JY-` 可有可无）。**先匹配 EVNK 再匹配 NNK**。结果在同目录 `results_flow/`，比较标签是 **JY-NNK / JY-EVNK**。

## JZ（第四套实验，免疫亚群降维，比较 JZ-AB / JZ-EVB）

数据在 `E:/R/fuction of cell-wjz`。分析与上面的免疫亚群降维一致（同一套 P1/P2/P3 染色、圈门、图件），只换目录和组别。

组别：**JZ-EVB**（生物学重复 EVB-1、EVB-2、EVB-3；每个生物学重复两个技术重复，如 EVB1-1、EVB1-2）与 **JZ-AB**（AB-1、AB-2、AB-3；技术重复如 AB1-1、AB1-2）。比较是 **JZ-AB / JZ-EVB**（AB 相对 EVB）。统计同样是技术重复先平均，再去掉 1 个极端生物学重复，n=2。tSNE/UMAP 仍用全部管子。

四套方案完全独立。把下面文件**整套**拷到 `E:/R/fuction of cell-wjz`（不要去 `E:/R/fuction of cell` 或 `cell-ljy` 找原流程，也不要 `source` `Flow_dimred_pipeline.R`、`ICI_*` 或 `JY_*`）：

- `JZ_Flow_dimred_pipeline.R`（入口）
- `JZ_flow_engine.R`（本方案自己的分析函数库）
- `JZ_flow_panel_map.json`
- `JZ_Flow_dimred_all_subsets.R`
- `JZ_Flow_dimred_trajectory.R`
- `JZ_Flow_dimred_functional_state.R`（总结果已出时，只补 NKT/B 活化与耗竭）

```r
setwd("E:/R/fuction of cell-wjz")
source("JZ_Flow_dimred_pipeline.R")
```

文件名：`EVB1-1_P1_unmixed.fcs`、`JZ-AB-2-2_P3_unmixed.fcs`、`AB-3-1_P2_unmixed.fcs` 都可以（`JZ-` 可有可无）。**先匹配 EVB 再匹配 AB**。结果在同目录 `results_flow/`，比较标签是 **JZ-AB / JZ-EVB**。

若日志出现「没有匹配到任何分析通道」，多半是 Cytek 通道名带 `-A`（如 `BUV496-A`）或 desc 为空，不是文件没找到。用最新脚本再跑；仍失败时日志会列出通道名，并在 `results_flow/00_logs/` 写 `*_channels.csv`。

读真实 FCS 需要 `flowCore`。若日志出现「读 FCS 需要 flowCore」，先在 R 里装好再重新 `source`：

```r
install.packages("BiocManager")
BiocManager::install("flowCore")
```

UMAP/tSNE 需要 `install.packages(c("uwot", "Rtsne"))`。缺这些包时脚本会退回 PCA，图标题会标明 fallback。

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
