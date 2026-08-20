# TG BRCA 细胞 RNA-seq 分析

针对 NTC、TG_sh1、TG_sh5 三组人乳腺癌细胞 RNA-seq 的差异表达与富集分析流程。

## 数据位置

默认读取 `E:/R/TG_BRCA/TG`。若该路径不存在，则使用当前工作目录。也可在运行前设置：

```r
Sys.setenv(TG_RNASEQ_DIR = "E:/R/TG_BRCA/TG")
```

优先读取同目录下 Cuffdiff 文件（已不再使用 `shTG(20220723).xlsx`）：

- `genes.read_group_tracking`（首选，含重复）
- `genes.count_tracking`
- `genes.fpkm_tracking`

Cuffdiff 的 `condition` 名需能对应到 `NTC`、`TG_sh1`、`TG_sh5`（例如 `NTC`、`shTG1`、`shTG5` 均可）。若缺少 TG_sh5，脚本会跳过 sh5 相关比较并写入日志。

## 运行

在 R 或 RStudio 中：

```r
setwd("E:/R/TG_BRCA/TG")   # 或本仓库目录
source("TG_RNAseq_pipeline.R")
```

首次运行会自动安装 CRAN / Bioconductor 依赖，需要网络。KEGG / Reactome / MSigDB 查询也需要网络。

## 分析步骤

1. 过滤低表达基因，再标准化（counts 用 DESeq2 size factor；FPKM/TPM 用 log2(x+1) + quantile）。
2. 三种比较：
   - `TG_sh1 vs NTC`
   - `TG_sh5 vs NTC`
   - `(TG_sh1 + TG_sh5)/2 vs NTC`（两个 knockdown 等权平均）
   - `common_up`：两次单独比较中的共同上调基因
3. 每个比较再按两套策略取上调基因：
   - FoldChange：FC ≥ 1、1.25、1.5、2
   - 上调排名：top 50 / 75 / 100 / 150 / 200 / 250 / 300
4. **每个比较 × 每个阈值都必须出图**，不能只出表格：
   - 差异基因表 + log2FC 柱状图
   - 火山图
   - 热图
   - GO 图（BP/MF/CC：barplot、dotplot、emapplot、cnetplot）
   - 通路富集图（Reactome + MSigDB Hallmark）
   - KEGG 图（barplot、dotplot、emapplot；并尽量输出 pathview 通路图）
   - GSEA 图（dotplot、ridgeplot、gseaplot2、Hallmark enrichment curve）
   默认 `padj < 0.05`；若无显著条目，会放宽阈值并在图标题标明 `relaxed cutoff`。

## 结果目录

每个比较下都有完整的 4 个 FC 文件夹和 7 个 topN 文件夹：

```
results/
  00_logs/
  TG_sh1_vs_NTC/
  TG_sh5_vs_NTC/
  TGsh_mean_vs_NTC/
  common_up/
    FoldChange/FC_1|FC_1.25|FC_1.5|FC_2/
      DE_selected_genes.xlsx
      DE_log2FC_barplot.pdf
      volcano.pdf
      heatmap.pdf
      GO/
      Pathway/
      KEGG/
      GSEA/
    TopRank/top50|top75|top100|top150|top200|top250|top300/
      （同上全套图）
```

后续修改本仓库中的 R 分析时，请遵循 `.cursor/rules/tg-brca-rnaseq.mdc`。
