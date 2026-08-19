# TG BRCA 细胞 RNA-seq 分析

针对 NTC、TG_sh1、TG_sh5 三组人乳腺癌细胞 RNA-seq 的差异表达与富集分析流程。

## 数据位置

默认读取 `E:/R/TG_BRCA/TG`。若该路径不存在，则使用当前工作目录。也可在运行前设置：

```r
Sys.setenv(TG_RNASEQ_DIR = "E:/R/TG_BRCA/TG")
```

优先读取表达矩阵 `shTG(20220723).xlsx`。若 Excel 列名无法识别分组，会回退到同目录下 Cuffdiff 文件：

- `genes.read_group_tracking`
- `genes.count_tracking`
- `genes.fpkm_tracking`

Excel 中样本列名需能对应到 `NTC`、`TG_sh1`、`TG_sh5`（例如 `NTC_1`、`TG_sh1-2`、`shTG5` 均可）。若缺少 TG_sh5，脚本会跳过 sh5 相关比较并写入日志。

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
4. 每个非空基因集输出：差异表、火山图、热图、GO（BP/MF/CC）、Reactome、MSigDB Hallmark、KEGG、GSEA。默认 `padj < 0.05`；重复不足无法估计 P 值时改为仅按 FC/排名筛选。

## 结果目录

```
results/
  00_logs/
  00_QC_PCA.*
  normalized_log_matrix.csv
  TG_sh1_vs_NTC/
  TG_sh5_vs_NTC/
  TGsh_mean_vs_NTC/
  common_up/
    FoldChange/FC_1 ... FC_2/
    TopRank/top50 ... top300/
```

后续修改本仓库中的 R 分析时，请遵循 `.cursor/rules/tg-brca-rnaseq.mdc`。
