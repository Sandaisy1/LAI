# LAI

TCGA-BRCA 神经相关 GO 通路分析：**每一个 GO 通路单独**与临床相关，并 **对每一个 GO 通路单独** 寻找负相关基因（禁止把 17 条通路合并成一条签名）。

不依赖 `GSVA`、`AnnotationDbi`、`org.Hs.eg.db`。通路活性为 combined z-score；GO 基因来自 QuickGO 或 NCBI `gene2go`。

## 输入数据

将下列文件放到 `data/`（文件名含这些前缀即可）：

1. `TCGA-BRCA.clinical`
2. `TCGA-BRCA.protein`
3. `TCGA-BRCA.star_fpkm`
4. `TCGA-BRCA.survival`
5. `gencode.v36.annotation.gtf.gene`

可选：`data/go_genesets.tsv`（列 `GO_ID`, `symbol`）或 `.gmt`，则不再联网取 GO 基因。

## 运行

```bash
Rscript scripts/install_r_packages.R
Rscript scripts/go_pathway_per_term_analysis.R data results
```

只需 CRAN 包：`data.table`、`survival`。

## 输出

- `results/pathway_zscores_per_GO.tsv`：每个 GO 一行的样本活性
- `results/clinical/<GO_ID>_clinical.tsv`：该 GO 与临床/生存的相关
- `results/negcorr/<GO_ID>_negative_genes.tsv`：与该 GO 活性负相关的基因（rho < 0 且 FDR < 0.05）
