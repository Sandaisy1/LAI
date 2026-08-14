# LAI

TCGA-BRCA 神经相关 GO 通路分析：**每一个 GO 通路单独**与临床相关，并 **对每一个 GO 通路单独** 寻找负相关基因（禁止把 17 条通路合并成一条签名）。

## 输入数据

将下列文件放到 `data/`（文件名含这些前缀即可）：

1. `TCGA-BRCA.clinical`
2. `TCGA-BRCA.protein`
3. `TCGA-BRCA.star_fpkm`
4. `TCGA-BRCA.survival`
5. `gencode.v36.annotation.gtf.gene`

## 运行

```bash
Rscript scripts/install_r_packages.R
Rscript scripts/go_pathway_per_term_analysis.R data results
```

## 输出

- `results/gsva_scores_per_GO.tsv`：每个 GO 一行的样本活性
- `results/clinical/<GO_ID>_clinical.tsv`：该 GO 与临床/生存的相关（每个通路一张表）
- `results/negcorr/<GO_ID>_negative_genes.tsv`：与该 GO 活性负相关的基因（rho < 0 且 FDR < 0.05）

Cursor 规则见 `.cursor/rules/r-go-pathway-tcga-brca.mdc`，编写 R 代码时会强制“每个 GO，不是所有 GO 合并”。
