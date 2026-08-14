---
name: tcga-brca-go-per-pathway
description: Analyzes each listed GO pathway independently in TCGA-BRCA (never pooling terms). Computes per-GO clinical correlation and per-GO negatively correlated genes from star_fpkm, clinical, survival, protein, and gencode annotation. Use when working with TCGA-BRCA, GO:0023041/GO:0007411-style IDs, 每个通路, 负相关基因, ssGSEA, axon guidance, or axonogenesis.
---

# TCGA-BRCA 逐个 GO 通路分析

## 硬性约束（必须遵守）

1. **每个 GO 通路单独分析，禁止合并。** 禁止把所列 GO 基因并成一个大基因集，禁止对通路评分取平均后再做相关，禁止把所有通路的相关 p 值放到一起做一次 FDR。
2. 临床相关与负相关基因筛选都必须在 **单个 GO ID** 的循环内完成。
3. 用户给出的 GO ID 与注释必须 **逐字保留**，不得改写、删减或替换。

当前通路列表（verbatim）：

```r
go_ids <- c(
  "GO:0023041",  # neuronal signal transduction
  "GO:1904457",  # positive regulation of neuronal action potential
  "GO:1904340",  # positive regulation of dopaminergic neuron differentiation
  "GO:2001224",  # positive regulation of neuron migration
  "GO:2001222",  # regulation of neuron migration
  "GO:0019227",  # neuronal action potential propagation
  "GO:0019228",  # neuronal action potential
  "GO:1902847",  # regulation of neuronal signal transduction
  "GO:0031102",  # neuron projection regeneration
  "GO:0097492",  # sympathetic neuron axon guidance
  "GO:0097491",  # sympathetic neuron projection guidance
  "GO:0097374",  # sensory neuron axon guidance
  "GO:0007158",  # neuron cell-cell adhesion
  "GO:1902667",  # regulation of axon guidance
  "GO:0031103",  # axon regeneration
  "GO:0007411",  # axon guidance
  "GO:0007409"   # axonogenesis
)
```

## 数据文件

默认放在 `data/`，按前缀匹配（允许 `.tsv` / `.txt` / `.gz`）：

| 前缀 | 用途 |
|------|------|
| `TCGA-BRCA.clinical` | 临床表型 |
| `TCGA-BRCA.protein` | RPPA 蛋白 |
| `TCGA-BRCA.star_fpkm` | STAR FPKM 表达矩阵 |
| `TCGA-BRCA.survival` | 生存（OS / OS.time） |
| `gencode.v36.annotation.gtf.gene` | ENSG ↔ SYMBOL |

表达矩阵：基因为行、样本为列；Ensembl ID 去掉版本号（`ENSGxxxx.15` → `ENSGxxxx`）。默认只保留原发瘤（barcode 第 14–15 位为 `01`）。FPKM 使用 `log2(FPKM + 1)`。

## 分析流程

入口脚本：`R/analyze_go_per_pathway.R`（配置在 `R/config.R`）。

对 **每一个** GO ID：

1. 用 `org.Hs.eg.db` 的 `GOALL` 取该 term 及其子 term 基因；**不要**与其他 GO 取并集。
2. 仅用该通路基因计算样本通路评分（优先 ssGSEA，否则通路基因 z-score 均值）。
3. 通路评分 vs 临床变量（连续 Spearman；分类 Wilcoxon / Kruskal-Wallis；生存 Cox / KM）。
4. 通路评分 vs 全基因组表达，Spearman；保留 **rho < 0** 且通过该通路内 BH-FDR 的基因。
5. 结果写到 `results/<GO_ID>/`，互不覆盖。

基因相关的多重检验只在 **该通路内部** 做 BH。蛋白数据若存在，同样按每个 GO 单独相关。

## 输出

每个 `results/<GO_ID>/` 至少包含：

- `geneset.tsv` 通路基因
- `pathway_scores.tsv` 样本评分
- `clinical_correlation.tsv` 临床相关
- `survival_cox.tsv` 与 KM 图
- `negative_gene_correlation.tsv` 负相关基因
- 可选 `protein_correlation.tsv`

汇总表（仍按行分通路，不是合并分析）：`results/summary_clinical_correlation.tsv`、`results/summary_negative_gene_counts.tsv`。

## 修改代码时

- 改 GO 列表只动 `R/config.R`。
- 新增统计必须放进 `for (go_id in go_ids)`，以 `go_id` 为 key 存结果。
- 某个 GO 映射到 0 个表达基因时记录警告并跳过，不要用其他通路填补。

数据格式细节见 [reference.md](reference.md)。
