# TCGA-BRCA 输入格式与样本匹配

## 表达矩阵（star_fpkm）

典型 UCSC Xena / GDC 长表：

```
Ensembl_ID          TCGA-3C-AAAU-01A    TCGA-3C-AALJ-01A
ENSG00000000003.15  12.3                8.1
```

- 第一列为基因 ID（Ensembl 或 Symbol）。
- 其余列为样本。
- 分析前：去 Ensembl 版本号、`log2(x + 1)`、按基因去重（保留均值）。

## 临床（clinical）

行为样本或患者，列为表型。常见列名（大小写/下划线可变）：

- `sample` / `submitter_id` / `_PATIENT`
- `age_at_diagnosis` / `age_at_initial_pathologic_diagnosis`
- `ajcc_pathologic_stage` / `tumor_stage`
- ER / PR / HER2 status
- `gender` / `race` / `histological_type`

连续变量用 Spearman；二分类 Wilcoxon；多分类 Kruskal-Wallis。

## 生存（survival）

```
sample    OS    OS.time
```

`OS` 为 0/1 事件，`OS.time` 为随访时间。每个 GO 单独做：

- `coxph(Surv(OS.time, OS) ~ score)`
- 按该通路评分中位数分组的 KM 曲线

## 蛋白（protein）

RPPA：蛋白为行、样本为列。每个 GO 的通路评分与每个蛋白 Spearman 相关；负相关蛋白单独输出。

## 注释（gencode.v36.annotation.gtf.gene）

至少包含 `gene_id` 与 `gene_name`。用于把表达矩阵的 ENSG 映射为 SYMBOL，并标注负相关基因的 `gene_type`。

## 样本 ID 对齐

1. `.` 换成 `-`，转大写。
2. 表达/蛋白用 15 位 aliquot（`TCGA-XX-XXXX-01A`）。
3. 临床/生存若只有患者 ID，用 12 位（`TCGA-XX-XXXX`）匹配。
4. 默认丢弃非 `01` 原发瘤；同一患者多份原发瘤时保留一份。

## 通路评分

- **ssGSEA**（GSVA）：每个 gene set 单独作为一个通路；一次调用可传入 named list，但 list 的每个元素仍是单个 GO，禁止 `Reduce(union, gene_sets)`。
- **回退**：对该 GO 在表达矩阵中的基因做样本内 z-score 后取均值。
- 表达矩阵中该 GO 基因数 `< min_genes`（默认 5）则跳过该通路。

## 负相关基因

对每个 GO、每个基因：

```
spearman(pathway_score, gene_expression)
```

保留 `rho < 0` 且 `p.adjust(p, "BH") < fdr_cutoff`（默认 0.05）。FDR 只在该 GO 的全基因组相关结果内计算。`in_geneset` 列标记该基因是否属于当前 GO。
