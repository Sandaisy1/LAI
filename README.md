# LAI

TCGA-BRCA 中 **每一个** 神经元/轴突相关 GO 通路单独分析：与临床的关联，以及与该通路得分呈负相关的基因。所列 GO **不会**被合并成一个基因集或一个综合得分。

## 你需要准备的数据

把下列文件放在同一目录（扩展名可以是 `.tsv` / `.txt` / `.csv`，可 `.gz`），文件名包含这些片段即可：

| 文件 | 用途 |
|------|------|
| `TCGA-BRCA.clinical` | 临床表型 |
| `TCGA-BRCA.protein` | RPPA 蛋白（可选；缺失则跳过蛋白负相关） |
| `TCGA-BRCA.star_fpkm` | STAR FPKM 表达 |
| `TCGA-BRCA.survival` | 生存（OS/DSS/PFI 等） |
| `gencode.v36.annotation.gtf.gene` | Ensembl ID ↔ gene symbol / biotype |

## 分析的 17 条 GO（逐条独立）

通路列表在 `R/00_config.R` 的 `GO_TERMS`。亲代/子代（例如 `GO:0007411` 轴突导向 与 `GO:0097374` 感觉神经元轴突导向）也会分开计算。

对 **每一个** GO ID，脚本会：

1. 只用该 GO 的基因（默认 `GOALL` 含子术语注释）计算每个病人的通路活性（默认：成员基因跨样本 z-score 再取均值；可改 `--method ssgsea`）
2. 将该通路得分与临床变量做 Spearman / Wilcoxon / Kruskal-Wallis，并与生存做 Cox 与 KM
3. 将全基因组表达（及蛋白）与 **该通路得分** 做 Spearman，筛 `rho ≤ -0.2` 且 FDR ≤ 0.05 的负相关特征

## 运行

```bash
# 依赖（真实数据需要 Bioconductor：org.Hs.eg.db）
Rscript R/run_go_tcga_brca.R --install-deps

# 正式分析
Rscript R/run_go_tcga_brca.R --data-dir /path/to/your/data --out-dir ./results

# 无 org.Hs.eg.db 时，可用 MSigDB 风格 GMT（第一列必须是 GO:xxxxxxx）
Rscript R/run_go_tcga_brca.R --data-dir /path/to/your/data --gmt /path/to/c5.go.v2023.Hs.symbols.gmt

# 不读真实数据，验证“每个 GO 独立打分”
Rscript R/run_go_tcga_brca.R --demo --out-dir ./results_demo
```

## 输出

`results/` 下 **每个 GO 一个文件夹**（如 `GO_0007411/`）：

- `go_member_genes.tsv` — 该通路映射到表达矩阵的基因
- `pathway_score.tsv` — 每个样本的该通路得分
- `clinical_association.tsv` — 该通路 vs 临床
- `survival_cox.tsv` / `os_km.pdf` — 该通路 vs 生存
- `negative_correlated_genes.tsv` — 与该通路负相关的基因
- `negative_correlated_proteins.tsv` — 与该通路负相关的蛋白（若有蛋白数据）

根目录的 `summary_*_long.tsv` 只是各 GO 结果的长表拼接，**不是**把通路合并后再统计。

Cursor 规则见 `.cursor/rules/r-go-per-term-tcga.mdc`：后续改 R 代码时会约束模型不要把这些 GO 合并分析。
