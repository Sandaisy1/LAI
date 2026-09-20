# 乳腺癌转移小鼠模型：原位 + 多器官组学文献与公开数据

检索范围：小鼠（或人源细胞种到小鼠）乳腺癌，同时或分别覆盖**原位乳腺/脂肪垫肿瘤**与**转移灶**（肺、骨、肝、脑、淋巴结）。数据类型：bulk / 单细胞 / 空间 RNA，以及蛋白组（组织、分泌组、细胞外囊泡）。

**结论先说：** 目前没有一份公开数据同时包含「原位组织 + 肺 + 骨 + 肝 + 脑」的配对 bulk RNA-seq。要用多器官比较，需要把下面几套数据拼起来，或优先用最接近的多器官 scRNA-seq。

与本仓库 TG BRCA 细胞 RNA-seq 的对接方式：把本实验上调基因，和下列「原位 vs 转移」或「器官嗜性」签名做 overlap / GSEA，而不是把小鼠组织直接并进现有 Cuffdiff 流程。

---

## 怎么选数据

| 目的 | 首选 | 注意 |
| --- | --- | --- |
| 一次看脑、肺、肝、骨微环境 | E-MTAB-16621（VO-PyMT 心内注射 scRNA-seq） | 实验性转移，**没有乳腺原位瘤** |
| 原位 + 多个远端器官（肺、淋巴结、骨髓） | GSE165393（MMTV-PyMT bulk + scRNA-seq） | **没有肝、脑组织** |
| 原位 vs 肺（RNA，4T1） | GSE146012（分选肿瘤细胞 bulk RNA-seq） | 样品数少（n=8） |
| 原位 vs 骨（经典、配对） | GSE37975（4T1.2 芯片） | 不是 RNA-seq |
| 原位 vs 肝（组织切片） | GSE238214（LCM 芯片 + 代谢） | 脾内注射成肝转移，不是完全自发 |
| 原位 vs 肺 vs 脑（同一套 4T1 衍生系） | GSE54773（芯片） | 测的是**体外再培养的细胞系**，不是新鲜组织 |
| 肺 + 肝转移龛（肿瘤细胞 + 龛细胞） | GSE252507（scRNA-seq） | 预印本；原位覆盖有限 |
| 四个器官嗜性的**蛋白** | PXD055261（4T1 器官嗜性系 EV 蛋白组） | 培养上清 EV，不是转移组织本身 |
| BRCA1 小鼠 | GSE148614 / GSE122075 | 几乎都是**原发瘤**，不是转移灶 |

注射方式会改变生物学含义：

- 脂肪垫原位（orthotopic / MFP）：最接近「原位 → 转移」
- 心内注射（intracardiac）：跳过原发瘤，直接多种子器官
- 脾内 / 尾静脉：分别偏向肝、肺定植

---

## 1. 优先推荐（最接近「原位 + 多器官」）

### 1.1 原位 + 肺 + 淋巴结 + 骨髓（RNA-seq）

| 项目 | 内容 |
| --- | --- |
| 文章 | Transcriptome analysis of heterogeneity in mouse model of metastatic breast cancer |
| 期刊 | *Breast Cancer Research* 2021 |
| DOI | https://doi.org/10.1186/s13058-021-01468-x |
| 数据 | **GSE165393**（GEO） |
| 模型 | MMTV-PyMT，自发转移 |
| 组织 | 乳腺脂肪垫原发瘤、肺、淋巴结、骨髓来源的肿瘤细胞系 |
| 类型 | bulk RNA-seq（NextSeq 500）+ 部分样本 scRNA-seq（ddSeq） |
| 样品规模 | 17 个 GSM：14 个 bulk + 3 个 scRNA-seq |
| 用途 | 本仓库上调基因 vs PyMT 器官来源转录组；看器官嗜性和 CD44 高低亚群 |

#### GSE165393 该下哪些

GEO 只给了处理后的表，没有 count 矩阵。和本仓库 TG 上调基因做 overlap / 器官比较时：

**先下这两个（一共不到 1 MB）：**

1. `GSE165393_AllTissues_TPM.csv.gz` — bulk TPM 矩阵（必下）
2. `GSE165393_series_matrix.txt.gz` — 样品注释（很小；表达表是空的，不要当表达矩阵用）

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165393/suppl/GSE165393_AllTissues_TPM.csv.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165393/matrix/GSE165393_series_matrix.txt.gz
```

TPM 表有 55450 个 `ENSMUSG` 基因 × **14 个 bulk 样品**（不要把 3 个 scRNA GSM 并进去）：

| 组织 | 样品列 | 生物学含义 |
| --- | --- | --- |
| 原位（乳腺脂肪垫） | `MFP1.TPM`, `MFP2.TPM` | 原发瘤来源，n=2，不分 CD44 |
| 肺转移 | `LUlow1/2`, `LUhigh1/2` | CD44low / CD44high 各 2 |
| 淋巴结 | `LN1low/2low`, `LN1high/2high` | 同上 |
| 骨髓 | `BM1low/2low`, `BM1high/2high` | 同上；不是皮质骨转移灶 |

比较建议：先做「器官主效应」（MFP vs 肺 / LN / BM，CD44 高低先合并或作为协变量），再单独看 CD44high vs low。GEO **没有 raw count**，不要拿 TPM 去跑 DESeq2；用 `log2(TPM+1)` + limma。基因 ID 是小鼠 Ensembl，和人 TG 结果比较前要做同源转换。

**默认不要下：**

- 三个 scRNA 矩阵（`LungH_rep1` / `LungL_rep1` / `LymphNodeH`，各约 2800–3600 细胞、raw counts）— 只有肺和淋巴结、没有原位和骨髓，适合看亚群，不适合先做器官 FC
- SRA FASTQ（SRP302977）— 只有需要自己重定量出 count 时才下
- Series matrix 里的表达表 — `data_row_count = 0`，是空的

注意：这些是器官取出后 FACS（CD44/EpCAM）再培养的**细胞系**，不是新鲜肿瘤块。没有肝、脑。

### 1.2 脑 + 肺 + 肝 + 骨同一实验（scRNA-seq，无原位瘤）

| 项目 | 内容 |
| --- | --- |
| 文章 | Single-cell profiling of synchronous multi-organ metastasis reveals a systemic CD74+ lipid-associated macrophage niche driving polymetastatic breast cancer |
| 来源 | bioRxiv 2026（预印本） |
| DOI | https://doi.org/10.64898/2026.01.31.701004 |
| 数据 | **E-MTAB-16621**（ArrayExpress / ENA ERP188936） |
| 模型 | VO-PyMT，左心室注射，同步多器官转移 |
| 组织 | 脑、肺、肝、骨；FACS 分成肿瘤细胞（GFP+）、近端龛（mCherry+）、远端间质 |
| 类型 | 10x Chromium scRNA-seq，约 73,051 细胞 / 11 个组织 |
| 用途 | 四器官转移龛比较的目前最完整小鼠图谱；不能回答「原位 vs 转移」 |

### 1.3 肺 vs 肝转移龛（scRNA-seq）

| 项目 | 内容 |
| --- | --- |
| 文章 | Charting the liver and lung metastatic niche in breast cancer |
| 来源 | bioRxiv 2024（预印本） |
| DOI | https://doi.org/10.1101/2024.10.04.616491 |
| 数据 | **GSE252507** |
| 模型 | TNBC 小鼠（含 4T1 / MVT1 等），自发内脏转移 + 条形码示踪 + Cherry-niche |
| 组织 | 肺转移、肝转移（肿瘤细胞 + 邻近龛细胞） |
| 类型 | scRNA-seq；部分实验有 hashing |
| 机制线索 | 肝内皮 BMP2 促进再播散；肺巨噬细胞 Granulin 抑制转移表型 |

---

## 2. 按转移器官：原位 vs 转移

### 2.1 肺

| 数据 | 文章 | 模型 | 材料 | 类型 | 备注 |
| --- | --- | --- | --- | --- | --- |
| **GSE146012** | So et al. *Cancer Research* 2020. Induction of DNMT3B by PGE2 and IL6 at distant metastatic sites. | 4T1 / BALB/c | GFP 分选：原发瘤 vs 肺转移（含 MFP 与尾静脉两条路径） | bulk mRNA-seq，8 个样本 | 配套 **GSE146010** 为肺转移 DNMT3B ChIP-seq |
| **GSE273439** | Bertolazzi et al. *J Exp Clin Cancer Res* 2025. Extraction of a stromal metastatic gene signature. DOI: 10.1186/s13046-025-03353-3. PMID: 40065328 | 4T1 脂肪垫接种 28 天 | 原发乳腺瘤 + 肺，FFPE | Visium 空间转录组 | 原位与肺在同一实验；**该下哪些见 §4.1** |
| **GSE131508** | Ombrato et al. *Nature* 2019. Metastatic niche labelling reveals tissue parenchyma stem cell features. DOI: 10.1038/s41586-019-1487-6 | Labelling-4T1 | 肺转移龛（mCherry+）vs 远端肺（mCherry−） | bulk / 后续 scRNA | 经典 Cherry-niche；偏龛细胞不是肿瘤细胞 |
| **GSE318532** | Targeting FASN/GPAM in AT2 cells decreases lung metastasis | 4T1 脂肪垫 | 转移肺 vs 对照肺 | Visium | 看肺泡 II 型细胞脂质支持，不是分选肿瘤细胞 |
| PXD005860 | *J Cancer* 2017. Proteomic analysis of lung metastases… CTSB/CTSL | MMTV-PyMT 肺转移灶 | 肺转移组织 | LC-MS/MS | 蛋白组；比较转基因 vs WT，不是原位 vs 肺 |
| — | *PLOS One* 2015. Differential proteome… TGF-β in 4T1. DOI: 10.1371/journal.pone.0126483 | 4T1 肺转移 | 肺转移组织 ± TGF-β 抑制剂 | Orbitrap 定量蛋白组 | 6694 蛋白；治疗对照，不是原位配对 |

肺空间多组学补充：*Cell Death & Disease* 2024, Metabolic shifts in lipid utilization… DOI: 10.1038/s41419-024-07205-4（4T1 / PyMT 肺转移 Visium + 蛋白相关分析）。

### 2.2 骨

| 数据 | 文章 | 模型 | 材料 | 类型 | 备注 |
| --- | --- | --- | --- | --- | --- |
| **GSE37975** | Bidwell et al. *Nature Medicine* 2012. Silencing of Irf7 pathways… bone metastasis. DOI: 10.1038/nm.2830 | 4T1.2 自发骨转移 | **同鼠配对**原发瘤 vs 骨转移，4 对 | Affymetrix 430 2.0 芯片 | 骨转移研究最常引用的配对转录组；Irf7 / I 型干扰素 |
| **GSE241165** | Lymphotoxin-β promotes bone colonization… bioRxiv 2023. DOI: 10.1101/2023.08.15.553179 | 4T1 vs 4T1.2 心内注射 | 体外、脂肪垫原发、骨 D4/D10/D16 分选肿瘤细胞 | SMART-seq2 单细胞（1248 个样本记录） | 覆盖原位对照 + 骨定植时程 |
| **GSE160101 / GSE160102 / GSE160100** | 骨髓内皮重塑相关 RNA-seq（4T1 家族） | 4T1.2 心内注射 | 细胞系 67NR / 66cl4 / 4T1 / 4T1.2；转移 vs 邻近 vs 正常骨髓内皮 | bulk RNA-seq | 看骨微环境内皮，不是肿瘤细胞全组织 |

骨转移分泌组（蛋白，细胞系）：Blanco et al. Global secretome analysis identifies novel mediators of bone metastasis. PMC3434351（MDA-MB-231、4T1/4T1.2、膀胱癌系 SILAC 分泌组）。

### 2.3 肝

| 数据 | 文章 | 模型 | 材料 | 类型 | 备注 |
| --- | --- | --- | --- | --- | --- |
| **GSE238214** | Metastatic breast cancer cells are metabolically reprogrammed to maintain redox homeostasis. PMC11321393 | 4T1-2776 / 4T1-2792 亲肝系 | 脂肪垫原发瘤（LCM core/margin）vs 脾内肝转移（core/margin/adjacent/distant，多时间点） | 芯片，64 个样本 | 肝转移 GSH / ROS 解毒；有代谢示踪 |
| **GSE252507** | 见 1.3 | 自发内脏转移 | 肝转移肿瘤细胞 + 龛 | scRNA-seq | 与肺对照最合适 |

4T1 亲肝衍生系来源：Tabariès / Siegel 实验室对 4T1-2776、4T1-2792 的体内筛选（后续被 GSE238214 与 EV 蛋白组沿用）。

### 2.4 脑

| 数据 | 文章 | 模型 | 材料 | 类型 | 备注 |
| --- | --- | --- | --- | --- | --- |
| **GSE54773** | *Science Translational Medicine* 2020. Connexins orchestrate progression of breast cancer metastasis to the brain. DOI: 10.1126/scitranslmed.aax8933 | 4T1 原位筛选：4T1-T2（原发）、4T1-LM2（肺）、4T1-BM2（脑） | 各 3 只鼠来源细胞系，体外平行培养后抽 RNA | 芯片，9 个样本 | **同一亲本、三个解剖部位**；测的是细胞系不是新鲜脑组织 |
| E-MTAB-16621 | 见 1.2 | VO-PyMT 心内 | 脑转移灶 + 近/远端间质 | scRNA-seq | 新鲜组织单细胞，无配对原发瘤 |

脑转移 TMT 蛋白组（机制文，未必有完整组织库）：IL6/CCL2 from M2-polarized microglia… *Front Pharmacol* 2025. DOI: 10.3389/fphar.2025.1547333。

### 2.5 淋巴结

| 数据 | 文章 | 模型 | 类型 |
| --- | --- | --- | --- |
| **GSE165393** | 见 1.1 | PyMT 淋巴结来源细胞 | bulk + scRNA |
| **GSE168181** | Single-cell sequencing reveals cancer cell heterogeneity in a murine breast cancer lymph node metastasis model | 自发淋巴结转移：原发瘤 vs 引流淋巴结 | scRNA-seq，6 个样本 |

---

## 3. 蛋白组（肺 / 骨 / 肝 / 脑都能碰到的）

组织蛋白组几乎都是**单器官转移灶**，没有四器官配对组织 MS。四器官比较目前最完整的是**器官嗜性细胞系的 EV 蛋白组**。

| 数据 | 文章 | 材料 | 覆盖器官 | 类型 |
| --- | --- | --- | --- | --- |
| **PXD055261** | Characterization of extracellular vesicle-associated DNA and proteins derived from organotropic metastatic breast cancer cells. *J Exp Clin Cancer Res* 2025. DOI: 10.1186/s13046-025-03418-3 | NMuMG、67NR、4T1 亲本，以及肺（4T1-533/537）、骨（592/593）、肝（2776/2792）、脑（BP/LM）衍生系的小 EV | 肺、骨、肝、脑 | LC-MS/MS，约 698 蛋白 |
| PXD005860 | PyMT ± CTSB/CTSL 肺转移 | 肺转移组织 | 肺 | 二甲基标记定量 |
| — | *PLOS One* 2015 TGF-β / 4T1 | 肺转移组织 | 肺 | Orbitrap |
| — | Blanco et al. 骨转移分泌组 | 4T1 vs 4T1.2 等培养上清 | 骨 | SILAC |
| — | *PLOS One* 2011 MDSC proteome, 67NR vs 4T1 | 荷瘤小鼠 MDSC，不是肿瘤组织 | 间接反映肺转移负荷 | shotgun MS |

PXD055261 的细胞系与 GSE238214（肝）、GSE54773（脑/肺）同属 4T1 体内筛选家族，适合做 **RNA（芯片/RNA-seq）+ EV 蛋白** 的跨组学对照。

---

## 4. 空间转录组（原位和/或转移灶）

| 数据 | 组织 | 平台 | 模型 |
| --- | --- | --- | --- |
| **GSE273439** | 原发乳腺瘤 + 肺 | Visium FFPE | 4T1 |
| **GSE318532** | 转移肺 vs 正常肺 | Visium 新鲜冰冻 | 4T1 |
| **GSE300613** | 原发瘤 Visium；配套 scRNA-seq 覆盖淋巴结、肝、肺转移 | Visium + scRNA | **MDA-MB-231 人源异种移植**（不是纯小鼠瘤） |

GSE300613 是人源细胞在小鼠体内，基因表达是人/鼠混合；用来验证人 BRCA 细胞签名比纯小鼠瘤更直接，但免疫微环境是 NSG/免疫缺陷背景。

### 4.1 GSE273439 该下哪些

这不是 bulk TPM，而是 **10x Visium FFPE**（Space Ranger → mm10）。两只 BALB/c（Mouse 627 / 628），各一块原发乳腺瘤（TM）和一块肺（LUNG），共 4 个 visium 切片。文章里的 **Nanostring DSP**（人原发瘤 + 肝转移 24 个 ROI、stromal Met 签名）**不在 GEO**，在期刊 Supplementary Table 1–7。

配对分析脚本：`GSE273439_Mouse_breast.R`（627 原位对 627 肺，628 对 628；本套没有骨）。

**必下这一个包（约 295 MB），四个样品全在里面：**

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273439/suppl/GSE273439_RAW.tar
```

可选、很小：`GSE273439_series_matrix.txt.gz`（3 KB）只当样品注释用。`data_row_count = 0`，**不要当表达矩阵**。

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273439/matrix/GSE273439_series_matrix.txt.gz
```

解压后每个 GSM 有 9 个文件。Seurat `Load10X_Spatial()` / SpatialExperiment 真正要用的是 MTX + 坐标 + 缩放因子 + PNG；大 TIFF 可以后下。

| GSM | 标题 | 组织 | 分析必用 | 可后下 |
| --- | --- | --- | --- | --- |
| GSM8428408 | Mouse_627_TM | 原发乳腺瘤 | `matrix.mtx.gz`（21 MB）、`features.tsv.gz`、`barcodes.tsv.gz`、`tissue_positions_list.csv.gz`、`scalefactors_json.json.gz`、`tissue_hires/lowres_image.png.gz` | `aligned_fiducials.jpg.gz`；`.tif.gz` 11 MB |
| GSM8428409 | Mouse_627_LUNG | 肺 | 同上（mtx 18 MB） | `.tif.gz` 47 MB |
| GSM8428410 | Mouse_628_TM | 原发乳腺瘤 | 同上（mtx 17 MB） | `.tif.gz` 61 MB |
| GSM8428411 | Mouse_628_LUNG | 肺 | 同上（mtx 39 MB） | `.tif.gz` 72 MB |

文件名带前缀，例如 `GSM8428408_09.9_627_TM_V11Y03-083_A_matrix.mtx.gz`。读进 Seurat 前按样品建目录并改成 10x 标准名：

```
Mouse_627_TM/
  filtered_feature_bc_matrix/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}
  spatial/{tissue_positions_list.csv,scalefactors_json.json,tissue_hires_image.png,tissue_lowres_image.png}
```

`spatial/` 里的 png / json / csv 需要 **gunzip**（Seurat 不认 `.png.gz`）。GEO 说明里写了 `detected_tissue_image.jpg.gz`，实际 **没有这个文件**。

**默认不要下：**

- SRA FASTQ（PRJNA1141928 / SRX25515152–155，NextSeq 2000）— 只有要自己重跑 Space Ranger 时才下
- 四张全分辨率 `.tif.gz`（合计约 190 MB）— 热图/聚类用 PNG 就够；重新配准或重跑 spaceranger 再下
- Series matrix 里的表达表 — 空的
- 期刊 DSP / stromal Met 基因表 — 去论文补充材料，不在 `RAW.tar`

注意：只有原位和肺，没有骨、肝、脑。小鼠 Visium 是验证集；签名本身来自人 DSP。和 TG 上调基因做 overlap 时，先把小鼠基因转人同源，再和 Supplementary Table 2（129 up / 99 down）比，不要把 visium spot count 直接并进 Cuffdiff。

---

## 5. 和本仓库 BRCA 背景相关、但不是转移灶

本仓库分析的是人 BRCA 细胞 TG 敲低 RNA-seq。下列 BRCA1 小鼠数据几乎只有**原发乳腺瘤**，不能当转移组织对照，只适合看 BRCA1 缺陷原发转录组。

| 数据 | 模型 | 类型 | 转移组织？ |
| --- | --- | --- | --- |
| GSE148614 / GSE148565 | MMTV-Cre; Brca1fl/fl ± Trp53 | bulk + scRNA-seq，23 个原发瘤 | 否 |
| GSE122075 | K14-Cre; p53F/F Brca1F/F | 芯片，原发瘤 | 否 |

---

## 6. 建议的最小下载清单

若下一步要和 TG 上调基因做 overlap，建议按这个顺序取矩阵（不必一次下完原始 FASTQ）：

1. **GSE165393** 处理后的 counts / RSEM（原位 + 肺 + LN + 骨髓）
2. **GSE146012** 补充文件 `GSE146012_RNA_seq_logtransformed_count.txt.gz`（4T1 原位 vs 肺）
3. **GSE37975** Series Matrix（4T1.2 原位 vs 骨，芯片）
4. **GSE238214** Series Matrix（原位 vs 肝，芯片）
5. **GSE54773** Series Matrix（原发 / 肺 / 脑衍生系，芯片）
6. **PXD055261** MaxQuant/蛋白表（四器官 EV 蛋白）
7. 若要单细胞多器官：E-MTAB-16621 或 GSE252507 的 h5ad / 10x 矩阵
8. 若要原位 vs 肺的空间转录组：**GSE273439** 只下 `GSE273439_RAW.tar`（Visium MTX + 坐标 + PNG；见 §4.1）

芯片与 RNA-seq 不要直接合并 counts；只在基因符号层做签名 overlap。小鼠基因映射到人用 `org.Mm.eg.db` / `org.Hs.eg.db` 同源转换后再和本仓库结果比较。

---

## 7. 常用模型速查（4T1 家族）

同一株 BALB/c 乳腺瘤拆出的同基因系列（Aslakson & Miller 1992）：

| 细胞系 | 转移能力 |
| --- | --- |
| 67NR | 不离开原发部位 |
| 168FARN | 淋巴结 |
| 4T07 | 微转移，肉眼灶少 |
| 66cl4 | 肺 |
| 4T1 | 肺、肝、骨、脑（自发，免疫健全） |
| 4T1.2 | 4T1 亚克隆，骨转移强、肺相对弱 |

转基因自发：MMTV-PyMT（肺为主，骨髓/LN 可见）、C3(1)-Tag、MMTV-Neu。实验性多器官：心内注射 VO-PyMT 或 4T1。

---

## 8. 检索日期与缺口

检索日期：2026-09-20。主要来源：GEO、ArrayExpress、PRIDE、PubMed/PMC。

仍缺的公开资源：

- 同一只鼠、同一平台的「乳腺原位 + 肺 + 骨 + 肝 + 脑」bulk RNA-seq
- 四器官**转移组织**（不是 EV / 不是细胞系）的配对蛋白组
- BRCA1 条件敲除小鼠的配对转移灶转录组

机器可读表见同目录 `mouse_breast_cancer_metastasis_datasets.csv`。
