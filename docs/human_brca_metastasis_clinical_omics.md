# 人源乳腺癌转移：临床组织 RNA / 蛋白组文献目录

筛选目标：**临床样本为主**，组织同时覆盖**人源原发（原位）组织**和**人源转移组织**（优先肺、骨、肝、脑），分子层为 **RNA 或蛋白**。细胞系、小鼠、PDX 不当病人；仅原发灶加随访结局、没有转移灶组织的队列单独列出，不当作「转移组织」。

详细字段见 `docs/human_brca_metastasis_clinical_omics.csv`。

## 怎么选数据

| 研究问题 | 优先队列 |
| --- | --- |
| 同一患者原发 vs 远处转移（全转录组，含肺/肝/脑/骨） | AURORA US + UNC RAP（GSE209998 / GSE193103 / GSE110590） |
| 最大配对原发–转移 RNA-seq | AURORA EU（152 对，需申请） |
| 脑转移配对 RNA | Cosgrove / Varešlija（GSE184869 / GSE173661） |
| 肺 / 胸膜配对 RNA | Sinn 等（GSE145752，NanoString 269 基因） |
| 多器官转移灶部位比较（转移灶活检为主） | Brasó-Maristany（GSE175692）；Zhang（GSE14020） |
| 原发 vs 淋巴结蛋白组 | Pozniak / Geiger（PXD000815） |
| 单细胞 / 空间、临床活检 | Klughammer HTAN（Nat Med 2024） |

**不要**把下列数据当成「人源转移组织」：`GSE2603`、`GSE5327`（原发灶 + 肺转移**临床结局**，没有转移灶活检）、TCGA-BRCA / CPTAC-BRCA（几乎全是原发）、MDA-MB-231 衍生肺/骨/脑亚系。

---

## A. 首选：人源原发 + 远处转移组织（肺 / 骨 / 肝 / 脑）

### A1. AURORA US 多组学（目前最完整的公开 bulk RNA 资源）

- **Garcia-Recio et al.**, *Nature Cancer* (2023). DOI: [10.1038/s43018-022-00491-x](https://www.nature.com/articles/s43018-022-00491-x)
- **样本**：55 例患者；51 个原发 + 102 个转移灶。转移部位：肝 28、肺 13、淋巴结 12、脑 11，另有骨等共约 20 个部位。20 例为尸检，可在同一患者内比较多个转移器官。
- **数据**：rRNA depletion RNA-seq、WES、low-pass WGS、甲基化芯片。无蛋白组。
- **下载**：RNA-seq GEO [GSE209998](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE209998)（SuperSeries [GSE212375](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212375)）；原始数据 dbGaP [phs002622](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002622.v1.p1)；切片 TCIA `AURORA-Metastatic-Breast-Multiomics`。
- **器官比较**：该文把 AURORA 与 RAP、GEICAM RNA-seq 合并后比较肝 / 肺 / 脑 vs 原发；肝转移免疫签名偏低，脑转移免疫/基质偏低，肺 vs 原发差异最小。

### A2. UNC Rapid Autopsy（同一患者多个转移器官）

- **Siegel et al.**, *J Clin Invest* (2018). DOI: [10.1172/JCI96153](https://doi.org/10.1172/jci96153)
- **样本**：16 例患者，原发 + 67 个配对转移灶（每例 2–7 个），含肝、肺、脑等。
- **数据**：RNA-seq + WES。处理后表达矩阵公开。
- **下载**：GEO [GSE110590](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110590)（RSEM 上分位数标准化）；后续 RAP 扩展与正常组织：[GSE193103](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193103)；原始 FASTQ dbGaP [phs000676](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000676.v2.p1)。
- **要点**：表达更像「同一患者」而不是「同一器官」；转移灶相对原发上调迁移/代谢，下调核酸加工。

### A3. AURORA EU / BIG（最大配对 RNA-seq，需申请）

- **Aftimos et al.**, *Cancer Discovery* (2021). DOI: [10.1158/2159-8290.CD-20-1647](https://doi.org/10.1158/2159-8290.cd-20-1647)
- **样本**：381 例转移性乳腺癌；配对原发–转移：TGS 252 对、RNA-seq 152 对、SNP array 67 对。转移活检部位有注释（肝、骨、肺、淋巴结、皮肤等）。
- **数据**：靶向基因测序 + Illumina RNA-seq。无蛋白组。
- **下载**：申请入口 [aurora-ctc.hubruxelles.be](https://aurora-ctc.hubruxelles.be/)；说明页列出 RNA.zip 原始 counts。
- **要点**：36% 发生 PAM50 亚型转换；转移灶免疫评分更低；`ESR1`/`PTEN`/`RB1` 等在转移灶富集。

### A4. 配对原发 vs 脑转移 RNA-seq

- **Cosgrove, Varešlija et al.**, *Nature Communications* (2022). DOI: [10.1038/s41467-022-27987-5](https://www.nature.com/articles/s41467-022-27987-5)
  - 45 例患者、90 个肿瘤（配对原发 + 脑转移）；exome-capture RNA-seq；部分病例另有 WES。
  - GEO：[GSE184869](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184869)；同队列处理后矩阵也见于 [GSE173661](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173661)（log2 TMM-CPM，16,714 个蛋白编码基因）。
- **Varešlija et al.**, *JNCI* (2019). DOI: [10.1093/jnci/djy110](https://doi.org/10.1093/jnci/djy110)
  - 更早的 21 对 TruSeq RNA-seq；脑转移相对原发 1,314 上调 / 1,702 下调（FC > 1.5, padj < 0.05）；突出 RET / HER2。

### A5. 配对原发 vs 肺 / 胸膜转移（靶向 RNA）

- **Sinn et al.**, *JCO Precision Oncology* (2020). DOI: [10.1200/PO.19.00337](https://doi.org/10.1200/PO.19.00337)
- **样本**：57 对 FFPE，原发 vs 活检证实的肺或胸膜转移。
- **数据**：NanoString 269 个乳腺癌基因 + PAM50，**不是**全转录组。
- **下载**：GEO [GSE145752](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145752)
- **要点**：Luminal A 最容易亚型转换（常转为 Luminal B）。

### A6. 配对原发–转移 NanoString（多器官，含肝/骨/肺）

- **Cejalvo et al.**, *Cancer Research* (2017). DOI: [10.1158/0008-5472.CAN-16-2479](https://aacrjournals.org/cancerres/article/77/9/2213/625048)
- **样本**：123 对原发 + 转移；转移部位：皮肤 35、淋巴结 24、肝 20、骨 16、肺 7、卵巢/腹膜 7、胸膜 6，另有脑等。含 GEICAM/2009-03 ConvertHER。
- **数据**：NanoString PAM50 / 乳腺癌基因，不是全转录组。AURORA US 后来对其中一部分做了 RNA-seq 并入合并队列。
- **下载**：GEO [GSE92977](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92977)
- **要点**：肝转移亚型转换率最高，肺最低。

### A7. 配对原发 vs 脑转移蛋白组（队列小于 RNA）

远处器官的**临床组织蛋白组**明显少于转录组。

- **JCO 2025 摘要 1050**，「Landscape analysis of proteins in the development of breast cancer brain metastasis」。DOI: [10.1200/jco.2025.43.16_suppl.1050](https://doi.org/10.1200/jco.2025.43.16_suppl.1050)
  - 21 例手术切除的配对原发 + 脑转移；DIA 质谱，9,430 个蛋白组，692 个差异蛋白（CRYAB、GFAP、ECM/胶原等）。会议摘要，全文/PRIDE 需跟进。
- **PXD073600**：标题为 primary TNBC and paired brain metastases 的蛋白组（iProX / ProteomeXchange）。适合作为脑转移蛋白组补充，使用前核对样本表是否均为人源配对组织。

---

## B. 人源原发 + 淋巴结转移（蛋白组相对完整；不是肺/骨/肝/脑）

淋巴结是区域转移，不是远处器官；蛋白组文献却主要在这一层。

- **Pozniak, Geiger et al.**, *Cell Systems* (2016). DOI: [10.1016/j.cels.2016.02.001](https://doi.org/10.1016/j.cels.2016.02.001)
  - 88 例临床样本：luminal 原发（LN− / LN+）、配对淋巴结转移、正常乳腺上皮。Super-SILAC，>10,000 蛋白。
  - 原发 vs 配对 LN 几乎相同（仅约 10 个显著蛋白）；原发 vs 正常差异很大。
  - PRIDE：[PXD000815](https://www.ebi.ac.uk/pride/archive/projects/PXD000815)
- **Franco et al.**, *Data in Brief* (2019). DOI: [10.1016/j.dib.2019.104125](https://doi.org/10.1016/j.dib.2019.104125)
  - 原发、腋窝淋巴结转移、对侧/癌旁非肿瘤乳腺；label-free nLC-MS/MS。
  - PRIDE：[PXD012431](https://www.ebi.ac.uk/pride/archive/projects/PXD012431)
- **Milioli et al.**, *OMICS* (2015)。7 对原发 vs 淋巴结转移，2DE-MS，队列小，作历史对照即可。

---

## C. 临床转移灶为主、原发很少或不配对

这些有肺/骨/肝/脑组织，但**本文本身缺少配对原发**。可与 TCGA / AURORA 原发对照，分析时必须写明「非配对」。

- **Brasó-Maristany et al.**, *Molecular Oncology* (2022). DOI: [10.1002/1878-0261.13021](https://doi.org/10.1002/1878-0261.13021)
  - 176 例患者、184 个转移灶、11 个器官（含肝、肺、脑、骨）；允许 7 例新发转移时的原发活检。
  - NanoString Breast Cancer 360（771 基因）。GEO [GSE175692](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175692)
  - 表达变异主要来自 PAM50 亚型而非器官；肺的免疫签名高于脑/肝。
- **Zhang, Massagué et al.**, *Cancer Cell* (2009). PMID: 19573813。GEO SuperSeries [GSE14020](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14020)（[GSE14017](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14017) U133Plus2；[GSE14018](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14018) U133A）
  - 58 个转移灶芯片：肺、脑、骨、肝。**没有原发灶**。
- **Robinson et al.**, MET500, *Nature* (2017). DOI: [10.1038/nature23306](https://doi.org/10.1038/nature23306)
  - 500 例转移实体瘤 WES + RNA-seq，含乳腺癌；肝相对多，肺/骨偏少。dbGaP [phs000673](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000673.v2.p1)。转移灶为主，不是乳腺癌配对原发。

---

## D. 单细胞 / 空间（临床活检）

- **Klughammer et al.**, HTAN, *Nature Medicine* (2024). DOI: [10.1038/s41591-024-03215-z](https://www.nature.com/articles/s41591-024-03215-z)
  - 60 例患者、67 个转移性乳腺癌活检：scRNA-seq 30、snRNA-seq 37；15 例匹配空间（Slide-seq / MERFISH / ExSeq / CODEX）。
  - 部位：肝 37、腋窝 9、乳腺 7（转移诊断后的原发部位）、骨 5、胸壁 3、颈 3、脑 1、肺 1、皮肤 1。
  - 数据走 HTAN / CELLxGENE，不在本仓库 Cuffdiff 流程里分析。
- **Labrie et al.**, SMMART / HTAN, *Cell Reports Medicine* (2022). DOI: [10.1016/j.xcrm.2022.100525](https://doi.org/10.1016/j.xcrm.2022.100525)
  - **1 例患者** 3.5 年：原发、多次肝活检、骨活检；DNA / RNA / 蛋白 / 多重空间成像。适合个案，不适合队列统计。
- **GSE225600**：4 例患者原发 + 配对转移淋巴结的 scRNA + Visium。淋巴结，不是肺/骨/肝/脑。
- **GSE190772**：2 例骨转移 scRNA（转移灶为主）。

---

## E. 易混入、但不符合「原发组织 + 转移组织」

| 数据集 | 实际是什么 | 为什么不够 |
| --- | --- | --- |
| GSE2603（Minn / MSK） | 原发灶芯片 + 肺/骨转移**随访结局** | 没有转移灶组织 |
| GSE5327 | ER− 原发 + 肺转移结局 | 同上 |
| GSE12276（Bos 脑转移签名） | 原发灶 + 脑转移风险 | 同上 |
| Cimas *ecancer* 2019 | 51 例 TNBC **原发** PRM 蛋白 + 日后是否 CNS 转移 | 没有脑转移组织 |
| TCGA-BRCA、CPTAC-BRCA | 手术原发为主 | 几乎没有配对远处转移组织 |
| MDA-MB-231 肺/骨/脑亚系、多数 PDX | 实验模型 | 禁止当病人 |

GSE2603 / GSE5327 只适用于「原发灶预测哪一类转移结局」，不能回答「转移灶相对原发灶的 RNA/蛋白差异」。

---

## 数据类型小结

1. **同时有人源原发和远处转移、且 RNA 可下载或可申请**：AURORA US、RAP、AURORA EU、Cosgrove/Varešlija 脑转移、Cejalvo/ConvertHER、Sinn 肺/胸膜。
2. **肺、骨、肝、脑在同一研究里都出现**：AURORA US、RAP、Cejalvo（骨/肝/肺明确，脑很少）、Brasó-Maristany（转移灶）、Zhang GSE14020（仅转移灶）、Klughammer（肝多、骨/肺/脑少）。
3. **蛋白组**：淋巴结层最完整（Pozniak PXD000815）。远处器官几乎只有配对脑转移的小队列（JCO 2025 摘要、PXD073600）。公开的临床「原发 vs 肺/骨/肝」深度蛋白组仍然很少。
4. **单细胞/空间**：HTAN Klughammer 是目前最大的临床转移活检图谱，但原发配对有限，肝活检占多数。
