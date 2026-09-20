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
| 原发 vs 淋巴结蛋白组 | Pozniak / Geiger（PXD000815）；复旦余科达/邵志敏 65 对（Protein & Cell 2026） |
| 单细胞 / 空间、临床活检 | Klughammer HTAN（Nat Med 2024） |
| 复旦邵志敏 / 江一舟组 | 见 **F 节**：原发多组学很大；远处转移公开数据几乎全是 DNA panel / IHC |
| 中山大学宋尔卫组 | 见 **G 节**：机制强；临床组织主要是 IHC/IF（含肝转移），没有公开配对 RNA/蛋白组 |

**不要**把下列数据当成「人源转移组织」：`GSE2603`、`GSE5327`（原发灶 + 肺转移**临床结局**，没有转移灶活检）、TCGA-BRCA / CPTAC-BRCA / **CBCGA 与 FUSCC 原发 TNBC**（几乎全是原发）、MDA-MB-231 衍生肺/骨/脑亚系。

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
- **Liu, Yu, Shao et al.**, *Advanced Science* (2023)；GEO [GSE225600](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225600)
  - 复旦肿瘤 **4 对**原发 + 腋窝转移淋巴结，scRNA + Visium。通讯余科达，邵志敏为共同作者。淋巴结，不是肺/骨/肝/脑。
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

## F. 复旦大学附属肿瘤医院：邵志敏 / 江一舟课题组

按原筛选标准核对：**临床组织 + 人源原发 + 人源转移（优先肺/骨/肝/脑）+ RNA 或蛋白**。

**结论：该组目前没有公开的「配对原发 vs 肺/骨/肝/脑转移组织」全转录组或质谱蛋白组。** 组学主力是原发灶；转移相关公开队列主要是 484 基因 DNA panel 和 IHC。淋巴结层的 RNA/蛋白在同科室余科达组（邵志敏参与，江一舟不是通讯）。脑转移 RNA+蛋白在同医院内科胡夕春/张剑组，不是邵/江组。

### F1. 最接近标准：原发 + 转移组织，但是淋巴结，且通讯多为余科达

| 文章 | 组织 | 数据 | 作者关系 | 入口 |
| --- | --- | --- | --- | --- |
| **Zhai, Yu, Shao et al.**, *Protein & Cell* (2026). DOI: [10.1093/procel/pwag002](https://doi.org/10.1093/procel/pwag002) | 65 对治疗naive：**原发 + 腋窝转移淋巴结 + 癌旁**（195 份） | WES；RNA-seq（约 55 原发 / 58 淋巴结 / 51 癌旁）；TMT 蛋白组 65+65；磷酸化组 65 原发 / 55 淋巴结（8,588 蛋白，24,178 磷酸化位点） | 邵志敏作者；通讯余科达、谭敏佳、金怡婷。**江一舟不在作者中** | 原文补充表；ProteomeXchange 号需跟全文 Data availability |
| **Liu, Yu, Shao et al.**, *Adv Sci* (2023) | 4 对原发 + 腋窝转移淋巴结 | scRNA + Visium | 通讯余科达；邵志敏共同作者 | GEO `GSE225600` |

淋巴结是区域转移，不要写成肺/骨/肝/脑。该文结论是：淋巴结转移差异主要在**蛋白/磷酸化**，基因组和转录组差别较小（ANGPTL4、HMGB1、MARCKSL1-S104、FKBP15-S320、PRKCB）。

### F2. 邵志敏 + 江一舟：有转移组织，但是 DNA / IHC，不是 RNA 或蛋白组

这些文章**有**人源转移活检（含肝、肺、淋巴结、胸壁；临床事件含骨），也常有配对原发，但分子层是靶向测序或免疫组化。

- **Zhu, Jiang, Shao, Wang et al.**, *JCI* (2025/2026). DOI: [10.1172/JCI188989](https://www.jci.org/articles/view/188989)
  - 296 例转移性 TNBC 的**转移灶** 484 基因 panel；临床转移事件：淋巴结 213、肺 142、骨 129、肝 99、胸壁 93。活检部位以淋巴结、肝、肺、胸壁为主。另 105 对原发–转移验证 `PKD1`。
  - RNA-seq 出现在机制实验（细胞系/小鼠），不是病人转移灶转录组。
  - 数据：NODE `OEZ00021764`、`OEZ00021765`。
- **Jiang, Shao et al.**, FUTURE 试验, *Cell Research* (2021). DOI: [10.1038/s41422-020-0375-9](https://www.nature.com/articles/s41422-020-0375-9)
  - 难治转移性 TNBC：对可及转移灶再活检，做 IHC 亚型（AR / CD8 / FOXC1）+ 484 基因 panel。部位含淋巴结、肺、肝、骨、胸壁、乳腺。**没有**转移灶 RNA-seq/蛋白组。
- **Zhu, Jiang, Shao et al.**, *Cancer Biology & Medicine* (2024). DOI: [10.20892/j.issn.2095-3941.2024.0009](https://doi.org/10.20892/j.issn.2095-3941.2024.0009)
  - 880 例（465 早期 + 415 转移）；40 对原发–转移。亚型用 IHC，不是转录组。MES 亚型容易转换并偏脑转移。

### F3. 邵 / 江组大规模组学：几乎全是原发灶（有转移随访，没有转移组织组学）

| 队列 | 文章 | 样本 | 数据 | 注意 |
| --- | --- | --- | --- | --- |
| FUSCC TNBC | Jiang, Shao et al., *Cancer Cell* (2019). DOI: [10.1016/j.ccell.2019.02.001](https://doi.org/10.1016/j.ccell.2019.02.001) | 465 例**原发** TNBC；随访中 65 例复发/转移 | WES 279、CNA 401、RNA-seq 360 | 原发 + 结局，不是转移组织 |
| CBCGA | Jiang, Shao et al., *Nature Cancer* (2024). DOI: [10.1038/s43018-024-00725-0](https://www.nature.com/articles/s43018-024-00725-0) | 773 例中国乳腺癌**原发** | 基因组、转录组 752、蛋白组 278、代谢组等 | 亚洲最大原发多组学，不是转移图谱 |
| FUSCC-BRCA | Ma, Jiang, Shao et al., *Cancer Cell* (2024). DOI: [10.1016/j.ccell.2024.03.006](https://doi.org/10.1016/j.ccell.2024.03.006) | 873 例亚洲乳腺癌**原发** | WES+CNA 873、RNA-seq 842、TMT 蛋白 261、代谢 509 | 同上 |
| TNBC 蛋白组 | Gong, Jiang, Shao et al., *Cell Reports* (2022). DOI: [10.1016/j.celrep.2022.110460](https://doi.org/10.1016/j.celrep.2022.110460) | 90 例 TNBC **原发** | 蛋白组 / 磷酸化 / 转录因子占用 | 原发 |
| HER2-low 蛋白 | Dai, Jiang, Shao et al., *Nat Commun* (2023) | 中国乳腺癌原发 | TMT 蛋白组 | iProX `PXD042886` |

这些可以做「原发灶预测转移风险」，**不能**做「转移灶相对原发灶的 RNA/蛋白差异」。

### F4. 同医院、但不是邵志敏 / 江一舟组（避免混进）

- **Lin, Zhang Jian, Hu Xichun et al.**, *Nature Communications* (2026). DOI: [10.1038/s41467-026-72927-2](https://www.nature.com/articles/s41467-026-72927-2)
  - 复旦肿瘤内科：TNBC **原发 vs 脑转移** 的 RNA-seq + 蛋白组 + 代谢组；另有配对脑脊液/血浆。NDUFB9。
  - 数据：NGDC [PRJCA032856](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA032856)、OMIX 代谢组。这是 FUSCC 里最接近「原发 + 脑转移 RNA/蛋白」的公开队列，课题组是张剑/胡夕春，不是邵/江。

---

## G. 中山大学孙逸仙纪念医院：宋尔卫课题组

按同一标准核对。宋组以**转移机制**见长（GM-CSF–CCL18、NET–CCDC25、外泌体 HISLA），临床组织主要用于 IHC / IF / 血清，而不是 bulk RNA-seq 或组织质谱。

**结论：没有公开的「配对原发 vs 肺/骨/肝/脑」全转录组或组织蛋白组。** 最接近的是人源原发 + 肝转移灶的免疫荧光（有组织，不是 RNA/蛋白矩阵）。

### G1. 有人源原发和肝转移组织，但是 IHC / IF，不是 RNA 或质谱

- **Yang, Su, Song et al.**, *Nature* (2020). DOI: [10.1038/s41586-020-2394-6](https://doi.org/10.1038/s41586-020-2394-6)
  - 对人源**原发灶和转移灶**做 MPO / H3Cit 免疫荧光；**肝转移** NET 浸润最多。早期乳腺癌血清 MPO–DNA 可预测日后肝转移。原发灶 CCDC25 IHC 与预后相关。
  - 机制：NET-DNA 经癌细胞受体 CCDC25 → ILK–β-parvin。小鼠肝/肺模型。
  - **没有**病人转移灶 RNA-seq / 蛋白组矩阵。
- **Zhang, Song, Yang et al.**, *Nature Communications* (2026). DOI: [10.1038/s41467-026-76459-7](https://www.nature.com/articles/s41467-026-76459-7)
  - 临床肝转移标本验证 NET 与 NK 功能障碍；纵向 scRNA-seq 来自 **4T1 小鼠**肝转移，不是病人组织转录组。
- **Su, Liu, Song et al.**, *Cancer Cell* (2014). DOI: [10.1016/j.ccr.2014.03.021](https://doi.org/10.1016/j.ccr.2014.03.021)
  - 约 1015 例**原发** IHC（GM-CSF / CCL18 / EMT）+ 151 例血清。GEO `GSE51938` 是细胞系条件培养基细胞因子芯片，不是转移组织。
- **Chen, Song et al.**, *Cancer Cell* (2011)：CCL18–PITPNM3，同样是原发组织 IHC 为主。
- **Chen, Su, Song et al.**, *Nature Cell Biology* (2019)：TAM 外泌体 HISLA，原发切片 IHC。

### G2. 转移性乳腺癌临床队列，但是血浆蛋白 / 原发或可及病灶 IHC

- **Liu Jieqiong, Song et al.**, *Nature Communications* (2022). DOI: [10.1038/s41467-022-30569-0](https://www.nature.com/articles/s41467-022-30569-0)
  - NCT04303741：46 例晚期 TNBC（局部晚期或转移）。Olink 做的是**血浆**免疫肿瘤 panel；组织侧是 IHC / 多重荧光（TLS、PML、PLOD3），不是配对原发–转移转录组或组织质谱。

### G3. 不要和同校其他医院、或其他单位的单细胞队列搞混

| 容易混进来的 | 实际是谁 | 数据 | 为什么不是宋组 |
| --- | --- | --- | --- |
| Zou, Tang Hailin et al., *Advanced Science* (2023). DOI: [10.1002/advs.202203699](https://doi.org/10.1002/advs.202203699) | **中山大学肿瘤防治中心**唐海林，不是孙逸仙纪念医院 | 6 例乳腺癌**肝或脑转移** scRNA，44,473 细胞 | 同大学不同医院、不同课题组 |
| Xu et al., *Oncogenesis* (2021)；GEO `GSE180286` | 南京医大一附院管晓翔等 | 5 例原发 + 10 个配对淋巴结 scRNA | 不是中大宋组 |
| Liu, Yu, Shao；`GSE225600` | 复旦余科达 | 4 对原发 + 腋窝淋巴结 scRNA+Visium | 复旦，不是宋组 |

宋组若要做「原发 vs 肝转移 RNA/蛋白」，公开数据不够，需要向组内要原始组织组学，或用 AURORA/RAP / 唐海林肝脑 scRNA 作为外源队列，并写明来源。

---

## 数据类型小结

1. **同时有人源原发和远处转移、且 RNA 可下载或可申请**：AURORA US、RAP、AURORA EU、Cosgrove/Varešlija 脑转移、Cejalvo/ConvertHER、Sinn 肺/胸膜。复旦邵/江组、中大宋尔卫组目前都没有对等的公开队列。
2. **肺、骨、肝、脑在同一研究里都出现**：AURORA US、RAP、Cejalvo（骨/肝/肺明确，脑很少）、Brasó-Maristany（转移灶）、Zhang GSE14020（仅转移灶）、Klughammer（肝多、骨/肺/脑少）。邵/江组 JCI 296 例转移性 TNBC 覆盖这些器官，但是 **DNA panel**。
3. **蛋白组**：淋巴结层最完整（Pozniak PXD000815；复旦余科达/邵志敏 65 对 TMT+磷酸化）。远处器官几乎只有配对脑转移的小队列（JCO 2025 摘要、PXD073600、FUSCC 张剑/胡夕春 NDUFB9）。公开的临床「原发 vs 肺/骨/肝」深度蛋白组仍然很少。
4. **单细胞/空间**：HTAN Klughammer 是目前最大的临床转移活检图谱，但原发配对有限，肝活检占多数。复旦 `GSE225600` 是原发+淋巴结。
5. **复旦邵/江组怎么用**：要 RNA/蛋白请用他们的**原发**队列（CBCGA / FUSCC TNBC）或余科达组淋巴结多组学；要肺/肝/骨转移组织目前只能用他们的 **DNA/IHC** 转移队列，或回到 AURORA/RAP。
6. **宋尔卫组怎么用**：机制和肝转移 **IHC/IF**（CCDC25 / NET）可用；不要把该组写成有公开配对 RNA/蛋白图谱。同校肝/脑转移 scRNA 属于中肿唐海林，不是宋组。
