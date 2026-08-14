# 用户可改配置。GO ID 与注释保持原文，禁止把通路合并成一个基因集。

data_dir <- "data"
output_dir <- "results"

# 文件前缀（允许 .tsv / .txt / .gz 等后缀）
file_prefixes <- list(
  clinical   = "TCGA-BRCA.clinical",
  protein    = "TCGA-BRCA.protein",
  expression = "TCGA-BRCA.star_fpkm",
  survival   = "TCGA-BRCA.survival",
  gtf        = "gencode.v36.annotation.gtf.gene"
)

# verbatim
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

go_names <- c(
  "GO:0023041" = "neuronal signal transduction",
  "GO:1904457" = "positive regulation of neuronal action potential",
  "GO:1904340" = "positive regulation of dopaminergic neuron differentiation",
  "GO:2001224" = "positive regulation of neuron migration",
  "GO:2001222" = "regulation of neuron migration",
  "GO:0019227" = "neuronal action potential propagation",
  "GO:0019228" = "neuronal action potential",
  "GO:1902847" = "regulation of neuronal signal transduction",
  "GO:0031102" = "neuron projection regeneration",
  "GO:0097492" = "sympathetic neuron axon guidance",
  "GO:0097491" = "sympathetic neuron projection guidance",
  "GO:0097374" = "sensory neuron axon guidance",
  "GO:0007158" = "neuron cell-cell adhesion",
  "GO:1902667" = "regulation of axon guidance",
  "GO:0031103" = "axon regeneration",
  "GO:0007411" = "axon guidance",
  "GO:0007409" = "axonogenesis"
)

# 通路基因映射：GOALL 含该 term 的子 term
go_keytype <- "GOALL"

# 评分：ssgsea 优先；GSVA 不可用时用 mean_z
score_method <- "ssgsea"
min_genes_per_go <- 5L
log2_fpkm <- TRUE
auto_install_packages <- TRUE

# 样本：默认只保留原发瘤 01
keep_sample_types <- c("01")
barcode_sample_chars <- 15L
barcode_patient_chars <- 12L

# 负相关基因：每个 GO 内部单独 BH
fdr_cutoff <- 0.05
rho_cutoff <- 0
min_paired_samples <- 20L

# 生存 KM 按该通路评分中位数分组
km_split <- "median"

# 临床列：NULL 表示按关键词自动挑选
clinical_columns <- NULL
clinical_keywords <- c(
  "age", "stage", "grade", "er_status", "pr_status", "her2",
  "subtype", "pam50", "race", "gender", "sex", "laterality",
  "histolog", "pathologic_t", "pathologic_n", "pathologic_m",
  "tumor_status", "menopause", "lymph_node"
)
