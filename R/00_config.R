# TCGA-BRCA × 神经元/轴突 GO：通路列表与分析参数
# 增删通路只改 GO_TERMS。每个 ID 必须被独立分析，禁止合并。

GO_TERMS <- c(
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

# 输入文件名片段（在 --data-dir 下自动匹配，支持 .tsv/.txt/.csv/.gz）
DATA_FILE_KEYS <- list(
  clinical = "TCGA-BRCA.clinical",
  protein  = "TCGA-BRCA.protein",
  expr     = "TCGA-BRCA.star_fpkm",
  survival = "TCGA-BRCA.survival",
  gencode  = "gencode.v36.annotation.gtf.gene"
)

ANALYSIS_PARAMS <- list(
  sample_type_prefix     = "01",   # 原发肿瘤 01；设为 NULL 则不过滤样本类型
  min_go_genes           = 3L,     # 映射到表达矩阵后少于此数则跳过该 GO
  min_expr_frac          = 0.20,   # 基因在至少 20% 样本中表达才进入全基因组相关
  protein_coding_only    = TRUE,
  go_keytype             = "GOALL", # 包含子术语注释；每个 GO 仍单独取基因
  min_spearman_rho_neg   = -0.20,  # 负相关强度下限（同时要求 rho < 0）
  fdr_cutoff             = 0.05,
  km_quantile            = 0.50,   # KM 高低分组：通路得分中位数
  seed                   = 1L
)
