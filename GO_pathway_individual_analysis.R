# ============================================================
# TCGA-BRCA：每个 GO 通路单独的临床相关性 & 负相关基因分析
# 核心原则：每个 GO ID 单独取基因、单独打分、单独统计，绝不把通路合并
#
# 用法：将本脚本放在数据目录（默认 E:/R/BRCA）后，RStudio 中 Source 整份脚本
# 或：setwd("E:/R/BRCA"); source("GO_pathway_individual_analysis.R")
# 结果目录：results_GO_individual/  （每个 GO 一个子文件夹）
#
# ★★★ 请重新 Source 整份脚本，不要只接着跑旧的 score[...] 行 ★★★
# 已修改处均用  “★★★ 已修改”  标记，全文搜索即可定位。
#   1) pathway_zmean：不再使用 scale()/t()
#   2) 通路分数变量由 score 改为 go_score（score 会撞上已有函数）
#   3) 临床箱线图：不再只打印 null device 1，改为报告保存了几张图
#   4) 生存分析 7.2：同样禁止 score[...]；打分成功后会把 go_score 同步到 score 以兼容旧行
#   5) 汇总热图：禁止 t(scale(score_mat))，改用手写 z-score
#   6) 热图绘制：必须用 ComplexHeatmap::Heatmap，避免被 heatmaps::Heatmap 覆盖
#   7) OS 森林图：ggplot2 4.0 用 geom_errorbar(orientation="y")，不要漏掉行末 +
#   8) 多个 GO：已写在 go_list；用 go_to_run 控制本次跑哪些
#   9) 设完 go_to_run 后：先 source("run_go_individual_analysis.R")，再 run_go_individual_analysis()
#      两个 .R 文件必须都放在 E:/R/BRCA
# ============================================================

## 0. 可调参数 ------------------------------------------------
work_dir <- "E:/R/BRCA"   # 与当前工程一致；若不存在则使用当前目录
use_primary_tumor_only <- TRUE
min_pathway_genes <- 2
neg_fdr_cutoff <- 0.05
neg_r_cutoff <- 0          # 负相关：Spearman r < 0
strict_r_cutoff <- -0.3    # 高置信负相关附加阈值
out_dir <- "results_GO_individual"

# 用户指定的 GO 通路（每个都单独分析）
go_list <- c(
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

go_name_map <- c(
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

# ★★★ 已修改开始：一次跑多个 GO，不用在控制台逐个输入 ★★★
# go_list 里已经有你列的 17 个通路。默认全部单独分析。
# 若只想补跑某几个（例如还没有文件夹的），改成：
#   go_to_run <- c("GO:0007411", "GO:0007409")
# 若要加新通路：先写进 go_list 和 go_name_map，再保持 go_to_run <- go_list
go_to_run <- go_list
# ★★★ 已修改结束 ★★★

## 1. 加载包 --------------------------------------------------
# 数据处理
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(stringr)
library(purrr)

# 生信 / TCGA
library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# 生存分析
library(survival)
library(survminer)

# 可视化
library(ggplot2)
library(ggpubr)
library(enrichplot)
library(ComplexHeatmap)
library(EnhancedVolcano)
# 若会话里加载过 heatmaps 包，Heatmap() 会被盖掉，后面一律写 ComplexHeatmap::Heatmap

## 2. 工作目录与数据读取 --------------------------------------
if (dir.exists(work_dir)) {
  setwd(work_dir)
} else {
  message("未找到 ", work_dir, " ，改用当前工作目录：", getwd())
}

need_files <- c(
  "TCGA-BRCA.star_fpkm.tsv",
  "gencode.v36.annotation.gtf.gene.probemap",
  "TCGA-BRCA.clinical.tsv",
  "TCGA-BRCA.survival.tsv"
)
miss <- need_files[!file.exists(need_files)]
if (length(miss) > 0) {
  stop("以下文件不在工作目录中：\n  ", paste(miss, collapse = "\n  "))
}

fpkm_data     <- fread("TCGA-BRCA.star_fpkm.tsv")
probe_annot   <- fread("gencode.v36.annotation.gtf.gene.probemap")
clinical_data <- fread("TCGA-BRCA.clinical.tsv")
survival_data <- fread("TCGA-BRCA.survival.tsv")
protein_data  <- if (file.exists("TCGA-BRCA.protein.tsv")) fread("TCGA-BRCA.protein.tsv") else NULL
if (is.null(protein_data)) message("未找到 TCGA-BRCA.protein.tsv，将跳过蛋白相关性")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "per_GO"), showWarnings = FALSE)

## 3. 工具函数 ------------------------------------------------
normalize_barcode <- function(x) {
  x <- toupper(gsub("\\.", "-", as.character(x)))
  x <- sub("A$", "", x)
  ifelse(nchar(x) >= 15, substr(x, 1, 15), x)
}

sample_type_code <- function(ids) {
  substr(normalize_barcode(ids), 14, 15)
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

first_present <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) NA_character_ else hit[1]
}

get_go_genes <- function(go_id) {
  pick <- function(keytype) {
    tryCatch(
      AnnotationDbi::select(
        org.Hs.eg.db,
        keys = go_id,
        keytype = keytype,
        columns = c("SYMBOL", "ENSEMBL", "ENTREZID")
      ),
      error = function(e) NULL
    )
  }
  res <- suppressMessages(pick("GOALL"))
  if (is.null(res) || nrow(res) == 0) res <- suppressMessages(pick("GO"))
  if (is.null(res) || nrow(res) == 0) {
    return(data.table(GO = go_id, SYMBOL = character(), ENSEMBL = character(), ENTREZID = character()))
  }
  res <- as.data.table(res)
  if ("GOALL" %in% names(res)) setnames(res, "GOALL", "GO", skip_absent = TRUE)
  if (!"GO" %in% names(res)) res[, GO := go_id]
  res <- unique(res[!is.na(SYMBOL) & SYMBOL != "", .(GO, SYMBOL, ENSEMBL, ENTREZID)])
  res
}

# ★★★ 已修改开始：pathway_zmean 不再调用 scale()/t() ★★★
pathway_zmean <- function(expr_mat, genes) {
  genes <- unique(intersect(as.character(genes), rownames(expr_mat)))
  if (length(genes) < min_pathway_genes) return(NULL)
  sub <- as.matrix(expr_mat[genes, , drop = FALSE])
  storage.mode(sub) <- "double"
  gene_mean <- rowMeans(sub, na.rm = TRUE)
  gene_sd <- sqrt(rowMeans((sub - gene_mean)^2, na.rm = TRUE))
  gene_sd[!is.finite(gene_sd) | gene_sd < 1e-12] <- 1
  z <- (sub - gene_mean) / gene_sd
  z[!is.finite(z)] <- 0
  sc <- colMeans(z, na.rm = TRUE)
  names(sc) <- colnames(sub)
  attr(sc, "n_genes") <- length(genes)
  attr(sc, "genes") <- genes
  sc
}
# ★★★ 已修改结束 ★★★

spearman_vs_score <- function(mat, score) {
  common <- intersect(colnames(mat), names(score))
  if (length(common) < 10) return(data.table())
  mat <- mat[, common, drop = FALSE]
  score <- score[common]
  keep <- apply(mat, 1, function(x) sd(x, na.rm = TRUE) > 0)
  mat <- mat[keep, , drop = FALSE]
  n <- ncol(mat)
  r <- as.numeric(cor(base::t(as.matrix(mat)), score, method = "spearman", use = "pairwise.complete.obs"))  # ★★★ 已修改：base::t
  names(r) <- rownames(mat)
  r <- pmin(pmax(r, -0.999999), 0.999999)
  tstat <- r * sqrt((n - 2) / pmax(1e-12, 1 - r^2))
  p <- 2 * pt(-abs(tstat), df = n - 2)
  data.table(
    feature = names(r),
    spearman_r = r,
    pvalue = p,
    fdr = p.adjust(p, method = "BH"),
    n = n
  )
}

assoc_clinical_feature <- function(score, feature, feature_name) {
  df <- data.frame(score = as.numeric(score), feature = feature, stringsAsFactors = FALSE)
  df <- df[is.finite(df$score) & !is.na(df$feature), , drop = FALSE]
  df$feature <- as.character(df$feature)
  df <- df[df$feature != "" & !tolower(df$feature) %in% c("na", "nan", "not reported", "unknown", "[not available]", "[unknown]"), ]
  if (nrow(df) < 10) return(NULL)

  feat_num <- suppressWarnings(as.numeric(df$feature))
  if (mean(is.finite(feat_num)) >= 0.8) {
    df$feature <- feat_num
    df <- df[is.finite(df$feature), ]
    if (nrow(df) < 10 || sd(df$feature) == 0) return(NULL)
    ct <- suppressWarnings(cor.test(df$score, df$feature, method = "spearman", exact = FALSE))
    return(data.table(
      feature = feature_name, type = "continuous", n = nrow(df),
      stat = unname(ct$estimate), pvalue = ct$p.value,
      method = "Spearman", detail = sprintf("rho=%.3f", unname(ct$estimate))
    ))
  }

  df$feature <- droplevels(factor(df$feature))
  tab <- table(df$feature)
  df <- df[df$feature %in% names(tab)[tab >= 3], ]
  df$feature <- droplevels(factor(df$feature))
  if (nlevels(df$feature) < 2 || nlevels(df$feature) > 15) return(NULL)

  if (nlevels(df$feature) == 2) {
    wt <- suppressWarnings(wilcox.test(score ~ feature, data = df))
    return(data.table(
      feature = feature_name, type = "categorical", n = nrow(df),
      stat = unname(wt$statistic), pvalue = wt$p.value,
      method = "Wilcoxon",
      detail = paste(paste0(names(table(df$feature)), "=", as.integer(table(df$feature))), collapse = "; ")
    ))
  }
  kt <- suppressWarnings(kruskal.test(score ~ feature, data = df))
  data.table(
    feature = feature_name, type = "categorical", n = nrow(df),
    stat = unname(kt$statistic), pvalue = kt$p.value,
    method = "Kruskal-Wallis",
    detail = paste(paste0(names(table(df$feature)), "=", as.integer(table(df$feature))), collapse = "; ")
  )
}

simplify_stage <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("IV|STAGE.?4|\\b4\\b", x)] <- "Stage IV"
  out[is.na(out) & grepl("III|STAGE.?3|\\b3\\b", x)] <- "Stage III"
  out[is.na(out) & grepl("II|STAGE.?2|\\b2\\b", x)] <- "Stage II"
  out[is.na(out) & grepl("I|STAGE.?1|\\b1\\b", x)] <- "Stage I"
  out[grepl("X|NA|NOT|UNKNOWN", x)] <- NA_character_
  out
}

detect_id_col <- function(dt, prefer) {
  nms <- names(dt)
  hit <- first_present(nms, prefer)
  if (!is.na(hit)) return(hit)
  nms[1]
}

as_sample_matrix <- function(dt, id_prefer) {
  id_col <- detect_id_col(dt, id_prefer)
  ids <- as.character(dt[[id_col]])
  mat <- as.matrix(dt[, !id_col, with = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- ids
  mat
}

## 4. 表达矩阵预处理 ------------------------------------------
gene_id_col <- names(fpkm_data)[1]
expr_ids <- as.character(fpkm_data[[gene_id_col]])
expr <- as.matrix(fpkm_data[, -1, with = FALSE])
storage.mode(expr) <- "numeric"
rownames(expr) <- expr_ids
colnames(expr) <- normalize_barcode(colnames(expr))

# 注释：probemap 通常为 id / gene
annot_id_col   <- first_present(names(probe_annot), c("id", "Id", "gene_id", names(probe_annot)[1]))
annot_gene_col <- first_present(names(probe_annot), c("gene", "Gene", "symbol", "Symbol", names(probe_annot)[2]))
annot <- unique(probe_annot[, .(ensembl = as.character(get(annot_id_col)),
                                symbol  = as.character(get(annot_gene_col)))])
annot[, ensembl_nv := sub("\\..*$", "", ensembl)]

map_dt <- data.table(
  ensembl = rownames(expr),
  ensembl_nv = sub("\\..*$", "", rownames(expr))
)
map_dt <- merge(map_dt, unique(annot[, .(ensembl, symbol)]), by = "ensembl", all.x = TRUE)
map_dt <- merge(
  map_dt,
  unique(annot[, .(ensembl_nv, symbol_nv = symbol)]),
  by = "ensembl_nv",
  all.x = TRUE
)
map_dt[, symbol := fifelse(is.na(symbol) | symbol == "", symbol_nv, symbol)]

# 重复基因符号保留平均表达最高者
gene_mean <- rowMeans(expr, na.rm = TRUE)
map_dt[, mean_exp := gene_mean[ensembl]]
map_dt <- map_dt[!is.na(symbol) & symbol != "" & !symbol %in% c(".", "-")]
setorder(map_dt, -mean_exp)
map_dt <- map_dt[!duplicated(symbol)]
expr <- expr[map_dt$ensembl, , drop = FALSE]
rownames(expr) <- map_dt$symbol

# 仅保留原发肿瘤（01）；如需全部样本可把 use_primary_tumor_only 设为 FALSE
if (use_primary_tumor_only) {
  keep_s <- sample_type_code(colnames(expr)) == "01"
  if (sum(keep_s) >= 20) {
    expr <- expr[, keep_s, drop = FALSE]
    message("仅保留原发肿瘤样本：", ncol(expr))
  } else {
    message("原发肿瘤样本过少，改用全部样本：", ncol(expr))
  }
}

# 去重复样本列
if (any(duplicated(colnames(expr)))) {
  expr <- expr[, !duplicated(colnames(expr)), drop = FALSE]
}

mx <- suppressWarnings(max(expr, na.rm = TRUE))
if (is.finite(mx) && mx > 50) {
  message("检测到原始 FPKM（max=", round(mx, 2), "），进行 log2(x+1) 转换")
  expr <- log2(expr + 1)
} else {
  message("表达值范围较小（max=", round(mx, 2), "），视为已 log 转换")
}

## 5. 临床 / 生存 / 蛋白对齐 ----------------------------------
clin_id <- detect_id_col(clinical_data, c("sampleID", "sample", "bcr_patient_barcode", "submitter_id"))
surv_id <- detect_id_col(survival_data, c("sample", "sampleID", "bcr_patient_barcode"))

clinical_data <- copy(clinical_data)
survival_data <- copy(survival_data)
clinical_data[, sample_std := normalize_barcode(get(clin_id))]
survival_data[, sample_std := normalize_barcode(get(surv_id))]
clinical_data <- clinical_data[!duplicated(sample_std)]
survival_data <- survival_data[!duplicated(sample_std)]

# 临床字段：优先常见 BRCA 变量，并纳入分期/受体/组织学等匹配列
preferred_clin <- c(
  "age_at_initial_pathologic_diagnosis", "age_at_diagnosis", "age",
  "ajcc_pathologic_tumor_stage", "pathologic_stage", "clinical_stage",
  "ajcc_tumor_pathologic_pt", "ajcc_nodes_pathologic_pn", "ajcc_metastasis_pathologic_pm",
  "er_status_by_ihc", "pr_status_by_ihc", "her2_status_by_ihc",
  "ER.Status", "PR.Status", "HER2.Final.Status",
  "histological_type", "histologic_diagnosis", "histology",
  "race", "ethnicity", "gender", "menopause_status", "tumor_status"
)
skip_clin <- unique(c(
  clin_id, "sample_std", "sample", "sampleID", "patient", "_PATIENT",
  "bcr_patient_barcode", "submitter_id", "project_id", "days_to_collection"
))
auto_clin <- grep(
  "age|stage|er_|pr_|her2|histolog|race|gender|menopause|subtype|pam50|tumor_status|pathologic|grade",
  names(clinical_data), ignore.case = TRUE, value = TRUE
)
clin_features <- unique(c(preferred_clin[preferred_clin %in% names(clinical_data)], auto_clin))
clin_features <- setdiff(clin_features, skip_clin)

# 分期简化列
stage_col <- first_present(names(clinical_data), c(
  "ajcc_pathologic_tumor_stage", "pathologic_stage", "clinical_stage", "ajcc_pathologic_stage"
))
if (!is.na(stage_col)) {
  clinical_data[, stage_simplified := simplify_stage(get(stage_col))]
  clin_features <- unique(c("stage_simplified", clin_features))
}

# 蛋白矩阵：自动判断样本在行还是在列
prot_mat <- NULL
if (!is.null(protein_data) && nrow(protein_data) > 0) {
  prot_mat <- tryCatch({
    prot_c1 <- as.character(protein_data[[1]])
    if (mean(grepl("^TCGA", prot_c1, ignore.case = TRUE), na.rm = TRUE) > 0.5) {
      tmp <- as_sample_matrix(protein_data, names(protein_data)[1])
      tmp <- t(tmp)
    } else {
      tmp <- as_sample_matrix(protein_data, names(protein_data)[1])
    }
    colnames(tmp) <- normalize_barcode(colnames(tmp))
    if (any(duplicated(colnames(tmp)))) {
      tmp <- tmp[, !duplicated(colnames(tmp)), drop = FALSE]
    }
    rownames(tmp) <- make.unique(as.character(rownames(tmp)))
    tmp
  }, error = function(e) {
    message("蛋白矩阵解析失败，跳过蛋白分析：", conditionMessage(e))
    NULL
  })
}

common_samples <- colnames(expr)
message("表达矩阵：", nrow(expr), " 基因 x ", ncol(expr), " 样本")
message("与临床重叠：", length(intersect(common_samples, clinical_data$sample_std)))
message("与生存重叠：", length(intersect(common_samples, survival_data$sample_std)))
message("与蛋白重叠：", if (is.null(prot_mat)) 0 else length(intersect(common_samples, colnames(prot_mat))))

## 6-8. 真正开始分析 ------------------------------------------
# ★★★ 已修改：函数写在 run_go_individual_analysis.R，必须先 source 再调用 ★★★
# 旧做法 go_to_run <- c(...) 或直接 run_go_individual_analysis() 会报“没有这个函数”
src_fun <- "run_go_individual_analysis.R"
if (!file.exists(src_fun)) {
  stop("请把 run_go_individual_analysis.R 和主脚本放在同一目录（当前：", getwd(), "）")
}
source(src_fun, encoding = "UTF-8")
run_go_individual_analysis()
