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
#   8) 多个 GO：已写在 go_list；用 go_to_run 控制本次跑哪些，必须整段 for 循环，不要只跑一个文件夹
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

## 6. 逐个 GO 取基因（绝不合并） ------------------------------
# ★★★ 已修改：按 go_to_run 逐个取基因，不是合并成一个基因集 ★★★
message("本次将单独分析 ", length(go_to_run), " 个 GO：", paste(go_to_run, collapse = ", "))
go_gene_list <- lapply(go_to_run, get_go_genes)
names(go_gene_list) <- go_to_run

go_set_summary <- rbindlist(lapply(go_to_run, function(g) {
  dt <- go_gene_list[[g]]
  in_expr <- intersect(unique(dt$SYMBOL), rownames(expr))
  data.table(
    GO = g,
    GO_name = unname(go_name_map[g]),
    n_annotated = uniqueN(dt$SYMBOL),
    n_in_expression = length(in_expr),
    genes_in_expression = paste(in_expr, collapse = ";")
  )
}), fill = TRUE)
fwrite(go_set_summary, file.path(out_dir, "00_GO_gene_sets.csv"))

## 7. 每个 GO 单独分析 ----------------------------------------
all_scores <- list()
all_clin   <- list()
all_surv   <- list()
all_neg    <- list()
all_prot   <- list()

surv_endpoints <- list(
  OS  = c("OS",  "OS.time"),
  DSS = c("DSS", "DSS.time"),
  PFI = c("PFI", "PFI.time"),
  DFI = c("DFI", "DFI.time")
)

for (go_id in go_to_run) {
  go_title <- unname(go_name_map[go_id])
  message("\n========== ", go_id, " | ", go_title, " ==========")

  go_dir <- file.path(out_dir, "per_GO", paste0(safe_name(go_id), "_", safe_name(go_title)))
  dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)

  genes_here <- unique(go_gene_list[[go_id]]$SYMBOL)
  fwrite(
    data.table(GO = go_id, GO_name = go_title, gene = genes_here,
               in_expression = genes_here %in% rownames(expr)),
    file.path(go_dir, "genes.csv")
  )

  # ★★★ 已修改开始：变量名必须用 go_score，禁止再用 score ★★★
  # 旧代码（会报 closure 不可取子集，请勿运行）：
  #   score <- pathway_zmean(expr, genes_here)
  #   clin_use <- clinical_data[sample_std %in% names(score)]
  #   sc_clin <- score[clin_use$sample_std]
  go_score <- tryCatch(pathway_zmean(expr, genes_here), error = function(e) {
    message("  通路打分失败：", conditionMessage(e))
    NULL
  })
  if (!is.numeric(go_score) || length(go_score) == 0) {
    message("  可用基因 < ", min_pathway_genes, " 或打分失败，跳过该通路")
    next
  }
  n_used <- attr(go_score, "n_genes")
  used_genes <- attr(go_score, "genes")
  message("  通路基因用于打分：", n_used)
  fwrite(data.table(sample = names(go_score), pathway_score = as.numeric(go_score)),
         file.path(go_dir, "pathway_score.csv"))
  all_scores[[go_id]] <- go_score
  score <- go_score  # ★★★ 已修改：覆盖函数 score，避免后面旧代码 score[...] 报 closure

  # ---- 7.1 临床相关性（仅本通路分数） ----
  clin_use <- clinical_data[sample_std %in% names(go_score)]
  setkey(clin_use, sample_std)
  sc_clin <- go_score[clin_use$sample_std]
  # ★★★ 已修改结束 ★★★
  clin_rows <- lapply(clin_features, function(ft) {
    if (!ft %in% names(clin_use)) return(NULL)
    assoc_clinical_feature(sc_clin, clin_use[[ft]], ft)
  })
  clin_tab <- rbindlist(clin_rows, fill = TRUE)
  if (nrow(clin_tab) > 0) {
    clin_tab[, `:=`(GO = go_id, GO_name = go_title, fdr = p.adjust(pvalue, method = "BH"))]
    setcolorder(clin_tab, c("GO", "GO_name"))
    fwrite(clin_tab, file.path(go_dir, "clinical_association.csv"))
    all_clin[[go_id]] <- clin_tab

    # ★★★ 已修改开始：临床箱线图 ★★★
    # 旧输出 `null device 1` 只是 dev.off() 关闭 PDF，不是分析结果，也看不出画了几张图。
    sig_cat <- clin_tab[type == "categorical" & pvalue < 0.05][order(pvalue)]
    if (nrow(sig_cat) == 0) {
      message("  本通路无显著分类临床变量（p < 0.05），不画箱线图")
    } else {
      pdf_clin <- file.path(go_dir, "clinical_boxplots.pdf")
      n_box <- 0L
      pdf(pdf_clin, width = 7, height = 5)
      for (ft in head(sig_cat$feature, 6)) {
        plot_df <- data.frame(
          score = as.numeric(sc_clin),
          group = as.character(clin_use[[ft]]),
          stringsAsFactors = FALSE
        )
        plot_df <- plot_df[is.finite(plot_df$score) & !is.na(plot_df$group) & plot_df$group != "", ]
        tab <- table(plot_df$group)
        plot_df <- plot_df[plot_df$group %in% names(tab)[tab >= 3], ]
        if (length(unique(plot_df$group)) < 2) next
        ok <- tryCatch({
          p <- ggboxplot(plot_df, x = "group", y = "score", fill = "group", outlier.shape = NA) +
            stat_compare_means() +
            labs(title = paste0(go_id, "\n", go_title),
                 subtitle = paste0(ft, "  (p=", signif(sig_cat$pvalue[sig_cat$feature == ft][1], 3), ")"),
                 x = NULL, y = "Pathway score") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
          print(p)
          TRUE
        }, error = function(e) {
          message("  箱线图失败（", ft, "）：", conditionMessage(e))
          FALSE
        })
        if (isTRUE(ok)) n_box <- n_box + 1L
      }
      invisible(dev.off())
      if (n_box == 0L) {
        if (file.exists(pdf_clin)) file.remove(pdf_clin)
        message("  箱线图均被跳过（分组不足），未保存 PDF")
      } else {
        message("  已保存临床箱线图 ", n_box, " 张：", pdf_clin)
      }
    }
    # ★★★ 已修改结束 ★★★
  }

  # ---- 7.2 生存分析（仅本通路分数） ----
  # ★★★ 已修改开始：7.2 生存分析，score 改成 go_score ★★★
  # 旧代码（会报 closure 不可取子集，请勿运行）：
  #   surv_use <- survival_data[sample_std %in% names(score)]
  #   sc_surv <- score[surv_use$sample_std]
  if (!exists("go_score") || !is.numeric(go_score)) {
    stop("go_score 不是数值。请先运行通路打分段（搜：变量名必须用 go_score），不要运行 score[...]")
  }
  surv_use <- survival_data[sample_std %in% names(go_score)]
  sc_surv <- go_score[surv_use$sample_std]
  # ★★★ 已修改结束 ★★★
  surv_rows <- list()

  for (ep in names(surv_endpoints)) {
    ev_col <- surv_endpoints[[ep]][1]
    tm_col <- surv_endpoints[[ep]][2]
    if (!all(c(ev_col, tm_col) %in% names(surv_use))) next

    d <- data.frame(
      sample = surv_use$sample_std,
      time = as.numeric(surv_use[[tm_col]]),
      event = as.numeric(surv_use[[ev_col]]),
      score = as.numeric(sc_surv),
      stringsAsFactors = FALSE
    )
    d <- d[is.finite(d$time) & d$time > 0 & d$event %in% c(0, 1) & is.finite(d$score), ]
    if (nrow(d) < 20 || sum(d$event) < 5) next

    d$group <- ifelse(d$score >= median(d$score, na.rm = TRUE), "High", "Low")
    d$group <- factor(d$group, levels = c("Low", "High"))

    cox_cont <- tryCatch(coxph(Surv(time, event) ~ score, data = d), error = function(e) NULL)
    cox_grp  <- tryCatch(coxph(Surv(time, event) ~ group, data = d), error = function(e) NULL)
    fit_km   <- tryCatch(survfit(Surv(time, event) ~ group, data = d), error = function(e) NULL)

    if (!is.null(cox_cont)) {
      s1 <- summary(cox_cont)
      surv_rows[[paste0(ep, "_cont")]] <- data.table(
        GO = go_id, GO_name = go_title, endpoint = ep, model = "continuous_score",
        n = nrow(d), events = sum(d$event),
        HR = s1$conf.int[1, 1],
        HR_low = s1$conf.int[1, 3],
        HR_high = s1$conf.int[1, 4],
        pvalue = s1$coefficients[1, "Pr(>|z|)"]
      )
    }
    if (!is.null(cox_grp)) {
      s2 <- summary(cox_grp)
      surv_rows[[paste0(ep, "_grp")]] <- data.table(
        GO = go_id, GO_name = go_title, endpoint = ep, model = "High_vs_Low",
        n = nrow(d), events = sum(d$event),
        HR = s2$conf.int[1, 1],
        HR_low = s2$conf.int[1, 3],
        HR_high = s2$conf.int[1, 4],
        pvalue = s2$coefficients[1, "Pr(>|z|)"]
      )
    }

    if (!is.null(fit_km) && ep == "OS") {
      tryCatch({
        d$time_month <- d$time / 30.44
        fit_m <- survfit(Surv(time_month, event) ~ group, data = d)
        pdf(file.path(go_dir, "KM_OS.pdf"), width = 7, height = 6)
        print(ggsurvplot(
          fit_m, data = d, pval = TRUE, risk.table = TRUE,
          legend.title = "Pathway score",
          legend.labs = c("Low", "High"),
          xlab = "Time (months)", ylab = "Overall survival",
          title = paste0(go_id, "\n", go_title),
          ggtheme = theme_bw()
        ))
        dev.off()
      }, error = function(e) {
        if (dev.cur() > 1) dev.off()
        message("  KM 作图失败：", conditionMessage(e))
      })
    }
  }

  surv_tab <- rbindlist(surv_rows, fill = TRUE)
  if (nrow(surv_tab) > 0) {
    fwrite(surv_tab, file.path(go_dir, "survival_cox.csv"))
    all_surv[[go_id]] <- surv_tab
  }

  # ---- 7.3 全基因组负相关基因（仅对本通路分数） ----
  cor_tab <- spearman_vs_score(expr, go_score)  # ★★★ 已修改：传入 go_score
  if (nrow(cor_tab) > 0) {
    setnames(cor_tab, "feature", "gene")
    cor_tab[, `:=`(
      GO = go_id,
      GO_name = go_title,
      in_this_GO = gene %in% used_genes
    )]
    fwrite(cor_tab[order(-spearman_r)], file.path(go_dir, "all_gene_correlation.csv"))

    neg_tab <- cor_tab[spearman_r < neg_r_cutoff & fdr < neg_fdr_cutoff]
    setorder(neg_tab, spearman_r)
    fwrite(neg_tab, file.path(go_dir, "negative_correlated_genes.csv"))
    fwrite(neg_tab[spearman_r <= strict_r_cutoff],
           file.path(go_dir, "negative_correlated_genes_strict_r0.3.csv"))
    all_neg[[go_id]] <- neg_tab
    message("  显著负相关基因：", nrow(neg_tab),
            " ；|r|>=0.3：", nrow(neg_tab[spearman_r <= strict_r_cutoff]))

    vol_df <- as.data.frame(cor_tab)
    rownames(vol_df) <- vol_df$gene
    key_neg <- head(neg_tab$gene, 15)
    tryCatch({
      pdf(file.path(go_dir, "volcano_gene_correlation.pdf"), width = 9, height = 7)
      print(EnhancedVolcano(
        vol_df,
        lab = vol_df$gene,
        x = "spearman_r",
        y = "fdr",
        xlab = "Spearman r (gene vs this GO score)",
        ylab = "-Log10 FDR",
        title = paste(go_id, go_title),
        pCutoff = neg_fdr_cutoff,
        FCcutoff = 0.3,
        xlim = c(-1, 1),
        selectLab = key_neg,
        drawConnectors = TRUE
      ))
      dev.off()
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      message("  火山图失败：", conditionMessage(e))
    })

    if (nrow(neg_tab) > 0) {
      topn <- head(neg_tab, 30)
      pbar <- ggplot(topn, aes(x = reorder(gene, spearman_r), y = spearman_r)) +
        geom_col(fill = "#3C5488") +
        coord_flip() +
        labs(title = paste0("Top negative genes: ", go_id),
             x = NULL, y = "Spearman r") +
        theme_bw()
      ggsave(file.path(go_dir, "top_negative_genes.pdf"), pbar, width = 7, height = 8)
    }
  }

  # ---- 7.4 蛋白水平相关性（补充，仍按本通路分数） ----
  if (!is.null(prot_mat)) {
    prot_common <- intersect(colnames(prot_mat), names(go_score))  # ★★★ 已修改：go_score
    if (length(prot_common) >= 20) {
      prot_cor <- spearman_vs_score(prot_mat, go_score)
      if (nrow(prot_cor) > 0) {
        setnames(prot_cor, "feature", "protein")
        prot_cor[, `:=`(GO = go_id, GO_name = go_title)]
        fwrite(prot_cor[order(spearman_r)], file.path(go_dir, "protein_correlation.csv"))
        all_prot[[go_id]] <- prot_cor[spearman_r < 0 & fdr < neg_fdr_cutoff][order(spearman_r)]
      }
    }
  }
}

## 8. 汇总输出（仍是“每个 GO 一行/一堆结果”，不是合并通路） ----
if (length(all_scores) > 0) {
  score_mat <- do.call(cbind, lapply(all_scores, function(x) {
    x[colnames(expr)]
  }))
  colnames(score_mat) <- names(all_scores)
  rownames(score_mat) <- colnames(expr)
  fwrite(data.table(sample = rownames(score_mat), as.data.table(score_mat)),
         file.path(out_dir, "01_pathway_scores_each_GO.csv"))

  # 各通路分数热图（每列仍是单独 GO 分数，不是合并基因集）
  anno_df <- data.frame(row.names = rownames(score_mat))
  clin_hm <- clinical_data[sample_std %in% rownames(score_mat)]
  if (nrow(clin_hm) > 0) {
    if ("stage_simplified" %in% names(clin_hm)) {
      anno_df$Stage <- clin_hm$stage_simplified[match(rownames(anno_df), clin_hm$sample_std)]
    }
    er_col <- first_present(names(clin_hm), c("er_status_by_ihc", "ER.Status"))
    if (!is.na(er_col)) {
      anno_df$ER <- as.character(clin_hm[[er_col]][match(rownames(anno_df), clin_hm$sample_std)])
    }
  }
  ha <- NULL
  if (ncol(anno_df) > 0) ha <- ComplexHeatmap::HeatmapAnnotation(df = anno_df)
  # ★★★ 已修改开始：热图缩放，禁止 t(scale(score_mat)) ★★★
  # 旧代码（会报 scale 找不到 matrix 方法，请勿运行）：
  #   z_score <- t(scale(score_mat))
  sm <- as.matrix(score_mat)
  storage.mode(sm) <- "double"
  z_score <- matrix(
    NA_real_,
    nrow = ncol(sm),
    ncol = nrow(sm),
    dimnames = list(colnames(sm), rownames(sm))
  )
  for (j in seq_len(ncol(sm))) {
    x <- sm[, j]
    sdx <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(sdx) || sdx < 1e-12) sdx <- 1
    z_score[j, ] <- (x - mean(x, na.rm = TRUE)) / sdx
  }
  z_score[!is.finite(z_score)] <- 0
  # ★★★ 已修改结束 ★★★
  # ★★★ 已修改开始：热图必须用 ComplexHeatmap::，不要直接 Heatmap()/draw() ★★★
  # 旧代码（heatmaps 包会抢走 Heatmap，报“参数没有用”）：
  #   draw(Heatmap(z_score, name = "z-score", ...))
  tryCatch({
    pdf_hm <- file.path(out_dir, "01_pathway_score_heatmap.pdf")
    pdf(pdf_hm, width = 12, height = 6)
    ht <- ComplexHeatmap::Heatmap(
      z_score,
      name = "z-score",
      top_annotation = ha,
      show_column_names = FALSE,
      row_names_gp = grid::gpar(fontsize = 8),
      column_title = "Each GO pathway score (not pooled)"
    )
    ComplexHeatmap::draw(ht)
    invisible(dev.off())
    message("  已保存通路分数热图：", pdf_hm)
  }, error = function(e) {
    if (dev.cur() > 1) invisible(dev.off())
    message("热图绘制失败：", conditionMessage(e))
  })
  # ★★★ 已修改结束 ★★★
}

if (length(all_clin) > 0) {
  clin_all <- rbindlist(all_clin, fill = TRUE)
  fwrite(clin_all, file.path(out_dir, "02_clinical_association_each_GO.csv"))
}

if (length(all_surv) > 0) {
  surv_all <- rbindlist(all_surv, fill = TRUE)
  fwrite(surv_all, file.path(out_dir, "03_survival_cox_each_GO.csv"))

  os_hl <- surv_all[endpoint == "OS" & model == "High_vs_Low"]
  if (nrow(os_hl) > 0) {
    os_hl[, lab := paste(GO, GO_name)]
    os_hl[, lab := factor(lab, levels = rev(lab))]
    # ★★★ 已修改开始：整段森林图请整块粘贴，geom_errorbar 行末尾必须有 + ★★★
    pfor <- ggplot(os_hl, aes(x = HR, y = lab)) +
      geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
      geom_errorbar(aes(xmin = HR_low, xmax = HR_high), orientation = "y", width = 0.2) +
      geom_point(aes(color = pvalue < 0.05), size = 3) +
      scale_x_log10() +
      labs(title = "OS Cox: High vs Low (each GO separately)",
           x = "Hazard ratio", y = NULL, color = "p < 0.05") +
      theme_bw()
    ggsave(file.path(out_dir, "03_OS_forest_each_GO.pdf"), pfor, width = 10, height = 6)
    message("  已保存 OS 森林图：", file.path(out_dir, "03_OS_forest_each_GO.pdf"))
    # ★★★ 已修改结束 ★★★
  }
}

if (length(all_neg) > 0) {
  neg_all <- rbindlist(all_neg, fill = TRUE)
  fwrite(neg_all, file.path(out_dir, "04_negative_genes_each_GO.csv"))
  neg_count <- neg_all[, .(n_negative_genes = .N,
                           n_strict_r0.3 = sum(spearman_r <= strict_r_cutoff)),
                       by = .(GO, GO_name)]
  fwrite(neg_count, file.path(out_dir, "04_negative_genes_count_each_GO.csv"))
}

if (length(all_prot) > 0) {
  fwrite(rbindlist(all_prot, fill = TRUE),
         file.path(out_dir, "05_negative_proteins_each_GO.csv"))
}

fwrite(go_set_summary, file.path(out_dir, "00_GO_gene_sets.csv"))
message("\n分析完成。每个 GO 的独立结果在：", normalizePath(out_dir, mustWork = FALSE))
message("请重点查看 per_GO/ 下各通路文件夹，以及 04_negative_genes_each_GO.csv")
