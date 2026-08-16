################################################################################
# TCGA乳腺癌（BRCA）多组学数据分析完整流程
# 包括：RNA-seq (FPKM)、临床数据、生存数据、蛋白表达数据
# 目标：筛选与乳腺癌进展相关的 GO 通路基因
#
# 核心原则：每个 GO ID 单独取基因、单独打分、单独统计，绝不把通路合并
#
# ============================ 请按这个跑 ============================
# 1. 用本文件完整覆盖：E:/R/BRCA/GO_pathway_individual_analysis.R
# 2. RStudio：Session -> Restart R（清掉旧的空函数 / score 函数冲突）
# 3. 打开本文件，点 Source（整份运行）。不要从中间逐段粘贴。
# 4. 不要在读数据之前调用 run_go_individual_analysis()
# 5. 不要 library(heatmaps)；热图只用 ggplot，不要 Heatmap()/draw()
# 6. 向量用小写 c()，不要写成 C()（C 是 contrasts）
# 7. 函数定义之后才会分析；本脚本末尾会自动调用一次
# 8. 汇总后会画两张气泡图：预后（OS/DSS）、转移（PFI/DFI + M/N/分期）
# 9. 若 17 个 GO 已经跑完、只想补画气泡图：把下面 run_mode 改成 "bubbles_only" 再 Source
#    或在控制台运行：draw_prognosis_metastasis_bubbles()
# ====================================================================
################################################################################

# ==============================================================================
# 第一部分：加载所需 R 包
# ==============================================================================
# 数据处理
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(stringr)
library(purrr)

# 基因注释（每个 GO 单独取基因）
library(org.Hs.eg.db)
library(AnnotationDbi)

# 生存分析
library(survival)
library(survminer)

# 可视化（不要加载 heatmaps：会盖掉 ComplexHeatmap::Heatmap）
library(ggplot2)
library(ggpubr)
library(EnhancedVolcano)

# 不再加载 SummarizedExperiment / TCGAbiolinks / edgeR / limma / clusterProfiler /
# enrichplot / ComplexHeatmap / heatmaps：本流程读的是 Xena tsv，不需要它们，
# 而且 SummarizedExperiment 会把 scale()/t() 变成 S4 泛型，普通 matrix 会报错。

# ==============================================================================
# 第二部分：可调参数、工作目录与数据读取
# ==============================================================================
# 探索性阈值已放宽：乳腺肿瘤里神经/轴突 GO 基因少、死亡事件少。
# 原来 FDR<0.05 且 |r|>=0.3、生存至少 20 人/5 个事件，很多通路会空。
work_dir <- "E:/R/BRCA"   # 与当前工程一致；若不存在则使用当前目录
use_primary_tumor_only <- TRUE
min_pathway_genes <- 1     # 有 1 个基因也计算通路分数（原 2）
neg_pvalue_cutoff <- 0.05  # 负相关主列表：名义 p（不再先卡 FDR）
neg_fdr_cutoff <- 0.25     # 另存一份较宽的 FDR 列表
neg_r_cutoff <- 0          # 负相关：Spearman r < 0
strict_r_cutoff <- -0.15   # 高置信负相关（原 -0.3 对 BRCA 过严）
min_surv_n <- 10           # 生存分析最少样本（原 20）
min_surv_events <- 3       # 最少事件数（原 5）
min_clin_n <- 6            # 临床关联最少样本（原 10）
min_group_n <- 2           # 每组最少人数（原 3）
clin_plot_p_cutoff <- 0.10 # 箱线图：p<0.10 也画（原 0.05）
out_dir <- "results_GO_individual"
# "full" = 读数据并做完全部分析（含气泡图）
# "bubbles_only" = 不重跑全基因组相关，只用已有 results_GO_individual 画两张气泡图
run_mode <- "full"

# 指定的 GO 通路（每个都单独分析；必须用小写 c()）
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

# 默认分析 go_list 里全部 17 个通路。若只想补跑某几个，改成：
#   go_to_run <- c("GO:0007411", "GO:0007409")
go_to_run <- go_list

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
if (!identical(run_mode, "bubbles_only")) {
  miss <- need_files[!file.exists(need_files)]
  if (length(miss) > 0) {
    stop("以下文件不在工作目录中：\n  ", paste(miss, collapse = "\n  "))
  }

  # 读取数据（不要在这之前调用 run_go_individual_analysis）
  fpkm_data     <- fread("TCGA-BRCA.star_fpkm.tsv")
  probe_annot   <- fread("gencode.v36.annotation.gtf.gene.probemap")
  clinical_data <- fread("TCGA-BRCA.clinical.tsv")
  survival_data <- fread("TCGA-BRCA.survival.tsv")
  protein_data  <- if (file.exists("TCGA-BRCA.protein.tsv")) fread("TCGA-BRCA.protein.tsv") else NULL
  if (is.null(protein_data)) message("未找到 TCGA-BRCA.protein.tsv，将跳过蛋白相关性")
} else {
  message("run_mode = bubbles_only：跳过 FPKM 读取，只用已有结果画气泡图")
  if (file.exists("TCGA-BRCA.clinical.tsv")) {
    clinical_data <- fread("TCGA-BRCA.clinical.tsv")
  }
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "per_GO"), showWarnings = FALSE)

# ==============================================================================
# 第三部分：工具函数
# ==============================================================================
# 标准化条形码：转大写、点号变横杠、去掉末尾 A、截取前 15 位
normalize_barcode <- function(x) {
  x <- toupper(gsub("\\.", "-", as.character(x)))
  x <- sub("A$", "", x)
  ifelse(nchar(x) >= 15, substr(x, 1, 15), x)
}

# 从标准化条形码提取第 14–15 位作为样本类型编码（01=原发肿瘤）
sample_type_code <- function(ids) {
  substr(normalize_barcode(ids), 14, 15)
}

# 只保留字母数字，其余换成下划线（用于文件夹名）
safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# 按优先级查找第一个存在的列名
first_present <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) NA_character_ else hit[1]
}

# 根据单个 GO 编号，从 org.Hs.eg.db 取基因（SYMBOL / ENSEMBL / ENTREZID）
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
  extra <- tryCatch({
    go2eg <- as.list(org.Hs.eg.db::org.Hs.egGO2ALLEGS)
    eg <- unique(as.character(go2eg[[go_id]]))
    if (length(eg) == 0) {
      NULL
    } else {
      AnnotationDbi::select(
        org.Hs.eg.db,
        keys = eg,
        keytype = "ENTREZID",
        columns = c("SYMBOL", "ENSEMBL")
      )
    }
  }, error = function(e) NULL)
  if (!is.null(extra) && nrow(extra) > 0) {
    extra <- as.data.table(extra)
    extra[, GO := go_id]
    res <- if (is.null(res) || nrow(res) == 0) extra else rbind(as.data.table(res), extra, fill = TRUE)
  }
  if (is.null(res) || nrow(res) == 0) {
    return(data.table(GO = go_id, SYMBOL = character(), ENSEMBL = character(), ENTREZID = character()))
  }
  res <- as.data.table(res)
  if ("GOALL" %in% names(res)) setnames(res, "GOALL", "GO", skip_absent = TRUE)
  if (!"GO" %in% names(res)) res[, GO := go_id]
  res <- unique(res[!is.na(SYMBOL) & SYMBOL != "", .(GO, SYMBOL, ENSEMBL, ENTREZID)])
  res
}

# 通路活性：基因 z-score 后对样本取均值。不要用 scale()/t()（S4 泛型会报错）
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

# 每个基因/蛋白与通路分数的 Spearman 相关
spearman_vs_score <- function(mat, score_vec) {
  common <- intersect(colnames(mat), names(score_vec))
  if (length(common) < 5) return(data.table())
  mat <- mat[, common, drop = FALSE]
  score_vec <- score_vec[common]
  keep <- apply(mat, 1, function(x) sd(x, na.rm = TRUE) > 0)
  mat <- mat[keep, , drop = FALSE]
  n <- ncol(mat)
  r <- as.numeric(cor(
    base::t(as.matrix(mat)),
    score_vec,
    method = "spearman",
    use = "pairwise.complete.obs"
  ))
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

# 通路分数 vs 临床特征（连续=Spearman，两组=Wilcoxon，多组=Kruskal-Wallis）
assoc_clinical_feature <- function(score_vec, feature, feature_name) {
  df <- data.frame(score = as.numeric(score_vec), feature = feature, stringsAsFactors = FALSE)
  df <- df[is.finite(df$score) & !is.na(df$feature), , drop = FALSE]
  df$feature <- as.character(df$feature)
  df <- df[df$feature != "" & !tolower(df$feature) %in% c("na", "nan", "not reported", "unknown", "[not available]", "[unknown]"), ]
  if (nrow(df) < min_clin_n) return(NULL)

  feat_num <- suppressWarnings(as.numeric(df$feature))
  if (mean(is.finite(feat_num)) >= 0.8) {
    df$feature <- feat_num
    df <- df[is.finite(df$feature), ]
    if (nrow(df) < min_clin_n || sd(df$feature) == 0) return(NULL)
    ct <- suppressWarnings(cor.test(df$score, df$feature, method = "spearman", exact = FALSE))
    return(data.table(
      feature = feature_name, type = "continuous", n = nrow(df),
      stat = unname(ct$estimate), pvalue = ct$p.value,
      method = "Spearman", detail = sprintf("rho=%.3f", unname(ct$estimate))
    ))
  }

  df$feature <- droplevels(factor(df$feature))
  tab <- table(df$feature)
  df <- df[df$feature %in% names(tab)[tab >= min_group_n], ]
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

# 分期文本标准化为 Stage I–IV
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

# 远处转移 M0 / M1（MX、未知记为 NA）
classify_m <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("\\bM1\\b|M1[ABC]?", x)] <- "M1"
  out[is.na(out) & grepl("\\bM0\\b", x)] <- "M0"
  out[grepl("MX|NOT AVAILABLE|UNKNOWN|NOT REPORTED", x)] <- NA_character_
  out
}

# 淋巴结 N0 / Nplus
classify_n <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("\\bN[1-3]", x)] <- "Nplus"
  out[is.na(out) & grepl("\\bN0\\b", x)] <- "N0"
  out[grepl("NX|NOT AVAILABLE|UNKNOWN|NOT REPORTED", x)] <- NA_character_
  out
}

classify_early_late <- function(stage_simplified) {
  x <- as.character(stage_simplified)
  out <- rep(NA_character_, length(x))
  out[x %in% c("Stage I", "Stage II")] <- "Early_I_II"
  out[x %in% c("Stage III", "Stage IV")] <- "Late_III_IV"
  out
}

# 在临床表上补转移相关二分类列（可重复调用）
ensure_metastasis_clin_columns <- function(clin) {
  if (is.null(clin) || nrow(clin) == 0) return(clin)
  clin <- copy(clin)
  if (!"sample_std" %in% names(clin)) {
    cid <- detect_id_col(clin, c("sampleID", "sample", "bcr_patient_barcode", "submitter_id"))
    clin[, sample_std := normalize_barcode(get(cid))]
    clin <- clin[!duplicated(sample_std)]
  }
  if (!"stage_simplified" %in% names(clin)) {
    stage_col <- first_present(names(clin), c(
      "ajcc_pathologic_tumor_stage", "pathologic_stage", "clinical_stage", "ajcc_pathologic_stage"
    ))
    if (!is.na(stage_col)) clin[, stage_simplified := simplify_stage(get(stage_col))]
  }
  m_col <- first_present(names(clin), c(
    "ajcc_metastasis_pathologic_pm", "pathologic_M", "pathologic_m",
    "ajcc_pathologic_m", "M", "metastasis"
  ))
  n_col <- first_present(names(clin), c(
    "ajcc_nodes_pathologic_pn", "pathologic_N", "pathologic_n",
    "ajcc_pathologic_n", "N"
  ))
  if (!is.na(m_col)) clin[, meta_M := classify_m(get(m_col))] else if (!"meta_M" %in% names(clin)) clin[, meta_M := NA_character_]
  if (!is.na(n_col)) clin[, meta_N := classify_n(get(n_col))] else if (!"meta_N" %in% names(clin)) clin[, meta_N := NA_character_]
  if ("stage_simplified" %in% names(clin)) {
    clin[, meta_stage := classify_early_late(stage_simplified)]
  } else if (!"meta_stage" %in% names(clin)) {
    clin[, meta_stage := NA_character_]
  }
  clin
}

# 两组 Wilcoxon：pos 为“更差/转移相关”组，effect = median(pos) - median(neg)
two_group_score_assoc <- function(score_vec, group, pos_level, neg_level, feature_name) {
  df <- data.frame(
    score = as.numeric(score_vec),
    group = as.character(group),
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$score) & df$group %in% c(pos_level, neg_level), , drop = FALSE]
  n_pos <- sum(df$group == pos_level)
  n_neg <- sum(df$group == neg_level)
  if (n_pos < min_group_n || n_neg < min_group_n) return(NULL)
  wt <- suppressWarnings(wilcox.test(score ~ group, data = df))
  med_pos <- stats::median(df$score[df$group == pos_level], na.rm = TRUE)
  med_neg <- stats::median(df$score[df$group == neg_level], na.rm = TRUE)
  data.table(
    feature = feature_name,
    n = nrow(df),
    n_pos = n_pos,
    n_neg = n_neg,
    effect = med_pos - med_neg,
    pvalue = wt$p.value,
    method = "Wilcoxon",
    effect_type = "delta_median_score",
    detail = sprintf(
      "%s median=%.3f (n=%d); %s median=%.3f (n=%d)",
      pos_level, med_pos, n_pos, neg_level, med_neg, n_neg
    )
  )
}

collect_metastasis_clinical <- function(go_id, go_title, score_vec, clin_use) {
  rows <- list()
  if ("meta_M" %in% names(clin_use)) {
    rows[["M"]] <- two_group_score_assoc(score_vec, clin_use$meta_M, "M1", "M0", "Distant_M_M1_vs_M0")
  }
  if ("meta_N" %in% names(clin_use)) {
    rows[["N"]] <- two_group_score_assoc(score_vec, clin_use$meta_N, "Nplus", "N0", "Node_Nplus_vs_N0")
  }
  if ("meta_stage" %in% names(clin_use)) {
    rows[["S"]] <- two_group_score_assoc(
      score_vec, clin_use$meta_stage, "Late_III_IV", "Early_I_II", "Stage_Late_vs_Early"
    )
  }
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(data.table())
  dt <- rbindlist(rows, fill = TRUE)
  if (is.null(dt) || nrow(dt) == 0) return(data.table())
  dt[, `:=`(GO = go_id, GO_name = go_title)]
  dt
}

# 气泡图：x=指标，y=通路，点大小=-log10(p)，颜色=效应方向
make_go_bubble_plot <- function(plot_dt, title, subtitle, caption, color_name) {
  d <- copy(as.data.table(plot_dt))
  d <- d[is.finite(pvalue) & is.finite(effect)]
  if (nrow(d) == 0) return(NULL)
  d[, neglogp := pmin(10, -log10(pmax(pvalue, 1e-12)))]
  d[, go_lab := paste0(GO, "  ", GO_name)]
  go_rank <- d[, .(m = max(neglogp, na.rm = TRUE)), by = go_lab]
  setorder(go_rank, m)
  d[, go_lab := factor(go_lab, levels = go_rank$go_lab)]
  if (!is.factor(d$feature_lab)) {
    d[, feature_lab := factor(feature_lab, levels = unique(as.character(feature_lab)))]
  }
  lim <- max(0.4, stats::quantile(abs(d$effect), 0.98, na.rm = TRUE), na.rm = TRUE)
  ggplot(d, aes(x = feature_lab, y = go_lab, size = neglogp, color = effect)) +
    geom_point(alpha = 0.9) +
    scale_size_continuous(range = c(2.5, 12), name = expression(-log[10](p))) +
    scale_color_gradient2(
      low = "#3C5488", mid = "grey90", high = "#E64B35",
      midpoint = 0, limits = c(-lim, lim), name = color_name
    ) +
    labs(title = title, subtitle = subtitle, caption = caption, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.major.x = element_blank(),
      plot.caption = element_text(hjust = 0, size = 8, color = "grey30")
    )
}

draw_prognosis_metastasis_bubbles <- function(result_dir = NULL, surv_all = NULL, all_scores = NULL) {
  if (is.null(result_dir)) {
    result_dir <- if (exists("out_dir", inherits = TRUE)) {
      normalizePath(out_dir, mustWork = FALSE)
    } else {
      "results_GO_individual"
    }
  }
  dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(surv_all)) {
    surv_file <- file.path(result_dir, "03_survival_cox_each_GO.csv")
    if (!file.exists(surv_file)) {
      stop("找不到 ", surv_file, " 。请先跑完整分析，或等当前 17 个 GO 跑完再画气泡图。")
    }
    surv_all <- fread(surv_file)
  }
  surv_all <- as.data.table(surv_all)
  if (!"model" %in% names(surv_all)) stop("生存结果表缺少 model 列")
  hl <- surv_all[surv_all[["model"]] == "High_vs_Low"]
  if (nrow(hl) == 0) hl <- surv_all[surv_all[["model"]] == "continuous_score"]
  hl <- hl[is.finite(HR) & HR > 0]
  hl[, effect := log2(HR)]
  hl[, effect_type := "log2_HR"]
  hl[, feature_lab := as.character(endpoint)]

  prog_map <- c(
    OS = "OS 总生存",
    DSS = "DSS 疾病特异生存"
  )
  meta_surv_map <- c(
    PFI = "PFI 无进展间隔",
    DFI = "DFI 无病间隔"
  )
  prog <- hl[endpoint %in% names(prog_map)]
  if (nrow(prog) > 0) prog[, feature_lab := unname(prog_map[as.character(endpoint)])]
  meta_s <- hl[endpoint %in% names(meta_surv_map)]
  if (nrow(meta_s) > 0) meta_s[, feature_lab := unname(meta_surv_map[as.character(endpoint)])]

  # 转移临床：M1、N+、III-IV；优先用内存中的通路分数，否则读 CSV
  meta_c <- data.table()
  clin_use <- NULL
  if (exists("clinical_data", inherits = TRUE)) {
    clin_use <- ensure_metastasis_clin_columns(get("clinical_data", inherits = TRUE))
  } else if (file.exists("TCGA-BRCA.clinical.tsv")) {
    tmp_clin <- fread("TCGA-BRCA.clinical.tsv")
    cid <- detect_id_col(tmp_clin, c("sampleID", "sample", "bcr_patient_barcode", "submitter_id"))
    tmp_clin[, sample_std := normalize_barcode(get(cid))]
    tmp_clin <- tmp_clin[!duplicated(sample_std)]
    clin_use <- ensure_metastasis_clin_columns(tmp_clin)
  }

  score_list <- all_scores
  if (is.null(score_list) || length(score_list) == 0) {
    score_file <- file.path(result_dir, "01_pathway_scores_each_GO.csv")
    if (file.exists(score_file)) {
      sm <- fread(score_file)
      go_cols <- setdiff(names(sm), "sample")
      score_list <- lapply(go_cols, function(g) {
        v <- as.numeric(sm[[g]])
        names(v) <- as.character(sm$sample)
        v
      })
      names(score_list) <- go_cols
    }
  }

  if (!is.null(clin_use) && length(score_list) > 0) {
    meta_c <- rbindlist(lapply(names(score_list), function(g) {
      sc <- score_list[[g]]
      title <- if (exists("go_name_map", inherits = TRUE) && g %in% names(go_name_map)) {
        unname(go_name_map[g])
      } else {
        g
      }
      common <- intersect(clin_use$sample_std, names(sc))
      if (length(common) < 10) return(NULL)
      cu <- clin_use[sample_std %in% common]
      collect_metastasis_clinical(g, title, sc[cu$sample_std], cu)
    }), fill = TRUE)
    if (is.null(meta_c)) meta_c <- data.table()
  }

  meta_lab_map <- c(
    Distant_M_M1_vs_M0 = "远处转移 M1 vs M0",
    Node_Nplus_vs_N0 = "淋巴结 N+ vs N0",
    Stage_Late_vs_Early = "分期 III-IV vs I-II"
  )
  if (nrow(meta_c) > 0) {
    meta_c[, feature_lab := unname(meta_lab_map[as.character(feature)])]
    meta_c[is.na(feature_lab), feature_lab := feature]
  }

  meta_all <- rbindlist(list(
    if (nrow(meta_s) > 0) meta_s[, .(GO, GO_name, feature, feature_lab, effect, pvalue, n, effect_type, HR, HR_low, HR_high)] else NULL,
    if (nrow(meta_c) > 0) meta_c[, .(GO, GO_name, feature, feature_lab, effect, pvalue, n, effect_type, HR = NA_real_, HR_low = NA_real_, HR_high = NA_real_)] else NULL
  ), fill = TRUE)

  n_go <- uniqueN(c(prog$GO, meta_all$GO))
  ht <- max(6.5, min(14, 0.38 * max(1, n_go) + 3.2))

  if (nrow(prog) > 0) {
    prog[, feature_lab := factor(feature_lab, levels = unname(prog_map))]
    fwrite(prog, file.path(result_dir, "06_prognosis_bubble_data.csv"))
    p_prog <- make_go_bubble_plot(
      prog,
      title = "各 GO 通路与乳腺癌预后",
      subtitle = "每个 GO 单独打分；High vs Low Cox（未合并通路）",
      caption = "颜色：log2(HR)。红=通路高分组预后更差（HR>1）；蓝=高分组预后更好（HR<1）。点越大 p 越小。",
      color_name = "log2(HR)"
    )
    if (!is.null(p_prog)) {
      pdf_p <- file.path(result_dir, "06_prognosis_bubble_each_GO.pdf")
      ggsave(pdf_p, p_prog, width = 9.5, height = ht)
      ggsave(file.path(result_dir, "06_prognosis_bubble_each_GO.png"), p_prog, width = 9.5, height = ht, dpi = 150)
      message("  已保存预后气泡图：", pdf_p)
    }
  } else {
    message("  无 OS/DSS High vs Low 结果，跳过预后气泡图")
  }

  if (nrow(meta_all) > 0) {
    lv <- c(unname(meta_surv_map), unname(meta_lab_map))
    lv <- lv[lv %in% as.character(meta_all$feature_lab)]
    meta_all[, feature_lab := factor(as.character(feature_lab), levels = lv)]
    fwrite(meta_all, file.path(result_dir, "07_metastasis_bubble_data.csv"))
    p_meta <- make_go_bubble_plot(
      meta_all,
      title = "各 GO 通路与乳腺癌转移",
      subtitle = "PFI/DFI：High vs Low Cox；M/N/分期：转移或晚分期组 vs 对照的通路分数差",
      caption = paste0(
        "颜色：生存为 log2(HR)，临床为 Δscore（转移/晚分期组中位数 − 对照）。",
        "红=高通路分数偏向更差转移表型；蓝相反。点越大 p 越小。",
        "TCGA-BRCA 的 M1 样本通常很少，M 列可能缺失或很弱。"
      ),
      color_name = "effect"
    )
    if (!is.null(p_meta)) {
      pdf_m <- file.path(result_dir, "07_metastasis_bubble_each_GO.pdf")
      ggsave(pdf_m, p_meta, width = 10.5, height = ht)
      ggsave(file.path(result_dir, "07_metastasis_bubble_each_GO.png"), p_meta, width = 10.5, height = ht, dpi = 150)
      message("  已保存转移气泡图：", pdf_m)
    }
  } else {
    message("  无转移相关结果，跳过转移气泡图")
  }
  invisible(list(prognosis = prog, metastasis = meta_all))
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

# ==============================================================================
# 第四部分：表达矩阵预处理
# ==============================================================================
if (!identical(run_mode, "bubbles_only")) {
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

# ==============================================================================
# 第五部分：临床 / 生存 / 蛋白对齐
# ==============================================================================
clin_id <- detect_id_col(clinical_data, c("sampleID", "sample", "bcr_patient_barcode", "submitter_id"))
surv_id <- detect_id_col(survival_data, c("sample", "sampleID", "bcr_patient_barcode"))

clinical_data <- copy(clinical_data)
survival_data <- copy(survival_data)
clinical_data[, sample_std := normalize_barcode(get(clin_id))]
survival_data[, sample_std := normalize_barcode(get(surv_id))]
clinical_data <- clinical_data[!duplicated(sample_std)]
survival_data <- survival_data[!duplicated(sample_std)]

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

stage_col <- first_present(names(clinical_data), c(
  "ajcc_pathologic_tumor_stage", "pathologic_stage", "clinical_stage", "ajcc_pathologic_stage"
))
if (!is.na(stage_col)) {
  clinical_data[, stage_simplified := simplify_stage(get(stage_col))]
  clin_features <- unique(c("stage_simplified", clin_features))
}
clinical_data <- ensure_metastasis_clin_columns(clinical_data)
message(
  "转移分组样本数：M1=", sum(clinical_data$meta_M == "M1", na.rm = TRUE),
  "  N+=", sum(clinical_data$meta_N == "Nplus", na.rm = TRUE),
  "  III-IV=", sum(clinical_data$meta_stage == "Late_III_IV", na.rm = TRUE)
)

prot_mat <- NULL
if (!is.null(protein_data) && nrow(protein_data) > 0) {
  prot_mat <- tryCatch({
    prot_c1 <- as.character(protein_data[[1]])
    if (mean(grepl("^TCGA", prot_c1, ignore.case = TRUE), na.rm = TRUE) > 0.5) {
      tmp <- as_sample_matrix(protein_data, names(protein_data)[1])
      tmp <- base::t(as.matrix(tmp))
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
}  # 结束 if (!identical(run_mode, "bubbles_only")) 的数据预处理

# ==============================================================================
# 第六–八部分：每个 GO 单独分析（全部包在函数里，末尾调用一次）
# ==============================================================================
run_go_individual_analysis <- function(go_ids = NULL) {
  if (is.null(go_ids)) {
    go_ids <- get("go_to_run", envir = .GlobalEnv, inherits = FALSE)
  }
  go_to_run <- unique(as.character(go_ids))
  if (length(go_to_run) == 0) stop("go_to_run 是空的，没有要分析的 GO")
  if (!exists("expr", envir = .GlobalEnv, inherits = FALSE) || is.null(get("expr", .GlobalEnv))) {
    stop("还没有表达矩阵 expr。请从脚本开头 Source 整份，不要在读数据之前调用本函数。")
  }

  out_dir_abs <- normalizePath(out_dir, mustWork = FALSE)
  dir.create(file.path(out_dir_abs, "per_GO"), recursive = TRUE, showWarnings = FALSE)
  message("输出根目录：", out_dir_abs)
  fwrite(data.table(
    min_pathway_genes, neg_pvalue_cutoff, neg_fdr_cutoff, neg_r_cutoff, strict_r_cutoff,
    min_surv_n, min_surv_events, min_clin_n, min_group_n, clin_plot_p_cutoff
  ), file.path(out_dir_abs, "00_thresholds_used.csv"))

  ## 6. 逐个 GO 取基因（绝不合并）
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
  fwrite(go_set_summary, file.path(out_dir_abs, "00_GO_gene_sets.csv"))

  ## 7. 每个 GO 单独分析（for 必须包住 7.1–7.4 全部内容）
  all_scores <- list()
  all_clin   <- list()
  all_surv   <- list()
  all_neg    <- list()
  all_prot   <- list()
  all_meta   <- list()

  surv_endpoints <- list(
    OS  = c("OS",  "OS.time"),
    DSS = c("DSS", "DSS.time"),
    PFI = c("PFI", "PFI.time"),
    DFI = c("DFI", "DFI.time")
  )

  for (go_id in go_to_run) {
    go_title <- unname(go_name_map[go_id])
    if (is.na(go_title) || identical(go_title, "NA")) go_title <- go_id
    message("\n========== ", go_id, " | ", go_title, " ==========")

    go_dir <- file.path(out_dir_abs, "per_GO", safe_name(go_id))
    dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)
    message("  写入：", go_dir)

    genes_here <- unique(go_gene_list[[go_id]]$SYMBOL)
    fwrite(
      data.table(
        GO = go_id, GO_name = go_title, gene = genes_here,
        in_expression = genes_here %in% rownames(expr)
      ),
      file.path(go_dir, "genes.csv")
    )

    go_score <- tryCatch(pathway_zmean(expr, genes_here), error = function(e) {
      message("  通路打分失败：", conditionMessage(e))
      NULL
    })

    if (!is.numeric(go_score) || length(go_score) == 0) {
      message("  可用基因 < ", min_pathway_genes, " 或打分失败，跳过该通路：", go_id)
      writeLines(paste("SKIPPED", go_id, go_title), file.path(go_dir, "SKIPPED.txt"))
    } else {
      n_used <- attr(go_score, "n_genes")
      used_genes <- attr(go_score, "genes")
      message("  通路基因用于打分：", n_used)
      fwrite(
        data.table(sample = names(go_score), pathway_score = as.numeric(go_score)),
        file.path(go_dir, "pathway_score.csv")
      )
      all_scores[[go_id]] <- go_score

      # ---- 7.1 临床相关性（仅本通路分数） ----
      clin_use <- clinical_data[sample_std %in% names(go_score)]
      setkey(clin_use, sample_std)
      sc_clin <- go_score[clin_use$sample_std]
      clin_rows <- lapply(clin_features, function(ft) {
        if (!ft %in% names(clin_use)) return(NULL)
        assoc_clinical_feature(sc_clin, clin_use[[ft]], ft)
      })
      clin_tab <- rbindlist(clin_rows, fill = TRUE)
      meta_clin_tab <- collect_metastasis_clinical(go_id, go_title, sc_clin, clin_use)
      if (nrow(meta_clin_tab) > 0) {
        fwrite(meta_clin_tab, file.path(go_dir, "metastasis_clinical.csv"))
        all_meta[[go_id]] <- meta_clin_tab
      }
      if (nrow(clin_tab) > 0) {
        clin_tab[, `:=`(GO = go_id, GO_name = go_title, fdr = p.adjust(pvalue, method = "BH"))]
        setcolorder(clin_tab, c("GO", "GO_name"))
        fwrite(clin_tab, file.path(go_dir, "clinical_association.csv"))
        all_clin[[go_id]] <- clin_tab

        # 不要写 type == "categorical"：type 会撞上 BiocGenerics::type
        sig_cat <- clin_tab[clin_tab[["type"]] == "categorical" & clin_tab[["pvalue"]] < clin_plot_p_cutoff]
        if (nrow(sig_cat) > 0) sig_cat <- sig_cat[order(sig_cat[["pvalue"]])]
        if (nrow(sig_cat) == 0) {
          message("  本通路无分类临床变量达到 p < ", clin_plot_p_cutoff, "，不画箱线图")
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
            plot_df <- plot_df[plot_df$group %in% names(tab)[tab >= min_group_n], ]
            if (length(unique(plot_df$group)) < 2) next
            ok <- tryCatch({
              p <- ggboxplot(plot_df, x = "group", y = "score", fill = "group", outlier.shape = NA) +
                stat_compare_means() +
                labs(
                  title = paste0(go_id, "\n", go_title),
                  subtitle = paste0(ft, "  (p=", signif(sig_cat$pvalue[sig_cat$feature == ft][1], 3), ")"),
                  x = NULL, y = "Pathway score"
                ) +
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
      }

      # ---- 7.2 生存分析（仅本通路分数） ----
      surv_use <- survival_data[sample_std %in% names(go_score)]
      sc_surv <- go_score[surv_use$sample_std]
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
        if (nrow(d) < min_surv_n || sum(d$event) < min_surv_events) next

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
      cor_tab <- spearman_vs_score(expr, go_score)
      if (nrow(cor_tab) > 0) {
        setnames(cor_tab, "feature", "gene")
        cor_tab[, `:=`(
          GO = go_id,
          GO_name = go_title,
          in_this_GO = gene %in% used_genes
        )]
        fwrite(cor_tab[order(-spearman_r)], file.path(go_dir, "all_gene_correlation.csv"))

        neg_tab <- cor_tab[spearman_r < neg_r_cutoff & pvalue < neg_pvalue_cutoff]
        setorder(neg_tab, spearman_r)
        fwrite(neg_tab, file.path(go_dir, "negative_correlated_genes.csv"))
        fwrite(
          cor_tab[spearman_r < 0 & fdr < neg_fdr_cutoff][order(spearman_r)],
          file.path(go_dir, "negative_correlated_genes_fdr.csv")
        )
        fwrite(
          neg_tab[spearman_r <= strict_r_cutoff],
          file.path(go_dir, "negative_correlated_genes_strict.csv")
        )
        all_neg[[go_id]] <- neg_tab
        message(
          "  负相关基因（r<0 且 p<", neg_pvalue_cutoff, "）：", nrow(neg_tab),
          " ；r<=", strict_r_cutoff, "：", nrow(neg_tab[spearman_r <= strict_r_cutoff])
        )

        vol_df <- as.data.frame(cor_tab)
        rownames(vol_df) <- vol_df$gene
        key_neg <- head(neg_tab$gene, 15)
        tryCatch({
          pdf(file.path(go_dir, "volcano_gene_correlation.pdf"), width = 9, height = 7)
          print(EnhancedVolcano(
            vol_df,
            lab = vol_df$gene,
            x = "spearman_r",
            y = "pvalue",
            xlab = "Spearman r (gene vs this GO score)",
            ylab = "-Log10 p",
            title = paste(go_id, go_title),
            pCutoff = neg_pvalue_cutoff,
            FCcutoff = abs(strict_r_cutoff),
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

      # ---- 7.4 蛋白水平相关性（仍按本通路分数） ----
      if (!is.null(prot_mat)) {
        prot_common <- intersect(colnames(prot_mat), names(go_score))
        if (length(prot_common) >= 20) {
          prot_cor <- spearman_vs_score(prot_mat, go_score)
          if (nrow(prot_cor) > 0) {
            setnames(prot_cor, "feature", "protein")
            prot_cor[, `:=`(GO = go_id, GO_name = go_title)]
            fwrite(prot_cor[order(spearman_r)], file.path(go_dir, "protein_correlation.csv"))
            all_prot[[go_id]] <- prot_cor[spearman_r < 0 & pvalue < neg_pvalue_cutoff][order(spearman_r)]
          }
        }
      }
    }  # 结束 else（打分成功才做 7.1–7.4）
  }    # 结束 for (go_id in go_to_run)

  ## 8. 汇总输出（每个 GO 一行/一堆结果，不是合并通路）
  if (length(all_scores) > 0) {
    score_mat <- do.call(cbind, lapply(all_scores, function(x) {
      x[colnames(expr)]
    }))
    colnames(score_mat) <- names(all_scores)
    rownames(score_mat) <- colnames(expr)
    fwrite(
      data.table(sample = rownames(score_mat), as.data.table(score_mat)),
      file.path(out_dir_abs, "01_pathway_scores_each_GO.csv")
    )

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
    pdf_hm <- file.path(out_dir_abs, "01_pathway_score_heatmap.pdf")
    plot_dt <- as.data.table(as.table(z_score))
    setnames(plot_dt, c("GO", "sample", "z"))
    p_hm <- ggplot(plot_dt, aes(x = sample, y = GO, fill = z)) +
      geom_tile() +
      scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35") +
      labs(title = "Each GO pathway score (not pooled)", x = NULL, y = NULL, fill = "z-score") +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank()
      )
    ggsave(pdf_hm, p_hm, width = 12, height = 6)
    message("  已保存通路分数热图：", pdf_hm)
  }

  if (length(all_clin) > 0) {
    clin_all <- rbindlist(all_clin, fill = TRUE)
    fwrite(clin_all, file.path(out_dir_abs, "02_clinical_association_each_GO.csv"))
  }

  if (length(all_surv) > 0) {
    surv_all <- rbindlist(all_surv, fill = TRUE)
    fwrite(surv_all, file.path(out_dir_abs, "03_survival_cox_each_GO.csv"))

    os_hl <- surv_all[surv_all[["endpoint"]] == "OS" & surv_all[["model"]] == "High_vs_Low"]
    if (nrow(os_hl) > 0) {
      os_hl[, lab := paste(GO, GO_name)]
      os_hl[, lab := factor(lab, levels = rev(lab))]
      pfor <- ggplot(os_hl, aes(x = HR, y = lab)) +
        geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
        geom_errorbar(aes(xmin = HR_low, xmax = HR_high), orientation = "y", width = 0.2) +
        geom_point(aes(color = pvalue < 0.05), size = 3) +
        scale_x_log10() +
        labs(
          title = "OS Cox: High vs Low (each GO separately)",
          x = "Hazard ratio", y = NULL, color = "p < 0.05"
        ) +
        theme_bw()
      ggsave(file.path(out_dir_abs, "03_OS_forest_each_GO.pdf"), pfor, width = 10, height = 6)
      message("  已保存 OS 森林图：", file.path(out_dir_abs, "03_OS_forest_each_GO.pdf"))
    }
  }

  if (length(all_meta) > 0) {
    fwrite(rbindlist(all_meta, fill = TRUE), file.path(out_dir_abs, "07_metastasis_clinical_each_GO.csv"))
  }

  tryCatch({
    draw_prognosis_metastasis_bubbles(
      result_dir = out_dir_abs,
      surv_all = if (length(all_surv) > 0) rbindlist(all_surv, fill = TRUE) else NULL,
      all_scores = all_scores
    )
  }, error = function(e) {
    message("  气泡图失败：", conditionMessage(e))
  })

  if (length(all_neg) > 0) {
    neg_all <- rbindlist(all_neg, fill = TRUE)
    fwrite(neg_all, file.path(out_dir_abs, "04_negative_genes_each_GO.csv"))
    neg_count <- neg_all[, .(
      n_negative_genes = .N,
      n_strict = sum(spearman_r <= strict_r_cutoff)
    ), by = .(GO, GO_name)]
    fwrite(neg_count, file.path(out_dir_abs, "04_negative_genes_count_each_GO.csv"))
  }

  if (length(all_prot) > 0) {
    fwrite(
      rbindlist(all_prot, fill = TRUE),
      file.path(out_dir_abs, "05_negative_proteins_each_GO.csv")
    )
  }

  fwrite(go_set_summary, file.path(out_dir_abs, "00_GO_gene_sets.csv"))
  folders <- list.files(file.path(out_dir_abs, "per_GO"))
  fwrite(data.table(folder = folders), file.path(out_dir_abs, "00_per_GO_folder_list.csv"))
  message("\n分析完成。每个 GO 的独立结果在：", normalizePath(out_dir_abs, mustWork = FALSE))
  message("per_GO 现有 ", length(folders), " 个文件夹：")
  message(paste("  ", folders, collapse = "\n"))
  message("请重点查看 per_GO/ 下各通路文件夹，以及 04_negative_genes_each_GO.csv")
  message("预后气泡图：06_prognosis_bubble_each_GO.pdf ；转移气泡图：07_metastasis_bubble_each_GO.pdf")
  invisible(TRUE)
}

# 数据已经读完、expr 已建好之后，才调用一次。不要把这一行挪到脚本开头。
# 17 个 GO 若已经跑完、只想补预后/转移气泡图：把 run_mode 改成 "bubbles_only"
# 或在控制台直接运行（需已有 03_survival_cox_each_GO.csv）：
#   draw_prognosis_metastasis_bubbles()
if (identical(run_mode, "bubbles_only")) {
  draw_prognosis_metastasis_bubbles()
} else {
  go_to_run <- go_list
  run_go_individual_analysis()
}

# ==============================================================================
# 第九部分（追加，不改上面已有代码）：
# 5 个神经/轴突 GO 与肿瘤增殖、转移通路的相关性气泡图
# 横轴：肿瘤增殖 / 转移信号通路
# 纵轴：下面 5 个 GO 与各肿瘤通路的 Spearman 相关
#   GO:2001222  regulation of neuron migration
#   GO:1904340  positive regulation of dopaminergic neuron differentiation
#   GO:1902667  regulation of axon guidance
#   GO:0097374  sensory neuron axon guidance
#   GO:0036518  chemorepulsion of dopaminergic neuron axon
# 需要已经有表达矩阵 expr（上面分析跑完即可）。也可单独选中本段运行。
# ==============================================================================
focus5_go_map <- c(
  "GO:2001222" = "regulation of neuron migration",
  "GO:1904340" = "positive regulation of dopaminergic neuron differentiation",
  "GO:1902667" = "regulation of axon guidance",
  "GO:0097374" = "sensory neuron axon guidance",
  "GO:0036518" = "chemorepulsion of dopaminergic neuron axon"
)

tumor_path_dt <- data.table(
  GO = c(
    "GO:0008283", "GO:0008284", "GO:0007049", "GO:0051301",
    "GO:0006260", "GO:0000082", "GO:0045787",
    "GO:0016477", "GO:0030335", "GO:0001837", "GO:0001525",
    "GO:0045766", "GO:0030198", "GO:0007155", "GO:0001666",
    "GO:0000165", "GO:0016055", "GO:0007179", "GO:0048015"
  ),
  path_name = c(
    "cell proliferation",
    "pos. reg. of cell proliferation",
    "cell cycle",
    "cell division",
    "DNA replication",
    "G1/S transition",
    "pos. reg. of cell cycle",
    "cell migration",
    "pos. reg. of cell migration",
    "EMT",
    "angiogenesis",
    "pos. reg. of angiogenesis",
    "ECM organization",
    "cell adhesion",
    "response to hypoxia",
    "MAPK cascade",
    "Wnt signaling",
    "TGF-beta receptor signaling",
    "PI3K-mediated signaling"
  ),
  category = c(
    rep("增殖", 7),
    rep("转移", 12)
  )
)

plot_focus5_vs_tumor_pathways <- function() {
  if (!exists("expr", inherits = TRUE) || is.null(get("expr", inherits = TRUE))) {
    stop("还没有表达矩阵 expr。请先跑完上面的分析，或从脚本开头 Source 到数据预处理结束，再运行本段。")
  }
  expr_use <- get("expr", inherits = TRUE)
  res_dir <- if (exists("out_dir", inherits = TRUE)) {
    normalizePath(out_dir, mustWork = FALSE)
  } else {
    "results_GO_individual"
  }
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

  get_score_one <- function(go_id) {
    score_file <- file.path(res_dir, "01_pathway_scores_each_GO.csv")
    if (file.exists(score_file)) {
      sm <- fread(score_file)
      if (go_id %in% names(sm)) {
        v <- as.numeric(sm[[go_id]])
        names(v) <- as.character(sm$sample)
        return(v)
      }
    }
    genes <- unique(get_go_genes(go_id)$SYMBOL)
    sc <- tryCatch(pathway_zmean(expr_use, genes), error = function(e) NULL)
    sc
  }

  message("计算 5 个关注 GO 的通路分数...")
  focus_scores <- lapply(names(focus5_go_map), get_score_one)
  names(focus_scores) <- names(focus5_go_map)
  n_ok <- vapply(focus_scores, function(x) is.numeric(x) && length(x) > 0, logical(1))
  if (!any(n_ok)) stop("5 个关注 GO 都没有得到通路分数")
  for (g in names(focus_scores)[!n_ok]) {
    message("  跳过（无基因或打分失败）：", g, " | ", unname(focus5_go_map[g]))
  }

  message("计算肿瘤增殖/转移通路分数...")
  tumor_scores <- lapply(tumor_path_dt$GO, get_score_one)
  names(tumor_scores) <- tumor_path_dt$GO
  t_ok <- vapply(tumor_scores, function(x) is.numeric(x) && length(x) > 0, logical(1))
  if (!any(t_ok)) stop("肿瘤增殖/转移通路都没有得到通路分数")
  for (g in names(tumor_scores)[!t_ok]) {
    message("  跳过肿瘤通路（无基因或打分失败）：", g, " | ", tumor_path_dt[GO == g, path_name])
  }

  cor_two <- function(a, b) {
    common <- intersect(names(a), names(b))
    if (length(common) < 10) return(NULL)
    aa <- as.numeric(a[common])
    bb <- as.numeric(b[common])
    keep <- is.finite(aa) & is.finite(bb)
    if (sum(keep) < 10) return(NULL)
    ct <- suppressWarnings(cor.test(aa[keep], bb[keep], method = "spearman", exact = FALSE))
    data.table(spearman_r = unname(ct$estimate), pvalue = ct$p.value, n = sum(keep))
  }

  cor_rows <- list()
  for (fg in names(focus_scores)[n_ok]) {
    for (tg in names(tumor_scores)[t_ok]) {
      one <- cor_two(focus_scores[[fg]], tumor_scores[[tg]])
      if (is.null(one)) next
      one[, `:=`(
        focus_GO = fg,
        focus_name = unname(focus5_go_map[fg]),
        tumor_GO = tg,
        tumor_name = tumor_path_dt[GO == tg, path_name],
        category = tumor_path_dt[GO == tg, category]
      )]
      cor_rows[[paste(fg, tg, sep = "_")]] <- one
    }
  }
  cor_dt <- rbindlist(cor_rows, fill = TRUE)
  if (nrow(cor_dt) == 0) stop("没有算出任何通路-通路相关")
  cor_dt[, fdr := p.adjust(pvalue, method = "BH")]
  cor_dt[, neglogp := pmin(10, -log10(pmax(pvalue, 1e-12)))]
  cor_dt[, focus_lab := factor(
    paste0(focus_GO, "\n", focus_name),
    levels = rev(paste0(names(focus5_go_map), "\n", unname(focus5_go_map)))
  )]
  cor_dt[, tumor_lab := factor(tumor_name, levels = tumor_path_dt$path_name)]
  cor_dt[, category := factor(category, levels = c("增殖", "转移"))]
  fwrite(cor_dt, file.path(res_dir, "08_focus5GO_vs_prolif_met_correlation.csv"))

  lim <- max(0.3, as.numeric(stats::quantile(abs(cor_dt$spearman_r), 0.98, na.rm = TRUE)), na.rm = TRUE)
  p_bub <- ggplot(cor_dt, aes(x = tumor_lab, y = focus_lab, size = neglogp, color = spearman_r)) +
    geom_point(alpha = 0.9) +
    facet_grid(. ~ category, scales = "free_x", space = "free_x") +
    scale_size_continuous(range = c(2, 11), name = expression(-log[10](p))) +
    scale_color_gradient2(
      low = "#3C5488", mid = "grey90", high = "#E64B35",
      midpoint = 0, limits = c(-lim, lim), name = "Spearman r"
    ) +
    labs(
      title = "5 个 GO 通路与肿瘤增殖、转移通路的相关性",
      subtitle = "每个通路单独打分后做 Spearman；未把 GO 基因集合并",
      x = "肿瘤增殖 / 转移信号通路",
      y = "关注的 5 个 GO 通路"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      strip.background = element_rect(fill = "grey95"),
      panel.grid.major.x = element_blank()
    )
  pdf_out <- file.path(res_dir, "08_focus5GO_vs_prolif_met_bubble.pdf")
  ggsave(pdf_out, p_bub, width = 13, height = 6.2)
  ggsave(file.path(res_dir, "08_focus5GO_vs_prolif_met_bubble.png"), p_bub, width = 13, height = 6.2, dpi = 150)
  message("已保存：", pdf_out)
  print(p_bub)
  invisible(cor_dt)
}

if (exists("expr", inherits = TRUE) && !is.null(get("expr", inherits = TRUE))) {
  plot_focus5_vs_tumor_pathways()
} else {
  message("没有 expr，跳过 5 个 GO 与增殖/转移通路气泡图。有表达矩阵后运行：plot_focus5_vs_tumor_pathways()")
}

