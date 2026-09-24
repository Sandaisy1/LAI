## 整段复制到 RStudio 运行。不要 Source 整份分析脚本。
## 需要：
##   E:/R/BRCA/results_GO_individual/01_pathway_scores_each_GO.csv
##   E:/R/BRCA/TCGA-BRCA.clinical.tsv
##   E:/R/BRCA/TCGA-BRCA.survival.tsv（预后，可缺）
## 转移 = M1 或 Stage IV；未转移 = M0 且不是 IV。每个 GO 单独比较。

library(data.table)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)

setwd("E:/R/BRCA")
res_dir <- "results_GO_individual"
score_file <- file.path(res_dir, "01_pathway_scores_each_GO.csv")
if (!file.exists(score_file)) stop("找不到 ", score_file)
if (!file.exists("TCGA-BRCA.clinical.tsv")) stop("找不到 TCGA-BRCA.clinical.tsv")

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
  "GO:0007409" = "axonogenesis",
  "GO:0036518" = "chemorepulsion of dopaminergic neuron axon"
)

normalize_barcode <- function(x) {
  x <- toupper(gsub("\\.", "-", as.character(x)))
  x <- sub("A$", "", x)
  ifelse(nchar(x) >= 15, substr(x, 1, 15), x)
}
first_present <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) NA_character_ else hit[1]
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
classify_m <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("\\bM1\\b|M1[ABC]?", x)] <- "M1"
  out[is.na(out) & grepl("\\bM0\\b", x)] <- "M0"
  out
}
classify_n <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("\\bN[1-3]", x)] <- "Nplus"
  out[is.na(out) & grepl("\\bN0\\b", x)] <- "N0"
  out
}

sm <- fread(score_file)
go_cols <- setdiff(names(sm), "sample")
score_mat <- as.matrix(sm[, go_cols, with = FALSE])
storage.mode(score_mat) <- "double"
rownames(score_mat) <- as.character(sm$sample)
colnames(score_mat) <- go_cols

clin <- fread("TCGA-BRCA.clinical.tsv")
idc <- first_present(names(clin), c("sampleID", "sample", "bcr_patient_barcode", "submitter_id", names(clin)[1]))
clin[, sample_std := normalize_barcode(get(idc))]
clin <- clin[!duplicated(sample_std)]
st <- first_present(names(clin), c("ajcc_pathologic_tumor_stage", "pathologic_stage", "clinical_stage", "ajcc_pathologic_stage"))
mc <- first_present(names(clin), c("ajcc_metastasis_pathologic_pm", "pathologic_M", "pathologic_m", "ajcc_pathologic_m"))
nc <- first_present(names(clin), c("ajcc_nodes_pathologic_pn", "pathologic_N", "pathologic_n", "ajcc_pathologic_n"))
if (!is.na(st)) clin[, stage_simplified := simplify_stage(get(st))]
if (!is.na(mc)) clin[, meta_M := classify_m(get(mc))]
if (!is.na(nc)) clin[, meta_N := classify_n(get(nc))]

surv <- if (file.exists("TCGA-BRCA.survival.tsv")) fread("TCGA-BRCA.survival.tsv") else NULL
if (!is.null(surv)) {
  sid <- first_present(names(surv), c("sample", "sampleID", "bcr_patient_barcode", names(surv)[1]))
  surv[, sample_std := normalize_barcode(get(sid))]
  surv <- surv[!duplicated(sample_std)]
}

ann <- data.table(sample = rownames(score_mat))
ann <- merge(ann, clin, by.x = "sample", by.y = "sample_std", all.x = TRUE)
if (!is.null(surv)) {
  keep <- intersect(c("sample_std", "OS", "OS.time", "PFI", "PFI.time"), names(surv))
  ann <- merge(ann, surv[, keep, with = FALSE], by.x = "sample", by.y = "sample_std", all.x = TRUE)
}

ann[, distant_M := factor(meta_M, levels = c("M0", "M1"))]
ann[, node_N := factor(meta_N, levels = c("N0", "Nplus"))]
ann[, any_met := NA_character_]
ann[meta_M == "M1" | stage_simplified == "Stage IV", any_met := "转移"]
ann[is.na(any_met) & meta_M == "M0" & (is.na(stage_simplified) | stage_simplified != "Stage IV"), any_met := "未转移"]
ann[is.na(any_met) & stage_simplified %in% c("Stage I", "Stage II", "Stage III") & (is.na(meta_M) | meta_M != "M1"), any_met := "未转移"]
ann[, any_met := factor(any_met, levels = c("未转移", "转移"))]
if ("PFI" %in% names(ann)) {
  ann[, progressed := NA_character_]
  ann[as.numeric(PFI) == 0, progressed := "未进展"]
  ann[as.numeric(PFI) == 1, progressed := "进展"]
  ann[, progressed := factor(progressed, levels = c("未进展", "进展"))]
}

fwrite(ann[, intersect(c("sample", "any_met", "distant_M", "node_N", "progressed",
                         "stage_simplified", "OS", "OS.time", "PFI", "PFI.time"), names(ann)), with = FALSE],
       file.path(res_dir, "09_sample_metastasis_prognosis.csv"))
message("转移=", sum(ann$any_met == "转移", na.rm = TRUE),
        "  未转移=", sum(ann$any_met == "未转移", na.rm = TRUE),
        "  M1=", sum(ann$distant_M == "M1", na.rm = TRUE))

compare_one <- function(go_score_vec, group, pos, neg, grouping) {
  df <- data.frame(pathway_score = as.numeric(go_score_vec), group = as.character(group), stringsAsFactors = FALSE)
  df <- df[is.finite(df$pathway_score) & df$group %in% c(pos, neg), ]
  if (sum(df$group == pos) < 2 || sum(df$group == neg) < 2) return(NULL)
  wt <- suppressWarnings(wilcox.test(pathway_score ~ group, data = df))
  data.table(
    grouping = grouping, pos_level = pos, neg_level = neg,
    n_pos = sum(df$group == pos), n_neg = sum(df$group == neg),
    median_pos = stats::median(df$pathway_score[df$group == pos], na.rm = TRUE),
    median_neg = stats::median(df$pathway_score[df$group == neg], na.rm = TRUE),
    delta_median = stats::median(df$pathway_score[df$group == pos], na.rm = TRUE) -
      stats::median(df$pathway_score[df$group == neg], na.rm = TRUE),
    pvalue = wt$p.value
  )
}

stat_rows <- list(); long_rows <- list()
for (g in colnames(score_mat)) {
  gnm <- if (g %in% names(go_name_map)) unname(go_name_map[g]) else g
  sc <- as.numeric(score_mat[, g]); names(sc) <- rownames(score_mat); sc <- sc[ann$sample]
  one <- rbindlist(Filter(Negate(is.null), list(
    compare_one(sc, ann$any_met, "转移", "未转移", "any_met"),
    compare_one(sc, ann$distant_M, "M1", "M0", "distant_M"),
    compare_one(sc, ann$node_N, "Nplus", "N0", "node_N"),
    if ("progressed" %in% names(ann)) compare_one(sc, ann$progressed, "进展", "未进展", "PFI_progressed") else NULL
  )), fill = TRUE)
  if (nrow(one) > 0) { one[, `:=`(GO = g, GO_name = gnm)]; stat_rows[[g]] <- one }
  long_rows[[g]] <- data.table(sample = ann$sample, GO = g, GO_name = gnm,
                               pathway_score = as.numeric(sc), any_met = ann$any_met)
}
stat_dt <- rbindlist(stat_rows, fill = TRUE)
stat_dt[, fdr := p.adjust(pvalue, method = "BH"), by = grouping]
fwrite(stat_dt, file.path(res_dir, "09_nerve_GO_score_by_metastasis.csv"))

long_any <- rbindlist(long_rows, fill = TRUE)
long_any <- long_any[!is.na(any_met) & is.finite(pathway_score)]
long_any[, go_lab := factor(paste(GO, GO_name), levels = unique(paste(GO, GO_name)))]
p_box <- ggplot(long_any, aes(x = any_met, y = pathway_score, fill = any_met)) +
  geom_boxplot(outlier.size = 0.4, width = 0.65) +
  stat_compare_means(size = 2.6, label = "p.format") +
  facet_wrap(~ go_lab, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("未转移" = "#4DBBD5", "转移" = "#E64B35")) +
  labs(title = "神经相关 GO 在转移 vs 未转移中的表达",
       subtitle = "转移 = M1 或 Stage IV；每个 GO 单独打分",
       x = NULL, y = "Pathway score", fill = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(size = 7),
        axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(res_dir, "09_nerve_GO_boxplot_any_met.pdf"), p_box, width = 12, height = 10)
ggsave(file.path(res_dir, "09_nerve_GO_boxplot_any_met.png"), p_box, width = 12, height = 10, dpi = 150)
print(p_box)

any_stat <- stat_dt[grouping == "any_met"]
any_stat[, lab := factor(paste(GO, GO_name), levels = paste(GO, GO_name)[order(delta_median)])]
p_for <- ggplot(any_stat, aes(x = delta_median, y = lab)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_point(aes(color = pvalue < 0.05, size = -log10(pmax(pvalue, 1e-12)))) +
  scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "#E64B35"), name = "p < 0.05") +
  labs(title = "转移 − 未转移：神经 GO 通路分数差",
       x = "Δ median (转移 − 未转移)", y = NULL, size = expression(-log[10](p))) +
  theme_bw()
ggsave(file.path(res_dir, "09_nerve_GO_delta_any_met.pdf"), p_for, width = 10, height = 6)
print(p_for)

surv_rows <- list()
if (all(c("OS", "OS.time") %in% names(ann))) {
  pdf(file.path(res_dir, "09_nerve_GO_KM_OS_within_met.pdf"), width = 10, height = 5)
  n_km <- 0L
  for (g in colnames(score_mat)) {
    gnm <- if (g %in% names(go_name_map)) unname(go_name_map[g]) else g
    d0 <- data.frame(
      time = as.numeric(ann[["OS.time"]]), event = as.numeric(ann[["OS"]]),
      pathway_score = as.numeric(score_mat[ann$sample, g]),
      met = as.character(ann$any_met), stringsAsFactors = FALSE
    )
    d0 <- d0[is.finite(d0$time) & d0$time > 0 & d0$event %in% c(0, 1) & is.finite(d0$pathway_score) & !is.na(d0$met), ]
    for (mg in c("未转移", "转移")) {
      d <- d0[d0$met == mg, ]
      if (nrow(d) < 10 || sum(d$event) < 3) next
      d$group <- factor(ifelse(d$pathway_score >= median(d$pathway_score, na.rm = TRUE), "High", "Low"),
                        levels = c("Low", "High"))
      cox_g <- tryCatch(coxph(Surv(time, event) ~ group, data = d), error = function(e) NULL)
      if (!is.null(cox_g)) {
        s2 <- summary(cox_g)
        surv_rows[[paste(g, mg)]] <- data.table(
          GO = g, GO_name = gnm, met_group = mg, n = nrow(d), events = sum(d$event),
          HR = s2$conf.int[1, 1], HR_low = s2$conf.int[1, 3], HR_high = s2$conf.int[1, 4],
          pvalue = s2$coefficients[1, "Pr(>|z|)"]
        )
      }
      fit <- tryCatch(survfit(Surv(time / 30.44, event) ~ group, data = d), error = function(e) NULL)
      if (!is.null(fit)) {
        print(ggsurvplot(fit, data = d, pval = TRUE, risk.table = TRUE,
                         legend.labs = c("Low", "High"),
                         xlab = "Time (months)", ylab = "Overall survival",
                         title = paste0(g, " | ", gnm, "\n", mg, " 亚组 High vs Low"),
                         ggtheme = theme_bw()))
        n_km <- n_km + 1L
      }
    }
  }
  invisible(dev.off())
  message("KM 页数：", n_km)
}

if (length(surv_rows) > 0) {
  surv_dt <- rbindlist(surv_rows, fill = TRUE)
  fwrite(surv_dt, file.path(res_dir, "09_nerve_GO_OS_cox_within_met.csv"))
  p_hr <- ggplot(surv_dt, aes(x = HR, y = reorder(paste(GO, met_group), HR))) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
    geom_errorbar(aes(xmin = HR_low, xmax = HR_high), orientation = "y", width = 0.2) +
    geom_point(aes(color = pvalue < 0.05), size = 2.5) +
    scale_x_log10() +
    facet_wrap(~ met_group, scales = "free_y") +
    labs(title = "OS High vs Low（转移 / 未转移亚组内）", x = "Hazard ratio", y = NULL, color = "p < 0.05") +
    theme_bw()
  ggsave(file.path(res_dir, "09_nerve_GO_OS_forest_within_met.pdf"), p_hr, width = 11, height = 7)
  print(p_hr)
}

message("完成。主表：", file.path(res_dir, "09_nerve_GO_score_by_metastasis.csv"))
print(stat_dt[grouping == "any_met"][order(pvalue)])
