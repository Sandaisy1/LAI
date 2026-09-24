## 只画两张气泡图：预后 + 转移
## 整段复制到 RStudio 运行。需要已经有：
##   E:/R/BRCA/results_GO_individual/03_survival_cox_each_GO.csv
##   E:/R/BRCA/results_GO_individual/01_pathway_scores_each_GO.csv（转移临床部分用）
## 不要 Source 整份分析脚本，不要改 run_mode。

library(data.table)
library(ggplot2)

setwd("E:/R/BRCA")
res_dir <- "results_GO_individual"
surv_file <- file.path(res_dir, "03_survival_cox_each_GO.csv")
score_file <- file.path(res_dir, "01_pathway_scores_each_GO.csv")
if (!file.exists(surv_file)) stop("找不到 ", surv_file, " ，等 17 个 GO 跑完再贴这段")

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

make_bubble <- function(d, title, color_name, outfile) {
  d <- as.data.table(d)
  d <- d[is.finite(pvalue) & is.finite(effect)]
  if (nrow(d) == 0) {
    message("没有可画的点：", title)
    return(invisible(NULL))
  }
  d[, neglogp := pmin(10, -log10(pmax(pvalue, 1e-12)))]
  d[, go_lab := paste0(GO, "  ", GO_name)]
  go_rank <- d[, .(m = max(neglogp, na.rm = TRUE)), by = go_lab]
  setorder(go_rank, m)
  d[, go_lab := factor(go_lab, levels = go_rank$go_lab)]
  if (!is.factor(d$feature_lab)) d[, feature_lab := factor(feature_lab, levels = unique(as.character(feature_lab)))]
  lim <- max(0.4, as.numeric(stats::quantile(abs(d$effect), 0.98, na.rm = TRUE)), na.rm = TRUE)
  ht <- max(6.5, min(14, 0.38 * uniqueN(d$go_lab) + 3.2))
  p <- ggplot(d, aes(x = feature_lab, y = go_lab, size = neglogp, color = effect)) +
    geom_point(alpha = 0.9) +
    scale_size_continuous(range = c(2.5, 12), name = expression(-log[10](p))) +
    scale_color_gradient2(
      low = "#3C5488", mid = "grey90", high = "#E64B35",
      midpoint = 0, limits = c(-lim, lim), name = color_name
    ) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid.major.x = element_blank())
  ggsave(paste0(outfile, ".pdf"), p, width = 10, height = ht)
  ggsave(paste0(outfile, ".png"), p, width = 10, height = ht, dpi = 150)
  message("已保存：", outfile, ".pdf")
  print(p)
  invisible(p)
}

## ---- 预后：OS / DSS ----
surv_all <- fread(surv_file)
hl <- surv_all[surv_all[["model"]] == "High_vs_Low"]
if (nrow(hl) == 0) hl <- surv_all[surv_all[["model"]] == "continuous_score"]
hl <- hl[is.finite(HR) & HR > 0]
hl[, effect := log2(HR)]
hl[is.na(GO_name) | GO_name == "", GO_name := unname(go_name_map[GO])]

prog <- hl[endpoint %in% c("OS", "DSS")]
prog[, feature_lab := factor(
  fifelse(endpoint == "OS", "OS 总生存", "DSS 疾病特异生存"),
  levels = c("OS 总生存", "DSS 疾病特异生存")
)]
fwrite(prog, file.path(res_dir, "06_prognosis_bubble_data.csv"))
make_bubble(
  prog,
  "各 GO 通路与乳腺癌预后",
  "log2(HR)",
  file.path(res_dir, "06_prognosis_bubble_each_GO")
)

## ---- 转移：PFI / DFI + M / N / 分期 ----
meta_s <- hl[endpoint %in% c("PFI", "DFI")]
meta_s[, feature_lab := fifelse(endpoint == "PFI", "PFI 无进展间隔", "DFI 无病间隔")]

meta_c <- data.table()
if (file.exists(score_file) && file.exists("TCGA-BRCA.clinical.tsv")) {
  clin <- fread("TCGA-BRCA.clinical.tsv")
  idc <- first_present(names(clin), c("sampleID", "sample", "bcr_patient_barcode", "submitter_id", names(clin)[1]))
  clin[, sample_std := normalize_barcode(get(idc))]
  clin <- clin[!duplicated(sample_std)]
  st <- first_present(names(clin), c("ajcc_pathologic_tumor_stage", "pathologic_stage", "clinical_stage", "ajcc_pathologic_stage"))
  mc <- first_present(names(clin), c("ajcc_metastasis_pathologic_pm", "pathologic_M", "pathologic_m", "ajcc_pathologic_m"))
  nc <- first_present(names(clin), c("ajcc_nodes_pathologic_pn", "pathologic_N", "pathologic_n", "ajcc_pathologic_n"))
  if (!is.na(st)) clin[, meta_stage := ifelse(simplify_stage(get(st)) %in% c("Stage I", "Stage II"), "Early_I_II",
                                      ifelse(simplify_stage(get(st)) %in% c("Stage III", "Stage IV"), "Late_III_IV", NA_character_))]
  if (!is.na(mc)) clin[, meta_M := classify_m(get(mc))]
  if (!is.na(nc)) clin[, meta_N := classify_n(get(nc))]

  sm <- fread(score_file)
  go_cols <- setdiff(names(sm), "sample")
  two_group <- function(sc, grp, pos, neg, feat, g, gnm) {
    df <- data.frame(score = as.numeric(sc), group = as.character(grp), stringsAsFactors = FALSE)
    df <- df[is.finite(df$score) & df$group %in% c(pos, neg), ]
    if (sum(df$group == pos) < 2 || sum(df$group == neg) < 2) return(NULL)
    wt <- suppressWarnings(wilcox.test(score ~ group, data = df))
    data.table(
      GO = g, GO_name = gnm, feature = feat,
      effect = stats::median(df$score[df$group == pos], na.rm = TRUE) -
        stats::median(df$score[df$group == neg], na.rm = TRUE),
      pvalue = wt$p.value, n = nrow(df)
    )
  }
  meta_c <- rbindlist(lapply(go_cols, function(g) {
    sc <- as.numeric(sm[[g]])
    names(sc) <- as.character(sm$sample)
    gnm <- unname(go_name_map[g])
    if (is.na(gnm) || gnm == "") gnm <- g
    cu <- clin[sample_std %in% names(sc)]
    sc2 <- sc[cu$sample_std]
    rbindlist(list(
      if ("meta_M" %in% names(cu)) two_group(sc2, cu$meta_M, "M1", "M0", "远处转移 M1 vs M0", g, gnm) else NULL,
      if ("meta_N" %in% names(cu)) two_group(sc2, cu$meta_N, "Nplus", "N0", "淋巴结 N+ vs N0", g, gnm) else NULL,
      if ("meta_stage" %in% names(cu)) two_group(sc2, cu$meta_stage, "Late_III_IV", "Early_I_II", "分期 III-IV vs I-II", g, gnm) else NULL
    ), fill = TRUE)
  }), fill = TRUE)
}

meta_all <- rbindlist(list(
  if (nrow(meta_s) > 0) meta_s[, .(GO, GO_name, feature_lab, effect, pvalue, n)] else NULL,
  if (nrow(meta_c) > 0) {
    meta_c[, feature_lab := feature]
    meta_c[, .(GO, GO_name, feature_lab, effect, pvalue, n)]
  } else NULL
), fill = TRUE)
if (nrow(meta_all) > 0) {
  lv <- c("PFI 无进展间隔", "DFI 无病间隔", "远处转移 M1 vs M0", "淋巴结 N+ vs N0", "分期 III-IV vs I-II")
  lv <- lv[lv %in% as.character(meta_all$feature_lab)]
  meta_all[, feature_lab := factor(as.character(feature_lab), levels = lv)]
  fwrite(meta_all, file.path(res_dir, "07_metastasis_bubble_data.csv"))
  make_bubble(
    meta_all,
    "各 GO 通路与乳腺癌转移",
    "effect",
    file.path(res_dir, "07_metastasis_bubble_each_GO")
  )
} else {
  message("没有转移相关结果可画")
}
