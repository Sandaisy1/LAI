#!/usr/bin/env Rscript
################################################################################
# E:/R/FUSC：FUSCCTNBC Excel 临床分型 + FPKM，每个 GO 单独分析
#
# 输入文件（必须放在 E:/R/FUSC，扩展名 .xlsx 或 .xls）：
#   1) OEP00000155_样本_Human_1786804447300   病人信息（分型/预后）
#   2) FUSCCTNBC_Expression_RNAseqFPKM         RNA-seq FPKM
#
# 任务 1：各 GO 通路 vs 乳腺癌分型、预后（气泡图）
# 任务 2：寻找与各 GO 通路活性呈负相关的基因
#
# ============================ 请按这个跑 ============================
# 1. 两份 Excel 必须在 E:/R/FUSC（资源管理器里应能同时看到它们）
# 2. RStudio：Session -> Restart R
# 3. 控制台先运行：setwd("E:/R/FUSC"); getwd()
# 4. 打开本文件，点 Source（整份运行）。不要从中间逐行粘贴。
# 5. 不要把 17 个 GO 的基因合并后再打分
# ====================================================================
################################################################################

suppressPackageStartupMessages({
  cran_pkgs <- c("data.table", "readxl", "ggplot2", "survival")
  for (p in cran_pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
  library(data.table)
  library(readxl)
  library(ggplot2)
  library(survival)
})

# ---------------------------------------------------------------------------
# 可调参数
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
work_dir <- if (length(args) >= 1) args[[1]] else "E:/R/FUSC"
out_dir  <- if (length(args) >= 2) args[[2]] else "results_GO_excel"

# 强制切到 E:/R/FUSC；控制台单独跑 locate_fusc_excels(".") 时也要先 setwd
if (dir.exists(work_dir)) {
  setwd(work_dir)
} else {
  message("未找到 ", work_dir, " ，改用当前工作目录：", getwd())
  work_dir <- getwd()
}
message("当前工作目录：", getwd())

min_pathway_genes <- 2
min_clin_n <- 6
min_group_n <- 2
min_surv_n <- 10
min_surv_events <- 3
fdr_cutoff <- 0.05
rho_cutoff <- 0
neg_r_strict <- -0.15

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

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "per_GO"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "negcorr"), recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path(work_dir, "cache_GO")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 工具函数
# ---------------------------------------------------------------------------
first_present <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (!length(hit)) {
    low <- tolower(nms)
    for (cand in candidates) {
      i <- which(low == tolower(cand))
      if (length(i)) return(nms[[i[[1]]]])
    }
    return(NA_character_)
  }
  hit[[1]]
}

pick_col <- function(df, candidates) {
  hit <- first_present(names(df), candidates)
  if (is.na(hit)) NULL else hit
}

normalize_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\s+", "", x)
  x <- gsub("\\.", "-", x)
  x
}

safe_name <- function(x) gsub("[:/\\\\]+", "_", x)

is_missing_label <- function(x) {
  x <- trimws(tolower(as.character(x)))
  x %in% c("", "na", "nan", "null", "none", ".", "-", "unknown", "not reported",
           "[not available]", "[unknown]", "not available", "不确定", "未知")
}

looks_gene <- function(x) {
  x <- as.character(x)
  mean(grepl("^(ENSG|NM_|NR_|[A-Za-z][A-Za-z0-9._-]{1,})$", x)) > 0.5
}

looks_sample <- function(x) {
  x <- as.character(x)
  mean(grepl("TCGA|^GSM|^SAMN|^patient|^sample|^P[0-9]|[-_][0-9]", x, ignore.case = TRUE)) > 0.3 ||
    mean(grepl("^[A-Za-z0-9._-]{4,}$", x)) > 0.8
}

numeric_frac <- function(x) {
  if (is.numeric(x)) return(mean(is.finite(x)))
  y <- suppressWarnings(as.numeric(as.character(x)))
  mean(is.finite(y))
}

read_excel_dt <- function(path) {
  sheets <- tryCatch(excel_sheets(path), error = function(e) "Sheet1")
  dt <- as.data.table(read_excel(path, sheet = sheets[[1]], .name_repair = "minimal"))
  # 去掉全空列
  keep <- vapply(dt, function(col) any(!is.na(col) & as.character(col) != ""), logical(1))
  dt[, keep, with = FALSE]
}

list_excel <- function(dir) {
  hits <- list.files(dir, pattern = "\\.(xlsx|xls)$", full.names = TRUE, ignore.case = TRUE)
  hits <- hits[!grepl("(~\\$)|(^\\.)", basename(hits))]
  hits
}

# 按文件名前缀定位（用户给的名字可以没有扩展名）
find_named_excel <- function(dirs, stem, fuzzy = NULL) {
  dirs <- unique(dirs[dir.exists(dirs)])
  for (dir in dirs) {
    hits <- list.files(dir, full.names = TRUE, ignore.case = TRUE)
    hits <- hits[!startsWith(basename(hits), "~$")]
    if (!length(hits)) next
    bn <- basename(hits)
    stem_l <- tolower(stem)
    sans <- tolower(sub("\\.(xlsx|xls)$", "", bn, ignore.case = TRUE))
    exact <- which(sans == stem_l)
    if (length(exact)) return(hits[[exact[[1]]]])
    pref <- which(startsWith(sans, stem_l) & grepl("\\.(xlsx|xls)$", bn, ignore.case = TRUE))
    if (length(pref)) return(hits[[pref[[1]]]])
    if (!is.null(fuzzy)) {
      fz <- which(grepl(fuzzy, bn, ignore.case = TRUE) & grepl("\\.(xlsx|xls)$", bn, ignore.case = TRUE))
      if (length(fz)) return(hits[[fz[[1]]]])
    }
  }
  NULL
}

locate_fusc_excels <- function(dirs = c(getwd(), "E:/R/FUSC")) {
  dirs <- unique(c(dirs, "E:/R/FUSC", "E:\\R\\FUSC"))
  clin_path <- find_named_excel(dirs, "OEP00000155_样本_Human_1786804447300", fuzzy = "OEP00000155")
  fpkm_path <- find_named_excel(
    dirs,
    "FUSCCTNBC_Expression_RNAseqFPKM",
    fuzzy = "FUSCCTNBC.*FPKM|RNAseqFPKM"
  )
  if (is.null(clin_path) || is.null(fpkm_path)) {
    found <- unique(unlist(lapply(dirs[dir.exists(dirs)], function(d) {
      paste0(normalizePath(d, winslash = "/", mustWork = FALSE), ": ",
             paste(basename(list_excel(d)), collapse = ", "))
    })))
    stop(
      "未找到指定 Excel。R 当前工作目录是：", getwd(), "\n",
      "请先运行：setwd(\"E:/R/FUSC\")\n",
      "需要：OEP00000155_样本_Human_1786804447300.xlsx\n",
      "      FUSCCTNBC_Expression_RNAseqFPKM.xlsx\n",
      "各目录中的 Excel：\n  ", paste(found, collapse = "\n  ")
    )
  }
  list(clin = clin_path, fpkm = fpkm_path)
}

excel_to_expr <- function(dt) {
  is_num_col <- vapply(dt, function(x) numeric_frac(x) >= 0.6, logical(1))
  ann_idx <- which(!is_num_col)
  num_idx <- which(is_num_col)
  if (!length(num_idx)) stop("FPKM 表里没有数值列，请确认打开的是 FUSCCTNBC_Expression_RNAseqFPKM")

  # 前几列注释 + 后面样本：取 gene/symbol 列作行名
  if (length(ann_idx) && min(ann_idx) == 1 && min(num_idx) > 1 &&
      (looks_gene(as.character(dt[[ann_idx[1]]])) || nrow(dt) >= length(num_idx))) {
    gene_col <- pick_col(dt, c("gene", "Gene", "symbol", "Symbol", "gene_name", "GeneSymbol", "SYMBOL", names(dt)[ann_idx]))
    if (is.null(gene_col)) gene_col <- names(dt)[ann_idx[length(ann_idx)]]
    mat <- as.matrix(data.frame(lapply(dt[, num_idx, with = FALSE], function(x) {
      suppressWarnings(as.numeric(as.character(x)))
    })))
    rownames(mat) <- make.unique(as.character(dt[[gene_col]]))
    colnames(mat) <- normalize_id(names(dt)[num_idx])
    return(mat)
  }
  # 样本在行、基因在列
  first <- names(dt)[[1]]
  other <- dt[, -1, with = FALSE]
  mat <- as.matrix(data.frame(lapply(other, function(x) suppressWarnings(as.numeric(as.character(x))))))
  rownames(mat) <- make.unique(normalize_id(as.character(dt[[first]])))
  colnames(mat) <- as.character(names(other))
  t(mat)
}

drop_normal_samples <- function(mat) {
  cn <- colnames(mat)
  is_n <- grepl("normal|adjacent|(^|[._-])N$|(^|[._-])N[0-9]*$", cn, ignore.case = TRUE) &
    !grepl("TNBC|BLIS|LAR|IM|MES", cn, ignore.case = TRUE)
  if (sum(!is_n) >= min_clin_n && sum(is_n) > 0) {
    message("去掉疑似癌旁/正常样本：", sum(is_n), "；保留肿瘤：", sum(!is_n))
    mat <- mat[, !is_n, drop = FALSE]
  }
  mat
}

compact_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

match_expr_clin <- function(expr, clin) {
  exact <- intersect(colnames(expr), clin$sample_id)
  if (length(exact) >= min_clin_n) {
    return(list(expr = expr[, exact, drop = FALSE], clin = clin[match(exact, sample_id)]))
  }
  emap <- data.table(expr_id = colnames(expr), key_id = compact_id(colnames(expr)))
  cmap <- data.table(sample_id = clin$sample_id, key_id = compact_id(clin$sample_id))
  emap <- emap[key_id != "" & !duplicated(key_id)]
  cmap <- cmap[key_id != "" & !duplicated(key_id)]
  common_key <- intersect(emap$key_id, cmap$key_id)
  if (length(common_key) >= min_clin_n) {
    expr2 <- expr[, emap$expr_id[match(common_key, emap$key_id)], drop = FALSE]
    colnames(expr2) <- cmap$sample_id[match(common_key, cmap$key_id)]
    return(list(expr = expr2, clin = clin[match(colnames(expr2), sample_id)]))
  }
  list(expr = expr, clin = clin)
}

collapse_duplicate_genes <- function(mat) {
  rn <- sub("\\..*$", "", rownames(mat))
  if (!any(duplicated(rn))) {
    rownames(mat) <- rn
    return(mat)
  }
  means <- rowMeans(mat, na.rm = TRUE)
  ord <- order(means, decreasing = TRUE)
  mat <- mat[ord, , drop = FALSE]
  rn <- sub("\\..*$", "", rownames(mat))
  mat <- mat[!duplicated(rn), , drop = FALSE]
  rownames(mat) <- sub("\\..*$", "", rownames(mat))
  mat
}

ihc_pos <- function(x) {
  x <- tolower(as.character(x))
  grepl("pos|\\+|阳性|yes|true|^1$|3\\+|2\\+", x) & !grepl("neg|阴性|equivocal|不确定", x)
}

ihc_neg <- function(x) {
  x <- tolower(as.character(x))
  grepl("neg|阴性|\\b0\\b|^-$", x)
}

derive_ihc_subtype <- function(clin) {
  er_col <- pick_col(clin, c("ER", "er_status", "ER.Status", "er_status_by_ihc", "雌激素受体", "ER状态"))
  pr_col <- pick_col(clin, c("PR", "pr_status", "PR.Status", "pr_status_by_ihc", "孕激素受体", "PR状态"))
  her2_col <- pick_col(clin, c("HER2", "her2_status", "HER2.Final.Status", "her2_status_by_ihc", "HER2状态", "Her2"))
  if (is.null(her2_col) && is.null(er_col)) return(rep(NA_character_, nrow(clin)))
  er <- if (is.null(er_col)) rep(NA, nrow(clin)) else clin[[er_col]]
  pr <- if (is.null(pr_col)) rep(NA, nrow(clin)) else clin[[pr_col]]
  her2 <- if (is.null(her2_col)) rep(NA, nrow(clin)) else clin[[her2_col]]
  hr_pos <- ihc_pos(er) | ihc_pos(pr)
  hr_neg <- ihc_neg(er) & (is.null(pr_col) | ihc_neg(pr))
  h2_pos <- ihc_pos(her2)
  h2_neg <- ihc_neg(her2)
  out <- rep(NA_character_, nrow(clin))
  out[hr_pos & h2_neg] <- "Luminal (HR+/HER2-)"
  out[hr_pos & h2_pos] <- "Luminal-HER2 (HR+/HER2+)"
  out[hr_neg & h2_pos] <- "HER2+"
  out[(hr_neg | (!hr_pos & ihc_neg(er))) & h2_neg] <- "TNBC"
  out
}

normalize_subtype_label <- function(x) {
  raw <- trimws(as.character(x))
  raw[is_missing_label(raw)] <- NA_character_
  out <- rep(NA_character_, length(raw))
  out[grepl("\\bLAR\\b|luminal androgen|腔面雄激素", raw, ignore.case = TRUE)] <- "LAR"
  out[is.na(out) & grepl("\\bBLIS\\b|basal-like immune|基底样免疫抑制", raw, ignore.case = TRUE)] <- "BLIS"
  out[is.na(out) & grepl("\\bMES\\b|mesenchymal|间充质", raw, ignore.case = TRUE)] <- "MES"
  out[is.na(out) & grepl("\\bIM\\b|immunomodulatory|免疫调节", raw, ignore.case = TRUE)] <- "IM"
  out[is.na(out) & grepl("三阴|TNBC", raw, ignore.case = TRUE)] <- "TNBC"
  out[is.na(out) & grepl("HER2\\s*\\+|HER2阳性|Her2-?enr|HER2_pos", raw, ignore.case = TRUE) &
        !grepl("luminal|腔面|HR\\+|LAR", raw, ignore.case = TRUE)] <- "HER2+"
  out[is.na(out) & grepl("LumA|Luminal\\s*A|腔面\\s*A", raw, ignore.case = TRUE)] <- "Luminal A"
  out[is.na(out) & grepl("LumB|Luminal\\s*B|腔面\\s*B", raw, ignore.case = TRUE)] <- "Luminal B"
  out[is.na(out) & grepl("Luminal|腔面|HR\\+", raw, ignore.case = TRUE)] <- "Luminal"
  out[is.na(out) & !is.na(raw)] <- raw[is.na(out) & !is.na(raw)]
  out
}

guess_subtype_column <- function(clin) {
  named <- pick_col(clin, c(
    "subtype", "Subtype", "PAM50", "pam50", "BRCA_Subtype_PAM50",
    "molecular_subtype", "IHC_subtype", "Fudan_subtype", "TNBC_subtype",
    "mRNA_subtype", "分型", "亚型", "分子分型", "乳腺癌分型", "转录组分型"
  ))
  if (!is.null(named)) return(named)
  rx <- "LAR|BLIS|\\bIM\\b|\\bMES\\b|三阴|TNBC|HER2|Luminal|腔面"
  char_cols <- names(clin)[vapply(clin, function(x) is.character(x) || is.factor(x), logical(1))]
  char_cols <- setdiff(char_cols, c("sample_id"))
  best <- NULL
  best_n <- 0
  for (cc in char_cols) {
    n <- sum(grepl(rx, as.character(clin[[cc]])), na.rm = TRUE)
    if (n > best_n) {
      best_n <- n
      best <- cc
    }
  }
  if (!is.null(best) && best_n >= min_clin_n) best else NULL
}

parse_event <- function(ev) {
  if (is.numeric(ev)) return(as.integer(ev > 0))
  as.integer(grepl("dead|deceased|death|event|relapse|recur|1|是|死亡|复发|进展", as.character(ev), ignore.case = TRUE))
}

prepare_clinical <- function(clin) {
  id_col <- pick_col(clin, c(
    "sample", "sampleID", "Sample", "SampleName", "sample_name", "patient", "Patient",
    "ID", "id", "bcr_patient_barcode", "submitter_id",
    "样本", "样本名称", "样品名称", "样本编号", "病人", "病人编号", names(clin)[1]
  ))
  clin[, sample_id := normalize_id(get(id_col))]
  clin <- clin[!is.na(sample_id) & sample_id != "" & !duplicated(sample_id)]

  sub_col <- guess_subtype_column(clin)
  if (!is.null(sub_col)) {
    message("分型列：", sub_col)
    clin[, subtype := normalize_subtype_label(get(sub_col))]
  } else {
    clin[, subtype := NA_character_]
  }
  derived <- derive_ihc_subtype(clin)
  clin[, subtype := fifelse(is.na(subtype) | subtype == "", derived, subtype)]

  time_col <- pick_col(clin, c(
    "OS.time", "OS_time", "os_time", "_OS", "overall_survival_time", "OS_months",
    "生存时间", "总生存时间", "OS时间", "随访时间", "days_to_death", "days_to_last_followup"
  ))
  event_col <- pick_col(clin, c(
    "OS", "OS.event", "os_event", "_EVENT", "event", "status", "vital_status",
    "生存状态", "死亡", "预后", "结局", "死亡状态"
  ))
  if (!is.null(time_col)) clin[, os_time := as.numeric(get(time_col))]
  if (!is.null(event_col)) clin[, os_event := parse_event(get(event_col))]

  rfs_time_col <- pick_col(clin, c(
    "RFS.time", "RFS_time", "DFS.time", "DFS_time", "PFS.time", "DFI.time",
    "无病生存时间", "无复发生存时间", "复发时间"
  ))
  rfs_event_col <- pick_col(clin, c(
    "RFS", "RFS.event", "DFS", "DFS.event", "PFS", "复发", "复发生存状态", "无病生存状态"
  ))
  if (!is.null(rfs_time_col)) clin[, rfs_time := as.numeric(get(rfs_time_col))]
  if (!is.null(rfs_event_col)) clin[, rfs_event := parse_event(get(rfs_event_col))]
  clin
}

zscore_one_go <- function(expr_mat, genes) {
  g <- intersect(unique(genes), rownames(expr_mat))
  if (length(g) < min_pathway_genes) return(NULL)
  sub <- expr_mat[g, , drop = FALSE]
  z <- t(scale(t(sub)))
  z[!is.finite(z)] <- NA
  sc <- colMeans(z, na.rm = TRUE)
  sc[!is.finite(sc)] <- NA
  attr(sc, "n_genes") <- length(g)
  attr(sc, "genes") <- g
  sc
}

two_group_delta <- function(score, group, pos, neg) {
  ok <- is.finite(score) & group %in% c(pos, neg)
  if (sum(group[ok] == pos) < min_group_n || sum(group[ok] == neg) < min_group_n) return(NULL)
  y <- score[ok]
  g <- group[ok]
  wt <- suppressWarnings(wilcox.test(y[g == pos], y[g == neg]))
  data.table(
    n = sum(ok),
    n_pos = sum(g == pos),
    n_neg = sum(g == neg),
    effect = stats::median(y[g == pos], na.rm = TRUE) - stats::median(y[g == neg], na.rm = TRUE),
    pvalue = wt$p.value,
    method = "Wilcoxon"
  )
}

cox_high_low <- function(score, time, event) {
  dt <- data.table(score = as.numeric(score), time = as.numeric(time), event = as.integer(event))
  dt <- dt[is.finite(score) & is.finite(time) & is.finite(event) & time > 0]
  if (nrow(dt) < min_surv_n || sum(dt$event == 1) < min_surv_events) return(NULL)
  dt[, grp := factor(
    fifelse(score >= median(score, na.rm = TRUE), "High", "Low"),
    levels = c("Low", "High")
  )]
  if (uniqueN(dt$grp) < 2) return(NULL)
  fit <- tryCatch(coxph(Surv(time, event) ~ grp, data = dt), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  hr <- unname(s$coefficients[1, "exp(coef)"])
  p <- unname(s$coefficients[1, "Pr(>|z|)"])
  data.table(n = nrow(dt), n_event = sum(dt$event == 1), HR = hr, effect = log2(pmax(hr, 1e-8)), pvalue = p)
}

make_bubble <- function(d, title, subtitle, caption, color_name, outfile_stub) {
  d <- as.data.table(d)
  d <- d[is.finite(pvalue) & is.finite(effect)]
  if (!nrow(d)) {
    message("没有可画的点：", title)
    return(invisible(NULL))
  }
  d[, neglogp := pmin(10, -log10(pmax(pvalue, 1e-12)))]
  d[, go_lab := paste0(GO, "  ", GO_name)]
  go_rank <- d[, .(m = max(neglogp, na.rm = TRUE)), by = go_lab]
  setorder(go_rank, m)
  d[, go_lab := factor(go_lab, levels = go_rank$go_lab)]
  if (!is.factor(d$feature_lab)) {
    d[, feature_lab := factor(feature_lab, levels = unique(as.character(feature_lab)))]
  }
  lim <- max(0.4, as.numeric(stats::quantile(abs(d$effect), 0.98, na.rm = TRUE)), na.rm = TRUE)
  p <- ggplot(d, aes(x = feature_lab, y = go_lab, size = neglogp, color = effect)) +
    geom_point(alpha = 0.9) +
    scale_size_continuous(range = c(2.5, 12), name = expression(-log[10](p))) +
    scale_color_gradient2(
      low = "#3C5488", mid = "grey90", high = "#E64B35",
      midpoint = 0, limits = c(-lim, lim), name = color_name
    ) +
    labs(title = title, subtitle = subtitle, caption = caption, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      panel.grid.major.x = element_blank(),
      plot.caption = element_text(hjust = 0, size = 8, color = "grey30")
    )
  ht <- max(6.5, min(14, 0.38 * uniqueN(d$GO) + 3.2))
  ggsave(paste0(outfile_stub, ".pdf"), p, width = 10, height = ht)
  ggsave(paste0(outfile_stub, ".png"), p, width = 10, height = ht, dpi = 150)
  message("已保存气泡图：", outfile_stub)
  invisible(p)
}

read_local_genesets <- function(path, wanted) {
  if (grepl("\\.(xlsx|xls)$", path, ignore.case = TRUE)) {
    dt <- as.data.table(read_excel(path, .name_repair = "minimal"))
  } else {
    dt <- fread(path)
  }
  go_col <- pick_col(dt, c("GO_ID", "go_id", "GO", "term"))
  sym_col <- pick_col(dt, c("symbol", "gene_name", "Symbol", "gene", "SYMBOL"))
  if (is.null(go_col) || is.null(sym_col)) stop(path, " 需要 GO_ID 与 symbol 两列")
  raw <- split(as.character(dt[[sym_col]]), as.character(dt[[go_col]]))
  out <- lapply(wanted, function(g) unique(na.omit(as.character(raw[[g]]))))
  names(out) <- wanted
  out
}

fetch_quickgo_one <- function(go_id) {
  q <- paste0(
    "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch",
    "?goId=", utils::URLencode(go_id, reserved = TRUE),
    "&taxonId=9606&goUsage=descendants&geneProductType=protein"
  )
  tmp <- tempfile(fileext = ".tsv")
  ok <- FALSE
  if (nzchar(Sys.which("curl"))) {
    st <- suppressWarnings(system2("curl", c("-sL", "-H", "Accept: text/tsv", "--fail", "-o", tmp, q),
                                   stdout = FALSE, stderr = FALSE))
    ok <- identical(st, 0L) && file.exists(tmp) && isTRUE(file.info(tmp)$size > 20)
  }
  if (!ok) return(character(0))
  dt <- tryCatch(fread(tmp, sep = "\t"), error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(character(0))
  sym_col <- pick_col(dt, c("SYMBOL", "symbol", "GENE PRODUCT ID"))
  if (is.null(sym_col)) return(character(0))
  unique(na.omit(as.character(dt[[sym_col]])))
}

load_go_genesets <- function(wanted) {
  local_hits <- list.files(
    ".",
    pattern = "go_genesets\\.(tsv|txt|csv|xlsx|xls|gmt)$",
    full.names = TRUE, ignore.case = TRUE
  )
  if (length(local_hits)) {
    message("使用本地基因集：", local_hits[[1]])
    return(read_local_genesets(local_hits[[1]], wanted))
  }
  message("从 QuickGO 逐个下载 GO 基因（含子术语）")
  out <- lapply(wanted, function(x) character(0))
  names(out) <- wanted
  n_ok <- 0L
  for (go_id in wanted) {
    cache <- file.path(cache_dir, paste0(safe_name(go_id), "_genes.txt"))
    if (file.exists(cache) && isTRUE(file.info(cache)$size > 0)) {
      out[[go_id]] <- unique(readLines(cache, warn = FALSE))
    } else {
      syms <- tryCatch(fetch_quickgo_one(go_id), error = function(e) character(0))
      out[[go_id]] <- setdiff(unique(syms), c("", "-", "NA"))
      if (length(out[[go_id]])) writeLines(out[[go_id]], cache)
    }
    if (length(out[[go_id]])) n_ok <- n_ok + 1L
    message("  ", go_id, ": ", length(out[[go_id]]), " 个基因")
  }
  if (!n_ok) {
    stop("无法获取 GO 基因。请在 ", getwd(),
         " 放置 go_genesets.tsv（列：GO_ID, symbol）后重跑。")
  }
  out
}

spearman_vs_score <- function(expr_mat, score) {
  common <- intersect(colnames(expr_mat), names(score))
  mat <- expr_mat[, common, drop = FALSE]
  sc <- score[common]
  keep <- apply(mat, 1, function(v) sd(v, na.rm = TRUE) > 0)
  mat <- mat[keep, , drop = FALSE]
  n <- ncol(mat)
  rho <- as.numeric(cor(t(mat), sc, method = "spearman", use = "pairwise.complete.obs"))
  rho <- pmin(pmax(rho, -0.999999), 0.999999)
  tstat <- rho * sqrt((n - 2) / pmax(1e-12, 1 - rho^2))
  pval <- 2 * pt(-abs(tstat), df = n - 2)
  data.table(
    gene = rownames(mat),
    rho = rho,
    pvalue = pval,
    fdr = p.adjust(pval, method = "BH"),
    n = n
  )
}

# ---------------------------------------------------------------------------
# 读入指定的两份 Excel
# ---------------------------------------------------------------------------
pair <- locate_fusc_excels(c(getwd(), work_dir, "E:/R/FUSC"))
message("临床 Excel：", pair$clin)
message("FPKM Excel：", pair$fpkm)

clin_raw <- read_excel_dt(pair$clin)
fpkm_raw <- read_excel_dt(pair$fpkm)
clin <- prepare_clinical(clin_raw)
expr <- excel_to_expr(fpkm_raw)
expr <- collapse_duplicate_genes(expr)
colnames(expr) <- normalize_id(colnames(expr))
expr <- expr[, !duplicated(colnames(expr)), drop = FALSE]
expr <- drop_normal_samples(expr)

mx <- suppressWarnings(max(expr, na.rm = TRUE))
if (is.finite(mx) && mx > 50) {
  message("检测到原始 FPKM（max=", round(mx, 2), "），进行 log2(x+1)")
  expr <- log2(expr + 1)
}

aligned <- match_expr_clin(expr, clin)
expr <- aligned$expr
clin <- aligned$clin
common <- intersect(colnames(expr), clin$sample_id)

if (length(common) < min_clin_n) {
  stop("临床与 FPKM 匹配到的样本太少（n=", length(common),
       "）。请确认 OEP00000155 样本表与 FUSCCTNBC_Expression_RNAseqFPKM 使用同一套样本 ID。")
}
expr <- expr[, common, drop = FALSE]
clin <- clin[match(common, sample_id)]
message("匹配样本数：", length(common),
        "；分型非空：", sum(!is.na(clin$subtype)),
        "；有 OS：", if ("os_time" %in% names(clin)) sum(is.finite(clin$os_time)) else 0,
        "；有 RFS：", if ("rfs_time" %in% names(clin)) sum(is.finite(clin$rfs_time)) else 0)

# ---------------------------------------------------------------------------
# GO 基因集 + 每个通路单独打分
# ---------------------------------------------------------------------------
go_genesets <- load_go_genesets(go_ids)
go_sizes <- rbindlist(lapply(go_ids, function(g) {
  data.table(
    GO = g, GO_name = unname(go_name_map[g]),
    n_annotated = length(go_genesets[[g]]),
    n_in_expr = sum(go_genesets[[g]] %in% rownames(expr))
  )
}))
fwrite(go_sizes, file.path(out_dir, "00_GO_gene_set_sizes.csv"))

usable <- go_sizes$GO[go_sizes$n_in_expr >= min_pathway_genes]
skipped <- setdiff(go_ids, usable)
if (length(skipped)) {
  message("以下 GO 在表达矩阵中基因 < ", min_pathway_genes, "，已跳过（未合并到其他通路）：",
          paste(skipped, collapse = ", "))
}
if (!length(usable)) {
  stop("没有任何 GO 在表达矩阵里有足够基因。请检查基因名是否为 SYMBOL，或提供 go_genesets.xlsx/tsv。")
}

all_scores <- list()
subtype_rows <- list()
surv_rows <- list()
neg_summary <- list()

for (go_id in usable) {
  go_title <- unname(go_name_map[go_id])
  message("\n========== ", go_id, " | ", go_title, " ==========")
  go_dir <- file.path(out_dir, "per_GO", safe_name(go_id))
  dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)

  go_score <- zscore_one_go(expr, go_genesets[[go_id]])
  if (is.null(go_score)) {
    writeLines("SKIPPED: too few genes", file.path(go_dir, "SKIPPED.txt"))
    next
  }
  all_scores[[go_id]] <- go_score
  fwrite(
    data.table(sample = names(go_score), pathway_score = as.numeric(go_score)),
    file.path(go_dir, "pathway_score.csv")
  )

  # ---- 分型：每个亚型 vs 其余 ----
  if ("subtype" %in% names(clin) && sum(!is.na(clin$subtype)) >= min_clin_n) {
    grp <- factor(clin$subtype)
    y <- go_score[clin$sample_id]
    ok <- !is.na(grp) & is.finite(y)
    if (sum(ok) >= min_clin_n && nlevels(droplevels(grp[ok])) >= 2) {
      kt <- kruskal.test(y[ok] ~ grp[ok])
      fwrite(
        data.table(GO = go_id, GO_name = go_title, method = "Kruskal-Wallis",
                   pvalue = kt$p.value, n = sum(ok)),
        file.path(go_dir, "subtype_overall.csv")
      )
      for (lv in levels(droplevels(grp[ok]))) {
        flag <- ifelse(as.character(grp) == lv, lv, "Other")
        hit <- two_group_delta(y, flag, lv, "Other")
        if (is.null(hit)) next
        hit[, `:=`(GO = go_id, GO_name = go_title, feature = lv,
                   feature_lab = paste0(lv, " vs 其余"))]
        subtype_rows[[length(subtype_rows) + 1]] <- hit
      }
    }
  }

  # ---- 预后：OS / RFS High vs Low Cox ----
  if (all(c("os_time", "os_event") %in% names(clin))) {
    cox <- cox_high_low(go_score[clin$sample_id], clin$os_time, clin$os_event)
    if (!is.null(cox)) {
      cox[, `:=`(GO = go_id, GO_name = go_title, feature = "OS", feature_lab = "OS 总生存")]
      surv_rows[[length(surv_rows) + 1]] <- cox
      fwrite(cox, file.path(go_dir, "survival_OS.csv"))
    }
  }
  if (all(c("rfs_time", "rfs_event") %in% names(clin))) {
    cox <- cox_high_low(go_score[clin$sample_id], clin$rfs_time, clin$rfs_event)
    if (!is.null(cox)) {
      cox[, `:=`(GO = go_id, GO_name = go_title, feature = "RFS", feature_lab = "RFS 无复发生存")]
      surv_rows[[length(surv_rows) + 1]] <- cox
      fwrite(cox, file.path(go_dir, "survival_RFS.csv"))
    }
  }

  # ---- 负相关基因（仅本通路分数） ----
  tab <- spearman_vs_score(expr, go_score)
  tab[, `:=`(GO = go_id, in_this_GO = gene %in% go_genesets[[go_id]])]
  neg <- tab[is.finite(rho) & rho < rho_cutoff & fdr < fdr_cutoff][order(rho)]
  fwrite(tab, file.path(out_dir, "negcorr", paste0(safe_name(go_id), "_all_genes.csv")))
  fwrite(neg, file.path(out_dir, "negcorr", paste0(safe_name(go_id), "_negative_genes.csv")))
  fwrite(neg[rho <= neg_r_strict], file.path(go_dir, "negative_genes_strict.csv"))
  neg_summary[[go_id]] <- data.table(
    GO = go_id, GO_name = go_title,
    n_tested = nrow(tab), n_negative_fdr = nrow(neg),
    min_rho = if (nrow(neg)) min(neg$rho) else NA_real_
  )
  message("  负相关基因：", nrow(neg), "（FDR < ", fdr_cutoff, " 且 rho < 0）")
}

# ---------------------------------------------------------------------------
# 汇总 + 气泡图
# ---------------------------------------------------------------------------
if (length(all_scores)) {
  score_dt <- data.table(sample = names(all_scores[[1]]))
  for (g in names(all_scores)) score_dt[, (g) := all_scores[[g]][sample]]
  fwrite(score_dt, file.path(out_dir, "01_pathway_scores_each_GO.csv"))
}

if (length(subtype_rows)) {
  sub_all <- rbindlist(subtype_rows, fill = TRUE)
  fwrite(sub_all, file.path(out_dir, "02_subtype_each_GO.csv"))
  make_bubble(
    sub_all,
    title = "各 GO 通路与乳腺癌分型",
    subtitle = "每个 GO 单独打分；该亚型 vs 其余样本的通路分数差（Fudan TNBC 或 IHC 分型，未合并通路）",
    caption = "颜色：Δscore = 该亚型中位数 − 其余中位数。红=该亚型通路活性更高；蓝相反。点越大 p 越小。",
    color_name = "Δscore",
    outfile_stub = file.path(out_dir, "02_subtype_bubble_each_GO")
  )
} else {
  message("未检测到可用的分型列（subtype / 分型 / 三阴 / HER2 等），跳过分型气泡图")
}

if (length(surv_rows)) {
  surv_all <- rbindlist(surv_rows, fill = TRUE)
  fwrite(surv_all, file.path(out_dir, "03_prognosis_each_GO.csv"))
  surv_all[, feature_lab := factor(feature_lab, levels = intersect(c("OS 总生存", "RFS 无复发生存"), unique(feature_lab)))]
  make_bubble(
    surv_all,
    title = "各 GO 通路与乳腺癌预后",
    subtitle = "每个 GO 单独打分；High vs Low Cox（未合并通路）",
    caption = "颜色：log2(HR)。红=通路高分组预后更差（HR>1）；蓝=高分组预后更好。点越大 p 越小。",
    color_name = "log2(HR)",
    outfile_stub = file.path(out_dir, "03_prognosis_bubble_each_GO")
  )
} else {
  message("临床表中没有可用的生存时间/结局列，跳过预后气泡图。可在 Excel 中增加 OS.time/OS 或 RFS.time/RFS 后重跑。")
}

if (length(neg_summary)) {
  fwrite(rbindlist(neg_summary), file.path(out_dir, "negcorr", "04_negative_gene_counts_per_GO.csv"))
}

message("完成。结果目录：", normalizePath(out_dir, mustWork = FALSE))
message("分型气泡图：02_subtype_bubble_each_GO.pdf")
message("预后气泡图：03_prognosis_bubble_each_GO.pdf")
message("各通路负相关基因：negcorr/*_negative_genes.csv")
