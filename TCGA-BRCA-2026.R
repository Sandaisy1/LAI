################################################################################
# TCGA-BRCA-2026.R
# 单独脚本：放到 E:/R/TCGA-BRCA-2026，RStudio 打开后从第一行 Source
# 不要从中间粘贴，不要运行旧的 factor(meta_M) 行
#
# 数据目录（.tsv / .tsv.gz 均可；空的 .tsv 会自动跳过，改读 .gz）：
#   TCGA-BRCA.star_fpkm.tsv(.gz)
#   TCGA-BRCA.clinical.tsv(.gz)
#   TCGA-BRCA.survival.tsv(.gz)
#   gencode.v36.annotation.gtf.gene.probemap
#
# 神经浸润：只分析下面 5 个神经信号 GO，每个单独取基因、单独打分（不合并）
#   GO:0019227  neuronal action potential propagation
#   GO:1902847  regulation of neuronal signal transduction
#   GO:0097374  sensory neuron axon guidance
#   GO:1902667  regulation of axon guidance
#   GO:0007409  axonogenesis
# 原位肿瘤：Primary Tumor / 条形码 01
#
# 1) 神经浸润 vs 转移：主图是气泡图（纵轴=各神经 GO，横轴=转移定义）
#    a. 原位肿瘤 诊断时远处转移 M1 vs M0
#    b. 原位肿瘤 AJCC Stage IV vs I–III
#    c. 原位肿瘤 淋巴结 N+ vs N0
#    d. 转移组织 vs 原位肿瘤
# 2) 原位肿瘤内与转移负相关的基因（标准同 1a/1b/1c）
# 3) 原位肿瘤内与每个神经 GO 通路分数负相关的基因（按 GO 分开）
#
# 负相关：Spearman r < 0 且 p < 0.05；另存 r <= -0.15 且 p < 0.05
################################################################################

library(data.table)
library(ggplot2)
library(ggpubr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ==============================================================================
# 参数
# ==============================================================================
work_dir <- Sys.getenv("TCGA_BRCA_2026_DIR", unset = "E:/R/TCGA-BRCA-2026")
out_dir <- "results_TCGA-BRCA-2026_5GO"
min_pathway_genes <- 1
min_group_n <- 2
neg_pvalue_cutoff <- 0.05
neg_r_cutoff <- 0
strict_r_cutoff <- -0.15
min_expr_frac <- 0.20

go_list <- c(
  "GO:0019227",
  "GO:1902847",
  "GO:0097374",
  "GO:1902667",
  "GO:0007409"
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
  "GO:0007409" = "axonogenesis",
  "GO:0036518" = "chemorepulsion of dopaminergic neuron axon"
)

# ==============================================================================
# 工具函数（全部先定义，读完表达矩阵再分析）
# ==============================================================================
normalize_barcode <- function(x) {
  x <- toupper(gsub("\\.", "-", as.character(x)))
  x <- sub("A$", "", x)
  ifelse(nchar(x) >= 15, substr(x, 1, 15), x)
}
patient_id <- function(x) {
  x <- normalize_barcode(x)
  ifelse(nchar(x) >= 12, substr(x, 1, 12), x)
}
sample_type_code <- function(x) {
  x <- normalize_barcode(x)
  ifelse(nchar(x) >= 15, substr(x, 14, 15), NA_character_)
}
first_present <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) NA_character_ else hit[1]
}
safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}
go_title <- function(go_id) {
  go_id <- as.character(go_id)
  out <- unname(go_name_map[go_id])
  miss <- is.na(out) | !nzchar(out)
  out[miss] <- go_id[miss]
  out
}
go_lab <- function(go_id) paste0(as.character(go_id), "  ", go_title(go_id))

existing_files <- function(stems) {
  cands <- unique(unlist(lapply(stems, function(s) c(s, paste0(s, ".gz")))))
  cands[file.exists(cands) & !is.na(file.info(cands)$isdir) & !file.info(cands)$isdir]
}
probe_table_file <- function(f) {
  info <- file.info(f)
  if (is.null(info) || is.na(info$size) || info$size < 200) {
    return(data.table(file = f, size = 0, n_col = 0L, readable = FALSE))
  }
  hdr <- tryCatch(names(fread(f, nrows = 0, fill = TRUE, showProgress = FALSE)),
                  error = function(e) character())
  data.table(file = f, size = as.numeric(info$size), n_col = length(hdr),
             readable = length(hdr) >= 1L)
}
pick_best_table <- function(stems, must = TRUE, min_cols = 1L, label = "数据表") {
  hits <- existing_files(stems)
  if (length(hits) == 0) {
    if (must) stop("找不到", label, "：", paste(stems, collapse = " / "), "（.tsv 或 .tsv.gz）")
    return(NA_character_)
  }
  tab <- rbindlist(lapply(hits, probe_table_file), fill = TRUE)
  message(label, "候选：\n",
          paste(sprintf("  %s  大小=%.1fMB  列=%d  可读=%s",
                        tab$file, tab$size / 1024^2, tab$n_col, tab$readable),
                collapse = "\n"))
  ok <- tab[readable == TRUE & n_col >= min_cols]
  if (nrow(ok) == 0) {
    if (must) {
      stop(label, "存在但读不了。请删掉空的 .tsv，保留 .tsv.gz 后再 Source。")
    }
    return(NA_character_)
  }
  setorder(ok, -n_col, -size)
  ok$file[1]
}
safe_fread <- function(path, ...) {
  out <- tryCatch(fread(path, showProgress = TRUE, ...), error = function(e) e)
  if (!inherits(out, "error")) return(out)
  alt <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    sub("\\.gz$", "", path, ignore.case = TRUE)
  } else {
    paste0(path, ".gz")
  }
  if (file.exists(alt) && isTRUE(file.info(alt)$size > 200)) {
    message("读 ", path, " 失败（", conditionMessage(out), "），改读 ", alt)
    return(fread(alt, showProgress = TRUE, ...))
  }
  stop("读不了 ", path, "：", conditionMessage(out), "。空的 .tsv 请删除，改用 .tsv.gz。")
}
pick_best_clinical <- function() {
  f <- pick_best_table(
    c("TCGA-BRCA.clinical.tsv", "TCGA-BRCA.GDC_phenotype.tsv"),
    must = TRUE, min_cols = 5L, label = "临床表"
  )
  hdr <- names(fread(f, nrows = 0, fill = TRUE, showProgress = FALSE))
  n_hit <- sum(grepl("ajcc_pathologic_m|ajcc_pathologic_stage|pathologic_m|pathologic_stage",
                     hdr, ignore.case = TRUE))
  message("选用临床：", f, "  列=", length(hdr), "  分期/M列=", n_hit)
  f
}
pick_clin_col <- function(dt, patterns) {
  nms <- names(dt)
  skip <- grepl(
    "^(sample|patient|barcode|submitter|case_id|id|project|age_|days_|year_|uuid)",
    nms, ignore.case = TRUE
  )
  nms2 <- nms[!skip]
  for (p in patterns) {
    hit <- grep(p, nms2, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

simplify_stage <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("IV|STAGE.?4", x)] <- "Stage IV"
  out[is.na(out) & grepl("III|STAGE.?3", x)] <- "Stage III"
  out[is.na(out) & grepl("II|STAGE.?2", x)] <- "Stage II"
  out[is.na(out) & grepl("I|STAGE.?1", x)] <- "Stage I"
  out[grepl("\\bX\\b|NOT AVAILABLE|UNKNOWN|NOT REPORTED|NOT APPLICABLE", x)] <- NA_character_
  out
}
classify_m <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("M1[ABC]?|\\bM1\\b", x)] <- "M1"
  out[is.na(out) & grepl("\\bM0\\b|CM0", x)] <- "M0"
  out[grepl("MX|NOT AVAILABLE|UNKNOWN|NOT REPORTED", x)] <- NA_character_
  out
}
classify_n <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("N[1-3]", x)] <- "Nplus"
  out[is.na(out) & grepl("\\bN0\\b", x)] <- "N0"
  out[grepl("NX|NOT AVAILABLE|UNKNOWN|NOT REPORTED", x)] <- NA_character_
  out
}
classify_sample_type <- function(x, barcode) {
  out <- rep(NA_character_, length(x))
  xl <- toupper(as.character(x))
  out[grepl("PRIMARY", xl)] <- "原位肿瘤"
  out[grepl("METASTATIC", xl)] <- "转移组织"
  out[grepl("NORMAL", xl)] <- "正常"
  code <- sample_type_code(barcode)
  out[is.na(out) & code == "01"] <- "原位肿瘤"
  out[is.na(out) & code %in% c("06", "07")] <- "转移组织"
  out[is.na(out) & code == "11"] <- "正常"
  out
}

get_go_genes <- function(go_id) {
  pick <- function(keytype) {
    tryCatch(
      AnnotationDbi::select(
        org.Hs.eg.db, keys = go_id, keytype = keytype,
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
    if (length(eg) == 0) NULL else {
      AnnotationDbi::select(org.Hs.eg.db, keys = eg, keytype = "ENTREZID",
                            columns = c("SYMBOL", "ENSEMBL"))
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
  unique(res[!is.na(SYMBOL) & SYMBOL != "", .(GO, SYMBOL, ENSEMBL, ENTREZID)])
}

# 通路活性：基因 z-score 后对样本取均值。不要用 scale()/t()，变量不要叫 score
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
  go_score <- colMeans(z, na.rm = TRUE)
  names(go_score) <- colnames(sub)
  attr(go_score, "n_genes") <- length(genes)
  attr(go_score, "genes") <- genes
  go_score
}

spearman_vs_go_score <- function(mat, go_score_vec) {
  common <- intersect(colnames(mat), names(go_score_vec))
  if (length(common) < 5) return(data.table())
  mat <- mat[, common, drop = FALSE]
  go_score_vec <- go_score_vec[common]
  keep <- apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
  mat <- mat[keep, , drop = FALSE]
  n <- ncol(mat)
  r <- as.numeric(cor(
    base::t(as.matrix(mat)), go_score_vec,
    method = "spearman", use = "pairwise.complete.obs"
  ))
  names(r) <- rownames(mat)
  r <- pmin(pmax(r, -0.999999), 0.999999)
  tstat <- r * sqrt((n - 2) / pmax(1e-12, 1 - r^2))
  p <- 2 * stats::pt(-abs(tstat), df = n - 2)
  data.table(feature = names(r), spearman_r = r, pvalue = p,
             fdr = p.adjust(p, method = "BH"), n = n)
}

spearman_vs_binary <- function(mat, group, pos, neg) {
  g <- as.character(group)
  names(g) <- names(group)
  common <- intersect(colnames(mat), names(g))
  g <- g[common]
  keep_s <- g %in% c(pos, neg)
  if (sum(g[keep_s] == pos) < min_group_n || sum(g[keep_s] == neg) < min_group_n) {
    return(data.table())
  }
  y <- ifelse(g[keep_s] == pos, 1, 0)
  mat <- mat[, names(y), drop = FALSE]
  keep <- apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
  mat <- mat[keep, , drop = FALSE]
  n <- ncol(mat)
  r <- as.numeric(cor(
    base::t(as.matrix(mat)), y,
    method = "spearman", use = "pairwise.complete.obs"
  ))
  names(r) <- rownames(mat)
  r <- pmin(pmax(r, -0.999999), 0.999999)
  tstat <- r * sqrt((n - 2) / pmax(1e-12, 1 - r^2))
  p <- 2 * stats::pt(-abs(tstat), df = n - 2)
  data.table(
    feature = names(r), spearman_r = r, pvalue = p,
    fdr = p.adjust(p, method = "BH"), n = n,
    n_pos = sum(y == 1), n_neg = sum(y == 0),
    pos_level = pos, neg_level = neg
  )
}

compare_go_groups <- function(go_score_vec, group, pos, neg, grouping) {
  df <- data.frame(
    pathway_score = as.numeric(go_score_vec),
    group = as.character(group),
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$pathway_score) & df$group %in% c(pos, neg), ]
  n_pos <- sum(df$group == pos)
  n_neg <- sum(df$group == neg)
  if (n_pos < 1 && n_neg < 1) return(NULL)
  pval <- NA_real_
  if (n_pos >= min_group_n && n_neg >= min_group_n) {
    wt <- suppressWarnings(stats::wilcox.test(pathway_score ~ group, data = df))
    pval <- wt$p.value
  }
  med_pos <- if (n_pos > 0) stats::median(df$pathway_score[df$group == pos], na.rm = TRUE) else NA_real_
  med_neg <- if (n_neg > 0) stats::median(df$pathway_score[df$group == neg], na.rm = TRUE) else NA_real_
  data.table(
    grouping = grouping, pos_level = pos, neg_level = neg,
    n_pos = n_pos, n_neg = n_neg,
    median_pos = med_pos, median_neg = med_neg,
    delta_median = med_pos - med_neg,
    pvalue = pval
  )
}

as_symbol_matrix <- function(fpkm_data, probe_annot) {
  fpkm_data <- as.data.table(fpkm_data)
  id_col <- names(fpkm_data)[1]
  rid <- as.character(fpkm_data[[id_col]])
  ens <- sub("\\..*$", "", rid)
  mat <- as.matrix(fpkm_data[, -1, with = FALSE])
  storage.mode(mat) <- "double"
  colnames(mat) <- normalize_barcode(colnames(mat))
  if (anyDuplicated(colnames(mat))) {
    mat <- mat[, !duplicated(colnames(mat)), drop = FALSE]
  }
  probe_annot <- as.data.table(probe_annot)
  pid <- first_present(names(probe_annot), c("id", "ensembl", "Ensembl", "ENSEMBL", names(probe_annot)[1]))
  psym <- first_present(names(probe_annot), c("gene", "symbol", "Gene", "SYMBOL", names(probe_annot)[2]))
  map_ens <- sub("\\..*$", "", as.character(probe_annot[[pid]]))
  map_sym <- as.character(probe_annot[[psym]])
  if (mean(grepl("^ENSG", ens, ignore.case = TRUE), na.rm = TRUE) >= 0.5) {
    sym <- map_sym[match(ens, map_ens)]
    keep <- !is.na(sym) & sym != "" & !grepl("^XLOC_", sym)
    mat <- mat[keep, , drop = FALSE]
    sym <- sym[keep]
  } else {
    sym <- rid
  }
  if (anyDuplicated(sym)) {
    dt <- data.table(symbol = sym, as.data.table(mat))
    dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = symbol]
    mat <- as.matrix(dt[, -1, with = FALSE])
    rownames(mat) <- dt$symbol
  } else {
    rownames(mat) <- sym
  }
  mat
}

save_plot <- function(p, path_stub, width = 11, height = 8) {
  tryCatch(
    ggsave(paste0(path_stub, ".pdf"), p, width = width, height = height),
    error = function(e) message("保存 PDF 失败：", conditionMessage(e))
  )
  tryCatch(
    ggsave(paste0(path_stub, ".png"), p, width = width, height = height, dpi = 150),
    error = function(e) message("保存 PNG 失败：", conditionMessage(e))
  )
  if (interactive()) {
    tryCatch(print(p), error = function(e) {
      message("预览图失败（已跳过）：", conditionMessage(e))
    })
  }
}

head_dt <- function(dt, n) {
  if (is.null(dt) || nrow(dt) == 0L || n <= 0L) return(dt[0])
  dt[seq_len(min(as.integer(n), nrow(dt)))]
}

pick_volcano_labels <- function(plot_dt, n_neg = 12L, n_pos = 6L) {
  neg <- head_dt(plot_dt[significant_neg == TRUE], n_neg)
  pos <- head_dt(plot_dt[spearman_r > 0][order(-spearman_r)], n_pos)
  out <- unique(rbindlist(list(neg, pos), fill = TRUE))
  if (nrow(out) == 0L) return(out)
  out[is.finite(spearman_r) & is.finite(neglogp) & !is.na(feature) & nzchar(as.character(feature))]
}

# 不用 ggrepel：ggplot2 4.x 下会报 depth(NULL)
add_volcano_labels <- function(p, lab_dt) {
  if (is.null(lab_dt) || nrow(lab_dt) == 0L) return(p)
  p + ggplot2::geom_text(
    data = as.data.frame(lab_dt),
    aes(x = spearman_r, y = neglogp, label = feature),
    size = 2.3, vjust = -0.55, inherit.aes = FALSE, check_overlap = TRUE
  )
}

# 每个 GO 拆成两列：未转移 / 转移，分别画该组通路分数中位数
expand_met_two_cols <- function(stat_dt) {
  d <- copy(as.data.table(stat_dt))
  if (nrow(d) == 0) return(d)
  rbindlist(list(
    d[, .(
      GO, GO_name, grouping, panel,
      side = "Non-metastatic", x_lab = neg_lab,
      n = n_neg, median_score = median_neg, pvalue
    )],
    d[, .(
      GO, GO_name, grouping, panel,
      side = "Metastatic", x_lab = pos_lab,
      n = n_pos, median_score = median_pos, pvalue
    )]
  ), fill = TRUE)
}

plot_go_bubble_two_cols <- function(stat_dt, title, subtitle, path_stub, facet = FALSE) {
  long <- expand_met_two_cols(stat_dt)
  long <- long[is.finite(median_score)]
  if (nrow(long) == 0) return(invisible(NULL))
  go_lv <- unique(paste0(stat_dt$GO, "  ", stat_dt$GO_name))
  long[, y_lab := factor(paste0(GO, "  ", GO_name), levels = rev(go_lv))]
  long[, x_lab := factor(x_lab, levels = unique(c(stat_dt$neg_lab, stat_dt$pos_lab)))]
  long[, neglogp := ifelse(is.finite(pvalue), pmin(10, -log10(pmax(pvalue, 1e-12))), 0.5)]
  fill_lim <- max(abs(long$median_score), na.rm = TRUE)
  if (!is.finite(fill_lim) || fill_lim < 1e-8) fill_lim <- 0.1
  p <- ggplot(long, aes(x = x_lab, y = y_lab)) +
    geom_point(
      aes(size = neglogp, fill = median_score),
      shape = 21, color = "black", stroke = 0.5 / ggplot2::.pt
    ) +
    scale_fill_gradientn(
      colours = c("#3C5488", "#5B7FA6", "#FFFFFF", "#EE8A7A", "#E64B35"),
      values = c(0, 0.47, 0.50, 0.53, 1),
      limits = c(-fill_lim, fill_lim),
      oob = scales::squish,
      name = "Pathway score\nmedian"
    ) +
    scale_size_continuous(range = c(3, 11), name = expression(-log[10](p))) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1, size = 11),
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      strip.text = element_text(size = 10)
    )
  if (isTRUE(facet) && "panel" %in% names(long) && uniqueN(long$panel) > 1) {
    p <- p + facet_wrap(~ panel, nrow = 1, scales = "free_x")
  }
  n_panel <- if (isTRUE(facet)) max(1, uniqueN(long$panel)) else 1
  save_plot(p, path_stub, max(8, 4.2 * n_panel + 4), max(6, 0.38 * uniqueN(long$y_lab) + 2.4))
}

# ==============================================================================
# 读数据（分析函数在 expr 建好之后才调用）
# ==============================================================================
if (dir.exists(work_dir)) {
  setwd(work_dir)
} else {
  message("未找到 ", work_dir, " ，改用当前目录：", getwd())
}

fpkm_file <- pick_best_table(c("TCGA-BRCA.star_fpkm.tsv"), must = TRUE, min_cols = 20L, label = "FPKM")
clin_file <- pick_best_clinical()
probe_file <- pick_best_table(c(
  "gencode.v36.annotation.gtf.gene.probemap",
  "gencode.v36.annotation.gtf.gene.probemap.tsv"
), must = TRUE, min_cols = 2L, label = "基因注释")
surv_file <- pick_best_table(c("TCGA-BRCA.survival.tsv"), must = FALSE, min_cols = 2L, label = "生存表")

message("FPKM：", fpkm_file)
message("临床：", clin_file)
message("注释：", probe_file)
if (!is.na(surv_file)) message("生存：", surv_file)

fpkm_data <- safe_fread(fpkm_file)
probe_annot <- safe_fread(probe_file)
clinical_data <- safe_fread(clin_file)
survival_data <- if (!is.na(surv_file)) safe_fread(surv_file) else NULL
protein_data <- NULL

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("正在把 FPKM 映射为基因符号矩阵…")
expr_all <- as_symbol_matrix(fpkm_data, probe_annot)
mx <- suppressWarnings(max(expr_all, na.rm = TRUE))
if (is.finite(mx) && mx > 50) {
  message("检测到原始 FPKM（max=", round(mx, 2), "），做 log2(x+1)")
  expr_all <- log2(expr_all + 1)
} else {
  message("表达值范围较小（max=", round(mx, 2), "），视为已 log 转换")
}

build_annotation <- function(sample_ids, clin, surv) {
  ann <- data.table(sample = normalize_barcode(sample_ids))
  ann <- ann[!is.na(sample) & sample != ""]
  ann <- ann[!duplicated(sample)]
  ann[, patient := patient_id(sample)]
  ann[, barcode_type := sample_type_code(sample)]

  clin <- as.data.table(clin)
  idc <- first_present(names(clin), c(
    "sample", "sampleID", "sample_id", "submitter_id.samples",
    "bcr_sample_barcode", "bcr_patient_barcode", "submitter_id", names(clin)[1]
  ))
  clin[, sample_std := normalize_barcode(clin[[idc]])]
  clin[, patient_std := patient_id(clin[[idc]])]

  st <- pick_clin_col(clin, c(
    "ajcc_pathologic_stage", "ajcc_pathologic_tumor_stage",
    "pathologic_stage", "tumor_stage", "clinical_stage"
  ))
  mc <- pick_clin_col(clin, c(
    "ajcc_pathologic_m", "ajcc_metastasis_pathologic_pm", "pathologic_m", "clinical_m"
  ))
  nc <- pick_clin_col(clin, c(
    "ajcc_pathologic_n", "ajcc_nodes_pathologic_pn", "pathologic_n", "clinical_n"
  ))
  tc <- pick_clin_col(clin, c("sample_type\\.samples", "sample_type", "tumor_descriptor"))
  message("选用临床列：stage=", st, "  M=", mc, "  N=", nc, "  sample_type=", tc)

  extra <- clin[, unique(c("sample_std", "patient_std",
                           na.omit(c(st, mc, nc, tc)))), with = FALSE]
  extra_s <- extra[!duplicated(sample_std)]
  ann <- merge(ann, extra_s, by.x = "sample", by.y = "sample_std", all.x = TRUE)
  extra_p <- extra[!duplicated(patient_std)]
  for (nm in setdiff(names(extra_p), c("sample_std", "patient_std", names(ann)))) {
    src <- extra_p[[nm]][match(ann$patient, extra_p$patient_std)]
    cur <- if (nm %in% names(ann)) ann[[nm]] else rep(NA, nrow(ann))
    miss <- is.na(cur) | as.character(cur) %in% c("", "NA")
    cur[miss] <- src[miss]
    ann[, (nm) := cur]
  }

  n <- nrow(ann)
  st_vec <- if (!is.na(st) && st %in% names(ann)) simplify_stage(ann[[st]]) else rep(NA_character_, n)
  m_vec  <- if (!is.na(mc) && mc %in% names(ann)) classify_m(ann[[mc]]) else rep(NA_character_, n)
  n_vec  <- if (!is.na(nc) && nc %in% names(ann)) classify_n(ann[[nc]]) else rep(NA_character_, n)
  raw_type <- if (!is.na(tc) && tc %in% names(ann)) ann[[tc]] else rep(NA_character_, n)
  type_vec <- classify_sample_type(raw_type, ann$sample)
  stage_iv <- rep(NA_character_, n)
  stage_iv[!is.na(st_vec) & st_vec == "Stage IV"] <- "Stage IV"
  stage_iv[!is.na(st_vec) & st_vec %in% c("Stage I", "Stage II", "Stage III")] <- "Stage I-III"

  ann[, stage_simplified := st_vec]
  ann[, distant_M := factor(m_vec, levels = c("M0", "M1"))]
  ann[, node_N := factor(n_vec, levels = c("N0", "Nplus"))]
  ann[, stage_IV := factor(stage_iv, levels = c("Stage I-III", "Stage IV"))]
  ann[, sample_class := factor(type_vec, levels = c("原位肿瘤", "转移组织", "正常"))]
  ann
}

ann_all <- build_annotation(colnames(expr_all), clinical_data, survival_data)
ann_all <- ann_all[sample %in% colnames(expr_all)]
fwrite(ann_all, file.path(out_dir, "00_sample_annotation.csv"))
message(
  "样本：原位=", sum(ann_all$sample_class == "原位肿瘤", na.rm = TRUE),
  "  转移组织=", sum(ann_all$sample_class == "转移组织", na.rm = TRUE),
  "  正常=", sum(ann_all$sample_class == "正常", na.rm = TRUE),
  "  M1=", sum(ann_all$distant_M == "M1", na.rm = TRUE),
  "  Stage IV=", sum(ann_all$stage_simplified == "Stage IV", na.rm = TRUE),
  "  N+=", sum(ann_all$node_N == "Nplus", na.rm = TRUE)
)

primary_ids <- ann_all$sample[ann_all$sample_class == "原位肿瘤"]
met_ids <- ann_all$sample[ann_all$sample_class == "转移组织"]
tumor_ids <- unique(c(primary_ids, met_ids))
if (length(primary_ids) < 20) stop("原位肿瘤样本太少：", length(primary_ids))

expr_tumor <- expr_all[, intersect(tumor_ids, colnames(expr_all)), drop = FALSE]
keep_g <- rowMeans(is.finite(expr_tumor) & expr_tumor > 0, na.rm = TRUE) >= min_expr_frac
expr_tumor <- expr_tumor[keep_g, , drop = FALSE]
expr_primary <- expr_tumor[, intersect(primary_ids, colnames(expr_tumor)), drop = FALSE]
message("表达矩阵：原位 ", ncol(expr_primary), " 样本 x ", nrow(expr_primary), " 基因；肿瘤(含转移组织) ", ncol(expr_tumor))

# ==============================================================================
# 主分析（必须在 expr 建好之后）
# ==============================================================================
run_tcga_brca_2026 <- function() {
  if (!exists("expr_primary", inherits = TRUE)) stop("还没有 expr_primary，请从脚本开头 Source")

  score_one <- function(expr_mat, go_id) {
    genes <- unique(get_go_genes(go_id)$SYMBOL)
    sc <- pathway_zmean(expr_mat, genes)
    if (is.null(sc)) {
      message("  ", go_id, " 映射基因不足，跳过")
      return(NULL)
    }
    message("  ", go_id, "  ", go_title(go_id), "  基因数=", attr(sc, "n_genes"))
    sc
  }

  message("计算原位肿瘤各神经 GO 通路分数（每个 GO 单独，不合并）")
  score_primary_list <- lapply(go_list, function(g) score_one(expr_primary, g))
  names(score_primary_list) <- go_list
  score_primary_list <- Filter(Negate(is.null), score_primary_list)
  if (length(score_primary_list) == 0) stop("没有任何神经 GO 能打分")

  message("计算原位+转移组织各神经 GO 通路分数（供 1d）")
  score_tumor_list <- lapply(names(score_primary_list), function(g) score_one(expr_tumor, g))
  names(score_tumor_list) <- names(score_primary_list)
  score_tumor_list <- Filter(Negate(is.null), score_tumor_list)

  mat_from_list <- function(lst) {
    common <- Reduce(intersect, lapply(lst, names))
    mat <- do.call(cbind, lapply(lst, function(x) x[common]))
    colnames(mat) <- names(lst)
    rownames(mat) <- common
    mat
  }
  score_primary_mat <- mat_from_list(score_primary_list)
  score_tumor_mat <- mat_from_list(score_tumor_list)
  fwrite(data.table(sample = rownames(score_primary_mat), as.data.table(score_primary_mat)),
         file.path(out_dir, "01_pathway_scores_primary_each_GO.csv"))
  fwrite(data.table(sample = rownames(score_tumor_mat), as.data.table(score_tumor_mat)),
         file.path(out_dir, "01_pathway_scores_tumor_each_GO.csv"))

  gene_map_rows <- lapply(names(score_primary_list), function(g) {
    data.table(GO = g, GO_name = go_title(g), n_genes = attr(score_primary_list[[g]], "n_genes"),
               genes = paste(attr(score_primary_list[[g]], "genes"), collapse = ";"))
  })
  fwrite(rbindlist(gene_map_rows), file.path(out_dir, "00_GO_genes_used.csv"))

  ann_p <- ann_all[sample %in% rownames(score_primary_mat)]
  ann_t <- ann_all[sample %in% rownames(score_tumor_mat)]

  designs <- list(
    list(key = "a_distant_M", title = "1a Distant metastasis", panel = "1a Distant M",
         group = setNames(as.character(ann_p$distant_M), ann_p$sample),
         pos = "M1", neg = "M0",
         pos_lab = "M1", neg_lab = "M0",
         score_mat = score_primary_mat),
    list(key = "b_AJCC_stageIV", title = "1b AJCC stage", panel = "1b AJCC stage",
         group = setNames(as.character(ann_p$stage_IV), ann_p$sample),
         pos = "Stage IV", neg = "Stage I-III",
         pos_lab = "Stage IV", neg_lab = "Stage I-III",
         score_mat = score_primary_mat),
    list(key = "c_node_N", title = "1c Lymph node", panel = "1c Lymph node",
         group = setNames(as.character(ann_p$node_N), ann_p$sample),
         pos = "Nplus", neg = "N0",
         pos_lab = "N+", neg_lab = "N0",
         score_mat = score_primary_mat),
    list(key = "d_sample_type", title = "1d Sample type", panel = "1d Sample type",
         group = setNames(as.character(ann_t$sample_class), ann_t$sample),
         pos = "转移组织", neg = "原位肿瘤",
         pos_lab = "Metastatic tissue", neg_lab = "Primary tumor",
         score_mat = score_tumor_mat)
  )

  all_go_stats <- list()
  for (ds in designs) {
    message("气泡图：", ds$title)
    stat_rows <- list()
    sm <- ds$score_mat
    for (g in colnames(sm)) {
      sc <- as.numeric(sm[, g])
      names(sc) <- rownames(sm)
      one <- compare_go_groups(sc, ds$group[names(sc)], ds$pos, ds$neg, ds$key)
      if (is.null(one)) next
      one[, `:=`(
        GO = g, GO_name = go_title(g),
        panel = ds$panel, pos_lab = ds$pos_lab, neg_lab = ds$neg_lab
      )]
      stat_rows[[g]] <- one
    }
    stat_dt <- rbindlist(stat_rows, fill = TRUE)
    if (nrow(stat_dt) == 0) {
      message("  分组人数不足，跳过 ", ds$key)
      next
    }
    stat_dt[, fdr := p.adjust(pvalue, method = "BH")]
    fwrite(stat_dt, file.path(out_dir, paste0("02_", ds$key, "_GO_vs_metastasis.csv")))
    fwrite(expand_met_two_cols(stat_dt),
           file.path(out_dir, paste0("02_", ds$key, "_GO_vs_metastasis_two_cols.csv")))
    all_go_stats[[ds$key]] <- stat_dt
    plot_go_bubble_two_cols(
      stat_dt,
      title = paste0(ds$title, ": neural GO in non-metastatic vs metastatic"),
      subtitle = "Y = neuronal GO (scored separately); X = non-metastatic | metastatic; fill = median pathway score; size = -log10(Wilcoxon p)",
      path_stub = file.path(out_dir, paste0("02_", ds$key, "_bubble"))
    )
  }

  if (length(all_go_stats) > 0) {
    bubble <- rbindlist(all_go_stats, fill = TRUE)
    fwrite(bubble, file.path(out_dir, "02_summary_GO_vs_metastasis.csv"))
    fwrite(expand_met_two_cols(bubble),
           file.path(out_dir, "02_summary_GO_vs_metastasis_two_cols.csv"))
    plot_go_bubble_two_cols(
      bubble,
      title = "Neural GO activity versus breast cancer metastasis",
      subtitle = "Each panel: non-metastatic | metastatic; fill = median pathway score; size = -log10(p)",
      path_stub = file.path(out_dir, "02_summary_bubble_GO_vs_metastasis"),
      facet = TRUE
    )
    if (interactive()) {
      message("主气泡图已保存：02_summary_bubble_GO_vs_metastasis.png")
    }
  }

  # ---- 2) 原位肿瘤内，与转移负相关的基因 ----
  met_defs <- list(
    list(key = "a_distant_M", title = "2a Genes negatively correlated with distant M1 (primary tumors)",
         group = setNames(as.character(ann_p$distant_M), ann_p$sample),
         pos = "M1", neg = "M0"),
    list(key = "b_AJCC_stageIV", title = "2b Genes negatively correlated with Stage IV (primary tumors)",
         group = setNames(as.character(ann_p$stage_IV), ann_p$sample),
         pos = "Stage IV", neg = "Stage I-III"),
    list(key = "c_node_N", title = "2c Genes negatively correlated with N+ (primary tumors)",
         group = setNames(as.character(ann_p$node_N), ann_p$sample),
         pos = "Nplus", neg = "N0")
  )
  for (md in met_defs) {
    message("全基因组相关：", md$title)
    tab <- spearman_vs_binary(expr_primary, md$group, md$pos, md$neg)
    if (nrow(tab) == 0) next
    tab[, `:=`(
      significant_neg = spearman_r < neg_r_cutoff & pvalue < neg_pvalue_cutoff,
      strict_neg = spearman_r <= strict_r_cutoff & pvalue < neg_pvalue_cutoff
    )]
    setorder(tab, spearman_r)
    fwrite(tab, file.path(out_dir, paste0("03_", md$key, "_genes_vs_metastasis_all.csv")))
    fwrite(tab[significant_neg == TRUE], file.path(out_dir, paste0("03_", md$key, "_genes_NEG_vs_metastasis.csv")))
    fwrite(tab[strict_neg == TRUE], file.path(out_dir, paste0("03_", md$key, "_genes_NEG_strict_vs_metastasis.csv")))
    message("  负相关基因 ", sum(tab$significant_neg),
            "（严格 r<=", strict_r_cutoff, "：", sum(tab$strict_neg), "）")

    plot_dt <- copy(tab)
    plot_dt[, neglogp := pmin(12, -log10(pmax(pvalue, 1e-12)))]
    plot_dt[, col := ifelse(significant_neg, "Negative",
                            ifelse(spearman_r > 0 & pvalue < neg_pvalue_cutoff, "Positive", "NS"))]
    top_lab <- pick_volcano_labels(plot_dt, 12L, 6L)
    p_vol <- tryCatch({
      p <- ggplot(plot_dt, aes(x = spearman_r, y = neglogp, color = col)) +
        geom_point(alpha = 0.45, size = 0.7) +
        geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
        geom_hline(yintercept = -log10(neg_pvalue_cutoff), linetype = 2, color = "grey50") +
        scale_color_manual(values = c("Negative" = "#3C5488", "Positive" = "#E64B35", "NS" = "grey75")) +
        labs(title = md$title,
             subtitle = paste0("Spearman: gene vs ", md$pos, " (1) / ", md$neg, " (0); blue = negative vs metastasis"),
             x = "Spearman r", y = expression(-log[10](p)), color = NULL) +
        theme_bw()
      add_volcano_labels(p, top_lab)
    }, error = function(e) {
      message("火山图失败（已跳过）：", conditionMessage(e))
      NULL
    })
    if (!is.null(p_vol)) {
      save_plot(p_vol, file.path(out_dir, paste0("03_", md$key, "_volcano_genes_vs_metastasis")), 10, 7)
    }
  }

  if (length(met_ids) >= min_group_n) {
    common <- intersect(colnames(expr_tumor), c(primary_ids, met_ids))
    grp <- ifelse(common %in% met_ids, "转移组织", "原位肿瘤")
    names(grp) <- common
    tab <- spearman_vs_binary(expr_tumor[, common, drop = FALSE], grp, "转移组织", "原位肿瘤")
    if (nrow(tab) > 0) {
      tab[, `:=`(
        significant_neg = spearman_r < neg_r_cutoff & pvalue < neg_pvalue_cutoff,
        strict_neg = spearman_r <= strict_r_cutoff & pvalue < neg_pvalue_cutoff
      )]
      setorder(tab, spearman_r)
      fwrite(tab[significant_neg == TRUE],
             file.path(out_dir, "03_d_sample_type_genes_NEG_vs_metastatic_tissue.csv"))
      fwrite(tab, file.path(out_dir, "03_d_sample_type_genes_vs_metastatic_tissue_all.csv"))
    }
  }

  # ---- 3) 原位肿瘤内，与每个神经 GO 负相关的基因 ----
  dir.create(file.path(out_dir, "04_neg_vs_neural_per_GO"), showWarnings = FALSE)
  summary_neg <- list()
  for (g in names(score_primary_list)) {
    message("与神经浸润负相关：", g, " ", go_title(g))
    go_score_vec <- score_primary_list[[g]]
    tab <- spearman_vs_go_score(expr_primary, go_score_vec)
    if (nrow(tab) == 0) next
    tab[, `:=`(
      GO = g, GO_name = go_title(g),
      significant_neg = spearman_r < neg_r_cutoff & pvalue < neg_pvalue_cutoff,
      strict_neg = spearman_r <= strict_r_cutoff & pvalue < neg_pvalue_cutoff
    )]
    setorder(tab, spearman_r)
    gdir <- file.path(out_dir, "04_neg_vs_neural_per_GO", paste0("GO_", safe_name(sub("GO:", "", g))))
    dir.create(gdir, showWarnings = FALSE)
    fwrite(tab, file.path(gdir, "genes_vs_neural_GO_all.csv"))
    fwrite(tab[significant_neg == TRUE], file.path(gdir, "genes_NEG_vs_neural_GO.csv"))
    fwrite(tab[strict_neg == TRUE], file.path(gdir, "genes_NEG_strict_vs_neural_GO.csv"))
    summary_neg[[g]] <- data.table(
      GO = g, GO_name = go_title(g),
      n_pathway_genes = attr(go_score_vec, "n_genes"),
      n_tested = nrow(tab),
      n_neg = sum(tab$significant_neg),
      n_neg_strict = sum(tab$strict_neg)
    )

    plot_dt <- copy(tab)
    plot_dt[, neglogp := pmin(12, -log10(pmax(pvalue, 1e-12)))]
    plot_dt[, col := ifelse(significant_neg, "Negative",
                            ifelse(spearman_r > 0 & pvalue < neg_pvalue_cutoff, "Positive", "NS"))]
    top_lab <- pick_volcano_labels(plot_dt, 10L, 5L)
    p_vol <- tryCatch({
      p <- ggplot(plot_dt, aes(x = spearman_r, y = neglogp, color = col)) +
        geom_point(alpha = 0.4, size = 0.65) +
        geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
        geom_hline(yintercept = -log10(neg_pvalue_cutoff), linetype = 2, color = "grey50") +
        scale_color_manual(values = c("Negative" = "#3C5488", "Positive" = "#E64B35", "NS" = "grey75")) +
        labs(title = "Genes negatively correlated with neural invasion",
             subtitle = paste0(g, "  ", go_title(g), "; blue = negative vs this GO score"),
             x = "Spearman r (gene vs neural GO score)",
             y = expression(-log[10](p)), color = NULL) +
        theme_bw()
      add_volcano_labels(p, top_lab)
    }, error = function(e) {
      message("火山图失败（已跳过）：", conditionMessage(e))
      NULL
    })
    if (!is.null(p_vol)) {
      save_plot(p_vol, file.path(gdir, "volcano_neg_vs_neural_GO"), 9, 6.5)
    }
  }
  if (length(summary_neg) > 0) {
    sum_dt <- rbindlist(summary_neg, fill = TRUE)
    sum_dt[, y_lab := paste0(GO, "  ", GO_name)]
    fwrite(sum_dt, file.path(out_dir, "04_summary_neg_genes_vs_each_neural_GO.csv"))
    p_n <- ggplot(sum_dt, aes(x = n_neg, y = reorder(y_lab, n_neg))) +
      geom_col(fill = "#3C5488", width = 0.7) +
      labs(title = "Number of genes negatively correlated with each neural GO",
           subtitle = paste0("Spearman r < 0 and p < ", neg_pvalue_cutoff, "; GO sets not pooled"),
           x = "Number of negative genes", y = NULL) +
      theme_bw()
    save_plot(p_n, file.path(out_dir, "04_summary_neg_gene_counts"), 10, 6)
  }

  message("完成。结果目录：", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
  message("主气泡图：02_summary_bubble_GO_vs_metastasis.png")
  message("分面气泡图：02_a/b/c/d_*_bubble.png")
  message("转移负相关基因：03_*_genes_NEG_vs_metastasis.csv")
  message("神经 GO 负相关基因：04_neg_vs_neural_per_GO/")
  invisible(TRUE)
}

if (exists("expr_primary") && ncol(expr_primary) > 10) {
  run_tcga_brca_2026()
} else {
  stop("表达矩阵未建好，请从第一行完整 Source")
}
