################################################################################
# TCGA-BRCA-2026.R
# 单独脚本：放到 E:/R/TCGA-BRCA-2026，RStudio 打开后从第一行 Source
# 不要从中间粘贴，不要运行旧的 factor(meta_M) 行
#
# 数据目录（.tsv / .tsv.gz 均可；空的 .tsv 会自动跳过，改读 .gz）：
#   TCGA-BRCA.star_fpkm.tsv(.gz)         # 若 .tsv 是 WinRAR 留下的空文件，删掉即可
#   TCGA-BRCA.clinical.tsv(.gz)
#   TCGA-BRCA.survival.tsv(.gz)          # 可选
#   gencode.v36.annotation.gtf.gene.probemap
#
# 神经浸润标定（不用神经 GO，两套签名分开打分，不把基因合并成一套）：
#   A. 施旺细胞 marker
#   B. 神经营养因子（配体 + 常用受体）
# 原位肿瘤：Primary Tumor / 条形码 01
#
# 1) 神经浸润 vs 乳腺癌转移（施旺 / 神经营养 各做一遍）
#    a. 诊断时远处转移：原位肿瘤 M1 vs M0
#    b. AJCC 分期：原位肿瘤 Stage IV vs I–III，以及 I/II/III/IV
#    c. 淋巴结：原位肿瘤 N+ vs N0
#    d. 样本类型：转移组织 vs 原位肿瘤（转移组织很少，约 7 例）
# 2) 原位肿瘤内，与转移呈负相关的基因（转移标准与 1a/1b/1c 相同）
# 3) 原位肿瘤内，与神经浸润分数呈负相关的基因（施旺、神经营养分开）
################################################################################

library(data.table)
library(ggplot2)
library(ggpubr)

# ==============================================================================
# 参数
# ==============================================================================
work_dir <- Sys.getenv("TCGA_BRCA_2026_DIR", unset = "E:/R/TCGA-BRCA-2026")
out_dir <- "results_TCGA-BRCA-2026"
min_signature_genes <- 1
min_group_n <- 2
neg_pvalue_cutoff <- 0.05
neg_r_cutoff <- 0
strict_r_cutoff <- -0.15
min_expr_frac <- 0.20

# 施旺细胞 marker（髓鞘 / 施旺谱系 / 神经束膜常用基因，HGNC 符号）
schwann_markers <- c(
  "S100B", "SOX10", "MPZ", "MBP", "PMP22", "PLP1", "MAG", "PRX",
  "GFAP", "NGFR", "NCAM1", "L1CAM", "CDH19", "ERBB3", "NRG1",
  "GAP43", "SCN7A", "EGR2", "POU3F1", "MAL", "GJB1", "NES"
)

# 神经营养因子：配体与常用受体（不与施旺 marker 合并打分）
neurotrophin_markers <- c(
  "NGF", "BDNF", "NTF3", "NTF4", "GDNF", "NRTN", "ARTN", "PSPN",
  "CNTF", "CDNF", "MANF", "VGF",
  "NTRK1", "NTRK2", "NTRK3", "RET",
  "GFRA1", "GFRA2", "GFRA3", "GFRA4", "CNTFR"
)

# 图上单独展示的关键基因（必须属于上面两套之一）
highlight_genes <- c(
  "S100B", "SOX10", "MPZ", "MBP", "NGFR", "ERBB3",
  "NGF", "BDNF", "NTF3", "GDNF", "ARTN", "NTRK1", "NTRK2", "RET"
)

signature_list <- list(
  Schwann = schwann_markers,
  Neurotrophin = neurotrophin_markers
)
signature_title <- c(
  Schwann = "施旺细胞 marker",
  Neurotrophin = "神经营养因子"
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

existing_files <- function(stems) {
  cands <- unique(unlist(lapply(stems, function(s) c(s, paste0(s, ".gz")))))
  cands[file.exists(cands) & !is.na(file.info(cands)$isdir) & !file.info(cands)$isdir]
}

# WinRAR 没解压完时会留下空的 .tsv；空文件 / 0 列一律跳过，改用 .gz
probe_table_file <- function(f) {
  info <- file.info(f)
  if (is.null(info) || is.na(info$size) || info$size < 200) {
    return(data.table(file = f, size = if (is.null(info) || is.na(info$size)) 0 else info$size,
                      n_col = 0L, readable = FALSE, err = "empty_or_tiny"))
  }
  hdr <- tryCatch(names(fread(f, nrows = 0, fill = TRUE, showProgress = FALSE)),
                  error = function(e) character())
  data.table(
    file = f, size = as.numeric(info$size), n_col = length(hdr),
    readable = length(hdr) >= 1L,
    err = if (length(hdr) >= 1L) NA_character_ else "unreadable_header"
  )
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
      stop(
        label, "存在但读不了（常见原因：WinRAR 还没解压完，留下了空的 .tsv）。",
        "请删掉空的 .tsv，直接保留 .tsv.gz，或等解压完成后再 Source。\n",
        "已看到：", paste(tab$file, collapse = " ; ")
      )
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
  stop(
    "读不了 ", path, "：", conditionMessage(out),
    "。若这是空的 .tsv，请删除它并改用同名 .tsv.gz（data.table 可以直接读 gz）。"
  )
}

pick_one_file <- function(stems, must = TRUE) {
  pick_best_table(stems, must = must, min_cols = 1L, label = "文件")
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
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  xl <- toupper(x)
  out[grepl("PRIMARY", xl)] <- "原位肿瘤"
  out[grepl("METASTATIC", xl)] <- "转移组织"
  out[grepl("NORMAL", xl)] <- "正常"
  code <- sample_type_code(barcode)
  out[is.na(out) & code == "01"] <- "原位肿瘤"
  out[is.na(out) & code %in% c("06", "07")] <- "转移组织"
  out[is.na(out) & code == "11"] <- "正常"
  out
}

# 签名活性：基因 z-score 后对样本取均值。不要用 scale()/t()，变量不要叫 score
signature_zmean <- function(expr_mat, genes) {
  genes <- unique(intersect(as.character(genes), rownames(expr_mat)))
  if (length(genes) < min_signature_genes) return(NULL)
  sub <- as.matrix(expr_mat[genes, , drop = FALSE])
  storage.mode(sub) <- "double"
  gene_mean <- rowMeans(sub, na.rm = TRUE)
  gene_sd <- sqrt(rowMeans((sub - gene_mean)^2, na.rm = TRUE))
  gene_sd[!is.finite(gene_sd) | gene_sd < 1e-12] <- 1
  z <- (sub - gene_mean) / gene_sd
  z[!is.finite(z)] <- 0
  sig_score <- colMeans(z, na.rm = TRUE)
  names(sig_score) <- colnames(sub)
  attr(sig_score, "n_genes") <- length(genes)
  attr(sig_score, "genes") <- genes
  sig_score
}

spearman_vs_sig <- function(mat, sig_vec) {
  common <- intersect(colnames(mat), names(sig_vec))
  if (length(common) < 5) return(data.table())
  mat <- mat[, common, drop = FALSE]
  sig_vec <- sig_vec[common]
  keep <- apply(mat, 1, function(x) stats::sd(x, na.rm = TRUE) > 0)
  mat <- mat[keep, , drop = FALSE]
  n <- ncol(mat)
  r <- as.numeric(cor(
    base::t(as.matrix(mat)),
    sig_vec,
    method = "spearman",
    use = "pairwise.complete.obs"
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

compare_sig_groups <- function(sig_vec, group, pos, neg, grouping) {
  df <- data.frame(
    signature_score = as.numeric(sig_vec),
    group = as.character(group),
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$signature_score) & df$group %in% c(pos, neg), ]
  n_pos <- sum(df$group == pos)
  n_neg <- sum(df$group == neg)
  if (n_pos < min_group_n || n_neg < min_group_n) return(NULL)
  wt <- suppressWarnings(stats::wilcox.test(signature_score ~ group, data = df))
  data.table(
    grouping = grouping, pos_level = pos, neg_level = neg,
    n_pos = n_pos, n_neg = n_neg,
    median_pos = stats::median(df$signature_score[df$group == pos], na.rm = TRUE),
    median_neg = stats::median(df$signature_score[df$group == neg], na.rm = TRUE),
    delta_median = stats::median(df$signature_score[df$group == pos], na.rm = TRUE) -
      stats::median(df$signature_score[df$group == neg], na.rm = TRUE),
    pvalue = wt$p.value
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
  ggsave(paste0(path_stub, ".pdf"), p, width = width, height = height)
  ggsave(paste0(path_stub, ".png"), p, width = width, height = height, dpi = 150)
  if (interactive()) print(p)
}
theme_pub <- function() {
  theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 9),
      plot.title = element_text(size = 13, face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1)
    )
}

maybe_label <- function() {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    return(ggrepel::geom_text_repel)
  }
  message("未安装 ggrepel，火山图不标基因名。可运行 install.packages(\"ggrepel\")")
  function(...) ggplot2::geom_blank()
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

# ==============================================================================
# 临床分组（只用向量赋值）
# ==============================================================================
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
  label_fun <- maybe_label()

  used_rows <- list()
  score_one <- function(expr_mat, sig_id, genes) {
    present <- intersect(genes, rownames(expr_mat))
    missing <- setdiff(genes, rownames(expr_mat))
    used_rows[[length(used_rows) + 1]] <<- data.table(
      signature = sig_id, signature_name = unname(signature_title[sig_id]),
      gene = genes,
      in_matrix = genes %in% rownames(expr_mat)
    )
    sc <- signature_zmean(expr_mat, present)
    if (is.null(sc)) {
      message("  ", sig_id, " 在表达矩阵中找不到 marker，跳过")
      return(NULL)
    }
    message("  ", sig_id, "  ", signature_title[sig_id],
            "  用到 ", attr(sc, "n_genes"), " 个基因：",
            paste(attr(sc, "genes"), collapse = ", "))
    if (length(missing) > 0) message("    缺失：", paste(missing, collapse = ", "))
    sc
  }

  message("计算原位肿瘤神经浸润分数（施旺 / 神经营养 分开）")
  score_primary_list <- list()
  for (sid in names(signature_list)) {
    score_primary_list[[sid]] <- score_one(expr_primary, sid, signature_list[[sid]])
  }
  score_primary_list <- Filter(Negate(is.null), score_primary_list)
  if (length(score_primary_list) == 0) stop("施旺细胞 marker 和神经营养因子都无法打分")

  message("计算原位+转移组织神经浸润分数（供 1d）")
  score_tumor_list <- list()
  for (sid in names(score_primary_list)) {
    score_tumor_list[[sid]] <- score_one(expr_tumor, sid, signature_list[[sid]])
  }
  score_tumor_list <- Filter(Negate(is.null), score_tumor_list)

  if (length(used_rows) > 0) {
    used_dt <- unique(rbindlist(used_rows, fill = TRUE))
    fwrite(used_dt, file.path(out_dir, "00_neural_markers_used.csv"))
  }

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
         file.path(out_dir, "01_neural_invasion_scores_primary.csv"))
  fwrite(data.table(sample = rownames(score_tumor_mat), as.data.table(score_tumor_mat)),
         file.path(out_dir, "01_neural_invasion_scores_tumor.csv"))

  ann_p <- ann_all[sample %in% rownames(score_primary_mat)]
  ann_t <- ann_all[sample %in% rownames(score_tumor_mat)]

  # ---- 1) 神经浸润 vs 转移：四个思路 × 两套签名 ----
  designs <- list(
    list(key = "a_distant_M", title = "1a 诊断时远处转移（原位肿瘤）",
         subtitle = "神经浸润分数：原位肿瘤 M1 vs M0",
         group = setNames(as.character(ann_p$distant_M), ann_p$sample),
         pos = "M1", neg = "M0", score_mat = score_primary_mat,
         fill = c("M0" = "#4DBBD5", "M1" = "#E64B35")),
    list(key = "b_AJCC_stageIV", title = "1b AJCC 分期（原位肿瘤）",
         subtitle = "神经浸润分数：Stage IV vs Stage I–III",
         group = setNames(as.character(ann_p$stage_IV), ann_p$sample),
         pos = "Stage IV", neg = "Stage I-III", score_mat = score_primary_mat,
         fill = c("Stage I-III" = "#4DBBD5", "Stage IV" = "#E64B35")),
    list(key = "c_node_N", title = "1c 淋巴结（原位肿瘤）",
         subtitle = "神经浸润分数：N+ vs N0",
         group = setNames(as.character(ann_p$node_N), ann_p$sample),
         pos = "Nplus", neg = "N0", score_mat = score_primary_mat,
         fill = c("N0" = "#4DBBD5", "Nplus" = "#E64B35")),
    list(key = "d_sample_type", title = "1d 样本类型（原位 vs 转移组织）",
         subtitle = "神经浸润分数：转移组织 vs 原位肿瘤（转移组织例数很少）",
         group = setNames(as.character(ann_t$sample_class), ann_t$sample),
         pos = "转移组织", neg = "原位肿瘤", score_mat = score_tumor_mat,
         fill = c("原位肿瘤" = "#4DBBD5", "转移组织" = "#E64B35"))
  )

  all_sig_stats <- list()
  for (ds in designs) {
    message("作图：", ds$title)
    stat_rows <- list()
    long_rows <- list()
    sm <- ds$score_mat
    for (sid in colnames(sm)) {
      sc <- as.numeric(sm[, sid])
      names(sc) <- rownames(sm)
      grp <- ds$group[names(sc)]
      one <- compare_sig_groups(sc, grp, ds$pos, ds$neg, ds$key)
      if (!is.null(one)) {
        one[, `:=`(signature = sid, signature_name = unname(signature_title[sid]))]
        stat_rows[[sid]] <- one
      }
      long_rows[[sid]] <- data.table(
        sample = names(sc),
        signature = sid,
        signature_name = unname(signature_title[sid]),
        signature_score = as.numeric(sc),
        group = factor(as.character(grp), levels = c(ds$neg, ds$pos))
      )
    }
    stat_dt <- rbindlist(stat_rows, fill = TRUE)
    if (nrow(stat_dt) == 0) {
      message("  分组人数不足，跳过 ", ds$key)
      next
    }
    stat_dt[, fdr := p.adjust(pvalue, method = "BH")]
    fwrite(stat_dt, file.path(out_dir, paste0("02_", ds$key, "_neural_vs_group.csv")))
    all_sig_stats[[ds$key]] <- stat_dt

    long_dt <- rbindlist(long_rows, fill = TRUE)
    long_dt <- long_dt[!is.na(group) & is.finite(signature_score)]
    p_box <- ggplot(long_dt, aes(x = group, y = signature_score, fill = group)) +
      geom_boxplot(outlier.size = 0.4, width = 0.65) +
      stat_compare_means(size = 3.2, label = "p.format") +
      facet_wrap(~ signature_name, scales = "free_y", ncol = 2) +
      scale_fill_manual(values = ds$fill) +
      labs(title = ds$title, subtitle = ds$subtitle,
           x = NULL, y = "Neural invasion signature score", fill = NULL) +
      theme_pub()
    save_plot(p_box, file.path(out_dir, paste0("02_", ds$key, "_boxplot")), 10, 6)

    p_delta <- ggplot(stat_dt, aes(x = delta_median, y = reorder(signature_name, delta_median))) +
      geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
      geom_point(aes(color = pvalue < 0.05, size = -log10(pmax(pvalue, 1e-12)))) +
      scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "#E64B35"), name = "p < 0.05") +
      labs(title = paste0(ds$title, "：神经浸润分数差"),
           subtitle = paste0("Δ median (", ds$pos, " − ", ds$neg, ")"),
           x = "Δ median signature score", y = NULL, size = expression(-log[10](p))) +
      theme_bw()
    save_plot(p_delta, file.path(out_dir, paste0("02_", ds$key, "_delta")), 8, 4)

    # 关键单个 marker
    expr_use <- if (identical(ds$key, "d_sample_type")) expr_tumor else expr_primary
    genes_ok <- intersect(highlight_genes, rownames(expr_use))
    if (length(genes_ok) >= 1) {
      g_long <- rbindlist(lapply(genes_ok, function(gn) {
        v <- as.numeric(expr_use[gn, ])
        names(v) <- colnames(expr_use)
        data.table(
          gene = gn, sample = names(v), expression = v,
          group = factor(as.character(ds$group[names(v)]), levels = c(ds$neg, ds$pos))
        )
      }), fill = TRUE)
      g_long <- g_long[!is.na(group) & is.finite(expression)]
      p_g <- ggplot(g_long, aes(x = group, y = expression, fill = group)) +
        geom_boxplot(outlier.size = 0.3, width = 0.65) +
        stat_compare_means(size = 2.4, label = "p.format") +
        facet_wrap(~ gene, scales = "free_y", ncol = 4) +
        scale_fill_manual(values = ds$fill) +
        labs(title = paste0(ds$title, "：关键 marker 单基因表达"),
             subtitle = "施旺细胞 marker 与神经营养因子；未合并打分",
             x = NULL, y = "log expression", fill = NULL) +
        theme_pub()
      save_plot(p_g, file.path(out_dir, paste0("02_", ds$key, "_key_markers")), 12, 8)
    }
  }

  # 1b 附加：原位肿瘤按 Stage I–IV
  if (sum(!is.na(ann_p$stage_simplified)) >= 10) {
    long_st <- rbindlist(lapply(colnames(score_primary_mat), function(sid) {
      sc <- score_primary_mat[, sid]
      data.table(
        signature = sid, signature_name = unname(signature_title[sid]),
        sample = names(sc), signature_score = as.numeric(sc),
        stage = factor(as.character(ann_p$stage_simplified[match(names(sc), ann_p$sample)]),
                       levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
      )
    }), fill = TRUE)
    long_st <- long_st[!is.na(stage) & is.finite(signature_score)]
    p_st <- ggplot(long_st, aes(x = stage, y = signature_score, fill = stage)) +
      geom_boxplot(outlier.size = 0.35, width = 0.7) +
      facet_wrap(~ signature_name, scales = "free_y", ncol = 2) +
      scale_fill_manual(values = c("Stage I" = "#3C5488", "Stage II" = "#4DBBD5",
                                   "Stage III" = "#E64B35", "Stage IV" = "#F39B7F")) +
      labs(title = "1b AJCC 分期（原位肿瘤，I–IV）",
           subtitle = "施旺细胞 marker 与神经营养因子分开打分",
           x = NULL, y = "Neural invasion signature score", fill = NULL) +
      theme_pub()
    save_plot(p_st, file.path(out_dir, "02_b_AJCC_stage_I_to_IV_boxplot"), 10, 6)
  }

  if (length(all_sig_stats) > 0) {
    bubble <- rbindlist(all_sig_stats, fill = TRUE)
    bubble[, grp_lab := factor(grouping, levels = c(
      "a_distant_M", "b_AJCC_stageIV", "c_node_N", "d_sample_type"
    ), labels = c("1a 远处转移 M1 vs M0", "1b Stage IV vs I–III",
                  "1c N+ vs N0", "1d 转移组织 vs 原位"))]
    p_bub <- ggplot(bubble, aes(x = grp_lab, y = signature_name)) +
      geom_point(aes(size = -log10(pmax(pvalue, 1e-12)),
                     color = pmin(pmax(delta_median, -1), 1))) +
      scale_color_gradient2(low = "#3C5488", mid = "white", high = "#E64B35",
                            midpoint = 0, name = "Δ median") +
      scale_size_continuous(name = expression(-log[10](p))) +
      labs(title = "神经浸润 vs 乳腺癌转移（四个定义）",
           subtitle = "红=转移侧神经浸润更高，蓝=更低；施旺与神经营养未合并",
           x = NULL, y = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "right")
    save_plot(p_bub, file.path(out_dir, "02_summary_bubble_neural_vs_metastasis"), 10, 5)
    fwrite(bubble, file.path(out_dir, "02_summary_neural_vs_metastasis.csv"))
  }

  # ---- 2) 原位肿瘤内，与转移负相关的基因 ----
  met_defs <- list(
    list(key = "a_distant_M", title = "2a 与诊断时远处转移负相关的基因（原位肿瘤）",
         group = setNames(as.character(ann_p$distant_M), ann_p$sample),
         pos = "M1", neg = "M0"),
    list(key = "b_AJCC_stageIV", title = "2b 与 Stage IV 负相关的基因（原位肿瘤）",
         group = setNames(as.character(ann_p$stage_IV), ann_p$sample),
         pos = "Stage IV", neg = "Stage I-III"),
    list(key = "c_node_N", title = "2c 与淋巴结 N+ 负相关的基因（原位肿瘤）",
         group = setNames(as.character(ann_p$node_N), ann_p$sample),
         pos = "Nplus", neg = "N0")
  )

  for (md in met_defs) {
    message("全基因组相关：", md$title)
    tab <- spearman_vs_binary(expr_primary, md$group, md$pos, md$neg)
    if (nrow(tab) == 0) {
      message("  分组不足，跳过")
      next
    }
    tab[, `:=`(
      direction = ifelse(spearman_r < 0, "negative", "positive"),
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
    plot_dt[, col := ifelse(significant_neg, "负相关",
                            ifelse(spearman_r > 0 & pvalue < neg_pvalue_cutoff, "正相关", "不显著"))]
    top_lab <- rbind(
      plot_dt[significant_neg == TRUE][1:min(12L, .N)],
      plot_dt[spearman_r > 0][order(-spearman_r)][1:min(6L, .N)]
    )
    p_vol <- ggplot(plot_dt, aes(x = spearman_r, y = neglogp, color = col)) +
      geom_point(alpha = 0.45, size = 0.7) +
      geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
      geom_hline(yintercept = -log10(neg_pvalue_cutoff), linetype = 2, color = "grey50") +
      scale_color_manual(values = c("负相关" = "#3C5488", "正相关" = "#E64B35", "不显著" = "grey75")) +
      label_fun(data = top_lab, aes(label = feature), size = 2.4, max.overlaps = 30, show.legend = FALSE) +
      labs(title = md$title,
           subtitle = paste0("Spearman：基因 vs ", md$pos, "(1)/", md$neg, "(0)；蓝=与转移负相关"),
           x = "Spearman r", y = expression(-log[10](p)), color = NULL) +
      theme_bw()
    save_plot(p_vol, file.path(out_dir, paste0("03_", md$key, "_volcano_genes_vs_metastasis")), 10, 7)
  }

  if (length(met_ids) >= min_group_n) {
    message("全基因组：原位 vs 转移组织（对应 1d）")
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

  # ---- 3) 原位肿瘤内，与神经浸润负相关的基因（施旺 / 神经营养分开）----
  dir.create(file.path(out_dir, "04_neg_vs_neural_per_signature"), showWarnings = FALSE)
  summary_neg <- list()
  for (sid in names(score_primary_list)) {
    message("与神经浸润负相关：", sid, " ", signature_title[sid])
    sig_vec <- score_primary_list[[sid]]
    tab <- spearman_vs_sig(expr_primary, sig_vec)
    if (nrow(tab) == 0) next
    tab[, `:=`(
      signature = sid, signature_name = unname(signature_title[sid]),
      significant_neg = spearman_r < neg_r_cutoff & pvalue < neg_pvalue_cutoff,
      strict_neg = spearman_r <= strict_r_cutoff & pvalue < neg_pvalue_cutoff
    )]
    setorder(tab, spearman_r)
    sdir <- file.path(out_dir, "04_neg_vs_neural_per_signature", safe_name(sid))
    dir.create(sdir, showWarnings = FALSE)
    fwrite(tab, file.path(sdir, "genes_vs_neural_score_all.csv"))
    fwrite(tab[significant_neg == TRUE], file.path(sdir, "genes_NEG_vs_neural_score.csv"))
    fwrite(tab[strict_neg == TRUE], file.path(sdir, "genes_NEG_strict_vs_neural_score.csv"))
    summary_neg[[sid]] <- data.table(
      signature = sid, signature_name = unname(signature_title[sid]),
      n_signature_genes = attr(sig_vec, "n_genes"),
      n_tested = nrow(tab),
      n_neg = sum(tab$significant_neg),
      n_neg_strict = sum(tab$strict_neg)
    )

    plot_dt <- copy(tab)
    plot_dt[, neglogp := pmin(12, -log10(pmax(pvalue, 1e-12)))]
    plot_dt[, col := ifelse(significant_neg, "负相关",
                            ifelse(spearman_r > 0 & pvalue < neg_pvalue_cutoff, "正相关", "不显著"))]
    top_lab <- rbind(
      plot_dt[significant_neg == TRUE][1:min(10L, .N)],
      plot_dt[spearman_r > 0][order(-spearman_r)][1:min(5L, .N)]
    )
    p_vol <- ggplot(plot_dt, aes(x = spearman_r, y = neglogp, color = col)) +
      geom_point(alpha = 0.4, size = 0.65) +
      geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
      geom_hline(yintercept = -log10(neg_pvalue_cutoff), linetype = 2, color = "grey50") +
      scale_color_manual(values = c("负相关" = "#3C5488", "正相关" = "#E64B35", "不显著" = "grey75")) +
      label_fun(data = top_lab, aes(label = feature), size = 2.3, max.overlaps = 25, show.legend = FALSE) +
      labs(title = "3 原位肿瘤：与神经浸润负相关的基因",
           subtitle = paste0(signature_title[sid], "；蓝=与该签名分数负相关"),
           x = "Spearman r (gene vs neural invasion score)",
           y = expression(-log[10](p)), color = NULL) +
      theme_bw()
    save_plot(p_vol, file.path(sdir, "volcano_neg_vs_neural"), 9, 6.5)
  }
  if (length(summary_neg) > 0) {
    sum_dt <- rbindlist(summary_neg, fill = TRUE)
    fwrite(sum_dt, file.path(out_dir, "04_summary_neg_genes_vs_neural_signatures.csv"))
    p_n <- ggplot(sum_dt, aes(x = n_neg, y = reorder(signature_name, n_neg))) +
      geom_col(fill = "#3C5488", width = 0.65) +
      labs(title = "原位肿瘤中与神经浸润负相关的基因数",
           subtitle = paste0("Spearman r < 0 且 p < ", neg_pvalue_cutoff,
                             "；施旺与神经营养分开，未合并"),
           x = "负相关基因数", y = NULL) +
      theme_bw()
    save_plot(p_n, file.path(out_dir, "04_summary_neg_gene_counts"), 8, 4)
  }

  message("完成。结果目录：", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
  message("神经浸润分数：01_neural_invasion_scores_primary.csv")
  message("主图：02_summary_bubble_neural_vs_metastasis.png")
  message("转移负相关基因：03_*_genes_NEG_vs_metastasis.csv")
  message("神经浸润负相关基因：04_neg_vs_neural_per_signature/")
  invisible(TRUE)
}

if (exists("expr_primary") && ncol(expr_primary) > 10) {
  run_tcga_brca_2026()
} else {
  stop("表达矩阵未建好，请从第一行完整 Source")
}
