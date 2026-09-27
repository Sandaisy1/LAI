################################################################################
# AURORA_nerve.R
# 单独脚本，不要和 nerve_TCGA.R / TCGA-BRCA-2026.R 在同一会话 Source
# 数据目录：E:/R/cBioportal breast cancer/AURORA
# 也自动识别子目录 brca_aurora_2023（cBioPortal AURORA US, Nat Cancer 2023）
# RStudio 打开后从第一行 Source
#
# 需要的文件（cBioPortal 下载即可；.txt / .txt.gz）：
#   data_clinical_sample.txt
#   data_clinical_patient.txt
#   data_mrna_seq_v2_rsem.txt                  # 优先，原始 RSEM
#   或 data_mrna_seq_v2_rsem_zscores_*.txt     # 仅有 z-score 时也可用
#
# 神经浸润只按下面 5 个神经信号 GO，每个单独取基因、单独打分（不合并）：
#   GO:0019227  neuronal action potential propagation
#   GO:1902847  regulation of neuronal signal transduction
#   GO:0097374  sensory neuron axon guidance
#   GO:1902667  regulation of axon guidance
#   GO:0007409  axonogenesis
#
# 打分：主 z-mean；补充 z-median、ssGSEA
# 转移定义（人数不足的组会跳过并写日志）：
#   a 原发样本 诊断 M1 vs M0（PRIM_M / YPM）
#   b 原发样本 Stage IV vs I-III
#   c 原发样本 N+ vs N0（YPN / PRIM_N）
#   d 转移组织 vs 原发组织（SAMPLE_TYPE，AURORA 主比较）
#   e 同一患者配对：转移 vs 原发（有配对 RNA 才做）
#
# 结果目录：results_AURORA_nerve/
################################################################################

library(data.table)
library(ggplot2)

aurora_out_dir <- "results_AURORA_nerve"
min_set_genes <- 1
min_group_n <- 2
min_expr_frac <- 0.20

go_list <- c(
  "GO:0019227",
  "GO:1902847",
  "GO:0097374",
  "GO:1902667",
  "GO:0007409"
)

go_name_map <- c(
  "GO:0019227" = "neuronal action potential propagation",
  "GO:1902847" = "regulation of neuronal signal transduction",
  "GO:0097374" = "sensory neuron axon guidance",
  "GO:1902667" = "regulation of axon guidance",
  "GO:0007409" = "axonogenesis"
)

go_fallback <- list(
  "GO:0019227" = c(
    "CNTNAP1", "SCN1B", "SCN1A", "SCN2A", "SCN8A", "SCN2B", "SCN4B",
    "ANK3", "NFASC", "NRCAM", "CNTN2", "SPTBN4", "KCNA1", "KCNQ2", "KCNQ3"
  ),
  "GO:1902847" = c(
    "CLU", "ADCY1", "CAMK2A", "GRIN1", "GRIN2B", "DLG4",
    "HOMER1", "SYNGAP1", "RGS4", "GNAO1", "PRKCA", "MAPK1"
  ),
  "GO:0097374" = c(
    "NRP1", "NRP2", "SEMA3A", "SEMA3F", "PLXNA1", "PLXNA3", "PLXNA4",
    "NTN1", "DCC", "SLIT1", "ROBO1"
  ),
  "GO:1902667" = c(
    "ATOH7", "KIF21A", "MYCBP2", "NOVA2", "POU4F2", "PTPRO",
    "ROBO3", "SLIT2", "TUBB2B", "YTHDF1"
  ),
  "GO:0007409" = strsplit(paste(
    "ABL1,ADGRB1,ADNP,ALCAM,AMIGO1,ANK3,ANOS1,APBB1,APLP1,APOE,APP,ATL1,ATOH7,AUTS2,",
    "BAIAP2,BDNF,BMPR2,BRSK1,BRSK2,CDH2,CDK5,CDK5R1,CELSR1,CELSR2,CELSR3,CHN1,",
    "CNTN1,CNTN2,CNTNAP1,CTNNA2,CXCL12,DCC,DCLK1,DISC1,DOCK7,DPYSL5,DRAXIN,DSCAM,",
    "EFNA1,EFNA2,EFNA3,EFNA4,EFNA5,EFNB1,EFNB2,EFNB3,ENAH,EPHA3,EPHA4,EPHA5,EPHA6,",
    "EPHA7,EPHA8,EPHB1,EPHB2,EPHB3,FEZ1,FEZ2,FGF13,FYN,FZD3,GAP43,GDNF,GSK3B,",
    "ISL1,ISL2,KALRN,KIF21A,KIF5C,KLF7,L1CAM,LHX1,LIMK1,LRRC4C,MACF1,MAP1B,MAP2,",
    "MAPT,MYCBP2,NCAM1,NEFH,NEO1,NFASC,NRCAM,NRP1,NRP2,NRXN1,NTN1,NTN4,NTNG1,",
    "NTRK1,NTRK2,PAFAH1B1,PAK1,PLXNA3,PLXNA4,PLXNB1,POU4F1,POU4F2,PTEN,PTK2,",
    "PTPRO,RELN,RET,ROBO1,ROBO2,ROBO3,RTN4,RTN4R,SEMA3A,SEMA3C,SEMA3E,SEMA3F,",
    "SEMA4D,SEMA5A,SEMA6A,SEMA6D,SHH,SLIT1,SLIT2,SLIT3,SLITRK1,SPAST,SPTBN4,",
    "TENM1,TENM2,TENM3,TENM4,TUBB3,ULK1,UNC5A,UNC5B,UNC5C,UNC5D,WNT5A,WNT7A"
  ), ",", fixed = TRUE)[[1]]
)
go_fallback <- lapply(go_fallback, function(x) unique(trimws(x)))

go_title <- function(go_id) {
  go_id <- as.character(go_id)
  out <- unname(go_name_map[go_id])
  miss <- is.na(out) | !nzchar(out)
  out[miss] <- go_id[miss]
  out
}
go_lab <- function(go_id) paste0(as.character(go_id), "  ", go_title(go_id))

# ==============================================================================
# 工具函数（全部先定义，读完表达矩阵再分析）
# ==============================================================================
first_present <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) NA_character_ else hit[1]
}
norm_id <- function(x) {
  x <- toupper(gsub("[. ]", "-", as.character(x)))
  gsub("-+", "-", x)
}
looks_like_cbioportal <- function(dir) {
  if (!dir.exists(dir)) return(FALSE)
  hits <- list.files(dir, pattern = "data_clinical_sample|data_mrna", ignore.case = TRUE)
  length(hits) > 0
}
resolve_aurora_dir <- function() {
  env <- Sys.getenv("AURORA_NERVE_DIR", unset = "")
  cands <- unique(c(
    env,
    "E:/R/cBioportal breast cancer/AURORA",
    "E:/R/cBioportal breast cancer/AURORA/brca_aurora_2023",
    file.path(getwd(), "brca_aurora_2023"),
    getwd()
  ))
  cands <- cands[nzchar(cands)]
  extra <- unlist(lapply(cands, function(d) {
    if (!dir.exists(d)) return(character())
    file.path(d, list.files(d, pattern = "aurora", ignore.case = TRUE, include.dirs = TRUE))
  }))
  cands <- unique(c(cands, extra))
  hit <- cands[vapply(cands, looks_like_cbioportal, logical(1))]
  if (length(hit) == 0) {
    stop(
      "找不到 AURORA / cBioPortal 数据。请把脚本放到数据目录，或设置 AURORA_NERVE_DIR。\n",
      "已试：", paste(cands, collapse = " | ")
    )
  }
  normalizePath(hit[1], winslash = "/", mustWork = FALSE)
}

fread_cbioportal <- function(path, ...) {
  if (!file.exists(path)) stop("找不到文件：", path)
  npeek <- 30L
  peek <- tryCatch(readLines(path, n = npeek, warn = FALSE), error = function(e) character())
  n_hash <- sum(grepl("^#", peek))
  out <- tryCatch(
    fread(path, sep = "\t", header = TRUE, skip = n_hash, fill = TRUE, showProgress = TRUE, ...),
    error = function(e) e
  )
  if (!inherits(out, "error")) return(out)
  alt <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    sub("\\.gz$", "", path, ignore.case = TRUE)
  } else {
    paste0(path, ".gz")
  }
  if (file.exists(alt)) {
    message("读 ", path, " 失败，改读 ", alt)
    peek2 <- tryCatch(readLines(alt, n = npeek, warn = FALSE), error = function(e) character())
    n_hash2 <- sum(grepl("^#", peek2))
    return(fread(alt, sep = "\t", header = TRUE, skip = n_hash2, fill = TRUE, showProgress = TRUE, ...))
  }
  stop("读不了 ", path, "：", conditionMessage(out))
}

list_data_files <- function(dir, pattern) {
  hits <- list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  hits <- hits[!is.na(file.info(hits)$isdir) & !file.info(hits)$isdir]
  hits
}

pick_expression_file <- function(dir) {
  all_m <- list_data_files(dir, "^data_mrna.*\\.(txt|tsv)(\\.gz)?$")
  if (length(all_m) == 0) {
    stop("找不到表达矩阵：", dir, " 下应有 data_mrna_seq_*.txt")
  }
  rawish <- all_m[!grepl("zscore", basename(all_m), ignore.case = TRUE)]
  pick <- if (length(rawish) > 0) rawish else all_m
  info <- file.info(pick)
  pick <- pick[order(-info$size)]
  message("表达候选：\n", paste(sprintf("  %s  %.1fMB", pick, file.info(pick)$size / 1024^2), collapse = "\n"))
  pick[1]
}

pick_clin_file <- function(dir, stem) {
  hits <- list_data_files(dir, paste0("^", stem, "\\.(txt|tsv)(\\.gz)?$"))
  if (length(hits) == 0) return(NA_character_)
  hits[which.max(file.info(hits)$size)]
}

as_symbol_matrix_cbioportal <- function(expr_dt) {
  expr_dt <- as.data.table(expr_dt)
  nms <- names(expr_dt)
  gene_col <- first_present(nms, c(
    "Hugo_Symbol", "HUGO_SYMBOL", "hugo_symbol", "Gene_Symbol",
    "gene_symbol", "SYMBOL", "Gene", "gene"
  ))
  drop <- intersect(nms, c(
    gene_col, "Entrez_Gene_Id", "ENTREZ_GENE_ID", "entrez_gene_id",
    "Entrez", "ENSEMBL", "Ensembl", "ensembl"
  ))
  if (is.na(gene_col)) {
    gene_col <- nms[1]
    drop <- nms[1]
  }
  samp_cols <- setdiff(nms, drop)
  if (length(samp_cols) < 5) stop("表达表样本列太少：", length(samp_cols))
  sym <- as.character(expr_dt[[gene_col]])
  mat <- as.matrix(expr_dt[, samp_cols, with = FALSE])
  storage.mode(mat) <- "double"
  colnames(mat) <- norm_id(colnames(mat))
  keep <- !is.na(sym) & nzchar(sym) & !grepl("^NA$|^-$|^\\.$", sym)
  mat <- mat[keep, , drop = FALSE]
  sym <- sym[keep]
  if (anyDuplicated(sym)) {
    dt <- data.table(symbol = sym, as.data.table(mat))
    dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = symbol]
    mat <- as.matrix(dt[, -1, with = FALSE])
    rownames(mat) <- dt$symbol
  } else {
    rownames(mat) <- sym
  }
  if (anyDuplicated(colnames(mat))) {
    mat <- mat[, !duplicated(colnames(mat)), drop = FALSE]
  }
  mat
}

detect_already_z <- function(mat, path) {
  if (grepl("zscore", basename(path), ignore.case = TRUE)) return(TRUE)
  finite <- mat[is.finite(mat)]
  if (length(finite) < 100) return(FALSE)
  isTRUE(max(abs(finite), na.rm = TRUE) < 40 && abs(stats::median(finite, na.rm = TRUE)) < 1)
}

classify_sample_class <- function(sample_type, sample_id) {
  out <- rep(NA_character_, length(sample_type))
  xl <- toupper(as.character(sample_type))
  sid <- toupper(as.character(sample_id))
  out[grepl("METASTA", xl)] <- "Metastatic"
  out[grepl("PRIMARY|PRIMARY SOLID", xl)] <- "Primary"
  out[is.na(out) & grepl("-TTM|\\bTTM", sid)] <- "Metastatic"
  out[is.na(out) & grepl("-TTP|\\bTTP", sid)] <- "Primary"
  out
}
classify_m <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("M1", x)] <- "M1"
  out[is.na(out) & grepl("M0", x)] <- "M0"
  out[grepl("MX|UNKNOWN|NOT AVAILABLE|NOT REPORTED", x)] <- NA_character_
  out
}
classify_n <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("N[1-3]", x)] <- "Nplus"
  out[is.na(out) & grepl("N0", x)] <- "N0"
  out[grepl("NX|UNKNOWN|NOT AVAILABLE|NOT REPORTED", x)] <- NA_character_
  out
}
classify_stage <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("IV|STAGE.?4", x)] <- "Stage IV"
  out[is.na(out) & grepl("III|STAGE.?3|II|STAGE.?2|I\\b|STAGE.?1", x)] <- "Stage I-III"
  out[grepl("X\\b|UNKNOWN|NOT AVAILABLE|NOT REPORTED", x)] <- NA_character_
  out
}

# 通路活性。already_z=TRUE 时不再做 z（cBioPortal 已是全样本 z-score）。变量不要叫 score
pathway_zstat <- function(expr_mat, genes, how = "mean", already_z = FALSE) {
  genes <- unique(intersect(as.character(genes), rownames(expr_mat)))
  if (length(genes) < min_set_genes) return(NULL)
  sub <- as.matrix(expr_mat[genes, , drop = FALSE])
  storage.mode(sub) <- "double"
  if (isTRUE(already_z)) {
    z <- sub
    z[!is.finite(z)] <- 0
  } else {
    gene_mean <- rowMeans(sub, na.rm = TRUE)
    gene_sd <- sqrt(rowMeans((sub - gene_mean)^2, na.rm = TRUE))
    gene_sd[!is.finite(gene_sd) | gene_sd < 1e-12] <- 1
    z <- (sub - gene_mean) / gene_sd
    z[!is.finite(z)] <- 0
  }
  if (identical(how, "median")) {
    set_score <- apply(z, 2, stats::median, na.rm = TRUE)
  } else {
    set_score <- colMeans(z, na.rm = TRUE)
  }
  names(set_score) <- colnames(sub)
  attr(set_score, "n_genes") <- length(genes)
  attr(set_score, "genes") <- genes
  set_score
}
pathway_zmean <- function(expr_mat, genes, already_z = FALSE) {
  pathway_zstat(expr_mat, genes, "mean", already_z)
}
pathway_zmedian <- function(expr_mat, genes, already_z = FALSE) {
  pathway_zstat(expr_mat, genes, "median", already_z)
}

ssgsea_via_gsva <- function(expr_mat, gene_sets) {
  if (!requireNamespace("GSVA", quietly = TRUE)) return(NULL)
  gsets <- lapply(gene_sets, function(g) unique(intersect(as.character(g), rownames(expr_mat))))
  gsets <- Filter(function(g) length(g) >= min_set_genes, gsets)
  if (length(gsets) == 0) return(NULL)
  mat <- as.matrix(expr_mat)
  storage.mode(mat) <- "double"
  scored <- tryCatch({
    if (exists("ssgseaParam", envir = asNamespace("GSVA"), inherits = FALSE)) {
      param <- GSVA::ssgseaParam(mat, gsets, normalize = TRUE)
      GSVA::gsva(param, verbose = FALSE)
    } else {
      GSVA::gsva(mat, gsets, method = "ssgsea", ssgsea.norm = TRUE, verbose = FALSE)
    }
  }, error = function(e) {
    message("GSVA ssGSEA 失败，改用脚本内实现：", conditionMessage(e))
    NULL
  })
  if (is.null(scored)) return(NULL)
  scored <- as.matrix(scored)
  out <- lapply(colnames(scored), function(nm) {
    v <- as.numeric(scored[, nm])
    names(v) <- rownames(scored)
    attr(v, "n_genes") <- length(gsets[[nm]])
    attr(v, "genes") <- gsets[[nm]]
    v
  })
  names(out) <- colnames(scored)
  out
}

ssgsea_builtin <- function(expr_mat, gene_sets, tau = 0.25) {
  mat <- as.matrix(expr_mat)
  storage.mode(mat) <- "double"
  genes_all <- rownames(mat)
  ng <- nrow(mat)
  ns <- ncol(mat)
  set_idx <- lapply(gene_sets, function(g) unique(which(genes_all %in% unique(as.character(g)))))
  keep <- vapply(set_idx, length, integer(1)) >= min_set_genes
  set_idx <- set_idx[keep]
  if (length(set_idx) == 0) return(list())
  out_mat <- matrix(NA_real_, ns, length(set_idx),
                    dimnames = list(colnames(mat), names(set_idx)))
  pos_tau <- seq_len(ng)^tau
  for (j in seq_len(ns)) {
    o <- order(mat[, j], decreasing = TRUE, na.last = TRUE)
    for (k in seq_along(set_idx)) {
      hit <- o %in% set_idx[[k]]
      n_hit <- sum(hit)
      n_miss <- ng - n_hit
      if (n_hit < min_set_genes || n_miss < 1) next
      hit_w <- numeric(ng)
      hit_w[hit] <- pos_tau[hit]
      s <- sum(hit_w)
      if (!is.finite(s) || s <= 0) next
      walk <- cumsum(hit_w / s - ifelse(hit, 0, 1 / n_miss))
      out_mat[j, k] <- sum(walk)
    }
  }
  for (k in seq_len(ncol(out_mat))) {
    mx <- max(abs(out_mat[, k]), na.rm = TRUE)
    if (is.finite(mx) && mx > 0) out_mat[, k] <- out_mat[, k] / mx
  }
  out <- lapply(colnames(out_mat), function(nm) {
    v <- as.numeric(out_mat[, nm])
    names(v) <- rownames(out_mat)
    attr(v, "n_genes") <- length(set_idx[[nm]])
    attr(v, "genes") <- genes_all[set_idx[[nm]]]
    v
  })
  names(out) <- colnames(out_mat)
  out
}

score_ssgsea_sets <- function(expr_mat, gene_sets) {
  via <- ssgsea_via_gsva(expr_mat, gene_sets)
  if (!is.null(via) && length(via) > 0) {
    message("  ssGSEA 后端：GSVA")
    return(via)
  }
  message("  ssGSEA 后端：脚本内置（未安装 GSVA）")
  ssgsea_builtin(expr_mat, gene_sets)
}

get_go_genes <- function(go_id) {
  mapped <- NULL
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    pick <- function(keytype) {
      tryCatch(
        AnnotationDbi::select(
          org.Hs.eg.db::org.Hs.eg.db, keys = go_id, keytype = keytype,
          columns = c("SYMBOL", "ENSEMBL", "ENTREZID")
        ),
        error = function(e) NULL
      )
    }
    mapped <- suppressMessages(pick("GOALL"))
    if (is.null(mapped) || nrow(mapped) == 0) mapped <- suppressMessages(pick("GO"))
    extra <- tryCatch({
      go2eg <- as.list(org.Hs.eg.db::org.Hs.egGO2ALLEGS)
      eg <- unique(as.character(go2eg[[go_id]]))
      if (length(eg) == 0) NULL else {
        AnnotationDbi::select(
          org.Hs.eg.db::org.Hs.eg.db, keys = eg, keytype = "ENTREZID",
          columns = c("SYMBOL", "ENSEMBL")
        )
      }
    }, error = function(e) NULL)
    if (!is.null(extra) && nrow(extra) > 0) {
      extra <- as.data.table(extra)
      extra[, GO := go_id]
      mapped <- if (is.null(mapped) || nrow(mapped) == 0) extra else {
        rbind(as.data.table(mapped), extra, fill = TRUE)
      }
    }
  }
  if (!is.null(mapped) && nrow(mapped) > 0) {
    mapped <- as.data.table(mapped)
    if ("GOALL" %in% names(mapped)) setnames(mapped, "GOALL", "GO", skip_absent = TRUE)
    if (!"GO" %in% names(mapped)) mapped[, GO := go_id]
    syms <- unique(mapped[!is.na(SYMBOL) & SYMBOL != "", SYMBOL])
    if (length(syms) > 0) return(list(genes = syms, source = "org.Hs.eg.db"))
  }
  fb <- go_fallback[[go_id]]
  if (is.null(fb) || length(fb) == 0) return(list(genes = character(), source = "none"))
  list(genes = unique(fb), source = "fallback")
}

compare_groups <- function(value_vec, group, pos, neg, grouping) {
  df <- data.frame(
    nerve_value = as.numeric(value_vec),
    group = as.character(group),
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$nerve_value) & df$group %in% c(pos, neg), ]
  n_pos <- sum(df$group == pos)
  n_neg <- sum(df$group == neg)
  if (n_pos < 1 && n_neg < 1) return(NULL)
  pval <- NA_real_
  if (n_pos >= min_group_n && n_neg >= min_group_n) {
    wt <- suppressWarnings(stats::wilcox.test(nerve_value ~ group, data = df))
    pval <- wt$p.value
  }
  med_pos <- if (n_pos > 0) stats::median(df$nerve_value[df$group == pos], na.rm = TRUE) else NA_real_
  med_neg <- if (n_neg > 0) stats::median(df$nerve_value[df$group == neg], na.rm = TRUE) else NA_real_
  data.table(
    grouping = grouping, pos_level = pos, neg_level = neg,
    n_pos = n_pos, n_neg = n_neg,
    median_pos = med_pos, median_neg = med_neg,
    delta_median = med_pos - med_neg,
    pvalue = pval, test = "wilcoxon_unpaired"
  )
}

compare_paired_patients <- function(value_vec, ann, grouping) {
  dt <- data.table(
    sample = names(value_vec),
    nerve_value = as.numeric(value_vec)
  )
  dt <- merge(dt, ann[, .(sample, patient, sample_class)], by = "sample", all.x = TRUE)
  dt <- dt[is.finite(nerve_value) & sample_class %in% c("Primary", "Metastatic")]
  if (nrow(dt) == 0) return(NULL)
  pat <- dt[, .(
    primary_mean = mean(nerve_value[sample_class == "Primary"], na.rm = TRUE),
    met_mean = mean(nerve_value[sample_class == "Metastatic"], na.rm = TRUE)
  ), by = patient]
  pat <- pat[is.finite(primary_mean) & is.finite(met_mean)]
  n <- nrow(pat)
  if (n < min_group_n) return(NULL)
  pval <- suppressWarnings(stats::wilcox.test(pat$met_mean, pat$primary_mean, paired = TRUE)$p.value)
  data.table(
    grouping = grouping, pos_level = "Metastatic", neg_level = "Primary",
    n_pos = n, n_neg = n,
    median_pos = stats::median(pat$met_mean, na.rm = TRUE),
    median_neg = stats::median(pat$primary_mean, na.rm = TRUE),
    delta_median = stats::median(pat$met_mean - pat$primary_mean, na.rm = TRUE),
    pvalue = pval, test = "wilcoxon_paired"
  )
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

expand_two_cols <- function(stat_dt) {
  d <- copy(as.data.table(stat_dt))
  if (nrow(d) == 0) return(d)
  rbindlist(list(
    d[, .(
      GO, GO_name, grouping, panel,
      side = "Non-metastatic", x_lab = neg_lab,
      n = n_neg, median_value = median_neg, pvalue
    )],
    d[, .(
      GO, GO_name, grouping, panel,
      side = "Metastatic", x_lab = pos_lab,
      n = n_pos, median_value = median_pos, pvalue
    )]
  ), fill = TRUE)
}

plot_bubble_two_cols <- function(stat_dt, title, subtitle, path_stub, facet = FALSE) {
  long <- expand_two_cols(stat_dt)
  long <- long[is.finite(median_value)]
  if (nrow(long) == 0) return(invisible(NULL))
  go_lv <- unique(go_lab(stat_dt$GO))
  long[, y_lab := factor(go_lab(GO), levels = rev(go_lv))]
  long[, x_lab := factor(x_lab, levels = unique(c(stat_dt$neg_lab, stat_dt$pos_lab)))]
  long[, neglogp := ifelse(is.finite(pvalue), pmin(10, -log10(pmax(pvalue, 1e-12))), 0.5)]
  fill_lim <- max(abs(long$median_value), na.rm = TRUE)
  if (!is.finite(fill_lim) || fill_lim < 1e-8) fill_lim <- 0.1
  p <- ggplot(long, aes(x = x_lab, y = y_lab)) +
    geom_point(
      aes(size = neglogp, fill = median_value),
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
    scale_x_discrete(expand = expansion(add = 0.85)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1, size = 11),
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      strip.text = element_text(size = 10),
      panel.spacing.x = grid::unit(0.55, "lines")
    )
  if (isTRUE(facet) && "panel" %in% names(long) && uniqueN(long$panel) > 1) {
    p <- p + facet_wrap(~ panel, nrow = 1, scales = "free_x")
  }
  n_panel <- if (isTRUE(facet)) max(1, uniqueN(long$panel)) else 1
  n_y <- uniqueN(long$y_lab)
  fig_w <- if (isTRUE(facet)) max(9.5, 2.7 * n_panel + 3.2) else 6.8
  fig_h <- max(5.4, 0.48 * n_y + 2.4)
  save_plot(p, path_stub, fig_w, fig_h)
}

build_aurora_annotation <- function(sample_ids, clin_sample, clin_patient) {
  ann <- data.table(sample = norm_id(sample_ids))
  ann <- ann[!is.na(sample) & sample != ""]
  ann <- ann[!duplicated(sample)]

  cs <- as.data.table(clin_sample)
  names(cs) <- toupper(names(cs))
  sid <- first_present(names(cs), c("SAMPLE_ID", "SAMPLE", "SAMPLEID"))
  pid <- first_present(names(cs), c("PATIENT_ID", "PATIENT", "PATIENTID"))
  if (is.na(sid)) stop("临床样本表没有 SAMPLE_ID")
  cs[, sample := norm_id(cs[[sid]])]
  cs[, patient := if (!is.na(pid)) as.character(cs[[pid]]) else sub("-(TTP|TTM).*$", "", sample, ignore.case = TRUE)]
  st <- first_present(names(cs), c("SAMPLE_TYPE", "TUMOR_TYPE", "TISSUE_TYPE"))
  ms <- first_present(names(cs), c("METASTATIC_SITE", "METASTATIC.SITE"))
  keep_s <- unique(c("sample", "patient", na.omit(c(st, ms))))
  extra_s <- cs[, keep_s, with = FALSE]
  extra_s <- extra_s[!duplicated(sample)]
  ann <- merge(ann, extra_s, by = "sample", all.x = TRUE)

  if (!is.null(clin_patient) && nrow(clin_patient) > 0) {
    cp <- as.data.table(clin_patient)
    names(cp) <- toupper(names(cp))
    pp <- first_present(names(cp), c("PATIENT_ID", "PATIENT", "PATIENTID"))
    if (!is.na(pp)) {
      cp[, patient := as.character(cp[[pp]])]
      want <- c(
        "PRIM_M", "YPM", "PRIM_N", "YPN", "PRIM_STAGE_DX",
        "YSTAGE_AT_PSTAGING", "DFS_STATUS", "OS_STATUS"
      )
      keep_p <- unique(c("patient", intersect(want, names(cp))))
      extra_p <- cp[, keep_p, with = FALSE]
      extra_p <- extra_p[!duplicated(patient)]
      ann <- merge(ann, extra_p, by = "patient", all.x = TRUE)
    }
  }

  raw_type <- if (!is.na(st) && st %in% names(ann)) ann[[st]] else rep(NA_character_, nrow(ann))
  type_vec <- classify_sample_class(raw_type, ann$sample)
  m_raw <- if ("PRIM_M" %in% names(ann)) ann$PRIM_M else rep(NA_character_, nrow(ann))
  if ("YPM" %in% names(ann)) {
    miss <- is.na(classify_m(m_raw))
    m_raw[miss] <- ann$YPM[miss]
  }
  n_raw <- if ("YPN" %in% names(ann)) ann$YPN else rep(NA_character_, nrow(ann))
  if ("PRIM_N" %in% names(ann)) {
    miss <- is.na(classify_n(n_raw))
    n_raw[miss] <- ann$PRIM_N[miss]
  }
  st_raw <- if ("PRIM_STAGE_DX" %in% names(ann)) ann$PRIM_STAGE_DX else rep(NA_character_, nrow(ann))
  if ("YSTAGE_AT_PSTAGING" %in% names(ann)) {
    miss <- is.na(classify_stage(st_raw))
    st_raw[miss] <- ann$YSTAGE_AT_PSTAGING[miss]
  }

  ann[, sample_class := factor(type_vec, levels = c("Primary", "Metastatic"))]
  ann[, distant_M := factor(classify_m(m_raw), levels = c("M0", "M1"))]
  ann[, node_N := factor(classify_n(n_raw), levels = c("N0", "Nplus"))]
  ann[, stage_IV := factor(classify_stage(st_raw), levels = c("Stage I-III", "Stage IV"))]
  if (is.na(ann$patient[1]) || !nzchar(ann$patient[1])) {
    ann[, patient := sub("-(TTP|TTM).*$", "", sample, ignore.case = TRUE)]
  }
  message(
    "选用临床：SAMPLE_TYPE=", st,
    "  M=PRIM_M/YPM  N=YPN/PRIM_N  stage=PRIM_STAGE_DX/YSTAGE"
  )
  ann
}

# ==============================================================================
# 读数据（分析函数在 expr 建好之后才调用）
# ==============================================================================
aurora_work_dir <- resolve_aurora_dir()
setwd(aurora_work_dir)
message("数据目录：", aurora_work_dir)

expr_file <- pick_expression_file(aurora_work_dir)
samp_file <- pick_clin_file(aurora_work_dir, "data_clinical_sample")
pat_file <- pick_clin_file(aurora_work_dir, "data_clinical_patient")
if (is.na(samp_file)) stop("找不到 data_clinical_sample.txt")

message("表达：", expr_file)
message("样本临床：", samp_file)
if (!is.na(pat_file)) message("患者临床：", pat_file)

expr_dt <- fread_cbioportal(expr_file)
clin_sample <- fread_cbioportal(samp_file)
clin_patient <- if (!is.na(pat_file)) fread_cbioportal(pat_file) else NULL

dir.create(aurora_out_dir, showWarnings = FALSE, recursive = TRUE)

message("正在把表达表映射为基因符号矩阵…")
aurora_expr_all <- as_symbol_matrix_cbioportal(expr_dt)
aurora_already_z <- detect_already_z(aurora_expr_all, expr_file)
mx <- suppressWarnings(max(aurora_expr_all, na.rm = TRUE))
if (!isTRUE(aurora_already_z) && is.finite(mx) && mx > 50) {
  message("检测到原始表达（max=", round(mx, 2), "），做 log2(x+1)")
  aurora_expr_all <- log2(pmax(aurora_expr_all, 0) + 1)
} else if (isTRUE(aurora_already_z)) {
  message("表达已是 z-score（", basename(expr_file), "），z-mean/z-median 不再重复标准化")
} else {
  message("表达值范围较小（max=", round(mx, 2), "），视为已 log 转换")
}

aurora_ann <- build_aurora_annotation(colnames(aurora_expr_all), clin_sample, clin_patient)
aurora_ann <- aurora_ann[sample %in% colnames(aurora_expr_all)]
fwrite(aurora_ann, file.path(aurora_out_dir, "00_sample_annotation.csv"))

message(
  "样本：原发=", sum(aurora_ann$sample_class == "Primary", na.rm = TRUE),
  "  转移组织=", sum(aurora_ann$sample_class == "Metastatic", na.rm = TRUE),
  "  M1=", sum(aurora_ann$distant_M == "M1", na.rm = TRUE),
  "  Stage IV=", sum(aurora_ann$stage_IV == "Stage IV", na.rm = TRUE),
  "  N+=", sum(aurora_ann$node_N == "Nplus", na.rm = TRUE)
)

keep_g <- rowMeans(is.finite(aurora_expr_all), na.rm = TRUE) >= min_expr_frac
if (!isTRUE(aurora_already_z)) {
  keep_g <- keep_g & rowMeans(is.finite(aurora_expr_all) & aurora_expr_all > 0, na.rm = TRUE) >= min_expr_frac
}
aurora_expr <- aurora_expr_all[keep_g, , drop = FALSE]
primary_ids <- aurora_ann$sample[aurora_ann$sample_class == "Primary"]
met_ids <- aurora_ann$sample[aurora_ann$sample_class == "Metastatic"]
aurora_expr_primary <- aurora_expr[, intersect(primary_ids, colnames(aurora_expr)), drop = FALSE]
message(
  "表达矩阵：全部 ", ncol(aurora_expr), " 样本 x ", nrow(aurora_expr),
  " 基因；原发 ", ncol(aurora_expr_primary)
)
if (ncol(aurora_expr) < 10) stop("有表达的样本太少：", ncol(aurora_expr))

# ==============================================================================
# 主分析（必须在 expr 建好之后）
# ==============================================================================
run_aurora_nerve <- function() {
  if (!exists("aurora_expr", inherits = TRUE)) stop("还没有 aurora_expr，请从脚本开头 Source")

  mat_from_list <- function(lst) {
    if (length(lst) == 0) return(NULL)
    common <- Reduce(intersect, lapply(lst, names))
    if (length(common) == 0) return(NULL)
    mat <- do.call(cbind, lapply(lst, function(x) x[common]))
    colnames(mat) <- names(lst)
    rownames(mat) <- common
    mat
  }

  message("收集五个神经 GO 基因（每个 GO 单独，不合并）")
  go_map <- lapply(go_list, get_go_genes)
  names(go_map) <- go_list
  go_map <- Filter(function(x) length(x$genes) >= min_set_genes, go_map)
  if (length(go_map) == 0) stop("没有任何神经 GO 能打分")
  gene_sets <- lapply(go_map, function(x) x$genes)
  fwrite(rbindlist(lapply(names(go_map), function(g) {
    data.table(
      GO = g, GO_name = go_title(g), gene_source = go_map[[g]]$source,
      n_genes = length(go_map[[g]]$genes),
      genes = paste(go_map[[g]]$genes, collapse = ";")
    )
  }), fill = TRUE), file.path(aurora_out_dir, "00_GO_genes_used.csv"))
  fwrite(data.table(
    method = c("zmean", "zmedian", "ssgsea"),
    role = c("primary", "supplement", "supplement"),
    description = c(
      "Gene-wise z-score then mean (skipped if matrix is already z-scored)",
      "Gene-wise z-score then median (skipped if matrix is already z-scored)",
      "ssGSEA (Barbie 2009); GSVA if installed, otherwise built-in"
    )
  ), file.path(aurora_out_dir, "00_scoring_methods.csv"))

  score_z_method <- function(expr_mat, how) {
    lst <- lapply(names(gene_sets), function(g) {
      fn <- if (identical(how, "median")) pathway_zmedian else pathway_zmean
      set_score <- fn(expr_mat, gene_sets[[g]], already_z = aurora_already_z)
      if (is.null(set_score)) {
        message("  ", g, " 映射基因不足，跳过")
        return(NULL)
      }
      message("  ", g, "  ", go_title(g), "  基因数=", attr(set_score, "n_genes"))
      set_score
    })
    names(lst) <- names(gene_sets)
    Filter(Negate(is.null), lst)
  }

  plot_one_method <- function(method, sm_all) {
    if (is.null(sm_all)) {
      message("  ", method$title, " 没有可画的分数，跳过")
      return(invisible(NULL))
    }
    mdir <- file.path(aurora_out_dir, method$id)
    dir.create(mdir, showWarnings = FALSE, recursive = TRUE)
    fwrite(data.table(sample = rownames(sm_all), as.data.table(sm_all)),
           file.path(mdir, "01_pathway_scores.csv"))

    ann_use <- aurora_ann[sample %in% rownames(sm_all)]
    ann_p <- ann_use[sample_class == "Primary"]
    sm_p <- sm_all[intersect(rownames(sm_all), ann_p$sample), , drop = FALSE]

    designs <- list(
      list(key = "a_distant_M", title = "Distant M at first diagnosis", panel = "1a Distant M",
           group = setNames(as.character(ann_p$distant_M), ann_p$sample),
           pos = "M1", neg = "M0", pos_lab = "M1", neg_lab = "M0",
           score_mat = sm_p, paired = FALSE),
      list(key = "b_AJCC_stageIV", title = "AJCC stage at diagnosis", panel = "1b AJCC stage",
           group = setNames(as.character(ann_p$stage_IV), ann_p$sample),
           pos = "Stage IV", neg = "Stage I-III",
           pos_lab = "Stage IV", neg_lab = "Stage I-III",
           score_mat = sm_p, paired = FALSE),
      list(key = "c_node_N", title = "Lymph node", panel = "1c Lymph node",
           group = setNames(as.character(ann_p$node_N), ann_p$sample),
           pos = "Nplus", neg = "N0", pos_lab = "N+", neg_lab = "N0",
           score_mat = sm_p, paired = FALSE),
      list(key = "d_sample_type", title = "Sample type", panel = "1d Sample type",
           group = setNames(as.character(ann_use$sample_class), ann_use$sample),
           pos = "Metastatic", neg = "Primary",
           pos_lab = "Metastatic tissue", neg_lab = "Primary tumor",
           score_mat = sm_all, paired = FALSE),
      list(key = "e_paired", title = "Paired primary vs metastasis", panel = "1e Paired",
           group = NULL, pos = "Metastatic", neg = "Primary",
           pos_lab = "Metastatic tissue", neg_lab = "Primary tumor",
           score_mat = sm_all, paired = TRUE)
    )

    all_stat <- list()
    for (ds in designs) {
      message("气泡图：", method$title, " / ", ds$title)
      stat_rows <- list()
      sm <- ds$score_mat
      if (is.null(sm) || nrow(sm) == 0) {
        message("  无分数矩阵，跳过 ", ds$key)
        next
      }
      for (g in colnames(sm)) {
        value_vec <- as.numeric(sm[, g])
        names(value_vec) <- rownames(sm)
        one <- if (isTRUE(ds$paired)) {
          compare_paired_patients(value_vec, ann_use, ds$key)
        } else {
          compare_groups(value_vec, ds$group[names(value_vec)], ds$pos, ds$neg, ds$key)
        }
        if (is.null(one)) next
        if (!isTRUE(ds$paired) && (one$n_pos < min_group_n || one$n_neg < min_group_n)) {
          next
        }
        one[, `:=`(
          scoring = method$id, GO = g, GO_name = go_title(g),
          panel = ds$panel, pos_lab = ds$pos_lab, neg_lab = ds$neg_lab
        )]
        stat_rows[[g]] <- one
      }
      stat_dt <- rbindlist(stat_rows, fill = TRUE)
      if (nrow(stat_dt) == 0) {
        message("  分组人数不足，跳过 ", method$id, " / ", ds$key)
        next
      }
      stat_dt[, fdr := p.adjust(pvalue, method = "BH")]
      fwrite(stat_dt, file.path(mdir, paste0("02_", ds$key, "_GO_vs_metastasis.csv")))
      fwrite(expand_two_cols(stat_dt),
             file.path(mdir, paste0("02_", ds$key, "_GO_vs_metastasis_two_cols.csv")))
      all_stat[[ds$key]] <- stat_dt
      plot_bubble_two_cols(
        stat_dt,
        title = paste0(ds$title, ": neural GO (", method$title, ")"),
        subtitle = paste0(
          "AURORA US; scoring = ", method$title,
          "; Y = neuronal GO (scored separately); X = two groups; fill = median; size = -log10(p)"
        ),
        path_stub = file.path(mdir, paste0("02_", ds$key, "_bubble"))
      )
    }

    if (length(all_stat) == 0) return(invisible(NULL))
    bubble <- rbindlist(all_stat, fill = TRUE)
    fwrite(bubble, file.path(mdir, "02_summary_GO_vs_metastasis.csv"))
    fwrite(expand_two_cols(bubble),
           file.path(mdir, "02_summary_GO_vs_metastasis_two_cols.csv"))
    plot_bubble_two_cols(
      bubble,
      title = paste0("AURORA neural GO versus metastasis (", method$title, ")"),
      subtitle = paste0(
        "Scoring = ", method$title,
        "; five GOs scored separately and not pooled"
      ),
      path_stub = file.path(mdir, "02_summary_bubble_GO_vs_metastasis"),
      facet = TRUE
    )
    if (isTRUE(method$primary)) {
      fwrite(bubble, file.path(aurora_out_dir, "02_summary_GO_vs_metastasis.csv"))
      plot_bubble_two_cols(
        bubble,
        title = "AURORA neural GO versus metastasis (z-mean, primary)",
        subtitle = "Primary scoring = z-mean; five GOs not pooled",
        path_stub = file.path(aurora_out_dir, "02_summary_bubble_GO_vs_metastasis"),
        facet = TRUE
      )
    }
    bubble
  }

  score_methods <- list(
    list(id = "zmean", title = "z-mean", primary = TRUE),
    list(id = "zmedian", title = "z-median", primary = FALSE),
    list(id = "ssgsea", title = "ssGSEA", primary = FALSE)
  )
  all_method_stats <- list()
  for (method in score_methods) {
    message("打分：", method$title, if (isTRUE(method$primary)) "（主）" else "（补充）")
    if (identical(method$id, "ssgsea")) {
      lst <- score_ssgsea_sets(aurora_expr, gene_sets)
    } else {
      how <- if (identical(method$id, "zmedian")) "median" else "mean"
      lst <- score_z_method(aurora_expr, how)
    }
    sm_all <- mat_from_list(lst)
    all_method_stats[[method$id]] <- plot_one_method(method, sm_all)
  }

  keep <- Filter(function(x) is.data.table(x) && nrow(x) > 0, all_method_stats)
  if (length(keep) > 0) {
    fwrite(rbindlist(keep, fill = TRUE),
           file.path(aurora_out_dir, "03_all_scoring_vs_metastasis.csv"))
  }

  message("完成。结果目录：", normalizePath(aurora_out_dir, winslash = "/", mustWork = FALSE))
  message("主图（转移组织 vs 原发）：zmean/02_d_sample_type_bubble.png")
  message("配对图：zmean/02_e_paired_bubble.png")
  invisible(TRUE)
}

if (exists("aurora_expr") && ncol(aurora_expr) > 5) {
  run_aurora_nerve()
} else {
  stop("表达矩阵未建好，请从第一行完整 Source")
}
