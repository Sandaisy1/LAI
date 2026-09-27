################################################################################
# nerve_TCGA.R
# 单独脚本，不要和 TCGA-BRCA-2026.R 在同一会话 Source
# 放到 E:/R/TCGA-BRCA-2026，RStudio 打开后从第一行 Source
#
# 数据（.tsv / .tsv.gz 均可；空的 .tsv 会跳过，改读 .gz）：
#   TCGA-BRCA.star_fpkm.tsv(.gz)
#   TCGA-BRCA.clinical.tsv(.gz)
#   TCGA-BRCA.survival.tsv(.gz)          # 可选
#   gencode.v36.annotation.gtf.gene.probemap
#
# 神经浸润只按下面 5 个神经信号 GO，每个单独取基因、单独打分（不合并）：
#   GO:0019227  neuronal action potential propagation
#   GO:1902847  regulation of neuronal signal transduction
#   GO:0097374  sensory neuron axon guidance
#   GO:1902667  regulation of axon guidance
#   GO:0007409  axonogenesis
#
# 对 4 种转移定义做 Wilcoxon，并画两列气泡图（Non-metastatic | Metastatic）
#   a 原位肿瘤 M1 vs M0
#   b 原位肿瘤 Stage IV vs I-III
#   c 原位肿瘤 N+ vs N0
#   d 转移组织 vs 原位肿瘤
#
# 结果目录：results_nerve_TCGA/
################################################################################

library(data.table)
library(ggplot2)

nerve_work_dir <- Sys.getenv(
  "NERVE_TCGA_DIR",
  unset = Sys.getenv("TCGA_BRCA_2026_DIR", unset = "E:/R/TCGA-BRCA-2026")
)
nerve_out_dir <- "results_nerve_TCGA"
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
  ifelse(nchar(x) >= 15, substr(x, 14, 15), x)
}
first_present <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) NA_character_ else hit[1]
}

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
pick_best_table <- function(stems, must = TRUE, min_cols = 1L, label = "table") {
  hits <- existing_files(stems)
  if (length(hits) == 0) {
    if (must) stop("找不到", label, "：", paste(stems, collapse = " / "))
    return(NA_character_)
  }
  tab <- rbindlist(lapply(hits, probe_table_file), fill = TRUE)
  message(label, "候选：\n",
          paste(sprintf("  %s  大小=%.1fMB  列=%d  可读=%s",
                        tab$file, tab$size / 1024^2, tab$n_col, tab$readable),
                collapse = "\n"))
  ok <- tab[readable == TRUE & n_col >= min_cols]
  if (nrow(ok) == 0) {
    if (must) stop(label, "存在但读不了。请删掉空的 .tsv，保留 .tsv.gz 后再 Source。")
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
  stop("读不了 ", path, "：", conditionMessage(out))
}
pick_best_clinical <- function() {
  f <- pick_best_table(
    c("TCGA-BRCA.clinical.tsv", "TCGA-BRCA.GDC_phenotype.tsv"),
    must = TRUE, min_cols = 5L, label = "临床表"
  )
  hdr <- names(fread(f, nrows = 0, fill = TRUE, showProgress = FALSE))
  message("选用临床：", f, "  列=", length(hdr))
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
  out[grepl("PRIMARY", xl)] <- "PrimaryTumor"
  out[grepl("METASTATIC", xl)] <- "MetastaticTissue"
  out[grepl("NORMAL", xl)] <- "Normal"
  code <- sample_type_code(barcode)
  out[is.na(out) & code == "01"] <- "PrimaryTumor"
  out[is.na(out) & code %in% c("06", "07")] <- "MetastaticTissue"
  out[is.na(out) & code == "11"] <- "Normal"
  out
}

# 通路活性：基因 z 后对样本取均值。不要用 scale()/t()，变量不要叫 score
pathway_zmean <- function(expr_mat, genes) {
  genes <- unique(intersect(as.character(genes), rownames(expr_mat)))
  if (length(genes) < min_set_genes) return(NULL)
  sub <- as.matrix(expr_mat[genes, , drop = FALSE])
  storage.mode(sub) <- "double"
  gene_mean <- rowMeans(sub, na.rm = TRUE)
  gene_sd <- sqrt(rowMeans((sub - gene_mean)^2, na.rm = TRUE))
  gene_sd[!is.finite(gene_sd) | gene_sd < 1e-12] <- 1
  z <- (sub - gene_mean) / gene_sd
  z[!is.finite(z)] <- 0
  set_score <- colMeans(z, na.rm = TRUE)
  names(set_score) <- colnames(sub)
  attr(set_score, "n_genes") <- length(genes)
  attr(set_score, "genes") <- genes
  set_score
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
    if (length(syms) > 0) {
      return(list(genes = syms, source = "org.Hs.eg.db"))
    }
  }
  fb <- go_fallback[[go_id]]
  if (is.null(fb) || length(fb) == 0) {
    return(list(genes = character(), source = "none"))
  }
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

build_annotation <- function(sample_ids, clin) {
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

  keep_cols <- unique(c("sample_std", "patient_std", na.omit(c(st, mc, nc, tc))))
  extra <- clin[, keep_cols, with = FALSE]
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
  ann[, sample_class := factor(type_vec, levels = c("PrimaryTumor", "MetastaticTissue", "Normal"))]
  ann
}

# ==============================================================================
# 读数据（分析函数在 expr 建好之后才调用）
# ==============================================================================
if (dir.exists(nerve_work_dir)) {
  setwd(nerve_work_dir)
} else {
  message("未找到 ", nerve_work_dir, " ，改用当前目录：", getwd())
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

dir.create(nerve_out_dir, showWarnings = FALSE, recursive = TRUE)

message("正在把 FPKM 映射为基因符号矩阵…")
nerve_expr_all <- as_symbol_matrix(fpkm_data, probe_annot)
mx <- suppressWarnings(max(nerve_expr_all, na.rm = TRUE))
if (is.finite(mx) && mx > 50) {
  message("检测到原始 FPKM（max=", round(mx, 2), "），做 log2(x+1)")
  nerve_expr_all <- log2(nerve_expr_all + 1)
} else {
  message("表达值范围较小（max=", round(mx, 2), "），视为已 log 转换")
}

nerve_ann <- build_annotation(colnames(nerve_expr_all), clinical_data)
nerve_ann <- nerve_ann[sample %in% colnames(nerve_expr_all)]
fwrite(nerve_ann, file.path(nerve_out_dir, "00_sample_annotation.csv"))

message(
  "样本：原位=", sum(nerve_ann$sample_class == "PrimaryTumor", na.rm = TRUE),
  "  转移组织=", sum(nerve_ann$sample_class == "MetastaticTissue", na.rm = TRUE),
  "  正常=", sum(nerve_ann$sample_class == "Normal", na.rm = TRUE),
  "  M1=", sum(nerve_ann$distant_M == "M1", na.rm = TRUE),
  "  Stage IV=", sum(nerve_ann$stage_simplified == "Stage IV", na.rm = TRUE),
  "  N+=", sum(nerve_ann$node_N == "Nplus", na.rm = TRUE)
)

primary_ids <- nerve_ann$sample[nerve_ann$sample_class == "PrimaryTumor"]
met_ids <- nerve_ann$sample[nerve_ann$sample_class == "MetastaticTissue"]
tumor_ids <- unique(c(primary_ids, met_ids))
if (length(primary_ids) < 20) stop("原位肿瘤样本太少：", length(primary_ids))

nerve_expr_tumor <- nerve_expr_all[, intersect(tumor_ids, colnames(nerve_expr_all)), drop = FALSE]
keep_g <- rowMeans(is.finite(nerve_expr_tumor) & nerve_expr_tumor > 0, na.rm = TRUE) >= min_expr_frac
nerve_expr_tumor <- nerve_expr_tumor[keep_g, , drop = FALSE]
nerve_expr_primary <- nerve_expr_tumor[, intersect(primary_ids, colnames(nerve_expr_tumor)), drop = FALSE]
message(
  "表达矩阵：原位 ", ncol(nerve_expr_primary), " 样本 x ", nrow(nerve_expr_primary),
  " 基因；肿瘤(含转移组织) ", ncol(nerve_expr_tumor)
)

# ==============================================================================
# 主分析（必须在 expr 建好之后）
# ==============================================================================
run_nerve_tcga <- function() {
  if (!exists("nerve_expr_primary", inherits = TRUE)) {
    stop("还没有 nerve_expr_primary，请从脚本开头 Source")
  }

  score_one_go <- function(expr_mat, go_id) {
    gm <- get_go_genes(go_id)
    set_score <- pathway_zmean(expr_mat, gm$genes)
    if (is.null(set_score)) {
      message("  ", go_id, " 映射基因不足，跳过")
      return(NULL)
    }
    attr(set_score, "gene_source") <- gm$source
    message("  ", go_id, "  ", go_title(go_id), "  基因数=", attr(set_score, "n_genes"),
            "  (", gm$source, ")")
    set_score
  }

  mat_from_list <- function(lst) {
    if (length(lst) == 0) return(NULL)
    common <- Reduce(intersect, lapply(lst, names))
    if (length(common) == 0) return(NULL)
    mat <- do.call(cbind, lapply(lst, function(x) x[common]))
    colnames(mat) <- names(lst)
    rownames(mat) <- common
    mat
  }

  message("计算各神经 GO 通路分数（每个 GO 单独，不合并）")
  go_primary <- lapply(go_list, function(g) score_one_go(nerve_expr_primary, g))
  names(go_primary) <- go_list
  go_primary <- Filter(Negate(is.null), go_primary)
  if (length(go_primary) == 0) stop("没有任何神经 GO 能打分")

  go_tumor <- lapply(names(go_primary), function(g) score_one_go(nerve_expr_tumor, g))
  names(go_tumor) <- names(go_primary)
  go_tumor <- Filter(Negate(is.null), go_tumor)

  sm_p <- mat_from_list(go_primary)
  sm_t <- mat_from_list(go_tumor)
  fwrite(data.table(sample = rownames(sm_p), as.data.table(sm_p)),
         file.path(nerve_out_dir, "01_pathway_scores_primary.csv"))
  fwrite(data.table(sample = rownames(sm_t), as.data.table(sm_t)),
         file.path(nerve_out_dir, "01_pathway_scores_tumor.csv"))

  gene_dt <- rbindlist(lapply(names(go_primary), function(g) {
    data.table(
      GO = g, GO_name = go_title(g),
      gene_source = attr(go_primary[[g]], "gene_source"),
      n_genes = attr(go_primary[[g]], "n_genes"),
      genes = paste(attr(go_primary[[g]], "genes"), collapse = ";")
    )
  }), fill = TRUE)
  fwrite(gene_dt, file.path(nerve_out_dir, "00_GO_genes_used.csv"))

  ann_p <- nerve_ann[sample %in% rownames(sm_p)]
  ann_t <- nerve_ann[sample %in% rownames(sm_t)]
  designs <- list(
    list(key = "a_distant_M", title = "Distant metastasis", panel = "1a Distant M",
         group = setNames(as.character(ann_p$distant_M), ann_p$sample),
         pos = "M1", neg = "M0", pos_lab = "M1", neg_lab = "M0",
         score_mat = sm_p),
    list(key = "b_AJCC_stageIV", title = "AJCC stage", panel = "1b AJCC stage",
         group = setNames(as.character(ann_p$stage_IV), ann_p$sample),
         pos = "Stage IV", neg = "Stage I-III",
         pos_lab = "Stage IV", neg_lab = "Stage I-III",
         score_mat = sm_p),
    list(key = "c_node_N", title = "Lymph node", panel = "1c Lymph node",
         group = setNames(as.character(ann_p$node_N), ann_p$sample),
         pos = "Nplus", neg = "N0", pos_lab = "N+", neg_lab = "N0",
         score_mat = sm_p),
    list(key = "d_sample_type", title = "Sample type", panel = "1d Sample type",
         group = setNames(as.character(ann_t$sample_class), ann_t$sample),
         pos = "MetastaticTissue", neg = "PrimaryTumor",
         pos_lab = "Metastatic tissue", neg_lab = "Primary tumor",
         score_mat = sm_t)
  )

  all_stat <- list()
  for (ds in designs) {
    message("气泡图：", ds$title)
    stat_rows <- list()
    sm <- ds$score_mat
    for (g in colnames(sm)) {
      value_vec <- as.numeric(sm[, g])
      names(value_vec) <- rownames(sm)
      one <- compare_groups(value_vec, ds$group[names(value_vec)], ds$pos, ds$neg, ds$key)
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
    fwrite(stat_dt, file.path(nerve_out_dir, paste0("02_", ds$key, "_GO_vs_metastasis.csv")))
    fwrite(expand_two_cols(stat_dt),
           file.path(nerve_out_dir, paste0("02_", ds$key, "_GO_vs_metastasis_two_cols.csv")))
    all_stat[[ds$key]] <- stat_dt
    plot_bubble_two_cols(
      stat_dt,
      title = paste0(ds$title, ": neural GO in non-metastatic vs metastatic"),
      subtitle = "Y = neuronal GO (scored separately); X = non-metastatic | metastatic; fill = median pathway score; size = -log10(Wilcoxon p)",
      path_stub = file.path(nerve_out_dir, paste0("02_", ds$key, "_bubble"))
    )
  }

  if (length(all_stat) > 0) {
    bubble <- rbindlist(all_stat, fill = TRUE)
    fwrite(bubble, file.path(nerve_out_dir, "02_summary_GO_vs_metastasis.csv"))
    fwrite(expand_two_cols(bubble),
           file.path(nerve_out_dir, "02_summary_GO_vs_metastasis_two_cols.csv"))
    plot_bubble_two_cols(
      bubble,
      title = "Neural GO activity versus breast cancer metastasis",
      subtitle = "Each panel: non-metastatic | metastatic; five GOs scored separately and not pooled",
      path_stub = file.path(nerve_out_dir, "02_summary_bubble_GO_vs_metastasis"),
      facet = TRUE
    )
  }

  message("完成。结果目录：", normalizePath(nerve_out_dir, winslash = "/", mustWork = FALSE))
  message("主气泡图：02_summary_bubble_GO_vs_metastasis.png")
  invisible(TRUE)
}

if (exists("nerve_expr_primary") && ncol(nerve_expr_primary) > 10) {
  run_nerve_tcga()
} else {
  stop("表达矩阵未建好，请从第一行完整 Source")
}
