#!/usr/bin/env Rscript
# =============================================================================
# Wang et al. Nature Communications 2024  TNBC 空间转录组
# 神经浸润：方案 + 肿瘤细胞高表达、可能促使神经浸润的基因
#
# 本脚本独立于 TG_RNAseq_pipeline.R 的 1–4 组比较，不修改原流程。
# 数据默认：E:/R/Nerve（Wang Zenodo 解压目录）
#
# 用法（Windows R / RStudio）：
#   setwd("E:/R/Nerve")
#   source("E:/R/TG_BRCA/TG/Wang_ST_nerve_infiltration.R")   # 或把本文件拷到 E:/R/Nerve
#   # 也可：Sys.setenv(WANG_ST_DIR = "E:/R/Nerve")
#
# 比较（各自单独出表/出图；病人不偷偷合并）：
#   1) 病理 Nerve 邻域：每个有神经 spot 的病人单独 1-vs-1
#   2) 这些病人的共同上调
#   3) 全队列 Schwann 签名邻域：~ patient + group 的伪 bulk DESeq2
# 显著性：先 p < 0.05，再按上调 FC >= 1.25 / 1.5 分层（不做 FC=1、FC=2，也不做 topN）。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")

# -----------------------------------------------------------------------------
# 0. 路径
# -----------------------------------------------------------------------------
resolve_nerve_dir <- function() {
  env_dir <- Sys.getenv("WANG_ST_DIR", unset = "")
  candidates <- c(
    env_dir,
    "E:/R/Nerve",
    "E:\\R\\Nerve",
    file.path(getwd(), "data", "Wang_ST_TNBC"),
    getwd()
  )
  candidates <- unique(candidates[nzchar(candidates)])
  looks_like_wang <- function(d) {
    dir.exists(d) && (
      dir.exists(file.path(d, "Robjects", "annotsBySpot")) ||
        dir.exists(file.path(d, "annotsBySpot")) ||
        file.exists(file.path(d, "Robjects.tar")) ||
        file.exists(file.path(d, "Clinical", "ids.RDS"))
    )
  }
  for (d in candidates) {
    if (looks_like_wang(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  stop("找不到 Wang ST 数据。请把 Zenodo 文件放到 E:/R/Nerve，或设置 WANG_ST_DIR。")
}

nerve_dir <- resolve_nerve_dir()

ensure_tar_extracted <- function(root, dirname, tarname) {
  dest <- file.path(root, dirname)
  if (dir.exists(dest) && length(list.files(dest, recursive = TRUE)) > 0) return(invisible(dest))
  tarf <- file.path(root, tarname)
  if (!file.exists(tarf)) return(invisible(NULL))
  message("解压 ", tarname, " -> ", root)
  utils::untar(tarf, exdir = root)
  dest
}

ensure_tar_extracted(nerve_dir, "Clinical", "Clinical.tar")
ensure_tar_extracted(nerve_dir, "Robjects", "Robjects.tar")
ensure_tar_extracted(nerve_dir, "misc", "misc.tar")
ensure_tar_extracted(nerve_dir, "patches", "patches.tar")
ensure_tar_extracted(nerve_dir, "classification", "classification.tar")

robj_dir <- if (dir.exists(file.path(nerve_dir, "Robjects"))) {
  file.path(nerve_dir, "Robjects")
} else {
  nerve_dir
}

# -----------------------------------------------------------------------------
# 1. 可选：载入原 RNA-seq 脚本的作图/富集函数（不跑主流程）
# -----------------------------------------------------------------------------
this_script_dir <- function() {
  ofile <- NULL
  if (sys.nframe() > 0) {
    for (i in sys.nframe():1) {
      e <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
      if (!is.null(e)) {
        ofile <- e
        break
      }
    }
  }
  if (is.null(ofile)) getwd() else dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE))
}

load_pipeline_functions_only <- function() {
  cands <- c(
    file.path(this_script_dir(), "TG_RNAseq_pipeline.R"),
    "TG_RNAseq_pipeline.R",
    file.path(getwd(), "TG_RNAseq_pipeline.R"),
    "E:/R/TG_BRCA/TG/TG_RNAseq_pipeline.R",
    "E:\\R\\TG_BRCA\\TG\\TG_RNAseq_pipeline.R"
  )
  pipe <- cands[file.exists(cands)][1]
  if (is.na(pipe) || !nzchar(pipe) || !file.exists(pipe)) {
    message("未找到 TG_RNAseq_pipeline.R，将只出差异表和基础图，不出 ORA/GSEA 全套。")
    return(FALSE)
  }
  lines <- readLines(pipe, warn = FALSE)
  main_at <- grep("^# 10\\. 主流程", lines)[1]
  if (is.na(main_at) || main_at < 2) stop("无法从 TG_RNAseq_pipeline.R 切出函数定义")
  eval(parse(text = lines[seq_len(main_at - 1)]), envir = .GlobalEnv)
  TRUE
}

have_pipeline <- FALSE
skip_pipe <- identical(Sys.getenv("WANG_ST_SKIP_PIPELINE"), "1")
if (skip_pipe) {
  message("WANG_ST_SKIP_PIPELINE=1：跳过原流程富集函数")
} else if (!exists("analyze_one_comparison", mode = "function")) {
  have_pipeline <- tryCatch(load_pipeline_functions_only(), error = function(e) {
    message("载入原流程函数失败: ", e$message)
    FALSE
  })
} else {
  have_pipeline <- TRUE
}

if (!have_pipeline) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    tryCatch(BiocManager::install("limma", update = FALSE, ask = FALSE),
             error = function(e) message("limma install failed: ", e$message))
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    tryCatch(install.packages("ggplot2", repos = "https://cloud.r-project.org"),
             error = function(e) message("ggplot2 install failed: ", e$message))
  }
}

if (!exists("log_msg", mode = "function")) {
  log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
}

result_dir <- file.path(nerve_dir, "results")
log_dir <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("wang_st_nerve_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg_file <- log_msg
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

# 本空间分析：先滤 p < 0.05，再只按上调 FC >= 1.25 / 1.5 分层
# （不做 FC=1、FC=2，也不做 topN；不改原 Cuffdiff 脚本的档位）
if (exists("padj_cutoff")) padj_cutoff <<- 0.05
p_cutoff <- 0.05
fc_cutoffs_local <- c("FC_1.25" = 1.25, "FC_1.5" = 1.5)

log_msg("Wang ST nerve dir: ", nerve_dir)
log_msg("Results: ", result_dir)

# -----------------------------------------------------------------------------
# 2. 分子名单（文献 / 轴突导向 / Schwann；平台上没有的基因会记成未检出）
# -----------------------------------------------------------------------------
schwann_genes <- c("SOX10", "MPZ", "PMP22", "S100B", "PLP1", "NGFR", "NCAM1", "MBP", "L1CAM")
neuron_genes  <- c("PRPH", "TUBB3", "UCHL1", "NRXN1", "SNAP25", "TH", "CALCA", "RET", "NEFL", "GAP43")
tumor_epi_genes <- c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1", "KRT14")

# 肿瘤细胞 → 神经 的配体（促使轴突生长 / PNI 的候选）
tumor_ligands <- data.frame(
  gene = c(
    "NGF", "BDNF", "NTF3", "NTF4", "GDNF", "ARTN", "NRTN", "PSPN",
    "NTN1", "NTN4", "SLIT2", "SLIT3", "SEMA3A", "SEMA3C", "SEMA3F",
    "EFNA1", "EFNA5", "EFNB2", "NRG1", "CXCL12", "VEGFA", "MDK", "PTN",
    "TGFB1", "WNT5A", "IL6", "CALCA", "ADM"
  ),
  receptor_on_nerve = c(
    "NTRK1/NGFR", "NTRK2", "NTRK3", "NTRK2", "GFRA1/RET", "GFRA3/RET", "GFRA2/RET", "GFRA4/RET",
    "DCC/UNC5", "UNC5/ITGB4", "ROBO1/2", "ROBO", "NRP1/PLXNA", "NRP2/PLXNA", "NRP2",
    "EPHA", "EPHA", "EPHB", "ERBB2/3/4", "CXCR4", "VEGFR/NRP1", "PTPRZ1/ALK", "PTPRZ1",
    "TGFBR", "FZD/ROR2", "IL6R/IL6ST", "RAMP1/CALCRL", "CALCRL"
  ),
  rationale = c(
    "感觉神经营养因子，TNBC PNI / NGF–CGRP 轴",
    "TrkB 介导的神经浸润与存活",
    "TrkC 神经营养",
    "NT-4/TrkB",
    "GDNF 家族，促轴突与 PNI（胰腺癌证据强）",
    "ARTN–GFRA3，促神经浸润",
    "neurturin / GFRA2",
    "persephin / GFRA4",
    "netrin-1 轴突导向",
    "netrin-4",
    "SLIT2–ROBO 导向 / 部分促侵袭",
    "SLIT3",
    "SEMA3A 轴突导向（本阵列常未检出）",
    "SEMA3C 导向与侵袭",
    "SEMA3F 导向",
    "ephrin-A1",
    "ephrin-A5",
    "ephrin-B2 血管/神经",
    "NRG1–ERBB，Schwann / 轴突互作",
    "CXCL12–CXCR4 趋化，与神经亲和有关",
    "VEGF 同时促血管与神经",
    "midkine 神经突生长",
    "pleiotrophin 神经突生长",
    "TGF-β / EMT，利于沿神经侵袭",
    "WNT5A 非经典 Wnt 与运动",
    "IL6 神经损伤炎症（本阵列常未检出）",
    "CGRP，感觉神经肽（本阵列常未检出）",
    "肾上腺髓质素，与 CGRP 受体共用 RAMP"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 3. 读入
# -----------------------------------------------------------------------------
safe_read_rds <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) {
    log_msg("readRDS failed: ", path, " | ", e$message)
    NULL
  })
}

load_ids <- function() {
  f <- file.path(nerve_dir, "Clinical", "ids.RDS")
  ids <- safe_read_rds(f)
  if (is.null(ids)) stop("缺少 Clinical/ids.RDS")
  ids
}

load_counts <- function(pid) {
  f1 <- file.path(robj_dir, "counts", paste0("TNBC", pid, ".RDS"))
  f2 <- file.path(robj_dir, "countsNonCorrected", paste0("TNBC", pid, ".RDS"))
  obj <- safe_read_rds(f1)
  if (is.null(obj)) obj <- safe_read_rds(f2)
  if (is.null(obj) || is.null(obj$cnts) || is.null(obj$spots)) return(NULL)
  cnts <- as.matrix(obj$cnts)
  spots <- obj$spots
  if (nrow(cnts) != nrow(spots)) {
    log_msg("TNBC", pid, " cnts/spots 行数不一致，跳过")
    return(NULL)
  }
  sid <- if ("slide" %in% names(spots)) as.character(spots$slide) else "slide1"
  sid[is.na(sid) | sid == ""] <- "slide1"
  spot_id <- paste(sid, paste(spots$x, spots$y, sep = "x"), sep = ".")
  rownames(cnts) <- spot_id
  spots$spot_id <- spot_id
  spots$patient <- as.character(pid)
  spots$slide_id <- sid
  storage.mode(cnts) <- "double"
  cnts[is.na(cnts)] <- 0
  list(cnts = cnts, spots = spots)
}

load_annots <- function(pid) {
  f <- file.path(robj_dir, "annotsBySpot", paste0("TNBC", pid, ".RDS"))
  obj <- safe_read_rds(f)
  if (is.null(obj) || is.null(obj$annots)) return(NULL)
  ann <- as.matrix(obj$annots)
  sp <- obj$spots
  key <- if (!is.null(rownames(ann))) rownames(ann) else paste(sp$x, sp$y, sep = "x")
  list(annots = ann, spots = sp, key = key)
}

list_patient_ids <- function() {
  d <- file.path(robj_dir, "counts")
  if (!dir.exists(d)) d <- file.path(robj_dir, "countsNonCorrected")
  fs <- list.files(d, pattern = "^TNBC[0-9]+\\.RDS$", full.names = FALSE)
  sort(as.integer(sub("^TNBC([0-9]+)\\.RDS$", "\\1", fs)))
}

# -----------------------------------------------------------------------------
# 4. 分数、邻域、分类
# -----------------------------------------------------------------------------
present_genes <- function(genes, universe) intersect(unique(genes), universe)

log1p_cpm <- function(cnts) {
  lib <- pmax(rowSums(cnts), 1)
  log2(sweep(cnts, 1, lib, "/") * 1e4 + 1)
}

module_score <- function(logmat, genes) {
  g <- present_genes(genes, colnames(logmat))
  if (length(g) < 2) return(rep(NA_real_, nrow(logmat)))
  z <- scale(logmat[, g, drop = FALSE])
  z[is.na(z)] <- 0
  rowMeans(z)
}

spot_ann_frac <- function(ann) {
  tot <- pmax(rowSums(ann), 1)
  col <- function(nm) if (nm %in% colnames(ann)) ann[, nm] else 0
  data.frame(
    key = rownames(ann),
    nerve_px = col("Nerve"),
    tumor_px = col("Tumor") + col("Tumor region"),
    lymph_px = col("Lymphocyte") + 0.5 * col("High TIL stroma"),
    stroma_px = col("Stroma cell") + col("Acellular stroma") +
      0.5 * col("Low TIL stroma") + 0.5 * col("High TIL stroma"),
    frac_nerve = col("Nerve") / tot,
    frac_tumor = (col("Tumor") + col("Tumor region")) / tot,
    frac_lymph = (col("Lymphocyte") + 0.5 * col("High TIL stroma")) / tot,
    stringsAsFactors = FALSE
  )
}

min_dist_to <- function(xy, idx_from, idx_to) {
  if (length(idx_to) == 0) return(rep(Inf, length(idx_from)))
  src <- as.matrix(xy[idx_from, , drop = FALSE])
  tgt <- as.matrix(xy[idx_to, , drop = FALSE])
  apply(src, 1, function(p) min(sqrt((p[1] - tgt[, 1])^2 + (p[2] - tgt[, 2])^2)))
}

classify_one_patient <- function(pid, ids, near_d = 2, far_d = 4,
                                 tumor_frac = 0.25, nerve_px_min = 1) {
  cc <- load_counts(pid)
  if (is.null(cc)) return(NULL)
  spots <- cc$spots
  cnts <- cc$cnts
  keep <- rowSums(cnts) >= 500
  if (sum(keep) < 50) {
    log_msg("TNBC", pid, " 有效 spot 太少，跳过")
    return(NULL)
  }
  spots <- spots[keep, , drop = FALSE]
  cnts <- cnts[keep, , drop = FALSE]
  logmat <- log1p_cpm(cnts)

  spots$score_schwann <- module_score(logmat, schwann_genes)
  spots$score_neuron  <- module_score(logmat, neuron_genes)
  spots$score_tumor   <- module_score(logmat, tumor_epi_genes)

  spots$is_nerve_path <- FALSE
  spots$frac_tumor <- NA_real_
  spots$frac_nerve <- NA_real_
  spots$frac_lymph <- NA_real_

  aa <- load_annots(pid)
  annot_slide <- NA_character_
  if (!is.null(ids) && "hasAnnot" %in% names(ids) && "id" %in% names(ids)) {
    w <- which(as.character(ids$id) == as.character(pid) & ids$hasAnnot)
    if (length(w) > 0) annot_slide <- rownames(ids)[w[1]]
  }

  if (!is.null(aa)) {
    on_annot <- if (!is.na(annot_slide) && "slide_id" %in% names(spots)) {
      spots$slide_id == annot_slide
    } else {
      rep(TRUE, nrow(spots))
    }
    key_sp <- paste(spots$x, spots$y, sep = "x")
    fr <- spot_ann_frac(aa$annots)
    m <- match(key_sp, fr$key)
    hit <- on_annot & !is.na(m)
    spots$frac_tumor[hit] <- fr$frac_tumor[m[hit]]
    spots$frac_nerve[hit] <- fr$frac_nerve[m[hit]]
    spots$frac_lymph[hit] <- fr$frac_lymph[m[hit]]
    spots$is_nerve_path[hit] <- fr$nerve_px[m[hit]] >= nerve_px_min
  }

  # 肿瘤细胞：优先病理 Tumor 比例，否则上皮签名
  spots$is_tumor <- ifelse(
    !is.na(spots$frac_tumor),
    spots$frac_tumor >= tumor_frac,
    spots$score_tumor > 0.5 & (is.na(spots$score_schwann) | spots$score_schwann < 1)
  )
  spots$is_tumor[is.na(spots$is_tumor)] <- FALSE

  # Schwann 高：切片内标准化后 > 1 SD
  spots$is_schwann_high <- FALSE
  for (sl in unique(spots$slide_id)) {
    w <- which(spots$slide_id == sl)
    sc <- spots$score_schwann[w]
    if (all(is.na(sc))) next
    mu <- mean(sc, na.rm = TRUE)
    sdv <- stats::sd(sc, na.rm = TRUE)
    if (is.na(sdv) || sdv == 0) sdv <- 1
    spots$is_schwann_high[w] <- sc > (mu + 1 * sdv)
  }

  xy <- as.matrix(spots[, c("x", "y")])
  spots$dist_nerve <- Inf
  spots$dist_schwann <- Inf
  for (sl in unique(spots$slide_id)) {
    w <- which(spots$slide_id == sl)
    nerve_i <- w[spots$is_nerve_path[w]]
    schw_i  <- w[spots$is_schwann_high[w]]
    if (length(nerve_i) > 0) {
      spots$dist_nerve[w] <- min_dist_to(xy, w, nerve_i)
    }
    if (length(schw_i) > 0) {
      spots$dist_schwann[w] <- min_dist_to(xy, w, schw_i)
    }
  }

  spots$near_nerve <- spots$is_tumor & !spots$is_nerve_path & spots$dist_nerve <= near_d
  spots$far_nerve  <- spots$is_tumor & !spots$is_nerve_path & spots$dist_nerve > far_d
  spots$near_schwann <- spots$is_tumor & !spots$is_schwann_high & spots$dist_schwann <= near_d
  spots$far_schwann  <- spots$is_tumor & !spots$is_schwann_high & spots$dist_schwann > far_d

  list(patient = as.character(pid), cnts = cnts, logmat = logmat, spots = spots,
       annot_slide = annot_slide)
}

# -----------------------------------------------------------------------------
# 5. 差异表达
# -----------------------------------------------------------------------------
empty_de <- function() {
  data.frame(gene = character(), log2FC = numeric(), AveExpr = numeric(),
             pvalue = numeric(), padj = numeric(), padj_BH = numeric(),
             stringsAsFactors = FALSE)
}

filter_genes_mat <- function(cnts, min_total = 50, min_spots = 10) {
  keep <- colSums(cnts) >= min_total & colSums(cnts > 0) >= min_spots
  keep
}

limma_two_group <- function(logmat, group, patient = NULL, label = "") {
  if (is.null(logmat) || nrow(logmat) < 4 || ncol(logmat) < 20) return(empty_de())
  group <- factor(group, levels = c("far", "near"))
  if (nlevels(droplevels(group)) < 2) return(empty_de())
  mat <- t(logmat)
  design <- if (!is.null(patient) && length(unique(patient)) > 1) {
    stats::model.matrix(~ patient + group)
  } else {
    stats::model.matrix(~ group)
  }
  fit <- tryCatch({
    if (requireNamespace("limma", quietly = TRUE)) {
      fit0 <- limma::lmFit(mat, design)
      limma::eBayes(fit0, trend = TRUE, robust = TRUE)
    } else {
      NULL
    }
  }, error = function(e) {
    log_msg("limma failed ", label, ": ", e$message)
    NULL
  })
  if (is.null(fit)) return(empty_de())
  coefn <- grep("groupnear", colnames(design), value = TRUE)
  if (length(coefn) == 0) coefn <- colnames(design)[ncol(design)]
  tt <- limma::topTable(fit, coef = coefn, number = Inf, sort.by = "none")
  data.frame(
    gene = rownames(tt),
    log2FC = tt$logFC,
    AveExpr = tt$AveExpr,
    pvalue = tt$P.Value,
    padj = tt$P.Value,
    padj_BH = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
}

deseq2_pseudobulk <- function(pb_counts, sample_info, label = "") {
  if (is.null(pb_counts) || ncol(pb_counts) < 4) return(empty_de())
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    log_msg("DESeq2 未安装，伪 bulk 改用 limma")
    logm <- log2(sweep(t(pb_counts), 1, pmax(colSums(pb_counts), 1), "/") * 1e6 + 1)
    return(limma_two_group(logm, sample_info$group, sample_info$patient, label))
  }
  si <- sample_info
  si$group <- factor(si$group, levels = c("far", "near"))
  si$patient <- factor(si$patient)
  dds <- tryCatch({
    mode(pb_counts) <- "integer"
    dds <- DESeq2::DESeqDataSetFromMatrix(pb_counts, si, design = ~ patient + group)
    keep <- rowSums(DESeq2::counts(dds)) >= 20
    dds <- dds[keep, ]
    DESeq2::DESeq(dds, quiet = TRUE)
  }, error = function(e) {
    log_msg("DESeq2 failed ", label, ": ", e$message)
    NULL
  })
  if (is.null(dds)) {
    logm <- log2(sweep(t(pb_counts), 1, pmax(colSums(pb_counts), 1), "/") * 1e6 + 1)
    return(limma_two_group(logm, si$group, si$patient, paste(label, "limma_fallback")))
  }
  res <- DESeq2::results(dds, contrast = c("group", "near", "far"))
  data.frame(
    gene = rownames(res),
    log2FC = as.numeric(res$log2FoldChange),
    AveExpr = as.numeric(res$baseMean),
    pvalue = as.numeric(res$pvalue),
    padj = as.numeric(res$pvalue),
    padj_BH = as.numeric(res$padj),
    stringsAsFactors = FALSE
  )
}

mean_fc_two_group <- function(logmat, near_id, far_id) {
  if (length(near_id) < 3 || length(far_id) < 3) return(empty_de())
  log2FC <- colMeans(logmat[near_id, , drop = FALSE]) - colMeans(logmat[far_id, , drop = FALSE])
  data.frame(
    gene = colnames(logmat),
    log2FC = as.numeric(log2FC),
    AveExpr = colMeans(logmat[c(near_id, far_id), , drop = FALSE]),
    pvalue = NA_real_,
    padj = NA_real_,
    padj_BH = NA_real_,
    stringsAsFactors = FALSE
  )
}

prepare_de_for_pipeline <- function(de) {
  de <- de[!is.na(de$log2FC), , drop = FALSE]
  de$padj <- de$pvalue
  de
}

# -----------------------------------------------------------------------------
# 6. 作图与分层输出
# -----------------------------------------------------------------------------
basic_volcano <- function(de, title, outfile, fc_line = 1) {
  if (nrow(de) == 0) return(invisible(NULL))
  df <- de
  df$y <- if (any(!is.na(df$pvalue))) -log10(pmax(df$pvalue, 1e-300)) else abs(df$log2FC)
  df$col <- "ns"
  df$col[!is.na(df$pvalue) & df$pvalue < p_cutoff & df$log2FC >= log2(fc_line)] <- "up"
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(log2FC, y, color = col)) +
      ggplot2::geom_point(alpha = 0.5, size = 0.8) +
      ggplot2::scale_color_manual(values = c(ns = "grey70", up = "#D62828")) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(title = title, y = "-log10(p)", x = "log2FC near / far")
    ggplot2::ggsave(paste0(outfile, ".pdf"), p, width = 7, height = 6)
    ggplot2::ggsave(paste0(outfile, ".png"), p, width = 7, height = 6, dpi = 150)
  }
}

spot_map_plot <- function(spots, title, outfile) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  df <- spots
  df$cls <- "other"
  df$cls[df$is_tumor] <- "tumor_far"
  df$cls[df$near_nerve | df$near_schwann] <- "tumor_near"
  df$cls[df$is_nerve_path] <- "nerve_path"
  df$cls[df$is_schwann_high] <- "schwann_high"
  pal <- c(other = "grey85", tumor_far = "#4C78A8", tumor_near = "#D62828",
           nerve_path = "#2A9D8F", schwann_high = "#E9C46A")
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, color = cls)) +
    ggplot2::geom_point(size = 0.9) +
    ggplot2::facet_wrap(~slide_id) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title, color = NULL)
  ggplot2::ggsave(paste0(outfile, ".pdf"), p, width = 9, height = 6)
  ggplot2::ggsave(paste0(outfile, ".png"), p, width = 9, height = 6, dpi = 150)
}

emit_comparison <- function(comp_name, de, heat_mat, sample_info) {
  de <- prepare_de_for_pipeline(de)
  base <- file.path(result_dir, comp_name)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(de, file.path(base, "DE_full.csv"), row.names = FALSE)
  have_p <- any(!is.na(de$pvalue))
  log_msg(comp_name, " genes=", nrow(de), " have_p=", have_p,
          " up p<", p_cutoff, " n=", sum(!is.na(de$pvalue) & de$pvalue < p_cutoff & de$log2FC > 0))

  writeLines(
    c("本比较只做上调 FC >= 1.25 和 FC >= 1.5（先 p < 0.05）。",
      "没有 FC=1、FC=2，也没有 TopRank。",
      "分层图在 FoldChange/FC_1.25 和 FoldChange/FC_1.5。",
      "全基因 GSEA 在 00_GSEA_all_genes_NOT_FC_or_topN（不是分层图）。",
      "细胞骨架/线粒体专项在 Focused_cytoskeleton_mito/。"),
    file.path(base, "00_READ_ME_先看这里.txt")
  )
  invisible(lapply(file.path(base, "FoldChange", names(fc_cutoffs_local)),
                   dir.create, recursive = TRUE, showWarnings = FALSE))

  gsea_cache <- list()
  use_pipe_plots <- have_pipeline && exists("emit_subset_analysis", mode = "function")
  assign("result_dir", result_dir, envir = .GlobalEnv)
  assign("padj_cutoff", p_cutoff, envir = .GlobalEnv)

  for (nm in names(fc_cutoffs_local)) {
    fc <- unname(fc_cutoffs_local[[nm]])
    keep <- !is.na(de$log2FC) & de$log2FC > 0 & (2^de$log2FC >= fc)
    if (have_p) keep <- keep & !is.na(de$pvalue) & de$pvalue < p_cutoff
    sub <- de[keep, , drop = FALSE]
    if (nrow(sub) > 0) sub <- sub[order(sub$log2FC, decreasing = TRUE), , drop = FALSE]
    od <- file.path(base, "FoldChange", nm)
    if (use_pipe_plots) {
      tryCatch(
        emit_subset_analysis(
          comp_name, sub, nm, paste0(comp_name, " | up FC >= ", fc),
          od, de, heat_mat, sample_info, gsea_cache, fc_line = fc
        ),
        error = function(e) log_msg("ERROR subset ", comp_name, " ", nm, ": ", e$message)
      )
    } else {
      utils::write.csv(sub, file.path(od, paste0(nm, "_DE_selected_genes.csv")), row.names = FALSE)
      basic_volcano(de, paste(comp_name, nm), file.path(od, paste0(nm, "_volcano")), fc)
    }
  }

  if (use_pipe_plots && exists("build_gsea_cache", mode = "function")) {
    tryCatch({
      log_msg("Building full-list GSEA after subset plots: ", comp_name)
      gsea_cache <- build_gsea_cache(de)
      full_gsea_dir <- file.path(base, "00_GSEA_all_genes_NOT_FC_or_topN")
      dir.create(full_gsea_dir, recursive = TRUE, showWarnings = FALSE)
      writeLines("全基因 GSEA，不是 FC 分层结果。分层图在 FoldChange/FC_1.25 和 FC_1.5。",
                 file.path(full_gsea_dir, "00_README.txt"))
      for (nm in c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Hallmark")) {
        plot_gsea_object(gsea_cache[[nm]], file.path(full_gsea_dir, paste0("allGenes_GSEA_", nm)),
                         paste("GSEA", nm, "|", comp_name, "| ALL genes, NOT FC subset"))
      }
      plot_fgsea_hallmark(gsea_cache$stats, full_gsea_dir,
                          paste("GSEA Hallmark |", comp_name, "| ALL genes"), prefix = "allGenes_")
    }, error = function(e) log_msg("full-list GSEA failed for ", comp_name, ": ", e$message))
    tryCatch({
      focus_stats <- if (!is.null(gsea_cache$stats)) gsea_cache$stats else ranked_entrez(de)
      run_focused_gsea(
        focus_stats, de, heat_mat, sample_info,
        file.path(base, "Focused_cytoskeleton_mito"),
        paste(comp_name, "| all genes")
      )
    }, error = function(e) log_msg("focused GSEA failed for ", comp_name, ": ", e$message))
  }
  invisible(TRUE)
}

build_heat_from_pb <- function(pb, si) {
  logm <- log2(sweep(pb, 2, pmax(colSums(pb), 1), "/") * 1e6 + 1)
  list(mat = logm, sample_info = si)
}

# -----------------------------------------------------------------------------
# 7. 主分析
# -----------------------------------------------------------------------------
ids <- load_ids()
pids <- list_patient_ids()
if (length(pids) == 0) stop("Robjects/counts 下没有 TNBC*.RDS。请先解压 Robjects.tar")
log_msg("Found count files for ", length(pids), " patients")

max_n <- suppressWarnings(as.integer(Sys.getenv("WANG_ST_MAX_PATIENTS", unset = "")))
if (!is.na(max_n) && max_n > 0) {
  pids <- utils::head(pids, max_n)
  log_msg("WANG_ST_MAX_PATIENTS=", max_n)
}

qc_rows <- list()
nerve_de_list <- list()
pb_schwann <- NULL
si_schwann <- NULL
ligand_rows <- list()

for (pid in pids) {
  log_msg("Patient TNBC", pid)
  obj <- tryCatch(classify_one_patient(pid, ids), error = function(e) {
    log_msg("classify failed TNBC", pid, ": ", e$message)
    NULL
  })
  if (is.null(obj)) next
  sp <- obj$spots
  qc_rows[[as.character(pid)]] <- data.frame(
    patient = pid,
    n_spots = nrow(sp),
    n_tumor = sum(sp$is_tumor),
    n_nerve_path = sum(sp$is_nerve_path),
    n_schwann_high = sum(sp$is_schwann_high),
    n_tumor_near_nerve = sum(sp$near_nerve),
    n_tumor_far_nerve = sum(sp$far_nerve),
    n_tumor_near_schwann = sum(sp$near_schwann),
    n_tumor_far_schwann = sum(sp$far_schwann),
    annot_slide = ifelse(is.na(obj$annot_slide), "", obj$annot_slide),
    stringsAsFactors = FALSE
  )
  tryCatch(
    spot_map_plot(sp, paste("TNBC", pid, "spot classes"),
                  file.path(result_dir, "00_QC_maps", paste0("TNBC", pid, "_spot_map"))),
    error = function(e) log_msg("map plot failed: ", e$message)
  )

  gkeep <- filter_genes_mat(obj$cnts)
  logm <- obj$logmat[, gkeep, drop = FALSE]
  cnts <- obj$cnts[, gkeep, drop = FALSE]

  # 病理神经 1-vs-1（不合并病人）
  if (sum(sp$near_nerve) >= 5 && sum(sp$far_nerve) >= 10) {
    grp <- ifelse(sp$near_nerve, "near", ifelse(sp$far_nerve, "far", NA))
    keep <- !is.na(grp)
    de_i <- limma_two_group(logm[keep, , drop = FALSE], grp[keep],
                            label = paste0("TNBC", pid, "_nerve"))
    if (nrow(de_i) == 0) {
      de_i <- mean_fc_two_group(logm, which(sp$near_nerve), which(sp$far_nerve))
      log_msg("TNBC", pid, " nerve: limma empty, fallback mean FC (no p)")
    }
    heat <- t(logm[keep, , drop = FALSE])
    si <- data.frame(sample = rownames(logm)[keep], group = grp[keep],
                     stringsAsFactors = FALSE)
    # heatmap 列太多时改用伪 bulk
    pb <- cbind(
      near = colSums(cnts[sp$near_nerve, , drop = FALSE]),
      far  = colSums(cnts[sp$far_nerve, , drop = FALSE])
    )
    colnames(pb) <- paste0("TNBC", pid, "_", c("near", "far"))
    heat_pb <- log2(sweep(pb, 2, pmax(colSums(pb), 1), "/") * 1e6 + 1)
    si_pb <- data.frame(sample = colnames(heat_pb),
                        group = c("near", "far"), stringsAsFactors = FALSE)
    nm <- paste0("TNBC", pid, "_tumor_near_nerve_vs_far")
    nerve_de_list[[nm]] <- de_i
    emit_comparison(nm, de_i, heat_pb, si_pb)
  }

  # Schwann 邻域伪 bulk（全队列，病人作为协变量，不把对照混成一组去平均）
  if (sum(sp$near_schwann) >= 5 && sum(sp$far_schwann) >= 10) {
    near_sum <- colSums(cnts[sp$near_schwann, , drop = FALSE])
    far_sum  <- colSums(cnts[sp$far_schwann, , drop = FALSE])
    pb <- cbind(near_sum, far_sum)
    colnames(pb) <- paste0("TNBC", pid, c("_near", "_far"))
    si <- data.frame(
      sample = colnames(pb),
      group = c("near", "far"),
      patient = as.character(pid),
      stringsAsFactors = FALSE
    )
    if (is.null(pb_schwann)) {
      pb_schwann <- pb
      si_schwann <- si
    } else {
      gn <- intersect(rownames(pb_schwann), rownames(pb))
      pb_schwann <- cbind(pb_schwann[gn, , drop = FALSE], pb[gn, , drop = FALSE])
      si_schwann <- rbind(si_schwann, si)
    }
  }

  # 配体在肿瘤 spot 的表达 vs 邻域 Schwann
  lig <- present_genes(tumor_ligands$gene, colnames(logm))
  tum <- which(sp$is_tumor)
  if (length(lig) > 0 && length(tum) >= 20) {
    cors <- sapply(lig, function(g) {
      stats::cor(logm[tum, g], -pmin(sp$dist_schwann[tum], 10),
                 method = "spearman", use = "complete.obs")
    })
    ligand_rows[[as.character(pid)]] <- data.frame(
      patient = pid,
      gene = lig,
      mean_in_tumor = colMeans(logm[tum, lig, drop = FALSE]),
      spearman_vs_near_schwann = as.numeric(cors),
      stringsAsFactors = FALSE
    )
  }
  rm(obj, cnts, logm)
  gc(verbose = FALSE)
}

qc <- do.call(rbind, qc_rows)
utils::write.csv(qc, file.path(log_dir, "spot_class_counts.csv"), row.names = FALSE)
log_msg("Patients with pathologist nerve spots: ",
        paste(qc$patient[qc$n_nerve_path > 0], collapse = ", "))

# 共同上调：病理神经 1-vs-1 交集
if (length(nerve_de_list) >= 2) {
  have_p <- all(vapply(nerve_de_list, function(d) any(!is.na(d$pvalue)), logical(1)))
  get_up <- function(de, fc) {
    keep <- !is.na(de$log2FC) & de$log2FC > 0 & (2^de$log2FC >= fc)
    if (have_p) keep <- keep & !is.na(de$pvalue) & de$pvalue < p_cutoff
    de$gene[keep]
  }
  for (nm in names(fc_cutoffs_local)) {
    sets <- lapply(nerve_de_list, get_up, fc = unname(fc_cutoffs_local[[nm]]))
    inter <- Reduce(intersect, sets)
    log_msg("common up pathologist nerve ", nm, " n=", length(inter))
  }
  genes_all <- Reduce(intersect, lapply(nerve_de_list, function(d) {
    keep <- d$log2FC > 0
    if (have_p) keep <- keep & !is.na(d$pvalue) & d$pvalue < p_cutoff
    d$gene[keep]
  }))
  if (length(genes_all) > 0) {
    tab <- data.frame(gene = genes_all, stringsAsFactors = FALSE)
    for (nm in names(nerve_de_list)) {
      d <- nerve_de_list[[nm]]
      tab[[paste0("log2FC_", nm)]] <- d$log2FC[match(tab$gene, d$gene)]
      tab[[paste0("pvalue_", nm)]] <- d$pvalue[match(tab$gene, d$gene)]
    }
    tab$log2FC <- rowMeans(as.matrix(tab[, grep("^log2FC_", names(tab)), drop = FALSE]), na.rm = TRUE)
    tab$pvalue <- apply(as.matrix(tab[, grep("^pvalue_", names(tab)), drop = FALSE]), 1, function(x) {
      if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    })
    tab$AveExpr <- NA_real_
    tab$padj <- tab$pvalue
    tab$padj_BH <- NA_real_
    heat <- matrix(tab$log2FC, ncol = 1, dimnames = list(tab$gene, "mean_log2FC"))
    si <- data.frame(sample = "mean_log2FC", group = "near", stringsAsFactors = FALSE)
    emit_comparison("common_up_pathologist_nerve", tab, heat, si)
  } else {
    log_msg("common_up_pathologist_nerve: empty intersection")
  }
} else {
  log_msg("病理 Nerve 重叠 spot 的病人不足 2 例，跳过共同上调。请看全队列 Schwann 邻域。")
}

# 全队列 Schwann 邻域
if (!is.null(pb_schwann) && length(unique(si_schwann$patient)) >= 2) {
  log_msg("Schwann-neighborhood pseudobulk patients: ",
          paste(unique(si_schwann$patient), collapse = ","))
  rownames(si_schwann) <- si_schwann$sample
  de_sw <- deseq2_pseudobulk(round(pmax(pb_schwann, 0)), si_schwann,
                             "tumor_near_schwann_vs_far")
  if (nrow(de_sw) == 0) {
    logm <- log2(sweep(t(pb_schwann), 1, pmax(colSums(pb_schwann), 1), "/") * 1e6 + 1)
    de_sw <- limma_two_group(logm, si_schwann$group, si_schwann$patient, "schwann_pb")
  }
  heat_obj <- build_heat_from_pb(pb_schwann, si_schwann)
  emit_comparison("tumor_near_schwann_vs_far", de_sw, heat_obj$mat, heat_obj$sample_info)
} else {
  log_msg("Schwann 邻域伪 bulk 病人不足，跳过 tumor_near_schwann_vs_far")
  de_sw <- empty_de()
}

# -----------------------------------------------------------------------------
# 8. 候选分子总表
# -----------------------------------------------------------------------------
lig_all <- if (length(ligand_rows) > 0) do.call(rbind, ligand_rows) else NULL
if (!is.null(lig_all)) {
  agg <- do.call(rbind, lapply(split(lig_all, lig_all$gene), function(d) {
    data.frame(
      gene = d$gene[1],
      n_patients = nrow(d),
      mean_tumor_expr = mean(d$mean_in_tumor, na.rm = TRUE),
      mean_spearman_near_schwann = mean(d$spearman_vs_near_schwann, na.rm = TRUE),
      n_pos_cor = sum(d$spearman_vs_near_schwann > 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
} else {
  agg <- data.frame(gene = tumor_ligands$gene, n_patients = 0,
                    mean_tumor_expr = NA, mean_spearman_near_schwann = NA, n_pos_cor = 0)
}

universe_example <- character()
if (length(pids) > 0) {
  ex <- load_counts(pids[1])
  if (!is.null(ex)) universe_example <- colnames(ex$cnts)
}

cand <- tumor_ligands
cand$detected_on_array <- cand$gene %in% universe_example
cand <- merge(cand, agg, by = "gene", all.x = TRUE)
if (nrow(de_sw) > 0) {
  cand$schwann_near_log2FC <- de_sw$log2FC[match(cand$gene, de_sw$gene)]
  cand$schwann_near_pvalue <- de_sw$pvalue[match(cand$gene, de_sw$gene)]
} else {
  cand$schwann_near_log2FC <- NA_real_
  cand$schwann_near_pvalue <- NA_real_
}
if (length(nerve_de_list) > 0) {
  for (i in seq_along(nerve_de_list)) {
    d <- nerve_de_list[[i]]
    cand[[paste0("pathologist_", names(nerve_de_list)[i], "_log2FC")]] <-
      d$log2FC[match(cand$gene, d$gene)]
  }
}
cand$up_in_schwann_near <- !is.na(cand$schwann_near_pvalue) &
  cand$schwann_near_pvalue < p_cutoff & cand$schwann_near_log2FC > 0
cand$priority <- (cand$detected_on_array) +
  2 * (!is.na(cand$mean_spearman_near_schwann) & cand$mean_spearman_near_schwann > 0.1) +
  3 * cand$up_in_schwann_near
cand <- cand[order(-cand$priority, cand$schwann_near_pvalue, -cand$mean_spearman_near_schwann), ]
utils::write.csv(cand, file.path(result_dir, "01_CANDIDATE_MOLECULES_tumor_to_nerve.csv"),
                 row.names = FALSE)
tryCatch(writexl::write_xlsx(cand, file.path(result_dir, "01_CANDIDATE_MOLECULES_tumor_to_nerve.xlsx")),
         error = function(e) NULL)

# TG knockdown 上调基因重叠（若已跑过原流程）
tg_hits <- c(
  "E:/R/TG_BRCA/TG/results",
  file.path(getwd(), "results")
)
tg_files <- unlist(lapply(tg_hits, function(d) {
  list.files(d, pattern = "^DE_full.csv$", recursive = TRUE, full.names = TRUE)
}))
tg_files <- tg_files[!grepl("Wang|nerve|schwann|pathologist", tg_files, ignore.case = TRUE)]
if (length(tg_files) > 0 && nrow(de_sw) > 0) {
  log_msg("Overlap with TG DE tables: ", length(tg_files), " files")
  ov_dir <- file.path(result_dir, "overlap_TG_knockdown")
  dir.create(ov_dir, recursive = TRUE, showWarnings = FALSE)
  sw_up <- de_sw$gene[!is.na(de_sw$pvalue) & de_sw$pvalue < p_cutoff & de_sw$log2FC > 0]
  for (f in tg_files) {
    tg <- tryCatch(utils::read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(tg) || !"gene" %in% names(tg)) next
    keep <- tg$log2FC > 0
    if ("pvalue" %in% names(tg) && any(!is.na(tg$pvalue))) {
      keep <- keep & !is.na(tg$pvalue) & tg$pvalue < p_cutoff
    }
    inter <- intersect(unique(tg$gene[keep]), sw_up)
    nm <- basename(dirname(f))
    utils::write.csv(data.frame(gene = inter),
                     file.path(ov_dir, paste0("overlap_", nm, ".csv")), row.names = FALSE)
    log_msg("overlap ", nm, " n=", length(inter))
  }
}

# -----------------------------------------------------------------------------
# 9. 方案说明（问题 1）
# -----------------------------------------------------------------------------
protocol <- c(
  "============================================================",
  "乳腺癌空间转录组：神经浸润分析方案（Wang et al. Nat Commun 2024）",
  "数据目录: ", nerve_dir,
  "============================================================",
  "",
  "【问题 1】明确神经浸润相关分子和方案",
  "",
  "A. 为什么不能指望『神经元 cluster』",
  "  瘤内主要是神经纤维 + Schwann 细胞，神经元胞体在神经节。",
  "  本队列 94 例里，病理 Nerve 与 ST spot 重叠的只有极少数（见 00_logs/spot_class_counts.csv）。",
  "  因此方案是双轨：",
  "    轨 1  病理 Nerve 邻域（金标准，病人少，单独 1-vs-1，不合并）",
  "    轨 2  Schwann 签名邻域（全队列，病人作为协变量）",
  "",
  "B. Spot 定义（写进日志，不事后改阈值去凑显著）",
  "  神经（病理）：Nerve 像素 >= 1",
  "  肿瘤：Tumor+Tumor region 像素比例 >= 0.25；无标注则用 EPCAM/KRT 签名",
  "  近：阵列坐标欧氏距离 <= 2（约两个 spot）",
  "  远：距离 > 4（中间空一圈，避免模糊带）",
  "  Schwann 高：切片内 z-score > 1；基因 SOX10 MPZ PMP22 S100B PLP1 NGFR NCAM1 MBP L1CAM",
  "",
  "C. 统计（本空间分析，不是 Cuffdiff 那六组）",
  "  先 p < 0.05，再只看上调 FC >= 1.25 和 FC >= 1.5",
  "  不做 FC=1、FC=2，也不做 top 50–300",
  "  只看上调（近神经肿瘤 > 远神经肿瘤）",
  "  有病理神经的病人各自 1-vs-1；共同上调 = 这些 1-vs-1 的交集",
  "  全队列 Schwann：伪 bulk ~ patient + group，不用把病人平均掉再比",
  "  低表达过滤后再算 FC；DESeq2 size factor / limma-voom，不用原始 count 直接除",
  "",
  "D. 文献里已经比较清楚的『肿瘤→神经』分子（本阵列未必都测到）",
  "  NGF–NTRK1/NGFR，BDNF–NTRK2，ARTN–GFRA3/RET，GDNF 家族",
  "  NTN1、SLIT2、SEMA3C/F、NRG1、CXCL12–CXCR4、MDK/PTN、VEGFA",
  "  TNBC 感觉神经：NGF → CGRP(CALCA) → CAF RAMP1（Cell 2026）",
  "  详细候选表: 01_CANDIDATE_MOLECULES_tumor_to_nerve.csv",
  "",
  "E. 实验上怎么验证（分析之后）",
  "  1. 候选配体在肿瘤细胞里敲低/过表达，背根神经节共培养看轴突长度",
  "  2. 阻断受体：NTRK 抑制剂、NGFR、GFRA3、CXCR4、RAMP1（rimegepant）",
  "  3. 切片 IHC：候选配体 vs PGP9.5/S100/SOX10 空间相邻",
  "  4. 不要把免疫细胞高表达基因当成肿瘤细胞促浸润分子",
  "",
  "【问题 2】肿瘤细胞高表达哪些基因可以促使神经浸润",
  "",
  "  主结果文件夹（按这个顺序看）:",
  "    1. 01_CANDIDATE_MOLECULES_tumor_to_nerve.csv",
  "       priority 高 = 阵列检出 + 靠近 Schwann 的肿瘤里升高 + 与距离负相关",
  "    2. results/tumor_near_schwann_vs_far/",
  "       全队列：靠近 Schwann 的肿瘤 vs 远离 Schwann 的肿瘤",
  "       FoldChange/FC_1.25 和 FC_1.5 才是分层图；没有 TopRank",
  "       全基因 GSEA 在 00_GSEA_all_genes_NOT_FC_or_topN",
  "    3. results/TNBC*_tumor_near_nerve_vs_far/",
  "       病理神经金标准，病人很少，只作验证，不要和 Schwann 结果混成一张表",
  "    4. results/common_up_pathologist_nerve/  上述金标准的共同上调",
  "",
  "  解释时注意：邻域基因升高 = 与神经空间相邻，不等于已证明『引起』浸润。",
  "  促浸润要用配体功能（轴突生长、Schwann 趋化）再验证。",
  "",
  "ORA 在 GO/ Pathway/ KEGG/，文件名 ORA_ 开头。",
  "GSEA 在 GSEA/，文件名 GSEA_ 开头。",
  "细胞骨架/线粒体专项在 Focused_cytoskeleton_mito/，不改全库 p 值。"
)
writeLines(protocol, file.path(result_dir, "00_PROTOCOL_神经浸润分析方案.txt"))
log_msg("Wrote protocol and candidate molecule table")
log_msg("Done. Open: ", file.path(result_dir, "00_PROTOCOL_神经浸润分析方案.txt"))
log_msg("and: ", file.path(result_dir, "01_CANDIDATE_MOLECULES_tumor_to_nerve.csv"))
