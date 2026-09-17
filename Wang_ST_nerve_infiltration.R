#!/usr/bin/env Rscript
# =============================================================================
# Wang et al. Nature Communications 2024  TNBC 空间转录组
# 神经浸润：三个问题（不做 GO / KEGG / GSEA / ORA）
#
# 只回答三件事：
#   问题1  神经浸润相关分子 + 统计方案
#   问题2  肿瘤细胞高表达、可能促使神经浸润的基因
#   问题3  肿瘤细胞低表达、可能促使神经浸润的基因（丢失排斥/屏障）
# 独立于 TG_RNAseq_pipeline.R，不改原 Cuffdiff 1–4。
# 数据：E:/R/Nerve 的 Robjects/（counts + annotsBySpot）。
# ids.RDS 可选（Windows 常把 Clinical/ids.RDS 展成 Clinicalids.RDS）。
#
# 用法（Windows R / RStudio）：
#   setwd("E:/R/Nerve")
#   source("E:/R/TG_BRCA/TG/Wang_ST_nerve_infiltration.R")   # 或把本文件拷到 E:/R/Nerve
#   # 也可：Sys.setenv(WANG_ST_DIR = "E:/R/Nerve")
#
# 比较（各自单独出表/出图；病人不偷偷合并）：
#   1) 病理 Nerve 邻域：每个有神经 spot 的病人单独 1-vs-1
#   2) 这些病人的共同上调 / 共同下调
#   3) 全队列 Schwann 签名邻域：~ patient + group 的伪 bulk DESeq2
# 显著性：先 p < 0.05，再按 |FC| >= 1.25 / 1.5 分上调（问题2）和下调（问题3）。
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
        file.exists(file.path(d, "Clinical", "ids.RDS")) ||
        file.exists(file.path(d, "ids.RDS")) ||
        file.exists(file.path(d, "Clinicalids.RDS")) ||
        file.exists(file.path(d, "Clinical.RDS")) ||
        file.exists(file.path(d, "Clinical.xlsx"))
    )
  }
  for (d in candidates) {
    if (looks_like_wang(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  stop("找不到 Wang ST 数据。请把 Zenodo 文件放到 E:/R/Nerve，或设置 WANG_ST_DIR。")
}

nerve_dir <- resolve_nerve_dir()

find_under <- function(root, filename) {
  hits <- list.files(root, pattern = paste0("^", filename, "$"),
                     recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  hits <- hits[!grepl("[/\\\\]__MACOSX[/\\\\]", hits)]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

ensure_tar_extracted <- function(root, dirname, tarname) {
  dest <- file.path(root, dirname)
  if (dir.exists(dest) && length(list.files(dest, recursive = TRUE)) > 0) return(invisible(dest))
  tars <- c(
    file.path(root, tarname),
    list.files(root, pattern = paste0("^", tarname, "$"), recursive = TRUE, full.names = TRUE)
  )
  tars <- unique(tars[file.exists(tars)])
  if (length(tars) == 0) return(invisible(NULL))
  message("解压 ", basename(tars[1]), " -> ", root)
  tryCatch(utils::untar(tars[1], exdir = root, tar = "internal"),
           error = function(e) utils::untar(tars[1], exdir = root))
  dest
}

ensure_tar_extracted(nerve_dir, "Clinical", "Clinical.tar")
ensure_tar_extracted(nerve_dir, "Robjects", "Robjects.tar")
ensure_tar_extracted(nerve_dir, "misc", "misc.tar")
ensure_tar_extracted(nerve_dir, "patches", "patches.tar")
ensure_tar_extracted(nerve_dir, "classification", "classification.tar")

ids_path <- find_under(nerve_dir, "ids.RDS")
if (is.na(ids_path)) {
  flat <- list.files(nerve_dir, pattern = "ids\\.RDS$", full.names = TRUE, ignore.case = TRUE)
  if (length(flat) > 0) ids_path <- flat[1]
}
if (!is.na(ids_path) && nzchar(ids_path)) {
  clinical_dir <- dirname(ids_path)
  if (basename(clinical_dir) == "Clinical") {
    nerve_dir <- dirname(clinical_dir)
  }
}

robj_dir <- if (dir.exists(file.path(nerve_dir, "Robjects"))) {
  file.path(nerve_dir, "Robjects")
} else if (dir.exists(file.path(nerve_dir, "annotsBySpot"))) {
  nerve_dir
} else {
  hit <- find_under(nerve_dir, "TNBC50.RDS")
  if (!is.na(hit) && basename(dirname(hit)) %in% c("counts", "countsNonCorrected", "annotsBySpot")) {
    dirname(dirname(hit))
  } else {
    nerve_dir
  }
}

# -----------------------------------------------------------------------------
# 1. 包（只要差异分析；不要 GO/GSEA）
# -----------------------------------------------------------------------------
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

if (!exists("log_msg", mode = "function")) {
  log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
}

result_dir <- file.path(nerve_dir, "results")
log_dir <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result_dir, "00_QC_maps"), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("wang_st_nerve_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg_file <- log_msg
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

# 先 p < 0.05，再 |FC| >= 1.25 / 1.5（上调=问题2，下调=问题3；不做 GO、不做 topN）
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

# 肿瘤细胞低表达后可能放开神经浸润：轴突排斥、基底膜、上皮屏障
tumor_repellents <- data.frame(
  gene = c(
    "SEMA3A", "SEMA3B", "SEMA3F", "SLIT2", "SLIT3",
    "EFNA1", "EFNA5", "EPHA2", "DCC", "UNC5B", "PLXNA1",
    "CDH1", "TJP1", "OCLN", "CLDN1", "COL4A1", "LAMA2", "LAMB1"
  ),
  receptor_or_partner = c(
    "NRP1/PLXNA", "NRP", "NRP2", "ROBO1/2", "ROBO",
    "EPHA", "EPHA", "EFNA", "NTN1", "NTN1", "SEMA3",
    "CTNNB1", "TJP", "CLDN", "OCLN", "ITG", "ITG", "ITG"
  ),
  rationale = c(
    "3 类 semaphorin，轴突排斥；肿瘤里丢失利于神经长入",
    "SEMA3B 排斥 / 抑癌",
    "SEMA3F–NRP2 排斥",
    "SLIT2–ROBO 排斥；多种癌中丢失与 PNI 有关",
    "SLIT3 排斥",
    "ephrin-A 接触排斥（下调才算问题3）",
    "ephrin-A5 排斥",
    "EphA2 接触导向",
    "DCC：netrin 依赖的排斥/吸引取决于受体组合",
    "UNC5B：netrin 排斥受体",
    "Plexin-A1，semaphorin 排斥",
    "E-cadherin 上皮屏障；丢失 / EMT 利于沿神经侵袭",
    "紧密连接",
    "occludin 紧密连接",
    "claudin-1 紧密连接",
    "IV 型胶原，基底膜",
    "laminin-α2，基底膜 / 神经周围基质",
    "laminin-β1，基底膜"
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

is_ids_map <- function(obj) {
  if (is.null(obj)) return(FALSE)
  if (is.matrix(obj)) obj <- as.data.frame(obj, stringsAsFactors = FALSE)
  if (!is.data.frame(obj)) return(FALSE)
  all(c("id", "hasAnnot") %in% names(obj))
}

coerce_ids_map <- function(obj) {
  if (is.matrix(obj)) obj <- as.data.frame(obj, stringsAsFactors = FALSE)
  if (is.data.frame(obj) && is_ids_map(obj)) return(obj)
  if (is.list(obj) && !is.data.frame(obj)) {
    for (el in obj) {
      got <- coerce_ids_map(el)
      if (!is.null(got)) return(got)
    }
  }
  NULL
}

list_rds_shallow <- function(d) {
  if (!dir.exists(d)) return(character())
  list.files(d, pattern = "\\.RDS$", full.names = TRUE, ignore.case = TRUE)
}

load_ids <- function() {
  # Windows 解 Clinical.tar 时常丢掉子目录：Clinical/ids.RDS → Clinicalids.RDS
  # Clinical.RDS 是 94×50 临床表，不是片子对照，不能当 ids 用。
  cand <- unique(c(
    file.path(nerve_dir, "Clinical", "ids.RDS"),
    file.path(nerve_dir, "ids.RDS"),
    file.path(nerve_dir, "Clinicalids.RDS"),
    file.path(nerve_dir, "Clinical_ids.RDS"),
    list.files(nerve_dir, pattern = "ids\\.RDS$", full.names = TRUE, ignore.case = TRUE),
    list_rds_shallow(file.path(nerve_dir, "Clinical")),
    list_rds_shallow(nerve_dir)
  ))
  cand <- cand[file.exists(cand)]
  cand <- cand[!grepl("Robjects|counts|annotsBySpot|images|BatchCorrection|__MACOSX",
                      cand, ignore.case = TRUE)]
  for (f in cand) {
    obj <- coerce_ids_map(safe_read_rds(f))
    if (is.null(obj)) next
    log_msg("片子↔病人对照: ", f,
            if (grepl("Clinicalids\\.RDS$", basename(f), ignore.case = TRUE))
              "（Windows 把 Clinical/ids.RDS 展成了这个文件名）" else "")
    return(obj)
  }
  shown <- paste(list.files(nerve_dir), collapse = ", ")
  log_msg("没有 ids.RDS（片子↔病人对照）。根目录若只有 Clinical.RDS / Clinical.xlsx，那是临床表，不是对照。")
  log_msg("将用 counts 坐标去对病理标注，分析继续。当前文件: ", shown)
  NULL
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
    key_sp <- paste(spots$x, spots$y, sep = "x")
    fr <- spot_ann_frac(aa$annots)
    if (is.na(annot_slide) && "slide_id" %in% names(spots)) {
      ov <- tapply(seq_len(nrow(spots)), spots$slide_id, function(ii) {
        length(intersect(key_sp[ii], fr$key))
      })
      if (length(ov) > 0 && max(as.numeric(ov), na.rm = TRUE) > 0) {
        annot_slide <- names(ov)[which.max(ov)]
        log_msg("TNBC", pid, " 病理标注对齐片子: ", annot_slide)
      }
    }
    on_annot <- if (!is.na(annot_slide) && "slide_id" %in% names(spots)) {
      spots$slide_id == annot_slide
    } else {
      key_sp %in% fr$key
    }
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

select_up <- function(de, fc = 1) {
  keep <- !is.na(de$log2FC) & de$log2FC > 0 & (2^de$log2FC >= fc)
  if (any(!is.na(de$pvalue))) {
    keep <- keep & !is.na(de$pvalue) & de$pvalue < p_cutoff
  }
  sub <- de[keep, , drop = FALSE]
  if (nrow(sub) == 0) return(sub)
  sub$FC <- 2^sub$log2FC
  if (any(!is.na(sub$pvalue))) {
    sub[order(sub$pvalue, -sub$log2FC), , drop = FALSE]
  } else {
    sub[order(-sub$log2FC), , drop = FALSE]
  }
}

select_down <- function(de, fc = 1) {
  # 近/远 FC <= 1/fc，即下调倍数 >= fc（1.25 或 1.5）
  keep <- !is.na(de$log2FC) & de$log2FC < 0 & (2^de$log2FC <= 1 / fc)
  if (any(!is.na(de$pvalue))) {
    keep <- keep & !is.na(de$pvalue) & de$pvalue < p_cutoff
  }
  sub <- de[keep, , drop = FALSE]
  if (nrow(sub) == 0) return(sub)
  sub$FC <- 2^sub$log2FC
  sub$down_FC <- 1 / pmax(sub$FC, 1e-12)
  if (any(!is.na(sub$pvalue))) {
    sub[order(sub$pvalue, sub$log2FC), , drop = FALSE]
  } else {
    sub[order(sub$log2FC), , drop = FALSE]
  }
}

# -----------------------------------------------------------------------------
# 6. 差异表（只出 CSV + 火山图，不出 GO/GSEA）
# -----------------------------------------------------------------------------
basic_volcano <- function(de, title, outfile, fc_line = 1) {
  if (nrow(de) == 0) return(invisible(NULL))
  df <- de
  df$y <- if (any(!is.na(df$pvalue))) -log10(pmax(df$pvalue, 1e-300)) else abs(df$log2FC)
  df$col <- "ns"
  sig <- is.na(df$pvalue) | df$pvalue < p_cutoff
  df$col[sig & !is.na(df$log2FC) & df$log2FC >= log2(fc_line)] <- "up"
  df$col[sig & !is.na(df$log2FC) & df$log2FC <= -log2(fc_line)] <- "down"
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(log2FC, y, color = col)) +
      ggplot2::geom_point(alpha = 0.5, size = 0.8) +
      ggplot2::scale_color_manual(values = c(ns = "grey70", up = "#D62828", down = "#1D4E89")) +
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

emit_de_tables <- function(comp_name, de) {
  de <- de[!is.na(de$log2FC), , drop = FALSE]
  base <- file.path(result_dir, comp_name)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(de, file.path(base, "DE_full.csv"), row.names = FALSE)
  have_p <- any(!is.na(de$pvalue))
  n_up <- sum((is.na(de$pvalue) | de$pvalue < p_cutoff) & de$log2FC > 0, na.rm = TRUE)
  n_dn <- sum((is.na(de$pvalue) | de$pvalue < p_cutoff) & de$log2FC < 0, na.rm = TRUE)
  log_msg(comp_name, " genes=", nrow(de), " have_p=", have_p,
          " up p<", p_cutoff, " n=", n_up, " down n=", n_dn)
  for (nm in names(fc_cutoffs_local)) {
    fc <- unname(fc_cutoffs_local[[nm]])
    up <- select_up(de, fc)
    dn <- select_down(de, fc)
    od_up <- file.path(base, paste0("up_", nm))
    od_dn <- file.path(base, paste0("down_", nm))
    dir.create(od_up, recursive = TRUE, showWarnings = FALSE)
    dir.create(od_dn, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(up, file.path(od_up, paste0(nm, "_upregulated_tumor_genes.csv")),
                     row.names = FALSE)
    utils::write.csv(dn, file.path(od_dn, paste0(nm, "_downregulated_tumor_genes.csv")),
                     row.names = FALSE)
    basic_volcano(de, paste(comp_name, "| |FC| >=", fc),
                  file.path(base, paste0("volcano_", nm)), fc)
    log_msg("  ", nm, " up n=", nrow(up), " down n=", nrow(dn))
  }
  invisible(de)
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
    nm <- paste0("TNBC", pid, "_tumor_near_nerve_vs_far")
    nerve_de_list[[nm]] <- de_i
    emit_de_tables(nm, de_i)
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

# 共同上调 / 共同下调：病理神经 1-vs-1 交集
if (length(nerve_de_list) >= 2) {
  have_p <- all(vapply(nerve_de_list, function(d) any(!is.na(d$pvalue)), logical(1)))
  get_up <- function(de, fc) {
    keep <- !is.na(de$log2FC) & de$log2FC > 0 & (2^de$log2FC >= fc)
    if (have_p) keep <- keep & !is.na(de$pvalue) & de$pvalue < p_cutoff
    de$gene[keep]
  }
  get_down <- function(de, fc) {
    keep <- !is.na(de$log2FC) & de$log2FC < 0 & (2^de$log2FC <= 1 / fc)
    if (have_p) keep <- keep & !is.na(de$pvalue) & de$pvalue < p_cutoff
    de$gene[keep]
  }
  for (nm in names(fc_cutoffs_local)) {
    fc <- unname(fc_cutoffs_local[[nm]])
    log_msg("common up pathologist nerve ", nm, " n=",
            length(Reduce(intersect, lapply(nerve_de_list, get_up, fc = fc))))
    log_msg("common down pathologist nerve ", nm, " n=",
            length(Reduce(intersect, lapply(nerve_de_list, get_down, fc = fc))))
  }
  collapse_inter <- function(genes_all) {
    if (length(genes_all) == 0) return(NULL)
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
    tab
  }
  genes_up <- Reduce(intersect, lapply(nerve_de_list, function(d) {
    keep <- d$log2FC > 0
    if (have_p) keep <- keep & !is.na(d$pvalue) & d$pvalue < p_cutoff
    d$gene[keep]
  }))
  genes_dn <- Reduce(intersect, lapply(nerve_de_list, function(d) {
    keep <- d$log2FC < 0
    if (have_p) keep <- keep & !is.na(d$pvalue) & d$pvalue < p_cutoff
    d$gene[keep]
  }))
  tab_up <- collapse_inter(genes_up)
  tab_dn <- collapse_inter(genes_dn)
  if (!is.null(tab_up)) emit_de_tables("common_up_pathologist_nerve", tab_up) else
    log_msg("common_up_pathologist_nerve: empty intersection")
  if (!is.null(tab_dn)) emit_de_tables("common_down_pathologist_nerve", tab_dn) else
    log_msg("common_down_pathologist_nerve: empty intersection")
} else {
  log_msg("病理 Nerve 重叠 spot 的病人不足 2 例，跳过共同上调/下调。请看全队列 Schwann 邻域。")
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
  emit_de_tables("tumor_near_schwann_vs_far", de_sw)
} else {
  log_msg("Schwann 邻域伪 bulk 病人不足，跳过 tumor_near_schwann_vs_far")
  de_sw <- empty_de()
}

# -----------------------------------------------------------------------------
# 8. 问题 1：文献配体；问题 2：肿瘤细胞上调基因总表（主结果）
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
cand$up_in_schwann_near <- (!is.na(cand$schwann_near_log2FC) & cand$schwann_near_log2FC > 0) &
  (is.na(cand$schwann_near_pvalue) | cand$schwann_near_pvalue < p_cutoff)
cand$priority <- as.integer(cand$detected_on_array) +
  2L * as.integer(!is.na(cand$mean_spearman_near_schwann) & cand$mean_spearman_near_schwann > 0.1) +
  3L * as.integer(cand$up_in_schwann_near)
cand <- cand[order(-cand$priority, cand$schwann_near_pvalue, -cand$mean_spearman_near_schwann), ]
utils::write.csv(cand, file.path(result_dir, "01_Q1_literature_ligands.csv"), row.names = FALSE)

repel <- tumor_repellents
repel$detected_on_array <- repel$gene %in% universe_example
if (nrow(de_sw) > 0) {
  repel$schwann_near_log2FC <- de_sw$log2FC[match(repel$gene, de_sw$gene)]
  repel$schwann_near_pvalue <- de_sw$pvalue[match(repel$gene, de_sw$gene)]
} else {
  repel$schwann_near_log2FC <- NA_real_
  repel$schwann_near_pvalue <- NA_real_
}
repel$down_in_schwann_near <- (!is.na(repel$schwann_near_log2FC) & repel$schwann_near_log2FC < 0) &
  (is.na(repel$schwann_near_pvalue) | repel$schwann_near_pvalue < p_cutoff)
repel$priority <- as.integer(repel$detected_on_array) +
  3L * as.integer(repel$down_in_schwann_near)
repel <- repel[order(-repel$priority, repel$schwann_near_pvalue, repel$schwann_near_log2FC), ]
utils::write.csv(repel, file.path(result_dir, "01_Q3_literature_repellents.csv"), row.names = FALSE)

# 问题 2 主表：靠近神经/Schwann 的肿瘤细胞里升高的基因，按 p 再按 FC 排
q2_src <- if (nrow(de_sw) > 0) {
  de_sw
} else if (length(nerve_de_list) > 0) {
  log_msg("全队列 Schwann 表为空，问题2改用病理神经 1-vs-1 合并排名")
  Reduce(function(a, b) {
    g <- unique(c(a$gene, b$gene))
    data.frame(
      gene = g,
      log2FC = rowMeans(cbind(a$log2FC[match(g, a$gene)], b$log2FC[match(g, b$gene)]), na.rm = TRUE),
      AveExpr = NA_real_,
      pvalue = suppressWarnings(pmax(a$pvalue[match(g, a$gene)], b$pvalue[match(g, b$gene)], na.rm = TRUE)),
      padj = NA_real_,
      padj_BH = NA_real_,
      stringsAsFactors = FALSE
    )
  }, nerve_de_list)
} else {
  empty_de()
}

path_up_n <- function(g) {
  if (length(nerve_de_list) == 0) return(0L)
  as.integer(sum(vapply(nerve_de_list, function(d) {
    i <- match(g, d$gene)
    if (is.na(i)) return(FALSE)
    ok <- !is.na(d$log2FC[i]) && d$log2FC[i] > 0
    if (!is.na(d$pvalue[i])) ok <- ok && d$pvalue[i] < p_cutoff
    ok
  }, logical(1))))
}

if (nrow(q2_src) > 0) {
  q2 <- q2_src
  q2$FC <- 2^q2$log2FC
  q2$known_tumor_to_nerve_ligand <- q2$gene %in% tumor_ligands$gene
  q2$ligand_receptor <- tumor_ligands$receptor_on_nerve[match(q2$gene, tumor_ligands$gene)]
  q2$ligand_rationale <- tumor_ligands$rationale[match(q2$gene, tumor_ligands$gene)]
  q2$mean_spearman_near_schwann <- agg$mean_spearman_near_schwann[match(q2$gene, agg$gene)]
  q2$n_pathologist_patients_also_up <- vapply(q2$gene, path_up_n, integer(1))
  q2$pass_p005_FC_1.25 <- q2$gene %in% select_up(q2_src, 1.25)$gene
  q2$pass_p005_FC_1.5  <- q2$gene %in% select_up(q2_src, 1.5)$gene
  q2 <- q2[order(is.na(q2$pvalue), q2$pvalue, -q2$log2FC), ]
  q2$genome_wide_rank <- seq_len(nrow(q2))
  utils::write.csv(q2, file.path(result_dir, "02_Q2_ALL_tumor_genes_near_nerve.csv"),
                   row.names = FALSE)
  up125 <- q2[q2$pass_p005_FC_1.25, , drop = FALSE]
  up15  <- q2[q2$pass_p005_FC_1.5, , drop = FALSE]
  utils::write.csv(up125, file.path(result_dir, "02_Q2_up_p005_FC1.25.csv"), row.names = FALSE)
  utils::write.csv(up15,  file.path(result_dir, "02_Q2_up_p005_FC1.5.csv"), row.names = FALSE)
  short <- up125
  if (nrow(short) > 0) {
    short <- short[order(-short$known_tumor_to_nerve_ligand, short$pvalue, -short$log2FC), ]
  }
  utils::write.csv(short, file.path(result_dir, "02_Q2_READ_THIS_tumor_genes_may_promote_nerve.csv"),
                   row.names = FALSE)
  log_msg("问题2 上调 FC>=1.25 n=", nrow(up125), " FC>=1.5 n=", nrow(up15))
} else {
  log_msg("问题2：没有可用的差异表。请确认 Robjects/counts 已解压且病人数足够。")
}

# 问题 3：靠近神经的肿瘤细胞里降低的基因（丢失排斥/屏障，可能促使浸润）
path_down_n <- function(g) {
  if (length(nerve_de_list) == 0) return(0L)
  as.integer(sum(vapply(nerve_de_list, function(d) {
    i <- match(g, d$gene)
    if (is.na(i)) return(FALSE)
    ok <- !is.na(d$log2FC[i]) && d$log2FC[i] < 0
    if (!is.na(d$pvalue[i])) ok <- ok && d$pvalue[i] < p_cutoff
    ok
  }, logical(1))))
}

if (nrow(q2_src) > 0) {
  q3 <- q2_src
  q3$FC <- 2^q3$log2FC
  q3$down_FC <- 1 / pmax(q3$FC, 1e-12)
  q3$known_nerve_repellent_or_barrier <- q3$gene %in% tumor_repellents$gene
  q3$repellent_partner <- tumor_repellents$receptor_or_partner[match(q3$gene, tumor_repellents$gene)]
  q3$repellent_rationale <- tumor_repellents$rationale[match(q3$gene, tumor_repellents$gene)]
  q3$n_pathologist_patients_also_down <- vapply(q3$gene, path_down_n, integer(1))
  q3$pass_p005_down_FC_1.25 <- q3$gene %in% select_down(q2_src, 1.25)$gene
  q3$pass_p005_down_FC_1.5  <- q3$gene %in% select_down(q2_src, 1.5)$gene
  q3 <- q3[order(is.na(q3$pvalue), q3$pvalue, q3$log2FC), ]
  q3$genome_wide_rank <- seq_len(nrow(q3))
  utils::write.csv(q3, file.path(result_dir, "03_Q3_ALL_tumor_genes_near_nerve.csv"),
                   row.names = FALSE)
  dn125 <- q3[q3$pass_p005_down_FC_1.25, , drop = FALSE]
  dn15  <- q3[q3$pass_p005_down_FC_1.5, , drop = FALSE]
  utils::write.csv(dn125, file.path(result_dir, "03_Q3_down_p005_FC1.25.csv"), row.names = FALSE)
  utils::write.csv(dn15,  file.path(result_dir, "03_Q3_down_p005_FC1.5.csv"), row.names = FALSE)
  short3 <- dn125
  if (nrow(short3) > 0) {
    short3 <- short3[order(-short3$known_nerve_repellent_or_barrier, short3$pvalue, short3$log2FC), ]
  }
  utils::write.csv(short3, file.path(result_dir, "03_Q3_READ_THIS_tumor_genes_low_may_promote_nerve.csv"),
                   row.names = FALSE)
  log_msg("问题3 下调 FC>=1.25 n=", nrow(dn125), " FC>=1.5 n=", nrow(dn15))
} else {
  log_msg("问题3：没有可用的差异表。")
}

# -----------------------------------------------------------------------------
# 9. 方案说明（问题 1）
# -----------------------------------------------------------------------------
protocol <- c(
  "============================================================",
  "乳腺癌空间转录组：神经浸润（Wang et al. Nat Commun 2024）",
  "数据目录: ", nerve_dir,
  "不做 GO / KEGG / GSEA。",
  "============================================================",
  "",
  "【问题 1】神经浸润相关分子和方案",
  "",
  "瘤内几乎没有神经元胞体，主要是纤维 + Schwann。病理 Nerve 与 ST spot",
  "重叠的病人极少（见 00_logs/spot_class_counts.csv），所以双轨：",
  "  轨1 病理 Nerve 邻域：有神经的病人各自 1-vs-1，不合并",
  "  轨2 Schwann 签名邻域：全队列伪 bulk ~ patient + group（主分析）",
  "",
  "Spot：Nerve 像素>=1；肿瘤 = Tumor 比例>=0.25（否则 EPCAM/KRT）；",
  "近 <=2 个 spot，远 >4；Schwann 高 = 切片内 z>1",
  "（SOX10 MPZ PMP22 S100B PLP1 NGFR NCAM1 MBP L1CAM）。",
  "",
  "统计：先 p<0.05，再 |FC|>=1.25 和 1.5。",
  "  问题2 = 近神经肿瘤上调；问题3 = 近神经肿瘤下调。",
  "过滤低表达后再算 FC。",
  "",
  "文献配体（肿瘤高表达→神经，问题2）：NGF, BDNF, ARTN, GDNF, CXCL12, MDK, PTN, VEGFA。",
  "表：01_Q1_literature_ligands.csv",
  "文献排斥/屏障（肿瘤低表达→放开浸润，问题3）：SEMA3A/F, SLIT2/3, CDH1, TJP1, COL4A1, LAMA2。",
  "表：01_Q3_literature_repellents.csv",
  "",
  "【问题 2】肿瘤细胞高表达、可能促使神经浸润的基因",
  "",
  "  02_Q2_READ_THIS_tumor_genes_may_promote_nerve.csv",
  "    p<0.05 且上调 FC>=1.25。已知促神经配体排在前面。",
  "  02_Q2_up_p005_FC1.25.csv / 02_Q2_up_p005_FC1.5.csv",
  "",
  "【问题 3】肿瘤细胞低表达、可能促使神经浸润的基因",
  "",
  "  03_Q3_READ_THIS_tumor_genes_low_may_promote_nerve.csv",
  "    p<0.05 且下调倍数>=1.25（近/远 FC<=1/1.25）。",
  "    已知轴突排斥 / 上皮屏障基因排在前面。",
  "  03_Q3_down_p005_FC1.25.csv / 03_Q3_down_p005_FC1.5.csv",
  "  03_Q3_ALL_tumor_genes_near_nerve.csv  全基因（未改 p 值）",
  "",
  "比较定义：靠近 Schwann/神经的肿瘤 spot vs 远离的肿瘤 spot。",
  "邻域升高或降低 = 空间相邻，还不等于已经证明『引起』浸润。",
  "问题2 功能验证：肿瘤细胞过表达后与 DRG 共培养看轴突是否更长。",
  "问题3 功能验证：把该基因补回去（过表达）看轴突是否被挡住。",
  "",
  "病理神经各病人表在 TNBC*_tumor_near_nerve_vs_far/，只作验证。"
)
writeLines(protocol, file.path(result_dir, "00_Q1_PROTOCOL.txt"))
writeLines(c(
  "先看这些文件：",
  "  问题1  00_Q1_PROTOCOL.txt",
  "         01_Q1_literature_ligands.csv（高表达促浸润配体）",
  "         01_Q3_literature_repellents.csv（低表达放开浸润）",
  "  问题2  02_Q2_READ_THIS_tumor_genes_may_promote_nerve.csv",
  "  问题3  03_Q3_READ_THIS_tumor_genes_low_may_promote_nerve.csv",
  "没有 GO/GSEA。"
), file.path(result_dir, "00_READ_ME.txt"))
log_msg("Wrote Q1 protocol and Q2/Q3 ranked tumor-gene tables")
log_msg("问题1: ", file.path(result_dir, "00_Q1_PROTOCOL.txt"))
log_msg("问题2: ", file.path(result_dir, "02_Q2_READ_THIS_tumor_genes_may_promote_nerve.csv"))
log_msg("问题3: ", file.path(result_dir, "03_Q3_READ_THIS_tumor_genes_low_may_promote_nerve.csv"))

