#!/usr/bin/env Rscript
# =============================================================================
# Wang et al. Nature Communications 2024  TNBC 空间转录组
# 远神经肿瘤组织 vs 近神经肿瘤组织（上调 = 远神经肿瘤更高）
#
# 本脚本只做这件事，不做 GO / KEGG / GSEA / ORA，不做下调，不做 FC=2，不做 topN。
# 独立于 TG_RNAseq_pipeline.R，不改原 Cuffdiff 1–4。
# 数据：E:/R/Nerve 的 Robjects/（counts + annotsBySpot）。
# ids.RDS 可选（Windows 常把 Clinical/ids.RDS 展成 Clinicalids.RDS）。
#
# 用法（Windows R / RStudio）：
#   setwd("E:/R/Nerve")
#   Sys.setenv(WANG_ST_DIR = "E:/R/Nerve")
#   source("Wang_ST_nerve_infiltration.R")
#   # 不要 source Nerve_RNA_primary_lung_tropism.R（那是另一套 GPL96 分析）
#
# FC 只做 >1、1.25、1.5。
# 符合条件的单个病人：results/00_eligible_single_patients/
# 综合上调基因：      results/01_combined_far_tumor_vs_near_tumor/
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

# 远神经肿瘤 vs 近神经肿瘤：全队列脚本默认 FC>1、1.25、1.5
p_cutoff <- 0.05
fc_cutoffs_up <- c("FC_gt_1" = 1, "FC_1.25" = 1.25, "FC_1.5" = 1.5)

wang_st_skip_main <- isTRUE(getOption("wang.st.skip_main", FALSE))

if (!wang_st_skip_main) {
result_dir <- file.path(nerve_dir, "results")
log_dir <- file.path(result_dir, "00_logs")
expr_dir <- file.path(result_dir, "00_eligible_single_patients")
combined_dir <- file.path(result_dir, "01_combined_far_tumor_vs_near_tumor")

# 删掉旧版「近 vs 远 / FC=2 / TopRank」目录，避免和这次要求混在一起
wipe_stale_wang_results <- function(root) {
  if (!dir.exists(root)) return(invisible(NULL))
  stale_names <- c(
    "tumor_near_schwann_vs_far",
    "tumor_near_nerve_vs_far",
    "pathologist_far_tumor_vs_near_tumor",
    "01_Q3_literature_repellents.csv",
    "01_CANDIDATE_MOLECULES_tumor_to_nerve.csv",
    "00_PROTOCOL_神经浸润分析方案.txt",
    "00_Q1_PROTOCOL.txt",
    "01_Q1_literature_ligands.csv",
    "01_Q3_literature_repellents.csv"
  )
  for (nm in stale_names) {
    p <- file.path(root, nm)
    if (file.exists(p) || dir.exists(p)) unlink(p, recursive = TRUE, force = TRUE)
  }
  extra <- list.files(root, full.names = TRUE, include.dirs = TRUE)
  extra <- extra[grepl("near_(nerve|schwann)_vs_far|tumor_near_.*vs_far",
                       basename(extra), ignore.case = TRUE)]
  if (length(extra) > 0) unlink(extra, recursive = TRUE, force = TRUE)
  invisible(NULL)
}
wipe_stale_wang_results(result_dir)

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result_dir, "00_QC_maps"), recursive = TRUE, showWarnings = FALSE)
dir.create(expr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(c(
  paste0("符合条件的单个病人：", expr_dir),
  "  打开 符合条件病人_打开这个.csv",
  paste0("远神经肿瘤 vs 近神经肿瘤 上调基因：", combined_dir),
  "  upregulated_far_tumor_vs_near_tumor_FC_gt_1.csv",
  "  upregulated_far_tumor_vs_near_tumor_FC_1.25.csv",
  "  upregulated_far_tumor_vs_near_tumor_FC_1.5.csv",
  "FC 只有 >1、1.25、1.5。不要看 tumor_near_schwann_vs_far / FoldChange / TopRank。"
), file.path(result_dir, "00_请先看这里_单个病人和上调基因.txt"))
log_file <- file.path(log_dir, paste0("wang_st_nerve_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg_file <- log_msg
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_msg("Wang ST nerve dir: ", nerve_dir)
log_msg("Results: ", result_dir)
} # end if (!wang_st_skip_main) 路径/日志初始化

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
  group <- factor(group, levels = c("near", "far"))
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
  coefn <- grep("groupfar", colnames(design), value = TRUE)
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
  si$group <- factor(si$group, levels = c("near", "far"))
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
  res <- DESeq2::results(dds, contrast = c("group", "far", "near"))
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
  log2FC <- colMeans(logmat[far_id, , drop = FALSE]) - colMeans(logmat[near_id, , drop = FALSE])
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
  # fc=1 → FC>1；fc=1.25/1.5 → FC >= 该档。log2FC>0 = 远神经肿瘤高于近神经肿瘤
  if (isTRUE(all.equal(as.numeric(fc), 1))) {
    keep <- !is.na(de$log2FC) & de$log2FC > 0
  } else {
    keep <- !is.na(de$log2FC) & de$log2FC > 0 & (2^de$log2FC >= fc)
  }
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

# -----------------------------------------------------------------------------
# 6. 差异表（只出 CSV + 火山图，不出 GO/GSEA；只出远 vs 近上调）
# -----------------------------------------------------------------------------
basic_volcano <- function(de, title, outfile, fc_line = 1) {
  if (nrow(de) == 0) return(invisible(NULL))
  df <- de
  df$y <- if (any(!is.na(df$pvalue))) -log10(pmax(df$pvalue, 1e-300)) else abs(df$log2FC)
  df$col <- "ns"
  sig <- is.na(df$pvalue) | df$pvalue < p_cutoff
  df$col[sig & !is.na(df$log2FC) & df$log2FC > 0 & df$log2FC >= log2(pmax(fc_line, 1 + 1e-8))] <- "up"
  if (isTRUE(all.equal(as.numeric(fc_line), 1))) {
    df$col[sig & !is.na(df$log2FC) & df$log2FC > 0] <- "up"
  }
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(log2FC, y, color = col)) +
      ggplot2::geom_point(alpha = 0.5, size = 0.8) +
      ggplot2::scale_color_manual(values = c(ns = "grey70", up = "#D62828")) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(title = title, y = "-log10(p)",
                    x = "log2FC (far-nerve tumor / near-nerve tumor)")
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

emit_de_tables <- function(comp_name, de, out_dir = NULL, also_top_level = FALSE) {
  de <- de[!is.na(de$log2FC), , drop = FALSE]
  base <- if (is.null(out_dir)) file.path(result_dir, comp_name) else out_dir
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(de, file.path(base, "DE_full_far_tumor_vs_near_tumor.csv"), row.names = FALSE)
  n_up <- sum((is.na(de$pvalue) | de$pvalue < p_cutoff) & de$log2FC > 0, na.rm = TRUE)
  log_msg(comp_name, " 远神经肿瘤 vs 近神经肿瘤  genes=", nrow(de),
          " 上调 p<", p_cutoff, " 且 FC>1 n=", n_up)
  for (nm in names(fc_cutoffs_up)) {
    fc <- unname(fc_cutoffs_up[[nm]])
    up <- select_up(de, fc)
    od <- file.path(base, paste0("up_", nm))
    dir.create(od, recursive = TRUE, showWarnings = FALSE)
    fname <- paste0("upregulated_far_tumor_vs_near_tumor_", nm, ".csv")
    utils::write.csv(up, file.path(od, fname), row.names = FALSE)
    if (isTRUE(also_top_level)) {
      utils::write.csv(up, file.path(base, fname), row.names = FALSE)
    }
    lab <- if (isTRUE(all.equal(fc, 1))) "FC>1" else paste0("FC>=", fc)
    basic_volcano(de, paste(comp_name, "| far-nerve tumor vs near-nerve tumor upregulated", lab),
                  file.path(od, paste0(nm, "_volcano")), fc)
    log_msg("  上调 ", lab, " n=", nrow(up))
  }
  invisible(de)
}

write_gene_mat <- function(mat, path) {
  utils::write.csv(
    data.frame(gene = rownames(mat), mat, check.names = FALSE, stringsAsFactors = FALSE),
    path, row.names = FALSE
  )
}

bind_gene_col <- function(mat, vec, colname) {
  vec <- vec[is.finite(vec)]
  if (length(vec) == 0) return(mat)
  if (is.null(mat)) {
    return(matrix(as.numeric(vec), ncol = 1,
                  dimnames = list(names(vec), colname)))
  }
  gn <- intersect(rownames(mat), names(vec))
  if (length(gn) == 0) return(mat)
  out <- cbind(mat[gn, , drop = FALSE], as.numeric(vec[gn]))
  colnames(out)[ncol(out)] <- colname
  out
}

bind_pb <- function(acc, pb, si) {
  if (is.null(acc$mat)) return(list(mat = pb, si = si))
  gn <- intersect(rownames(acc$mat), rownames(pb))
  list(
    mat = cbind(acc$mat[gn, , drop = FALSE], pb[gn, , drop = FALSE]),
    si = rbind(acc$si, si)
  )
}

# 每个病人：近神经肿瘤 vs 远神经肿瘤的基因表达（RNA count / mean logCPM）
export_patient_near_far <- function(pid, track, logm, cnts, near_i, far_i, eligible_de) {
  if (length(near_i) < 1 || length(far_i) < 1) return(NULL)
  od <- file.path(expr_dir, track, paste0("TNBC", pid))
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  near_mean <- colMeans(logm[near_i, , drop = FALSE])
  far_mean  <- colMeans(logm[far_i, , drop = FALSE])
  near_pb <- colSums(cnts[near_i, , drop = FALSE])
  far_pb  <- colSums(cnts[far_i, , drop = FALSE])
  tab <- data.frame(
    gene = colnames(logm),
    mean_logCPM_near_nerve_tumor = as.numeric(near_mean),
    mean_logCPM_far_nerve_tumor = as.numeric(far_mean),
    pseudobulk_count_near = as.numeric(near_pb),
    pseudobulk_count_far = as.numeric(far_pb),
    n_spots_near = length(near_i),
    n_spots_far = length(far_i),
    log2FC_far_minus_near = as.numeric(far_mean - near_mean),
    eligible_for_combined_DE = eligible_de,
    stringsAsFactors = FALSE
  )
  utils::write.csv(tab, file.path(od, "near_and_far_tumor_RNA_expression.csv"),
                   row.names = FALSE)
  meta <- data.frame(
    patient = pid,
    track = track,
    n_spots_near = length(near_i),
    n_spots_far = length(far_i),
    eligible_for_combined_DE = eligible_de,
    note = "Wang ST 是基因表达（RNA），不是蛋白质组",
    stringsAsFactors = FALSE
  )
  utils::write.csv(meta, file.path(od, "sample_info.csv"), row.names = FALSE)
  list(
    mean_near = near_mean, mean_far = far_mean,
    pb_near = near_pb, pb_far = far_pb
  )
}

# -----------------------------------------------------------------------------
# 7. 主分析（单独子集脚本可设 options(wang.st.skip_main=TRUE) 后只加载函数）
# -----------------------------------------------------------------------------
if (!isTRUE(getOption("wang.st.skip_main", FALSE))) {
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
pb_nerve_acc <- list(mat = NULL, si = NULL)
mean_nerve_wide <- NULL
mean_schwann_wide <- NULL
pb_nerve_wide <- NULL
pb_schwann_wide <- NULL
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
  tag <- paste0("TNBC", pid)

  # 病理神经：只导出符合条件的病人（近>=5 且 远>=10）
  n_near_n <- which(sp$near_nerve)
  n_far_n  <- which(sp$far_nerve)
  elig_nerve <- length(n_near_n) >= 5 && length(n_far_n) >= 10
  if (elig_nerve) {
    ex <- export_patient_near_far(pid, "pathologist_nerve", logm, cnts,
                                  n_near_n, n_far_n, TRUE)
    if (!is.null(ex)) {
      mean_nerve_wide <- bind_gene_col(mean_nerve_wide, ex$mean_near, paste0(tag, "_near"))
      mean_nerve_wide <- bind_gene_col(mean_nerve_wide, ex$mean_far,  paste0(tag, "_far"))
      pb_n <- cbind(ex$pb_near, ex$pb_far)
      rownames(pb_n) <- names(ex$pb_near)
      colnames(pb_n) <- paste0(tag, c("_near", "_far"))
      pb_nerve_wide <- if (is.null(pb_nerve_wide)) pb_n else {
        gn <- intersect(rownames(pb_nerve_wide), rownames(pb_n))
        cbind(pb_nerve_wide[gn, , drop = FALSE], pb_n[gn, , drop = FALSE])
      }
      si_n <- data.frame(
        sample = colnames(pb_n), group = c("near", "far"),
        patient = as.character(pid), stringsAsFactors = FALSE
      )
      pb_nerve_acc <- bind_pb(pb_nerve_acc, pb_n, si_n)
    }
    grp <- ifelse(sp$near_nerve, "near", ifelse(sp$far_nerve, "far", NA))
    keep <- !is.na(grp)
    de_i <- limma_two_group(logm[keep, , drop = FALSE], grp[keep],
                            label = paste0("TNBC", pid, "_nerve"))
    if (nrow(de_i) == 0) {
      de_i <- mean_fc_two_group(logm, n_near_n, n_far_n)
      log_msg("TNBC", pid, " nerve: limma empty, fallback mean FC (no p)")
    }
    nm <- paste0("TNBC", pid, "_far_tumor_vs_near_tumor")
    nerve_de_list[[nm]] <- de_i
    emit_de_tables(nm, de_i,
                   out_dir = file.path(expr_dir, "pathologist_nerve", tag, "DE_far_tumor_vs_near_tumor"))
  }

  # Schwann 邻域：只导出符合条件的病人
  n_near_s <- which(sp$near_schwann)
  n_far_s  <- which(sp$far_schwann)
  elig_schw <- length(n_near_s) >= 5 && length(n_far_s) >= 10
  if (elig_schw) {
    exs <- export_patient_near_far(pid, "schwann_neighborhood", logm, cnts,
                                   n_near_s, n_far_s, TRUE)
    if (!is.null(exs)) {
      mean_schwann_wide <- bind_gene_col(mean_schwann_wide, exs$mean_near, paste0(tag, "_near"))
      mean_schwann_wide <- bind_gene_col(mean_schwann_wide, exs$mean_far,  paste0(tag, "_far"))
    }
    grp_s <- ifelse(sp$near_schwann, "near", ifelse(sp$far_schwann, "far", NA))
    keep_s <- !is.na(grp_s)
    de_s <- limma_two_group(logm[keep_s, , drop = FALSE], grp_s[keep_s],
                            label = paste0("TNBC", pid, "_schwann"))
    if (nrow(de_s) == 0) {
      de_s <- mean_fc_two_group(logm, n_near_s, n_far_s)
    }
    emit_de_tables(paste0("TNBC", pid, "_far_tumor_vs_near_tumor"), de_s,
                   out_dir = file.path(expr_dir, "schwann_neighborhood", tag,
                                       "DE_far_tumor_vs_near_tumor"))
    near_sum <- colSums(cnts[n_near_s, , drop = FALSE])
    far_sum  <- colSums(cnts[n_far_s, , drop = FALSE])
    pb <- cbind(near_sum, far_sum)
    colnames(pb) <- paste0(tag, c("_near", "_far"))
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
if (!is.null(qc) && nrow(qc) > 0) {
  qc$eligible_pathologist_nerve_DE <- qc$n_tumor_near_nerve >= 5 & qc$n_tumor_far_nerve >= 10
  qc$eligible_schwann_DE <- qc$n_tumor_near_schwann >= 5 & qc$n_tumor_far_schwann >= 10
  qc$folder_pathologist <- ifelse(
    qc$eligible_pathologist_nerve_DE,
    file.path(expr_dir, "pathologist_nerve", paste0("TNBC", qc$patient)),
    ""
  )
  qc$folder_schwann <- ifelse(
    qc$eligible_schwann_DE,
    file.path(expr_dir, "schwann_neighborhood", paste0("TNBC", qc$patient)),
    ""
  )
}
utils::write.csv(qc, file.path(log_dir, "spot_class_counts.csv"), row.names = FALSE)
utils::write.csv(qc, file.path(expr_dir, "INDEX_eligible_patients.csv"), row.names = FALSE)
if (!is.null(qc) && nrow(qc) > 0) {
  hit <- qc[qc$eligible_pathologist_nerve_DE | qc$eligible_schwann_DE, , drop = FALSE]
  utils::write.csv(hit, file.path(expr_dir, "符合条件病人_打开这个.csv"), row.names = FALSE)
  log_msg("Patients with pathologist nerve spots: ",
          paste(qc$patient[qc$n_nerve_path > 0], collapse = ", "))
  log_msg("符合条件（病理神经 远vs近）病人: ",
          paste(qc$patient[qc$eligible_pathologist_nerve_DE], collapse = ", "))
  log_msg("符合条件（Schwann 远vs近）病人: ",
          paste(qc$patient[qc$eligible_schwann_DE], collapse = ", "))
  log_msg("单个病人文件夹: ", expr_dir)
}

write_combined_expression <- function(track, mean_wide, pb_mat, si) {
  od <- file.path(expr_dir, track)
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(mean_wide) && ncol(mean_wide) > 0) {
    write_gene_mat(mean_wide, file.path(od, "COMBINED_mean_logCPM_near_and_far.csv"))
    tryCatch(writexl::write_xlsx(
      data.frame(gene = rownames(mean_wide), mean_wide, check.names = FALSE),
      file.path(od, "COMBINED_mean_logCPM_near_and_far.xlsx")
    ), error = function(e) NULL)
  }
  if (!is.null(pb_mat) && ncol(pb_mat) > 0) {
    write_gene_mat(pb_mat, file.path(od, "COMBINED_pseudobulk_counts_near_and_far.csv"))
    if (!is.null(si)) {
      utils::write.csv(si, file.path(od, "COMBINED_sample_info.csv"), row.names = FALSE)
    }
  }
}

write_combined_expression("pathologist_nerve", mean_nerve_wide, pb_nerve_wide, pb_nerve_acc$si)
si_sw_export <- si_schwann
write_combined_expression("schwann_neighborhood", mean_schwann_wide, pb_schwann, si_sw_export)
where_txt <- c(
  "比较：远神经肿瘤组织 vs 近神经肿瘤组织。",
  "上调 = 远神经肿瘤更高（log2FC = 远 - 近）。",
  "FC 只有三档：>1、1.25、1.5。没有 FC=2，没有 topN，没有下调。",
  "这是空间转录组 RNA，不是蛋白。",
  "",
  paste0("【符合条件的单个病人】", expr_dir),
  "  先打开：符合条件病人_打开这个.csv",
  "  或：INDEX_eligible_patients.csv（eligible_* = TRUE 才有文件夹）",
  "  条件：近神经肿瘤 spot >= 5 且 远神经肿瘤 spot >= 10。",
  "",
  "  每人近/远肿瘤表达：",
  "    schwann_neighborhood/TNBC数字/near_and_far_tumor_RNA_expression.csv",
  "  每人上调基因（远 vs 近）：",
  "    schwann_neighborhood/TNBC数字/DE_far_tumor_vs_near_tumor/up_FC_gt_1/",
  "    schwann_neighborhood/TNBC数字/DE_far_tumor_vs_near_tumor/up_FC_1.25/",
  "    schwann_neighborhood/TNBC数字/DE_far_tumor_vs_near_tumor/up_FC_1.5/",
  "  病理神经（病人很少，常见只有 TNBC50）：",
  "    pathologist_nerve/TNBC数字/  同样结构",
  "",
  paste0("【这些人合在一起的上调基因】", combined_dir),
  "  upregulated_far_tumor_vs_near_tumor_FC_gt_1.csv",
  "  upregulated_far_tumor_vs_near_tumor_FC_1.25.csv",
  "  upregulated_far_tumor_vs_near_tumor_FC_1.5.csv",
  "",
  "不要看旧文件夹 tumor_near_schwann_vs_far / FoldChange / TopRank，脚本会删掉它们。",
  "不要 source Nerve_RNA_primary_lung_tropism.R（那是另一套 GPL96 分析）。"
)
writeLines(where_txt, file.path(expr_dir, "00_READ_ME.txt"))
writeLines(where_txt, file.path(result_dir, "00_请先看这里_单个病人和上调基因.txt"))
writeLines(where_txt, file.path(combined_dir, "00_READ_ME.txt"))

# 共同上调（病理神经，仅当 >=2 个病人）：远 vs 近，只报 FC>1 / 1.25 / 1.5
if (length(nerve_de_list) >= 2) {
  for (nm in names(fc_cutoffs_up)) {
    fc <- unname(fc_cutoffs_up[[nm]])
    inter <- Reduce(intersect, lapply(nerve_de_list, function(d) select_up(d, fc)$gene))
    log_msg("病理神经共同上调 远vs近 ", nm, " n=", length(inter))
  }
} else {
  log_msg("病理神经可做远vs近的病人不足 2 例。单个病人表在 00_eligible_single_patients/pathologist_nerve/")
}

# 能做远 vs 近的病人合在一起（Schwann 邻域，~ patient + group）
if (!is.null(pb_schwann) && length(unique(si_schwann$patient)) >= 2) {
  log_msg("综合 远神经肿瘤 vs 近神经肿瘤，病人: ",
          paste(unique(si_schwann$patient), collapse = ","))
  rownames(si_schwann) <- si_schwann$sample
  de_sw <- deseq2_pseudobulk(round(pmax(pb_schwann, 0)), si_schwann,
                             "far_tumor_vs_near_tumor")
  if (nrow(de_sw) == 0) {
    logm <- log2(sweep(t(pb_schwann), 1, pmax(colSums(pb_schwann), 1), "/") * 1e6 + 1)
    de_sw <- limma_two_group(logm, si_schwann$group, si_schwann$patient, "schwann_pb")
  }
  emit_de_tables("far_tumor_vs_near_tumor", de_sw,
                 out_dir = combined_dir, also_top_level = TRUE)
} else {
  log_msg("能做远 vs 近的病人不足，跳过综合分析")
  de_sw <- empty_de()
}

de_nerve_combined <- empty_de()
if (!is.null(pb_nerve_acc$mat) && length(unique(pb_nerve_acc$si$patient)) >= 2) {
  si_n <- pb_nerve_acc$si
  rownames(si_n) <- si_n$sample
  de_nerve_combined <- deseq2_pseudobulk(round(pmax(pb_nerve_acc$mat, 0)), si_n,
                                         "pathologist_far_vs_near")
  if (nrow(de_nerve_combined) == 0) {
    logm <- log2(sweep(t(pb_nerve_acc$mat), 1, pmax(colSums(pb_nerve_acc$mat), 1), "/") * 1e6 + 1)
    de_nerve_combined <- limma_two_group(logm, si_n$group, si_n$patient, "nerve_pb")
  }
  emit_de_tables("pathologist_far_tumor_vs_near_tumor", de_nerve_combined,
                 out_dir = file.path(combined_dir, "pathologist_if_ge2_patients"))
} else {
  log_msg("病理远vs近可合并病人不足 2 例，综合分析用 Schwann 邻域。")
}

# -----------------------------------------------------------------------------
# 文献配体对照表（可选；主结果仍是上面的远 vs 近上调基因）
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
  cand$log2FC_far_vs_near <- de_sw$log2FC[match(cand$gene, de_sw$gene)]
  cand$pvalue_far_vs_near <- de_sw$pvalue[match(cand$gene, de_sw$gene)]
} else {
  cand$log2FC_far_vs_near <- NA_real_
  cand$pvalue_far_vs_near <- NA_real_
}
cand$up_in_far_vs_near <- (!is.na(cand$log2FC_far_vs_near) & cand$log2FC_far_vs_near > 0) &
  (is.na(cand$pvalue_far_vs_near) | cand$pvalue_far_vs_near < p_cutoff)
utils::write.csv(cand, file.path(result_dir, "02_optional_literature_ligands.csv"),
                 row.names = FALSE)

log_msg("单个病人: ", expr_dir)
log_msg("综合上调: ", combined_dir)
log_msg("先打开: ", file.path(result_dir, "00_请先看这里_单个病人和上调基因.txt"))
} # end if (!wang.st.skip_main) 全队列主分析

