#!/usr/bin/env Rscript
# =============================================================================
# 单独运行：指定病人综合  远神经肿瘤组织 vs 近神经肿瘤组织
# 上调只做 FC>1、1.25。不做 FC=1.5、不做 FC=2、不做 topN、不做下调、不做 GO。
# 近/远：近 ≤2 spot，远 ≥10 spot（这批病人统一；旧版远>4 间隔太小）。
# 不改 Wang_ST_nerve_infiltration.R 的全队列结果（00_ / 01_ 文件夹不动）。
#
# 用法（Windows R / RStudio）：
#   setwd("E:/R/Nerve")
#   Sys.setenv(WANG_ST_DIR = "E:/R/Nerve")
#   source("Wang_ST_far_vs_near_selected_patients.R")
#
# 本脚本和 Wang_ST_nerve_infiltration.R 放同一目录。
# 结果：results/02_selected32_far_tumor_vs_near_tumor/
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)

# 病人 2,3,4,6,11,13,14,15,20,22,27,28,30,31,37,39,51,52,53,56,61,62,67,69,74,79,81,85,86,90,93,94
selected_patients <- c(
  2, 3, 4, 6, 11, 13, 14, 15, 20, 22, 27, 28, 30, 31, 37, 39,
  51, 52, 53, 56, 61, 62, 67, 69, 74, 79, 81, 85, 86, 90, 93, 94
)
ov <- Sys.getenv("WANG_ST_SELECTED_PIDS", unset = "")
if (nzchar(ov)) {
  selected_patients <- as.integer(unlist(strsplit(ov, "[,[:space:]]+")))
  selected_patients <- selected_patients[is.finite(selected_patients) & selected_patients > 0]
}

locate_wang_main <- function() {
  ofile <- NULL
  args <- commandArgs(trailingOnly = FALSE)
  farg <- grep("^--file=", args, value = TRUE)
  if (length(farg) > 0) ofile <- sub("^--file=", "", farg[1])
  if (is.null(ofile) || !nzchar(ofile)) {
    n <- sys.nframe()
    if (n >= 1) {
      for (i in n:1) {
        of <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
        if (!is.null(of) && nzchar(of)) {
          ofile <- of
          break
        }
      }
    }
  }
  here <- if (!is.null(ofile) && nzchar(ofile)) {
    dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE))
  } else {
    getwd()
  }
  cands <- unique(c(
    file.path(here, "Wang_ST_nerve_infiltration.R"),
    file.path(getwd(), "Wang_ST_nerve_infiltration.R"),
    "E:/R/Nerve/Wang_ST_nerve_infiltration.R",
    "E:/R/TG_BRCA/TG/Wang_ST_nerve_infiltration.R"
  ))
  hit <- cands[file.exists(cands)]
  if (length(hit) == 0) {
    stop("找不到 Wang_ST_nerve_infiltration.R。请和本脚本放在同一目录后再 source。")
  }
  hit[1]
}

old_opt <- getOption("wang.st.skip_main")
options(wang.st.skip_main = TRUE)
main_r <- locate_wang_main()
source(main_r, local = FALSE)
options(wang.st.skip_main = old_opt)

# 本脚本只要上调 FC>1、1.25；距离图只对 FC>1.25 的每个上调基因单独画
p_cutoff <- 0.05
fc_cutoffs_up <- c("FC_gt_1" = 1, "FC_1.25" = 1.25)
# Wang ST spot 中心距约 100 μm；距离用 spot 网格换算成 mm
spot_pitch_mm <- 0.1
# 指定病人统一用更大间隔：近 ≤2 spot，远 ≥10 spot（旧版近≤2、远>4 间隔太小）
# 可用环境变量覆盖：WANG_ST_NEAR_SPOTS / WANG_ST_FAR_SPOTS
near_spots <- 2
far_spots <- 10
env_near <- suppressWarnings(as.numeric(Sys.getenv("WANG_ST_NEAR_SPOTS", unset = "")))
env_far <- suppressWarnings(as.numeric(Sys.getenv("WANG_ST_FAR_SPOTS", unset = "")))
if (is.finite(env_near) && env_near > 0) near_spots <- env_near
if (is.finite(env_far) && env_far > near_spots) far_spots <- env_far
# 距离参照：施旺细胞签名
schwann_distance_markers <- schwann_genes

apply_selected_near_far <- function(sp, near_d, far_d) {
  d_s <- as.numeric(sp$dist_schwann)
  d_n <- as.numeric(sp$dist_nerve)
  sp$near_schwann <- sp$is_tumor & !sp$is_schwann_high & is.finite(d_s) & d_s <= near_d
  sp$far_schwann  <- sp$is_tumor & !sp$is_schwann_high & is.finite(d_s) & d_s >= far_d
  sp$near_nerve <- sp$is_tumor & !sp$is_nerve_path & is.finite(d_n) & d_n <= near_d
  sp$far_nerve  <- sp$is_tumor & !sp$is_nerve_path & is.finite(d_n) & d_n >= far_d
  sp
}

out_root <- file.path(nerve_dir, "results", "02_selected32_far_tumor_vs_near_tumor")
expr_dir <- file.path(out_root, "eligible_single_patients")
log_dir <- file.path(out_root, "00_logs")
dir.create(expr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
flag_no_de <- file.path(out_root, "NO_COMBINED_DE.txt")
if (file.exists(flag_no_de)) unlink(flag_no_de)

log_file <- file.path(log_dir, paste0("selected32_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

label_high_low <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  if (length(x) < 8) return(rep(NA_character_, length(x)))
  if (mean(x <= 0) >= 0.5) {
    ifelse(x > 0, "high", "low")
  } else {
    ifelse(x >= stats::median(x), "high", "low")
  }
}

plot_high_low_vs_distance <- function(df, outfile, gene_lab, xlab, write_spots = FALSE) {
  df <- df[is.finite(df$dist_mm) & df$grp %in% c("high", "low"), , drop = FALSE]
  n_hi <- sum(df$grp == "high")
  n_lo <- sum(df$grp == "low")
  if (nrow(df) < 40 || n_hi < 15 || n_lo < 15) {
    return(list(plotted = FALSE, n_high = n_hi, n_low = n_lo, note = "high/low spot 太少"))
  }
  xmax <- max(df$dist_mm, na.rm = TRUE)
  if (!is.finite(xmax) || xmax <= 0) {
    return(list(plotted = FALSE, n_high = n_hi, n_low = n_lo, note = "距离无效"))
  }
  df$grp <- factor(df$grp, levels = c("high", "low"))
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(list(plotted = FALSE, n_high = n_hi, n_low = n_lo, note = "无 ggplot2"))
  }
  pal <- c(high = "#E31A1C", low = "#377EB8")
  labs_grp <- c(
    high = paste0(gene_lab, "+ cells"),
    low = paste0(gene_lab, "- cells")
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(x = dist_mm, colour = grp)) +
    ggplot2::stat_density(geom = "line", position = "identity",
                          linewidth = 1.2, adjust = 1.15, bounds = c(0, Inf)) +
    ggplot2::scale_colour_manual(values = pal, labels = labs_grp, name = "Cell type") +
    ggplot2::coord_cartesian(xlim = c(0, xmax), expand = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position = c(0.72, 0.84),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    ) +
    ggplot2::labs(
      x = xlab,
      y = "Density"
    )
  ggplot2::ggsave(paste0(outfile, ".pdf"), p, width = 6.2, height = 5.0)
  ggplot2::ggsave(paste0(outfile, ".png"), p, width = 6.2, height = 5.0, dpi = 150)
  if (isTRUE(write_spots)) {
    utils::write.csv(df, paste0(outfile, "_spots.csv"), row.names = FALSE)
  }
  list(plotted = TRUE, n_high = n_hi, n_low = n_lo, note = "")
}

safe_gene_filename <- function(g) {
  gsub("[^A-Za-z0-9._-]+", "_", g)
}

dist_dir <- file.path(out_root, "03_distance_to_Schwann")
dir.create(dist_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(c(
  "比较：远神经肿瘤组织 vs 近神经肿瘤组织。",
  "上调 = 远神经肿瘤更高。FC 只有 >1、1.25。没有 1.5、没有 2、没有 topN、没有下调。",
  "这是指定 32 个病人的综合分析，不覆盖 00_eligible_single_patients / 01_combined。",
  "",
  paste0("近/远统一按这 32 个病人一套标准：近 ≤", near_spots,
         " spot，远 ≥", far_spots, " spot（中间 ", near_spots,
         "–", far_spots, " 不进近/远）。旧版远>4 间隔太小。"),
  "1 spot ≈ 0.1 mm。看 DISTANCE_CUTOFFS.txt 和 DISTANCE_per_patient.csv。",
  "",
  "先打开：INDEX_requested_vs_used.csv",
  "综合上调基因：",
  "  upregulated_far_tumor_vs_near_tumor_FC_gt_1.csv",
  "  upregulated_far_tumor_vs_near_tumor_FC_1.25.csv",
  "",
  "  03_distance_to_Schwann/",
  "  不要打开旧的 genes/（ARTN/NGF 不是这批上调基因，重跑会删掉）",
  "  本次要画的上调基因.txt",
  "  INDEX_FC_1.25_each_gene.csv",
  "  FC_1.25_each_gene/每个上调基因一张图",
  "  神经 marker：SOX10, MPZ, PMP22, S100B, PLP1, NGFR, NCAM1, MBP, L1CAM",
  "  这是空间转录组 RNA，不是蛋白质组。"
), file.path(out_root, "00_请先看这里.txt"))

log_msg("Wang ST dir: ", nerve_dir)
log_msg("Main functions from: ", main_r)
log_msg("Output: ", out_root)
log_msg("Requested patients n=", length(selected_patients), ": ",
        paste(selected_patients, collapse = ","))
log_msg("Near/far cutoffs (selected-patient average analysis): near<=",
        near_spots, " spot, far>=", far_spots, " spot (1 spot ~ 0.1 mm)")

available <- list_patient_ids()
ids <- load_ids()

status_rows <- list()
qc_rows <- list()
dist_rows <- list()
pb_schwann <- NULL
si_schwann <- NULL
mean_schwann_wide <- NULL
sig_cache <- list()

for (pid in selected_patients) {
  has_counts <- pid %in% available
  row <- data.frame(
    patient = pid,
    requested = TRUE,
    has_counts = has_counts,
    classified = FALSE,
    n_tumor_near_schwann = NA_integer_,
    n_tumor_far_schwann = NA_integer_,
    eligible_combined_DE = FALSE,
    folder = "",
    note = if (has_counts) "" else "Robjects/counts 下没有 TNBC*.RDS",
    stringsAsFactors = FALSE
  )
  if (!has_counts) {
    log_msg("TNBC", pid, " 无 counts，跳过")
    status_rows[[as.character(pid)]] <- row
    next
  }
  log_msg("Patient TNBC", pid)
  obj <- tryCatch(classify_one_patient(pid, ids), error = function(e) {
    log_msg("classify failed TNBC", pid, ": ", e$message)
    NULL
  })
  if (is.null(obj)) {
    row$note <- "classify 失败"
    status_rows[[as.character(pid)]] <- row
    next
  }
  row$classified <- TRUE
  obj$spots <- apply_selected_near_far(obj$spots, near_spots, far_spots)
  sp <- obj$spots
  td <- as.numeric(sp$dist_schwann[sp$is_tumor & !sp$is_schwann_high &
                                   is.finite(sp$dist_schwann)])
  dist_rows[[as.character(pid)]] <- data.frame(
    patient = pid,
    n_tumor_to_schwann = length(td),
    mean_dist_spot = if (length(td) > 0) mean(td) else NA_real_,
    median_dist_spot = if (length(td) > 0) stats::median(td) else NA_real_,
    q90_dist_spot = if (length(td) > 0) as.numeric(stats::quantile(td, 0.90, names = FALSE)) else NA_real_,
    max_dist_spot = if (length(td) > 0) max(td) else NA_real_,
    near_spots = near_spots,
    far_spots = far_spots,
    n_tumor_near_schwann = sum(sp$near_schwann),
    n_tumor_far_schwann = sum(sp$far_schwann),
    stringsAsFactors = FALSE
  )
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
    near_spots = near_spots,
    far_spots = far_spots,
    stringsAsFactors = FALSE
  )
  row$n_tumor_near_schwann <- sum(sp$near_schwann)
  row$n_tumor_far_schwann <- sum(sp$far_schwann)

  gkeep <- filter_genes_mat(obj$cnts)
  logm <- obj$logmat[, gkeep, drop = FALSE]
  cnts <- obj$cnts[, gkeep, drop = FALSE]
  tum <- which(sp$is_tumor & !sp$is_schwann_high & is.finite(sp$dist_schwann))
  if (length(tum) >= 20) {
    dist_mm <- as.numeric(sp$dist_schwann[tum]) * spot_pitch_mm
    sig_cache[[as.character(pid)]] <- list(
      patient = pid,
      dist_mm = dist_mm,
      logm = logm[tum, , drop = FALSE]
    )
  }

  n_near_s <- which(sp$near_schwann)
  n_far_s <- which(sp$far_schwann)
  elig <- length(n_near_s) >= 5 && length(n_far_s) >= 10
  row$eligible_combined_DE <- elig
  if (!elig) {
    row$note <- paste0("近(≤", near_spots, " spot)肿瘤=", length(n_near_s),
                       " 远(≥", far_spots, " spot)肿瘤=", length(n_far_s),
                       "（需要近>=5 且 远>=10 个spot）")
    status_rows[[as.character(pid)]] <- row
    rm(obj, cnts, logm)
    gc(verbose = FALSE)
    next
  }

  exs <- export_patient_near_far(pid, "schwann_neighborhood", logm, cnts,
                                 n_near_s, n_far_s, TRUE)
  if (!is.null(exs)) {
    tag <- paste0("TNBC", pid)
    mean_schwann_wide <- bind_gene_col(mean_schwann_wide, exs$mean_near, paste0(tag, "_near"))
    mean_schwann_wide <- bind_gene_col(mean_schwann_wide, exs$mean_far, paste0(tag, "_far"))
  }
  near_sum <- colSums(cnts[n_near_s, , drop = FALSE])
  far_sum <- colSums(cnts[n_far_s, , drop = FALSE])
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
  row$folder <- file.path(expr_dir, "schwann_neighborhood", paste0("TNBC", pid))
  row$note <- "进入综合 远 vs 近"
  status_rows[[as.character(pid)]] <- row
  rm(obj, cnts, logm)
  gc(verbose = FALSE)
}

st <- do.call(rbind, status_rows)
utils::write.csv(st, file.path(out_root, "INDEX_requested_vs_used.csv"), row.names = FALSE)
if (length(qc_rows) > 0) {
  utils::write.csv(do.call(rbind, qc_rows),
                   file.path(out_root, "spot_class_counts.csv"), row.names = FALSE)
}
if (length(dist_rows) > 0) {
  dist_tab <- do.call(rbind, dist_rows)
  utils::write.csv(dist_tab, file.path(out_root, "DISTANCE_per_patient.csv"),
                   row.names = FALSE)
  avg_mean <- mean(dist_tab$mean_dist_spot, na.rm = TRUE)
  avg_med <- mean(dist_tab$median_dist_spot, na.rm = TRUE)
  avg_q90 <- mean(dist_tab$q90_dist_spot, na.rm = TRUE)
  avg_max <- mean(dist_tab$max_dist_spot, na.rm = TRUE)
  writeLines(c(
    paste0("指定病人统一近/远（综合分析同一套）：近 ≤", near_spots,
           " spot，远 ≥", far_spots, " spot。"),
    "中间距离不进近组、也不进远组。1 spot ≈ 0.1 mm。",
    "",
    paste0("这批已分类病人肿瘤→施旺距离平均：mean=",
           signif(avg_mean, 4), "  median=", signif(avg_med, 4),
           "  q90=", signif(avg_q90, 4), "  max=", signif(avg_max, 4),
           " spot。"),
    "近=2、远=10 是按这批病人平均后定的更大间隔（旧版远>4 太近）。",
    "逐病人数字：DISTANCE_per_patient.csv"
  ), file.path(out_root, "DISTANCE_CUTOFFS.txt"))
  log_msg("Selected-patient avg tumor-Schwann dist (spot): mean=",
          signif(avg_mean, 4), " median=", signif(avg_med, 4),
          " q90=", signif(avg_q90, 4), " max=", signif(avg_max, 4))
}
used <- st$patient[st$eligible_combined_DE]
log_msg("Requested ", nrow(st),
        " | with counts ", sum(st$has_counts),
        " | eligible combined ", length(used),
        ": ", paste(used, collapse = ","))

if (!is.null(mean_schwann_wide) && ncol(mean_schwann_wide) > 0) {
  write_gene_mat(mean_schwann_wide,
                 file.path(out_root, "COMBINED_mean_logCPM_near_and_far.csv"))
}
if (!is.null(pb_schwann) && ncol(pb_schwann) > 0) {
  write_gene_mat(pb_schwann,
                 file.path(out_root, "COMBINED_pseudobulk_counts_near_and_far.csv"))
  utils::write.csv(si_schwann, file.path(out_root, "COMBINED_sample_info.csv"),
                   row.names = FALSE)
}

plot_each_fc125_gene <- function(up_tab) {
  gene_dir <- file.path(dist_dir, "FC_1.25_each_gene")
  if (dir.exists(gene_dir)) unlink(gene_dir, recursive = TRUE, force = TRUE)
  # 旧版误画的配体图（ARTN/NGF/CXCL12）在 genes/，必须删掉，避免当成上调基因
  old_ligand <- file.path(dist_dir, "genes")
  if (dir.exists(old_ligand)) {
    unlink(old_ligand, recursive = TRUE, force = TRUE)
    log_msg("已删除旧文件夹 genes/（那是配体 NGF/ARTN，不是这批上调基因）")
  }
  stale <- list.files(dist_dir, full.names = TRUE, include.dirs = TRUE)
  stale <- stale[grepl("^(up_FC_|ligand_|NGF_not|genes)$", basename(stale))]
  if (length(stale) > 0) unlink(stale, recursive = TRUE, force = TRUE)
  dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

  writeLines(c(
    "只画「远神经肿瘤 vs 近神经肿瘤、FC>1.25 上调」的那些基因，每个基因一张图。",
    "不要看旧的 genes/ 文件夹（ARTN、NGF、CXCL12 是文献配体，不是这批上调基因）。",
    "名单：INDEX_FC_1.25_each_gene.csv 和 本次要画的上调基因.txt",
    "施旺 marker：SOX10, MPZ, PMP22, S100B, PLP1, NGFR, NCAM1, MBP, L1CAM。",
    "红 = 该基因高表达，蓝 = 该基因低表达。1 spot = 0.1 mm。",
    paste0("近/远：近 ≤", near_spots, " spot，远 ≥", far_spots, " spot。"),
    "这是空间转录组 RNA，不是蛋白质组。"
  ), file.path(dist_dir, "00_READ_ME.txt"))

  if (is.null(up_tab) || nrow(up_tab) == 0 || length(sig_cache) == 0) {
    log_msg("没有 FC>1.25 上调基因，或没有肿瘤 spot，跳过逐基因距离图")
    return(invisible(NULL))
  }

  xlab <- "Distance to Schwann (mm)"
  idx <- list()
  genes <- unique(as.character(up_tab$gene))
  log_msg("FC>1.25 上调基因共 ", length(genes), " 个，只分析这些：")
  log_msg(paste(genes, collapse = ", "))
  writeLines(
    c(paste0("共 ", length(genes), " 个 FC>1.25 上调基因，逐个画高低表达 vs 施旺距离："),
      genes),
    file.path(dist_dir, "本次要画的上调基因.txt")
  )
  for (g in genes) {
    parts <- list()
    for (nm in names(sig_cache)) {
      ch <- sig_cache[[nm]]
      if (!g %in% colnames(ch$logm)) next
      expr <- as.numeric(ch$logm[, g])
      grp <- label_high_low(expr)
      parts[[nm]] <- data.frame(
        patient = ch$patient,
        dist_mm = ch$dist_mm,
        expr = expr,
        grp = grp,
        stringsAsFactors = FALSE
      )
    }
    row0 <- up_tab[match(g, up_tab$gene), , drop = FALSE]
    rec <- data.frame(
      gene = g,
      log2FC_far_vs_near = if ("log2FC" %in% names(row0)) row0$log2FC[1] else NA_real_,
      FC = if ("FC" %in% names(row0)) row0$FC[1] else NA_real_,
      pvalue = if ("pvalue" %in% names(row0)) row0$pvalue[1] else NA_real_,
      n_high = 0L, n_low = 0L, plotted = FALSE, note = "",
      pdf = "",
      stringsAsFactors = FALSE
    )
    if (length(parts) == 0) {
      rec$note <- "缓存里没有这个基因"
      idx[[g]] <- rec
      next
    }
    d <- do.call(rbind, parts)
    d <- d[!is.na(d$grp), , drop = FALSE]
    out_stub <- file.path(gene_dir, paste0(safe_gene_filename(g), "_high_vs_low_vs_Schwann_distance"))
    res <- plot_high_low_vs_distance(d, out_stub, g, xlab, write_spots = FALSE)
    rec$n_high <- as.integer(res$n_high)
    rec$n_low <- as.integer(res$n_low)
    rec$plotted <- isTRUE(res$plotted)
    rec$note <- res$note
    if (isTRUE(res$plotted)) rec$pdf <- paste0(out_stub, ".pdf")
    idx[[g]] <- rec
  }
  index <- do.call(rbind, idx)
  if ("FC" %in% names(index)) {
    index <- index[order(-index$FC, index$pvalue), , drop = FALSE]
  }
  utils::write.csv(index, file.path(dist_dir, "INDEX_FC_1.25_each_gene.csv"),
                   row.names = FALSE)
  log_msg("逐基因距离图完成：画出 ", sum(index$plotted), " / ", nrow(index))
}

up_fc125_tab <- NULL
if (is.null(pb_schwann) || length(unique(si_schwann$patient)) < 2) {
  log_msg("能进入综合分析的病人不足 2 例。看 INDEX_requested_vs_used.csv")
  writeLines(
    paste0("能进入综合分析的病人不足 2 例（需要近(≤", near_spots,
           " spot)>=5 且 远(≥", far_spots, " spot)>=10 个肿瘤spot）。"),
    flag_no_de
  )
} else {
  log_msg("综合 远神经肿瘤 vs 近神经肿瘤，病人: ",
          paste(unique(si_schwann$patient), collapse = ","))
  rownames(si_schwann) <- si_schwann$sample
  de_sw <- deseq2_pseudobulk(round(pmax(pb_schwann, 0)), si_schwann,
                             "selected32_far_tumor_vs_near_tumor")
  if (nrow(de_sw) == 0) {
    logm <- log2(sweep(t(pb_schwann), 1, pmax(colSums(pb_schwann), 1), "/") * 1e6 + 1)
    de_sw <- limma_two_group(logm, si_schwann$group, si_schwann$patient, "selected32_pb")
  }
  emit_de_tables("selected32_far_tumor_vs_near_tumor", de_sw,
                 out_dir = out_root, also_top_level = TRUE)
  up_fc125_tab <- select_up(de_sw, 1.25)
}

plot_each_fc125_gene(up_fc125_tab)

log_msg("先打开: ", file.path(out_root, "00_请先看这里.txt"))
log_msg("近/远: ≤", near_spots, " / ≥", far_spots, " spot  ",
        file.path(out_root, "DISTANCE_CUTOFFS.txt"))
log_msg("综合上调: ", out_root)
log_msg("FC>1.25 逐基因距离图: ", file.path(dist_dir, "FC_1.25_each_gene"))
