#!/usr/bin/env Rscript
# =============================================================================
# 单独运行：指定病人综合  远神经肿瘤组织 vs 近神经肿瘤组织
# 上调只做 FC>1、1.25。不做 FC=1.5、不做 FC=2、不做 topN、不做下调、不做 GO。
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

# 本脚本只要上调 FC>1、1.25
p_cutoff <- 0.05
fc_cutoffs_up <- c("FC_gt_1" = 1, "FC_1.25" = 1.25)
# Wang ST spot 中心距约 100 μm；距离用 spot 网格换算成 mm
spot_pitch_mm <- 0.1
# 高低表达要画距离图的基因（NGF 本阵列常没有，有则画）
highlow_genes <- c(
  "NGF", "BDNF", "NTF3", "NTF4", "GDNF", "ARTN", "NRTN",
  "NRG1", "CXCL12", "VEGFA", "MDK", "PTN", "TGFB1", "WNT5A", "NGFR"
)
# 距离参照：施旺细胞签名（与全队列脚本同一套 marker）
schwann_distance_markers <- schwann_genes

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

plot_high_low_vs_distance <- function(df, outfile, gene_lab, xlab) {
  df <- df[is.finite(df$dist_mm) & df$grp %in% c("high", "low"), , drop = FALSE]
  if (nrow(df) < 40) return(invisible(NULL))
  n_hi <- sum(df$grp == "high")
  n_lo <- sum(df$grp == "low")
  if (n_hi < 15 || n_lo < 15) return(invisible(NULL))
  xmax <- max(df$dist_mm, na.rm = TRUE)
  if (!is.finite(xmax) || xmax <= 0) return(invisible(NULL))
  df$grp <- factor(df$grp, levels = c("high", "low"))
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
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
  utils::write.csv(df, paste0(outfile, "_spots.csv"), row.names = FALSE)
  invisible(p)
}

score_gene_set <- function(mat, genes) {
  g <- intersect(unique(genes), colnames(mat))
  if (length(g) < 3) return(rep(NA_real_, nrow(mat)))
  z <- scale(mat[, g, drop = FALSE])
  z[is.na(z)] <- 0
  as.numeric(rowMeans(z))
}

dist_dir <- file.path(out_root, "03_distance_to_Schwann")
dir.create(dist_dir, recursive = TRUE, showWarnings = FALSE)

writeLines(c(
  "比较：远神经肿瘤组织 vs 近神经肿瘤组织。",
  "上调 = 远神经肿瘤更高。FC 只有 >1、1.25。没有 1.5、没有 2、没有 topN、没有下调。",
  "这是指定 32 个病人的综合分析，不覆盖 00_eligible_single_patients / 01_combined。",
  "",
  "先打开：INDEX_requested_vs_used.csv",
  "综合上调基因：",
  "  upregulated_far_tumor_vs_near_tumor_FC_gt_1.csv",
  "  upregulated_far_tumor_vs_near_tumor_FC_1.25.csv",
  "",
  "高低表达 vs 施旺细胞距离（类似 NGF+ / NGF- 密度图）：",
  "  03_distance_to_Schwann/",
  "  神经 marker：SOX10, MPZ, PMP22, S100B, PLP1, NGFR, NCAM1, MBP, L1CAM",
  "  主图：up_FC_gt_1_high_vs_low_vs_Schwann_distance.pdf"
), file.path(out_root, "00_请先看这里.txt"))

log_msg("Wang ST dir: ", nerve_dir)
log_msg("Main functions from: ", main_r)
log_msg("Output: ", out_root)
log_msg("Requested patients n=", length(selected_patients), ": ",
        paste(selected_patients, collapse = ","))

available <- list_patient_ids()
ids <- load_ids()

status_rows <- list()
qc_rows <- list()
pb_schwann <- NULL
si_schwann <- NULL
mean_schwann_wide <- NULL
hl_gene_rows <- list()
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
    dist_pni_mm <- as.numeric(sp$dist_nerve[tum]) * spot_pitch_mm
    g_here <- intersect(highlow_genes, colnames(logm))
    if (length(g_here) > 0) {
      for (g in g_here) {
        expr <- as.numeric(logm[tum, g])
        hl_gene_rows[[paste(pid, g, sep = "_")]] <- data.frame(
          patient = pid, gene = g, dist_mm = dist_mm,
          dist_pni_mm = dist_pni_mm, expr = expr,
          grp = label_high_low(expr),
          stringsAsFactors = FALSE
        )
      }
    }
    sig_cache[[as.character(pid)]] <- list(
      patient = pid,
      dist_mm = dist_mm,
      dist_pni_mm = dist_pni_mm,
      logm = logm[tum, , drop = FALSE]
    )
  }

  n_near_s <- which(sp$near_schwann)
  n_far_s <- which(sp$far_schwann)
  elig <- length(n_near_s) >= 5 && length(n_far_s) >= 10
  row$eligible_combined_DE <- elig
  if (!elig) {
    row$note <- paste0("近神经肿瘤 spot=", length(n_near_s),
                       " 远神经肿瘤 spot=", length(n_far_s),
                       "（需要近>=5 且 远>=10）")
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

emit_highlow_plots <- function(up_fc1, up_fc125) {
  writeLines(c(
    "Tumor spots split into gene-high vs gene-low, then density vs distance to Schwann-high spots.",
    "Schwann markers: SOX10, MPZ, PMP22, S100B, PLP1, NGFR, NCAM1, MBP, L1CAM.",
    "1 spot = 0.1 mm. Red = high / positive, blue = low / negative.",
    "NGF is often absent on this array; if missing it is listed below and NGFR / ligands are plotted instead."
  ), file.path(dist_dir, "00_READ_ME.txt"))

  missing_ngf <- !("NGF" %in% unique(if (length(hl_gene_rows) > 0) {
    do.call(rbind, hl_gene_rows)$gene
  } else character()))
  if (isTRUE(missing_ngf)) {
    writeLines("NGF was not detected on this array. Plotted detected ligands and the upregulated-gene signature instead.",
               file.path(dist_dir, "NGF_not_on_array.txt"))
    log_msg("NGF 本阵列未检出，改画可检出配体 + 上调基因签名 vs 施旺距离")
  }

  if (length(hl_gene_rows) > 0) {
    hl <- do.call(rbind, hl_gene_rows)
    hl <- hl[!is.na(hl$grp), , drop = FALSE]
    utils::write.csv(hl, file.path(dist_dir, "ligand_spots_distance_to_Schwann.csv"),
                     row.names = FALSE)
    for (g in unique(hl$gene)) {
      d <- hl[hl$gene == g, , drop = FALSE]
      plot_high_low_vs_distance(
        d, file.path(dist_dir, "genes", paste0(g, "_high_vs_low_vs_Schwann_distance")),
        g, "Distance to Schwann (mm)"
      )
    }
    log_msg("Per-gene distance plots: ", paste(unique(hl$gene), collapse = ","))
  }

  plot_signature <- function(up_genes, tag, lab) {
    if (length(sig_cache) == 0 || length(up_genes) < 3) return(invisible(NULL))
    parts <- list()
    for (nm in names(sig_cache)) {
      ch <- sig_cache[[nm]]
      sc <- score_gene_set(ch$logm, up_genes)
      if (all(!is.finite(sc))) next
      parts[[nm]] <- data.frame(
        patient = ch$patient,
        dist_mm = ch$dist_mm,
        expr = sc,
        grp = label_high_low(sc),
        stringsAsFactors = FALSE
      )
    }
    if (length(parts) == 0) return(invisible(NULL))
    d <- do.call(rbind, parts)
    plot_high_low_vs_distance(
      d,
      file.path(dist_dir, paste0(tag, "_high_vs_low_vs_Schwann_distance")),
      lab,
      "Distance to Schwann (mm)"
    )
    log_msg("Signature distance plot ", tag, " n_spots=", nrow(d),
            " n_genes=", length(intersect(up_genes, colnames(sig_cache[[1]]$logm))))
  }

  if (length(up_fc1) >= 3) plot_signature(up_fc1, "up_FC_gt_1", "Up FC>1")
  if (length(up_fc125) >= 3) plot_signature(up_fc125, "up_FC_1.25", "Up FC>=1.25")
}

up_fc1 <- character()
up_fc125 <- character()
if (is.null(pb_schwann) || length(unique(si_schwann$patient)) < 2) {
  log_msg("能进入综合分析的病人不足 2 例。看 INDEX_requested_vs_used.csv")
  writeLines(
    "能进入综合分析的病人不足 2 例（需要近神经肿瘤 spot>=5 且 远>=10）。",
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
  up_fc1 <- select_up(de_sw, 1)$gene
  up_fc125 <- select_up(de_sw, 1.25)$gene
}

emit_highlow_plots(up_fc1, up_fc125)

log_msg("先打开: ", file.path(out_root, "00_请先看这里.txt"))
log_msg("综合上调: ", out_root)
log_msg("距离图: ", dist_dir)
