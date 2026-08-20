#!/usr/bin/env Rscript
# =============================================================================
# 设计 5 / 6（也可单独跑；主流程 TG_RNAseq_pipeline.R 已包含这两组）
#   5) mean(TG_sh1, TG_sh5) vs NTC_rep0
#   6) mean(TG_sh1, TG_sh5) vs NTC_rep1
# 下调网格：
#   FoldChange: treat/control <= 1 / 0.8 / 0.6 / 0.5
#   TopRank: 下调前 50 / 75 / 100 / 150 / 200 / 250 / 300
#   火山图、热图、GO、通路、KEGG、GSEA
#
# 用法：
#   setwd("E:/R/TG_BRCA/TG")
#   source("TG_RNAseq_pipeline.R")                 # 六种设计（含下面两组）
#   source("TG_RNAseq_TGsh_mean_vs_NTC_reps.R")    # 若主流程已跑过则跳过
# =============================================================================

load_pipeline_functions_only <- function() {
  pipe <- "TG_RNAseq_pipeline.R"
  if (!file.exists(pipe)) {
    alt <- file.path(getwd(), "TG_RNAseq_pipeline.R")
    if (file.exists(alt)) pipe <- alt
  }
  if (!file.exists(pipe)) stop("找不到 TG_RNAseq_pipeline.R")
  lines <- readLines(pipe, warn = FALSE)
  main_at <- grep("^# 10\\. 主流程", lines)[1]
  if (is.na(main_at) || main_at < 2) stop("无法从原脚本切出函数定义")
  eval(parse(text = lines[seq_len(main_at - 1)]), envir = .GlobalEnv)
}

if (!exists("analyze_one_comparison", mode = "function")) {
  load_pipeline_functions_only()
}

if (!exists("mean_kd_vs_one_ntc_de", mode = "function")) {
  mean_kd_vs_one_ntc_de <- function(log_mat, sample_info, ntc_id, comp_name) {
    sh1 <- find_sample(sample_info, "TG_sh1")
    sh5 <- find_sample(sample_info, "TG_sh5")
    ntc <- find_sample(sample_info, "NTC", ntc_id)
    if (is.na(sh1) || is.na(sh5) || is.na(ntc)) return(NULL)
    if (!all(c(sh1, sh5, ntc) %in% colnames(log_mat))) return(NULL)
    log_msg(comp_name, " : mean(", sh1, ", ", sh5, ") vs ", ntc, " (KD mean vs one NTC, FC only)")
    sh_mean <- (log_mat[, sh1] + log_mat[, sh5]) / 2
    ntc_val <- log_mat[, ntc]
    data.frame(
      gene = rownames(log_mat),
      log2FC = as.numeric(sh_mean - ntc_val),
      AveExpr = as.numeric((sh_mean + ntc_val) / 2),
      pvalue = NA_real_,
      padj = NA_real_,
      treat_sample = paste0("mean(", sh1, ",", sh5, ")"),
      ntc_sample = ntc,
      stringsAsFactors = FALSE
    )
  }
}

already_in_main <- exists("de_list", envir = .GlobalEnv) &&
  all(c("TGsh_mean_vs_NTC_rep0", "TGsh_mean_vs_NTC_rep1") %in% names(get("de_list", envir = .GlobalEnv)))

if (isTRUE(already_in_main)) {
  log_msg("Designs 5-6 already ran in TG_RNAseq_pipeline.R; skip extra script")
} else {
  prepare_extra_expression <- function() {
    if (exists("log_mat", envir = .GlobalEnv) &&
        exists("si", envir = .GlobalEnv) &&
        exists("norm", envir = .GlobalEnv)) {
      log_msg("Reuse normalized matrix already in memory from TG_RNAseq_pipeline.R")
      return(list(
        log_mat = get("log_mat", envir = .GlobalEnv),
        heat_mat = get("norm", envir = .GlobalEnv)$heat_mat,
        sample_info = get("si", envir = .GlobalEnv)
      ))
    }

    csv_mat <- file.path(result_dir, "normalized_log_matrix.csv")
    csv_si <- file.path(log_dir, "sample_info.csv")
    if (file.exists(csv_mat) && file.exists(csv_si)) {
      log_msg("Reuse saved matrix: ", csv_mat)
      df <- utils::read.csv(csv_mat, check.names = FALSE, stringsAsFactors = FALSE)
      genes <- df[[1]]
      mat <- as.matrix(df[, -1, drop = FALSE])
      storage.mode(mat) <- "double"
      rownames(mat) <- genes
      si_df <- utils::read.csv(csv_si, check.names = FALSE, stringsAsFactors = FALSE)
      si_df <- si_df[match(colnames(mat), si_df$sample), ]
      if (!"ntc_id" %in% names(si_df)) si_df <- add_ntc_ids(si_df)
      return(list(log_mat = mat, heat_mat = mat, sample_info = si_df))
    }

    log_msg("Load and normalize expression for extra KD-mean vs single NTC comparisons")
    expr <- load_expression(project_dir)
    expr$sample_info <- add_ntc_ids(expr$sample_info)
    value_type <- detect_value_type(expr$mat)
    filt <- filter_low_expression(expr$mat, expr$sample_info, value_type)
    norm_obj <- normalize_expression(filt, expr$sample_info, value_type)
    si_df <- expr$sample_info[match(colnames(norm_obj$log_mat), expr$sample_info$sample), ]
    si_df <- add_ntc_ids(si_df)
    list(log_mat = norm_obj$log_mat, heat_mat = norm_obj$heat_mat, sample_info = si_df)
  }

  log_msg("Extra comparisons: mean(TG_sh1, TG_sh5) vs each NTC (downregulation grid)")
  prep <- prepare_extra_expression()
  log_mat_extra <- prep$log_mat
  heat_mat_extra <- prep$heat_mat
  si_extra <- prep$sample_info

  extra_list <- list(
    TGsh_mean_vs_NTC_rep0 = mean_kd_vs_one_ntc_de(
      log_mat_extra, si_extra, "NTC_rep0", "TGsh_mean_vs_NTC_rep0"
    ),
    TGsh_mean_vs_NTC_rep1 = mean_kd_vs_one_ntc_de(
      log_mat_extra, si_extra, "NTC_rep1", "TGsh_mean_vs_NTC_rep1"
    )
  )
  extra_list <- extra_list[!vapply(extra_list, is.null, logical(1))]
  if (length(extra_list) == 0) {
    stop("无法构建 mean(TG_sh1, TG_sh5) vs NTC_rep0 / NTC_rep1，请检查样本名")
  }

  for (nm in names(extra_list)) {
    have_p <- any(!is.na(extra_list[[nm]]$padj))
    tryCatch(
      analyze_one_comparison(
        nm, extra_list[[nm]], extra_list[[nm]],
        heat_mat_extra, si_extra, have_p, extra_list[[nm]]
      ),
      error = function(e) log_msg("ERROR in extra comparison ", nm, ": ", e$message)
    )
  }

  log_msg("Extra comparisons done: ", paste(names(extra_list), collapse = ", "))
  log_msg("Folders: ", file.path(result_dir, "TGsh_mean_vs_NTC_rep0"), " ; ",
          file.path(result_dir, "TGsh_mean_vs_NTC_rep1"))
}
