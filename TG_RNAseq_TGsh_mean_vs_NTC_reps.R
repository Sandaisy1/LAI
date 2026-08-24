#!/usr/bin/env Rscript
# =============================================================================
# 额外比较（不修改 TG_RNAseq_pipeline.R 的比较 1–4）
#   5) mean(TG_sh1, TG_sh5) vs NTC_rep0
#   6) mean(TG_sh1, TG_sh5) vs NTC_rep1
#   7) mean(common(TG_sh1 & TG_sh5)) vs mean(common(NTC_rep0 & NTC_rep1))
#      先找两 KD 都检测到的基因、两 NTC 都检测到的基因，取交集后再均值比较
# 分层不用 p，只按 FoldChange 与上调排名。
# 对列出 GO 出表、火山图、热图、气泡图与 BP/CC/MF 柱状图。
#
# 用法：
#   setwd("E:/R/TG_BRCA/TG")
#   source("TG_RNAseq_pipeline.R")                 # 比较 1–4
#   source("TG_RNAseq_TGsh_mean_vs_NTC_reps.R")    # 比较 5–7
# =============================================================================

if (!exists("analyze_one_comparison", mode = "function") ||
    !exists("mean_kd_vs_one_ntc_de", mode = "function") ||
    !exists("mean_common_kd_vs_mean_common_ntc_de", mode = "function")) {
  options(tg.rnaseq.functions_only = TRUE)
  pipe <- "TG_RNAseq_pipeline.R"
  if (!file.exists(pipe)) pipe <- file.path(getwd(), "TG_RNAseq_pipeline.R")
  if (!file.exists(pipe)) stop("找不到 TG_RNAseq_pipeline.R")
  sys.source(pipe, envir = .GlobalEnv, keep.source = TRUE)
  options(tg.rnaseq.functions_only = FALSE)
}

prepare_extra_expression <- function() {
  if (exists("log_mat", envir = .GlobalEnv) &&
      exists("si", envir = .GlobalEnv) &&
      exists("norm", envir = .GlobalEnv) &&
      exists("go_tab", envir = .GlobalEnv) &&
      exists("go_sets", envir = .GlobalEnv)) {
    log_msg("Reuse normalized matrix and GO sets from TG_RNAseq_pipeline.R")
    return(list(
      log_mat = get("log_mat", envir = .GlobalEnv),
      heat_mat = get("norm", envir = .GlobalEnv)$heat_mat,
      sample_info = get("si", envir = .GlobalEnv),
      go_tab = get("go_tab", envir = .GlobalEnv),
      go_sets = get("go_sets", envir = .GlobalEnv),
      pathway_genes = get("pathway_genes", envir = .GlobalEnv)
    ))
  }

  go_path <- find_custom_go_file(project_dir)
  go_tab <- parse_custom_go_file(go_path)
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
    go_sets <- map_go_to_symbols(go_tab$go_id, rownames(mat))
    names(go_sets) <- go_tab$go_id
    return(list(
      log_mat = mat, heat_mat = mat, sample_info = si_df,
      go_tab = go_tab, go_sets = go_sets,
      pathway_genes = all_pathway_genes(go_sets)
    ))
  }

  log_msg("Load Cuffdiff FPKM for extra comparisons 5-7")
  expr <- load_expression(project_dir)
  expr$sample_info <- add_ntc_ids(expr$sample_info)
  value_type <- infer_value_type(expr)
  filt <- filter_low_expression(expr$mat, expr$sample_info, value_type)
  norm_obj <- normalize_expression(filt, expr$sample_info, value_type)
  si_df <- expr$sample_info[match(colnames(norm_obj$log_mat), expr$sample_info$sample), ]
  si_df <- add_ntc_ids(si_df)
  go_sets <- map_go_to_symbols(go_tab$go_id, rownames(norm_obj$log_mat))
  names(go_sets) <- go_tab$go_id
  list(
    log_mat = norm_obj$log_mat, heat_mat = norm_obj$heat_mat, sample_info = si_df,
    go_tab = go_tab, go_sets = go_sets,
    pathway_genes = all_pathway_genes(go_sets)
  )
}

log_msg("Extra comparisons 5-7: KD-mean vs each NTC, plus common-gene means")
if (isTRUE(getOption("tg.rnaseq.restyle_only", FALSE))) {
  restyle_ora_bubbles()
} else {
prep <- prepare_extra_expression()

extra_list <- list(
  TGsh_mean_vs_NTC_rep0 = mean_kd_vs_one_ntc_de(
    prep$log_mat, prep$sample_info, "NTC_rep0", "TGsh_mean_vs_NTC_rep0"
  ),
  TGsh_mean_vs_NTC_rep1 = mean_kd_vs_one_ntc_de(
    prep$log_mat, prep$sample_info, "NTC_rep1", "TGsh_mean_vs_NTC_rep1"
  ),
  mean_common_TGsh_vs_mean_common_NTC = mean_common_kd_vs_mean_common_ntc_de(
    prep$log_mat, prep$sample_info
  )
)
extra_list <- extra_list[!vapply(extra_list, is.null, logical(1))]
if (length(extra_list) == 0) {
  stop("无法构建比较 5–7，请检查样本名 TG_sh1 / TG_sh5 / NTC_rep0 / NTC_rep1")
}

for (nm in names(extra_list)) {
  if (identical(nm, "mean_common_TGsh_vs_mean_common_NTC")) {
    tryCatch(
      write_common_mean_gene_sets(extra_list[[nm]], file.path(result_dir, nm)),
      error = function(e) log_msg("write common gene sets failed: ", e$message)
    )
  }
  tryCatch(
    analyze_one_comparison(
      nm, extra_list[[nm]], extra_list[[nm]],
      prep$heat_mat, prep$sample_info, prep$go_tab, prep$go_sets, prep$pathway_genes
    ),
    error = function(e) log_msg("ERROR in extra comparison ", nm, ": ", e$message)
  )
}

log_msg("Extra comparisons done: ", paste(names(extra_list), collapse = ", "))
}
