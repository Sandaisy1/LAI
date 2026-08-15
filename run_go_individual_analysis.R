# ★★★ 已修改：把本文件和主脚本放在同一目录（E:/R/BRCA）★★★
# 在已经读完数据、设好 go_to_run 的 R 会话里运行：
#   source("run_go_individual_analysis.R", encoding = "UTF-8")
#   run_go_individual_analysis()
#
# 本文件只定义函数，不会自动开始分析。

run_go_individual_analysis <- function(go_ids = NULL) {
  if (is.null(go_ids)) {
    if (!exists("go_to_run", inherits = TRUE)) {
      stop("请先设置 go_to_run，例如：go_to_run <- go_list")
    }
    go_ids <- go_to_run
  }
  go_to_run <- unique(as.character(go_ids))
  if (length(go_to_run) == 0) stop("go_to_run 是空的，没有要分析的 GO")
  if (!exists("expr", inherits = TRUE) || is.null(expr)) {
    stop("还没有表达矩阵 expr。请先从脚本开头 Source 到数据预处理结束，再调用本函数。")
  }

## 6. 逐个 GO 取基因（绝不合并） ------------------------------
# ★★★ 已修改：按 go_to_run 逐个取基因，不是合并成一个基因集 ★★★
message("本次将单独分析 ", length(go_to_run), " 个 GO：", paste(go_to_run, collapse = ", "))
go_gene_list <- lapply(go_to_run, get_go_genes)
names(go_gene_list) <- go_to_run

go_set_summary <- rbindlist(lapply(go_to_run, function(g) {
  dt <- go_gene_list[[g]]
  in_expr <- intersect(unique(dt$SYMBOL), rownames(expr))
  data.table(
    GO = g,
    GO_name = unname(go_name_map[g]),
    n_annotated = uniqueN(dt$SYMBOL),
    n_in_expression = length(in_expr),
    genes_in_expression = paste(in_expr, collapse = ";")
  )
}), fill = TRUE)
fwrite(go_set_summary, file.path(out_dir, "00_GO_gene_sets.csv"))

## 7. 每个 GO 单独分析 ----------------------------------------
all_scores <- list()
all_clin   <- list()
all_surv   <- list()
all_neg    <- list()
all_prot   <- list()

surv_endpoints <- list(
  OS  = c("OS",  "OS.time"),
  DSS = c("DSS", "DSS.time"),
  PFI = c("PFI", "PFI.time"),
  DFI = c("DFI", "DFI.time")
)

for (go_id in go_to_run) {
  go_title <- unname(go_name_map[go_id])
  message("\n========== ", go_id, " | ", go_title, " ==========")

  go_dir <- file.path(out_dir, "per_GO", paste0(safe_name(go_id), "_", safe_name(go_title)))
  dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)

  genes_here <- unique(go_gene_list[[go_id]]$SYMBOL)
  fwrite(
    data.table(GO = go_id, GO_name = go_title, gene = genes_here,
               in_expression = genes_here %in% rownames(expr)),
    file.path(go_dir, "genes.csv")
  )

  # ★★★ 已修改开始：变量名必须用 go_score，禁止再用 score ★★★
  # 旧代码（会报 closure 不可取子集，请勿运行）：
  #   score <- pathway_zmean(expr, genes_here)
  #   clin_use <- clinical_data[sample_std %in% names(score)]
  #   sc_clin <- score[clin_use$sample_std]
  go_score <- tryCatch(pathway_zmean(expr, genes_here), error = function(e) {
    message("  通路打分失败：", conditionMessage(e))
    NULL
  })
  if (!is.numeric(go_score) || length(go_score) == 0) {
    message("  可用基因 < ", min_pathway_genes, " 或打分失败，跳过该通路")
    next
  }
  n_used <- attr(go_score, "n_genes")
  used_genes <- attr(go_score, "genes")
  message("  通路基因用于打分：", n_used)
  fwrite(data.table(sample = names(go_score), pathway_score = as.numeric(go_score)),
         file.path(go_dir, "pathway_score.csv"))
  all_scores[[go_id]] <- go_score
  score <- go_score  # ★★★ 已修改：覆盖函数 score，避免后面旧代码 score[...] 报 closure

  # ---- 7.1 临床相关性（仅本通路分数） ----
  clin_use <- clinical_data[sample_std %in% names(go_score)]
  setkey(clin_use, sample_std)
  sc_clin <- go_score[clin_use$sample_std]
  # ★★★ 已修改结束 ★★★
  clin_rows <- lapply(clin_features, function(ft) {
    if (!ft %in% names(clin_use)) return(NULL)
    assoc_clinical_feature(sc_clin, clin_use[[ft]], ft)
  })
  clin_tab <- rbindlist(clin_rows, fill = TRUE)
  if (nrow(clin_tab) > 0) {
    clin_tab[, `:=`(GO = go_id, GO_name = go_title, fdr = p.adjust(pvalue, method = "BH"))]
    setcolorder(clin_tab, c("GO", "GO_name"))
    fwrite(clin_tab, file.path(go_dir, "clinical_association.csv"))
    all_clin[[go_id]] <- clin_tab

    # ★★★ 已修改开始：临床箱线图 ★★★
    # 旧输出 `null device 1` 只是 dev.off() 关闭 PDF，不是分析结果，也看不出画了几张图。
    sig_cat <- clin_tab[type == "categorical" & pvalue < 0.05][order(pvalue)]
    if (nrow(sig_cat) == 0) {
      message("  本通路无显著分类临床变量（p < 0.05），不画箱线图")
    } else {
      pdf_clin <- file.path(go_dir, "clinical_boxplots.pdf")
      n_box <- 0L
      pdf(pdf_clin, width = 7, height = 5)
      for (ft in head(sig_cat$feature, 6)) {
        plot_df <- data.frame(
          score = as.numeric(sc_clin),
          group = as.character(clin_use[[ft]]),
          stringsAsFactors = FALSE
        )
        plot_df <- plot_df[is.finite(plot_df$score) & !is.na(plot_df$group) & plot_df$group != "", ]
        tab <- table(plot_df$group)
        plot_df <- plot_df[plot_df$group %in% names(tab)[tab >= 3], ]
        if (length(unique(plot_df$group)) < 2) next
        ok <- tryCatch({
          p <- ggboxplot(plot_df, x = "group", y = "score", fill = "group", outlier.shape = NA) +
            stat_compare_means() +
            labs(title = paste0(go_id, "\n", go_title),
                 subtitle = paste0(ft, "  (p=", signif(sig_cat$pvalue[sig_cat$feature == ft][1], 3), ")"),
                 x = NULL, y = "Pathway score") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
          print(p)
          TRUE
        }, error = function(e) {
          message("  箱线图失败（", ft, "）：", conditionMessage(e))
          FALSE
        })
        if (isTRUE(ok)) n_box <- n_box + 1L
      }
      invisible(dev.off())
      if (n_box == 0L) {
        if (file.exists(pdf_clin)) file.remove(pdf_clin)
        message("  箱线图均被跳过（分组不足），未保存 PDF")
      } else {
        message("  已保存临床箱线图 ", n_box, " 张：", pdf_clin)
      }
    }
    # ★★★ 已修改结束 ★★★
  }

  # ---- 7.2 生存分析（仅本通路分数） ----
  # ★★★ 已修改开始：7.2 生存分析，score 改成 go_score ★★★
  # 旧代码（会报 closure 不可取子集，请勿运行）：
  #   surv_use <- survival_data[sample_std %in% names(score)]
  #   sc_surv <- score[surv_use$sample_std]
  if (!exists("go_score") || !is.numeric(go_score)) {
    stop("go_score 不是数值。请先运行通路打分段（搜：变量名必须用 go_score），不要运行 score[...]")
  }
  surv_use <- survival_data[sample_std %in% names(go_score)]
  sc_surv <- go_score[surv_use$sample_std]
  # ★★★ 已修改结束 ★★★
  surv_rows <- list()

  for (ep in names(surv_endpoints)) {
    ev_col <- surv_endpoints[[ep]][1]
    tm_col <- surv_endpoints[[ep]][2]
    if (!all(c(ev_col, tm_col) %in% names(surv_use))) next

    d <- data.frame(
      sample = surv_use$sample_std,
      time = as.numeric(surv_use[[tm_col]]),
      event = as.numeric(surv_use[[ev_col]]),
      score = as.numeric(sc_surv),
      stringsAsFactors = FALSE
    )
    d <- d[is.finite(d$time) & d$time > 0 & d$event %in% c(0, 1) & is.finite(d$score), ]
    if (nrow(d) < 20 || sum(d$event) < 5) next

    d$group <- ifelse(d$score >= median(d$score, na.rm = TRUE), "High", "Low")
    d$group <- factor(d$group, levels = c("Low", "High"))

    cox_cont <- tryCatch(coxph(Surv(time, event) ~ score, data = d), error = function(e) NULL)
    cox_grp  <- tryCatch(coxph(Surv(time, event) ~ group, data = d), error = function(e) NULL)
    fit_km   <- tryCatch(survfit(Surv(time, event) ~ group, data = d), error = function(e) NULL)

    if (!is.null(cox_cont)) {
      s1 <- summary(cox_cont)
      surv_rows[[paste0(ep, "_cont")]] <- data.table(
        GO = go_id, GO_name = go_title, endpoint = ep, model = "continuous_score",
        n = nrow(d), events = sum(d$event),
        HR = s1$conf.int[1, 1],
        HR_low = s1$conf.int[1, 3],
        HR_high = s1$conf.int[1, 4],
        pvalue = s1$coefficients[1, "Pr(>|z|)"]
      )
    }
    if (!is.null(cox_grp)) {
      s2 <- summary(cox_grp)
      surv_rows[[paste0(ep, "_grp")]] <- data.table(
        GO = go_id, GO_name = go_title, endpoint = ep, model = "High_vs_Low",
        n = nrow(d), events = sum(d$event),
        HR = s2$conf.int[1, 1],
        HR_low = s2$conf.int[1, 3],
        HR_high = s2$conf.int[1, 4],
        pvalue = s2$coefficients[1, "Pr(>|z|)"]
      )
    }

    if (!is.null(fit_km) && ep == "OS") {
      tryCatch({
        d$time_month <- d$time / 30.44
        fit_m <- survfit(Surv(time_month, event) ~ group, data = d)
        pdf(file.path(go_dir, "KM_OS.pdf"), width = 7, height = 6)
        print(ggsurvplot(
          fit_m, data = d, pval = TRUE, risk.table = TRUE,
          legend.title = "Pathway score",
          legend.labs = c("Low", "High"),
          xlab = "Time (months)", ylab = "Overall survival",
          title = paste0(go_id, "\n", go_title),
          ggtheme = theme_bw()
        ))
        dev.off()
      }, error = function(e) {
        if (dev.cur() > 1) dev.off()
        message("  KM 作图失败：", conditionMessage(e))
      })
    }
  }

  surv_tab <- rbindlist(surv_rows, fill = TRUE)
  if (nrow(surv_tab) > 0) {
    fwrite(surv_tab, file.path(go_dir, "survival_cox.csv"))
    all_surv[[go_id]] <- surv_tab
  }

  # ---- 7.3 全基因组负相关基因（仅对本通路分数） ----
  cor_tab <- spearman_vs_score(expr, go_score)  # ★★★ 已修改：传入 go_score
  if (nrow(cor_tab) > 0) {
    setnames(cor_tab, "feature", "gene")
    cor_tab[, `:=`(
      GO = go_id,
      GO_name = go_title,
      in_this_GO = gene %in% used_genes
    )]
    fwrite(cor_tab[order(-spearman_r)], file.path(go_dir, "all_gene_correlation.csv"))

    neg_tab <- cor_tab[spearman_r < neg_r_cutoff & fdr < neg_fdr_cutoff]
    setorder(neg_tab, spearman_r)
    fwrite(neg_tab, file.path(go_dir, "negative_correlated_genes.csv"))
    fwrite(neg_tab[spearman_r <= strict_r_cutoff],
           file.path(go_dir, "negative_correlated_genes_strict_r0.3.csv"))
    all_neg[[go_id]] <- neg_tab
    message("  显著负相关基因：", nrow(neg_tab),
            " ；|r|>=0.3：", nrow(neg_tab[spearman_r <= strict_r_cutoff]))

    vol_df <- as.data.frame(cor_tab)
    rownames(vol_df) <- vol_df$gene
    key_neg <- head(neg_tab$gene, 15)
    tryCatch({
      pdf(file.path(go_dir, "volcano_gene_correlation.pdf"), width = 9, height = 7)
      print(EnhancedVolcano(
        vol_df,
        lab = vol_df$gene,
        x = "spearman_r",
        y = "fdr",
        xlab = "Spearman r (gene vs this GO score)",
        ylab = "-Log10 FDR",
        title = paste(go_id, go_title),
        pCutoff = neg_fdr_cutoff,
        FCcutoff = 0.3,
        xlim = c(-1, 1),
        selectLab = key_neg,
        drawConnectors = TRUE
      ))
      dev.off()
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      message("  火山图失败：", conditionMessage(e))
    })

    if (nrow(neg_tab) > 0) {
      topn <- head(neg_tab, 30)
      pbar <- ggplot(topn, aes(x = reorder(gene, spearman_r), y = spearman_r)) +
        geom_col(fill = "#3C5488") +
        coord_flip() +
        labs(title = paste0("Top negative genes: ", go_id),
             x = NULL, y = "Spearman r") +
        theme_bw()
      ggsave(file.path(go_dir, "top_negative_genes.pdf"), pbar, width = 7, height = 8)
    }
  }

  # ---- 7.4 蛋白水平相关性（补充，仍按本通路分数） ----
  if (!is.null(prot_mat)) {
    prot_common <- intersect(colnames(prot_mat), names(go_score))  # ★★★ 已修改：go_score
    if (length(prot_common) >= 20) {
      prot_cor <- spearman_vs_score(prot_mat, go_score)
      if (nrow(prot_cor) > 0) {
        setnames(prot_cor, "feature", "protein")
        prot_cor[, `:=`(GO = go_id, GO_name = go_title)]
        fwrite(prot_cor[order(spearman_r)], file.path(go_dir, "protein_correlation.csv"))
        all_prot[[go_id]] <- prot_cor[spearman_r < 0 & fdr < neg_fdr_cutoff][order(spearman_r)]
      }
    }
  }
}

## 8. 汇总输出（仍是“每个 GO 一行/一堆结果”，不是合并通路） ----
if (length(all_scores) > 0) {
  score_mat <- do.call(cbind, lapply(all_scores, function(x) {
    x[colnames(expr)]
  }))
  colnames(score_mat) <- names(all_scores)
  rownames(score_mat) <- colnames(expr)
  fwrite(data.table(sample = rownames(score_mat), as.data.table(score_mat)),
         file.path(out_dir, "01_pathway_scores_each_GO.csv"))

  # 各通路分数热图（每列仍是单独 GO 分数，不是合并基因集）
  anno_df <- data.frame(row.names = rownames(score_mat))
  clin_hm <- clinical_data[sample_std %in% rownames(score_mat)]
  if (nrow(clin_hm) > 0) {
    if ("stage_simplified" %in% names(clin_hm)) {
      anno_df$Stage <- clin_hm$stage_simplified[match(rownames(anno_df), clin_hm$sample_std)]
    }
    er_col <- first_present(names(clin_hm), c("er_status_by_ihc", "ER.Status"))
    if (!is.na(er_col)) {
      anno_df$ER <- as.character(clin_hm[[er_col]][match(rownames(anno_df), clin_hm$sample_std)])
    }
  }
  ha <- NULL
  if (ncol(anno_df) > 0) ha <- ComplexHeatmap::HeatmapAnnotation(df = anno_df)
  # ★★★ 已修改开始：热图缩放，禁止 t(scale(score_mat)) ★★★
  # 旧代码（会报 scale 找不到 matrix 方法，请勿运行）：
  #   z_score <- t(scale(score_mat))
  sm <- as.matrix(score_mat)
  storage.mode(sm) <- "double"
  z_score <- matrix(
    NA_real_,
    nrow = ncol(sm),
    ncol = nrow(sm),
    dimnames = list(colnames(sm), rownames(sm))
  )
  for (j in seq_len(ncol(sm))) {
    x <- sm[, j]
    sdx <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(sdx) || sdx < 1e-12) sdx <- 1
    z_score[j, ] <- (x - mean(x, na.rm = TRUE)) / sdx
  }
  z_score[!is.finite(z_score)] <- 0
  # ★★★ 已修改结束 ★★★
  # ★★★ 已修改开始：热图必须用 ComplexHeatmap::，不要直接 Heatmap()/draw() ★★★
  # 旧代码（heatmaps 包会抢走 Heatmap，报“参数没有用”）：
  #   draw(Heatmap(z_score, name = "z-score", ...))
  tryCatch({
    pdf_hm <- file.path(out_dir, "01_pathway_score_heatmap.pdf")
    pdf(pdf_hm, width = 12, height = 6)
    ht <- ComplexHeatmap::Heatmap(
      z_score,
      name = "z-score",
      top_annotation = ha,
      show_column_names = FALSE,
      row_names_gp = grid::gpar(fontsize = 8),
      column_title = "Each GO pathway score (not pooled)"
    )
    ComplexHeatmap::draw(ht)
    invisible(dev.off())
    message("  已保存通路分数热图：", pdf_hm)
  }, error = function(e) {
    if (dev.cur() > 1) invisible(dev.off())
    message("热图绘制失败：", conditionMessage(e))
  })
  # ★★★ 已修改结束 ★★★
}

if (length(all_clin) > 0) {
  clin_all <- rbindlist(all_clin, fill = TRUE)
  fwrite(clin_all, file.path(out_dir, "02_clinical_association_each_GO.csv"))
}

if (length(all_surv) > 0) {
  surv_all <- rbindlist(all_surv, fill = TRUE)
  fwrite(surv_all, file.path(out_dir, "03_survival_cox_each_GO.csv"))

  os_hl <- surv_all[endpoint == "OS" & model == "High_vs_Low"]
  if (nrow(os_hl) > 0) {
    os_hl[, lab := paste(GO, GO_name)]
    os_hl[, lab := factor(lab, levels = rev(lab))]
    # ★★★ 已修改开始：整段森林图请整块粘贴，geom_errorbar 行末尾必须有 + ★★★
    pfor <- ggplot(os_hl, aes(x = HR, y = lab)) +
      geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
      geom_errorbar(aes(xmin = HR_low, xmax = HR_high), orientation = "y", width = 0.2) +
      geom_point(aes(color = pvalue < 0.05), size = 3) +
      scale_x_log10() +
      labs(title = "OS Cox: High vs Low (each GO separately)",
           x = "Hazard ratio", y = NULL, color = "p < 0.05") +
      theme_bw()
    ggsave(file.path(out_dir, "03_OS_forest_each_GO.pdf"), pfor, width = 10, height = 6)
    message("  已保存 OS 森林图：", file.path(out_dir, "03_OS_forest_each_GO.pdf"))
    # ★★★ 已修改结束 ★★★
  }
}

if (length(all_neg) > 0) {
  neg_all <- rbindlist(all_neg, fill = TRUE)
  fwrite(neg_all, file.path(out_dir, "04_negative_genes_each_GO.csv"))
  neg_count <- neg_all[, .(n_negative_genes = .N,
                           n_strict_r0.3 = sum(spearman_r <= strict_r_cutoff)),
                       by = .(GO, GO_name)]
  fwrite(neg_count, file.path(out_dir, "04_negative_genes_count_each_GO.csv"))
}

if (length(all_prot) > 0) {
  fwrite(rbindlist(all_prot, fill = TRUE),
         file.path(out_dir, "05_negative_proteins_each_GO.csv"))
}

fwrite(go_set_summary, file.path(out_dir, "00_GO_gene_sets.csv"))
message("\n分析完成。每个 GO 的独立结果在：", normalizePath(out_dir, mustWork = FALSE))
message("请重点查看 per_GO/ 下各通路文件夹，以及 04_negative_genes_each_GO.csv")
  invisible(TRUE)
}
