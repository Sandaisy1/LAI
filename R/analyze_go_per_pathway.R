#!/usr/bin/env Rscript
# 每个 GO 通路单独：1) 与临床病人相关；2) 寻找负相关基因。禁止合并通路。
# 在项目根目录运行: Rscript R/analyze_go_per_pathway.R

local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[[1]])
    project_root <- dirname(dirname(normalizePath(script_path)))
    setwd(project_root)
  }
})

source(file.path("R", "config.R"))

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

ensure_packages <- function() {
  cran <- c("data.table", "matrixStats", "survival")
  bioc <- c("org.Hs.eg.db", "GO.db", "AnnotationDbi")
  optional_bioc <- c("GSVA")

  if (!isTRUE(auto_install_packages)) {
    return(invisible(NULL))
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  for (pkg in cran) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  for (pkg in bioc) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
  }
  for (pkg in optional_bioc) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      try(BiocManager::install(pkg, ask = FALSE, update = FALSE), silent = TRUE)
    }
  }
}

find_data_file <- function(dir, prefix) {
  hits <- list.files(dir, full.names = TRUE)
  hits <- hits[grepl(prefix, basename(hits), fixed = TRUE)]
  hits <- hits[!grepl("\\.(md|R|r)$", hits)]
  if (length(hits) == 0) {
    return(NA_character_)
  }
  hits[[1]]
}

normalize_barcode <- function(x, nchar_keep = barcode_sample_chars) {
  x <- toupper(gsub("\\.", "-", as.character(x)))
  x <- sub("\\s+$", "", x)
  ifelse(nchar(x) >= nchar_keep, substr(x, 1, nchar_keep), x)
}

sample_type_code <- function(barcode15) {
  ifelse(nchar(barcode15) >= 15, substr(barcode15, 14, 15), NA_character_)
}

strip_ensembl_version <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

read_table_auto <- function(path) {
  log_msg("读取", path)
  header <- readLines(path, n = 1L, warn = FALSE)
  sep <- if (grepl("\t", header, fixed = TRUE)) "\t" else "auto"
  data.table::fread(path, sep = sep, header = TRUE, data.table = TRUE, check.names = FALSE)
}

pick_id_column <- function(dt, candidates) {
  nms <- names(dt)
  hit <- candidates[candidates %in% nms]
  if (length(hit)) {
    return(hit[[1]])
  }
  hit <- nms[tolower(nms) %in% tolower(candidates)]
  if (length(hit)) {
    return(hit[[1]])
  }
  names(dt)[[1]]
}

matrix_from_feature_table <- function(dt) {
  id_col <- names(dt)[[1]]
  ids <- strip_ensembl_version(dt[[id_col]])
  mat <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(mat) <- "numeric"
  colnames(mat) <- normalize_barcode(colnames(mat), barcode_sample_chars)
  rownames(mat) <- make.unique(ids)
  if (anyDuplicated(ids)) {
    log_msg("基因 ID 有重复，按去版本号后的 ID 取均值")
    split_idx <- split(seq_len(nrow(mat)), ids)
    mat <- do.call(rbind, lapply(split_idx, function(i) {
      colMeans(mat[i, , drop = FALSE], na.rm = TRUE)
    }))
    rownames(mat) <- names(split_idx)
  }
  mat
}

filter_primary_tumors <- function(mat) {
  codes <- sample_type_code(colnames(mat))
  keep <- codes %in% keep_sample_types
  if (!any(keep)) {
    log_msg("未检测到样本类型代码，保留全部样本")
    return(mat)
  }
  mat <- mat[, keep, drop = FALSE]
  patients <- normalize_barcode(colnames(mat), barcode_patient_chars)
  dup <- duplicated(patients)
  if (any(dup)) {
    log_msg("同一患者多份原发瘤，各保留第一份")
    mat <- mat[, !dup, drop = FALSE]
  }
  mat
}

load_gtf_map <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(NULL)
  }
  gtf <- read_table_auto(path)
  id_col <- pick_id_column(gtf, c("gene_id", "Geneid", "ensembl_id", "Ensembl_ID"))
  name_col <- pick_id_column(gtf, c("gene_name", "GeneSymbol", "symbol", "gene_symbol"))
  type_col <- names(gtf)[tolower(names(gtf)) %in% c("gene_type", "gene_biotype", "biotype")]
  map <- data.frame(
    ensembl = strip_ensembl_version(gtf[[id_col]]),
    symbol = as.character(gtf[[name_col]]),
    gene_type = if (length(type_col)) as.character(gtf[[type_col[[1]]]]) else NA_character_,
    stringsAsFactors = FALSE
  )
  map <- map[!duplicated(map$ensembl), ]
  rownames(map) <- map$ensembl
  map
}

go_term_name <- function(go_id) {
  if (!is.null(go_names) && go_id %in% names(go_names)) {
    return(unname(go_names[[go_id]]))
  }
  if (requireNamespace("GO.db", quietly = TRUE)) {
    nm <- tryCatch(AnnotationDbi::Term(GO.db::GOTERM[[go_id]]), error = function(e) NA_character_)
    if (!is.na(nm)) return(nm)
  }
  go_id
}

genes_for_one_go <- function(go_id) {
  # 只取这一个 GO，绝不与其他通路合并
  raw <- tryCatch(
    AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = go_id,
      keytype = go_keytype,
      columns = c("ENSEMBL", "SYMBOL", "ENTREZID")
    ),
    error = function(e) {
      log_msg("GO 映射失败", go_id, conditionMessage(e))
      NULL
    }
  )
  if (is.null(raw) || nrow(raw) == 0) {
    return(data.frame(
      go_id = character(), ensembl = character(), symbol = character(),
      entrez = character(), stringsAsFactors = FALSE
    ))
  }
  data.frame(
    go_id = go_id,
    ensembl = strip_ensembl_version(raw$ENSEMBL),
    symbol = as.character(raw$SYMBOL),
    entrez = as.character(raw$ENTREZID),
    stringsAsFactors = FALSE
  )
}

match_geneset_to_expr <- function(gs, expr) {
  expr_id <- rownames(expr)
  by_ensembl <- which(gs$ensembl %in% expr_id)
  by_symbol <- which(gs$symbol %in% expr_id)
  keys <- unique(c(gs$ensembl[by_ensembl], gs$symbol[by_symbol]))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  intersect(keys, expr_id)
}

mean_z_score <- function(expr, genes) {
  sub <- expr[genes, , drop = FALSE]
  z <- matrixStats::rowZscores(sub, na.rm = TRUE)
  z[!is.finite(z)] <- NA_real_
  colMeans(z, na.rm = TRUE)
}

ssgsea_scores <- function(expr, geneset_list) {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    return(NULL)
  }
  geneset_list <- lapply(geneset_list, intersect, rownames(expr))
  geneset_list <- geneset_list[lengths(geneset_list) >= min_genes_per_go]
  if (!length(geneset_list)) {
    return(NULL)
  }
  ns <- asNamespace("GSVA")
  scores <- tryCatch({
    if (exists("ssgseaParam", envir = ns, inherits = FALSE)) {
      param <- GSVA::ssgseaParam(expr, geneset_list)
      as.matrix(GSVA::gsva(param, verbose = FALSE))
    } else {
      GSVA::gsva(expr, geneset_list, method = "ssgsea", kcdf = "Gaussian", verbose = FALSE)
    }
  }, error = function(e) {
    log_msg("ssGSEA 失败，改用 mean_z:", conditionMessage(e))
    NULL
  })
  scores
}

compute_pathway_scores <- function(expr, geneset_list) {
  # geneset_list 必须是 named list，每个名字一个 GO ID
  if (identical(score_method, "ssgsea")) {
    ssg <- ssgsea_scores(expr, geneset_list)
    if (!is.null(ssg)) {
      return(ssg)
    }
    log_msg("GSVA 不可用或 ssGSEA 失败，每个通路改用 mean_z")
  }
  out <- lapply(names(geneset_list), function(go_id) {
    genes <- intersect(geneset_list[[go_id]], rownames(expr))
    if (length(genes) < min_genes_per_go) {
      return(NULL)
    }
    mean_z_score(expr, genes)
  })
  names(out) <- names(geneset_list)
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) {
    stop("没有任何 GO 通路达到 min_genes_per_go")
  }
  do.call(rbind, out)
}

is_binary <- function(x) {
  ux <- unique(x[!is.na(x)])
  length(ux) == 2L
}

is_numeric_like <- function(x) {
  if (is.numeric(x)) return(TRUE)
  y <- suppressWarnings(as.numeric(as.character(x)))
  mean(!is.na(y) | is.na(x)) > 0.8 && length(unique(y[!is.na(y)])) > 5
}

guess_clinical_columns <- function(clin) {
  if (!is.null(clinical_columns)) {
    return(intersect(clinical_columns, names(clin)))
  }
  nms <- names(clin)
  keep <- vapply(nms, function(nm) {
    any(grepl(paste(clinical_keywords, collapse = "|"), nm, ignore.case = TRUE))
  }, logical(1))
  skip <- grepl("id$|barcode|uuid|submitter|patient|sample$|^_", nms, ignore.case = TRUE)
  nms[keep & !skip]
}

align_vector_to_samples <- function(ids, values, sample_ids) {
  names(values) <- ids
  patient <- normalize_barcode(sample_ids, barcode_patient_chars)
  out <- values[sample_ids]
  miss <- is.na(out)
  if (any(miss)) {
    out[miss] <- values[patient[miss]]
  }
  out
}

correlate_one_go_clinical <- function(score, clin) {
  sample_ids <- names(score)
  clin_id <- pick_id_column(clin, c("sample", "Sample", "submitter_id", "_PATIENT", "bcr_patient_barcode"))
  clin_ids15 <- normalize_barcode(clin[[clin_id]], barcode_sample_chars)
  clin_ids12 <- normalize_barcode(clin[[clin_id]], barcode_patient_chars)
  cols <- guess_clinical_columns(clin)
  if (!length(cols)) {
    log_msg("未自动识别到临床列，请在 config.R 设置 clinical_columns")
    return(data.frame())
  }
  rows <- lapply(cols, function(col) {
    raw <- clin[[col]]
    vals <- align_vector_to_samples(
      ifelse(nchar(clin_ids15) >= barcode_sample_chars, clin_ids15, clin_ids12),
      raw,
      sample_ids
    )
    ok <- !is.na(score) & !is.na(vals) & vals != "" & vals != "not reported" &
      !grepl("^\\[|unknown|not evaluated", vals, ignore.case = TRUE)
    n <- sum(ok)
    if (n < min_paired_samples) {
      return(NULL)
    }
    x <- score[ok]
    y <- vals[ok]
    if (is_numeric_like(y)) {
      ynum <- as.numeric(as.character(y))
      ct <- suppressWarnings(cor.test(x, ynum, method = "spearman", exact = FALSE))
      return(data.frame(
        variable = col, type = "continuous", n = n,
        method = "spearman", estimate = unname(ct$estimate),
        p_value = ct$p.value, extra = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
    y <- droplevels(factor(as.character(y)))
    k <- nlevels(y)
    if (k < 2) return(NULL)
    if (k == 2) {
      wt <- wilcox.test(x ~ y)
      return(data.frame(
        variable = col, type = "binary", n = n,
        method = "wilcoxon", estimate = unname(diff(tapply(x, y, median))),
        p_value = wt$p.value, extra = paste(levels(y), collapse = " vs "),
        stringsAsFactors = FALSE
      ))
    }
    kt <- kruskal.test(x ~ y)
    data.frame(
      variable = col, type = "categorical", n = n,
      method = "kruskal-wallis", estimate = unname(kt$statistic),
      p_value = kt$p.value, extra = paste(k, "levels"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0) {
    return(data.frame())
  }
  out$fdr <- p.adjust(out$p_value, method = "BH")
  out[order(out$p_value), ]
}

correlate_one_go_survival <- function(score, surv) {
  sid <- pick_id_column(surv, c("sample", "Sample", "submitter_id", "_PATIENT"))
  nms <- names(surv)
  os_col <- nms[grepl("^os$", nms, ignore.case = TRUE) |
                grepl("^overall_survival$", nms, ignore.case = TRUE)]
  time_col <- nms[grepl("^os[._]?time$", nms, ignore.case = TRUE) |
                  grepl("^overall_survival_time$", nms, ignore.case = TRUE)]
  if (!length(os_col) || !length(time_col)) {
    log_msg("survival 表缺少 OS / OS.time")
    return(NULL)
  }
  ids <- normalize_barcode(surv[[sid]], barcode_sample_chars)
  ids12 <- normalize_barcode(surv[[sid]], barcode_patient_chars)
  os <- align_vector_to_samples(ifelse(nchar(ids) >= barcode_sample_chars, ids, ids12), surv[[os_col[[1]]]], names(score))
  ost <- align_vector_to_samples(ifelse(nchar(ids) >= barcode_sample_chars, ids, ids12), surv[[time_col[[1]]]], names(score))
  os <- as.numeric(as.character(os))
  ost <- as.numeric(as.character(ost))
  ok <- !is.na(score) & !is.na(os) & !is.na(ost) & ost > 0
  if (sum(ok) < min_paired_samples) {
    return(NULL)
  }
  df <- data.frame(score = as.numeric(score[ok]), OS = os[ok], OS.time = ost[ok])
  fit <- survival::coxph(survival::Surv(OS.time, OS) ~ score, data = df)
  s <- summary(fit)
  list(
    table = data.frame(
      n = nrow(df),
      events = sum(df$OS == 1, na.rm = TRUE),
      hr = unname(s$coefficients[1, "exp(coef)"]),
      hr_l95 = s$conf.int[1, "lower .95"],
      hr_u95 = s$conf.int[1, "upper .95"],
      p_value = s$coefficients[1, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    ),
    df = df
  )
}

plot_km <- function(surv_df, go_id, outfile) {
  grp <- ifelse(surv_df$score >= median(surv_df$score, na.rm = TRUE), "High", "Low")
  surv_df$grp <- factor(grp, levels = c("Low", "High"))
  fit <- survival::survfit(survival::Surv(OS.time, OS) ~ grp, data = surv_df)
  pdf(outfile, width = 6, height = 5)
  plot(
    fit, col = c("#3C5488", "#E64B35"), lwd = 2, xlab = "Time", ylab = "Overall survival",
    main = paste0(go_id, " (", go_term_name(go_id), ")")
  )
  legend("bottomleft", legend = levels(surv_df$grp), col = c("#3C5488", "#E64B35"), lwd = 2, bty = "n")
  sd <- survival::survdiff(survival::Surv(OS.time, OS) ~ grp, data = surv_df)
  p <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  legend("topright", legend = sprintf("log-rank p = %.3g", p), bty = "n")
  dev.off()
}

spearman_with_p <- function(mat, score) {
  rho <- as.numeric(cor(t(mat), score, method = "spearman", use = "pairwise.complete.obs"))
  finite_score <- is.finite(score)
  n <- as.numeric(rowSums(is.finite(mat) & matrix(finite_score, nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE)))
  tstat <- rho * sqrt((n - 2) / pmax(1e-16, 1 - rho^2))
  pval <- 2 * pt(-abs(tstat), df = pmax(n - 2, 1))
  pval[n < min_paired_samples] <- NA_real_
  rho[n < min_paired_samples] <- NA_real_
  list(rho = rho, p_value = pval, n = n)
}

negative_gene_correlation_one_go <- function(score, expr, geneset_ids, gene_map) {
  common <- intersect(names(score), colnames(expr))
  score <- score[common]
  expr <- expr[, common, drop = FALSE]
  if (length(common) < min_paired_samples) {
    return(data.frame())
  }
  sp <- spearman_with_p(expr, score)
  tab <- data.frame(
    gene_id = rownames(expr),
    rho = sp$rho,
    p_value = sp$p_value,
    n = sp$n,
    stringsAsFactors = FALSE
  )
  tab <- tab[is.finite(tab$rho) & is.finite(tab$p_value), ]
  tab$fdr <- p.adjust(tab$p_value, method = "BH")
  tab$in_geneset <- tab$gene_id %in% geneset_ids
  if (!is.null(gene_map)) {
    tab$symbol <- gene_map[tab$gene_id, "symbol"]
    tab$gene_type <- gene_map[tab$gene_id, "gene_type"]
    miss <- is.na(tab$symbol) & tab$gene_id %in% gene_map$symbol
    if (any(miss)) {
      tab$symbol[miss] <- tab$gene_id[miss]
    }
  } else {
    tab$symbol <- tab$gene_id
    tab$gene_type <- NA_character_
  }
  tab <- tab[tab$rho < rho_cutoff & tab$fdr < fdr_cutoff, ]
  tab[order(tab$rho), ]
}

protein_correlation_one_go <- function(score, protein) {
  common <- intersect(names(score), colnames(protein))
  if (length(common) < min_paired_samples) {
    return(data.frame())
  }
  score <- score[common]
  protein <- protein[, common, drop = FALSE]
  sp <- spearman_with_p(protein, score)
  tab <- data.frame(
    protein = rownames(protein),
    rho = sp$rho,
    p_value = sp$p_value,
    n = sp$n,
    stringsAsFactors = FALSE
  )
  tab <- tab[is.finite(tab$rho) & is.finite(tab$p_value), ]
  if (nrow(tab) == 0) return(data.frame())
  tab$fdr <- p.adjust(tab$p_value, method = "BH")
  tab[order(tab$p_value), ]
}

write_tsv <- function(x, path) {
  if (is.null(x) || (is.data.frame(x) && nrow(x) == 0)) {
    x <- data.frame(note = "no_results")
  }
  data.table::fwrite(x, path, sep = "\t")
}

safe_go_dirname <- function(go_id) {
  gsub(":", "_", go_id, fixed = TRUE)
}

run_pipeline <- function() {
  ensure_packages()
  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(AnnotationDbi)
  })
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  expr_path <- find_data_file(data_dir, file_prefixes$expression)
  clin_path <- find_data_file(data_dir, file_prefixes$clinical)
  surv_path <- find_data_file(data_dir, file_prefixes$survival)
  prot_path <- find_data_file(data_dir, file_prefixes$protein)
  gtf_path  <- find_data_file(data_dir, file_prefixes$gtf)

  if (is.na(expr_path)) stop("找不到表达矩阵，前缀: ", file_prefixes$expression)
  if (is.na(clin_path)) stop("找不到临床表，前缀: ", file_prefixes$clinical)

  expr <- filter_primary_tumors(matrix_from_feature_table(read_table_auto(expr_path)))
  if (isTRUE(log2_fpkm)) {
    expr <- log2(pmax(expr, 0) + 1)
  }
  log_msg("表达矩阵", nrow(expr), "基因 x", ncol(expr), "样本")

  clin <- read_table_auto(clin_path)
  surv <- if (!is.na(surv_path)) read_table_auto(surv_path) else NULL
  prot <- if (!is.na(prot_path)) filter_primary_tumors(matrix_from_feature_table(read_table_auto(prot_path))) else NULL
  gene_map <- load_gtf_map(gtf_path)

  log_msg("逐个构建 GO 基因集（不合并）")
  genesets <- lapply(go_ids, genes_for_one_go)
  names(genesets) <- go_ids

  matched <- lapply(go_ids, function(id) match_geneset_to_expr(genesets[[id]], expr))
  names(matched) <- go_ids
  n_matched <- lengths(matched)
  log_msg("各通路在表达矩阵中的基因数:")
  print(n_matched)

  usable <- names(matched)[n_matched >= min_genes_per_go]
  skipped <- setdiff(go_ids, usable)
  if (length(skipped)) {
    log_msg("基因数不足，跳过:", paste(skipped, collapse = ", "))
  }
  if (!length(usable)) {
    stop("没有任何 GO 通路达到 min_genes_per_go")
  }

  log_msg("计算每个 GO 的通路评分")
  scores <- compute_pathway_scores(expr, matched[usable])

  summary_clin <- list()
  summary_neg <- list()

  scored_ids <- intersect(usable, rownames(scores))
  write_tsv(
    data.frame(
      go_id = go_ids,
      go_name = vapply(go_ids, go_term_name, character(1)),
      n_in_expression = as.integer(n_matched[go_ids]),
      used_in_analysis = go_ids %in% scored_ids,
      stringsAsFactors = FALSE
    ),
    file.path(output_dir, "go_geneset_sizes.tsv")
  )
  for (go_id in scored_ids) {
    # 每个通路独立：临床相关 + 负相关基因
    log_msg("====", go_id, go_term_name(go_id), "====")
    out_go <- file.path(output_dir, safe_go_dirname(go_id))
    dir.create(out_go, showWarnings = FALSE, recursive = TRUE)

    gs <- genesets[[go_id]]
    gs$in_expression <- gs$ensembl %in% rownames(expr) | gs$symbol %in% rownames(expr)
    write_tsv(gs, file.path(out_go, "geneset.tsv"))

    score <- scores[go_id, ]
    names(score) <- colnames(scores)
    write_tsv(
      data.frame(sample = names(score), pathway_score = as.numeric(score), stringsAsFactors = FALSE),
      file.path(out_go, "pathway_scores.tsv")
    )

    clin_tab <- correlate_one_go_clinical(score, clin)
    if (nrow(clin_tab)) {
      clin_tab$go_id <- go_id
      clin_tab$go_name <- go_term_name(go_id)
      clin_tab <- clin_tab[, c("go_id", "go_name", setdiff(names(clin_tab), c("go_id", "go_name")))]
    }
    write_tsv(clin_tab, file.path(out_go, "clinical_correlation.tsv"))
    summary_clin[[go_id]] <- clin_tab

    if (!is.null(surv)) {
      surv_res <- correlate_one_go_survival(score, surv)
      if (!is.null(surv_res)) {
        surv_tab <- surv_res$table
        surv_tab$go_id <- go_id
        surv_tab$go_name <- go_term_name(go_id)
        write_tsv(surv_tab, file.path(out_go, "survival_cox.tsv"))
        plot_km(surv_res$df, go_id, file.path(out_go, "survival_km.pdf"))
      }
    }

    neg <- negative_gene_correlation_one_go(score, expr, matched[[go_id]], gene_map)
    if (nrow(neg)) {
      neg$go_id <- go_id
      neg$go_name <- go_term_name(go_id)
      neg <- neg[, c("go_id", "go_name", setdiff(names(neg), c("go_id", "go_name")))]
    }
    write_tsv(neg, file.path(out_go, "negative_gene_correlation.tsv"))
    summary_neg[[go_id]] <- data.frame(
      go_id = go_id,
      go_name = go_term_name(go_id),
      n_geneset = length(matched[[go_id]]),
      n_negative_genes = nrow(neg),
      min_rho = if (nrow(neg)) min(neg$rho) else NA_real_,
      stringsAsFactors = FALSE
    )

    if (!is.null(prot)) {
      prot_tab <- protein_correlation_one_go(score, prot)
      if (nrow(prot_tab)) {
        prot_tab$go_id <- go_id
        prot_neg <- prot_tab[prot_tab$rho < rho_cutoff & prot_tab$fdr < fdr_cutoff, ]
        write_tsv(prot_tab, file.path(out_go, "protein_correlation.tsv"))
        write_tsv(prot_neg, file.path(out_go, "negative_protein_correlation.tsv"))
      }
    }
  }

  summary_clin <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, summary_clin)
  write_tsv(if (length(summary_clin)) do.call(rbind, summary_clin) else NULL,
            file.path(output_dir, "summary_clinical_correlation.tsv"))
  write_tsv(do.call(rbind, summary_neg), file.path(output_dir, "summary_negative_gene_counts.tsv"))

  skipped_df <- data.frame(
    go_id = skipped,
    go_name = vapply(skipped, go_term_name, character(1)),
    n_matched_genes = unname(n_matched[skipped]),
    stringsAsFactors = FALSE
  )
  write_tsv(skipped_df, file.path(output_dir, "skipped_go_terms.tsv"))

  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  log_msg("完成。每个 GO 的结果在", output_dir, "/<GO_ID>/")
}

if (!interactive() || identical(Sys.getenv("RUN_GO_PIPELINE"), "1")) {
  run_pipeline()
}
