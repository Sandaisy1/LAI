#!/usr/bin/env Rscript
# =============================================================================
# 病人转移 / 预后 × 神经相关信号通路
# 独立脚本：不修改 TG_RNAseq_pipeline.R 的六组细胞系比较。
#
# 比较：有转移 vs 无转移（能估计 p，先 p < 0.01，再 FC / topN，只分析上调）
# 预后：生存曲线 + Cox；若有 good/poor 列则再比较通路分数
# 神经通路：专项分数、箱线、热图、Focused_neural GSEA/ORA
#           不改全库 GO/KEGG 的 p 值
#
# 用法：
#   setwd("E:/R/TG_BRCA/TG")
#   source("TG_RNAseq_patient_metastasis_neural.R")
#
# 输入（CSV/TSV）：
#   patient_expression.csv   第一列基因符号，其余列样品
#   patient_clinical.csv     样品、转移、生存/预后
# 也可用环境变量 TG_PATIENT_EXPR / TG_PATIENT_CLINICAL
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)

if (!exists("log_msg", mode = "function")) {
  log_msg <- function(...) {
    cat(paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = "")), "\n")
  }
}
if (!exists("has_pkg", mode = "function")) {
  has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
}

p_cutoff_patient <- 0.01
fc_cutoffs_patient <- c("FC_1" = 1, "FC_1.25" = 1.25, "FC_1.5" = 1.5, "FC_2" = 2)
top_ns_patient <- c(50, 75, 100, 150, 200, 250, 300)

# -----------------------------------------------------------------------------
# 神经相关 GO / KEGG（官方现行条目）+ 符号回退
# -----------------------------------------------------------------------------
neural_go_sets <- function() {
  c(
    NEURAL_GO_synaptic_signaling = "GO:0099536",
    NEURAL_GO_chemical_synaptic_transmission = "GO:0007268",
    NEURAL_GO_neurotransmitter_secretion = "GO:0007269",
    NEURAL_GO_neurotransmitter_transport = "GO:0006836",
    NEURAL_GO_neuropeptide_signaling = "GO:0007218",
    NEURAL_GO_neurotrophin_signaling = "GO:0038179",
    NEURAL_GO_neurotrophin_TRK_signaling = "GO:0048011",
    NEURAL_GO_glutamate_receptor_signaling = "GO:0007215",
    NEURAL_GO_GABA_signaling = "GO:0007214",
    NEURAL_GO_dopamine_receptor_signaling = "GO:0007212",
    NEURAL_GO_axon_guidance = "GO:0007411",
    NEURAL_GO_semaphorin_plexin_signaling = "GO:0071526",
    NEURAL_GO_ephrin_receptor_signaling = "GO:0048013",
    NEURAL_GO_neuron_projection_development = "GO:0031175",
    NEURAL_GO_neuron_development = "GO:0048666",
    NEURAL_GO_neuron_differentiation = "GO:0030182",
    NEURAL_GO_neuron_migration = "GO:0001764",
    NEURAL_GO_synapse_organization = "GO:0050808",
    NEURAL_GO_synaptic_plasticity = "GO:0048167",
    NEURAL_GO_nervous_system_development = "GO:0007399",
    NEURAL_GO_neurogenesis = "GO:0022008"
  )
}

neural_kegg_sets <- function() {
  c(
    NEURAL_KEGG_neurotrophin_signaling = "04722",
    NEURAL_KEGG_axon_guidance = "04360",
    NEURAL_KEGG_synaptic_vesicle_cycle = "04721",
    NEURAL_KEGG_glutamatergic_synapse = "04724",
    NEURAL_KEGG_cholinergic_synapse = "04725",
    NEURAL_KEGG_serotonergic_synapse = "04726",
    NEURAL_KEGG_GABAergic_synapse = "04727",
    NEURAL_KEGG_dopaminergic_synapse = "04728",
    NEURAL_KEGG_long_term_potentiation = "04720",
    NEURAL_KEGG_neuroactive_ligand_receptor = "04080"
  )
}

neural_symbol_fallback <- function() {
  list(
    NEURAL_GO_neurotrophin_TRK_signaling = c(
      "NTRK1", "NTRK2", "NTRK3", "NGF", "BDNF", "NTF3", "NTF4", "NGFR",
      "NTRK2", "SHC1", "GRB2", "SOS1", "MAPK1", "MAPK3", "AKT1", "PLCG1"
    ),
    NEURAL_GO_axon_guidance = c(
      "SEMA3A", "SEMA3F", "NRP1", "NRP2", "PLXNA1", "ROBO1", "SLIT2",
      "DCC", "NTN1", "EPHA4", "EFNB2", "UNC5B", "CXCR4", "MET"
    ),
    NEURAL_GO_synaptic_signaling = c(
      "SNAP25", "SYP", "DLG4", "SYN1", "STX1A", "VAMP2", "GRIN1", "GRIA1",
      "GABRA1", "SLC17A7", "SLC32A1"
    ),
    NEURAL_GO_glutamate_receptor_signaling = c(
      "GRIN1", "GRIN2A", "GRIN2B", "GRIA1", "GRIA2", "GRM1", "GRM5", "SLC1A2"
    ),
    NEURAL_GO_GABA_signaling = c(
      "GABRA1", "GABRB2", "GABRG2", "GABBR1", "GAD1", "GAD2", "SLC6A1"
    ),
    NEURAL_GO_dopamine_receptor_signaling = c(
      "DRD1", "DRD2", "TH", "SLC6A3", "DDC", "SLC18A2"
    ),
    NEURAL_KEGG_cholinergic_synapse = c(
      "CHAT", "ACHE", "SLC18A3", "CHRNA4", "CHRNB2", "CHRM1"
    ),
    NEURAL_ALL_core = c(
      "NTRK1", "NTRK2", "NTRK3", "NGF", "BDNF", "NTF3", "NGFR",
      "SNAP25", "SYP", "DLG4", "GRIN1", "GRIA1", "GABRA1", "TH", "SLC6A3",
      "CHAT", "SLC18A3", "SEMA3A", "ROBO1", "DCC", "EPHA4", "EFNB2", "NTN1"
    )
  )
}

neural_keyword_pattern <- function() {
  paste(
    "neuro", "axon", "synap", "neurotrophin", "neuropeptide",
    "neurotransmitter", "glutamat", "gaba", "dopamine", "serotonin",
    "cholinergic", "adrenergic", "semaphorin", "ephrin", "netrin",
    "neuron", "glia", "myelin", "nerve growth",
    sep = "|"
  )
}

classify_neural_term <- function(id, desc) {
  id <- as.character(id)[1]
  desc <- as.character(desc)[1]
  if (is.na(id)) id <- ""
  if (is.na(desc)) desc <- ""
  if (grepl("^NEURAL_", id)) return("neural")
  kegg <- unname(neural_kegg_sets())
  go <- unname(neural_go_sets())
  id_short <- sub("^hsa", "", id)
  if (id %in% go || id %in% kegg || id_short %in% kegg ||
      paste0("hsa", id_short) %in% paste0("hsa", kegg)) {
    return("neural")
  }
  txt <- tolower(paste(id, desc))
  if (grepl(neural_keyword_pattern(), txt, perl = TRUE, ignore.case = TRUE)) {
    return("neural")
  }
  NA_character_
}

# -----------------------------------------------------------------------------
# 临床表解析
# -----------------------------------------------------------------------------
normalize_key <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", as.character(x)))
}

find_column <- function(df, patterns) {
  nms <- names(df)
  key <- normalize_key(nms)
  pats <- normalize_key(patterns)
  for (p in pats) {
    hit <- which(key == p)
    if (length(hit) > 0) return(nms[hit[1]])
  }
  for (p in pats) {
    hit <- which(grepl(p, key, fixed = TRUE))
    if (length(hit) > 0) return(nms[hit[1]])
  }
  NA_character_
}

read_table_flex <- function(path) {
  if (!file.exists(path)) stop("找不到文件: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls") && has_pkg("readxl")) {
    return(as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE))
  }
  sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
  first <- readLines(path, n = 1, warn = FALSE)
  if (grepl("\t", first) && !grepl(",", first)) sep <- "\t"
  utils::read.table(
    path, header = TRUE, sep = sep, check.names = FALSE,
    stringsAsFactors = FALSE, comment.char = "", quote = "\""
  )
}

label_metastasis <- function(clinical) {
  out <- clinical
  met_col <- find_column(
    clinical,
    c("metastasis", "metastatic", "distant_metastasis", "distantmetastasis",
      "pathologic_M", "pathologicm", "ajcc_pathologic_m", "ajccpathologicm",
      "m_stage", "mstage", "path_m", "ajcc_m")
  )
  raw <- if (!is.na(met_col)) as.character(clinical[[met_col]]) else rep(NA_character_, nrow(clinical))
  raw <- trimws(raw)
  lab <- rep(NA_character_, length(raw))
  up <- toupper(gsub("[^A-Za-z0-9]", "", raw))
  lab[grepl("^(M1|YES|Y|TRUE|1|MET|METASTASIS|METASTATIC|DISTANT)$", up)] <- "Metastasis"
  lab[grepl("^(M0|NO|N|FALSE|0|NONMET|NOMETASTASIS|ABSENT|NONE)$", up)] <- "NoMetastasis"
  lab[grepl("MX|UNKNOWN|NA|NOTASSESSED", up)] <- NA_character_
  # 中文
  lab[grepl("无转移|未转移|没有转移", raw)] <- "NoMetastasis"
  lab[grepl("转移", raw) & is.na(lab)] <- "Metastasis"
  out$metastasis_group <- lab
  out$metastasis_source <- if (is.na(met_col)) NA_character_ else met_col
  out
}

label_survival <- function(clinical) {
  out <- clinical
  time_col <- find_column(
    clinical,
    c("os_time", "ostime", "OS.time", "days_to_death", "daystodeath",
      "days_to_last_follow_up", "daystolastfollowup", "survival_days",
      "survival_time", "futime", "osdays")
  )
  death_col <- find_column(clinical, c("days_to_death", "daystodeath"))
  follow_col <- find_column(clinical, c("days_to_last_follow_up", "daystolastfollowup"))
  event_col <- find_column(
    clinical,
    c("os_event", "osevent", "OS", "OS.event", "vital_status", "vitalstatus",
      "fustat", "event", "status", "os_status")
  )
  time <- if (!is.na(time_col)) suppressWarnings(as.numeric(clinical[[time_col]])) else rep(NA_real_, nrow(clinical))
  if ((all(!is.finite(time)) || is.na(time_col)) && !is.na(death_col) && !is.na(follow_col)) {
    dth <- suppressWarnings(as.numeric(clinical[[death_col]]))
    fol <- suppressWarnings(as.numeric(clinical[[follow_col]]))
    time <- ifelse(is.finite(dth), dth, fol)
  }
  event <- rep(NA_real_, nrow(clinical))
  if (!is.na(event_col)) {
    ev <- toupper(trimws(as.character(clinical[[event_col]])))
    event[grepl("^(1|DEAD|DECEASED|DEATH|TRUE|YES|EVENT)$", ev)] <- 1
    event[grepl("^(0|ALIVE|LIVING|FALSE|NO|CENSOR|CENSORED)$", ev)] <- 0
  }
  out$os_time <- time
  out$os_event <- event
  out
}

label_prognosis <- function(clinical) {
  out <- clinical
  col <- find_column(clinical, c("prognosis", "risk", "risk_group", "outcome_group"))
  lab <- rep(NA_character_, nrow(clinical))
  if (!is.na(col)) {
    raw <- tolower(trimws(as.character(clinical[[col]])))
    lab[grepl("poor|high|bad|unfav", raw)] <- "PoorPrognosis"
    lab[grepl("good|low|fav", raw)] <- "GoodPrognosis"
  }
  out$prognosis_group <- lab
  out
}

prepare_clinical <- function(clinical) {
  id_col <- find_column(
    clinical,
    c("sample", "sample_id", "sampleid", "barcode", "patient", "patient_id",
      "submitter_id", "bcr_patient_barcode")
  )
  if (is.na(id_col)) {
    stop("临床表找不到样品列（sample / sample_id / barcode / patient）")
  }
  out <- clinical
  out$sample <- as.character(out[[id_col]])
  out <- label_metastasis(out)
  out <- label_survival(out)
  out <- label_prognosis(out)
  out
}

read_expression_matrix <- function(path) {
  df <- read_table_flex(path)
  if (ncol(df) < 2) stop("表达矩阵至少需要一列基因和一列样品: ", path)
  genes <- as.character(df[[1]])
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  genes <- ifelse(is.na(genes) | !nzchar(genes), paste0("row", seq_along(genes)), genes)
  if (exists("pick_official_symbol", mode = "function")) {
    genes <- vapply(genes, pick_official_symbol, character(1), USE.NAMES = FALSE)
    genes[is.na(genes)] <- as.character(df[[1]])[is.na(genes)]
  }
  keep <- !is.na(genes) & nzchar(genes)
  mat <- mat[keep, , drop = FALSE]
  genes <- genes[keep]
  # 重复符号取均值
  if (anyDuplicated(genes)) {
    log_msg("Duplicate gene symbols: averaging ", sum(duplicated(genes)), " extras")
    split_idx <- split(seq_along(genes), genes)
    um <- vapply(split_idx, function(ii) {
      if (length(ii) == 1) return(ii[1])
      0
    }, numeric(1))
    # rebuild
    ug <- names(split_idx)
    newm <- matrix(NA_real_, nrow = length(ug), ncol = ncol(mat), dimnames = list(ug, colnames(mat)))
    for (g in ug) {
      ii <- split_idx[[g]]
      newm[g, ] <- if (length(ii) == 1) mat[ii, ] else colMeans(mat[ii, , drop = FALSE], na.rm = TRUE)
    }
    mat <- newm
  } else {
    rownames(mat) <- genes
  }
  colnames(mat) <- as.character(colnames(mat))
  mat
}

align_patient_data <- function(mat, clinical) {
  clin <- prepare_clinical(clinical)
  expr_ids <- colnames(mat)
  clin_ids <- clin$sample
  exact <- intersect(expr_ids, clin_ids)
  if (length(exact) >= 4) {
    mat2 <- mat[, exact, drop = FALSE]
    clin2 <- clin[match(exact, clin$sample), , drop = FALSE]
  } else {
    short <- function(x) toupper(gsub("[^A-Z0-9]", "", substr(x, 1, 15)))
    em <- short(expr_ids)
    cm <- short(clin_ids)
    map <- match(em, cm)
    keep <- which(!is.na(map))
    if (length(keep) < 4) {
      stop("表达矩阵列名与临床样品对不上。表达列示例: ",
           paste(utils::head(expr_ids, 3), collapse = ", "),
           " ；临床示例: ", paste(utils::head(clin_ids, 3), collapse = ", "))
    }
    mat2 <- mat[, keep, drop = FALSE]
    clin2 <- clin[map[keep], , drop = FALSE]
    clin2$sample <- colnames(mat2)
  }
  clin2$group <- clin2$metastasis_group
  list(mat = mat2, clinical = clin2)
}

find_patient_files <- function(project_dir) {
  expr_env <- Sys.getenv("TG_PATIENT_EXPR", unset = "")
  clin_env <- Sys.getenv("TG_PATIENT_CLINICAL", unset = "")
  if (nzchar(expr_env) && nzchar(clin_env)) {
    return(list(expr = expr_env, clinical = clin_env))
  }
  roots <- unique(c(
    project_dir,
    file.path(project_dir, "patient"),
    file.path(project_dir, "patient_metastasis_neural")
  ))
  if (identical(Sys.getenv("TG_PATIENT_USE_EXAMPLE"), "1")) {
    roots <- c(roots, file.path(getwd(), "examples", "patient_metastasis_neural"))
  }
  expr_names <- c(
    "patient_expression.csv", "patient_expression.tsv", "patient_expression.txt",
    "patient_expr.csv", "expression.csv"
  )
  clin_names <- c(
    "patient_clinical.csv", "patient_clinical.tsv", "patient_clinical.txt",
    "clinical.csv", "patient_info.csv"
  )
  for (root in roots) {
    ehit <- expr_names[file.exists(file.path(root, expr_names))]
    chit <- clin_names[file.exists(file.path(root, clin_names))]
    if (length(ehit) > 0 && length(chit) > 0) {
      return(list(expr = file.path(root, ehit[1]), clinical = file.path(root, chit[1])))
    }
  }
  NULL
}

# -----------------------------------------------------------------------------
# 过滤 / 标准化 / 差异（病人两组，可估计 p）
# -----------------------------------------------------------------------------
patient_filter_low <- function(mat, sample_info, value_type) {
  n_group <- table(sample_info$group, useNA = "no")
  min_n <- max(2, min(n_group[n_group > 0], na.rm = TRUE))
  if (value_type == "counts") {
    keep <- rowSums(mat >= 10, na.rm = TRUE) >= min_n
  } else {
    keep <- rowSums(mat > 1, na.rm = TRUE) >= min_n
  }
  if (sum(keep) < 50) {
    keep <- rowSums(is.finite(mat) & mat > 0, na.rm = TRUE) >= min_n
    log_msg("Patient filter fallback: keep expressed-in-", min_n, "-samples")
  }
  log_msg("Patient low-expression filter: keep ", sum(keep), " / ", nrow(mat))
  mat[keep, , drop = FALSE]
}

patient_normalize <- function(mat, sample_info, value_type) {
  grp <- factor(sample_info$group)
  if (value_type == "counts" && has_pkg("DESeq2")) {
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(pmax(mat, 0)),
      colData = data.frame(row.names = colnames(mat), group = grp),
      design = ~ group
    )
    dds <- DESeq2::estimateSizeFactors(dds)
    log_msg("Patient normalization: DESeq2 size factor")
    log_mat <- log2(DESeq2::counts(dds, normalized = TRUE) + 1)
    heat <- tryCatch(
      SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE)),
      error = function(e) log_mat
    )
    return(list(log_mat = log_mat, heat_mat = heat, method = "DESeq2"))
  }
  log_msg("Patient normalization: log2 + quantile")
  log_mat <- log2(pmax(mat, 0) + 1)
  if (has_pkg("limma")) {
    log_mat <- limma::normalizeQuantiles(log_mat)
  } else {
    # 分位数标准化（无 limma）
    rnk <- apply(log_mat, 2, rank, ties.method = "average")
    target <- sort(rowMeans(apply(log_mat, 2, sort)))
    log_mat <- apply(rnk, 2, function(r) {
      o <- order(r)
      out <- numeric(length(r))
      out[o] <- stats::approx(seq_along(target), target, xout = r[o], rule = 2)$y
      out
    })
    rownames(log_mat) <- rownames(mat)
    colnames(log_mat) <- colnames(mat)
  }
  list(log_mat = log_mat, heat_mat = log_mat, method = "log2_quantile")
}

patient_de_metastasis <- function(log_mat, sample_info, raw_mat = NULL, value_type = "fpkm") {
  ok <- sample_info$group %in% c("Metastasis", "NoMetastasis")
  si <- sample_info[ok, , drop = FALSE]
  if (sum(si$group == "Metastasis") < 2 || sum(si$group == "NoMetastasis") < 2) {
    stop("转移 / 无转移两组都至少需要 2 个样品才能估计 p 值")
  }
  mat <- log_mat[, si$sample, drop = FALSE]
  log_msg(
    "DE Metastasis (n=", sum(si$group == "Metastasis"),
    ") vs NoMetastasis (n=", sum(si$group == "NoMetastasis"), ")"
  )
  if (value_type == "counts" && !is.null(raw_mat) && has_pkg("DESeq2")) {
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(pmax(raw_mat[, si$sample, drop = FALSE], 0)),
      colData = data.frame(
        row.names = si$sample,
        group = factor(si$group, levels = c("NoMetastasis", "Metastasis"))
      ),
      design = ~ group
    )
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    res <- DESeq2::results(dds, contrast = c("group", "Metastasis", "NoMetastasis"))
    return(data.frame(
      gene = rownames(res),
      log2FC = as.numeric(res$log2FoldChange),
      AveExpr = as.numeric(res$baseMean),
      pvalue = as.numeric(res$pvalue),
      padj = as.numeric(res$padj),
      stringsAsFactors = FALSE
    ))
  }
  if (has_pkg("limma")) {
    grp <- factor(si$group, levels = c("NoMetastasis", "Metastasis"))
    design <- stats::model.matrix(~ grp)
    fit <- limma::eBayes(limma::lmFit(mat, design))
    tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
    return(data.frame(
      gene = rownames(mat),
      log2FC = tt$logFC,
      AveExpr = tt$AveExpr,
      pvalue = tt$P.Value,
      padj = tt$adj.P.Val,
      stringsAsFactors = FALSE
    ))
  }
  log_msg("limma/DESeq2 不可用，改用每基因 Wilcoxon（仍报告真实 p，不伪造）")
  met <- si$sample[si$group == "Metastasis"]
  nom <- si$sample[si$group == "NoMetastasis"]
  res <- lapply(seq_len(nrow(mat)), function(i) {
    a <- as.numeric(mat[i, met])
    b <- as.numeric(mat[i, nom])
    pv <- tryCatch(stats::wilcox.test(a, b, exact = FALSE)$p.value, error = function(e) NA_real_)
    data.frame(
      gene = rownames(mat)[i],
      log2FC = mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE),
      AveExpr = mean(c(a, b), na.rm = TRUE),
      pvalue = pv,
      stringsAsFactors = FALSE
    )
  })
  de <- do.call(rbind, res)
  de$padj <- stats::p.adjust(de$pvalue, method = "BH")
  de
}

passes_p01 <- function(pvalue) {
  !is.na(pvalue) & pvalue < p_cutoff_patient
}

select_up_fc_p01 <- function(de, fc) {
  keep <- passes_p01(de$pvalue) & !is.na(de$log2FC) & (2^de$log2FC >= fc)
  de[keep, , drop = FALSE]
}

select_up_topn_p01 <- function(de, n) {
  x <- de[passes_p01(de$pvalue) & !is.na(de$log2FC) & de$log2FC > 0, , drop = FALSE]
  x <- x[order(x$log2FC, decreasing = TRUE), , drop = FALSE]
  utils::head(x, n)
}

# -----------------------------------------------------------------------------
# 通路分数
# -----------------------------------------------------------------------------
build_neural_gene_sets_symbols <- function(universe) {
  sets <- neural_symbol_fallback()
  if (exists("entrez_for_go", mode = "function") && has_pkg("org.Hs.eg.db")) {
    log_msg("Building neural gene sets from GO/KEGG via org.Hs.eg.db")
    go <- neural_go_sets()
    for (nm in names(go)) {
      ez <- tryCatch(entrez_for_go(go[[nm]]), error = function(e) character())
      if (length(ez) >= 8 && exists("map_to_entrez", mode = "function")) {
        # map entrez back through AnnotationDbi
        sym <- tryCatch(
          AnnotationDbi::mapIds(org.Hs.eg.db, keys = ez, column = "SYMBOL",
                                keytype = "ENTREZID", multiVals = "first"),
          error = function(e) NULL
        )
        if (!is.null(sym)) sets[[nm]] <- unique(c(sets[[nm]], unname(sym)))
      } else if (length(ez) >= 8) {
        sets[[nm]] <- unique(c(sets[[nm]], ez))
      }
    }
    kegg <- neural_kegg_sets()
    if (exists("entrez_for_kegg_path", mode = "function")) {
      for (nm in names(kegg)) {
        ez <- tryCatch(entrez_for_kegg_path(kegg[[nm]]), error = function(e) character())
        sym <- tryCatch(
          AnnotationDbi::mapIds(org.Hs.eg.db, keys = ez, column = "SYMBOL",
                                keytype = "ENTREZID", multiVals = "first"),
          error = function(e) NULL
        )
        if (!is.null(sym)) sets[[nm]] <- unique(c(sets[[nm]], unname(sym)))
      }
    }
  } else {
    log_msg("org.Hs.eg.db 不可用，使用内置神经通路符号回退集")
  }
  sets <- lapply(sets, function(g) intersect(unique(as.character(g)), universe))
  sets[vapply(sets, length, integer(1)) >= 3]
}

pathway_mean_z <- function(log_mat, gene_sets) {
  z <- t(scale(t(log_mat)))
  z[!is.finite(z)] <- 0
  sc <- sapply(gene_sets, function(genes) {
    genes <- intersect(genes, rownames(z))
    if (length(genes) == 0) return(rep(NA_real_, ncol(z)))
    colMeans(z[genes, , drop = FALSE], na.rm = TRUE)
  })
  if (is.null(dim(sc))) {
    sc <- matrix(sc, ncol = 1, dimnames = list(colnames(log_mat), names(gene_sets)))
  }
  rownames(sc) <- colnames(log_mat)
  as.data.frame(sc, check.names = FALSE)
}

compare_scores_by_group <- function(scores, group, group_a = "Metastasis", group_b = "NoMetastasis") {
  paths <- colnames(scores)
  rows <- lapply(paths, function(p) {
    a <- as.numeric(scores[group == group_a, p])
    b <- as.numeric(scores[group == group_b, p])
    a <- a[is.finite(a)]
    b <- b[is.finite(b)]
    if (length(a) < 2 || length(b) < 2) {
      return(data.frame(
        pathway = p, n_a = length(a), n_b = length(b),
        mean_a = mean(a), mean_b = mean(b), delta = mean(a) - mean(b),
        wilcox_p = NA_real_, stringsAsFactors = FALSE
      ))
    }
    pv <- tryCatch(stats::wilcox.test(a, b, exact = FALSE)$p.value, error = function(e) NA_real_)
    data.frame(
      pathway = p,
      n_a = length(a), n_b = length(b),
      mean_a = mean(a), mean_b = mean(b),
      delta = mean(a) - mean(b),
      wilcox_p = pv,
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  names(tab)[names(tab) == "mean_a"] <- paste0("mean_", group_a)
  names(tab)[names(tab) == "mean_b"] <- paste0("mean_", group_b)
  tab$wilcox_padj <- stats::p.adjust(tab$wilcox_p, method = "BH")
  tab <- tab[order(tab$wilcox_p, na.last = TRUE), ]
  tab
}

# -----------------------------------------------------------------------------
# 作图
# -----------------------------------------------------------------------------
patient_save_gg <- function(plot, stub, width = 8, height = 6) {
  dir.create(dirname(stub), recursive = TRUE, showWarnings = FALSE)
  if (exists("save_gg", mode = "function")) {
    save_gg(plot, stub, width = width, height = height)
    return(invisible(TRUE))
  }
  if (!has_pkg("ggplot2")) return(invisible(FALSE))
  tryCatch(ggplot2::ggsave(paste0(stub, ".pdf"), plot, width = width, height = height),
           error = function(e) log_msg("pdf failed: ", e$message))
  tryCatch(ggplot2::ggsave(paste0(stub, ".png"), plot, width = width, height = height, dpi = 300),
           error = function(e) log_msg("png failed: ", e$message))
  invisible(TRUE)
}

plot_score_boxplot <- function(scores, group, title, outfile, group_levels = NULL) {
  df <- data.frame(group = group, scores, check.names = FALSE, stringsAsFactors = FALSE)
  long <- tidyr_like_pivot(df, "group")
  if (!is.null(group_levels)) {
    long <- long[long$group %in% group_levels, , drop = FALSE]
    long$group <- factor(long$group, levels = group_levels)
  } else {
    long <- long[!is.na(long$group), , drop = FALSE]
  }
  if (nrow(long) == 0) return(invisible(NULL))
  if (has_pkg("ggplot2")) {
    pal <- c(Metastasis = "#D62828", NoMetastasis = "#4C78A8",
             PoorPrognosis = "#F58518", GoodPrognosis = "#54A24B")
    p <- ggplot2::ggplot(long, ggplot2::aes(x = group, y = score, fill = group)) +
      ggplot2::geom_boxplot(outlier.size = 0.8, width = 0.65) +
      ggplot2::geom_jitter(width = 0.12, alpha = 0.35, size = 0.7) +
      ggplot2::facet_wrap(~ pathway, scales = "free_y") +
      ggplot2::scale_fill_manual(values = pal, na.value = "grey70") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1), legend.position = "none") +
      ggplot2::labs(title = title, x = NULL, y = "mean z-score")
    patient_save_gg(p, outfile, width = 12, height = 9)
  } else {
    grDevices::pdf(paste0(outfile, ".pdf"), width = 12, height = 9)
    paths <- unique(long$pathway)
    nc <- min(3, length(paths))
    nr <- ceiling(length(paths) / nc)
    graphics::par(mfrow = c(nr, nc), mar = c(6, 4, 2, 1))
    for (pw in paths) {
      sub <- long[long$pathway == pw, ]
      graphics::boxplot(score ~ group, data = sub, main = pw, las = 2, ylab = "mean z-score")
    }
    grDevices::dev.off()
  }
}

tidyr_like_pivot <- function(df, id_col) {
  paths <- setdiff(names(df), id_col)
  do.call(rbind, lapply(paths, function(p) {
    data.frame(
      group = df[[id_col]],
      pathway = p,
      score = as.numeric(df[[p]]),
      stringsAsFactors = FALSE
    )
  }))
}

plot_patient_heatmap <- function(heat_mat, sample_info, genes, title, outfile, max_genes = 80) {
  genes <- intersect(genes, rownames(heat_mat))
  if (length(genes) < 2) {
    writeLines("fewer than 2 genes", paste0(outfile, "_EMPTY.txt"))
    return(invisible(NULL))
  }
  if (length(genes) > max_genes) genes <- genes[seq_len(max_genes)]
  sub <- heat_mat[genes, sample_info$sample, drop = FALSE]
  ann <- data.frame(
    Metastasis = sample_info$metastasis_group,
    row.names = sample_info$sample
  )
  if ("prognosis_group" %in% names(sample_info) && any(!is.na(sample_info$prognosis_group))) {
    ann$Prognosis <- sample_info$prognosis_group
  }
  pal <- list(
    Metastasis = c(Metastasis = "#D62828", NoMetastasis = "#4C78A8"),
    Prognosis = c(PoorPrognosis = "#F58518", GoodPrognosis = "#54A24B")
  )
  if (has_pkg("pheatmap")) {
    hm_col <- if (has_pkg("RColorBrewer")) {
      grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100)
    } else {
      grDevices::colorRampPalette(c("blue", "white", "red"))(100)
    }
    grDevices::pdf(paste0(outfile, ".pdf"), width = max(8, min(16, 6 + 0.04 * ncol(sub))),
                   height = max(6, min(16, 0.16 * nrow(sub) + 3)))
    pheatmap::pheatmap(
      sub, scale = "row", annotation_col = ann, annotation_colors = pal,
      show_rownames = nrow(sub) <= 80, fontsize_row = 6, main = title, color = hm_col
    )
    grDevices::dev.off()
    grDevices::png(paste0(outfile, ".png"), width = 2400,
                   height = max(1800, 40 * nrow(sub) + 400), res = 300)
    pheatmap::pheatmap(
      sub, scale = "row", annotation_col = ann, annotation_colors = pal,
      show_rownames = nrow(sub) <= 80, fontsize_row = 6, main = title, color = hm_col
    )
    grDevices::dev.off()
  } else {
    grDevices::pdf(paste0(outfile, ".pdf"), width = 8, height = 6)
    graphics::image(t(scale(t(sub))), main = title, xlab = "samples", ylab = "genes")
    grDevices::dev.off()
  }
}

plot_survival_km <- function(time, event, group, title, outfile) {
  keep <- is.finite(time) & time > 0 & !is.na(event) & !is.na(group)
  if (sum(keep) < 6 || length(unique(group[keep])) < 2) {
    writeLines("not enough survival samples", paste0(outfile, "_EMPTY.txt"))
    return(invisible(NULL))
  }
  time <- time[keep]
  event <- event[keep]
  group <- droplevels(factor(group[keep]))
  if (!has_pkg("survival")) {
    writeLines("survival package not installed", paste0(outfile, "_EMPTY.txt"))
    return(invisible(NULL))
  }
  sv <- survival::Surv(time, event)
  fit <- survival::survfit(sv ~ group)
  diff <- tryCatch(survival::survdiff(sv ~ group), error = function(e) NULL)
  pval <- if (!is.null(diff)) {
    1 - stats::pchisq(diff$chisq, length(diff$n) - 1)
  } else NA_real_
  grDevices::pdf(paste0(outfile, ".pdf"), width = 7.5, height = 6)
  graphics::plot(
    fit, col = seq_along(levels(group)), lwd = 2, xlab = "Time",
    ylab = "Survival probability", main = paste0(title, "\nlog-rank p = ",
                                                 format.pval(pval, digits = 3))
  )
  graphics::legend("bottomleft", legend = levels(group), col = seq_along(levels(group)), lwd = 2, bty = "n")
  grDevices::dev.off()
  grDevices::png(paste0(outfile, ".png"), width = 1800, height = 1400, res = 200)
  graphics::plot(
    fit, col = seq_along(levels(group)), lwd = 2, xlab = "Time",
    ylab = "Survival probability", main = paste0(title, "\nlog-rank p = ",
                                                 format.pval(pval, digits = 3))
  )
  graphics::legend("bottomleft", legend = levels(group), col = seq_along(levels(group)), lwd = 2, bty = "n")
  grDevices::dev.off()
  data.frame(title = title, logrank_p = pval, n = length(time), stringsAsFactors = FALSE)
}

cox_pathway_table <- function(scores, time, event) {
  if (!has_pkg("survival")) return(NULL)
  keep <- is.finite(time) & time > 0 & !is.na(event)
  if (sum(keep) < 8) return(NULL)
  rows <- lapply(colnames(scores), function(p) {
    x <- as.numeric(scores[keep, p])
    if (length(unique(x[is.finite(x)])) < 3) return(NULL)
    fit <- tryCatch(
      survival::coxph(survival::Surv(time[keep], event[keep]) ~ x),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    s <- summary(fit)
    data.frame(
      pathway = p,
      HR = s$coefficients[1, "exp(coef)"],
      HR_low = s$conf.int[1, "lower .95"],
      HR_high = s$conf.int[1, "upper .95"],
      cox_p = s$coefficients[1, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  if (is.null(tab) || nrow(tab) == 0) return(NULL)
  tab$cox_padj <- stats::p.adjust(tab$cox_p, method = "BH")
  tab[order(tab$cox_p), ]
}

# -----------------------------------------------------------------------------
# 可选：加载细胞系脚本里的作图 / 富集函数（不跑主流程）
# -----------------------------------------------------------------------------
load_pipeline_functions_only <- function() {
  if (exists("analyze_one_comparison", mode = "function") &&
      exists("emit_subset_analysis", mode = "function")) {
    return(TRUE)
  }
  need <- c("clusterProfiler", "org.Hs.eg.db", "ggplot2")
  if (!all(vapply(need, requireNamespace, logical(1), quietly = TRUE))) {
    log_msg("Skip loading pipeline helpers: enrichment packages not installed")
    return(FALSE)
  }
  pipe <- "TG_RNAseq_pipeline.R"
  if (!file.exists(pipe)) {
    alt <- file.path(getwd(), "TG_RNAseq_pipeline.R")
    if (file.exists(alt)) pipe <- alt
  }
  if (!file.exists(pipe)) {
    log_msg("未找到 TG_RNAseq_pipeline.R，跳过 ORA/GSEA 复用")
    return(FALSE)
  }
  lines <- readLines(pipe, warn = FALSE)
  main_at <- grep("^# 10\\. 主流程", lines)[1]
  if (is.na(main_at) || main_at < 2) return(FALSE)
  tryCatch({
    eval(parse(text = lines[seq_len(main_at - 1)]), envir = .GlobalEnv)
    TRUE
  }, error = function(e) {
    log_msg("加载 pipeline 函数失败（细胞系主流程未改）: ", e$message)
    FALSE
  })
}

export_neural_focus <- function(x, stub, title) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) return(invisible(NULL))
  df <- as.data.frame(x)
  desc <- if ("Description" %in% names(df)) df$Description else df$ID
  df$genome_wide_rank <- seq_len(nrow(df))
  df$focus_class <- vapply(seq_len(nrow(df)), function(i) {
    classify_neural_term(df$ID[i], desc[i])
  }, character(1))
  hit <- df[!is.na(df$focus_class), , drop = FALSE]
  if (nrow(hit) == 0) {
    writeLines("no neural terms in this genome-wide result", paste0(stub, "_FOCUS_neural_EMPTY.txt"))
    return(invisible(NULL))
  }
  utils::write.csv(hit, paste0(stub, "_FOCUS_neural.csv"), row.names = FALSE)
}

get_neural_term2gene <- function() {
  rows <- list()
  append_local <- function(lst, name, entrez) {
    entrez <- unique(as.character(entrez))
    entrez <- entrez[nzchar(entrez) & !is.na(entrez)]
    if (length(entrez) < 5) return(lst)
    lst[[length(lst) + 1]] <- data.frame(gs_name = name, entrez = entrez, stringsAsFactors = FALSE)
    lst
  }
  if (exists("entrez_for_go", mode = "function") && has_pkg("org.Hs.eg.db")) {
    go <- neural_go_sets()
    for (nm in names(go)) rows <- append_local(rows, nm, entrez_for_go(go[[nm]]))
    kegg <- neural_kegg_sets()
    if (exists("entrez_for_kegg_path", mode = "function")) {
      for (nm in names(kegg)) rows <- append_local(rows, nm, entrez_for_kegg_path(kegg[[nm]]))
    }
  }
  if (length(rows) == 0 && exists("map_to_entrez", mode = "function")) {
    fb <- neural_symbol_fallback()
    for (nm in names(fb)) {
      mp <- map_to_entrez(fb[[nm]])
      rows <- append_local(rows, nm, mp$entrez)
    }
  }
  if (length(rows) == 0) return(data.frame(gs_name = character(), entrez = character()))
  unique(do.call(rbind, rows))
}

run_focused_neural <- function(stats, de, heat_mat, sample_info, outdir, label) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c("这不是改全库 GO/KEGG 的 p 值或排名。",
      "本文件夹只检验神经相关信号通路基因集。",
      "全库结果旁边的 *_FOCUS_neural.csv 保留原始 p 值和 genome_wide_rank。"),
    file.path(outdir, "00_README.txt")
  )
  t2g <- get_neural_term2gene()
  if (nrow(t2g) < 8 || length(stats) < 10 || !has_pkg("clusterProfiler")) {
    writeLines("too few genes / empty neural sets / no clusterProfiler",
               file.path(outdir, "GSEA_focused_neural_EMPTY.txt"))
    return(invisible(NULL))
  }
  obj <- tryCatch(
    clusterProfiler::GSEA(
      geneList = stats, TERM2GENE = t2g, minGSSize = 5, maxGSSize = 2500,
      pvalueCutoff = 1, eps = 0, verbose = FALSE
    ),
    error = function(e) {
      log_msg("focused neural GSEA failed: ", e$message)
      NULL
    }
  )
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    obj <- tryCatch(
      clusterProfiler::setReadable(obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
      error = function(e) obj
    )
    if (exists("plot_gsea_object", mode = "function")) {
      plot_gsea_object(obj, file.path(outdir, "GSEA_focused_neural"),
                       paste(label, "| GSEA focused neural"), also_export_focus = FALSE)
    } else {
      utils::write.csv(as.data.frame(obj), file.path(outdir, "GSEA_focused_neural.csv"), row.names = FALSE)
    }
  }
  neural_genes <- unique(de$gene[de$gene %in% unlist(neural_symbol_fallback())])
  plot_patient_heatmap(
    heat_mat, sample_info, neural_genes,
    paste(label, "| neural genes"),
    file.path(outdir, "heatmap_neural_genes")
  )
}

emit_patient_subset <- function(comp_name, sub, tag, title, outdir, full_de, heat_mat, sample_info, have_helpers) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(paste("comparison:", comp_name),
      paste("subset:", tag),
      paste("p_cutoff:", p_cutoff_patient),
      paste("n_genes:", nrow(sub))),
    file.path(outdir, paste0("00_", tag, "_THIS_FOLDER.txt"))
  )
  utils::write.csv(sub, file.path(outdir, paste0(tag, "_DE_selected_genes.csv")), row.names = FALSE)
  if (nrow(sub) == 0) {
    writeLines("no genes after p < 0.01 and this cutoff", file.path(outdir, paste0(tag, "_EMPTY.txt")))
    return(invisible(NULL))
  }
  if (have_helpers && exists("emit_subset_analysis", mode = "function")) {
    gsea_cache <- list()
    tryCatch(
      emit_subset_analysis(
        comp_name, sub, tag, title, outdir, full_de,
        heat_mat, sample_info, gsea_cache, fc_line = 1
      ),
      error = function(e) log_msg("emit_subset_analysis failed: ", e$message)
    )
    return(invisible(NULL))
  }
  if (exists("plot_volcano", mode = "function")) {
    tryCatch(plot_volcano(full_de, sub$gene, title, file.path(outdir, paste0(tag, "_volcano")), fc_line = 1),
             error = function(e) log_msg("volcano failed: ", e$message))
  }
  plot_patient_heatmap(heat_mat, sample_info, sub$gene, title, file.path(outdir, paste0(tag, "_heatmap")))
}

# -----------------------------------------------------------------------------
# 主分析
# -----------------------------------------------------------------------------
run_patient_metastasis_neural <- function(project_dir = NULL, result_dir = NULL) {
  if (is.null(project_dir)) {
    if (exists("resolve_project_dir", mode = "function")) {
      project_dir <- resolve_project_dir()
    } else {
      project_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    }
  }
  if (is.null(result_dir)) {
    result_dir <- file.path(project_dir, "results", "patient_metastasis_neural")
  }
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  log_msg("Patient neural / metastasis analysis. Project: ", project_dir)

  files <- find_patient_files(project_dir)
  if (is.null(files)) {
    stop(
      "未找到病人表达/临床表。请提供 patient_expression.csv 与 patient_clinical.csv，",
      "或设置 TG_PATIENT_EXPR / TG_PATIENT_CLINICAL。",
      "格式见 examples/patient_metastasis_neural/"
    )
  }
  log_msg("Expression: ", files$expr)
  log_msg("Clinical: ", files$clinical)

  have_helpers <- load_pipeline_functions_only()
  raw <- read_expression_matrix(files$expr)
  clin <- read_table_flex(files$clinical)
  aligned <- align_patient_data(raw, clin)
  mat <- aligned$mat
  si <- aligned$clinical
  si <- si[match(colnames(mat), si$sample), , drop = FALSE]
  rownames(si) <- si$sample

  n_met <- sum(si$metastasis_group == "Metastasis", na.rm = TRUE)
  n_nom <- sum(si$metastasis_group == "NoMetastasis", na.rm = TRUE)
  log_msg("Samples: ", ncol(mat), " | Metastasis=", n_met, " NoMetastasis=", n_nom,
          " excluded/unknown=", sum(is.na(si$metastasis_group)))
  utils::write.csv(si, file.path(result_dir, "sample_clinical_used.csv"), row.names = FALSE)

  if (exists("detect_value_type", mode = "function")) {
    value_type <- detect_value_type(mat)
  } else {
    x <- as.numeric(mat)
    x <- x[is.finite(x)]
    value_type <- if (mean(abs(x - round(x)) < 1e-6) > 0.85 && stats::quantile(x, 0.95) > 50) "counts" else "fpkm"
  }
  log_msg("Value type: ", value_type)

  si$group <- si$metastasis_group
  filt <- patient_filter_low(mat, si[!is.na(si$group), , drop = FALSE], value_type)
  # 过滤用有转移标签的样品，但标准化用全部对齐样品
  si_de <- si[!is.na(si$metastasis_group), , drop = FALSE]
  filt_de <- filt[, si_de$sample, drop = FALSE]
  norm <- patient_normalize(filt_de, si_de, value_type)
  log_mat <- norm$log_mat
  heat_mat <- norm$heat_mat
  utils::write.csv(
    cbind(gene = rownames(log_mat), as.data.frame(log_mat, check.names = FALSE)),
    file.path(result_dir, "normalized_log_matrix.csv"),
    row.names = FALSE
  )

  de <- patient_de_metastasis(log_mat, si_de, raw_mat = filt_de, value_type = value_type)
  de <- de[order(de$pvalue, na.last = TRUE), ]
  utils::write.csv(de, file.path(result_dir, "DE_full_metastasis_vs_nometastasis.csv"), row.names = FALSE)
  log_msg("DE genes: ", nrow(de), " | p < ", p_cutoff_patient, " : ",
          sum(passes_p01(de$pvalue)), " | up among them: ",
          sum(passes_p01(de$pvalue) & de$log2FC > 0, na.rm = TRUE))

  gene_sets <- build_neural_gene_sets_symbols(rownames(log_mat))
  if (length(gene_sets) == 0) {
    log_msg("WARNING: no neural genes mapped in this matrix")
  } else {
    scores <- pathway_mean_z(log_mat, gene_sets)
    scores$sample <- rownames(scores)
    scores$metastasis_group <- si_de$metastasis_group[match(scores$sample, si_de$sample)]
    scores$prognosis_group <- si_de$prognosis_group[match(scores$sample, si_de$sample)]
    scores$os_time <- si_de$os_time[match(scores$sample, si_de$sample)]
    scores$os_event <- si_de$os_event[match(scores$sample, si_de$sample)]
    utils::write.csv(scores, file.path(result_dir, "neural_pathway_scores.csv"), row.names = FALSE)

    score_only <- scores[, names(gene_sets), drop = FALSE]
    met_tab <- compare_scores_by_group(score_only, scores$metastasis_group, "Metastasis", "NoMetastasis")
    utils::write.csv(met_tab, file.path(result_dir, "neural_pathway_Metastasis_vs_NoMetastasis.csv"), row.names = FALSE)
    plot_score_boxplot(
      score_only, scores$metastasis_group,
      "Neural pathway scores | Metastasis vs NoMetastasis",
      file.path(result_dir, "boxplot_neural_by_metastasis"),
      c("NoMetastasis", "Metastasis")
    )
    if (any(!is.na(scores$prognosis_group))) {
      prog_tab <- compare_scores_by_group(score_only, scores$prognosis_group, "PoorPrognosis", "GoodPrognosis")
      utils::write.csv(prog_tab, file.path(result_dir, "neural_pathway_Poor_vs_Good_prognosis.csv"), row.names = FALSE)
      plot_score_boxplot(
        score_only, scores$prognosis_group,
        "Neural pathway scores | Poor vs Good prognosis",
        file.path(result_dir, "boxplot_neural_by_prognosis"),
        c("GoodPrognosis", "PoorPrognosis")
      )
    }

    all_neural <- unique(unlist(gene_sets, use.names = FALSE))
    plot_patient_heatmap(
      heat_mat, si_de, all_neural,
      "Neural pathway genes | metastasis annotation",
      file.path(result_dir, "heatmap_neural_genes_all_samples")
    )

    km_dir <- file.path(result_dir, "Survival")
    dir.create(km_dir, recursive = TRUE, showWarnings = FALSE)
    km_rows <- list()
    km_rows[[1]] <- plot_survival_km(
      scores$os_time, scores$os_event, scores$metastasis_group,
      "OS by metastasis", file.path(km_dir, "KM_by_metastasis")
    )
    for (p in names(gene_sets)) {
      sc <- as.numeric(score_only[[p]])
      hi <- ifelse(sc >= stats::median(sc, na.rm = TRUE), "HighScore", "LowScore")
      km_rows[[length(km_rows) + 1]] <- plot_survival_km(
        scores$os_time, scores$os_event, hi,
        paste("OS by", p), file.path(km_dir, paste0("KM_", gsub("[^A-Za-z0-9]+", "_", p)))
      )
    }
    km_tab <- do.call(rbind, km_rows[!vapply(km_rows, is.null, logical(1))])
    if (!is.null(km_tab)) {
      utils::write.csv(km_tab, file.path(km_dir, "logrank_summary.csv"), row.names = FALSE)
    }
    cox <- cox_pathway_table(score_only, scores$os_time, scores$os_event)
    if (!is.null(cox)) {
      utils::write.csv(cox, file.path(km_dir, "cox_neural_pathway_scores.csv"), row.names = FALSE)
    }
  }

  # 差异分层：p < 0.01 后 FC / topN，只分析上调
  comp_name <- "Metastasis_vs_NoMetastasis"
  base <- file.path(result_dir, comp_name)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c("病人转移比较：Metastasis vs NoMetastasis。",
      "显著性：先 p < 0.01，再按上调 FC 与 topN 分层。",
      "不要只看 GSEA。FoldChange/ 与 TopRank/ 里有火山图、热图、ORA、GSEA。",
      "神经专项在 Focused_neural/；全库表旁 *_FOCUS_neural.csv 保留原始 p 与 genome_wide_rank。",
      "00_GSEA_all_genes_NOT_FC_or_topN 是全基因 GSEA，不是分层图。"),
    file.path(base, "00_READ_ME_先看这里.txt")
  )
  utils::write.csv(de, file.path(base, "DE_full.csv"), row.names = FALSE)

  heat_si <- si_de
  heat_si$group <- heat_si$metastasis_group
  for (nm in names(fc_cutoffs_patient)) {
    fc <- unname(fc_cutoffs_patient[[nm]])
    sub <- select_up_fc_p01(de, fc)
    if (nrow(sub) > 0) sub <- sub[order(sub$log2FC, decreasing = TRUE), , drop = FALSE]
    emit_patient_subset(
      comp_name, sub, nm, paste0(comp_name, " | p<0.01 up FC >= ", fc),
      file.path(base, "FoldChange", nm), de, heat_mat, heat_si, have_helpers
    )
  }
  for (n in top_ns_patient) {
    tag <- paste0("top", n)
    sub <- select_up_topn_p01(de, n)
    emit_patient_subset(
      comp_name, sub, tag, paste0(comp_name, " | p<0.01 upregulated top ", n),
      file.path(base, "TopRank", tag), de, heat_mat, heat_si, have_helpers
    )
  }

  if (have_helpers && exists("build_gsea_cache", mode = "function") &&
      exists("ranked_entrez", mode = "function")) {
    tryCatch({
      gsea_cache <- build_gsea_cache(de)
      full_gsea_dir <- file.path(base, "00_GSEA_all_genes_NOT_FC_or_topN")
      dir.create(full_gsea_dir, recursive = TRUE, showWarnings = FALSE)
      writeLines("全基因 GSEA，不是 FC/topN 分层结果。", file.path(full_gsea_dir, "00_README.txt"))
      for (nm in c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Hallmark")) {
        if (exists("plot_gsea_object", mode = "function")) {
          plot_gsea_object(
            gsea_cache[[nm]], file.path(full_gsea_dir, paste0("allGenes_GSEA_", nm)),
            paste("GSEA", nm, "|", comp_name, "| ALL genes, NOT FC/topN")
          )
        }
        export_neural_focus(
          gsea_cache[[nm]], file.path(full_gsea_dir, paste0("allGenes_GSEA_", nm)),
          paste(nm, comp_name)
        )
      }
      run_focused_neural(
        gsea_cache$stats, de, heat_mat, heat_si,
        file.path(base, "Focused_neural"),
        paste(comp_name, "| all genes")
      )
    }, error = function(e) log_msg("patient GSEA/focus failed: ", e$message))
  } else {
    log_msg("Skip genome-wide GSEA (pipeline helpers or clusterProfiler unavailable)")
    run_focused_neural(
      setNames(de$log2FC, de$gene), de, heat_mat, heat_si,
      file.path(base, "Focused_neural"),
      paste(comp_name, "| all genes")
    )
  }

  writeLines(
    c("病人转移 / 预后 × 神经通路分析完成。",
      paste("结果目录:", result_dir),
      "先看: neural_pathway_Metastasis_vs_NoMetastasis.csv",
      "     boxplot_neural_by_metastasis.pdf",
      "     Survival/KM_by_metastasis.pdf",
      "     Metastasis_vs_NoMetastasis/Focused_neural/"),
    file.path(result_dir, "00_README.txt")
  )
  log_msg("Patient neural / metastasis analysis done: ", result_dir)
  invisible(list(de = de, clinical = si_de, result_dir = result_dir))
}

if (!isTRUE(getOption("tg.patient.skip_main", FALSE))) {
  run_patient_metastasis_neural()
}
