#!/usr/bin/env Rscript
# =============================================================================
# 乳腺癌多器官转移 bulk 转录组/蛋白组
# 问题1：肿瘤高表达哪些基因可能促神经浸润
# 问题2：神经浸润更倾向哪个继发部位
#
# 主数据：GSE175692（184 例，11 器官）；补充 GSE14020 芯片
# 独立于 TG_RNAseq_*.R / ST / scRNA 脚本
#
#   setwd("E:/R/BRCA_met_bulk_nerve")
#   source("BRCA_met_bulk_nerve_infiltration.R")
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")

need_cran <- c("dplyr", "tidyr", "tibble", "stringr", "ggplot2", "matrixStats")
need_bioc <- c("limma", "org.Hs.eg.db")
opt_cran <- c("writexl", "pheatmap", "preprocessCore")
opt_bioc <- c("GEOquery", "Biobase", "clusterProfiler", "enrichplot", "fgsea",
              "msigdbr", "DESeq2")

install_if_missing <- function(pkgs, bioc = FALSE, required = FALSE) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) == 0) return(invisible(TRUE))
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      tryCatch(install.packages("BiocManager", repos = "https://cloud.r-project.org"),
               error = function(e) NULL)
    }
    tryCatch(BiocManager::install(miss, update = FALSE, ask = FALSE),
             error = function(e) message("Bioc install failed: ", e$message))
  } else {
    tryCatch(install.packages(miss, repos = "https://cloud.r-project.org"),
             error = function(e) message("CRAN install failed: ", e$message))
  }
  still <- miss[!vapply(miss, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still) > 0 && required) stop("缺少 R 包: ", paste(still, collapse = ", "))
  if (length(still) > 0) message("可选包未装: ", paste(still, collapse = ", "))
  invisible(TRUE)
}
install_if_missing(need_cran, FALSE, TRUE)
install_if_missing(opt_cran, FALSE, FALSE)
install_if_missing(need_bioc, TRUE, FALSE)
install_if_missing(opt_bioc, TRUE, FALSE)
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
for (p in c(need_cran, opt_cran, need_bioc, opt_bioc)) {
  if (has_pkg(p)) suppressPackageStartupMessages(library(p, character.only = TRUE))
}

resolve_dir <- function() {
  env_dir <- Sys.getenv("BRCA_MET_BULK_DIR", unset = "")
  cands <- unique(c(env_dir, "E:/R/BRCA_met_bulk_nerve", "E:\\R\\BRCA_met_bulk_nerve", getwd()))
  for (d in cands[nzchar(cands)]) {
    if (dir.exists(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

met_dir <- resolve_dir()
result_dir <- file.path(met_dir, "results")
log_dir <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
log_msg("Bulk data directory: ", met_dir)

p_cutoff <- 0.01
fc_cutoffs_local <- c("FC_1.25" = 1.25, "FC_1.5" = 1.5, "FC_2" = 2)

ligand_genes <- c("NGF", "BDNF", "NTF3", "NTF4", "GDNF", "ARTN", "NRTN", "PSPN",
                  "NTN1", "NTN4", "SLIT2", "SLIT3", "SEMA3A", "SEMA3C", "SEMA3F",
                  "EFNA1", "EFNA5", "EFNB2", "NRG1", "CXCL12", "VEGFA", "MDK", "PTN",
                  "TGFB1", "WNT5A", "L1CAM", "NCAM1", "CCL2")
stroma_genes <- c("GFAP", "MBP", "MPZ", "PMP22", "S100B", "SOX10", "PLP1", "MOG",
                  "AQP4", "OLIG2", "RBFOX3", "MAP2", "SYP", "MAG")
nerve_term_pat <- paste("axon", "neur", "nerve", "schwann", "synap", "neurotroph",
                        "semaphorin", "ephrin", "glia", "oligodend", "astrocyte", sep = "|")

this_script_dir <- function() {
  ofile <- NULL
  if (sys.nframe() > 0) {
    for (i in sys.nframe():1) {
      e <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
      if (!is.null(e)) { ofile <- e; break }
    }
  }
  if (is.null(ofile)) getwd() else dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE))
}
load_pipeline_functions_only <- function() {
  cands <- c(file.path(this_script_dir(), "TG_RNAseq_pipeline.R"), "TG_RNAseq_pipeline.R",
             "E:/R/TG_BRCA/TG/TG_RNAseq_pipeline.R")
  pipe <- cands[file.exists(cands)][1]
  if (is.na(pipe) || !nzchar(pipe)) return(FALSE)
  lines <- readLines(pipe, warn = FALSE)
  main_at <- grep("^# 10\\. 主流程", lines)[1]
  if (is.na(main_at) || main_at < 2) return(FALSE)
  eval(parse(text = lines[seq_len(main_at - 1)]), envir = .GlobalEnv)
  TRUE
}
have_pipeline <- FALSE
if (!identical(Sys.getenv("BRCA_MET_BULK_SKIP_PIPELINE"), "1")) {
  have_pipeline <- tryCatch(load_pipeline_functions_only(), error = function(e) {
    log_msg("未载入 TG 作图函数: ", e$message); FALSE
  })
}

norm_organ <- function(x) {
  x <- tolower(paste(x, collapse = " "))
  if (grepl("brain|cns|cerebr", x)) return("brain")
  if (grepl("lung|pulmon", x)) return("lung")
  if (grepl("bone|osseo|marrow|femur|spine|vertebra|osseous", x)) return("bone")
  if (grepl("liver|hepat", x)) return("liver")
  if (grepl("skin|derm|cutan", x)) return("skin")
  if (grepl("ovary|ovarian", x)) return("ovary")
  if (grepl("pleura", x)) return("pleura")
  if (grepl("node|lymph|axilla", x)) return("lymph_node")
  if (grepl("breast|chest wall|chestwall", x)) return("breast")
  if (grepl("soft tissue|soft_tissue", x)) return("soft_tissue")
  "other"
}

present_genes <- function(genes, universe) intersect(unique(genes), universe)

zmean_score <- function(mat, genes) {
  g <- present_genes(genes, rownames(mat))
  if (length(g) < 2) return(rep(NA_real_, ncol(mat)))
  z <- t(scale(t(mat[g, , drop = FALSE])))
  z[is.na(z)] <- 0
  colMeans(z)
}

quantile_norm <- function(mat) {
  if (has_pkg("preprocessCore")) {
    out <- preprocessCore::normalize.quantiles(mat)
    dimnames(out) <- dimnames(mat)
    return(out)
  }
  # 无 preprocessCore 时：各列秩次再反 qnorm
  apply(mat, 2, function(v) {
    r <- rank(v, ties.method = "average", na.last = "keep")
    stats::qnorm((r - 0.5) / sum(is.finite(r)))
  })
}

collapse_by_symbol <- function(mat, symbols) {
  symbols <- as.character(symbols)
  symbols[is.na(symbols) | symbols %in% c("", "---", "NA")] <- NA
  keep <- !is.na(symbols)
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  if (any(duplicated(symbols))) {
    ord <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
    mat <- mat[ord, , drop = FALSE]
    symbols <- symbols[ord]
    keep_u <- !duplicated(symbols)
    mat <- mat[keep_u, , drop = FALSE]
    symbols <- symbols[keep_u]
  }
  rownames(mat) <- symbols
  mat
}

# -----------------------------------------------------------------------------
# 读入：本地矩阵 / GEO / DIA-NN
# -----------------------------------------------------------------------------
read_local_matrix <- function() {
  pg <- list.files(met_dir, pattern = "report\\.pg_matrix", recursive = TRUE, full.names = TRUE)
  if (length(pg) > 0) {
    log_msg("DIA-NN: ", pg[1])
    df <- utils::read.delim(pg[1], check.names = FALSE, stringsAsFactors = FALSE)
    gene_col <- intersect(c("Genes", "Gene", "Gene.Names", "Protein.Names"), names(df))[1]
    num <- vapply(df, is.numeric, logical(1))
    mat <- as.matrix(df[, num, drop = FALSE])
    mat[mat <= 0] <- NA
    mat <- log2(mat)
    symbols <- if (!is.na(gene_col)) df[[gene_col]] else rownames(df)
    symbols <- sub(";.*", "", as.character(symbols))
    mat <- collapse_by_symbol(mat, symbols)
    meta_f <- list.files(met_dir, pattern = "metadata|sample.*organ|pdata", ignore.case = TRUE,
                         recursive = TRUE, full.names = TRUE)
    si <- data.frame(sample = colnames(mat), organ = "other", stringsAsFactors = FALSE)
    if (length(meta_f) > 0) {
      md <- tryCatch(utils::read.csv(meta_f[1], stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(md)) {
        sc <- intersect(c("sample", "Sample", "id"), names(md))[1]
        oc <- intersect(c("organ", "Organ", "site", "Site"), names(md))[1]
        if (!is.na(sc) && !is.na(oc)) {
          si$organ <- vapply(as.character(md[[oc]][match(si$sample, md[[sc]])]),
                             norm_organ, character(1))
        }
      }
    }
    return(list(name = "DIA_NN_protein", mat = mat, si = si, platform = "protein"))
  }
  NULL
}

pdata_blob <- function(pd) {
  cols <- names(pd)[vapply(pd, function(x) is.character(x) || is.factor(x), logical(1))]
  apply(pd[, cols, drop = FALSE], 1, function(r) paste(r, collapse = " | "))
}

eset_to_dataset <- function(eset, name) {
  if (!has_pkg("Biobase")) return(NULL)
  mat <- as.matrix(Biobase::exprs(eset))
  pd <- as.data.frame(Biobase::pData(eset))
  fd <- tryCatch(as.data.frame(Biobase::fData(eset)), error = function(e) NULL)
  symbols <- rownames(mat)
  if (!is.null(fd) && nrow(fd) == nrow(mat)) {
    for (cc in c("GENE_SYMBOL", "Gene Symbol", "gene_assignment", "Symbol", "gene_symbol", "ID")) {
      if (cc %in% names(fd)) {
        symbols <- as.character(fd[[cc]])
        symbols <- sub(" // .*", "", symbols)
        symbols <- sub(";.*", "", symbols)
        break
      }
    }
  }
  storage.mode(mat) <- "double"
  if (max(mat, na.rm = TRUE) > 100) mat <- log2(mat + 1)
  mat <- collapse_by_symbol(mat, symbols)
  blob <- pdata_blob(pd)
  organ <- vapply(blob, norm_organ, character(1))
  si <- data.frame(sample = colnames(mat), organ = organ,
                   geo = if ("geo_accession" %in% names(pd)) as.character(pd$geo_accession) else colnames(mat),
                   stringsAsFactors = FALSE)
  rownames(si) <- si$sample
  list(name = name, mat = mat, si = si, platform = "array_or_panel")
}

download_geo <- function(acc) {
  if (!has_pkg("GEOquery")) {
    log_msg("未安装 GEOquery，跳过 ", acc)
    return(NULL)
  }
  dest <- file.path(met_dir, "GEO")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  log_msg("GEOquery ", acc)
  g <- tryCatch(GEOquery::getGEO(acc, destdir = dest, getGPL = TRUE),
                error = function(e) { log_msg("getGEO failed ", acc, ": ", e$message); NULL })
  if (is.null(g) || length(g) == 0) return(NULL)
  out <- list()
  for (i in seq_along(g)) {
    nm <- paste0(acc, if (length(g) > 1) paste0("_", i) else "")
    ds <- eset_to_dataset(g[[i]], nm)
    if (!is.null(ds)) out[[nm]] <- ds
  }
  out
}

datasets <- list()
local <- read_local_matrix()
if (!is.null(local)) datasets[[local$name]] <- local

geo_acc <- c("GSE175692", "GSE14017", "GSE14018")
if (!identical(Sys.getenv("BRCA_MET_BULK_SKIP_GEO"), "1")) {
  for (acc in geo_acc) {
    got <- download_geo(acc)
    for (nm in names(got)) datasets[[nm]] <- got[[nm]]
  }
}

if (length(datasets) == 0) {
  log_msg("没有读到表达矩阵。请看 BRCA_met_bulk_nerve_DOWNLOAD.txt，或安装 GEOquery 后联网运行。")
}

# -----------------------------------------------------------------------------
# 预处理 + 分数 + DEG
# -----------------------------------------------------------------------------
empty_de <- function() {
  data.frame(gene = character(), log2FC = numeric(), AveExpr = numeric(),
             pvalue = numeric(), padj = numeric(), padj_BH = numeric(),
             stringsAsFactors = FALSE)
}

limma_near_far <- function(mat, group, label = "") {
  group <- factor(group, levels = c("far", "near"))
  if (length(unique(group)) < 2 || ncol(mat) < 4) return(empty_de())
  design <- stats::model.matrix(~ group)
  fit <- tryCatch({
    if (!has_pkg("limma")) return(NULL)
    limma::eBayes(limma::lmFit(mat, design), trend = TRUE, robust = TRUE)
  }, error = function(e) { log_msg("limma failed ", label, ": ", e$message); NULL })
  if (is.null(fit)) return(empty_de())
  coefn <- grep("groupnear", colnames(design), value = TRUE)
  tt <- limma::topTable(fit, coef = coefn, number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), log2FC = tt$logFC, AveExpr = tt$AveExpr,
             pvalue = tt$P.Value, padj = tt$P.Value, padj_BH = tt$adj.P.Val,
             stringsAsFactors = FALSE)
}

basic_volcano <- function(de, title, outfile, fc_line = 1.25) {
  if (nrow(de) == 0) return(invisible(NULL))
  df <- de
  df$y <- if (any(!is.na(df$pvalue))) -log10(pmax(df$pvalue, 1e-300)) else abs(df$log2FC)
  df$col <- "ns"
  df$col[!is.na(df$pvalue) & df$pvalue < p_cutoff & df$log2FC >= log2(fc_line)] <- "up"
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  if (has_pkg("ggplot2")) {
    p <- ggplot2::ggplot(df, ggplot2::aes(log2FC, y, color = col)) +
      ggplot2::geom_point(alpha = 0.45, size = 0.8) +
      ggplot2::scale_color_manual(values = c(ns = "grey70", up = "#D62828")) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(title = title, y = "-log10(p)")
    ggplot2::ggsave(paste0(outfile, ".pdf"), p, width = 7, height = 6)
    ggplot2::ggsave(paste0(outfile, ".png"), p, width = 7, height = 6, dpi = 140)
  }
}

extract_focus_nerve <- function(outdir, tag) {
  csvs <- list.files(outdir, pattern = "ORA_.*\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(csvs) == 0) return(invisible(NULL))
  fdir <- file.path(outdir, "Focused_nerve_invasion")
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  chunks <- list()
  for (f in csvs) {
    df <- tryCatch(utils::read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    desc <- if ("Description" %in% names(df)) df$Description else df[[1]]
    hit <- grepl(nerve_term_pat, desc, ignore.case = TRUE)
    if (!any(hit)) next
    sub <- df[hit, , drop = FALSE]
    pcol <- intersect(c("pvalue", "p.adjust", "pval", "padj"), names(df))[1]
    if (length(pcol) == 1 && !is.na(pcol)) {
      ord <- order(df[[pcol]])
      sub$genome_wide_rank <- match(which(hit), ord)
    } else sub$genome_wide_rank <- seq_len(nrow(sub))
    sub$source_file <- basename(f)
    chunks[[length(chunks) + 1]] <- sub
  }
  if (length(chunks) == 0) return(invisible(NULL))
  common <- Reduce(intersect, lapply(chunks, names))
  hit <- do.call(rbind, lapply(chunks, function(d) d[, common, drop = FALSE]))
  utils::write.csv(hit, file.path(fdir, paste0(tag, "_FOCUS_nerve_invasion.csv")), row.names = FALSE)
}

emit_comparison <- function(comp_name, de, heat_mat, sample_info) {
  de <- de[!is.na(de$log2FC), , drop = FALSE]
  de$padj <- de$pvalue
  base <- file.path(result_dir, comp_name)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(de, file.path(base, "DE_full.csv"), row.names = FALSE)
  have_p <- any(!is.na(de$pvalue))
  log_msg(comp_name, " genes=", nrow(de), " have_p=", have_p)
  writeLines(
    c("bulk 转移灶：上调 FC >= 1.25 / 1.5 / 2（先 p < 0.01）。",
      "不做 FC=1，不做 TopRank。脑高 GFAP 不等于 PNI。"),
    file.path(base, "00_READ_ME_先看这里.txt")
  )
  invisible(lapply(file.path(base, "FoldChange", names(fc_cutoffs_local)),
                   dir.create, recursive = TRUE, showWarnings = FALSE))
  gsea_cache <- list()
  use_pipe <- have_pipeline && exists("emit_subset_analysis", mode = "function")
  assign("result_dir", result_dir, envir = .GlobalEnv)
  assign("padj_cutoff", p_cutoff, envir = .GlobalEnv)
  for (nm in names(fc_cutoffs_local)) {
    fc <- unname(fc_cutoffs_local[[nm]])
    keep <- !is.na(de$log2FC) & de$log2FC > 0 & (2^de$log2FC >= fc)
    if (have_p) keep <- keep & !is.na(de$pvalue) & de$pvalue < p_cutoff
    sub <- de[keep, , drop = FALSE]
    if (nrow(sub) > 0) sub <- sub[order(sub$log2FC, decreasing = TRUE), , drop = FALSE]
    od <- file.path(base, "FoldChange", nm)
    if (use_pipe) {
      tryCatch(
        emit_subset_analysis(comp_name, sub, nm, paste0(comp_name, " | up FC >= ", fc),
                             od, de, heat_mat, sample_info, gsea_cache, fc_line = fc),
        error = function(e) log_msg("ERROR subset ", e$message)
      )
    } else {
      utils::write.csv(sub, file.path(od, paste0(nm, "_DE_selected_genes.csv")), row.names = FALSE)
      basic_volcano(de, paste(comp_name, nm), file.path(od, paste0(nm, "_volcano")), fc)
    }
    tryCatch(extract_focus_nerve(od, nm), error = function(e) NULL)
  }
  invisible(de)
}

all_site <- list()
all_cand_de <- list()

process_dataset <- function(ds) {
  mat0 <- ds$mat
  si <- ds$si
  common <- intersect(colnames(mat0), si$sample)
  if (length(common) < 6) {
    log_msg(ds$name, " 样本太少"); return(invisible(NULL))
  }
  mat0 <- mat0[, common, drop = FALSE]
  si <- si[match(common, si$sample), , drop = FALSE]
  keep <- rowSums(is.finite(mat0)) >= max(4, floor(0.3 * ncol(mat0)))
  mat0 <- mat0[keep, , drop = FALSE]
  mat0[!is.finite(mat0)] <- min(mat0[is.finite(mat0)], na.rm = TRUE)
  mat <- quantile_norm(mat0)
  log_msg(ds$name, " genes=", nrow(mat), " samples=", ncol(mat),
          " organs=", paste(sprintf("%s=%s", names(table(si$organ)), as.integer(table(si$organ))),
                            collapse = ", "))
  log_msg(ds$name, " ligand genes present: ",
          paste(present_genes(ligand_genes, rownames(mat)), collapse = ","))
  log_msg(ds$name, " stroma genes present: ",
          paste(present_genes(stroma_genes, rownames(mat)), collapse = ","))

  si$score_ligand <- zmean_score(mat, ligand_genes)
  si$score_stroma <- zmean_score(mat, stroma_genes)
  ok <- is.finite(si$score_ligand) & is.finite(si$score_stroma)
  si$score_ligand_resid <- NA_real_
  if (sum(ok) >= 8) {
    fit <- stats::lm(score_ligand ~ score_stroma, data = si[ok, ])
    si$score_ligand_resid[ok] <- stats::resid(fit)
  }
  utils::write.csv(si, file.path(log_dir, paste0(ds$name, "_sample_scores.csv")), row.names = FALSE)

  organs <- names(which(table(si$organ) >= 3))
  rank_rows <- lapply(organs, function(org) {
    in_o <- si$organ == org
    rest <- !in_o
    wl <- tryCatch(stats::wilcox.test(si$score_ligand[in_o], si$score_ligand[rest])$p.value,
                   error = function(e) NA_real_)
    ws <- tryCatch(stats::wilcox.test(si$score_stroma[in_o], si$score_stroma[rest])$p.value,
                   error = function(e) NA_real_)
    wr <- tryCatch(stats::wilcox.test(si$score_ligand_resid[in_o], si$score_ligand_resid[rest])$p.value,
                   error = function(e) NA_real_)
    data.frame(
      dataset = ds$name, organ = org, n = sum(in_o),
      mean_ligand = mean(si$score_ligand[in_o], na.rm = TRUE),
      mean_stroma = mean(si$score_stroma[in_o], na.rm = TRUE),
      mean_ligand_resid = mean(si$score_ligand_resid[in_o], na.rm = TRUE),
      p_ligand_vs_rest = wl, p_stroma_vs_rest = ws, p_resid_vs_rest = wr,
      stringsAsFactors = FALSE
    )
  })
  site <- do.call(rbind, rank_rows)
  if (!is.null(site) && nrow(site) > 0) {
    site <- site[order(-site$mean_ligand_resid, -site$mean_ligand), ]
    site$rank_by_ligand_resid <- seq_len(nrow(site))
    all_site[[ds$name]] <<- site
    utils::write.csv(site, file.path(result_dir, paste0("02_SITE_RANK_", ds$name, ".csv")),
                     row.names = FALSE)
    if (has_pkg("ggplot2")) {
      p <- ggplot2::ggplot(si[si$organ %in% organs, ],
                           ggplot2::aes(x = organ, y = score_ligand_resid, fill = organ)) +
        ggplot2::geom_boxplot(outlier.size = 0.6, na.rm = TRUE) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none") +
        ggplot2::labs(title = paste(ds$name, "| neurotropic ligand score residualized on neural stroma"),
                      y = "ligand | stroma (higher = more tumor neurotropic program)",
                      subtitle = "Brain high GFAP/MBP is stroma, not automatically PNI")
      ggplot2::ggsave(file.path(result_dir, paste0("02_SITE_RANK_", ds$name, "_ligand_resid_boxplot.pdf")),
                      p, width = 8, height = 5)
      ggplot2::ggsave(file.path(result_dir, paste0("02_SITE_RANK_", ds$name, "_ligand_resid_boxplot.png")),
                      p, width = 8, height = 5, dpi = 140)
    }
  }

  heat <- mat
  si_heat <- data.frame(sample = si$sample, group = si$organ, stringsAsFactors = FALSE)
  # 各器官 vs 其余
  for (org in intersect(c("brain", "lung", "bone", "liver"), unique(si$organ))) {
    if (sum(si$organ == org) < 3 || sum(si$organ != org) < 3) next
    grp <- ifelse(si$organ == org, "near", "far")
    de <- limma_near_far(mat, grp, paste(ds$name, org, "vs rest"))
    emit_comparison(paste0(ds$name, "_", org, "_vs_other_mets"), de, heat, si_heat)
    all_cand_de[[paste(ds$name, org)]] <<- list(de = de, organ = org, dataset = ds$name)
  }
  pairs <- list(c("brain", "lung"), c("brain", "bone"), c("lung", "bone"))
  for (pr in pairs) {
    a <- pr[1]; b <- pr[2]
    keep <- si$organ %in% c(a, b)
    if (sum(si$organ == a) < 3 || sum(si$organ == b) < 3) next
    grp <- ifelse(si$organ[keep] == a, "near", "far")
    de <- limma_near_far(mat[, keep, drop = FALSE], grp, paste(ds$name, a, "vs", b))
    si2 <- si_heat[keep, , drop = FALSE]
    emit_comparison(paste0(ds$name, "_", a, "_vs_", b), de, heat[, keep, drop = FALSE], si2)
  }
}

for (nm in names(datasets)) {
  tryCatch(process_dataset(datasets[[nm]]),
           error = function(e) log_msg("dataset failed ", nm, ": ", e$message))
}

# -----------------------------------------------------------------------------
# 总表
# -----------------------------------------------------------------------------
if (length(all_site) > 0) {
  site_all <- do.call(rbind, all_site)
  utils::write.csv(site_all, file.path(result_dir, "02_SITE_RANK_neural_invasion.csv"), row.names = FALSE)
}

cand <- data.frame(gene = ligand_genes, class = "neurotropic_ligand", stringsAsFactors = FALSE)
for (k in names(all_cand_de)) {
  x <- all_cand_de[[k]]
  tag <- paste(x$dataset, x$organ, sep = "__")
  cand[[paste0(tag, "_log2FC")]] <- x$de$log2FC[match(cand$gene, x$de$gene)]
  cand[[paste0(tag, "_pvalue")]] <- x$de$pvalue[match(cand$gene, x$de$gene)]
}
fc_cols <- grep("_log2FC$", names(cand), value = TRUE)
pv_cols <- grep("_pvalue$", names(cand), value = TRUE)
if (length(fc_cols) > 0) {
  cand$n_up_p01 <- 0
  for (i in seq_along(fc_cols)) {
    fc <- cand[[fc_cols[i]]]
    pv <- if (i <= length(pv_cols)) cand[[pv_cols[i]]] else rep(NA_real_, nrow(cand))
    cand$n_up_p01 <- cand$n_up_p01 + as.integer(!is.na(fc) & fc > 0 & !is.na(pv) & pv < p_cutoff)
  }
  cand <- cand[order(-cand$n_up_p01), ]
}
utils::write.csv(cand, file.path(result_dir, "01_CANDIDATE_MOLECULES_tumor_to_nerve.csv"), row.names = FALSE)
if (has_pkg("writexl")) {
  tryCatch(writexl::write_xlsx(cand, file.path(result_dir, "01_CANDIDATE_MOLECULES_tumor_to_nerve.xlsx")),
           error = function(e) NULL)
}

protocol <- c(
  "============================================================",
  "乳腺癌转移 bulk：神经浸润基因 + 继发部位倾向",
  "数据目录: ", met_dir,
  "============================================================",
  "",
  "主队列 GSE175692（184 例，11 器官，771 基因）。全转录组补充 GSE14020。",
  "蛋白组需自备 DIA-NN；公开多器官转移蛋白组几乎没有。",
  "两个平台不要拼成一张矩阵。小鼠 PDX 不当病人。",
  "",
  "【问题 1】肿瘤高表达哪些基因可能促神经浸润",
  "  看 01_CANDIDATE_MOLECULES_tumor_to_nerve.csv",
  "  以及各 GSE*_brain|lung|bone_vs_other_mets / FoldChange/FC_1.25 等",
  "  先 p < 0.01，再上调 FC >= 1.25 / 1.5 / 2",
  "",
  "【问题 2】更倾向哪个继发部位",
  "  看 02_SITE_RANK_neural_invasion.csv",
  "  排序请用 mean_ligand_resid（配体分数对神经基质回归后的残差）",
  "  不要只用 mean_stroma：脑转移 GFAP/MBP 高常常是脑组织，不是外周 PNI",
  "  脑实质转移 ≠ 沿神经浸润（PNI）",
  "",
  "Focused_nerve_invasion/ 保留原始 p 与 genome_wide_rank。"
)
writeLines(protocol, file.path(result_dir, "00_PROTOCOL_神经浸润与继发部位.txt"))
log_msg("Done. ", file.path(result_dir, "02_SITE_RANK_neural_invasion.csv"))
