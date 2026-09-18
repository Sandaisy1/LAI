#!/usr/bin/env Rscript
# =============================================================================
# E:\R\Nerve RNA —— 原位（原发）乳腺癌：肺转移倾向 vs 其他器官转移倾向
#
# 问题1：倾向肺转移的原位乳腺癌组织 vs 倾向其他器官转移的原位乳腺癌组织
#         → 下调基因（仅 FC < 1；不做 1.25/1.5/2 分层）+ 火山图 + Excel
# 问题2：在「倾向肺转移」的原位乳腺癌组织中，与神经浸润负相关的蛋白有哪些？
#         （公开队列几乎无匹配蛋白组；默认用转录本丰度作蛋白代理，并支持自备蛋白矩阵）
#
# 【对应关系硬性要求】
#   肺转移 ↔ 原位灶 必须是「同一患者临床随访结局」对应：
#   - lung_tropic：该原发灶患者随访中确实发生肺转移（lm event / lung met）
#   - other_met：该原发灶患者发生转移但未发生肺转移
#   禁止用 LMS 签名预测、禁止用转移灶活检冒充原位灶、禁止无结局的样本入组。
#
# 主数据（原发灶 + 肺/骨转移结局，GPL96）：
#   GSE2603（MSK 原发灶，推荐）
#   GSE5327（Erasmus ER- 原发灶，验证）
# 可选蛋白组：proteomics_pg_matrix.csv + proteomics_sample_meta.csv
#
#   setwd("E:/R/Nerve RNA")
#   source("Nerve_RNA_primary_lung_tropism.R")
#
# 独立于 TG_RNAseq_*.R / GSE175692 转移灶脚本，不修改原流程。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")

need_cran <- c("ggplot2", "matrixStats")
opt_cran  <- c("writexl", "preprocessCore", "pheatmap")
need_bioc <- c("limma")
opt_bioc  <- c("GEOquery", "Biobase")

install_if_missing <- function(pkgs, bioc = FALSE, required = FALSE) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) == 0) return(invisible(TRUE))
  writable <- Filter(function(p) dir.exists(p) && file.access(p, 2) == 0, .libPaths())
  if (length(writable) == 0) {
    user_lib <- path.expand(file.path("~", "R", "library"))
    dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
    .libPaths(c(user_lib, .libPaths()))
    writable <- user_lib
  }
  lib <- writable[1]
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      tryCatch(install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = lib),
               error = function(e) NULL)
    }
    tryCatch(BiocManager::install(miss, update = FALSE, ask = FALSE, lib = lib),
             error = function(e) message("Bioc install failed: ", e$message))
  } else {
    tryCatch(install.packages(miss, repos = "https://cloud.r-project.org", lib = lib),
             error = function(e) message("CRAN install failed: ", e$message))
  }
  still <- miss[!vapply(miss, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still) > 0 && required) stop("缺少 R 包: ", paste(still, collapse = ", "))
  if (length(still) > 0) message("可选包未装，相关步骤跳过: ", paste(still, collapse = ", "))
  invisible(TRUE)
}
install_if_missing(need_cran, FALSE, TRUE)
install_if_missing(opt_cran, FALSE, FALSE)
install_if_missing(need_bioc, TRUE, TRUE)
install_if_missing(opt_bioc, TRUE, FALSE)
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
for (p in c(need_cran, opt_cran, need_bioc, opt_bioc)) {
  if (has_pkg(p)) suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# -----------------------------------------------------------------------------
# 0. 路径
# -----------------------------------------------------------------------------
resolve_dir <- function() {
  env_dir <- Sys.getenv("NERVE_RNA_DIR", unset = "")
  cands <- unique(c(env_dir, "E:/R/Nerve RNA", "E:\\R\\Nerve RNA",
                    "E:/R/Nerve_RNA", getwd()))
  for (d in cands[nzchar(cands)]) {
    if (!dir.exists(d)) next
    hits <- list.files(d, pattern = "GSE2603|GSE5327|GSE175692|proteomics",
                       ignore.case = TRUE)
    if (length(hits) > 0 || identical(normalizePath(d, winslash = "/", mustWork = FALSE),
                                      normalizePath(getwd(), winslash = "/", mustWork = FALSE))) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  stop("找不到数据目录。请创建 E:/R/Nerve RNA 并放入 GSE2603/GSE5327 series matrix。")
}

data_dir   <- resolve_dir()
result_dir <- file.path(data_dir, "results_primary_lung_tropism")
log_dir    <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("primary_lung_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
log_msg("Data directory: ", data_dir)

p_cutoff <- 0.01
# 问题1：只讨论下调 FC < 1（即 log2FC < 0）；不做 1.25/1.5/2 分层
fc_down_max <- 1  # fold_change < 1

# 神经浸润签名（用于问题2 分数；转录本 / 蛋白共用符号）
nerve_ligand_genes <- c(
  "BDNF", "NGFR", "NTRK1", "NTRK2", "NTRK3", "CXCL12", "VEGFA", "VEGFD",
  "TGFB1", "TGFB2", "TGFB3", "WNT5A", "WNT5B", "NCAM1", "HGF", "MET",
  "ERBB2", "CCL2", "CDH2", "VIM", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2",
  "MMP9", "MMP14", "ITGB1", "L1CAM", "NRG1", "ARTN", "GDNF", "NGF"
)

# -----------------------------------------------------------------------------
# 1. 工具
# -----------------------------------------------------------------------------
find_file <- function(patterns) {
  hits <- character()
  for (pat in patterns) {
    hits <- c(hits, list.files(data_dir, pattern = pat, recursive = TRUE,
                               full.names = TRUE, ignore.case = TRUE))
  }
  unique(hits[file.exists(hits)])
}

split_quoted <- function(line) {
  parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
  gsub('^"|"$', "", parts[-1])
}

quantile_norm <- function(mat) {
  if (has_pkg("preprocessCore")) {
    out <- preprocessCore::normalize.quantiles(as.matrix(mat))
    dimnames(out) <- dimnames(mat)
    return(out)
  }
  apply(mat, 2, function(v) {
    r <- rank(v, ties.method = "average", na.last = "keep")
    stats::qnorm((r - 0.5) / sum(is.finite(r)))
  })
}

zmean_score <- function(mat, genes) {
  g <- intersect(unique(genes), rownames(mat))
  if (length(g) == 0) return(rep(NA_real_, ncol(mat)))
  if (length(g) == 1) {
    z <- scale(as.numeric(mat[g, ]))[, 1]
    z[is.na(z)] <- 0
    return(as.numeric(z))
  }
  z <- t(scale(t(mat[g, , drop = FALSE])))
  z[is.na(z)] <- 0
  colMeans(z)
}

collapse_probes <- function(mat, probe2gene) {
  gene <- probe2gene[rownames(mat)]
  keep <- !is.na(gene) & nzchar(gene) & gene != "---"
  mat <- mat[keep, , drop = FALSE]
  gene <- gene[keep]
  # 多基因符号取第一个
  gene <- sub(" /// .*", "", gene)
  gene <- sub(";.*", "", gene)
  gene <- trimws(gene)
  ok <- nzchar(gene)
  mat <- mat[ok, , drop = FALSE]
  gene <- gene[ok]
  # 每基因取探针中位
  spl <- split(seq_len(nrow(mat)), gene)
  out <- matrix(NA_real_, nrow = length(spl), ncol = ncol(mat),
                dimnames = list(names(spl), colnames(mat)))
  for (gn in names(spl)) {
    idx <- spl[[gn]]
    if (length(idx) == 1L) out[gn, ] <- mat[idx, ]
    else out[gn, ] <- matrixStats::colMedians(mat[idx, , drop = FALSE], na.rm = TRUE)
  }
  out
}

write_xlsx_safe <- function(obj, path) {
  if (!has_pkg("writexl")) return(invisible(FALSE))
  tryCatch({ writexl::write_xlsx(obj, path); TRUE }, error = function(e) {
    log_msg("xlsx 写出失败: ", e$message); FALSE
  })
}

save_gg <- function(p, stub, w = 7, h = 5) {
  dir.create(dirname(stub), recursive = TRUE, showWarnings = FALSE)
  if (!has_pkg("ggplot2") || is.null(p)) return(invisible(NULL))
  ggplot2::ggsave(paste0(stub, ".pdf"), p, width = w, height = h)
  ggplot2::ggsave(paste0(stub, ".png"), p, width = w, height = h, dpi = 140)
}

# 火山图：高亮 p < 0.01 且 FC < 1（下调）
volcano_down <- function(de, title, stub) {
  if (nrow(de) == 0) return(invisible(NULL))
  df <- de
  df$y <- -log10(pmax(df$pvalue, 1e-300))
  df$col <- "ns"
  df$col[!is.na(df$pvalue) & df$pvalue < p_cutoff & df$fold_change < fc_down_max] <- "down"
  df$col[!is.na(df$pvalue) & df$pvalue < p_cutoff & df$fold_change > 1] <- "up"
  p <- ggplot2::ggplot(df, ggplot2::aes(log2FC, y, color = col)) +
    ggplot2::geom_point(alpha = 0.55, size = 1.1) +
    ggplot2::scale_color_manual(values = c(ns = "grey70", down = "#1D3557", up = "#D62828")) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = 2, color = "grey40") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("DOWN highlight: p < ", p_cutoff, " & FC < ", fc_down_max,
                        " (lung-tropic vs other-met primary)"),
      y = "-log10(p)", x = "log2FC (lung_tropic / other_organ_met)"
    )
  save_gg(p, stub, 7.5, 6)
}

limma_two_group <- function(mat, group, near_level = "lung_tropic", far_level = "other_met") {
  group <- factor(as.character(group), levels = c(far_level, near_level))
  if (nlevels(droplevels(group)) < 2 || ncol(mat) < 4) {
    return(data.frame(gene = character(), log2FC = numeric(), AveExpr = numeric(),
                      pvalue = numeric(), padj_BH = numeric(), stringsAsFactors = FALSE))
  }
  design <- stats::model.matrix(~ group)
  fit <- limma::eBayes(limma::lmFit(mat, design), trend = TRUE, robust = TRUE)
  coefn <- grep(paste0("group", near_level), colnames(design), value = TRUE)
  tt <- limma::topTable(fit, coef = coefn, number = Inf, sort.by = "none")
  data.frame(
    gene = rownames(tt), log2FC = tt$logFC, AveExpr = tt$AveExpr,
    pvalue = tt$P.Value, padj_BH = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# 2. GPL96 注释（本地 → NCBI FTP 自动下载 → GEOquery → Bioconductor）
# -----------------------------------------------------------------------------
read_gpl96_annot_file <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  hdr_i <- which(grepl("^ID\\t", lines))[1]
  if (is.na(hdr_i)) return(NULL)
  tab <- utils::read.delim(textConnection(lines[hdr_i:length(lines)]),
                           check.names = FALSE, stringsAsFactors = FALSE, quote = "")
  if (!("ID" %in% names(tab))) return(NULL)
  sym_col <- grep("Gene [Ss]ymbol|gene_assignment|Gene symbol", names(tab), value = TRUE)[1]
  if (is.na(sym_col)) return(NULL)
  setNames(as.character(tab[[sym_col]]), as.character(tab[["ID"]]))
}

download_gpl96_annot <- function(dest) {
  urls <- c(
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz",
    "http://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz"
  )
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  for (u in urls) {
    log_msg("下载 GPL96.annot.gz: ", u)
    ok <- tryCatch({
      utils::download.file(u, destfile = dest, mode = "wb", quiet = TRUE)
      file.exists(dest) && file.info(dest)$size > 1000
    }, error = function(e) {
      log_msg("下载失败: ", e$message)
      FALSE
    })
    if (isTRUE(ok)) return(dest)
  }
  if (file.exists(dest) && isTRUE(file.info(dest)$size < 1000)) unlink(dest)
  NULL
}

load_gpl96 <- function() {
  local <- find_file(c("GPL96\\.annot(\\.gz)?$", "GPL96\\.soft(\\.gz)?$"))
  map <- NULL
  if (length(local) > 0) {
    log_msg("Reading GPL96 annotation: ", local[1])
    map <- tryCatch(read_gpl96_annot_file(local[1]), error = function(e) {
      log_msg("本地 GPL96 解析失败: ", e$message); NULL
    })
  }

  # 自动从 NCBI FTP 拉取（不依赖 Rtools / GEOquery）
  if (is.null(map)) {
    dest <- file.path(data_dir, "GPL96.annot.gz")
    got <- download_gpl96_annot(dest)
    if (!is.null(got)) {
      map <- tryCatch(read_gpl96_annot_file(got), error = function(e) {
        log_msg("下载的 GPL96 解析失败: ", e$message); NULL
      })
    }
  }

  if (is.null(map) && has_pkg("GEOquery")) {
    log_msg("尝试 GEOquery::getGEO(GPL96)")
    geo_dir <- file.path(data_dir, "GEO")
    dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)
    gpl <- tryCatch(GEOquery::getGEO("GPL96", destdir = geo_dir),
                    error = function(e) { log_msg(e$message); NULL })
    if (!is.null(gpl)) {
      td <- tryCatch(GEOquery::Table(gpl), error = function(e) NULL)
      if (!is.null(td)) {
        sym_col <- grep("Gene Symbol|gene_assignment", names(td), value = TRUE)[1]
        id_col <- if ("ID" %in% names(td)) "ID" else names(td)[1]
        if (!is.na(sym_col)) map <- setNames(as.character(td[[sym_col]]), as.character(td[[id_col]]))
      }
    }
  }

  if (is.null(map) && requireNamespace("hgu133a.db", quietly = TRUE)) {
    log_msg("使用 Bioconductor hgu133a.db 作探针注释")
    map <- tryCatch({
      suppressPackageStartupMessages(library(hgu133a.db))
      prb <- AnnotationDbi::keys(hgu133a.db::hgu133a.db, keytype = "PROBEID")
      ann <- AnnotationDbi::select(hgu133a.db::hgu133a.db, keys = prb,
                                  columns = "SYMBOL", keytype = "PROBEID")
      ann <- ann[!is.na(ann$SYMBOL) & nzchar(ann$SYMBOL), , drop = FALSE]
      # 一探针多符号时取第一个
      ann <- ann[!duplicated(ann$PROBEID), , drop = FALSE]
      setNames(as.character(ann$SYMBOL), as.character(ann$PROBEID))
    }, error = function(e) { log_msg(e$message); NULL })
  }

  if (is.null(map) || length(map) < 100) {
    stop(
      "无法获得 GPL96 探针→基因注释。\n",
      "请任选其一：\n",
      "  1) 浏览器下载后放到 E:/R/Nerve RNA/GPL96.annot.gz\n",
      "     https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz\n",
      "  2) 确认本机可访问 ftp.ncbi.nlm.nih.gov 后重新 source 本脚本（会自动下载）\n",
      "  3) install.packages 无法用时，也可 BiocManager::install('hgu133a.db')"
    )
  }
  log_msg("GPL96 注释探针数: ", length(map))
  map
}

# -----------------------------------------------------------------------------
# 3. 读 series matrix + 划分肺倾向 / 其他器官倾向
# -----------------------------------------------------------------------------
parse_series_matrix_geo <- function(path) {
  log_msg("Reading: ", path)
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  titles <- split_quoted(lines[grepl("^!Sample_title", lines)][1])
  geos <- {
    gl <- lines[grepl("^!Sample_geo_accession", lines)][1]
    if (is.na(gl)) paste0("S", seq_along(titles)) else split_quoted(gl)
  }
  char_lines <- lines[grepl("^!Sample_characteristics_ch1", lines)]
  char_mat <- lapply(char_lines, split_quoted)
  # 列对齐
  n <- length(titles)
  char_mat <- lapply(char_mat, function(v) { length(v) <- n; v })

  # 必须匹配行首特征名，避免 "met event" 误匹配 "binary 5y met event"
  get_char <- function(key_pat) {
    out <- rep(NA_character_, n)
    anchor <- paste0("(?i)^\\s*", key_pat, "\\s*:")
    for (row in char_mat) {
      hit <- grepl(anchor, row, perl = TRUE)
      if (!any(hit, na.rm = TRUE)) next
      if (mean(hit, na.rm = TRUE) > 0.3) {
        out <- sub(paste0("(?i)^\\s*", key_pat, "\\s*:\\s*"), "", row, perl = TRUE)
        out <- trimws(out)
        out[out %in% c("", "--", "NA", "null")] <- NA_character_
        return(out)
      }
    }
    out
  }

  begin <- which(grepl("series_matrix_table_begin", lines))[1]
  end <- which(grepl("series_matrix_table_end", lines))[1]
  tab <- utils::read.delim(textConnection(lines[(begin + 1):(end - 1)]),
                           check.names = FALSE, stringsAsFactors = FALSE, quote = "\"")
  probes <- as.character(tab[[1]])
  mat <- as.matrix(tab[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- probes
  if (ncol(mat) == n) colnames(mat) <- geos

  list(mat = mat, title = titles, gsm = geos,
       chars = list(
         lm_event = get_char("lm event"),
         bm_event = get_char("bm event"),
         met_event = get_char("met event"),
         lung_met_all = get_char("lung met all"),
         lung_met_first = get_char("lung met first event"),
         metastasis = get_char("metastasis"),
         lms_status = get_char("lms status"),
         tissue_type = get_char("tissue type")
       ),
       source = basename(path))
}

# 临床结局取值规范化：只接受 0/1；-- / NA / 空 = 未知（不可入组）
norm_event01 <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "--", "NA", "na", "NULL", "null")] <- NA_character_
  out <- rep(NA_character_, length(x))
  out[!is.na(x) & x %in% c("0", "1")] <- x[!is.na(x) & x %in% c("0", "1")]
  out
}

assign_tropism_gse2603 <- function(ds) {
  # 只要原发灶：标题 B###-T；排除细胞系 / 无随访结局
  is_primary <- grepl("^B[0-9]+-T$", ds$title)
  lm <- norm_event01(ds$chars$lm_event)
  met <- norm_event01(ds$chars$met_event)
  bm <- norm_event01(ds$chars$bm_event)

  trop <- rep(NA_character_, length(ds$title))
  correspond <- rep(NA_character_, length(ds$title))
  exclude_reason <- rep(NA_character_, length(ds$title))

  # 必须同一患者临床结局对应：lm event=1 → 该原位灶对应肺转移
  lung_ok <- is_primary & !is.na(lm) & lm == "1"
  # other：有转移、明确无肺转移（lm=0）；优先有骨转移记录，否则「非肺其他转移」
  other_ok <- is_primary & !is.na(met) & met == "1" & !is.na(lm) & lm == "0"
  no_met_ok <- is_primary & !is.na(met) & met == "0" & !is.na(lm) & lm == "0"

  trop[lung_ok] <- "lung_tropic"
  correspond[lung_ok] <- "same_patient_primary_with_documented_lung_met(lm_event=1)"
  trop[other_ok] <- "other_met"
  correspond[other_ok] <- ifelse(
    !is.na(bm[other_ok]) & bm[other_ok] == "1",
    "same_patient_primary_with_documented_bone_met_no_lung(bm=1,lm=0)",
    "same_patient_primary_with_nonlung_met(met=1,lm=0)"
  )
  trop[no_met_ok] <- "no_met"
  correspond[no_met_ok] <- "same_patient_primary_no_metastasis(not_used_in_Q1)"

  exclude_reason[!is_primary] <- "not_primary_or_cell_line"
  miss <- is_primary & is.na(trop)
  exclude_reason[miss] <- "incomplete_clinical_met_outcome(lm/met missing)"

  data.frame(
    sample = colnames(ds$mat), gsm = ds$gsm, title = ds$title,
    cohort = "GSE2603", is_primary = is_primary, tropism = trop,
    correspondence = correspond,
    exclude_reason = exclude_reason,
    lm_event = lm, bm_event = bm, met_event = met,
    lung_met_first = NA_character_,
    used_lms_signature = FALSE,
    stringsAsFactors = FALSE
  )
}

assign_tropism_gse5327 <- function(ds) {
  # 全部为 ER- 原发灶；禁止用 lms status 签名代替真实肺转移结局
  lung_all <- norm_event01(ds$chars$lung_met_all)
  lung_first <- norm_event01(ds$chars$lung_met_first)
  met <- norm_event01(ds$chars$metastasis)

  trop <- rep(NA_character_, length(ds$title))
  correspond <- rep(NA_character_, length(ds$title))
  exclude_reason <- rep(NA_character_, length(ds$title))

  # 同一患者随访：lung_met_all=1 表示该原位灶对应发生过肺转移
  lung_ok <- !is.na(lung_all) & lung_all == "1"
  other_ok <- !is.na(met) & met == "1" & !is.na(lung_all) & lung_all == "0"
  no_met_ok <- !is.na(met) & met == "0" & !is.na(lung_all) & lung_all == "0"

  trop[lung_ok] <- "lung_tropic"
  correspond[lung_ok] <- ifelse(
    !is.na(lung_first[lung_ok]) & lung_first[lung_ok] == "1",
    "same_patient_primary_with_lung_as_first_met(lung_met_first=1)",
    "same_patient_primary_with_documented_lung_met(lung_met_all=1)"
  )
  trop[other_ok] <- "other_met"
  correspond[other_ok] <- "same_patient_primary_with_nonlung_met(met=1,lung_met_all=0)"
  trop[no_met_ok] <- "no_met"
  correspond[no_met_ok] <- "same_patient_primary_no_metastasis(not_used_in_Q1)"

  miss <- is.na(trop)
  exclude_reason[miss] <- "incomplete_clinical_met_outcome(lung_met/metastasis missing)"

  data.frame(
    sample = colnames(ds$mat), gsm = ds$gsm, title = ds$title,
    cohort = "GSE5327", is_primary = TRUE, tropism = trop,
    correspondence = correspond,
    exclude_reason = exclude_reason,
    lm_event = lung_all, bm_event = NA_character_, met_event = met,
    lung_met_first = lung_first,
    used_lms_signature = FALSE,
    stringsAsFactors = FALSE
  )
}

probe2gene <- load_gpl96()

download_series_matrix <- function(acc, dest_file) {
  # GSE2603 -> GSE2nnn; GSE5327 -> GSE5nnn
  num <- as.integer(sub("^GSE", "", acc))
  bucket <- paste0("GSE", floor(num / 1000), "nnn")
  urls <- c(
    sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz",
            bucket, acc, acc),
    sprintf("http://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz",
            bucket, acc, acc)
  )
  dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
  for (u in urls) {
    log_msg("下载 ", acc, " series matrix: ", u)
    ok <- tryCatch({
      utils::download.file(u, destfile = dest_file, mode = "wb", quiet = TRUE)
      file.exists(dest_file) && file.info(dest_file)$size > 10000
    }, error = function(e) {
      log_msg("下载失败: ", e$message)
      FALSE
    })
    if (isTRUE(ok)) return(dest_file)
  }
  if (file.exists(dest_file) && isTRUE(file.info(dest_file)$size < 10000)) unlink(dest_file)
  NULL
}

ensure_series_matrix <- function(acc) {
  pat <- paste0(acc, "_series_matrix\\.txt(\\.gz)?$")
  hits <- find_file(pat)
  if (length(hits) > 0) return(hits[1])
  dest <- file.path(data_dir, paste0(acc, "_series_matrix.txt.gz"))
  got <- download_series_matrix(acc, dest)
  if (!is.null(got)) return(got)
  # GEOquery 备用：下载后在 GEO/ 与 data_dir 再找一次
  if (has_pkg("GEOquery")) {
    geo_dir <- file.path(data_dir, "GEO")
    dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)
    log_msg("GEOquery download ", acc)
    tryCatch(GEOquery::getGEO(acc, destdir = geo_dir, getGPL = FALSE),
             error = function(e) { log_msg(e$message); NULL })
    hits <- unique(c(
      list.files(geo_dir, pattern = pat, full.names = TRUE, ignore.case = TRUE, recursive = TRUE),
      find_file(pat)
    ))
    hits <- hits[file.exists(hits)]
    if (length(hits) > 0) {
      # 拷一份到数据根目录，方便下次直接找到
      if (!file.exists(dest)) {
        tryCatch(file.copy(hits[1], dest, overwrite = FALSE), error = function(e) NULL)
      }
      return(if (file.exists(dest)) dest else hits[1])
    }
  }
  NULL
}

load_dataset <- function(acc) {
  path <- ensure_series_matrix(acc)
  if (is.null(path)) {
    log_msg("未找到且无法下载 ", acc)
    return(NULL)
  }
  raw <- tryCatch(parse_series_matrix_geo(path), error = function(e) {
    log_msg(acc, " 解析失败: ", e$message); NULL
  })
  if (is.null(raw)) return(NULL)
  si <- if (identical(acc, "GSE2603")) assign_tropism_gse2603(raw) else assign_tropism_gse5327(raw)
  mat_g <- collapse_probes(raw$mat, probe2gene)
  log_msg(acc, " 载入成功: genes=", nrow(mat_g), " samples=", ncol(mat_g),
          " lung_tropic=", sum(si$tropism == "lung_tropic", na.rm = TRUE),
          " other_met=", sum(si$tropism == "other_met", na.rm = TRUE))
  list(mat = mat_g, si = si, source = raw$source)
}

datasets <- list()
for (acc in c("GSE2603", "GSE5327")) {
  ds <- load_dataset(acc)
  if (!is.null(ds)) datasets[[acc]] <- ds
}

if (length(datasets) == 0) {
  stop(
    "本地缺少 GSE2603/GSE5327 series matrix，且自动下载失败。\n",
    "请用浏览器下载后放到 E:/R/Nerve RNA/ ：\n",
    "  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE2nnn/GSE2603/matrix/GSE2603_series_matrix.txt.gz\n",
    "  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5327/matrix/GSE5327_series_matrix.txt.gz\n",
    "详见 Nerve_RNA_DOWNLOAD.txt"
  )
}

# 写出「原位灶 ↔ 转移器官」对应表（强制临床结局对应，禁止签名冒充）
si_all <- do.call(rbind, lapply(datasets, function(d) d$si))
utils::write.csv(si_all, file.path(result_dir, "00_SAMPLE_CORRESPONDENCE_primary_to_met_organ.csv"),
                 row.names = FALSE)
if (any(si_all$used_lms_signature %in% TRUE)) {
  stop("内部错误：检测到用 LMS 签名分组，已禁止。")
}
n_lung <- sum(si_all$tropism == "lung_tropic", na.rm = TRUE)
n_other <- sum(si_all$tropism == "other_met", na.rm = TRUE)
log_msg("对应入组合计: lung_tropic=", n_lung, " other_met=", n_other,
        " （均为同一患者临床随访结局；详见 00_SAMPLE_CORRESPONDENCE_*.csv）")
if (n_lung < 3) stop("有临床肺转移对应的原位灶不足（lung_tropic < 3）")
if (n_other < 3) stop("有临床非肺转移对应的原位灶不足（other_met < 3）")

# -----------------------------------------------------------------------------
# 4. 问题1：lung_tropic vs other_met —— 下调基因（仅 FC < 1）
# -----------------------------------------------------------------------------
emit_q1 <- function(tag, mat, si) {
  # 只保留「有临床对应」的 lung / other；排除无结局、签名预测、细胞系
  keep <- which(
    !is.na(si$tropism) &
      si$tropism %in% c("lung_tropic", "other_met") &
      !is.na(si$correspondence) &
      !grepl("not_used|signature|LMS", si$correspondence, ignore.case = TRUE) &
      !(si$used_lms_signature %in% TRUE)
  )
  if (sum(si$tropism[keep] == "lung_tropic", na.rm = TRUE) < 3 ||
      sum(si$tropism[keep] == "other_met", na.rm = TRUE) < 3) {
    log_msg("跳过问题1 ", tag, "：有临床对应的 lung/other 样本不足")
    return(NULL)
  }
  mat2 <- mat[, keep, drop = FALSE]
  si2 <- si[keep, , drop = FALSE]
  log_msg(tag, " Q1 入组对应: ",
          paste(unique(si2$correspondence[si2$tropism == "lung_tropic"]), collapse = " | "))
  keep_g <- rowSums(is.finite(mat2)) >= max(4, floor(0.7 * ncol(mat2)))
  mat2 <- mat2[keep_g, , drop = FALSE]
  mat2[!is.finite(mat2)] <- stats::median(mat2[is.finite(mat2)], na.rm = TRUE)
  if (max(mat2, na.rm = TRUE) > 100) mat2 <- log2(pmax(mat2, 1))
  mat_n <- quantile_norm(mat2)

  base <- file.path(result_dir, "Q1_lung_tropic_vs_other_met", tag)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(si2, file.path(base, "sample_groups.csv"), row.names = FALSE)

  de <- limma_two_group(mat_n, si2$tropism)
  de$fold_change <- 2^de$log2FC
  de <- de[order(de$pvalue, de$log2FC), ]
  utils::write.csv(de, file.path(base, "DE_full_lung_vs_other.csv"), row.names = FALSE)

  # 仅：p < 0.01 且 FC < 1（下调）；不做其他 FC 分层
  down <- de[!is.na(de$pvalue) & de$pvalue < p_cutoff &
               !is.na(de$fold_change) & de$fold_change < fc_down_max, , drop = FALSE]
  down <- down[order(down$fold_change, decreasing = FALSE), , drop = FALSE]
  utils::write.csv(down, file.path(base, "DOWN_genes_p01_FC_lt1.csv"), row.names = FALSE)
  volcano_down(de, paste0(tag, " | lung-tropic vs other-met | DOWN FC<1"),
               file.path(base, "volcano_DOWN_FC_lt1"))

  sheets <- list(
    DOWN_FC_lt1 = down,
    DE_full = de,
    sample_groups = si2
  )
  write_xlsx_safe(sheets, file.path(base, paste0(tag, "_Q1_DOWN_FC_lt1.xlsx")))
  utils::write.csv(down, file.path(result_dir, paste0("01_", tag, "_DOWN_genes_p01_FC_lt1.csv")),
                   row.names = FALSE)
  write_xlsx_safe(list(DOWN_FC_lt1 = down, DE_full = de),
                  file.path(result_dir, paste0("01_", tag, "_DOWN_FC_lt1.xlsx")))

  log_msg(tag, " Q1: lung_tropic n=", sum(si2$tropism == "lung_tropic"),
          " other_met n=", sum(si2$tropism == "other_met"),
          " DOWN(p<0.01, FC<1)=", nrow(down))
  list(de = de, down = down, mat_n = mat_n, si = si2, base = base)
}

q1_res <- list()
for (nm in names(datasets)) {
  q1_res[[nm]] <- emit_q1(nm, datasets[[nm]]$mat, datasets[[nm]]$si)
}

# 合并两队列（同平台），cohort 作协变量
if (length(datasets) >= 2 &&
    !is.null(q1_res[["GSE2603"]]) && !is.null(q1_res[["GSE5327"]])) {
  common <- Reduce(intersect, lapply(datasets, function(d) rownames(d$mat)))
  mats <- list(); sis <- list()
  for (nm in names(datasets)) {
    si <- datasets[[nm]]$si
    keep <- which(
      !is.na(si$tropism) & si$tropism %in% c("lung_tropic", "other_met") &
        !is.na(si$correspondence) & !(si$used_lms_signature %in% TRUE)
    )
    m <- datasets[[nm]]$mat[common, keep, drop = FALSE]
    s <- si[keep, , drop = FALSE]
    if (max(m, na.rm = TRUE) > 100) m <- log2(pmax(m, 1))
    mats[[nm]] <- m
    sis[[nm]] <- s
  }
  mat_c <- do.call(cbind, mats)
  si_c <- do.call(rbind, sis)
  rownames(si_c) <- colnames(mat_c)
  mat_c <- quantile_norm(mat_c)
  trop <- factor(si_c$tropism, levels = c("other_met", "lung_tropic"))
  cohort <- factor(si_c$cohort)
  design <- stats::model.matrix(~ trop + cohort)
  fit <- limma::eBayes(limma::lmFit(mat_c, design), trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit, coef = "troplung_tropic", number = Inf, sort.by = "none")
  de <- data.frame(gene = rownames(tt), log2FC = tt$logFC, AveExpr = tt$AveExpr,
                   pvalue = tt$P.Value, padj_BH = tt$adj.P.Val,
                   fold_change = 2^tt$logFC, stringsAsFactors = FALSE)
  de <- de[order(de$pvalue, de$log2FC), ]
  down <- de[!is.na(de$pvalue) & de$pvalue < p_cutoff &
               !is.na(de$fold_change) & de$fold_change < fc_down_max, , drop = FALSE]
  down <- down[order(down$fold_change), , drop = FALSE]

  base <- file.path(result_dir, "Q1_lung_tropic_vs_other_met", "GSE2603_plus_GSE5327")
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(si_c, file.path(base, "sample_groups.csv"), row.names = FALSE)
  utils::write.csv(de, file.path(base, "DE_full_lung_vs_other.csv"), row.names = FALSE)
  utils::write.csv(down, file.path(base, "DOWN_genes_p01_FC_lt1.csv"), row.names = FALSE)
  volcano_down(de, "Pooled | lung-tropic vs other-met | DOWN FC<1",
               file.path(base, "volcano_DOWN_FC_lt1"))
  write_xlsx_safe(list(DOWN_FC_lt1 = down, DE_full = de, sample_groups = si_c),
                  file.path(base, "Pooled_Q1_DOWN_FC_lt1.xlsx"))
  utils::write.csv(down, file.path(result_dir, "01_POOLED_DOWN_genes_p01_FC_lt1.csv"),
                   row.names = FALSE)
  write_xlsx_safe(list(DOWN_FC_lt1 = down, DE_full = de),
                  file.path(result_dir, "01_POOLED_DOWN_FC_lt1.xlsx"))
  q1_res[["POOLED"]] <- list(de = de, down = down, mat_n = mat_c, si = si_c, base = base)
  log_msg("POOLED Q1: lung=", sum(si_c$tropism == "lung_tropic"),
          " other=", sum(si_c$tropism == "other_met"),
          " DOWN(p<0.01, FC<1)=", nrow(down))
}

# -----------------------------------------------------------------------------
# 5. 问题2：肺倾向原位组织中，与神经浸润负相关的「蛋白」
# -----------------------------------------------------------------------------
# 5a. 可选真实蛋白组
load_user_proteomics <- function() {
  mats <- find_file(c("proteomics_pg_matrix\\.csv$", "report\\.pg_matrix\\.tsv$",
                      "proteomics.*\\.(csv|tsv)$"))
  metas <- find_file(c("proteomics_sample_meta\\.csv$", "proteomics.*meta.*\\.(csv|tsv)$"))
  if (length(mats) == 0) return(NULL)
  path <- mats[1]
  log_msg("读取用户蛋白组: ", path)
  sep <- if (grepl("\\.tsv$", path, ignore.case = TRUE)) "\t" else ","
  tab <- utils::read.delim(path, sep = sep, check.names = FALSE, stringsAsFactors = FALSE)
  # 第一列基因/蛋白符号
  id_col <- names(tab)[1]
  mat <- as.matrix(tab[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- make.unique(as.character(tab[[id_col]]))
  si <- data.frame(sample = colnames(mat), tropism = "lung_tropic",
                   stringsAsFactors = FALSE)
  if (length(metas) > 0) {
    meta <- utils::read.delim(metas[1], sep = if (grepl("\\.tsv$", metas[1], ignore.case = TRUE)) "\t" else ",",
                              check.names = FALSE, stringsAsFactors = FALSE)
    sn <- intersect(c("sample", "Sample", "sample_id"), names(meta))[1]
    tn <- intersect(c("tropism", "group", "organ_tropism"), names(meta))[1]
    if (!is.na(sn)) {
      m <- match(si$sample, meta[[sn]])
      if (!is.na(tn)) si$tropism[!is.na(m)] <- as.character(meta[[tn]][m[!is.na(m)]])
    }
  }
  list(mat = mat, si = si, source = basename(path), is_protein = TRUE)
}

cor_neg_with_nerve <- function(mat, score, label, is_protein = FALSE) {
  ok <- is.finite(score)
  if (sum(ok) < 5) {
    log_msg("跳过问题2 ", label, "：肺倾向样本不足")
    return(NULL)
  }
  mat2 <- mat[, ok, drop = FALSE]
  sc <- score[ok]
  rho <- apply(mat2, 1, function(v) {
    if (sum(is.finite(v)) < 5) return(NA_real_)
    suppressWarnings(stats::cor(v, sc, method = "spearman", use = "pairwise.complete.obs"))
  })
  # Spearman 近似 p（大样本）；小样本用 cor.test
  pvals <- apply(mat2, 1, function(v) {
    if (sum(is.finite(v) & is.finite(sc)) < 5) return(NA_real_)
    tryCatch(stats::cor.test(v, sc, method = "spearman", exact = FALSE)$p.value,
             error = function(e) NA_real_)
  })
  out <- data.frame(
    feature = rownames(mat2),
    molecule_type = if (is_protein) "protein" else "RNA_proxy_for_protein",
    spearman_vs_nerve_score = as.numeric(rho),
    pvalue = as.numeric(pvals),
    mean_abundance = rowMeans(mat2, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  out$padj_BH <- p.adjust(out$pvalue, method = "BH")
  out$negatively_correlated <- !is.na(out$spearman_vs_nerve_score) &
    out$spearman_vs_nerve_score < 0 &
    !is.na(out$pvalue) & out$pvalue < p_cutoff
  out <- out[order(out$spearman_vs_nerve_score, out$pvalue), ]
  out
}

emit_q2 <- function(tag, mat, si, is_protein = FALSE, lung_sample_ids = NULL) {
  # 问题2必须用「有肺转移临床对应」的原位灶；若传入 Q1 的 lung 样本 ID，则严格取交集
  lung_idx <- which(
    !is.na(si$tropism) & si$tropism == "lung_tropic" &
      !is.na(si$correspondence) &
      grepl("documented_lung|lung_as_first", si$correspondence, ignore.case = TRUE) &
      !(si$used_lms_signature %in% TRUE)
  )
  if (!is.null(lung_sample_ids)) {
    lung_idx <- lung_idx[si$sample[lung_idx] %in% lung_sample_ids]
  }
  if (length(lung_idx) < 5) {
    log_msg("跳过问题2 ", tag, "：有肺转移临床对应的原位灶 n < 5")
    return(NULL)
  }
  mat_l <- mat[, lung_idx, drop = FALSE]
  si_l <- si[lung_idx, , drop = FALSE]
  log_msg(tag, " Q2 仅用临床对应肺转移原位灶 n=", nrow(si_l))
  if (max(mat_l, na.rm = TRUE) > 100 && !is_protein) mat_l <- log2(pmax(mat_l, 1))
  mat_l <- quantile_norm(mat_l)

  lig <- intersect(nerve_ligand_genes, rownames(mat_l))
  score <- zmean_score(mat_l, lig)
  si_l$nerve_invasion_score <- score
  log_msg(tag, " Q2 nerve signature genes used: ", paste(lig, collapse = ", "))

  tab <- cor_neg_with_nerve(mat_l, score, tag, is_protein = is_protein)
  if (is.null(tab)) return(NULL)
  base <- file.path(result_dir, "Q2_neg_cor_nerve_in_lung_tropic", tag)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(si_l, file.path(base, "lung_tropic_samples_with_nerve_score.csv"),
                   row.names = FALSE)
  utils::write.csv(tab, file.path(base, "ALL_features_vs_nerve_score.csv"), row.names = FALSE)
  neg <- tab[tab$negatively_correlated %in% TRUE, , drop = FALSE]
  utils::write.csv(neg, file.path(base, "NEG_cor_p01_vs_nerve_score.csv"), row.names = FALSE)

  # 图：分数分布 + top 负相关散点
  if (has_pkg("ggplot2")) {
    p1 <- ggplot2::ggplot(si_l, ggplot2::aes(x = "lung_tropic", y = nerve_invasion_score)) +
      ggplot2::geom_boxplot(fill = "#A8DADC", width = 0.4) +
      ggplot2::geom_jitter(width = 0.08, size = 1.5, alpha = 0.7) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(title = paste(tag, "nerve invasion score in lung-tropic primaries"),
                    x = NULL, y = "nerve_invasion_score")
    save_gg(p1, file.path(base, "nerve_score_boxplot"), 4.5, 5)

    topn <- utils::head(neg, 12)
    if (nrow(topn) > 0) {
      plot_df <- data.frame(
        feature = factor(rep(topn$feature, each = ncol(mat_l)), levels = topn$feature),
        abundance = as.numeric(t(mat_l[topn$feature, , drop = FALSE])),
        nerve = rep(score, times = nrow(topn))
      )
      p2 <- ggplot2::ggplot(plot_df, ggplot2::aes(abundance, nerve)) +
        ggplot2::geom_point(alpha = 0.7, size = 1.2, color = "#1D3557") +
        ggplot2::geom_smooth(method = "lm", se = FALSE, color = "#E63946", linewidth = 0.6) +
        ggplot2::facet_wrap(~ feature, scales = "free_x") +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::labs(title = paste(tag, "Top negative correlations vs nerve score"),
                      x = if (is_protein) "protein abundance" else "RNA abundance (protein proxy)",
                      y = "nerve_invasion_score")
      save_gg(p2, file.path(base, "top_neg_cor_scatter"), 10, 8)
    }
  }

  sheets <- list(
    NEG_cor_p01 = neg,
    ALL_correlations = tab,
    lung_samples = si_l,
    nerve_genes_used = data.frame(gene = lig, stringsAsFactors = FALSE)
  )
  write_xlsx_safe(sheets, file.path(base, paste0(tag, "_Q2_NEG_cor_nerve.xlsx")))
  utils::write.csv(neg, file.path(result_dir, paste0("02_", tag, "_NEG_cor_nerve_p01.csv")),
                   row.names = FALSE)
  write_xlsx_safe(list(NEG_cor_p01 = neg, ALL = tab),
                  file.path(result_dir, paste0("02_", tag, "_NEG_cor_nerve.xlsx")))
  log_msg(tag, " Q2: lung_tropic n=", ncol(mat_l),
          " neg-cor(p<0.01)=", nrow(neg),
          if (is_protein) " [protein]" else " [RNA proxy for protein]")
  invisible(list(tab = tab, neg = neg, si = si_l, lig = lig))
}

q2_res <- list()
for (nm in names(datasets)) {
  lung_ids <- NULL
  if (!is.null(q1_res[[nm]]) && !is.null(q1_res[[nm]]$si)) {
    lung_ids <- q1_res[[nm]]$si$sample[q1_res[[nm]]$si$tropism == "lung_tropic"]
  }
  q2_res[[nm]] <- emit_q2(nm, datasets[[nm]]$mat, datasets[[nm]]$si,
                          is_protein = FALSE, lung_sample_ids = lung_ids)
}
# 合并肺倾向样本做问题2：必须与 POOLED Q1 的 lung 样本严格对应
if (!is.null(q1_res[["POOLED"]])) {
  si_p <- q1_res[["POOLED"]]$si
  mat_p <- q1_res[["POOLED"]]$mat_n
  lung_ids <- si_p$sample[si_p$tropism == "lung_tropic"]
  q2_res[["POOLED"]] <- emit_q2("POOLED", mat_p, si_p, is_protein = FALSE,
                                lung_sample_ids = lung_ids)
}

prot <- load_user_proteomics()
if (!is.null(prot)) {
  # 仅肺倾向（或用户 meta 指定）
  q2_res[["PROTEOMICS"]] <- emit_q2("USER_PROTEOMICS", prot$mat, prot$si, is_protein = TRUE)
}

# -----------------------------------------------------------------------------
# 6. 方案说明
# -----------------------------------------------------------------------------
n_summary <- lapply(datasets, function(d) {
  paste0(d$si$cohort[1], ": lung_tropic=",
         sum(d$si$tropism == "lung_tropic", na.rm = TRUE),
         " other_met=", sum(d$si$tropism == "other_met", na.rm = TRUE),
         " no_met=", sum(d$si$tropism == "no_met", na.rm = TRUE))
})

protocol <- c(
  "============================================================",
  "E:/R/Nerve RNA | 原位乳腺癌：肺转移倾向 vs 其他器官倾向",
  "============================================================",
  "",
  "【问题1】倾向肺转移的原位组织 vs 倾向其他器官转移的原位组织 —— 下调基因",
  "  【对应硬性要求】原位灶必须与转移器官结局一一对应（同一患者临床随访）：",
  "    lung_tropic = 该患者原发灶 + 随访确实发生肺转移 (GSE2603 lm event=1; GSE5327 lung met all=1)",
  "    other_met   = 该患者原发灶 + 发生转移但未发生肺转移",
  "    禁止：LMS 签名预测、无结局样本、细胞系、转移灶活检冒充原位",
  "    样本对应表: 00_SAMPLE_CORRESPONDENCE_primary_to_met_organ.csv",
  "  规则: 先 p < 0.01，再只取 FC < 1（下调）；不做 1.25/1.5/2 分层",
  "  输出:",
  "    results_primary_lung_tropism/Q1_lung_tropic_vs_other_met/<队列>/",
  "      DOWN_genes_p01_FC_lt1.csv",
  "      volcano_DOWN_FC_lt1.png / .pdf",
  "      <队列>_Q1_DOWN_FC_lt1.xlsx",
  "    总表: 01_POOLED_DOWN_FC_lt1.xlsx / 01_*_DOWN_genes_p01_FC_lt1.csv",
  "",
  "【问题2】在倾向肺转移的原位组织中，与神经浸润负相关的蛋白",
  "  仅使用问题1中有肺转移临床对应的同一批 lung_tropic 原位灶样本",
  "  神经浸润分数: 神经亲和配体/侵袭签名基因的 z-mean",
  "  在 lung_tropic 样本内做 Spearman 相关；取 rho<0 且 p<0.01",
  "  默认无公开匹配蛋白组 → 用转录本丰度作蛋白代理 (molecule_type=RNA_proxy_for_protein)",
  "  若自备蛋白组: 放入 proteomics_pg_matrix.csv (+ proteomics_sample_meta.csv)",
  "  输出:",
  "    results_primary_lung_tropism/Q2_neg_cor_nerve_in_lung_tropic/<队列>/",
  "    总表: 02_*_NEG_cor_nerve.xlsx",
  "",
  "样本数: ", paste(unlist(n_summary), collapse = " | "),
  "",
  "数据文件:",
  "  GSE2603_series_matrix.txt.gz  (必须优先)",
  "  GSE5327_series_matrix.txt.gz  (验证/合并)",
  "  GPL96.annot.gz                (探针注释，推荐)",
  "",
  "运行:",
  "  setwd(\"E:/R/Nerve RNA\")",
  "  source(\"Nerve_RNA_primary_lung_tropism.R\")",
  "",
  "注意: 关联 ≠ 因果；脑实质转移 ≠ 外周神经浸润(PNI)。",
  "旧脚本 Nerve_RNA_nerve_infiltration.R 针对转移灶 GSE175692，与本问题不同。"
)
writeLines(protocol, file.path(result_dir, "00_PROTOCOL_原位肺倾向两个问题.txt"))
log_msg("Done. Results -> ", result_dir)
log_msg("Read -> ", file.path(result_dir, "00_PROTOCOL_原位肺倾向两个问题.txt"))
