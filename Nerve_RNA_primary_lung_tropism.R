#!/usr/bin/env Rscript
# =============================================================================
# E:\R\Nerve RNA —— 原位（原发）乳腺癌：肺转移倾向 vs 其他器官转移倾向
#
# 问题1：倾向肺转移的原位乳腺癌组织 vs 倾向其他器官转移的原位乳腺癌组织
#         → 低表达基因有哪些？火山图 + Excel
# 问题2：在「倾向肺转移」的原位乳腺癌组织中，与神经浸润负相关的蛋白有哪些？
#         （公开队列几乎无匹配蛋白组；默认用转录本丰度作蛋白代理，并支持自备蛋白矩阵）
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
fc_cutoffs <- c("FC_1.25" = 1.25, "FC_1.5" = 1.5, "FC_2" = 2)

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

volcano_down <- function(de, title, stub, fc_line = 1.25) {
  if (nrow(de) == 0) return(invisible(NULL))
  df <- de
  df$y <- -log10(pmax(df$pvalue, 1e-300))
  df$col <- "ns"
  df$col[!is.na(df$pvalue) & df$pvalue < p_cutoff & df$log2FC <= -log2(fc_line)] <- "down"
  df$col[!is.na(df$pvalue) & df$pvalue < p_cutoff & df$log2FC >= log2(fc_line)] <- "up"
  p <- ggplot2::ggplot(df, ggplot2::aes(log2FC, y, color = col)) +
    ggplot2::geom_point(alpha = 0.55, size = 1.1) +
    ggplot2::scale_color_manual(values = c(ns = "grey70", down = "#1D3557", up = "#D62828")) +
    ggplot2::geom_vline(xintercept = c(-log2(fc_line), log2(fc_line)), linetype = 2, color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = 2, color = "grey40") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, subtitle = paste0("Highlight DOWN: p<", p_cutoff,
                                                   " & FC>=", fc_line, " (lung-tropic vs other-met)"),
                  y = "-log10(p)", x = "log2FC (lung_tropic / other_organ_met)")
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

assign_tropism_gse2603 <- function(ds) {
  # 只要原发灶：标题 B###-T；排除细胞系
  is_primary <- grepl("^B[0-9]+-T$", ds$title)
  lm <- ds$chars$lm_event
  met <- ds$chars$met_event
  trop <- rep(NA_character_, length(ds$title))
  trop[is_primary & lm == "1"] <- "lung_tropic"
  trop[is_primary & met == "1" & lm == "0"] <- "other_met"
  trop[is_primary & met == "0"] <- "no_met"
  data.frame(
    sample = colnames(ds$mat), gsm = ds$gsm, title = ds$title,
    cohort = "GSE2603", is_primary = is_primary, tropism = trop,
    lm_event = lm, bm_event = ds$chars$bm_event, met_event = met,
    stringsAsFactors = FALSE
  )
}

assign_tropism_gse5327 <- function(ds) {
  # 全部为 ER- 原发灶
  lung <- ds$chars$lung_met_all
  met <- ds$chars$metastasis
  trop <- rep(NA_character_, length(ds$title))
  trop[lung == "1"] <- "lung_tropic"
  trop[met == "1" & lung == "0"] <- "other_met"
  trop[met == "0"] <- "no_met"
  data.frame(
    sample = colnames(ds$mat), gsm = ds$gsm, title = ds$title,
    cohort = "GSE5327", is_primary = TRUE, tropism = trop,
    lm_event = lung, bm_event = NA_character_, met_event = met,
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

# -----------------------------------------------------------------------------
# 4. 问题1：lung_tropic vs other_met —— 低表达基因
# -----------------------------------------------------------------------------
emit_q1 <- function(tag, mat, si) {
  keep <- which(!is.na(si$tropism) & si$tropism %in% c("lung_tropic", "other_met"))
  if (sum(si$tropism == "lung_tropic", na.rm = TRUE) < 3 ||
      sum(si$tropism == "other_met", na.rm = TRUE) < 3) {
    log_msg("跳过问题1 ", tag, "：lung/other 样本不足")
    return(NULL)
  }
  mat2 <- mat[, keep, drop = FALSE]
  si2 <- si[keep, , drop = FALSE]
  # 过滤 + 标准化
  keep_g <- rowSums(is.finite(mat2)) >= max(4, floor(0.7 * ncol(mat2)))
  mat2 <- mat2[keep_g, , drop = FALSE]
  mat2[!is.finite(mat2)] <- stats::median(mat2[is.finite(mat2)], na.rm = TRUE)
  # Affy series matrix 多为已 log2；若值域大则 log2
  if (max(mat2, na.rm = TRUE) > 100) mat2 <- log2(pmax(mat2, 1))
  mat_n <- quantile_norm(mat2)

  base <- file.path(result_dir, "Q1_lung_tropic_vs_other_met", tag)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(si2, file.path(base, "sample_groups.csv"), row.names = FALSE)

  de <- limma_two_group(mat_n, si2$tropism)
  de$fold_change <- 2^de$log2FC
  de$down_FC <- 2^(-de$log2FC)  # 下调倍数（lung 相对 other）
  de <- de[order(de$pvalue, de$log2FC), ]
  utils::write.csv(de, file.path(base, "DE_full_lung_vs_other.csv"), row.names = FALSE)

  sheets <- list(DE_full = de, sample_groups = si2)
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    od <- file.path(base, "FoldChange", nm)
    dir.create(od, recursive = TRUE, showWarnings = FALSE)
    down <- de[!is.na(de$pvalue) & de$pvalue < p_cutoff &
                 de$log2FC < 0 & de$down_FC >= fc, , drop = FALSE]
    down <- down[order(down$log2FC, decreasing = FALSE), , drop = FALSE]
    utils::write.csv(down, file.path(od, paste0(nm, "_DOWN_genes_lung_vs_other.csv")),
                     row.names = FALSE)
    volcano_down(de, paste0(tag, " | ", nm, " | lung-tropic vs other-met PRIMARY"),
                 file.path(od, paste0(nm, "_volcano")), fc)
    sheets[[paste0(nm, "_DOWN")]] <- down
  }
  xlsx_path <- file.path(base, paste0(tag, "_Q1_lung_vs_other_DOWN.xlsx"))
  write_xlsx_safe(sheets, xlsx_path)
  # 总表副本
  utils::write.csv(sheets[["FC_1.25_DOWN"]],
                   file.path(result_dir, paste0("01_", tag, "_DOWN_genes_p01_FC1.25.csv")),
                   row.names = FALSE)

  log_msg(tag, " Q1: lung_tropic n=", sum(si2$tropism == "lung_tropic"),
          " other_met n=", sum(si2$tropism == "other_met"),
          " DOWN(p<0.01,FC>=1.25)=", nrow(sheets[["FC_1.25_DOWN"]]))
  list(de = de, mat_n = mat_n, si = si2, sheets = sheets, base = base)
}

q1_res <- list()
for (nm in names(datasets)) {
  q1_res[[nm]] <- emit_q1(nm, datasets[[nm]]$mat, datasets[[nm]]$si)
}

# 合并两队列（同平台），cohort 作协变量
if (length(datasets) >= 2 &&
    !is.null(q1_res[["GSE2603"]]) && !is.null(q1_res[["GSE5327"]])) {
  # 用原始 gene 矩阵，重新取 tropism 样本
  common <- Reduce(intersect, lapply(datasets, function(d) rownames(d$mat)))
  mats <- list(); sis <- list()
  for (nm in names(datasets)) {
    si <- datasets[[nm]]$si
    keep <- si$tropism %in% c("lung_tropic", "other_met")
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
  # limma ~ tropism + cohort
  trop <- factor(si_c$tropism, levels = c("other_met", "lung_tropic"))
  cohort <- factor(si_c$cohort)
  design <- stats::model.matrix(~ trop + cohort)
  fit <- limma::eBayes(limma::lmFit(mat_c, design), trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit, coef = "troplung_tropic", number = Inf, sort.by = "none")
  de <- data.frame(gene = rownames(tt), log2FC = tt$logFC, AveExpr = tt$AveExpr,
                   pvalue = tt$P.Value, padj_BH = tt$adj.P.Val,
                   fold_change = 2^tt$logFC, down_FC = 2^(-tt$logFC),
                   stringsAsFactors = FALSE)
  de <- de[order(de$pvalue, de$log2FC), ]
  base <- file.path(result_dir, "Q1_lung_tropic_vs_other_met", "GSE2603_plus_GSE5327")
  dir.create(file.path(base, "FoldChange"), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(si_c, file.path(base, "sample_groups.csv"), row.names = FALSE)
  utils::write.csv(de, file.path(base, "DE_full_lung_vs_other.csv"), row.names = FALSE)
  sheets <- list(DE_full = de, sample_groups = si_c)
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    od <- file.path(base, "FoldChange", nm)
    dir.create(od, recursive = TRUE, showWarnings = FALSE)
    down <- de[!is.na(de$pvalue) & de$pvalue < p_cutoff &
                 de$log2FC < 0 & de$down_FC >= fc, , drop = FALSE]
    down <- down[order(down$log2FC), , drop = FALSE]
    utils::write.csv(down, file.path(od, paste0(nm, "_DOWN_genes_lung_vs_other.csv")),
                     row.names = FALSE)
    volcano_down(de, paste0("Pooled | ", nm, " | lung-tropic vs other-met PRIMARY"),
                 file.path(od, paste0(nm, "_volcano")), fc)
    sheets[[paste0(nm, "_DOWN")]] <- down
  }
  write_xlsx_safe(sheets, file.path(base, "Pooled_Q1_lung_vs_other_DOWN.xlsx"))
  utils::write.csv(sheets[["FC_1.25_DOWN"]],
                   file.path(result_dir, "01_POOLED_DOWN_genes_p01_FC1.25.csv"),
                   row.names = FALSE)
  write_xlsx_safe(list(FC_1.25_DOWN = sheets[["FC_1.25_DOWN"]], DE_full = de),
                  file.path(result_dir, "01_POOLED_lung_vs_other_DOWN.xlsx"))
  q1_res[["POOLED"]] <- list(de = de, mat_n = mat_c, si = si_c, sheets = sheets, base = base)
  log_msg("POOLED Q1: lung=", sum(si_c$tropism == "lung_tropic"),
          " other=", sum(si_c$tropism == "other_met"),
          " DOWN(FC1.25)=", nrow(sheets[["FC_1.25_DOWN"]]))
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

emit_q2 <- function(tag, mat, si, is_protein = FALSE) {
  # 禁止用带 NA 的逻辑下标（会掺进 NA 行，虚增样本）
  lung_idx <- which(!is.na(si$tropism) & si$tropism == "lung_tropic")
  if (length(lung_idx) < 5) {
    log_msg("跳过问题2 ", tag, "：lung_tropic n < 5")
    return(NULL)
  }
  mat_l <- mat[, lung_idx, drop = FALSE]
  si_l <- si[lung_idx, , drop = FALSE]
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
  q2_res[[nm]] <- emit_q2(nm, datasets[[nm]]$mat, datasets[[nm]]$si, is_protein = FALSE)
}
# 合并肺倾向样本做问题2
if (!is.null(q1_res[["POOLED"]])) {
  si_p <- q1_res[["POOLED"]]$si
  mat_p <- q1_res[["POOLED"]]$mat_n
  q2_res[["POOLED"]] <- emit_q2("POOLED", mat_p, si_p, is_protein = FALSE)
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
  "【问题1】倾向肺转移的原位组织 vs 倾向其他器官转移的原位组织 —— 低表达基因",
  "  定义:",
  "    lung_tropic = 原发灶且发生肺转移 (GSE2603: lm event=1; GSE5327: lung met all=1)",
  "    other_met   = 原发灶发生转移但未发生肺转移 (多为骨等其他部位)",
  "  规则: 先 p < 0.01，再下调 FC >= 1.25 / 1.5 / 2（lung 相对 other）",
  "  输出:",
  "    results_primary_lung_tropism/Q1_lung_tropic_vs_other_met/<队列>/",
  "      FoldChange/FC_*/FC_*_DOWN_genes_lung_vs_other.csv",
  "      FoldChange/FC_*/FC_*_volcano.png",
  "      <队列>_Q1_lung_vs_other_DOWN.xlsx",
  "    总表: 01_POOLED_lung_vs_other_DOWN.xlsx / 01_*_DOWN_genes_p01_FC1.25.csv",
  "",
  "【问题2】在倾向肺转移的原位组织中，与神经浸润负相关的蛋白",
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
