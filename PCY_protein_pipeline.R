#!/usr/bin/env Rscript
# =============================================================================
# PCY IP 蛋白质谱：PCY vs EV
# 对照 EV1-3，实验 PCY1-3。每个重复单独保留，不把重复合并成一个样品。
#
# 1. limma 比较 PCY vs EV，输出差异表和火山图
# 2. 上调且 p < 0.05、FC > 1 / 1.25 / 1.5 的蛋白，分别做 GO、KEGG、通路富集，画气泡图
#
# 运行（把定量表放在数据目录后）：
#   setwd("E:/R/PCY_protein")
#   source("PCY_protein_pipeline.R")
# 也可指定目录：Sys.setenv(PCY_PROTEIN_DIR = "E:/R/PCY_protein")
#
# 接受的输入（二选一）：
#   1. 六个样品各一个文件，文件名就是样品名：EV1.txt、EV2.txt、EV3.txt、PCY1.txt、PCY2.txt、PCY3.txt
#   2. 一张宽表（xlsx / csv / tsv / MaxQuant proteinGroups.txt），列名里能识别这六个样品
# 定量列优先使用 LFQ，其次 iBAQ，再次 Intensity / Abundance / 强度。
# 物种默认人类（org.Hs.eg.db，KEGG hsa）。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")
options(clusterProfiler.download.method = "auto")

# -----------------------------------------------------------------------------
# 0. 依赖包
# -----------------------------------------------------------------------------
cran_required <- c("ggplot2", "ggrepel")
cran_optional <- c("readxl", "writexl", "msigdbr")
bioc_required <- c("limma", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi")
bioc_optional <- c("ReactomePA")

install_if_missing <- function(pkgs, bioc = FALSE, required = TRUE) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) == 0) return(invisible(TRUE))
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    tryCatch(
      BiocManager::install(miss, update = FALSE, ask = FALSE),
      error = function(e) message("Bioconductor install failed: ", e$message)
    )
  } else {
    tryCatch(
      install.packages(miss, repos = "https://cloud.r-project.org"),
      error = function(e) message("CRAN install failed: ", e$message)
    )
  }
  still <- miss[!vapply(miss, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still) > 0 && required) {
    stop("缺少必需 R 包: ", paste(still, collapse = ", "))
  }
  if (length(still) > 0) message("可选包未安装，相关分析将跳过: ", paste(still, collapse = ", "))
  invisible(TRUE)
}

install_if_missing(cran_required, bioc = FALSE, required = TRUE)
install_if_missing(cran_optional, bioc = FALSE, required = FALSE)
install_if_missing(bioc_required, bioc = TRUE, required = TRUE)
install_if_missing(bioc_optional, bioc = TRUE, required = FALSE)

safe_library <- function(pkgs) {
  for (p in pkgs) {
    if (requireNamespace(p, quietly = TRUE)) {
      suppressPackageStartupMessages(library(p, character.only = TRUE))
    }
  }
}
safe_library(c(cran_required, cran_optional, bioc_required, bioc_optional))

has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

# -----------------------------------------------------------------------------
# 1. 路径与参数
# -----------------------------------------------------------------------------
# 差异蛋白：未校正 p < pvalue_cutoff。n = 3 时 BH 往往过严，padj 仍写入结果表。
# 富集名单在此之上再要求上调 FC 大于对应档位。只分析 PCY 相对 EV 的上调。
pvalue_cutoff <- 0.05
fc_cutoffs <- c("FC_1" = 1, "FC_1.25" = 1.25, "FC_1.5" = 1.5)
wiki_unavailable <- FALSE
min_valid_in_group <- 2L
orgdb_name <- "org.Hs.eg.db"
kegg_organism <- "hsa"
wp_organism <- "Homo sapiens"
reactome_organism <- "human"

resolve_project_dir <- function() {
  env_dir <- Sys.getenv("PCY_PROTEIN_DIR", unset = "")
  candidates <- c(env_dir, "E:/R/PCY_protein", "E:\\R\\PCY_protein", getwd())
  candidates <- unique(candidates[nzchar(candidates)])
  existing <- candidates[dir.exists(candidates)]
  if (length(existing) == 0) {
    stop(
      "找不到数据目录。请把质谱定量表放到 E:/R/PCY_protein，",
      "或设置环境变量 PCY_PROTEIN_DIR。"
    )
  }
  normalizePath(existing[1], winslash = "/", mustWork = FALSE)
}

project_dir <- resolve_project_dir()
result_dir <- file.path(project_dir, "results", "PCY_vs_EV")
log_dir <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

# -----------------------------------------------------------------------------
# 2. 读入定量表
# -----------------------------------------------------------------------------
sample_levels <- c("EV1", "EV2", "EV3", "PCY1", "PCY2", "PCY3")

parse_sample_id <- function(text) {
  x <- toupper(as.character(text))
  if (length(x) != 1 || is.na(x) || !nzchar(x)) return(NA_character_)
  if (grepl("PCY", x) && grepl("EV", x)) return(NA_character_)
  m <- regexec("(PCY|EV)[^A-Z0-9]*0*([123])([^0-9]|$)", x)
  r <- regmatches(x, m)[[1]]
  if (length(r) < 3) return(NA_character_)
  paste0(r[2], r[3])
}

quant_class <- function(nm) {
  if (grepl("肽段|肽数|覆盖率|覆盖度|得分|评分", nm)) return("ignore")
  if (grepl("强度|丰度|峰面积|定量", nm)) return("INTENSITY")
  u <- toupper(nm)
  if (grepl("PEPTIDE|\\bPSM\\b|SCORE|PROBABILITY|QVALUE|Q\\.VALUE|PVALUE|P\\.VALUE|FOLD|RATIO|COVERAGE|UNIQUE|SEQUENCE", u, perl = TRUE)) {
    return("ignore")
  }
  if (grepl("LFQ", u)) return("LFQ")
  if (grepl("IBAQ", u)) return("IBAQ")
  if (grepl("INTENSITY|ABUNDANCE|AREA|QUANTITY|\\bQUANT\\b|MS1", u, perl = TRUE)) return("INTENSITY")
  if (grepl("SPECTRAL|MS/MS COUNT|MSMS", u)) return("COUNT")
  "OTHER"
}

norm_header <- function(x) {
  tolower(gsub("[^a-z0-9\u4e00-\u9fff]", "", x, ignore.case = TRUE, perl = TRUE))
}

pick_column <- function(nms, aliases) {
  nn <- norm_header(nms)
  for (a in aliases) {
    hit <- which(nn == a)
    if (length(hit) > 0) return(nms[hit[1]])
  }
  NA_character_
}

guess_gene_column <- function(nms) {
  pick_column(nms, c(
    "genenames", "genename", "genesymbol", "genesymbols", "symbol",
    "pggenes", "geneid", "基因名", "基因符号", "基因名称", "gene", "基因"
  ))
}

guess_protein_column <- function(nms) {
  pick_column(nms, c(
    "majorityproteinids", "proteinids", "proteingroups", "proteingroup",
    "pgproteinaccessions", "uniprotids", "uniprot", "accession",
    "蛋白登录号", "蛋白质登录号", "登录号", "蛋白编号", "蛋白质编号",
    "proteinid", "protein", "entry"
  ))
}

flag_column <- function(nms, aliases) {
  pick_column(nms, aliases)
}

is_flagged <- function(x) {
  u <- toupper(trimws(as.character(x)))
  u %in% c("+", "TRUE", "T", "1", "YES", "Y")
}

to_numeric <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
  x <- gsub(",", "", trimws(as.character(x)), fixed = TRUE)
  x[x %in% c("", "NA", "NaN", "N/A", "n.a.", "Filtered", "Na")] <- NA
  suppressWarnings(as.numeric(x))
}

sample_hits <- function(nms) {
  ids <- vapply(nms, parse_sample_id, character(1))
  cls <- vapply(nms, quant_class, character(1))
  ok <- !is.na(ids) & cls != "ignore"
  data.frame(column = nms[ok], sample = ids[ok], class = cls[ok], stringsAsFactors = FALSE)
}

choose_quant_columns <- function(df) {
  info <- sample_hits(names(df))
  if (nrow(info) == 0) return(NULL)
  pref <- c("LFQ", "IBAQ", "INTENSITY", "COUNT", "OTHER")
  pick <- NULL
  for (cl in pref) {
    sub <- info[info$class == cl, , drop = FALSE]
    if (all(sample_levels %in% sub$sample)) {
      pick <- sub
      break
    }
  }
  if (is.null(pick)) {
    counts <- tapply(info$sample, info$class, function(s) length(unique(s)))
    cl <- names(which.max(counts))
    pick <- info[info$class == cl, , drop = FALSE]
  }
  pick <- pick[order(match(pick$sample, sample_levels), nchar(pick$column)), , drop = FALSE]
  pick <- pick[!duplicated(pick$sample), , drop = FALSE]
  missing <- setdiff(sample_levels, pick$sample)
  attr(pick, "missing") <- missing
  attr(pick, "class") <- unique(as.character(pick$class))
  pick
}

looks_like_header_cell <- function(x) {
  x <- as.character(x)
  if (length(x) != 1 || is.na(x) || !nzchar(trimws(x))) return(FALSE)
  if (grepl("基因|蛋白|登录|强度|丰度|峰面积|描述|名称|肽段|得分|分子量", x)) return(TRUE)
  nx <- norm_header(x)
  keys <- c(
    "gene", "protein", "accession", "uniprot", "intensity", "lfq", "ibaq",
    "abundance", "area", "peptide", "description", "symbol", "score",
    "coverage", "sequence", "entry", "name"
  )
  any(vapply(keys, function(k) grepl(k, nx, fixed = TRUE), logical(1)))
}

score_header_row <- function(cells) {
  cells <- trimws(as.character(cells))
  key <- sum(vapply(cells, looks_like_header_cell, logical(1)))
  n_samp <- length(unique(stats::na.omit(vapply(cells, parse_sample_id, character(1)))))
  key + n_samp * 5
}

promote_header <- function(raw) {
  raw <- as.data.frame(raw, stringsAsFactors = FALSE)
  if (nrow(raw) < 1) return(NULL)
  nscan <- min(20L, nrow(raw))
  scores <- vapply(seq_len(nscan), function(i) {
    score_header_row(unlist(raw[i, , drop = TRUE]))
  }, numeric(1))
  if (max(scores) > 0) {
    best_i <- which.max(scores)
    hdr <- trimws(as.character(unlist(raw[best_i, , drop = TRUE])))
    hdr <- sub("^\ufeff", "", hdr)
    empty <- !nzchar(hdr) | is.na(hdr)
    hdr[empty] <- paste0("V", which(empty))
    hdr <- make.unique(hdr)
    if (best_i >= nrow(raw)) return(NULL)
    df <- raw[(best_i + 1L):nrow(raw), , drop = FALSE]
    names(df) <- hdr
    return(df)
  }
  frac_num <- mean(!is.na(to_numeric(as.character(unlist(raw[1, , drop = TRUE])))))
  if (frac_num > 0.5) {
    names(raw) <- paste0("V", seq_len(ncol(raw)))
    return(raw)
  }
  hdr <- trimws(as.character(unlist(raw[1, , drop = TRUE])))
  hdr <- sub("^\ufeff", "", hdr)
  empty <- !nzchar(hdr) | is.na(hdr)
  hdr[empty] <- paste0("V", which(empty))
  if (nrow(raw) < 2) return(NULL)
  df <- raw[-1, , drop = FALSE]
  names(df) <- make.unique(hdr)
  df
}

read_text_flexible <- function(path) {
  first <- readLines(path, n = 8, warn = FALSE, encoding = "UTF-8")
  first <- first[nzchar(first)]
  sep <- "\t"
  if (length(first) > 0) {
    scores <- c(tab = sum(grepl("\t", first)), comma = sum(grepl(",", first)), semi = sum(grepl(";", first)))
    sep <- switch(names(which.max(scores)), tab = "\t", comma = ",", semi = ";", "\t")
  }
  read_one <- function(enc) {
    utils::read.delim(
      path, sep = sep, header = FALSE, check.names = FALSE, quote = "\"",
      comment.char = "", stringsAsFactors = FALSE, fileEncoding = enc,
      blank.lines.skip = FALSE, fill = TRUE
    )
  }
  raw <- tryCatch(suppressWarnings(read_one("UTF-8")), error = function(e) NULL)
  df <- if (is.null(raw)) NULL else promote_header(raw)
  raw2 <- tryCatch(suppressWarnings(read_one("GB18030")), error = function(e) NULL)
  df2 <- if (is.null(raw2)) NULL else promote_header(raw2)
  q1 <- if (is.null(df)) -1 else score_header_row(names(df))
  q2 <- if (is.null(df2)) -1 else score_header_row(names(df2))
  if (q2 > q1) df2 else df
}

read_excel_flexible <- function(path) {
  if (!has_pkg("readxl")) {
    log_msg("未安装 readxl，跳过 Excel: ", path)
    return(NULL)
  }
  sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) character())
  best <- NULL
  best_n <- -1L
  best_sheet <- NA_character_
  for (sh in sheets) {
    raw <- tryCatch(
      readxl::read_excel(path, sheet = sh, col_names = FALSE, .name_repair = "minimal"),
      error = function(e) NULL
    )
    if (is.null(raw)) next
    df <- promote_header(as.data.frame(raw, stringsAsFactors = FALSE))
    if (is.null(df)) next
    n <- score_header_row(names(df)) + nrow(df) / 1e6
    if (n > best_n) {
      best <- df
      best_n <- n
      best_sheet <- sh
    }
  }
  if (!is.null(best)) attr(best, "sheet") <- best_sheet
  best
}

read_input_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  df <- if (ext %in% c("xlsx", "xls")) read_excel_flexible(path) else read_text_flexible(path)
  if (is.null(df) || ncol(df) < 4 || nrow(df) < 1) return(NULL)
  quant <- choose_quant_columns(df)
  if (is.null(quant)) return(NULL)
  attr(df, "quant") <- quant
  attr(df, "n_samples") <- length(unique(quant$sample))
  df
}

list_input_files <- function(dir) {
  files <- list.files(dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[grepl("\\.(csv|tsv|txt|xlsx|xls|tab)$", files, ignore.case = TRUE)]
  files <- files[!grepl("(^|/)results(/|$)", files)]
  files <- files[!grepl("normalized_|_DE_|sessionInfo|pipeline_", basename(files))]
  files
}

filename_sample_id <- function(path) {
  base <- toupper(gsub("[^A-Za-z0-9]", "", tools::file_path_sans_ext(basename(path))))
  if (!grepl("^(PCY|EV)0*[123]$", base)) return(NA_character_)
  parse_sample_id(base)
}

read_table_generic <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) return(read_excel_flexible(path))
  read_text_flexible(path)
}

pick_abundance_column <- function(df, exclude = character()) {
  nms <- setdiff(names(df), exclude[!is.na(exclude)])
  if (length(nms) == 0) return(NA_character_)
  cls <- setNames(vapply(nms, quant_class, character(1)), nms)
  is_num <- setNames(vapply(nms, function(nm) {
    sum(is.finite(to_numeric(df[[nm]]))) >= 5
  }, logical(1)), nms)
  pref <- c("LFQ", "IBAQ", "INTENSITY", "COUNT", "OTHER")
  for (cl in pref) {
    cols <- nms[cls[nms] == cl & is_num[nms]]
    if (length(cols) == 0) next
    if (length(cols) == 1) return(cols)
    med <- vapply(cols, function(nm) {
      x <- to_numeric(df[[nm]])
      x <- x[is.finite(x) & x > 0]
      if (length(x) == 0) -Inf else stats::median(x)
    }, numeric(1))
    return(cols[which.max(med)])
  }
  NA_character_
}

extract_sample_abundance <- function(df, sample_id) {
  gene_col <- guess_gene_column(names(df))
  prot_col <- guess_protein_column(names(df))
  val_col <- pick_abundance_column(df, exclude = c(gene_col, prot_col))
  if (is.na(val_col)) {
    stop(sample_id, " 里没有找到强度列。列名: ", paste(names(df), collapse = ", "))
  }
  log_msg(
    sample_id, " | gene: ", ifelse(is.na(gene_col), "(none)", gene_col),
    " | protein: ", ifelse(is.na(prot_col), "(none)", prot_col),
    " | abundance: ", val_col
  )
  gene <- if (!is.na(gene_col)) clean_symbol(df[[gene_col]]) else rep(NA_character_, nrow(df))
  protein <- if (!is.na(prot_col)) extract_accession(df[[prot_col]]) else rep(NA_character_, nrow(df))
  value <- to_numeric(df[[val_col]])
  drop <- rep(FALSE, nrow(df))
  for (col in c(
    flag_column(names(df), c("reverse", "reversed")),
    flag_column(names(df), c("potentialcontaminant", "contaminant")),
    flag_column(names(df), c("onlyidentifiedbysite"))
  )) {
    if (!is.na(col)) drop <- drop | is_flagged(df[[col]])
  }
  id <- protein
  empty_id <- is.na(id) | !nzchar(id)
  id[empty_id] <- gene[empty_id]
  ok <- !drop & !is.na(id) & nzchar(id) & is.finite(value)
  id <- id[ok]
  gene <- gene[ok]
  protein <- protein[ok]
  value <- value[ok]
  if (length(id) == 0) stop(sample_id, " 没有可用的蛋白定量行。")
  if (any(duplicated(id))) {
    ord <- order(value, decreasing = TRUE)
    id <- id[ord]
    gene <- gene[ord]
    protein <- protein[ord]
    value <- value[ord]
    keep <- !duplicated(id)
    log_msg(sample_id, " | duplicate IDs collapsed: ", sum(!keep))
    id <- id[keep]
    gene <- gene[keep]
    protein <- protein[keep]
    value <- value[keep]
  }
  data.frame(id = id, gene = gene, protein = protein, value = value, stringsAsFactors = FALSE)
}

merge_sample_files <- function(files) {
  pieces <- lapply(files, function(f) {
    sid <- filename_sample_id(f)
    df <- read_table_generic(f)
    if (is.null(df) || nrow(df) < 1) stop("无法读取样品文件: ", f)
    extract_sample_abundance(df, sid)
  })
  names(pieces) <- vapply(files, filename_sample_id, character(1))
  ids <- unique(unlist(lapply(pieces, function(p) p$id), use.names = FALSE))
  gene <- setNames(rep(NA_character_, length(ids)), ids)
  protein <- gene
  mat <- matrix(NA_real_, nrow = length(ids), ncol = length(sample_levels),
                dimnames = list(ids, sample_levels))
  for (sid in names(pieces)) {
    p <- pieces[[sid]]
    hit <- match(p$id, ids)
    mat[hit, sid] <- p$value
    fill_gene <- is.na(gene[p$id]) & !is.na(p$gene) & nzchar(p$gene)
    gene[p$id[fill_gene]] <- p$gene[fill_gene]
    fill_pro <- is.na(protein[p$id]) & !is.na(p$protein) & nzchar(p$protein)
    protein[p$id[fill_pro]] <- p$protein[fill_pro]
  }
  out <- data.frame(
    Gene.names = unname(gene[ids]),
    Majority.protein.IDs = unname(protein[ids]),
    mat,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  quant <- data.frame(column = sample_levels, sample = sample_levels, class = "INTENSITY", stringsAsFactors = FALSE)
  attr(quant, "missing") <- character()
  attr(quant, "class") <- "per-sample file"
  attr(out, "quant") <- quant
  out
}

load_per_sample_files <- function(files) {
  ids <- vapply(files, filename_sample_id, character(1))
  hit <- files[!is.na(ids)]
  ids <- ids[!is.na(ids)]
  if (length(hit) == 0) return(NULL)
  if (any(duplicated(ids))) {
    log_msg("多个文件对应同一样品，保留较大的文件: ", paste(basename(hit[duplicated(ids)]), collapse = ", "))
    ord <- order(file.info(hit)$size, decreasing = TRUE)
    hit <- hit[ord]
    ids <- ids[ord]
    hit <- hit[!duplicated(ids)]
    ids <- ids[!duplicated(ids)]
  }
  missing <- setdiff(sample_levels, ids)
  if (length(missing) > 0) {
    log_msg("按文件名识别到的样品还不齐，缺少: ", paste(missing, collapse = ", "))
    return(NULL)
  }
  hit <- hit[match(sample_levels, ids)]
  log_msg("Merging one file per sample: ", paste(basename(hit), collapse = ", "))
  df <- merge_sample_files(hit)
  list(df = df, path = paste(hit, collapse = "; "), quant = attr(df, "quant"))
}

load_wide_table <- function(files) {
  best <- NULL
  best_score <- -1
  best_path <- NA_character_
  for (f in files) {
    if (!is.na(filename_sample_id(f))) next
    log_msg("Scanning ", f)
    df <- tryCatch(read_input_file(f), error = function(e) {
      log_msg("  read failed: ", e$message)
      NULL
    })
    if (is.null(df)) next
    score <- attr(df, "n_samples")
    if (grepl("proteingroups", basename(f), ignore.case = TRUE)) score <- score + 0.1
    log_msg("  recognized samples: ", score)
    if (score > best_score) {
      best <- df
      best_score <- score
      best_path <- f
    }
  }
  if (is.null(best)) return(NULL)
  quant <- attr(best, "quant")
  if (length(attr(quant, "missing")) > 0) {
    log_msg(
      "宽表 ", best_path, " 缺少样品列: ", paste(attr(quant, "missing"), collapse = ", ")
    )
    return(NULL)
  }
  log_msg("Using wide table: ", best_path)
  if (!is.null(attr(best, "sheet"))) log_msg("Excel sheet: ", attr(best, "sheet"))
  log_msg("Quantification class: ", paste(attr(quant, "class"), collapse = ", "))
  log_msg(paste(quant$sample, quant$column, sep = " <- ", collapse = " | "))
  list(df = best, path = best_path, quant = quant)
}

load_quant_table <- function(dir) {
  files <- list_input_files(dir)
  if (length(files) == 0) {
    stop("目录中没有 csv/tsv/txt/xlsx 定量表: ", dir)
  }
  wide <- load_wide_table(files)
  if (!is.null(wide)) return(wide)
  per <- load_per_sample_files(files)
  if (!is.null(per)) return(per)
  stop(
    "没有读到六个样品的定量。请提供 EV1、EV2、EV3、PCY1、PCY2、PCY3 六个文件，",
    "或一张列名包含这六个样品的宽表。"
  )
}

# -----------------------------------------------------------------------------
# 3. 基因名、过滤、标准化、缺失值填补
# -----------------------------------------------------------------------------
clean_symbol <- function(x) {
  vapply(as.character(x), function(one) {
    if (length(one) != 1 || is.na(one)) return(NA_character_)
    one <- trimws(one)
    if (!nzchar(one) || one %in% c("-", ".", "NA", "NaN")) return(NA_character_)
    if (grepl("|", one, fixed = TRUE)) {
      bits <- strsplit(one, "|", fixed = TRUE)[[1]]
      last <- bits[length(bits)]
      gm <- regmatches(last, regexec("^([A-Za-z0-9_.-]+)_", last))[[1]]
      if (length(gm) >= 2) return(gm[2])
    }
    parts <- trimws(unlist(strsplit(one, "[;,/]+")))
    parts <- parts[nzchar(parts) & !parts %in% c("-", ".", "NA")]
    if (length(parts) == 0) return(NA_character_)
    parts[1]
  }, character(1), USE.NAMES = FALSE)
}

extract_accession <- function(x) {
  vapply(as.character(x), function(one) {
    if (length(one) != 1 || is.na(one)) return(NA_character_)
    one <- trimws(strsplit(one, "[;,]")[[1]][1])
    if (!nzchar(one) || is.na(one)) return(NA_character_)
    if (grepl("|", one, fixed = TRUE)) {
      bits <- strsplit(one, "|", fixed = TRUE)[[1]]
      if (length(bits) >= 2) return(bits[2])
    }
    one
  }, character(1), USE.NAMES = FALSE)
}

is_uniprot <- function(x) {
  acc <- sub("-\\d+$", "", x)
  grepl("^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])$", acc, perl = TRUE) & !is.na(x)
}

map_accessions_to_symbols <- function(ids) {
  out <- ids
  acc <- sub("-\\d+$", "", ids)
  up <- is_uniprot(ids)
  if (any(up)) {
    mp <- tryCatch(
      suppressMessages(AnnotationDbi::select(
        org.Hs.eg.db, keys = unique(acc[up]), keytype = "UNIPROT", columns = "SYMBOL"
      )),
      error = function(e) NULL
    )
    if (!is.null(mp) && nrow(mp) > 0) {
      mp <- mp[!is.na(mp$SYMBOL) & !duplicated(mp$UNIPROT), , drop = FALSE]
      hit <- mp$SYMBOL[match(acc, mp$UNIPROT)]
      replace <- up & !is.na(hit)
      out[replace] <- hit[replace]
    }
  }
  ensp <- grepl("^ENSP[0-9]", ids) & !is.na(ids)
  if (any(ensp)) {
    mp <- tryCatch(
      suppressMessages(clusterProfiler::bitr(
        unique(ids[ensp]), fromType = "ENSEMBLPROT", toType = "SYMBOL", OrgDb = org.Hs.eg.db
      )),
      error = function(e) NULL
    )
    if (!is.null(mp) && nrow(mp) > 0) {
      mp <- mp[!duplicated(mp$ENSEMBLPROT), , drop = FALSE]
      hit <- mp$SYMBOL[match(ids, mp$ENSEMBLPROT)]
      replace <- ensp & !is.na(hit)
      out[replace] <- hit[replace]
    }
  }
  out
}

collapse_by_gene <- function(gene, protein_id, mat) {
  ok <- !is.na(gene) & nzchar(gene)
  gene <- gene[ok]
  protein_id <- protein_id[ok]
  mat <- mat[ok, , drop = FALSE]
  if (length(gene) == 0) stop("没有可用的基因名或蛋白编号。")
  if (!any(duplicated(gene))) {
    rownames(mat) <- gene
    return(list(mat = mat, protein_id = stats::setNames(protein_id, gene), n_collapsed = 0L))
  }
  spl <- split(seq_along(gene), gene)
  keep <- vapply(spl, function(idx) {
    if (length(idx) == 1) return(idx)
    sub <- mat[idx, , drop = FALSE]
    nval <- rowSums(is.finite(sub))
    mu <- rowMeans(sub, na.rm = TRUE)
    mu[!is.finite(mu)] <- -Inf
    idx[order(nval, mu, decreasing = TRUE)[1]]
  }, integer(1))
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- gene[keep]
  list(
    mat = mat,
    protein_id = stats::setNames(protein_id[keep], gene[keep]),
    n_collapsed = length(gene) - length(keep)
  )
}

filter_missing <- function(mat) {
  n_ev <- rowSums(is.finite(mat[, c("EV1", "EV2", "EV3"), drop = FALSE]))
  n_pcy <- rowSums(is.finite(mat[, c("PCY1", "PCY2", "PCY3"), drop = FALSE]))
  keep <- n_ev >= min_valid_in_group | n_pcy >= min_valid_in_group
  list(mat = mat[keep, , drop = FALSE], n_ev = n_ev[keep], n_pcy = n_pcy[keep])
}

looks_like_log2 <- function(mat) {
  x <- as.numeric(mat)
  x <- x[is.finite(x)]
  if (length(x) < 20) return(FALSE)
  if (any(x < 0)) return(TRUE)
  if (max(x) > 80 || stats::median(x) > 40) return(FALSE)
  prop_integer <- mean(abs(x - round(x)) < 1e-8)
  if (prop_integer > 0.9) return(FALSE)
  TRUE
}

median_center <- function(mat) {
  med <- apply(mat, 2, stats::median, na.rm = TRUE)
  if (any(!is.finite(med))) stop("某个样品在过滤后没有可用定量值，无法做中位数标准化。")
  sweep(mat, 2, med, "-")
}

impute_perseus <- function(mat, width = 0.3, shift = 1.8) {
  set.seed(20260927)
  n_imp <- 0L
  for (j in seq_len(ncol(mat))) {
    x <- mat[, j]
    ok <- is.finite(x)
    if (sum(ok) < 5 || all(ok)) next
    mu <- mean(x[ok])
    sdv <- stats::sd(x[ok])
    if (!is.finite(sdv) || sdv == 0) sdv <- 1e-3
    draws <- stats::rnorm(sum(!ok), mean = mu - shift * sdv, sd = max(width * sdv, 1e-6))
    draws <- pmin(draws, min(x[ok]) - 1e-6)
    x[!ok] <- draws
    n_imp <- n_imp + sum(!ok)
    mat[, j] <- x
  }
  attr(mat, "n_imputed") <- n_imp
  mat
}

build_matrix <- function(df, quant) {
  gene_col <- guess_gene_column(names(df))
  prot_col <- guess_protein_column(names(df))
  if (is.na(gene_col) && is.na(prot_col)) {
    gene_col <- names(df)[1]
    log_msg("未找到基因名列，改用第一列: ", gene_col)
  }
  if (!is.na(gene_col)) log_msg("Gene column: ", gene_col)
  if (!is.na(prot_col)) log_msg("Protein column: ", prot_col)

  rev_col <- flag_column(names(df), c("reverse", "reversed"))
  con_col <- flag_column(names(df), c("potentialcontaminant", "contaminant"))
  site_col <- flag_column(names(df), c("onlyidentifiedbysite"))
  drop <- rep(FALSE, nrow(df))
  for (col in c(rev_col, con_col, site_col)) {
    if (!is.na(col)) drop <- drop | is_flagged(df[[col]])
  }
  if (any(drop)) log_msg("Removed contaminant/reverse/site-only rows: ", sum(drop))
  df <- df[!drop, , drop = FALSE]

  gene <- if (!is.na(gene_col)) clean_symbol(df[[gene_col]]) else rep(NA_character_, nrow(df))
  protein_id <- if (!is.na(prot_col)) extract_accession(df[[prot_col]]) else gene
  empty <- is.na(gene) | !nzchar(gene)
  gene[empty] <- protein_id[empty]
  gene <- map_accessions_to_symbols(gene)

  mat <- as.matrix(as.data.frame(lapply(quant$column, function(col) to_numeric(df[[col]])), check.names = FALSE))
  colnames(mat) <- quant$sample
  mat <- mat[, sample_levels, drop = FALSE]
  mat[!is.finite(mat) | mat <= 0] <- NA

  collapsed <- collapse_by_gene(gene, ifelse(is.na(protein_id), gene, protein_id), mat)
  log_msg("Rows after contaminant filter: ", nrow(df), " | unique genes: ", nrow(collapsed$mat),
          " | collapsed duplicates: ", collapsed$n_collapsed)

  filt <- filter_missing(collapsed$mat)
  log_msg(
    "Kept proteins with >= ", min_valid_in_group, " values in EV or in PCY: ",
    nrow(filt$mat), " (removed ", nrow(collapsed$mat) - nrow(filt$mat), ")"
  )
  if (nrow(filt$mat) < 10) stop("过滤后蛋白少于 10 个，请检查定量列是不是强度而不是肽段数。")

  log_mat <- filt$mat
  if (looks_like_log2(log_mat)) {
    log_msg("Values look like log2 already; skip log2 transform.")
  } else {
    log_msg("log2 transform of linear intensities.")
    log_mat <- log2(log_mat)
  }
  log_mat <- median_center(log_mat)
  n_ev <- filt$n_ev
  n_pcy <- filt$n_pcy
  names(n_ev) <- rownames(log_mat)
  names(n_pcy) <- rownames(log_mat)
  imputed <- impute_perseus(log_mat)
  if (any(!is.finite(imputed))) {
    for (j in seq_len(ncol(imputed))) {
      x <- imputed[, j]
      bad <- !is.finite(x)
      if (any(bad) && any(!bad)) x[bad] <- min(x[!bad]) - 1
      imputed[, j] <- x
    }
  }
  log_msg("Imputed missing log2 values (Perseus downshift): ", attr(imputed, "n_imputed"))
  list(
    log_mat = imputed,
    n_ev = n_ev,
    n_pcy = n_pcy,
    protein_id = collapsed$protein_id[rownames(imputed)],
    n_imputed = attr(imputed, "n_imputed")
  )
}

# -----------------------------------------------------------------------------
# 4. limma：PCY - EV
# -----------------------------------------------------------------------------
run_limma <- function(log_mat) {
  group <- factor(rep(c("EV", "PCY"), each = 3), levels = c("EV", "PCY"))
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- c("EV", "PCY")
  fit <- limma::lmFit(log_mat, design)
  cont <- limma::makeContrasts(PCY_vs_EV = PCY - EV, levels = design)
  fit2 <- limma::contrasts.fit(fit, cont)
  fit2 <- tryCatch(
    limma::eBayes(fit2, trend = TRUE, robust = TRUE),
    error = function(e) {
      log_msg("eBayes trend/robust failed, fallback to ordinary eBayes: ", e$message)
      limma::eBayes(fit2)
    }
  )
  limma::topTable(fit2, coef = 1, number = Inf, sort.by = "none")
}

assemble_de <- function(tt, prep) {
  gene <- rownames(tt)
  de <- data.frame(
    gene = gene,
    protein_id = unname(prep$protein_id[gene]),
    log2FC = tt$logFC,
    FC = 2^tt$logFC,
    AveExpr = tt$AveExpr,
    t = tt$t,
    pvalue = tt$P.Value,
    padj = tt$adj.P.Val,
    B = tt$B,
    n_EV = unname(prep$n_ev[gene]),
    n_PCY = unname(prep$n_pcy[gene]),
    stringsAsFactors = FALSE
  )
  both <- de$n_EV >= min_valid_in_group & de$n_PCY >= min_valid_in_group
  de$regulation <- "NS"
  de$regulation[both & de$pvalue < pvalue_cutoff & de$log2FC > 0] <- "Up"
  de$regulation[both & de$pvalue < pvalue_cutoff & de$log2FC < 0] <- "Down"
  de$regulation[de$n_EV == 0 & de$n_PCY >= min_valid_in_group] <- "PCY only"
  de$regulation[de$n_PCY == 0 & de$n_EV >= min_valid_in_group] <- "EV only"
  de <- de[order(de$pvalue, -abs(de$log2FC)), , drop = FALSE]
  rownames(de) <- NULL
  de
}

select_up <- function(de, fc) {
  keep <- (de$regulation == "Up" & de$FC > fc) | de$regulation == "PCY only"
  sub <- de[keep, , drop = FALSE]
  sub[order(sub$FC, decreasing = TRUE), , drop = FALSE]
}

# -----------------------------------------------------------------------------
# 5. 火山图与气泡图
# -----------------------------------------------------------------------------
save_gg <- function(plot, path_stub, width = 8, height = 6) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  tryCatch(
    ggplot2::ggsave(paste0(path_stub, ".pdf"), plot, width = width, height = height),
    error = function(e) log_msg("pdf ggsave failed: ", e$message)
  )
  tryCatch(
    ggplot2::ggsave(paste0(path_stub, ".png"), plot, width = width, height = height, dpi = 300),
    error = function(e) log_msg("png ggsave failed: ", e$message)
  )
}

plot_volcano <- function(de, outfile) {
  df <- de
  df$y <- -log10(pmax(df$pvalue, 1e-300))
  lev <- c("Up", "Down", "PCY only", "EV only", "NS")
  df$regulation <- factor(df$regulation, levels = lev)
  cols <- c(
    "Up" = "#D62828", "Down" = "#1D4E89", "PCY only" = "#F77F00",
    "EV only" = "#2A9D8F", "NS" = "grey75"
  )
  pcy_only <- df[df$regulation == "PCY only", , drop = FALSE]
  pcy_only <- pcy_only[order(pcy_only$log2FC, decreasing = TRUE), , drop = FALSE]
  label_genes <- unique(c(
    utils::head(df$gene[df$regulation == "Up"], 8),
    utils::head(df$gene[df$regulation == "Down"], 8),
    utils::head(pcy_only$gene, 6)
  ))
  df$label <- ifelse(df$gene %in% label_genes, df$gene, NA_character_)
  n_up <- sum(df$regulation == "Up")
  n_down <- sum(df$regulation == "Down")
  n_po <- sum(df$regulation == "PCY only")
  shapes <- c("Up" = 16, "Down" = 16, "PCY only" = 17, "EV only" = 15, "NS" = 16)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = y, color = regulation, shape = regulation)) +
    ggplot2::geom_point(alpha = 0.8, size = 1.8) +
    ggplot2::scale_color_manual(values = cols, drop = TRUE) +
    ggplot2::scale_shape_manual(values = shapes, drop = TRUE) +
    ggplot2::geom_vline(xintercept = 0, color = "grey40", linewidth = 0.3) +
    ggplot2::geom_vline(
      xintercept = c(-log2(1.5), -log2(1.25), log2(1.25), log2(1.5)),
      linetype = 2, color = "grey45", linewidth = 0.3
    ) +
    ggplot2::geom_hline(yintercept = -log10(pvalue_cutoff), linetype = 2, color = "grey45", linewidth = 0.3) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label), size = 3, max.overlaps = 40, force = 2,
      box.padding = 0.45, point.padding = 0.3, min.segment.length = 0,
      segment.size = 0.25, seed = 1,
      na.rm = TRUE, show.legend = FALSE, color = "grey15"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = "PCY vs EV",
      subtitle = sprintf(
        "Up %d, Down %d, PCY only %d (p < %s). Dashed lines: p = %s, FC = 1.25 and 1.5",
        n_up, n_down, n_po, pvalue_cutoff, pvalue_cutoff
      ),
      x = "log2 fold change (PCY / EV)",
      y = "-log10(p value)",
      color = NULL,
      shape = NULL
    )
  save_gg(p, outfile, width = 8, height = 6.5)
}

ratio_to_num <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(p) {
    if (length(p) != 2) return(NA_real_)
    as.numeric(p[1]) / as.numeric(p[2])
  }, numeric(1))
}

wrap_text <- function(x, width = 42) {
  vapply(as.character(x), function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

note_empty <- function(stub, msg) {
  dir.create(dirname(stub), recursive = TRUE, showWarnings = FALSE)
  writeLines(msg, paste0(stub, "_EMPTY.txt"))
}

plot_bubble <- function(df, title, stub) {
  if (is.null(df) || nrow(df) == 0 || !"Description" %in% names(df)) {
    note_empty(stub, "no enrichment terms")
    return(invisible(NULL))
  }
  df <- df[is.finite(df$pvalue), , drop = FALSE]
  if (nrow(df) == 0) {
    note_empty(stub, "no enrichment terms")
    return(invisible(NULL))
  }
  if (!"p.adjust" %in% names(df)) df$p.adjust <- df$pvalue
  df <- df[order(df$p.adjust, df$pvalue), , drop = FALSE]
  nshow <- min(15L, nrow(df))
  df <- df[seq_len(nshow), , drop = FALSE]
  df$gene_ratio <- if ("GeneRatio" %in% names(df)) ratio_to_num(df$GeneRatio) else NA_real_
  if (all(!is.finite(df$gene_ratio)) && "Count" %in% names(df)) df$gene_ratio <- df$Count
  if (!"Count" %in% names(df)) df$Count <- 1
  df$label <- wrap_text(df$Description)
  df$label <- factor(df$label, levels = rev(unique(df$label)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = gene_ratio, y = label, size = Count, color = p.adjust)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_color_gradient(low = "#B2182B", high = "#2166AC", name = "Adjusted p") +
    ggplot2::scale_size(range = c(2.5, 8), name = "Count") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9)) +
    ggplot2::labs(title = title, x = "Gene ratio", y = NULL)
  save_gg(p, stub, width = 10, height = max(5, 0.42 * nshow + 1.8))
}

# -----------------------------------------------------------------------------
# 6. GO / KEGG / 通路（ORA 气泡图）
# -----------------------------------------------------------------------------
map_symbols_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  empty <- data.frame(gene = character(), entrez = character(), stringsAsFactors = FALSE)
  if (length(symbols) == 0) return(empty)
  mapped <- empty
  hit <- tryCatch(
    suppressMessages(clusterProfiler::bitr(
      symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db
    )),
    error = function(e) NULL
  )
  if (!is.null(hit) && nrow(hit) > 0) {
    hit <- hit[!duplicated(hit$SYMBOL), , drop = FALSE]
    mapped <- data.frame(gene = hit$SYMBOL, entrez = as.character(hit$ENTREZID), stringsAsFactors = FALSE)
  }
  missed <- setdiff(symbols, mapped$gene)
  if (length(missed) > 0) {
    al <- tryCatch(
      suppressMessages(AnnotationDbi::select(
        org.Hs.eg.db, keys = missed, keytype = "ALIAS", columns = "ENTREZID"
      )),
      error = function(e) NULL
    )
    if (!is.null(al) && nrow(al) > 0) {
      al <- al[!is.na(al$ENTREZID) & !duplicated(al$ALIAS), , drop = FALSE]
      mapped <- rbind(
        mapped,
        data.frame(gene = al$ALIAS, entrez = as.character(al$ENTREZID), stringsAsFactors = FALSE)
      )
    }
  }
  mapped[!duplicated(mapped$gene), , drop = FALSE]
}

msig_hallmark_map <- function() {
  if (!has_pkg("msigdbr")) return(NULL)
  msig <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e) tryCatch(msigdbr::msigdbr(species = "Homo sapiens", category = "H"), error = function(e2) NULL)
  )
  if (is.null(msig) || nrow(msig) == 0) return(NULL)
  gene_col <- intersect(c("ncbi_gene", "entrez_gene"), names(msig))[1]
  if (is.na(gene_col)) return(NULL)
  out <- msig[, c("gs_name", gene_col)]
  names(out) <- c("term", "gene")
  out$gene <- as.character(out$gene)
  out
}

is_download_error <- function(msg) {
  !is.na(msg) && nzchar(msg) && grepl(
    "cannot open|HTTP status|timed out|Timeout|Could not resolve|network",
    msg, ignore.case = TRUE
  )
}

enrich_terms <- function(strict_fun, relax_fun, label) {
  err1 <- NA_character_
  obj <- tryCatch(strict_fun(), error = function(e) {
    err1 <<- conditionMessage(e)
    log_msg(label, " failed: ", err1)
    NULL
  })
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    return(list(obj = obj, relaxed = FALSE, error = NA_character_))
  }
  if (is_download_error(err1)) {
    return(list(obj = NULL, relaxed = TRUE, error = err1))
  }
  err2 <- NA_character_
  obj <- tryCatch(relax_fun(), error = function(e) {
    err2 <<- conditionMessage(e)
    log_msg(label, " relaxed failed: ", err2)
    NULL
  })
  list(obj = obj, relaxed = TRUE, error = ifelse(is.na(err2), err1, err2))
}

write_enrich_bubble <- function(res, stub, title) {
  obj <- res$obj
  df <- if (is.null(obj)) NULL else as.data.frame(obj)
  if (is.null(df) || nrow(df) == 0) {
    plot_bubble(NULL, title, paste0(stub, "_bubble"))
    return(invisible(NULL))
  }
  if (isTRUE(res$relaxed)) {
    df <- df[order(df$pvalue), , drop = FALSE]
    df <- utils::head(df, 50)
    title <- paste0(title, " (no term at adjusted p < 0.05; top terms)")
  }
  dir.create(dirname(stub), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, paste0(stub, ".csv"), row.names = FALSE)
  plot_bubble(df, title, paste0(stub, "_bubble"))
}

set_readable <- function(obj) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) return(obj)
  tryCatch(
    clusterProfiler::setReadable(obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
    error = function(e) obj
  )
}

run_fc_enrichment <- function(genes, universe, outdir, tag, label) {
  go_dir <- file.path(outdir, "GO")
  kg_dir <- file.path(outdir, "KEGG")
  pw_dir <- file.path(outdir, "Pathway")
  dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(kg_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pw_dir, recursive = TRUE, showWarnings = FALSE)
  pref <- paste0(tag, "_ORA_")

  gene_map <- map_symbols_to_entrez(genes)
  univ_map <- map_symbols_to_entrez(universe)
  universe_entrez <- unique(univ_map$entrez)
  entrez <- intersect(unique(gene_map$entrez), universe_entrez)
  log_msg(tag, " genes ", length(unique(genes)), " mapped ", length(entrez),
          " | universe mapped ", length(universe_entrez))
  if (length(entrez) < 3) {
    msg <- paste("mapped Entrez genes < 3:", length(entrez))
    note_empty(file.path(go_dir, paste0(pref, "GO")), msg)
    note_empty(file.path(kg_dir, paste0(pref, "KEGG")), msg)
    note_empty(file.path(pw_dir, paste0(pref, "Pathway")), msg)
    return(invisible(NULL))
  }

  for (ont in c("BP", "MF", "CC")) {
    res <- enrich_terms(
      function() clusterProfiler::enrichGO(
        gene = entrez, universe = universe_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
        ont = ont, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
        minGSSize = 3, maxGSSize = 500, readable = TRUE
      ),
      function() clusterProfiler::enrichGO(
        gene = entrez, universe = universe_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
        ont = ont, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
        minGSSize = 3, maxGSSize = 500, readable = TRUE
      ),
      paste("GO", ont, tag)
    )
    write_enrich_bubble(
      res,
      file.path(go_dir, paste0(pref, "GO_", ont)),
      paste(label, "| GO", ont)
    )
  }

  kegg <- enrich_terms(
    function() clusterProfiler::enrichKEGG(
      gene = entrez, universe = universe_entrez, organism = kegg_organism,
      pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 3, maxGSSize = 500
    ),
    function() clusterProfiler::enrichKEGG(
      gene = entrez, universe = universe_entrez, organism = kegg_organism,
      pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 3, maxGSSize = 500
    ),
    paste("KEGG", tag)
  )
  kegg$obj <- set_readable(kegg$obj)
  write_enrich_bubble(kegg, file.path(kg_dir, paste0(pref, "KEGG")), paste(label, "| KEGG"))

  if (has_pkg("ReactomePA")) {
    rea <- enrich_terms(
      function() ReactomePA::enrichPathway(
        gene = entrez, universe = universe_entrez, organism = reactome_organism,
        pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 3, maxGSSize = 500, readable = TRUE
      ),
      function() ReactomePA::enrichPathway(
        gene = entrez, universe = universe_entrez, organism = reactome_organism,
        pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 3, maxGSSize = 500, readable = TRUE
      ),
      paste("Reactome", tag)
    )
    write_enrich_bubble(
      rea, file.path(pw_dir, paste0(pref, "Reactome_pathway")), paste(label, "| Reactome pathway")
    )
  } else {
    note_empty(file.path(pw_dir, paste0(pref, "Reactome_pathway_bubble")), "ReactomePA not installed")
  }

  if (isTRUE(wiki_unavailable)) {
    note_empty(file.path(pw_dir, paste0(pref, "WikiPathways_pathway_bubble")), "WikiPathways download unavailable")
  } else {
    old_timeout <- options(timeout = 40)
    wp <- enrich_terms(
      function() clusterProfiler::enrichWP(
        gene = entrez, universe = universe_entrez, organism = wp_organism,
        pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 3, maxGSSize = 500
      ),
      function() clusterProfiler::enrichWP(
        gene = entrez, universe = universe_entrez, organism = wp_organism,
        pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 3, maxGSSize = 500
      ),
      paste("WikiPathways", tag)
    )
    options(old_timeout)
    if (is_download_error(wp$error)) wiki_unavailable <<- TRUE
    wp$obj <- set_readable(wp$obj)
    write_enrich_bubble(
      wp, file.path(pw_dir, paste0(pref, "WikiPathways_pathway")), paste(label, "| WikiPathways")
    )
  }

  hallmark <- msig_hallmark_map()
  if (!is.null(hallmark)) {
    hm <- enrich_terms(
      function() clusterProfiler::enricher(
        entrez, universe = universe_entrez, TERM2GENE = hallmark,
        pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 3, maxGSSize = 500
      ),
      function() clusterProfiler::enricher(
        entrez, universe = universe_entrez, TERM2GENE = hallmark,
        pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 3, maxGSSize = 500
      ),
      paste("Hallmark", tag)
    )
    hm$obj <- set_readable(hm$obj)
    write_enrich_bubble(
      hm, file.path(pw_dir, paste0(pref, "Hallmark_pathway")), paste(label, "| MSigDB Hallmark pathway")
    )
  } else {
    note_empty(file.path(pw_dir, paste0(pref, "Hallmark_pathway_bubble")), "msigdbr not installed or Hallmark map empty")
  }
}

# -----------------------------------------------------------------------------
# 7. 主流程
# -----------------------------------------------------------------------------
log_msg("Project dir: ", project_dir)
loaded <- load_quant_table(project_dir)
prep <- build_matrix(loaded$df, loaded$quant)
utils::write.csv(
  data.frame(sample = sample_levels, group = rep(c("EV", "PCY"), each = 3), column = loaded$quant$column[match(sample_levels, loaded$quant$sample)]),
  file.path(log_dir, "sample_columns.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(gene = rownames(prep$log_mat), protein_id = unname(prep$protein_id), prep$log_mat, check.names = FALSE),
  file.path(result_dir, "normalized_imputed_log2_matrix.csv"),
  row.names = FALSE
)

tt <- run_limma(prep$log_mat)
de <- assemble_de(tt, prep)
utils::write.csv(de, file.path(result_dir, "PCY_vs_EV_DE_all.csv"), row.names = FALSE)
if (has_pkg("writexl")) {
  tryCatch(
    writexl::write_xlsx(de, file.path(result_dir, "PCY_vs_EV_DE_all.xlsx")),
    error = function(e) log_msg("xlsx export failed: ", e$message)
  )
}
log_msg(
  "DE proteins: ", nrow(de),
  " | Up ", sum(de$regulation == "Up"),
  " | Down ", sum(de$regulation == "Down"),
  " | PCY only ", sum(de$regulation == "PCY only"),
  " | EV only ", sum(de$regulation == "EV only"),
  " | padj < ", pvalue_cutoff, ": ", sum(de$padj < pvalue_cutoff, na.rm = TRUE)
)
plot_volcano(de, file.path(result_dir, "PCY_vs_EV_volcano"))

writeLines(
  c(
    "PCY vs EV（IP 蛋白质谱）",
    "",
    "火山图: PCY_vs_EV_volcano.pdf / .png",
    "全部蛋白差异表: PCY_vs_EV_DE_all.csv",
    "  log2FC 与 FC 都是 PCY / EV。FC > 1 表示 PCY 更高。",
    "  regulation = Up / Down 要求两组都至少有 2 个有效值，且 p < 0.05。",
    "  PCY only：PCY 至少 2 个重复有值、EV 三个重复都缺失，缺失值按低丰度填补后再检验。",
    "",
    "富集输入（只取上调）:",
    "  FoldChange/FC_1/    p < 0.05 且 FC > 1 的 Up，加上 PCY only",
    "  FoldChange/FC_1.25/ p < 0.05 且 FC > 1.25 的 Up，加上 PCY only",
    "  FoldChange/FC_1.5/  p < 0.05 且 FC > 1.5 的 Up，加上 PCY only",
    "每个档位:",
    "  GO/       GO 生物过程、分子功能、细胞组分，气泡图",
    "  KEGG/     KEGG 气泡图",
    "  Pathway/  Reactome、WikiPathways、MSigDB Hallmark 气泡图",
    "背景基因集是本次通过过滤的全部蛋白，不是整个人类基因组。",
    "气泡图文件名以 _bubble 结尾，并带有 FC_1 / FC_1.25 / FC_1.5。"
  ),
  file.path(result_dir, "00_READ_ME.txt")
)

for (nm in names(fc_cutoffs)) {
  fc <- unname(fc_cutoffs[[nm]])
  sub <- select_up(de, fc)
  outdir <- file.path(result_dir, "FoldChange", nm)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(sub, file.path(outdir, paste0(nm, "_up_genes.csv")), row.names = FALSE)
  writeLines(
    c(
      paste("comparison: PCY vs EV"),
      paste("filter: (Up AND p <", pvalue_cutoff, "AND FC >", fc, ") OR PCY only"),
      paste("n_genes:", nrow(sub)),
      "Enrichment background: all proteins that passed the missing-value filter."
    ),
    file.path(outdir, paste0("00_", nm, "_filter.txt"))
  )
  log_msg(nm, " upregulated genes for enrichment: ", nrow(sub))
  tryCatch(
    run_fc_enrichment(
      sub$gene, de$gene, outdir, nm,
      sprintf("PCY vs EV | FC > %s", fc)
    ),
    error = function(e) log_msg("ERROR enrichment ", nm, ": ", e$message)
  )
}

writeLines(
  c(
    paste("pvalue_cutoff", pvalue_cutoff),
    paste("fc_cutoffs", paste(names(fc_cutoffs), fc_cutoffs, sep = "=", collapse = ", ")),
    paste("min_valid_in_group", min_valid_in_group),
    paste("input", loaded$path),
    paste("n_proteins", nrow(de)),
    paste("n_imputed_values", prep$n_imputed),
    "normalization: log2 then per-sample median centering",
    "imputation: Perseus normal downshift (shift 1.8, width 0.3), seed 20260927",
    "test: limma eBayes trend+robust, contrast PCY - EV",
    "species: Homo sapiens"
  ),
  file.path(log_dir, "analysis_parameters.txt")
)
base::writeLines(utils::capture.output(utils::sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
log_msg("All done. Results in: ", result_dir)
