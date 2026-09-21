#!/usr/bin/env Rscript
# =============================================================================
# FUSCC_Human_breast.R
# 复旦肿瘤医院 FUSCC 人源乳腺癌：原发灶低表达 vs 肺/骨转移
#
# 数据目录（脚本与结果默认同目录）：
#   E:/R/Human breast cancer/FUSCC
#
# 四个问题（只分析「原发相对更低」= FC > 1）：
#   1) 原发内哪些基因低表达促进肺转移、骨转移
#   2) 器官特异：只肺或只骨
#   3) 轴突导向 / 施旺细胞 / 神经营养 三套签名给原发打神经浸润分
#   4) 三种神经浸润分别与肺转移的关系
#
# 阈值：p < 0.05，FC > 1 与 FC > 1.25。不做 top50–300。
#
# 公开组学矩阵几乎全是**原发灶**；转移信息在临床/预后/随访表里（DFS/RFS/DMFS、
# 复发部位、肺转移、骨转移）。Patient N 的原发 对应 Patient N 预后里的肺/骨转移。
# 两种模式自动切换：
#   A) PAIRED_TISSUE：文件夹里若真有同一患者的原发 + 肺/骨转移表达，做组织一一对应
#      （FC = 该患者转移组织 / 该患者原发）
#   B) CLINICAL_OUTCOME（默认，FUSCC 预后转移走这条）：每位患者一条原发，按该患者
#      随访是否肺/骨转移分组（FC = 无该器官转移 / 有该器官转移）
#
# 不要改 TG_RNAseq_*.R。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")
options(clusterProfiler.download.method = "auto")

# -----------------------------------------------------------------------------
# 0. 包
# -----------------------------------------------------------------------------
cran_required <- c(
  "dplyr", "tidyr", "tibble", "stringr", "ggplot2", "ggrepel",
  "pheatmap", "RColorBrewer", "matrixStats", "writexl", "readxl"
)
cran_optional <- c("ggvenn", "GSVA")
bioc_required <- c(
  "DESeq2", "limma", "clusterProfiler", "org.Hs.eg.db",
  "enrichplot", "AnnotationDbi", "fgsea"
)
bioc_optional <- c("ReactomePA", "msigdbr", "pathview")

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
  if (length(still) > 0) {
    message("可选包未安装，相关分析将跳过: ", paste(still, collapse = ", "))
  }
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
# 1. 路径
# -----------------------------------------------------------------------------
resolve_project_dir <- function() {
  env_dir <- Sys.getenv("FUSCC_DIR", unset = "")
  script_dir <- tryCatch({
    of <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(of)) dirname(normalizePath(sub("^--file=", "", of[1]))) else NA_character_
  }, error = function(e) NA_character_)
  candidates <- c(
    env_dir,
    "E:/R/Human breast cancer/FUSCC",
    "E:\\R\\Human breast cancer\\FUSCC",
    script_dir,
    file.path(getwd(), "FUSCC"),
    getwd()
  )
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  hit_file <- function(d) {
    if (!dir.exists(d)) return(FALSE)
    ff <- list.files(d, full.names = FALSE, ignore.case = TRUE)
    any(grepl("FUSCC|OEP|count|tpm|fpkm|clinic|pheno|matrix|rna|PXD|OMIX|prognos|follow|survival|dfs|rfs|dmfs|预后|随访",
              ff, ignore.case = TRUE))
  }
  for (d in candidates) {
    if (hit_file(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

project_dir <- resolve_project_dir()
result_dir  <- file.path(project_dir, "results_FUSCC_Human_breast")
log_dir     <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("FUSCC_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_msg("Working directory: ", project_dir)
log_msg("Results: ", result_dir)

p_cutoff  <- 0.05
fc_cutoffs <- c("FC_1" = 1, "FC_1.25" = 1.25)

# -----------------------------------------------------------------------------
# 2. 扫描 FUSCC 文件夹：表达矩阵 + 临床表
# -----------------------------------------------------------------------------
list_data_files <- function() {
  ff <- list.files(project_dir, full.names = TRUE, recursive = TRUE)
  ff <- ff[!grepl("results_FUSCC", ff, ignore.case = TRUE)]
  ff <- ff[file.info(ff)$isdir == FALSE]
  ff
}

read_table_any <- function(path, n_max = Inf) {
  ext <- tolower(sub(".*\\.", "", sub("\\.gz$", "", basename(path), ignore.case = TRUE)))
  if (ext %in% c("xlsx", "xls")) {
    if (!has_pkg("readxl")) stop("需要 readxl 读 Excel: ", path)
    sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) "1")
    pick <- sheets[grepl("clinic|pheno|sample|patient|rna|count|expr|sheet1|ratio", sheets, ignore.case = TRUE)]
    if (length(pick) == 0) pick <- sheets[1]
    as.data.frame(readxl::read_excel(path, sheet = pick[1], guess_max = 8000), stringsAsFactors = FALSE)
  } else {
    open_text <- function() {
      if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
    }
    peek <- function() {
      con <- open_text()
      on.exit(close(con), add = TRUE)
      skip <- 0L
      line <- ""
      repeat {
        raw_line <- readLines(con, n = 1, warn = FALSE)
        if (!length(raw_line)) break
        raw_line <- sub("^\ufeff", "", raw_line)
        if (!nzchar(trimws(raw_line))) {
          skip <- skip + 1L
          next
        }
        if (grepl("^#", raw_line) && !grepl("[\t,;]", sub("^#\\s*", "", raw_line))) {
          skip <- skip + 1L
          next
        }
        line <- raw_line
        break
      }
      list(skip = skip, line = line)
    }
    pk <- peek()
    counts <- c(
      tab = length(strsplit(pk$line, "\t", fixed = TRUE)[[1]]),
      comma = length(strsplit(pk$line, ",", fixed = TRUE)[[1]]),
      semi = length(strsplit(pk$line, ";", fixed = TRUE)[[1]])
    )
    sep <- c(tab = "\t", comma = ",", semi = ";")[names(counts)[which.max(counts)]]
    if (max(counts) < 2) {
      sep <- if (identical(ext, "csv")) "," else "\t"
    }
    con <- open_text()
    on.exit(close(con), add = TRUE)
    utils::read.delim(
      con, sep = unname(sep), check.names = FALSE, stringsAsFactors = FALSE,
      quote = "\"", fill = TRUE, comment.char = "",
      skip = pk$skip,
      nrows = if (is.finite(n_max)) n_max else -1
    )
  }
}

is_junk_omics_file <- function(f) {
  b <- basename(f)
  grepl(paste(
    "oncoscan", "probelevel", "probe.?level", "idat", "\\.cel(\\.|$)",
    "fastq", "\\.bam$", "\\.vcf", "copynumber", "copy.?number",
    "\\bcna\\b", "\\bcnv\\b", "snp.?array", "GSE118527"
  , collapse = "|"), b, ignore.case = TRUE)
}

score_expr_file <- function(f) {
  b <- basename(f)
  s <- 0
  if (is_junk_omics_file(f)) return(-999)
  if (grepl("count", b, ignore.case = TRUE)) s <- s + 50
  if (grepl("tpm", b, ignore.case = TRUE)) s <- s + 40
  if (grepl("fpkm", b, ignore.case = TRUE)) s <- s + 35
  if (grepl("rna|express|oep|gene", b, ignore.case = TRUE)) s <- s + 20
  if (grepl("protein|proteome|ratio_matrix|tmt", b, ignore.case = TRUE)) s <- s + 25
  if (grepl("matrix", b, ignore.case = TRUE)) s <- s + 15
  if (grepl("\\.(csv|tsv|txt|xlsx)(\\.gz)?$", b, ignore.case = TRUE)) s <- s + 5
  if (grepl("clinic|pheno|annot|supplement|mmc|prognos|follow|survival|预后|随访", b, ignore.case = TRUE) &&
      !grepl("rna|count|tpm|fpkm|ratio_matrix", b, ignore.case = TRUE)) s <- s - 40
  s
}

score_clin_file <- function(f) {
  b <- basename(f)
  if (is_junk_omics_file(f)) return(-999)
  s <- 0
  if (grepl("prognos|follow.?up|surviv|dfs|rfs|dmfs|ddfs|drfs|os_event|预后|随访|生存",
            b, ignore.case = TRUE)) s <- s + 60
  if (grepl("clinic", b, ignore.case = TRUE)) s <- s + 50
  if (grepl("pheno|annot|patient|clinical_data|meta.?data", b, ignore.case = TRUE)) s <- s + 30
  if (grepl("table.?s1|mmc1|supplement", b, ignore.case = TRUE)) s <- s + 25
  if (s > 0 && grepl("\\.(csv|tsv|txt|xlsx|xls)(\\.gz)?$", b, ignore.case = TRUE)) s <- s + 5
  if (grepl("count|tpm|fpkm|fastq|bam|vcf|ratio_matrix", b, ignore.case = TRUE)) s <- s - 40
  s
}

as_expr_matrix <- function(df) {
  if (is.null(df) || ncol(df) < 2) stop("表达表列数不足")
  nms <- as.character(names(df))
  gene_col <- find_col(
    nms,
    c("^gene.names$", "^gene.?name", "^gene.?symbol", "^symbol$", "^hgnc",
      "^gene$", "gene.names", "Gene names")
  )
  if (is.na(gene_col)) gene_col <- nms[1]
  is_num_col <- function(x) {
    if (is.numeric(x)) {
      return(mean(is.finite(x) | is.na(x)) > 0.3)
    }
    ch <- trimws(as.character(x))
    ch[ch %in% c("", "NA", "NaN", "NAN", "Filtered", "filtered", "NULL", "-")] <- NA
    v <- suppressWarnings(as.numeric(gsub(",", "", ch)))
    mean(!is.na(v)) > 0.5
  }
  num_ok <- vapply(df, is_num_col, logical(1))
  num_ok[match(gene_col, nms)] <- FALSE
  keep_cols <- unique(c(gene_col, nms[num_ok]))
  if (length(keep_cols) < 2) stop("表达表列数不足")
  df <- df[, keep_cols, drop = FALSE]
  numish <- vapply(df, is_num_col, logical(1))
  # 若看起来是样本×基因（行很少、列很多基因名），转置
  if (nrow(df) < 80 && ncol(df) > 200 && mean(numish[-1]) < 0.3) {
    rn <- as.character(df[[1]])
    mat <- t(data.matrix(df[, -1, drop = FALSE]))
    colnames(mat) <- rn
    df <- data.frame(gene = rownames(mat), mat, check.names = FALSE, stringsAsFactors = FALSE)
  }
  genes <- as.character(df[[1]])
  genes <- sub(";.*$", "", genes)
  genes <- trimws(genes)
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  genes[is.na(genes) | genes == ""] <- paste0("row", seq_along(genes))
  genes <- sub("\\.[0-9]+$", "", genes)
  if (any(duplicated(genes))) {
    means <- rowMeans(mat, na.rm = TRUE)
    ord <- order(means, decreasing = TRUE, na.last = TRUE)
    mat <- mat[ord, , drop = FALSE]
    genes <- genes[ord]
    keep <- !duplicated(genes)
    mat <- mat[keep, , drop = FALSE]
    genes <- genes[keep]
  }
  rownames(mat) <- genes
  colnames(mat) <- as.character(colnames(mat))
  if (ncol(mat) < 3 || nrow(mat) < 20) {
    stop("表达矩阵太小: ", nrow(mat), " x ", ncol(mat))
  }
  mat
}

norm_key <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[^A-Z0-9]+", "", x)
  x
}

find_col <- function(nms, patterns) {
  nms <- as.character(nms)
  for (pat in patterns) {
    hit <- which(grepl(pat, nms, ignore.case = TRUE))
    if (length(hit) > 0) return(nms[hit[1]])
  }
  NA_character_
}

find_cols <- function(nms, patterns) {
  nms <- as.character(nms)
  unique(unlist(lapply(patterns, function(pat) nms[grepl(pat, nms, ignore.case = TRUE)]),
                use.names = FALSE))
}

is_time_col <- function(nm) {
  grepl("time|month|day|year|年限|月数|天数|duration|interval|follow.?up.?m",
        nm, ignore.case = TRUE) &&
    !grepl("event|status|flag|site|部位|器官", nm, ignore.case = TRUE)
}

yesish <- function(x) {
  u <- toupper(trimws(as.character(x)))
  u %in% c("1", "YES", "Y", "TRUE", "T", "POS", "POSITIVE", "EVENT", "MET",
           "METASTASIS", "PRESENT", "RECURRENCE", "RELAPSE", "DISTANT",
           "PROGRESSION", "PD", "是", "有", "阳性", "复发", "转移", "进展") |
    grepl("LUNG|BONE|PULMO|OSSEO|肺转移|骨转移", u)
}

parse_organs <- function(x) {
  u <- toupper(paste(as.character(x), collapse = ";"))
  u <- gsub("[,|/]|和|及|、|；", ";", u)
  u <- gsub("锁骨上|SUPRACLAV", " ", u)
  u <- gsub("肺门淋巴|HILAR.?LN|HILAR.?LYMPH", " ", u)
  local_only <- grepl("LOCAL|CHEST.?WALL|胸壁|局部|残腔|切口", u) &&
    !grepl("LUNG|PULMO|BONE|OSSEO|肺|骨", u)
  data.frame(
    lung = !local_only && grepl("LUNG|PULMO|PULMONARY|肺", u),
    bone = !local_only && (
      grepl("BONE|OSSEO|SKELET|RIB|SPIN|SKULL|VERTEBR", u) ||
        grepl("骨转移|骨质|肋骨|脊柱|颅骨|髂骨|股骨|胸骨|骨髓", u) ||
        grepl("(^|;|[[:space:]])骨(;|[[:space:]]|$)", u)
    ),
    any_distant = !local_only && grepl(
      "LUNG|PULMO|BONE|OSSEO|LIVER|BRAIN|DISTANT|METASTA|VISCERAL|肺|骨|肝|脑|远处",
      u
    ),
    stringsAsFactors = FALSE
  )[1, ]
}

infer_kind_tissue <- function(nm, extra = "") {
  blob <- toupper(paste(nm, extra, collapse = " "))
  if (grepl("NORMAL|NAT|ADJACENT|癌旁|正常", blob) && !grepl("TUMOR|PRIMARY|MET", blob)) {
    return(c(kind = "normal", tissue = "Normal"))
  }
  if (grepl("\\bMET\\b|METASTA|TTM|转移", blob) && !grepl("PRIMARY|原发|预后|随访|DFS|RFS", blob)) {
    tissue <- "Other"
    if (grepl("LUNG|PULMO|肺", blob)) tissue <- "Lung"
    if (grepl("BONE|OSSEO|RIB|SPIN|SKULL|骨", blob)) tissue <- "Bone"
    if (grepl("LIVER|肝", blob)) tissue <- "Liver"
    if (grepl("BRAIN|脑", blob)) tissue <- "Brain"
    if (grepl("LN|LYMPH|淋巴结", blob)) tissue <- "LN"
    return(c(kind = "met", tissue = tissue))
  }
  if (grepl("PRIMARY|PRIM|TTP|PT-|原发|TUMOR|TUMOUR", blob)) {
    return(c(kind = "primary", tissue = "Breast"))
  }
  c(kind = "primary", tissue = "Breast")
}

looks_like_expr_table <- function(df) {
  if (is.null(df) || ncol(df) < 3 || nrow(df) < 20) return(FALSE)
  numish <- vapply(df, function(x) mean(suppressWarnings(!is.na(as.numeric(as.character(x))))), numeric(1))
  nrow(df) > 500 && ncol(df) > 30 && mean(numish[-1], na.rm = TRUE) > 0.75
}

read_clin_tables <- function(path) {
  ext <- tolower(sub(".*\\.", "", sub("\\.gz$", "", basename(path), ignore.case = TRUE)))
  out <- list()
  if (ext %in% c("xlsx", "xls")) {
    if (!has_pkg("readxl")) {
      log_msg("需要 readxl 读 Excel: ", path)
      return(out)
    }
    sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) character(0))
    if (length(sheets) == 0) return(out)
    for (sh in sheets) {
      df <- tryCatch(
        as.data.frame(readxl::read_excel(path, sheet = sh, guess_max = 8000),
                      stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      if (is.null(df) || ncol(df) < 2 || nrow(df) < 2) next
      if (looks_like_expr_table(df)) {
        log_msg("Skip expression-like sheet: ", basename(path), " / ", sh)
        next
      }
      attr(df, "source") <- paste0(basename(path), "|", sh)
      out[[length(out) + 1]] <- df
    }
  } else {
    df <- tryCatch(read_table_any(path), error = function(e) {
      log_msg("Clinical read failed (", basename(path), "): ", e$message)
      NULL
    })
    if (!is.null(df) && ncol(df) >= 2 && nrow(df) >= 2 && !looks_like_expr_table(df)) {
      attr(df, "source") <- basename(path)
      out[[1]] <- df
    }
  }
  out
}

id_patterns <- c(
  "^patient", "^sample", "sample.?id", "patient.?id", "barcode", "fuscc",
  "case.?id", "标本", "患者", "住院号", "病理号", "^pid$", "^id$"
)

harmonize_clin <- function(df) {
  id <- find_col(names(df), id_patterns)
  if (is.na(id)) return(NULL)
  df$.id_raw <- as.character(df[[id]])
  df$.key <- norm_key(df$.id_raw)
  keep <- !is.na(df$.key) & nzchar(df$.key) & df$.key != "NA"
  df <- df[keep, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df <- df[!duplicated(df$.key), , drop = FALSE]
  df
}

merge_clin_by_key <- function(tabs) {
  tabs <- Filter(Negate(is.null), lapply(tabs, harmonize_clin))
  if (length(tabs) == 0) return(NULL)
  out <- tabs[[1]]
  if (length(tabs) == 1) return(out)
  for (k in seq_along(tabs)[-1]) {
    extra <- setdiff(names(tabs[[k]]), names(out))
    extra <- extra[extra != ".id_raw"]
    if (length(extra) == 0) next
    add <- tabs[[k]][, c(".key", extra), drop = FALSE]
    out <- merge(out, add, by = ".key", all = TRUE, sort = FALSE)
  }
  out$.id_raw <- ifelse(is.na(out$.id_raw) | !nzchar(as.character(out$.id_raw)),
                        out$.key, as.character(out$.id_raw))
  out
}

all_files <- list_data_files()
if (length(all_files) == 0) {
  stop("在 ", project_dir, " 里没有找到数据文件。请把 FUSCC 表达矩阵和临床/预后表放到该目录。")
}
expr_scores <- vapply(all_files, score_expr_file, numeric(1))
clin_scores <- vapply(all_files, score_clin_file, numeric(1))
expr_order <- all_files[order(expr_scores, decreasing = TRUE)]
expr_order <- expr_order[expr_scores[match(expr_order, all_files)] > 0]
clin_cands <- all_files[clin_scores > 0]
clin_cands <- clin_cands[order(clin_scores[match(clin_cands, all_files)], decreasing = TRUE)]

log_msg("Expression candidates: ",
        if (length(expr_order) == 0) "NONE" else paste(basename(expr_order), collapse = "; "))
log_msg("Clinical/prognosis candidates: ",
        if (length(clin_cands) == 0) "NONE" else paste(basename(utils::head(clin_cands, 15)), collapse = "; "))

raw <- NULL
expr_path <- NA_character_
tried <- character(0)
for (f in expr_order) {
  tried <- c(tried, basename(f))
  log_msg("Trying expression file: ", basename(f))
  df <- tryCatch(read_table_any(f), error = function(e) {
    log_msg("  read failed: ", e$message)
    NULL
  })
  if (is.null(df)) next
  log_msg("  raw table: ", nrow(df), " rows x ", ncol(df), " cols")
  mat <- tryCatch(as_expr_matrix(df), error = function(e) {
    log_msg("  not an expression matrix: ", e$message)
    NULL
  })
  if (is.null(mat)) next
  raw <- mat
  expr_path <- f
  break
}
clin_paths <- setdiff(clin_cands, expr_path)

writeLines(
  c(paste("project_dir:", project_dir),
    paste("expression:", expr_path),
    paste("tried_expression:", paste(tried, collapse = " | ")),
    paste("clinical/prognosis:", paste(clin_paths, collapse = " | ")),
    "All files:",
    paste0("  ", basename(all_files), collapse = "\n")),
  file.path(log_dir, "detected_files.txt")
)
if (is.null(raw) || is.na(expr_path)) {
  stop(
    "找不到能读入的表达矩阵。当前目录: ", project_dir, "\n",
    "已忽略 OncoScan/ProbeLevel 拷贝数文件。\n",
    "请放 RNA counts/TPM/FPKM，或 iProX 蛋白 ratio_matrix_original.csv。\n",
    "尝试过: ", paste(tried, collapse = ", "), "\n",
    "临床表应类似 FUSCC_BRCA_panel_4000_clinical_data.tsv。"
  )
}

log_msg("Using expression file: ", basename(expr_path))
log_msg("Expression: ", nrow(raw), " genes x ", ncol(raw), " samples")

clin_tabs <- list()
for (p in clin_paths) {
  clin_tabs <- c(clin_tabs, read_clin_tables(p))
}
clin <- merge_clin_by_key(clin_tabs)
if (!is.null(clin)) {
  log_msg("Merged clinical/prognosis: ", nrow(clin), " patients x ", ncol(clin),
          " cols from ", length(clin_tabs), " table(s); names: ",
          paste(utils::head(names(clin), 40), collapse = ", "))
} else if (length(clin_paths) > 0) {
  log_msg("Clinical files found but no patient/sample ID column could be merged.")
}

id_col <- if (is.null(clin)) NA_character_ else find_col(names(clin), id_patterns)
site_cols <- if (is.null(clin)) character(0) else find_cols(
  names(clin),
  c("metastat.*site", "relapse.?site", "recur.*site", "dm.?site", "first.?site",
    "organ", "转移部位", "转移灶", "复发部位", "首发", "site of", "visceral",
    "distant.?site")
)
lung_cols <- if (is.null(clin)) character(0) else find_cols(
  names(clin),
  c("lung.?met", "pulmo", "lung.?event", "lung.?relapse", "肺转移", "肺复发")
)
bone_cols <- if (is.null(clin)) character(0) else find_cols(
  names(clin),
  c("bone.?met", "osseo", "bone.?event", "bone.?relapse", "骨转移", "骨复发")
)
any_met_cols <- if (is.null(clin)) character(0) else find_cols(
  names(clin),
  c("distant.?met", "dmfs.?event", "dmfs.?status", "ddfs.?event", "drfs.?event",
    "rfs.?event", "rfs.?status", "dfs.?event", "dfs.?status", "metastasis",
    "relapse", "recurrence", "dm_event",
    "复发", "远处转移", "转移事件", "预后事件")
)
any_met_cols <- any_met_cols[!vapply(any_met_cols, is_time_col, logical(1))]
type_col <- if (is.null(clin)) NA_character_ else find_col(
  names(clin),
  c("sample.?type", "tissue", "specimen", "取材", "组织类型")
)

# 所有列里只要单元格写了肺/骨，也算预后转移部位（不依赖列名）
organ_text_cols <- character(0)
if (!is.null(clin)) {
  organ_text_cols <- names(clin)[vapply(clin, function(x) {
    u <- paste(as.character(utils::head(x, 2000)), collapse = " ")
    grepl("LUNG|PULMO|BONE|OSSEO|肺转移|骨转移|肺,|骨,", u, ignore.case = TRUE) ||
      grepl("肺|骨转移", u)
  }, logical(1))]
  organ_text_cols <- setdiff(organ_text_cols, c(".key", ".id_raw"))
}

cn <- colnames(raw)
nk <- norm_key(cn)
pheno <- data.frame(
  sample = cn,
  patient = cn,
  kind = "primary",
  tissue = "Breast",
  lung_event = FALSE,
  bone_event = FALSE,
  any_met_event = FALSE,
  met_site_text = "",
  stringsAsFactors = FALSE
)
kt <- t(vapply(cn, function(s) infer_kind_tissue(s), character(2)))
pheno$kind <- kt[, 1]
pheno$tissue <- kt[, 2]

match_clin_idx <- function(sample_keys, clin_keys) {
  idx <- match(sample_keys, clin_keys)
  if (mean(is.na(idx)) <= 0.5) return(idx)
  relax <- function(x) {
    x <- gsub("^(FUSCC|FDR|BRCA|P|PT|SAMPLE)", "", x)
    x <- gsub("^0+", "", x)
    x
  }
  idx2 <- match(relax(sample_keys), relax(clin_keys))
  idx[is.na(idx)] <- idx2[is.na(idx)]
  if (mean(is.na(idx)) <= 0.5) return(idx)
  for (i in seq_along(sample_keys)) {
    if (!is.na(idx[i])) next
    sk <- sample_keys[i]
    if (!nzchar(sk) || nchar(sk) < 4) next
    hit <- which(clin_keys == sk | startsWith(sk, clin_keys) | startsWith(clin_keys, sk))
    if (length(hit) == 1) idx[i] <- hit[1]
  }
  idx
}

if (!is.null(clin) && ".key" %in% names(clin)) {
  idx <- match_clin_idx(nk, clin$.key)
  log_msg("Matched clinical/prognosis rows to samples: ", sum(!is.na(idx)), " / ", length(idx))
  pheno$patient <- ifelse(!is.na(idx), as.character(clin$.id_raw[idx]), pheno$patient)
  if (!is.na(type_col)) {
    extra <- ifelse(!is.na(idx), as.character(clin[[type_col]][idx]), "")
    kt2 <- t(mapply(infer_kind_tissue, pheno$sample, extra, USE.NAMES = FALSE))
    pheno$kind <- kt2[, 1]
    pheno$tissue <- kt2[, 2]
  }
  take <- function(cols) {
    if (length(cols) == 0) {
      return(rep("", length(idx)))
    }
    bits <- lapply(cols, function(cl) ifelse(!is.na(idx), as.character(clin[[cl]][idx]), ""))
    do.call(paste, c(bits, sep = ";"))
  }
  site_txt <- take(unique(c(site_cols, organ_text_cols)))
  pheno$met_site_text <- site_txt
  org <- lapply(site_txt, parse_organs)
  pheno$lung_event <- vapply(org, function(z) isTRUE(z$lung), logical(1))
  pheno$bone_event <- vapply(org, function(z) isTRUE(z$bone), logical(1))
  pheno$any_met_event <- vapply(org, function(z) isTRUE(z$any_distant), logical(1))
  for (cl in lung_cols) {
    pheno$lung_event <- pheno$lung_event | (!is.na(idx) & yesish(clin[[cl]][idx]))
  }
  for (cl in bone_cols) {
    pheno$bone_event <- pheno$bone_event | (!is.na(idx) & yesish(clin[[cl]][idx]))
  }
  for (cl in any_met_cols) {
    pheno$any_met_event <- pheno$any_met_event | (!is.na(idx) & yesish(clin[[cl]][idx]))
  }
}

pheno$any_met_event <- pheno$any_met_event | pheno$lung_event | pheno$bone_event | (pheno$kind == "met")
pheno$lung_event <- pheno$lung_event | (pheno$kind == "met" & pheno$tissue == "Lung")
pheno$bone_event <- pheno$bone_event | (pheno$kind == "met" & pheno$tissue == "Bone")
rownames(pheno) <- pheno$sample

col_used <- data.frame(
  role = c(
    rep("id", length(na.omit(id_col))),
    rep("met_site", length(site_cols)),
    rep("lung_flag", length(lung_cols)),
    rep("bone_flag", length(bone_cols)),
    rep("prognosis_event", length(any_met_cols)),
    rep("organ_text_any_col", length(organ_text_cols))
  ),
  column = c(
    if (is.na(id_col)) character(0) else id_col,
    site_cols, lung_cols, bone_cols, any_met_cols, organ_text_cols
  ),
  stringsAsFactors = FALSE
)
if (nrow(col_used) > 0) {
  utils::write.csv(col_used, file.path(log_dir, "prognosis_metastasis_columns.csv"),
                   row.names = FALSE)
}
log_msg("Prognosis lung events: ", sum(pheno$lung_event),
        " bone: ", sum(pheno$bone_event),
        " any_met: ", sum(pheno$any_met_event))

# -----------------------------------------------------------------------------
# 3. 过滤 + 标准化（四个问题共用一次）
# -----------------------------------------------------------------------------
keep_samp <- pheno$kind %in% c("primary", "met")
count_use <- raw[, keep_samp, drop = FALSE]
pheno_use <- pheno[keep_samp, , drop = FALSE]
mx <- stats::quantile(as.numeric(count_use), 0.99, na.rm = TRUE)
frac_int <- mean(abs(count_use - round(count_use)) < 1e-6, na.rm = TRUE)
use_deseq <- is.finite(mx) && mx > 50 && frac_int > 0.7 && has_pkg("DESeq2")
if (use_deseq) {
  log_msg("Treat as counts (q99=", signif(mx, 3), "); DESeq2 size factors")
  keep_gene <- rowSums(count_use, na.rm = TRUE) >= 10 & rowSums(count_use > 0, na.rm = TRUE) >= 3
  count_use <- count_use[keep_gene, , drop = FALSE]
  count_int <- round(pmax(count_use, 0))
  count_int[is.na(count_int)] <- 0
  dds <- DESeq2::DESeqDataSetFromMatrix(count_int, colData = pheno_use, design = ~ 1)
  dds <- DESeq2::estimateSizeFactors(dds)
  logmat <- log2(DESeq2::counts(dds, normalized = TRUE) + 1)
} else {
  log_msg("Treat as normalized / log-scale (q99=", signif(mx, 3), ")")
  keep_gene <- rowSums(is.finite(count_use)) >= 3 &
    matrixStats::rowSds(count_use, na.rm = TRUE) > 1e-6
  count_use <- count_use[keep_gene, , drop = FALSE]
  if (is.finite(mx) && mx > 20) {
    logmat <- log2(pmax(count_use, 0) + 1)
    if (has_pkg("limma")) logmat <- limma::normalizeQuantiles(logmat)
  } else {
    logmat <- count_use
  }
}
log_msg("After filter: ", nrow(logmat), " genes x ", ncol(logmat), " samples")
heat_mat <- t(scale(t(logmat)))
heat_mat[!is.finite(heat_mat)] <- 0
pheno_use <- pheno_use[colnames(logmat), , drop = FALSE]
utils::write.csv(pheno_use, file.path(result_dir, "00_logs", "sample_annotation.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. 一一对应
# -----------------------------------------------------------------------------
pick_primary <- function(rows) {
  if (nrow(rows) == 1) return(rows$sample[1])
  rows$sample[1]
}

build_pairs <- function(organ, organ_labels) {
  mets <- pheno_use[pheno_use$kind == "met" & pheno_use$tissue %in% organ_labels, , drop = FALSE]
  prim <- pheno_use[pheno_use$kind == "primary", , drop = FALSE]
  out <- list()
  for (pid in intersect(unique(mets$patient), unique(prim$patient))) {
    p_rows <- prim[prim$patient == pid, , drop = FALSE]
    m_rows <- mets[mets$patient == pid, , drop = FALSE]
    p_id <- pick_primary(p_rows)
    m_ids <- m_rows$sample
    out[[length(out) + 1]] <- data.frame(
      pair_id = pid, patient = pid, organ = organ,
      primary = p_id,
      metastasis = if (length(m_ids) == 1) m_ids else paste(m_ids, collapse = ";"),
      n_met_samples = length(m_ids),
      met_tissue = paste(unique(m_rows$tissue), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) {
    return(data.frame(
      pair_id = character(), patient = character(), organ = character(),
      primary = character(), metastasis = character(), n_met_samples = integer(),
      met_tissue = character(), stringsAsFactors = FALSE
    ))
  }
  dplyr::bind_rows(out)
}

empty_pairs <- function() {
  data.frame(
    pair_id = character(), patient = character(), organ = character(),
    primary = character(), metastasis = character(), n_met_samples = integer(),
    met_tissue = character(), stringsAsFactors = FALSE
  )
}

outcome_pairs <- function(organ, event_flag) {
  prim <- pheno_use[pheno_use$kind == "primary" & event_flag, , drop = FALSE]
  if (nrow(prim) == 0) return(empty_pairs())
  prim <- prim[!duplicated(prim$patient), ]
  data.frame(
    pair_id = prim$patient, patient = prim$patient, organ = organ,
    primary = prim$sample, metastasis = NA_character_,
    n_met_samples = 0L, met_tissue = "clinical_outcome",
    stringsAsFactors = FALSE
  )
}

pairs_lung_tissue <- build_pairs("Lung", "Lung")
pairs_bone_tissue <- build_pairs("Bone", c("Bone", "Rib", "Spine", "Spinal", "Skull"))
analysis_mode <- if (nrow(pairs_lung_tissue) >= 2 || nrow(pairs_bone_tissue) >= 2) {
  "PAIRED_TISSUE"
} else {
  "CLINICAL_OUTCOME"
}
pair_dir <- file.path(result_dir, "00_sample_pairing")
dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)

if (identical(analysis_mode, "PAIRED_TISSUE")) {
  pairs_lung <- pairs_lung_tissue
  pairs_bone <- pairs_bone_tissue
  writeLines(
    c("MODE: PAIRED_TISSUE",
      "FC = matched metastasis / matched primary (same patient)."),
    file.path(pair_dir, "00_MODE.txt")
  )
} else {
  lung_flag <- pheno_use$kind == "primary" & pheno_use$lung_event
  bone_flag <- pheno_use$kind == "primary" & pheno_use$bone_event
  pairs_lung <- outcome_pairs("Lung", lung_flag)
  pairs_bone <- outcome_pairs("Bone", bone_flag)
  writeLines(
    c("MODE: CLINICAL_OUTCOME",
      "FUSCC 表达矩阵是原发灶；转移来自同一患者的预后/随访（Patient N 原发 ↔ Patient N 肺/骨转移结局）。",
      "DFS/RFS/DMFS 事件、复发部位、肺转移、骨转移列都会读。",
      "FC = (无该器官转移) / (有该器官转移)；选出后来发生该器官转移的原发里更低的基因。",
      paste("lung_event n=", sum(lung_flag), " bone_event n=", sum(bone_flag),
            " any_met n=", sum(pheno_use$kind == "primary" & pheno_use$any_met_event)),
      paste("id_col=", id_col),
      paste("site_cols=", paste(site_cols, collapse = ",")),
      paste("lung_cols=", paste(lung_cols, collapse = ",")),
      paste("bone_cols=", paste(bone_cols, collapse = ",")),
      paste("prognosis_event_cols=", paste(any_met_cols, collapse = ",")),
      paste("organ_text_cols=", paste(utils::head(organ_text_cols, 20), collapse = ","))),
    file.path(pair_dir, "00_MODE.txt")
  )
  if (sum(lung_flag) < 2 && sum(bone_flag) < 2) {
    log_msg("WARNING: 预后表对上了患者，但肺/骨转移人数不足。",
            "请确认随访表里有转移部位或肺转移/骨转移列（不是只有 OS）。")
  }
}

utils::write.csv(pairs_lung, file.path(pair_dir, "pairs_primary_vs_lung_1to1.csv"), row.names = FALSE)
utils::write.csv(pairs_bone, file.path(pair_dir, "pairs_primary_vs_bone_1to1.csv"), row.names = FALSE)
log_msg("MODE: ", analysis_mode)
log_msg("Lung n=", nrow(pairs_lung), if (nrow(pairs_lung) > 0) paste(":", paste(utils::head(pairs_lung$patient, 20), collapse = ", ")) else "")
log_msg("Bone n=", nrow(pairs_bone), if (nrow(pairs_bone) > 0) paste(":", paste(utils::head(pairs_bone$patient, 20), collapse = ", ")) else "")

# 每个配对患者：原发向量 vs 该器官转移（多灶则均值）
pair_vectors <- function(pairs) {
  prim_mat <- matrix(NA_real_, nrow = nrow(logmat), ncol = nrow(pairs),
                     dimnames = list(rownames(logmat), pairs$patient))
  met_mat <- prim_mat
  for (i in seq_len(nrow(pairs))) {
    p <- pairs$primary[i]
    m_raw <- pairs$metastasis[i]
    m <- if (is.na(m_raw) || !nzchar(m_raw)) character(0) else strsplit(m_raw, ";", fixed = TRUE)[[1]]
    m <- intersect(m, colnames(logmat))
    prim_mat[, i] <- logmat[, p]
    met_mat[, i] <- if (length(m) == 0) {
      NA_real_
    } else if (length(m) == 1) {
      logmat[, m]
    } else {
      rowMeans(logmat[, m, drop = FALSE])
    }
  }
  list(primary = prim_mat, met = met_mat, patients = pairs$patient)
}

# -----------------------------------------------------------------------------
# 5. 配对 limma：转移 vs 原发；FC = 转移/原发（>1 表示原发更低）
# -----------------------------------------------------------------------------
paired_limma <- function(pairs, label) {
  if (nrow(pairs) < 2) {
    log_msg(label, ": fewer than 2 pairs, skip limma")
    return(NULL)
  }
  if (nrow(pairs) < 3) {
    log_msg("WARNING ", label, ": only ", nrow(pairs),
            " pairs. p-values have very few residual df; interpret cautiously.")
  }
  vec <- pair_vectors(pairs)
  expr <- cbind(vec$primary, vec$met)
  colnames(expr) <- c(paste0(vec$patients, "_Primary"), paste0(vec$patients, "_Met"))
  patient <- factor(c(vec$patients, vec$patients))
  tissue  <- factor(c(rep("Primary", ncol(vec$primary)), rep("Met", ncol(vec$met))),
                    levels = c("Primary", "Met"))
  design <- stats::model.matrix(~ patient + tissue)
  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit, coef = "tissueMet", number = Inf, sort.by = "none")
  fc_pairs <- 2^(vec$met - vec$primary)
  colnames(fc_pairs) <- paste0("FC_", vec$patients)
  de <- data.frame(
    gene = rownames(tt),
    log2FC = tt$logFC,
    FC = 2^tt$logFC,
    AveExpr = tt$AveExpr,
    t = tt$t,
    pvalue = tt$P.Value,
    padj = tt$adj.P.Val,
    n_pairs = nrow(pairs),
    n_pairs_FC_gt_1 = rowSums(fc_pairs > 1, na.rm = TRUE),
    n_pairs_FC_gt_1.25 = rowSums(fc_pairs > 1.25, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  de <- cbind(de, as.data.frame(fc_pairs))
  de <- de[order(de$pvalue, -de$log2FC), ]
  rownames(de) <- NULL
  attr(de, "pair_expr") <- vec
  attr(de, "pairs") <- pairs
  attr(de, "label") <- label
  attr(de, "fc_meaning") <- paste(
    "PAIRED_TISSUE: FC = matched metastasis / matched primary (same patient).",
    "Selected genes are LOW in the primary relative to that patient's metastasis."
  )
  attr(de, "log2fc_xlab") <- "log2FC (matched metastasis / primary); >0 = lower in primary"
  de
}

# 临床结局：Patient N 原发 vs Patient N 是否发生该器官转移
# FC = 无该器官转移 / 有该器官转移；>1 表示发生该器官转移的原发里更低
grouped_limma_outcome <- function(event_ids, nonevent_ids, label) {
  event_ids <- unique(intersect(event_ids, colnames(logmat)))
  nonevent_ids <- unique(intersect(nonevent_ids, colnames(logmat)))
  nonevent_ids <- setdiff(nonevent_ids, event_ids)
  if (length(event_ids) < 2 || length(nonevent_ids) < 2) {
    log_msg(label, ": need >=2 event and >=2 non-event primaries (got ",
            length(event_ids), " vs ", length(nonevent_ids), ")")
    return(NULL)
  }
  ids <- c(nonevent_ids, event_ids)
  expr <- logmat[, ids, drop = FALSE]
  grp <- factor(
    c(rep("no_event", length(nonevent_ids)), rep("event", length(event_ids))),
    levels = c("event", "no_event")
  )
  design <- stats::model.matrix(~ grp)
  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit, coef = "grpno_event", number = Inf, sort.by = "none")
  de <- data.frame(
    gene = rownames(tt),
    log2FC = tt$logFC,
    FC = 2^tt$logFC,
    AveExpr = tt$AveExpr,
    t = tt$t,
    pvalue = tt$P.Value,
    padj = tt$adj.P.Val,
    n_event = length(event_ids),
    n_no_event = length(nonevent_ids),
    stringsAsFactors = FALSE
  )
  de <- de[order(de$pvalue, -de$log2FC), ]
  rownames(de) <- NULL
  attr(de, "heatmap_mode") <- "group"
  attr(de, "heat_mat") <- logmat
  attr(de, "heat_samples") <- ids
  attr(de, "heat_groups") <- stats::setNames(as.character(grp), ids)
  attr(de, "fc_meaning") <- paste(
    "CLINICAL_OUTCOME: FC = (no that-organ met) / (has that-organ met).",
    "Patient N primary is paired to Patient N lung/bone outcome.",
    "Selected genes are LOW in primaries that later metastasized to that organ."
  )
  attr(de, "log2fc_xlab") <- "log2FC (no-organ-met / organ-met); >0 = lower in primaries with that-organ met"
  attr(de, "label") <- label
  de
}

select_up_in_met <- function(de, fc_min, p_min = p_cutoff) {
  if (is.null(de) || nrow(de) == 0) {
    return(de[0, , drop = FALSE])
  }
  de[!is.na(de$pvalue) & de$pvalue < p_min & de$FC > fc_min, , drop = FALSE]
}

# -----------------------------------------------------------------------------
# 6. 绘图 / ORA / GSEA
# -----------------------------------------------------------------------------
save_gg <- function(plot, path_stub, width = 8, height = 6) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".pdf"), plot, width = width, height = height),
           error = function(e) log_msg("pdf ggsave failed: ", e$message))
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".png"), plot, width = width, height = height, dpi = 300),
           error = function(e) log_msg("png ggsave failed: ", e$message))
}

note_empty <- function(stub, msg) writeLines(msg, paste0(stub, "_EMPTY.txt"))

try_save_plot <- function(fun, stub, width = 9, height = 7) {
  p <- tryCatch(fun(), error = function(e) {
    log_msg("Plot failed (", basename(stub), "): ", e$message)
    NULL
  })
  if (is.null(p)) return(invisible(FALSE))
  save_gg(p, stub, width = width, height = height)
}

plot_volcano <- function(de, highlight, title, outfile, fc_line = 1) {
  df <- de
  df$y <- -log10(pmax(df$pvalue, 1e-300))
  df$set <- ifelse(df$gene %in% highlight, "selected", "other")
  df$label <- ifelse(df$gene %in% utils::head(highlight, 15), df$gene, NA)
  lfc_line <- log2(fc_line)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = y, color = set)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.4) +
    ggplot2::scale_color_manual(values = c(other = "grey70", selected = "#D62828")) +
    ggplot2::geom_vline(xintercept = c(-lfc_line, lfc_line), linetype = 2, color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = 2, color = "grey40") +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = 30, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 12) +
    xlab <- attr(de, "log2fc_xlab")
    if (is.null(xlab) || !nzchar(xlab)) {
      xlab <- "log2FC; >0 = lower in the index group (primary / event / high PNI)"
    }
    ggplot2::labs(
      title = title,
      x = xlab,
      y = "-log10(p value)", color = NULL
    )
  save_gg(p, outfile)
}

plot_de_bar <- function(sub, title, outfile) {
  if (nrow(sub) == 0) return(invisible(NULL))
  df <- sub[order(sub$log2FC, decreasing = TRUE), , drop = FALSE]
  if (nrow(df) > 60) df <- rbind(utils::head(df, 30), utils::tail(df, 30))
  df$gene <- factor(df$gene, levels = rev(unique(df$gene)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = gene, y = log2FC)) +
    ggplot2::geom_col(fill = "#D62828", width = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title, x = NULL, y = "log2FC (selected; >0 = lower in primary/event/high PNI)")
  save_gg(p, outfile, width = 8, height = max(5, min(16, 0.22 * nrow(df) + 2)))
}

group_heatmap <- function(mat, samples, groups, genes, title, outfile) {
  genes <- intersect(genes, rownames(mat))
  samples <- intersect(samples, colnames(mat))
  if (length(genes) > 200) genes <- genes[seq_len(200)]
  if (length(genes) < 2 || length(samples) < 2) {
    log_msg("Group heatmap skipped: ", title)
    return(invisible(NULL))
  }
  sub <- mat[genes, samples, drop = FALSE]
  ann <- data.frame(Group = groups[match(samples, names(groups))], row.names = samples)
  pal <- list(Group = c(
    low = "#4C78A8", high = "#D62828",
    Primary = "#4C78A8", Metastasis = "#D62828",
    no_event = "#4C78A8", event = "#D62828",
    lung_only = "#4C78A8", bone_only = "#D62828"
  ))
  pal$Group <- pal$Group[names(pal$Group) %in% unique(ann$Group)]
  draw <- function() {
    pheatmap::pheatmap(
      sub, scale = "row", annotation_col = ann, annotation_colors = pal,
      cluster_cols = TRUE, cluster_rows = TRUE,
      show_rownames = nrow(sub) <= 80, fontsize_row = 6, main = title,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100)
    )
  }
  grDevices::pdf(paste0(outfile, ".pdf"), width = 8, height = max(6, min(18, 0.18 * nrow(sub) + 3)))
  on.exit({
    while (grDevices::dev.cur() > 1) grDevices::dev.off()
  }, add = TRUE)
  draw()
  grDevices::dev.off()
  grDevices::png(paste0(outfile, ".png"), width = 2400, height = max(1800, 40 * nrow(sub) + 400), res = 300)
  draw()
  grDevices::dev.off()
}

paired_heatmap <- function(de_obj, genes, title, outfile) {
  vec <- attr(de_obj, "pair_expr")
  genes <- intersect(genes, rownames(vec$primary))
  if (length(genes) > 200) genes <- genes[seq_len(200)]
  if (length(genes) < 2) {
    log_msg("Heatmap skipped (<2 genes): ", title)
    return(invisible(NULL))
  }
  mat <- matrix(NA_real_, nrow = length(genes), ncol = 2 * length(vec$patients),
                dimnames = list(genes, NULL))
  cn <- character(0)
  ann <- data.frame(Patient = character(0), Tissue = character(0), stringsAsFactors = FALSE)
  col_i <- 1
  for (j in seq_along(vec$patients)) {
    mat[, col_i] <- vec$primary[genes, j]
    mat[, col_i + 1] <- vec$met[genes, j]
    cn <- c(cn, paste0(vec$patients[j], "_P"), paste0(vec$patients[j], "_M"))
    ann <- rbind(ann, data.frame(
      Patient = vec$patients[j], Tissue = "Primary", stringsAsFactors = FALSE
    ))
    ann <- rbind(ann, data.frame(
      Patient = vec$patients[j], Tissue = "Metastasis", stringsAsFactors = FALSE
    ))
    col_i <- col_i + 2
  }
  colnames(mat) <- cn
  rownames(ann) <- cn
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  pal <- list(
    Tissue = c(Primary = "#4C78A8", Metastasis = "#D62828")
  )
  draw <- function() {
    pheatmap::pheatmap(
      z, scale = "none", annotation_col = ann, annotation_colors = pal,
      cluster_cols = FALSE, cluster_rows = TRUE,
      show_rownames = nrow(z) <= 80, fontsize_row = 6, main = title,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100)
    )
  }
  grDevices::pdf(paste0(outfile, ".pdf"), width = max(8, 0.45 * ncol(z) + 4),
                 height = max(6, min(18, 0.18 * nrow(z) + 3)))
  on.exit({
    while (grDevices::dev.cur() > 1) grDevices::dev.off()
  }, add = TRUE)
  draw()
  grDevices::dev.off()
  grDevices::png(paste0(outfile, ".png"), width = max(1800, 80 * ncol(z)),
                 height = max(1800, 40 * nrow(z) + 400), res = 300)
  draw()
  grDevices::dev.off()
}

pick_official_symbol <- function(x) {
  x <- trimws(as.character(x))
  if (length(x) != 1 || is.na(x) || x %in% c("", "-", ".", "NA")) return(NA_character_)
  parts <- unlist(strsplit(x, "[,;|/]+"))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) return(NA_character_)
  parts[1]
}

map_to_entrez <- function(symbols) {
  symbols <- unique(vapply(symbols, pick_official_symbol, character(1), USE.NAMES = FALSE))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  if (length(symbols) == 0) return(data.frame(gene = character(), entrez = character()))
  m <- tryCatch(
    clusterProfiler::bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
  )
  if (nrow(m) == 0) return(data.frame(gene = character(), entrez = character()))
  m <- m[!duplicated(m[[1]]), ]
  data.frame(gene = m[[1]], entrez = m[[2]], stringsAsFactors = FALSE)
}

ranked_entrez <- function(de) {
  mp <- map_to_entrez(de$gene)
  de2 <- merge(de, mp, by = "gene")
  de2 <- de2[!is.na(de2$entrez) & !is.na(de2$log2FC), ]
  de2 <- de2[order(abs(de2$log2FC), decreasing = TRUE), ]
  de2 <- de2[!duplicated(de2$entrez), ]
  stats <- de2$log2FC
  names(stats) <- de2$entrez
  sort(stats, decreasing = TRUE)
}

enrich_or_relax <- function(strict_fun, relax_fun, label) {
  obj <- tryCatch(strict_fun(), error = function(e) {
    log_msg(label, " strict failed: ", e$message)
    NULL
  })
  if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
    attr(obj, "relaxed") <- FALSE
    return(obj)
  }
  obj2 <- tryCatch(relax_fun(), error = function(e) {
    log_msg(label, " relaxed failed: ", e$message)
    NULL
  })
  if (!is.null(obj2)) attr(obj2, "relaxed") <- TRUE
  obj2
}

plot_ora_object <- function(x, stub, title, fold_change = NULL) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) {
    note_empty(stub, "no enrichment terms")
    return(invisible(NULL))
  }
  df <- as.data.frame(x)
  utils::write.csv(df, paste0(stub, ".csv"), row.names = FALSE)
  nshow <- min(15, nrow(df))
  try_save_plot(function() enrichplot::dotplot(x, showCategory = nshow) + ggplot2::ggtitle(title),
                paste0(stub, "_dotplot"), 9, 7)
  try_save_plot(function() enrichplot::barplot(x, showCategory = nshow) + ggplot2::ggtitle(title),
                paste0(stub, "_barplot"), 9, 7)
}

run_ora_plots <- function(genes, de_sub, outdir, label, tag) {
  go_dir <- file.path(outdir, "GO")
  pw_dir <- file.path(outdir, "Pathway")
  kg_dir <- file.path(outdir, "KEGG")
  dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(kg_dir, recursive = TRUE, showWarnings = FALSE)
  pref <- paste0(tag, "_")
  mp <- map_to_entrez(genes)
  entrez <- unique(mp$entrez)
  fc_sym <- setNames(de_sub$log2FC, de_sub$gene)
  if (length(entrez) < 3) {
    writeLines(paste("mapped_entrez", length(entrez)), file.path(outdir, paste0(pref, "ORA_skipped.txt")))
    note_empty(file.path(go_dir, paste0(pref, "ORA_GO")), "too few mapped genes")
    note_empty(file.path(pw_dir, paste0(pref, "ORA_Pathway")), "too few mapped genes")
    note_empty(file.path(kg_dir, paste0(pref, "ORA_KEGG")), "too few mapped genes")
    return(invisible(NULL))
  }
  for (ont in c("BP", "MF", "CC")) {
    ego <- enrich_or_relax(
      function() clusterProfiler::enrichGO(
        gene = entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
        pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
      ),
      function() clusterProfiler::enrichGO(
        gene = entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
        pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
      ),
      paste("enrichGO", ont)
    )
    plot_ora_object(ego, file.path(go_dir, paste0(pref, "ORA_GO_", ont)),
                    paste(label, "| ORA GO", ont), fold_change = fc_sym)
  }
  ek <- enrich_or_relax(
    function() clusterProfiler::enrichKEGG(
      gene = entrez, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.2
    ),
    function() clusterProfiler::enrichKEGG(
      gene = entrez, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1
    ),
    "enrichKEGG"
  )
  if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) {
    ek <- tryCatch(clusterProfiler::setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
                   error = function(e) ek)
  }
  plot_ora_object(ek, file.path(kg_dir, paste0(pref, "ORA_KEGG")),
                  paste(label, "| ORA KEGG"), fold_change = fc_sym)
  if (has_pkg("ReactomePA")) {
    er <- enrich_or_relax(
      function() ReactomePA::enrichPathway(
        gene = entrez, organism = "human", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
      ),
      function() ReactomePA::enrichPathway(
        gene = entrez, organism = "human", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
      ),
      "enrichPathway"
    )
    plot_ora_object(er, file.path(pw_dir, paste0(pref, "ORA_Reactome_pathway")),
                    paste(label, "| ORA Reactome"), fold_change = fc_sym)
  } else {
    note_empty(file.path(pw_dir, paste0(pref, "ORA_Reactome_pathway")), "ReactomePA not installed")
  }
  writeLines(
    c("GO/Pathway/KEGG 是 ORA，不是 GSEA。",
      "GSEA 在同级 GSEA/ 目录，文件名以 GSEA_ 开头。"),
    file.path(outdir, paste0(pref, "00_ORA_is_not_GSEA.txt"))
  )
}

run_gsea_plots <- function(full_de, sub, outdir, tag, label) {
  gsea_dir <- file.path(outdir, "GSEA")
  dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
  pref <- paste0(tag, "_")
  stats <- ranked_entrez(full_de)
  if (length(stats) < 8) {
    note_empty(file.path(gsea_dir, paste0(pref, "GSEA")), "too few ranked genes")
    return(invisible(NULL))
  }
  kegg <- enrich_or_relax(
    function() clusterProfiler::gseKEGG(
      geneList = stats, organism = "hsa", minGSSize = 5, maxGSSize = 500,
      pvalueCutoff = 0.05, verbose = FALSE, eps = 0
    ),
    function() clusterProfiler::gseKEGG(
      geneList = stats, organism = "hsa", minGSSize = 3, maxGSSize = 500,
      pvalueCutoff = 1, verbose = FALSE, eps = 0
    ),
    paste("gseKEGG", tag)
  )
  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
    kegg <- tryCatch(clusterProfiler::setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
                     error = function(e) kegg)
  }
  plot_ora_object(kegg, file.path(gsea_dir, paste0(pref, "GSEA_KEGG")),
                  paste(label, "| GSEA KEGG"))
  if (has_pkg("msigdbr") && has_pkg("fgsea")) {
    hm <- tryCatch({
      m <- tryCatch(
        msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
        error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H")
      )
      gs <- split(as.character(m$entrez_gene), m$gs_name)
      fgsea::fgsea(pathways = gs, stats = stats, minSize = 5, maxSize = 500)
    }, error = function(e) {
      log_msg("fgsea hallmark failed: ", e$message)
      NULL
    })
    if (!is.null(hm) && nrow(hm) > 0) {
      hm <- hm[order(hm$pval), ]
      utils::write.csv(as.data.frame(hm), file.path(gsea_dir, paste0(pref, "GSEA_Hallmark.csv")),
                       row.names = FALSE)
    }
  }
}

emit_subset <- function(comp_name, de_full, sub, tag, title, outdir, fc_line) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  fc_note <- attr(de_full, "fc_meaning")
  if (is.null(fc_note) || !nzchar(fc_note)) {
    if (identical(attr(de_full, "heatmap_mode"), "group")) {
      fc_note <- "FC direction as in the comparison title; selected genes are LOW in the high-risk / event / high-PNI group."
    } else {
      fc_note <- paste(
        "FC = matched metastasis / matched primary (same patient).",
        "Selected genes are LOW in the primary relative to the paired metastasis."
      )
    }
  }
  writeLines(
    c(paste("comparison:", comp_name),
      paste("subset:", tag),
      paste("title:", title),
      paste("n_genes:", if (is.null(sub)) 0 else nrow(sub)),
      paste("mode:", analysis_mode),
      fc_note),
    file.path(outdir, paste0("00_", tag, "_THIS_FOLDER.txt"))
  )
  if (is.null(sub) || nrow(sub) == 0) {
    writeLines("no genes", file.path(outdir, paste0(tag, "_EMPTY.txt")))
    log_msg(comp_name, " ", tag, ": EMPTY")
    return(invisible(NULL))
  }
  utils::write.csv(sub, file.path(outdir, paste0(tag, "_DE_selected_genes.csv")), row.names = FALSE)
  tryCatch(writexl::write_xlsx(sub, file.path(outdir, paste0(tag, "_DE_selected_genes.xlsx"))),
           error = function(e) log_msg("xlsx write failed: ", e$message))
  log_msg(comp_name, " ", tag, ": n = ", nrow(sub))
  tryCatch(plot_de_bar(sub, paste(title, "| genes low in primary"), file.path(outdir, paste0(tag, "_DE_log2FC_barplot"))),
           error = function(e) log_msg("barplot failed: ", e$message))
  tryCatch(plot_volcano(de_full, sub$gene, title, file.path(outdir, paste0(tag, "_volcano")), fc_line = fc_line),
           error = function(e) log_msg("volcano failed: ", e$message))
  hm_mode <- attr(de_full, "heatmap_mode")
  if (identical(hm_mode, "group")) {
    tryCatch(group_heatmap(
      attr(de_full, "heat_mat"), attr(de_full, "heat_samples"),
      attr(de_full, "heat_groups"), sub$gene, title,
      file.path(outdir, paste0(tag, "_heatmap_high_vs_low"))
    ), error = function(e) {
      while (grDevices::dev.cur() > 1) grDevices::dev.off()
      log_msg("group heatmap failed: ", e$message)
    })
  } else {
    tryCatch(paired_heatmap(de_full, sub$gene, title, file.path(outdir, paste0(tag, "_heatmap_paired_1to1"))),
             error = function(e) {
               while (grDevices::dev.cur() > 1) grDevices::dev.off()
               log_msg("heatmap failed: ", e$message)
             })
  }
  tryCatch(run_ora_plots(sub$gene, sub, outdir, title, tag),
           error = function(e) log_msg("ORA failed: ", e$message))
  tryCatch(run_gsea_plots(de_full, sub, outdir, tag, title),
           error = function(e) log_msg("GSEA failed: ", e$message))
}

analyze_paired <- function(comp_name, pairs, folder) {
  base <- file.path(result_dir, folder)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  de <- paired_limma(pairs, comp_name)
  if (is.null(de)) {
    writeLines("not enough pairs", file.path(base, "SKIPPED_not_enough_pairs.txt"))
    return(NULL)
  }
  utils::write.csv(de, file.path(base, paste0(comp_name, "_full_genome_paired_DE.csv")), row.names = FALSE)
  tryCatch(writexl::write_xlsx(de, file.path(base, paste0(comp_name, "_full_genome_paired_DE.xlsx"))),
           error = function(e) log_msg(e$message))
  gsea_all <- file.path(base, "00_GSEA_all_genes_NOT_FC_or_topN")
  dir.create(gsea_all, recursive = TRUE, showWarnings = FALSE)
  tryCatch(run_gsea_plots(de, de, gsea_all, "ALLgenes", paste(comp_name, "all genes")),
           error = function(e) log_msg("all-gene GSEA failed: ", e$message))
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    sub <- select_up_in_met(de, fc_min = fc)
    outdir <- file.path(base, "FoldChange", nm)
    title <- paste0(comp_name, " | p<", p_cutoff, " & FC>", fc, " (low in primary)")
    emit_subset(comp_name, de, sub, nm, title, outdir, fc_line = fc)
  }
  de
}

analyze_grouped <- function(comp_name, event_ids, nonevent_ids, folder) {
  base <- file.path(result_dir, folder)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  de <- grouped_limma_outcome(event_ids, nonevent_ids, comp_name)
  if (is.null(de)) {
    writeLines(
      c("not enough event vs non-event primaries",
        paste("n_event_in=", length(unique(event_ids))),
        paste("n_nonevent_in=", length(unique(nonevent_ids)))),
      file.path(base, "SKIPPED_not_enough_groups.txt")
    )
    return(NULL)
  }
  utils::write.csv(de, file.path(base, paste0(comp_name, "_full_genome_clinical_outcome_DE.csv")), row.names = FALSE)
  tryCatch(writexl::write_xlsx(de, file.path(base, paste0(comp_name, "_full_genome_clinical_outcome_DE.xlsx"))),
           error = function(e) log_msg(e$message))
  gsea_all <- file.path(base, "00_GSEA_all_genes_NOT_FC_or_topN")
  dir.create(gsea_all, recursive = TRUE, showWarnings = FALSE)
  tryCatch(run_gsea_plots(de, de, gsea_all, "ALLgenes", paste(comp_name, "all genes")),
           error = function(e) log_msg("all-gene GSEA failed: ", e$message))
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    sub <- select_up_in_met(de, fc_min = fc)
    outdir <- file.path(base, "FoldChange", nm)
    title <- paste0(comp_name, " | p<", p_cutoff, " & FC>", fc,
                    " (low in primaries with that-organ met)")
    emit_subset(comp_name, de, sub, nm, title, outdir, fc_line = fc)
  }
  de
}

primaries_one <- pheno_use[pheno_use$kind == "primary", , drop = FALSE]
primaries_one <- primaries_one[!duplicated(primaries_one$patient), ]
utils::write.csv(
  primaries_one[, intersect(
    c("sample", "patient", "lung_event", "bone_event", "any_met_event", "met_site_text"),
    names(primaries_one)
  )],
  file.path(pair_dir, "patient_primary_to_organ_outcome_1to1.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Q1 肺 / 骨：组织配对或临床结局（Patient N 原发 ↔ Patient N 该器官）
# -----------------------------------------------------------------------------
if (identical(analysis_mode, "PAIRED_TISSUE")) {
  de_lung <- analyze_paired("paired_primary_vs_lung", pairs_lung, "01_lung_paired_met_vs_primary")
  de_bone <- analyze_paired("paired_primary_vs_bone", pairs_bone, "02_bone_paired_met_vs_primary")
} else {
  de_lung <- analyze_grouped(
    "clinical_primary_lung_outcome",
    primaries_one$sample[primaries_one$lung_event],
    primaries_one$sample[!primaries_one$lung_event],
    "01_lung_clinical_outcome_vs_primary"
  )
  de_bone <- analyze_grouped(
    "clinical_primary_bone_outcome",
    primaries_one$sample[primaries_one$bone_event],
    primaries_one$sample[!primaries_one$bone_event],
    "02_bone_clinical_outcome_vs_primary"
  )
}

# -----------------------------------------------------------------------------
# Q2 器官特异：肺对过阈值但骨对不过，或反过来
# -----------------------------------------------------------------------------
organ_dir <- file.path(result_dir, "03_organ_specific")
dir.create(organ_dir, recursive = TRUE, showWarnings = FALSE)

write_setdiff <- function(a, b, name_a, name_b, out_stub, de_ref, folder_tag) {
  genes_a <- if (is.null(a) || nrow(a) == 0) character(0) else a$gene
  genes_b <- if (is.null(b) || nrow(b) == 0) character(0) else b$gene
  only_a <- setdiff(genes_a, genes_b)
  sub <- if (is.null(a) || length(only_a) == 0) a[0, ] else a[a$gene %in% only_a, , drop = FALSE]
  title <- paste0(name_a, " specific vs ", name_b)
  emit_subset(title, de_ref, sub, folder_tag, title, out_stub, fc_line = 1)
  invisible(sub)
}

if (!is.null(de_lung) && !is.null(de_bone)) {
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    lung_sub <- select_up_in_met(de_lung, fc)
    bone_sub <- select_up_in_met(de_bone, fc)
    write_setdiff(
      lung_sub, bone_sub, "Lung", "Bone",
      file.path(organ_dir, nm, "lung_specific_not_bone"),
      de_lung, paste0(nm, "_lung_specific")
    )
    write_setdiff(
      bone_sub, lung_sub, "Bone", "Lung",
      file.path(organ_dir, nm, "bone_specific_not_lung"),
      de_bone, paste0(nm, "_bone_specific")
    )
    both <- intersect(lung_sub$gene, bone_sub$gene)
    overlap <- if (length(both) == 0) lung_sub[0, ] else lung_sub[lung_sub$gene %in% both, ]
    emit_subset(
      "lung_and_bone_shared", de_lung, overlap, paste0(nm, "_shared"),
      paste0("Shared lung & bone | p<", p_cutoff, " FC>", fc),
      file.path(organ_dir, nm, "shared_lung_and_bone"), fc
    )
    if (has_pkg("ggvenn")) {
      vdf <- list(Lung = lung_sub$gene, Bone = bone_sub$gene)
      p <- ggvenn::ggvenn(vdf, fill_color = c("#4C78A8", "#D62828")) +
        ggplot2::ggtitle(paste("Organ-specific DE", nm, "p<", p_cutoff))
      save_gg(p, file.path(organ_dir, nm, paste0(nm, "_venn_lung_vs_bone")))
    }
  }
} else {
  writeLines(
    "Need both lung and bone DE (paired tissue or clinical outcome) to call organ-specific genes.",
    file.path(organ_dir, "SKIPPED.txt")
  )
}

# 原发水平：只转肺 vs 只转骨（CLINICAL_OUTCOME 用随访器官；PAIRED_TISSUE 用配对组织）
lung_only_pid <- setdiff(pairs_lung$patient, pairs_bone$patient)
bone_only_pid <- setdiff(pairs_bone$patient, pairs_lung$patient)
tropism <- data.frame(
  patient = c(lung_only_pid, bone_only_pid),
  tropism = c(rep("lung_only", length(lung_only_pid)),
              rep("bone_only", length(bone_only_pid))),
  stringsAsFactors = FALSE
)
utils::write.csv(tropism, file.path(organ_dir, "patients_lung_only_vs_bone_only.csv"), row.names = FALSE)
log_msg("Lung-only paired patients: ", paste(lung_only_pid, collapse = ", "))
log_msg("Bone-only paired patients: ", paste(bone_only_pid, collapse = ", "))

if (length(lung_only_pid) >= 2 && length(bone_only_pid) >= 2) {
  prim_ids <- c(
    pairs_lung$primary[match(lung_only_pid, pairs_lung$patient)],
    pairs_bone$primary[match(bone_only_pid, pairs_bone$patient)]
  )
  grp <- factor(c(rep("lung_only", length(lung_only_pid)),
                  rep("bone_only", length(bone_only_pid))),
                levels = c("bone_only", "lung_only"))
  design <- stats::model.matrix(~ grp)
  fit <- limma::lmFit(logmat[, prim_ids, drop = FALSE], design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit, coef = "grplung_only", number = Inf)
  trop_de <- data.frame(
    gene = rownames(tt), log2FC = tt$logFC, FC = 2^tt$logFC,
    pvalue = tt$P.Value, padj = tt$adj.P.Val, AveExpr = tt$AveExpr,
    stringsAsFactors = FALSE
  )
  trop_de <- trop_de[order(trop_de$pvalue), ]
  utils::write.csv(trop_de, file.path(organ_dir, "primary_lung_only_vs_bone_only_full_DE.csv"),
                   row.names = FALSE)
  attr(trop_de, "heatmap_mode") <- "group"
  attr(trop_de, "heat_mat") <- logmat
  attr(trop_de, "heat_samples") <- prim_ids
  attr(trop_de, "heat_groups") <- stats::setNames(as.character(grp), prim_ids)
  attr(trop_de, "fc_meaning") <- "FC = lung_only primary / bone_only primary; genes with FC < 1/threshold are lower in lung-only primaries."
  attr(trop_de, "log2fc_xlab") <- "log2FC (lung-only / bone-only primary)"
  # 原发里更低、且偏向只转肺：log2FC(lung_only / bone_only) < 0
  for (nm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[nm]])
    sub <- trop_de[!is.na(trop_de$pvalue) & trop_de$pvalue < p_cutoff & trop_de$FC < (1 / fc), ]
    if (nrow(sub) > 0) sub$FC_low_in_lung_primary <- 1 / sub$FC
    emit_subset(
      "primary_low_in_lung_only", trop_de, sub, nm,
      paste0("Primary lower in lung-only vs bone-only | p<", p_cutoff, " FC>", fc),
      file.path(organ_dir, "primary_unpaired_lung_only_low", nm), fc
    )
  }
}

# -----------------------------------------------------------------------------
# Q3 神经浸润签名（原发）
# 公开表通常无病理 PNI；用轴突导向 / 施旺 / 神经营养 三套基因在原发打分
# -----------------------------------------------------------------------------
pni_curated <- list(
  axon_guidance = c(
    "SEMA3A", "SEMA3B", "SEMA3C", "SEMA3D", "SEMA3E", "SEMA3F", "SEMA3G",
    "SEMA4A", "SEMA4B", "SEMA4C", "SEMA4D", "SEMA4F", "SEMA4G",
    "SEMA5A", "SEMA5B", "SEMA6A", "SEMA6B", "SEMA6C", "SEMA6D", "SEMA7A",
    "PLXNA1", "PLXNA2", "PLXNA3", "PLXNA4", "PLXNB1", "PLXNB2", "PLXNB3",
    "PLXNC1", "PLXND1", "NRP1", "NRP2",
    "EPHA1", "EPHA2", "EPHA3", "EPHA4", "EPHA5", "EPHA6", "EPHA7", "EPHA8",
    "EPHB1", "EPHB2", "EPHB3", "EPHB4", "EPHB6",
    "EFNA1", "EFNA2", "EFNA3", "EFNA4", "EFNA5", "EFNB1", "EFNB2", "EFNB3",
    "NTN1", "NTN4", "DCC", "UNC5A", "UNC5B", "UNC5C", "UNC5D",
    "SLIT1", "SLIT2", "SLIT3", "ROBO1", "ROBO2", "ROBO3", "ROBO4",
    "NTNG1", "NTNG2", "RHOA", "RAC1", "CDC42", "PAK1", "GSK3B", "MAPK1"
  ),
  schwann = c(
    "SOX10", "S100B", "MPZ", "MBP", "PMP22", "MAG", "PRX", "EGR2", "POU3F1",
    "ERBB2", "ERBB3", "NRG1", "NCAM1", "L1CAM", "NGFR", "GFAP", "GAP43",
    "MAL", "GJB1", "PMP2", "DRP2", "PLP1", "MPZL1", "CDH19", "FOXD3",
    "SOX2", "OCT6", "POU3F2", "ITGA4", "ITGB8", "LAMA2", "LAMB2"
  ),
  neurotrophic = c(
    "NGF", "BDNF", "NTF3", "NTF4", "NTRK1", "NTRK2", "NTRK3", "NGFR",
    "GDNF", "NRTN", "ARTN", "PSPN", "GFRA1", "GFRA2", "GFRA3", "GFRA4",
    "RET", "CNTF", "CNTFR", "LIF", "LIFR", "OSM", "OSMR", "IGF1", "IGF1R",
    "VEGFA", "NTRK2", "SORT1", "NTRK1"
  )
)

expand_msig <- function(key_words) {
  if (!has_pkg("msigdbr")) return(character(0))
  m <- tryCatch({
    tryCatch(
      msigdbr::msigdbr(species = "Homo sapiens", collection = "C5"),
      error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "C5")
    )
  }, error = function(e) NULL)
  if (is.null(m) || !("gs_name" %in% names(m))) return(character(0))
  hit <- m[grepl(key_words, m$gs_name, ignore.case = TRUE), ]
  unique(as.character(hit$gene_symbol))
}

pni_sets <- list(
  axon_guidance = unique(c(pni_curated$axon_guidance, expand_msig("AXON_GUIDANCE"))),
  schwann = unique(c(pni_curated$schwann, expand_msig("SCHWANN"))),
  neurotrophic = unique(c(
    pni_curated$neurotrophic,
    expand_msig("NEUROTROPH")
  ))
)

zscore_sig <- function(mat, genes) {
  genes <- intersect(unique(genes), rownames(mat))
  if (length(genes) < 3) {
    return(list(score = setNames(rep(NA_real_, ncol(mat)), colnames(mat)), genes = genes))
  }
  z <- t(scale(t(mat[genes, , drop = FALSE])))
  z[!is.finite(z)] <- NA
  list(score = colMeans(z, na.rm = TRUE), genes = genes)
}

prim_samples <- pheno_use$sample[pheno_use$kind == "primary"]
pni_dir <- file.path(result_dir, "04_PNI_primary")
dir.create(pni_dir, recursive = TRUE, showWarnings = FALSE)

score_df <- data.frame(
  sample = prim_samples,
  patient = pheno_use$patient[match(prim_samples, pheno_use$sample)],
  has_paired_lung = pheno_use$patient[match(prim_samples, pheno_use$sample)] %in% pairs_lung$patient,
  has_paired_bone = pheno_use$patient[match(prim_samples, pheno_use$sample)] %in% pairs_bone$patient,
  stringsAsFactors = FALSE
)
# CLINICAL_OUTCOME 时 has_paired_lung = 该患者后来发生肺转移（不是转移组织）
used_genes <- list()
for (nm in names(pni_sets)) {
  sc <- zscore_sig(logmat[, prim_samples, drop = FALSE], pni_sets[[nm]])
  score_df[[paste0("score_", nm)]] <- sc$score[prim_samples]
  used_genes[[nm]] <- sc$genes
  writeLines(sc$genes, file.path(pni_dir, paste0("signature_genes_", nm, ".txt")))
  log_msg("PNI signature ", nm, ": ", length(sc$genes), " genes in matrix")
}
utils::write.csv(score_df, file.path(pni_dir, "primary_PNI_signature_scores.csv"), row.names = FALSE)

unpaired_limma <- function(expr, group, coef_name) {
  group <- factor(group)
  design <- stats::model.matrix(~ group)
  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  coef <- grep(coef_name, colnames(design), value = TRUE)
  if (length(coef) == 0) coef <- colnames(design)[ncol(design)]
  tt <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
  data.frame(
    gene = rownames(tt), log2FC = tt$logFC, FC = 2^tt$logFC,
    AveExpr = tt$AveExpr, t = tt$t, pvalue = tt$P.Value, padj = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
}

# 高神经浸润 vs 低：基因在高分组更低 => log2FC(high/low) < 0，FC_low/high > 阈值
pni_de_list <- list()
for (nm in names(pni_sets)) {
  sc <- score_df[[paste0("score_", nm)]]
  ok <- prim_samples[!is.na(sc)]
  sc2 <- sc[match(ok, prim_samples)]
  hi <- ok[sc2 >= stats::median(sc2, na.rm = TRUE)]
  lo <- ok[sc2 < stats::median(sc2, na.rm = TRUE)]
  log_msg("PNI ", nm, " high n=", length(hi), " low n=", length(lo))
  expr <- logmat[, c(lo, hi), drop = FALSE]
  grp <- factor(c(rep("low", length(lo)), rep("high", length(hi))), levels = c("low", "high"))
  de <- unpaired_limma(expr, grp, "high")
  de <- de[order(de$pvalue), ]
  # 高 PNI 原发里更低：FC_high/low < 1
  de$FC_low_over_high <- 1 / de$FC
  pni_de_list[[nm]] <- de
  base <- file.path(pni_dir, nm)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(de, file.path(base, paste0(nm, "_high_vs_low_full_DE.csv")), row.names = FALSE)
  grp_ids <- c(lo, hi)
  grp_lab <- setNames(c(rep("low", length(lo)), rep("high", length(hi))), grp_ids)
  for (fcnm in names(fc_cutoffs)) {
    fc <- unname(fc_cutoffs[[fcnm]])
    sub <- de[!is.na(de$pvalue) & de$pvalue < p_cutoff & de$FC_low_over_high > fc, ]
    if (nrow(sub) > 0) {
      sub$log2FC <- log2(sub$FC_low_over_high)
      sub$FC <- sub$FC_low_over_high
    }
    full_for_volcano <- de
    full_for_volcano$log2FC <- log2(pmax(full_for_volcano$FC_low_over_high, 1e-8))
    full_for_volcano$FC <- full_for_volcano$FC_low_over_high
    attr(full_for_volcano, "heatmap_mode") <- "group"
    attr(full_for_volcano, "heat_mat") <- logmat
    attr(full_for_volcano, "heat_samples") <- grp_ids
    attr(full_for_volcano, "heat_groups") <- grp_lab
    attr(full_for_volcano, "fc_meaning") <- paste(
      "PNI high vs low (median split on primary signature).",
      "FC = low-PNI / high-PNI; selected genes are LOW in high-score primaries."
    )
    attr(full_for_volcano, "log2fc_xlab") <- "log2FC (low-PNI / high-PNI); >0 = lower in high-PNI primaries"
    emit_subset(
      paste0("PNI_", nm), full_for_volcano, sub, fcnm,
      paste0("PNI ", nm, " | genes lower in high-score primaries | p<", p_cutoff, " FC>", fc),
      file.path(base, "FoldChange", fcnm), fc
    )
  }
}

# -----------------------------------------------------------------------------
# Q4 三种神经浸润 vs 肺转移（分别做）
# -----------------------------------------------------------------------------
q4_dir <- file.path(result_dir, "05_PNI_vs_lung_met")
dir.create(q4_dir, recursive = TRUE, showWarnings = FALSE)

score_long <- tidyr::pivot_longer(
  score_df,
  cols = dplyr::starts_with("score_"),
  names_to = "signature",
  values_to = "score"
)
score_long$signature <- sub("^score_", "", score_long$signature)
if (identical(analysis_mode, "PAIRED_TISSUE")) {
  score_long$lung_status <- ifelse(score_long$has_paired_lung,
                                   "primary_with_paired_lung",
                                   "primary_without_paired_lung")
  q4_box_title <- "Primary PNI signature scores vs paired lung metastasis"
} else {
  score_long$lung_status <- ifelse(score_long$has_paired_lung,
                                   "primary_with_lung_met_outcome",
                                   "primary_without_lung_met_outcome")
  q4_box_title <- "Primary PNI signature scores vs Patient N lung-met outcome"
}

p <- ggplot2::ggplot(score_long, ggplot2::aes(x = lung_status, y = score, fill = signature)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  ggplot2::geom_jitter(width = 0.12, size = 1.4, alpha = 0.8) +
  ggplot2::facet_wrap(~ signature, scales = "free_y") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)) +
  ggplot2::labs(
    title = q4_box_title,
    x = NULL, y = "z-score (primary tumor)"
  )
save_gg(p, file.path(q4_dir, "PNI_scores_by_lung_pair_status"), width = 11, height = 5)

wilcox_rows <- list()
for (nm in names(pni_sets)) {
  sub <- score_df[!is.na(score_df[[paste0("score_", nm)]]), ]
  a <- sub[[paste0("score_", nm)]][sub$has_paired_lung]
  b <- sub[[paste0("score_", nm)]][!sub$has_paired_lung]
  wt <- tryCatch(stats::wilcox.test(a, b, exact = FALSE), error = function(e) NULL)
  wilcox_rows[[nm]] <- data.frame(
    signature = nm,
    n_primary_with_lung = length(a),
    n_primary_without_lung = length(b),
    median_with_lung = stats::median(a, na.rm = TRUE),
    median_without_lung = stats::median(b, na.rm = TRUE),
    comparison = if (identical(analysis_mode, "PAIRED_TISSUE")) {
      "primary_with_vs_without_paired_lung_tissue"
    } else {
      "primary_with_vs_without_lung_met_outcome"
    },
    wilcox_p = if (is.null(wt)) NA_real_ else unname(wt$p.value),
    stringsAsFactors = FALSE
  )
}
wilcox_tab <- dplyr::bind_rows(wilcox_rows)
utils::write.csv(wilcox_tab, file.path(q4_dir, "PNI_score_wilcox_paired_lung_vs_not.csv"), row.names = FALSE)

# 配对肺组织：原发 vs 肺转移 的三种签名分数（同一患者一一对应）
# CLINICAL_OUTCOME 没有转移组织，跳过这段 Wilcoxon
tissue_pairs_ok <- identical(analysis_mode, "PAIRED_TISSUE") &&
  nrow(pairs_lung) >= 2 &&
  !any(is.na(pairs_lung$metastasis))
if (tissue_pairs_ok) {
  vec <- pair_vectors(pairs_lung)
  paired_score <- list()
  for (nm in names(pni_sets)) {
    sc_p <- zscore_sig(vec$primary, pni_sets[[nm]])$score
    sc_m <- zscore_sig(vec$met, pni_sets[[nm]])$score
    df <- data.frame(
      patient = names(sc_p),
      signature = nm,
      score_primary = unname(sc_p),
      score_lung_met = unname(sc_m),
      delta_met_minus_primary = unname(sc_m - sc_p),
      stringsAsFactors = FALSE
    )
    paired_score[[nm]] <- df
    wt <- tryCatch(stats::wilcox.test(df$score_lung_met, df$score_primary, paired = TRUE, exact = FALSE),
                   error = function(e) NULL)
    log_msg("Paired PNI ", nm, " lung vs primary Wilcoxon p=",
            if (is.null(wt)) "NA" else signif(wt$p.value, 3))
    pdf <- tidyr::pivot_longer(df, cols = c("score_primary", "score_lung_met"),
                               names_to = "tissue", values_to = "score")
    pdf$tissue <- ifelse(pdf$tissue == "score_primary", "Primary", "Lung met")
    pp <- ggplot2::ggplot(pdf, ggplot2::aes(x = tissue, y = score, group = patient)) +
      ggplot2::geom_line(color = "grey60") +
      ggplot2::geom_point(ggplot2::aes(color = tissue), size = 2.5) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(
        title = paste0("1-to-1 PNI ", nm, " : primary vs paired lung"),
        x = NULL, y = "signature z-score"
      )
    save_gg(pp, file.path(q4_dir, nm, paste0(nm, "_paired_primary_vs_lung_score")))
  }
  utils::write.csv(dplyr::bind_rows(paired_score),
                   file.path(q4_dir, "paired_PNI_scores_primary_vs_lung.csv"),
                   row.names = FALSE)
} else {
  writeLines(
    c("SKIPPED: no matched lung-metastasis tissue RNA/protein.",
      paste("MODE:", analysis_mode),
      "Q4 still compares PNI scores in primaries with vs without Patient N lung-met outcome,",
      "and intersects PNI-low genes with Q1 lung-low genes."),
    file.path(q4_dir, "SKIPPED_no_matched_lung_tissue.txt")
  )
}

# 基因交集：原发相对肺转移（或肺结局）更低 ∩ 高神经浸润原发里更低
if (!is.null(de_lung)) {
  for (nm in names(pni_de_list)) {
    for (fcnm in names(fc_cutoffs)) {
      fc <- unname(fc_cutoffs[[fcnm]])
      lung_sub <- select_up_in_met(de_lung, fc)
      pni_de <- pni_de_list[[nm]]
      pni_sub <- pni_de[!is.na(pni_de$pvalue) & pni_de$pvalue < p_cutoff & pni_de$FC_low_over_high > fc, ]
      both <- intersect(lung_sub$gene, pni_sub$gene)
      ov <- if (length(both) == 0) lung_sub[0, ] else lung_sub[lung_sub$gene %in% both, ]
      emit_subset(
        paste0("PNI_", nm, "_AND_lung"),
        de_lung, ov, fcnm,
        paste0(nm, " PNI-low genes also low in primary vs lung (tissue or outcome) | ", fcnm),
        file.path(q4_dir, nm, "overlap_with_lung_paired", fcnm), fc
      )
    }
  }
}

# -----------------------------------------------------------------------------
# 总表
# -----------------------------------------------------------------------------
summarize_n <- function(de, fc) {
  if (is.null(de) || nrow(de) == 0) return(0L)
  nrow(select_up_in_met(de, fc))
}
q1_lung_lab <- if (identical(analysis_mode, "PAIRED_TISSUE")) "Q1_lung_paired" else "Q1_lung_clinical_outcome"
q1_bone_lab <- if (identical(analysis_mode, "PAIRED_TISSUE")) "Q1_bone_paired" else "Q1_bone_clinical_outcome"
summary_tab <- data.frame(
  question = c(
    q1_lung_lab, q1_bone_lab,
    "Q2_lung_specific_FC1.25", "Q2_bone_specific_FC1.25",
    "Q3_axon_guidance", "Q3_schwann", "Q3_neurotrophic"
  ),
  n_pairs_or_primaries = c(
    nrow(pairs_lung), nrow(pairs_bone),
    nrow(pairs_lung), nrow(pairs_bone),
    sum(!is.na(score_df$score_axon_guidance)),
    sum(!is.na(score_df$score_schwann)),
    sum(!is.na(score_df$score_neurotrophic))
  ),
  n_genes_p0.05_FC_gt_1 = c(
    summarize_n(de_lung, 1), summarize_n(de_bone, 1),
    NA, NA, NA, NA, NA
  ),
  n_genes_p0.05_FC_gt_1.25 = c(
    summarize_n(de_lung, 1.25), summarize_n(de_bone, 1.25),
    NA, NA, NA, NA, NA
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(summary_tab, file.path(result_dir, "00_run_summary.csv"), row.names = FALSE)
writeLines(
  c(
    "FUSCC_Human_breast.R finished.",
    paste("MODE:", analysis_mode),
    paste("Lung n:", nrow(pairs_lung), paste(utils::head(pairs_lung$patient, 30), collapse = ", ")),
    paste("Bone n:", nrow(pairs_bone), paste(utils::head(pairs_bone$patient, 30), collapse = ", ")),
    if (identical(analysis_mode, "PAIRED_TISSUE")) {
      "Heatmaps keep patient order: Primary_i next to Met_i (cluster_cols = FALSE)."
    } else {
      "Patient N primary is paired to Patient N prognosis/follow-up lung/bone metastasis (DFS/RFS/DMFS and met site)."
    },
    "No Top50-300. Thresholds: p < 0.05 and FC > 1 / 1.25 only.",
    "PNI is signature-based (axon guidance / Schwann / neurotrophic); pathology PNI is not required."
  ),
  file.path(result_dir, "00_README_results.txt")
)
log_msg("Done. Results in ", result_dir)
