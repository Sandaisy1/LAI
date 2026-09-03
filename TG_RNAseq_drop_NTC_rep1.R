#!/usr/bin/env Rscript
# =============================================================================
# 从 Cuffdiff 原始结果中删除 NTC_rep1，其余样品按原格式另存。
# 不覆盖 E:/R/TG_BRCA/TG 里的原文件。
#
# 输出目录：<数据目录>/without_NTC_rep1/
# 保留：NTC_rep0、TG_sh1、TG_sh5
#
# 用法：
#   setwd("E:/R/TG_BRCA/TG")
#   source("TG_RNAseq_drop_NTC_rep1.R")
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

resolve_project_dir <- function() {
  env_dir <- Sys.getenv("TG_RNASEQ_DIR", unset = "")
  candidates <- c(
    env_dir,
    "E:/R/TG_BRCA/TG",
    "E:\\R\\TG_BRCA\\TG",
    getwd()
  )
  candidates <- unique(candidates[nzchar(candidates)])
  for (d in candidates) {
    if (dir.exists(d) && (
      file.exists(file.path(d, "genes.read_group_tracking")) ||
      file.exists(file.path(d, "genes.fpkm_tracking")) ||
      file.exists(file.path(d, "genes.count_tracking")) ||
      file.exists(file.path(d, "read_groups.info"))
    )) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  stop("未找到 Cuffdiff 原始数据目录。请把工作目录设为 E:/R/TG_BRCA/TG")
}

norm_key <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

is_ntc_name <- function(x) {
  n <- norm_key(x)
  grepl("NTC|SHNC|NEGCTRL", n) ||
    grepl("CTRL|CONTROL", n) ||
    grepl("^NC[0-9]*$", n)
}

is_rep1_token <- function(x) {
  r <- norm_key(x)
  r %in% c("1", "REP1")
}

is_rep0_token <- function(x) {
  r <- norm_key(x)
  r %in% c("0", "REP0")
}

# 条件名本身带 NTC_rep1，或 NTC 条件下的 replicate=1
row_is_ntc_rep1 <- function(condition, replicate = NA) {
  cond <- as.character(condition)
  n <- norm_key(cond)
  if (is_ntc_name(cond) && grepl("REP1", n)) return(TRUE)
  if (is_ntc_name(cond) && grepl("REP0", n)) return(FALSE)
  if (is_ntc_name(cond) && is_rep1_token(replicate)) return(TRUE)
  FALSE
}

row_is_ntc_rep0 <- function(condition, replicate = NA) {
  cond <- as.character(condition)
  n <- norm_key(cond)
  if (is_ntc_name(cond) && grepl("REP0", n)) return(TRUE)
  if (is_ntc_name(cond) && grepl("REP1", n)) return(FALSE)
  if (is_ntc_name(cond) && is_rep0_token(replicate)) return(TRUE)
  FALSE
}

read_cuff <- function(path) {
  utils::read.delim(
    path, check.names = FALSE, stringsAsFactors = FALSE,
    colClasses = "character", quote = "", comment.char = "",
    na.strings = NULL
  )
}

write_cuff <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

replicate_col <- function(df) {
  nms <- names(df)
  hit <- intersect(c("replicate_num", "replicate", "rep", "replicate_id"), nms)
  if (length(hit) > 0) return(hit[1])
  NA_character_
}

condition_col <- function(df) {
  nms <- names(df)
  hit <- intersect(c("condition", "sample", "group"), nms)
  if (length(hit) > 0) return(hit[1])
  NA_character_
}

log_msg <- function(...) {
  cat(paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = "")), "\n")
}

is_cuffdiff_file <- function(fname) {
  grepl("\\.(info|diff)$", fname, ignore.case = TRUE) ||
    grepl("tracking$", fname, ignore.case = TRUE)
}

q_parts <- function(nm) {
  m <- regexec("^q([0-9]+)(_.*)$", nm)
  g <- regmatches(nm, m)[[1]]
  if (length(g) < 3) return(NULL)
  list(idx = as.integer(g[2]), suffix = g[3])
}

drop_and_reindex_q <- function(df, drop_idx, n_cond) {
  if (length(drop_idx) == 0) return(df)
  keep_old <- setdiff(seq_len(n_cond) - 1L, drop_idx)
  new_of_old <- setNames(seq_along(keep_old) - 1L, as.character(keep_old))
  keep <- rep(TRUE, ncol(df))
  new_names <- names(df)
  for (i in seq_along(new_names)) {
    p <- q_parts(new_names[i])
    if (is.null(p)) next
    if (p$idx %in% drop_idx) {
      keep[i] <- FALSE
    } else if (as.character(p$idx) %in% names(new_of_old)) {
      new_names[i] <- paste0("q", unname(new_of_old[[as.character(p$idx)]]), p$suffix)
    }
  }
  out <- df[, keep, drop = FALSE]
  names(out) <- new_names[keep]
  out
}

safe_log2fc <- function(v1, v2) {
  a <- suppressWarnings(as.numeric(v1))
  b <- suppressWarnings(as.numeric(v2))
  out <- rep("-", length(a))
  ok <- is.finite(a) & is.finite(b) & a > 0 & b > 0
  out[ok] <- format(log2(b[ok] / a[ok]), scientific = FALSE, digits = 12, trim = TRUE)
  out[is.finite(a) & is.finite(b) & a <= 0 & b > 0] <- "inf"
  out[is.finite(a) & is.finite(b) & a > 0 & b <= 0] <- "-inf"
  out[is.finite(a) & is.finite(b) & a <= 0 & b <= 0] <- "0"
  out
}

filter_long_by_ntc_rep1 <- function(df) {
  cc <- condition_col(df)
  rc <- replicate_col(df)
  if (is.na(cc)) return(df)
  cond <- df[[cc]]
  repv <- if (is.na(rc)) rep(NA_character_, nrow(df)) else df[[rc]]
  drop <- mapply(row_is_ntc_rep1, cond, repv, USE.NAMES = FALSE)
  df[!drop, , drop = FALSE]
}

name_is_dropped_condition <- function(x, dropped_conditions) {
  z <- as.character(x)
  z %in% dropped_conditions ||
    norm_key(z) %in% norm_key(dropped_conditions) ||
    (is_ntc_name(z) && grepl("REP1", norm_key(z)))
}

filter_diff_by_dropped_conditions <- function(df, dropped_conditions) {
  if (length(dropped_conditions) == 0) return(df)
  pair_cols <- intersect(c("sample_1", "sample_2", "sample1", "sample2",
                           "condition_1", "condition_2"), names(df))
  if (length(pair_cols) == 0) {
    cc <- condition_col(df)
    if (is.na(cc)) return(df)
    keep <- !vapply(df[[cc]], name_is_dropped_condition, logical(1),
                    dropped_conditions = dropped_conditions)
    return(df[keep, , drop = FALSE])
  }
  keep <- rep(TRUE, nrow(df))
  for (col in pair_cols) {
    keep <- keep & !vapply(df[[col]], name_is_dropped_condition, logical(1),
                           dropped_conditions = dropped_conditions)
  }
  df[keep, , drop = FALSE]
}

replace_q_with_remaining_rep <- function(wide, rg, cond_name, q_idx) {
  if (is.null(rg) || !("tracking_id" %in% names(wide)) || !("tracking_id" %in% names(rg))) {
    return(wide)
  }
  sub <- rg[as.character(rg$condition) == as.character(cond_name), , drop = FALSE]
  if (nrow(sub) == 0) return(wide)
  m <- match(wide$tracking_id, sub$tracking_id)
  hit <- !is.na(m)
  if (!any(hit)) return(wide)
  fpkm_col <- paste0("q", q_idx, "_FPKM")
  count_col <- paste0("q", q_idx, "_count")
  if (fpkm_col %in% names(wide) && "FPKM" %in% names(sub)) {
    wide[[fpkm_col]][hit] <- sub$FPKM[m[hit]]
  }
  if (count_col %in% names(wide) && "raw_frags" %in% names(sub)) {
    wide[[count_col]][hit] <- sub$raw_frags[m[hit]]
  }
  wide
}

update_diff_values_from_fpkm <- function(diff_df, fpkm_df, cond_order) {
  if (is.null(fpkm_df) || !("tracking_id" %in% names(fpkm_df))) return(diff_df)
  id_col <- intersect(c("test_id", "tracking_id", "gene_id"), names(diff_df))[1]
  if (is.na(id_col)) return(diff_df)
  s1 <- intersect(c("sample_1", "sample1"), names(diff_df))[1]
  s2 <- intersect(c("sample_2", "sample2"), names(diff_df))[1]
  v1 <- intersect(c("value_1", "value1"), names(diff_df))[1]
  v2 <- intersect(c("value_2", "value2"), names(diff_df))[1]
  lfc <- grep("log2", names(diff_df), value = TRUE, ignore.case = TRUE)[1]
  if (anyNA(c(s1, s2, v1, v2))) return(diff_df)
  q_of <- function(cond) {
    i <- match(as.character(cond), as.character(cond_order))
    if (is.na(i)) return(NA_integer_)
    i - 1L
  }
  for (i in seq_len(nrow(diff_df))) {
    id <- diff_df[[id_col]][i]
    row <- match(id, fpkm_df$tracking_id)
    if (is.na(row) && "gene_id" %in% names(fpkm_df)) {
      row <- match(id, fpkm_df$gene_id)
    }
    if (is.na(row)) next
    q1 <- q_of(diff_df[[s1]][i])
    q2 <- q_of(diff_df[[s2]][i])
    c1 <- if (!is.na(q1)) paste0("q", q1, "_FPKM") else NA_character_
    c2 <- if (!is.na(q2)) paste0("q", q2, "_FPKM") else NA_character_
    if (!is.na(c1) && c1 %in% names(fpkm_df)) diff_df[[v1]][i] <- fpkm_df[[c1]][row]
    if (!is.na(c2) && c2 %in% names(fpkm_df)) diff_df[[v2]][i] <- fpkm_df[[c2]][row]
    if (!is.na(lfc)) diff_df[[lfc]][i] <- safe_log2fc(diff_df[[v1]][i], diff_df[[v2]][i])
    for (stat_col in intersect(c("test_stat", "p_value", "q_value"), names(diff_df))) {
      diff_df[[stat_col]][i] <- "-"
    }
    if ("significant" %in% names(diff_df)) diff_df[["significant"]][i] <- "no"
  }
  diff_df
}

rg_file_for <- function(fname) {
  if (grepl("^genes\\.", fname)) return("genes.read_group_tracking")
  if (grepl("^isoforms\\.", fname)) return("isoforms.read_group_tracking")
  if (grepl("^cds", fname)) return("cds.read_group_tracking")
  if (grepl("^tss_groups\\.", fname)) return("tss_groups.read_group_tracking")
  "genes.read_group_tracking"
}

fpkm_file_for_diff <- function(fname) {
  if (grepl("isoform", fname, ignore.case = TRUE)) return("isoforms.fpkm_tracking")
  if (grepl("^cds", fname, ignore.case = TRUE)) return("cds.fpkm_tracking")
  if (grepl("tss", fname, ignore.case = TRUE)) return("tss_groups.fpkm_tracking")
  "genes.fpkm_tracking"
}

# -----------------------------------------------------------------------------
# 主流程
# -----------------------------------------------------------------------------
project_dir <- resolve_project_dir()
out_dir <- file.path(project_dir, "without_NTC_rep1")
if (identical(normalizePath(out_dir, winslash = "/", mustWork = FALSE),
              normalizePath(project_dir, winslash = "/", mustWork = FALSE))) {
  stop("输出目录不能与输入目录相同")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_msg("Input:  ", project_dir)
log_msg("Output: ", out_dir)

info_path <- file.path(project_dir, "read_groups.info")
rg_path <- file.path(project_dir, "genes.read_group_tracking")
if (file.exists(info_path)) {
  sample_tab <- read_cuff(info_path)
} else if (file.exists(rg_path)) {
  raw_rg <- read_cuff(rg_path)
  sample_tab <- unique(raw_rg[, intersect(c("condition", "replicate"), names(raw_rg)), drop = FALSE])
} else {
  stop("需要 read_groups.info 或 genes.read_group_tracking 才能识别 NTC_rep1")
}

cc <- condition_col(sample_tab)
rc <- replicate_col(sample_tab)
if (is.na(cc)) stop("样本表没有 condition 列")
repv <- if (is.na(rc)) rep("0", nrow(sample_tab)) else sample_tab[[rc]]
drop_mask <- mapply(row_is_ntc_rep1, sample_tab[[cc]], repv, USE.NAMES = FALSE)
if (!any(drop_mask)) {
  stop("没有识别到 NTC_rep1。请检查 read_groups.info / genes.read_group_tracking 的 condition 与 replicate")
}

cond_order <- unique(as.character(sample_tab[[cc]]))
n_cond <- length(cond_order)
keep_tab <- sample_tab[!drop_mask, , drop = FALSE]
keep_conditions <- unique(as.character(keep_tab[[cc]]))
dropped_conditions <- setdiff(cond_order, keep_conditions)
drop_q_idx <- which(cond_order %in% dropped_conditions) - 1L
remaining_cond_order <- cond_order[!(cond_order %in% dropped_conditions)]

log_msg("Conditions in:  ", paste(cond_order, collapse = ", "))
log_msg("Drop NTC_rep1 rows: ", sum(drop_mask), " / ", nrow(sample_tab))
if (length(dropped_conditions) > 0) {
  log_msg("Dropped whole condition(s): ", paste(dropped_conditions, collapse = ", "))
  log_msg("q indices dropped: ", paste(drop_q_idx, collapse = ", "))
} else {
  log_msg("NTC_rep1 is a replicate of remaining NTC condition; condition-level files kept, NTC values replaced by NTC_rep0")
}
log_msg("Conditions out: ", paste(remaining_cond_order, collapse = ", "))

files <- list.files(project_dir, full.names = FALSE)
files <- files[is_cuffdiff_file(files)]
files <- files[!grepl("^without_NTC_rep1$", files)]
file_rank <- function(fname) {
  if (identical(fname, "read_groups.info")) return(1L)
  if (grepl("read_group_tracking$", fname)) return(2L)
  if (grepl("fpkm_tracking$|count_tracking$", fname)) return(3L)
  if (grepl("\\.diff$", fname)) return(4L)
  5L
}
files <- files[order(file_rank(files), files)]

# 先处理 long 表，供后续用剩余 replicate 替换 NTC 条件值
written <- character()
rg_cache <- list()

for (fname in files) {
  src <- file.path(project_dir, fname)
  if (!file.exists(src) || dir.exists(src)) next
  df <- tryCatch(read_cuff(src), error = function(e) {
    log_msg("Skip (read failed) ", fname, ": ", e$message)
    NULL
  })
  if (is.null(df) || ncol(df) == 0) next

  is_rg_long <- grepl("read_group_tracking$", fname) || identical(fname, "read_groups.info")
  is_wide_q <- any(grepl("^q[0-9]+_", names(df)))
  is_diff <- grepl("\\.diff$", fname) && !is_wide_q

  if (is_rg_long) {
    out_df <- filter_long_by_ntc_rep1(df)
    if (grepl("read_group_tracking$", fname)) rg_cache[[fname]] <- out_df
  } else if (is_wide_q) {
    out_df <- drop_and_reindex_q(df, drop_q_idx, n_cond)
    if (length(dropped_conditions) == 0) {
      ntc_cond <- unique(as.character(keep_tab[[cc]][mapply(
        row_is_ntc_rep0, keep_tab[[cc]],
        if (is.na(rc)) rep("0", nrow(keep_tab)) else keep_tab[[rc]],
        USE.NAMES = FALSE
      )]))
      if (length(ntc_cond) == 1) {
        q_ntc <- match(ntc_cond, remaining_cond_order) - 1L
        rg_name <- rg_file_for(fname)
        rg_use <- rg_cache[[rg_name]]
        if (is.null(rg_use) && file.exists(file.path(out_dir, rg_name))) {
          rg_use <- read_cuff(file.path(out_dir, rg_name))
        }
        if (is.null(rg_use) && file.exists(file.path(project_dir, rg_name))) {
          rg_use <- filter_long_by_ntc_rep1(read_cuff(file.path(project_dir, rg_name)))
        }
        if (!is.null(rg_use) && !is.na(q_ntc)) {
          out_df <- replace_q_with_remaining_rep(out_df, rg_use, ntc_cond, q_ntc)
        }
      }
    }
  } else if (is_diff) {
    out_df <- filter_diff_by_dropped_conditions(df, dropped_conditions)
    if (length(dropped_conditions) == 0) {
      fpkm_name <- fpkm_file_for_diff(fname)
      fpkm_use <- NULL
      if (file.exists(file.path(out_dir, fpkm_name))) {
        fpkm_use <- read_cuff(file.path(out_dir, fpkm_name))
      }
      if (!is.null(fpkm_use)) {
        out_df <- update_diff_values_from_fpkm(out_df, fpkm_use, remaining_cond_order)
      }
    }
  } else {
    cc2 <- condition_col(df)
    rc2 <- replicate_col(df)
    if (!is.na(cc2)) {
      out_df <- filter_long_by_ntc_rep1(df)
      if (length(dropped_conditions) > 0) {
        out_df <- out_df[!(out_df[[cc2]] %in% dropped_conditions), , drop = FALSE]
      }
    } else {
      out_df <- df
    }
  }

  write_cuff(out_df, file.path(out_dir, fname))
  written <- c(written, fname)
  log_msg(fname, " : ", nrow(df), " -> ", nrow(out_df), " rows, ", ncol(df), " -> ", ncol(out_df), " cols")
}

readme <- c(
  "Cuffdiff 原始结果已去掉 NTC_rep1，格式与原目录相同。",
  paste("输入:", project_dir),
  paste("输出:", out_dir),
  paste("原条件:", paste(cond_order, collapse = ", ")),
  paste("保留条件:", paste(remaining_cond_order, collapse = ", ")),
  if (length(dropped_conditions) > 0) {
    paste("删除的整个 condition:", paste(dropped_conditions, collapse = ", "))
  } else {
    "NTC_rep1 是 NTC 的一个 replicate：已从 *.read_group_tracking 和 read_groups.info 删除该行；条件水平 FPKM/count 已换成 NTC_rep0。gene_exp.diff 的 value/log2FC 已按剩余样品重算，p/q 记为 -（不再是 Cuffdiff 原检验）。"
  },
  "原目录未被修改。分析请改指向本文件夹，或把 TG_RNASEQ_DIR 设为本文件夹。"
)
writeLines(readme, file.path(out_dir, "00_REMOVED_NTC_rep1.txt"))
log_msg("Wrote ", length(written), " files")
log_msg("Done. New raw data: ", out_dir)
