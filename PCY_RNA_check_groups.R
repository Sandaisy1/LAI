#!/usr/bin/env Rscript
# 只查看 Cuffdiff 里的组别名称，不做差异分析。
#
# source("E:/R/PCY_RNA/PCY_RNA_check_groups.R", encoding = "UTF-8")
# 或
# Rscript PCY_RNA_check_groups.R "E:/R/PCY_RNA"

options(stringsAsFactors = FALSE)

find_data_dir <- function() {
  cmd <- commandArgs(trailingOnly = TRUE)
  candidates <- c(
    Sys.getenv("PCY_RNA_DIR", unset = ""),
    if (length(cmd) >= 1) cmd[[1]] else "",
    "E:/R/PCY_RNA",
    "E:\\R\\PCY_RNA",
    getwd()
  )
  candidates <- unique(candidates[nzchar(candidates)])
  for (d in candidates) {
    if (dir.exists(d) && any(file.exists(file.path(d, c(
      "genes.read_group_tracking", "genes.count_tracking", "genes.fpkm_tracking", "read_groups.info"
    ))))) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  stop("未找到 Cuffdiff 文件。请把脚本放到 E:/R/PCY_RNA 再运行。", call. = FALSE)
}

read_table <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
}

say <- function(...) {
  message(paste0(...))
}

print_block <- function(title, x) {
  say("")
  say("---- ", title, " ----")
  if (is.null(x)) {
    say("（没有）")
    return(invisible(NULL))
  }
  print(x)
  invisible(x)
}

report_read_group <- function(path) {
  if (!file.exists(path)) {
    say("没有文件: ", path)
    return(invisible(NULL))
  }
  rg <- read_table(path)
  say("文件: ", path)
  say("列名: ", paste(names(rg), collapse = ", "))
  if (!all(c("condition", "replicate") %in% names(rg))) {
    say("这个表没有 condition / replicate 列，下面只列出列名。")
    return(invisible(rg))
  }
  cond <- as.character(rg$condition)
  rep <- as.character(rg$replicate)
  say("")
  say("condition 的全部名称（原始写法）:")
  cond_n <- sort(unique(cond))
  for (nm in cond_n) {
    say("  ", nm, "    行数=", sum(cond == nm, na.rm = TRUE))
  }
  sample_id <- paste(cond, rep, sep = "_rep")
  tab <- as.data.frame(table(condition = cond, replicate = rep), stringsAsFactors = FALSE)
  tab <- tab[tab$Freq > 0, , drop = FALSE]
  tab$sample_id <- paste(tab$condition, tab$replicate, sep = "_rep")
  names(tab)[names(tab) == "Freq"] <- "n_rows"
  tab <- tab[order(tab$condition, tab$replicate), c("condition", "replicate", "sample_id", "n_rows")]
  rownames(tab) <- NULL
  print_block("每个 condition + replicate（这就是样品）", tab)
  say("不同样品数: ", length(unique(sample_id)))
  invisible(tab)
}

report_read_groups_info <- function(path) {
  if (!file.exists(path)) {
    say("没有文件: ", basename(path))
    return(invisible(NULL))
  }
  info <- read_table(path)
  say("文件: ", path)
  say("列名: ", paste(names(info), collapse = ", "))
  show_n <- min(30L, nrow(info))
  print_block(paste0("read_groups.info 前 ", show_n, " 行"), utils::head(info, show_n))
  invisible(info)
}

report_wide_columns <- function(path, pattern) {
  if (!file.exists(path)) {
    say("没有文件: ", basename(path))
    return(invisible(NULL))
  }
  hdr <- names(read_table(path))
  hit <- grep(pattern, hdr, value = TRUE, ignore.case = TRUE)
  hit <- hit[!grepl("variance|conf|status|dispersion|uncertainty", hit, ignore.case = TRUE)]
  say("文件: ", path)
  if (length(hit) == 0) {
    say("没有匹配到样品列。全部列名: ", paste(hdr, collapse = ", "))
    return(invisible(hdr))
  }
  say("样品列:")
  for (nm in hit) say("  ", nm)
  invisible(hit)
}

check_groups <- function(dir = NULL) {
  if (is.null(dir)) dir <- find_data_dir()
  say("数据目录: ", dir)
  say("下面打印的是文件里的原始组别名，没有改名。")
  say("")
  say("======== genes.read_group_tracking ========")
  report_read_group(file.path(dir, "genes.read_group_tracking"))
  say("")
  say("======== read_groups.info ========")
  report_read_groups_info(file.path(dir, "read_groups.info"))
  say("")
  say("======== genes.count_tracking 的样品列 ========")
  report_wide_columns(file.path(dir, "genes.count_tracking"), "_count$|^q[0-9]+_count$")
  say("")
  say("======== genes.fpkm_tracking 的样品列 ========")
  report_wide_columns(file.path(dir, "genes.fpkm_tracking"), "_FPKM$|^q[0-9]+_FPKM$")
  invisible(dir)
}

check_groups()
