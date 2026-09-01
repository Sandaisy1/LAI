#!/usr/bin/env Rscript
# =============================================================================
# JZ 免疫细胞亚群降维：函数库（不要单独 source 本文件）
# 这是 E:/R/fuction of cell-wjz 的方案，比较 JZ-AB / JZ-EVB。
# 圈门、降维、去极端生物学重复与 E:/R/fuction of cell 的免疫亚群分析一致，
# 但文件、目录、组别完全独立：禁止 source Flow_dimred_pipeline.R、ICI_* 或 JY_*。
# 入口：source("JZ_Flow_dimred_pipeline.R")
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANGUAGE = "en")
user_lib <- path.expand("~/R/library")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

# -----------------------------------------------------------------------------
# 0. 依赖
# -----------------------------------------------------------------------------
cran_required <- c("ggplot2", "dplyr", "jsonlite", "RColorBrewer", "scales")
cran_optional <- c("uwot", "Rtsne", "pheatmap", "cowplot", "ggrepel", "matrixStats")
bioc_optional <- c("flowCore", "FlowSOM")

install_if_missing <- function(pkgs, bioc = FALSE, required = TRUE) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) == 0) return(invisible(TRUE))
  message("Installing: ", paste(miss, collapse = ", "))
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      tryCatch(install.packages("BiocManager", repos = "https://cloud.r-project.org"),
               error = function(e) message("BiocManager install failed: ", e$message))
    }
    if (requireNamespace("BiocManager", quietly = TRUE)) {
      tryCatch(BiocManager::install(miss, update = FALSE, ask = FALSE),
               error = function(e) message("Bioconductor install failed: ", e$message))
    }
  } else {
    tryCatch(install.packages(miss, repos = "https://cloud.r-project.org"),
             error = function(e) message("CRAN install failed: ", e$message))
  }
  still <- miss[!vapply(miss, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still) > 0 && required) {
    stop("缺少必需 R 包: ", paste(still, collapse = ", "))
  }
  if (length(still) > 0) message("可选包未安装，相关步骤将跳过或降级: ", paste(still, collapse = ", "))
  invisible(TRUE)
}

install_if_missing(cran_required, bioc = FALSE, required = TRUE)
install_if_missing(cran_optional, bioc = FALSE, required = FALSE)

ensure_flowcore <- function() {
  if (requireNamespace("flowCore", quietly = TRUE)) return(invisible(TRUE))
  message("正在安装 flowCore（读 FCS 必需），可能要几分钟…")
  install_if_missing("flowCore", bioc = TRUE, required = FALSE)
  if (requireNamespace("flowCore", quietly = TRUE)) return(invisible(TRUE))
  stop(
    "读 FCS 需要 Bioconductor 包 flowCore。请先在 R 控制台运行这两行，装好后再 source 脚本：\n",
    "  install.packages(\"BiocManager\")\n",
    "  BiocManager::install(\"flowCore\")"
  )
}

safe_library <- function(pkgs) {
  for (p in pkgs) {
    if (requireNamespace(p, quietly = TRUE)) {
      suppressPackageStartupMessages(library(p, character.only = TRUE))
    }
  }
}
safe_library(c(cran_required, cran_optional, bioc_optional))
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

# -----------------------------------------------------------------------------
# 1. 路径
# -----------------------------------------------------------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)))
  }
  ofile <- NULL
  n <- sys.nframe()
  if (n > 0) {
    for (i in n:1) {
      ofile <- sys.frame(i)$ofile
      if (!is.null(ofile)) break
    }
  }
  if (!is.null(ofile) && nzchar(ofile)) {
    return(dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir()

# 读写目录：只在 E:/R/fuction of cell-wjz；不要回退到 fuction of cell 或 Internation
flow_primary_data_dir <- "E:/R/fuction of cell-wjz"
flow_legacy_data_dir <- "E:/R/fuction of cell-wjz"
flow_legacy_data_dir2 <- "E:/R/fuction of cell-wjz"

jz_is_foreign_dir <- function(p) {
  s <- gsub("\\", "/", as.character(p), fixed = TRUE)
  if (grepl("cell-wjz", s, ignore.case = TRUE)) return(FALSE)
  if (grepl("cell-ljy", s, ignore.case = TRUE)) return(TRUE)
  if (grepl("Internation cell immune", s, ignore.case = TRUE)) return(TRUE)
  grepl("fuction of cell|function of cell", s, ignore.case = TRUE)
}

resolve_flow_dir <- function() {
  env_dir <- Sys.getenv("JZ_FLOW_DIR", unset = "")
  if (!nzchar(env_dir)) env_dir <- Sys.getenv("FLOW_DIR", unset = "")
  preferred <- c(
    env_dir,
    flow_primary_data_dir,
    "E:\\R\\fuction of cell-wjz",
    "E:/R/function of cell-wjz",
    "E:\\R\\function of cell-wjz"
  )
  preferred <- unique(preferred[nzchar(preferred)])
  preferred <- preferred[!vapply(preferred, jz_is_foreign_dir, logical(1))]
  for (d in preferred) {
    if (dir.exists(d)) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  candidates <- c(
    flow_primary_data_dir,
    "E:\\R\\fuction of cell-wjz",
    file.path(script_dir, "fuction of cell-wjz"),
    file.path(script_dir, "function of cell-wjz"),
    script_dir,
    getwd()
  )
  candidates <- unique(candidates[nzchar(candidates)])
  candidates <- candidates[!vapply(candidates, jz_is_foreign_dir, logical(1))]
  for (d in candidates) {
    if (!dir.exists(d)) next
    hits <- list.files(d, pattern = "(?i)_unmixed\\.fcs$", recursive = TRUE)
    if (length(hits) > 0) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  if (dir.exists(flow_primary_data_dir)) {
    return(normalizePath(flow_primary_data_dir, winslash = "/", mustWork = FALSE))
  }
  normalizePath(script_dir, winslash = "/", mustWork = FALSE)
}

project_dir <- resolve_flow_dir()
result_dir <- file.path(project_dir, "results_flow")
log_dir <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("flow_dimred_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

demo_flag <- toupper(Sys.getenv("FLOW_DEMO", unset = "0")) %in% c("1", "TRUE", "YES")
n_per_sample_dr <- as.integer(Sys.getenv("FLOW_CELLS_PER_SAMPLE", unset = "8000"))
if (is.na(n_per_sample_dr) || n_per_sample_dr < 200) n_per_sample_dr <- 8000
asinh_cofactor <- 150
seed_value <- 42

# -----------------------------------------------------------------------------
# 2. 文件名：JZ-EVB / JZ-AB
#    生物学重复 EVB-1/2/3 与 AB-1/2/3；每个生物学重复两个技术重复 EVB1-1、AB1-2
#    必须先匹配 EVB，再匹配 AB，避免 EVB 被拆成 EV 或 AB
# -----------------------------------------------------------------------------
flow_ctrl_group <- "JZ-EVB"
flow_trt_group <- "JZ-AB"
flow_group_levels <- c("JZ-EVB", "JZ-AB")
flow_cohort <- "JZ"

# 比较是 JZ-AB / JZ-EVB（处理 / 对照）。文件名里仍写 EVB、AB。
jz_canon_group <- function(grp) {
  g <- toupper(gsub("[^A-Za-z0-9]", "", as.character(grp)))
  if (!nzchar(g) || is.na(g)) return(NA_character_)
  if (grepl("EVB", g, fixed = TRUE)) return("JZ-EVB")
  if (identical(g, "AB") || identical(g, "JZAB")) return("JZ-AB")
  NA_character_
}

jz_short_arm <- function(grp) {
  g <- jz_canon_group(grp)
  if (identical(g, "JZ-EVB")) return("EVB")
  if (identical(g, "JZ-AB")) return("AB")
  NA_character_
}

# 带技术重复的名字优先，避免把 EVB1-1_P1 的 “1” 当成 panel
flow_re_tech <- "^(?:JZ[_-]?)?(EVB|AB)[-_ ]?([123])[-_ ]([12])[-_ ]+(?:PANEL[-_ ]?)?P?0?([123])[-_ ].*(unmixed|raw)\\.fcs$"
flow_re_bio <- "^(?:JZ[_-]?)?(EVB|AB)[-_ ]?([123])[-_ ]+(?:PANEL[-_ ]?)?P?0?([123])[-_ ].*(unmixed|raw)\\.fcs$"
flow_re_folder_tech <- "^([123])[-_ ]([12])[-_ ]+(?:PANEL[-_ ]?)?P?0?([123])[-_ ].*(unmixed|raw)\\.fcs$"

group_from_path <- function(path) {
  parts <- unlist(strsplit(gsub("\\", "/", as.character(path), fixed = TRUE), "/"))
  parts <- toupper(gsub("[^A-Za-z0-9]", "", parts))
  parts <- parts[nzchar(parts)]
  for (p in rev(parts)) {
    g <- jz_canon_group(p)
    if (!is.na(g) && nzchar(g)) return(g)
  }
  NA_character_
}

parse_fcs_filename <- function(path) {
  b <- basename(path)
  grp <- NA_character_
  bio <- NA_character_
  tech <- NA_character_
  panel_n <- NA_character_
  kind <- NA_character_
  m <- regexec(flow_re_tech, b, ignore.case = TRUE)
  hit <- regmatches(b, m)[[1]]
  if (length(hit) >= 6) {
    grp <- toupper(hit[2])
    bio <- hit[3]
    tech <- hit[4]
    panel_n <- hit[5]
    kind <- tolower(hit[6])
  } else {
    m <- regexec(flow_re_bio, b, ignore.case = TRUE)
    hit <- regmatches(b, m)[[1]]
    if (length(hit) >= 5) {
      grp <- toupper(hit[2])
      bio <- hit[3]
      panel_n <- hit[4]
      kind <- tolower(hit[5])
    } else {
      m <- regexec(flow_re_folder_tech, b, ignore.case = TRUE)
      hit <- regmatches(b, m)[[1]]
      if (length(hit) >= 5) {
        grp <- group_from_path(path)
        bio <- hit[2]
        tech <- hit[3]
        panel_n <- hit[4]
        kind <- tolower(hit[5])
      }
    }
  }
  if (is.na(grp) || !nzchar(grp) || is.na(bio) || !nzchar(bio) ||
      is.na(panel_n) || !nzchar(panel_n)) {
    return(NULL)
  }
  grp <- jz_canon_group(grp)
  if (is.na(grp) || !nzchar(grp)) return(NULL)
  arm <- jz_short_arm(grp)
  # 组别 JZ-EVB / JZ-AB；生物学重复 EVB-1；技术管 EVB1-1
  bio_sample <- paste0(arm, "-", bio)
  sample <- if (!is.na(tech) && nzchar(tech)) paste0(arm, bio, "-", tech) else bio_sample
  list(
    file = b,
    path = path,
    cohort = flow_cohort,
    group = grp,
    replicate = bio,
    tech_rep = if (!is.na(tech) && nzchar(tech)) tech else NA_character_,
    bio_sample = bio_sample,
    sample = sample,
    panel = paste0("P", panel_n),
    kind = kind
  )
}

list_unmixed_files <- function(root) {
  grab <- function(rec) {
    list.files(root, pattern = "(?i)unmixed\\.fcs$", full.names = TRUE, recursive = rec)
  }
  files <- unique(c(grab(FALSE), grab(TRUE)))
  files <- files[!grepl("(^|\\\\|/)results_flow(\\\\|/|$)", files, ignore.case = TRUE)]
  meta <- lapply(files, parse_fcs_filename)
  ok <- !vapply(meta, is.null, logical(1))
  if (any(!ok) && length(files) > 0) {
    log_msg("Skip unmatched unmixed names: ", paste(basename(files[!ok]), collapse = ", "))
  }
  if (!any(ok)) return(data.frame())
  rows <- lapply(meta[ok], function(x) {
    data.frame(
      file = x$file, path = x$path, cohort = x$cohort, group = x$group,
      replicate = x$replicate, tech_rep = x$tech_rep,
      bio_sample = x$bio_sample, sample = x$sample,
      panel = x$panel, kind = x$kind, stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  df <- df[!duplicated(df$path), , drop = FALSE]
  tech_ord <- ifelse(is.na(df$tech_rep) | !nzchar(df$tech_rep), "0", df$tech_rep)
  df <- df[order(df$panel, df$group, df$replicate, tech_ord), ]
  rownames(df) <- NULL
  df
}

ensure_bio_sample <- function(df) {
  if (is.null(df) || !nrow(df)) return(df)
  if (!"bio_sample" %in% names(df) || all(is.na(df$bio_sample) | !nzchar(as.character(df$bio_sample)))) {
    df$bio_sample <- as.character(df$sample)
  }
  df$bio_sample <- as.character(df$bio_sample)
  df
}

# 技术重复先在生物学重复内平均；免疫亚群统计再去掉 1 个极端生物学重复后 n=2。
# 不要把 EVB1-1/EVB1-2 当成两个生物学 n。JZ 统计 n=2。
aggregate_freq_by_bio <- function(freq_df, id_col) {
  freq_df <- ensure_bio_sample(freq_df)
  if (!id_col %in% names(freq_df)) stop("aggregate_freq_by_bio: missing ", id_col)
  keep <- unique(c("bio_sample", "group", id_col, "percent"))
  keep <- intersect(keep, names(freq_df))
  d <- freq_df[, keep, drop = FALSE]
  d <- d[is.finite(d$percent), , drop = FALSE]
  if (!nrow(d)) return(d)
  fmla <- stats::as.formula(paste("percent ~ bio_sample + group +", id_col))
  aggregate(fmla, data = d, FUN = mean, na.rm = TRUE)
}

# 免疫亚群降维：三个生物学重复里去掉离中位数更远的那个（最大值或最小值），
# 剩下的用于均值/SD/检验。降维仍用全部管子。n<3 时不去。
# JZ 与免疫亚群方案一样去极端值；ICI 才是 n=3。不要把本方案改成 ICI。
if (!exists("flow_trim_bio_extremes", inherits = TRUE)) {
  flow_trim_bio_extremes <- TRUE
}

flow_should_trim_bio <- function(trim_bio = NULL) {
  if (!is.null(trim_bio)) return(isTRUE(trim_bio))
  isTRUE(get0("flow_trim_bio_extremes", ifnotfound = TRUE, inherits = TRUE))
}

which_extreme_bio <- function(vals) {
  x <- as.numeric(vals)
  finite <- is.finite(x)
  if (sum(finite) < 3L) return(integer(0))
  idx <- which(finite)
  xx <- x[idx]
  med <- stats::median(xx)
  i_min <- which.min(xx)
  i_max <- which.max(xx)
  dmin <- med - xx[i_min]
  dmax <- xx[i_max] - med
  if (!is.finite(dmin) || !is.finite(dmax)) return(integer(0))
  drop_local <- if (dmax >= dmin) i_max else i_min
  idx[drop_local]
}

trim_bio_extremes <- function(freq_df, id_col, group_col = "group") {
  if (is.null(freq_df) || !nrow(freq_df)) return(freq_df)
  if (!id_col %in% names(freq_df)) stop("trim_bio_extremes: missing ", id_col)
  ids <- unique(as.character(freq_df[[id_col]]))
  grps <- unique(as.character(freq_df[[group_col]]))
  keep_rows <- rep(TRUE, nrow(freq_df))
  dropped <- list()
  for (id in ids) {
    for (g in grps) {
      ii <- which(as.character(freq_df[[id_col]]) == id &
                    as.character(freq_df[[group_col]]) == g)
      if (length(ii) < 3L) next
      drop_local <- which_extreme_bio(freq_df$percent[ii])
      if (!length(drop_local)) next
      row_i <- ii[drop_local]
      keep_rows[row_i] <- FALSE
      bio <- if ("bio_sample" %in% names(freq_df)) {
        as.character(freq_df$bio_sample[row_i])
      } else {
        as.character(row_i)
      }
      dropped[[length(dropped) + 1L]] <- data.frame(
        id = id,
        group = g,
        dropped_bio = bio,
        dropped_percent = freq_df$percent[row_i],
        stringsAsFactors = FALSE
      )
    }
  }
  out <- freq_df[keep_rows, , drop = FALSE]
  attr(out, "dropped") <- if (length(dropped)) do.call(rbind, dropped) else NULL
  out
}

maybe_trim_bio <- function(freq_df, id_col, trim_bio = NULL) {
  if (!flow_should_trim_bio(trim_bio)) return(freq_df)
  trim_bio_extremes(freq_df, id_col)
}

bio_percent_table <- function(samp, trim_bio = NULL) {
  s2 <- ensure_bio_sample(samp)
  s2 <- s2[is.finite(s2$percent), , drop = FALSE]
  if (!nrow(s2)) return(s2)
  out <- aggregate(percent ~ bio_sample + group, data = s2, FUN = mean, na.rm = TRUE)
  out$sample <- out$bio_sample
  if (!flow_should_trim_bio(trim_bio)) return(out)
  out$.trim_id <- "subset"
  out <- trim_bio_extremes(out, ".trim_id")
  dropped <- attr(out, "dropped")
  out$.trim_id <- NULL
  attr(out, "dropped") <- dropped
  out
}

# -----------------------------------------------------------------------------
# 3. Panel 地图
# Windows 常把“另存为”弄成 JZ_flow_panel_map.json.txt（资源管理器显示“文本文档”）
# -----------------------------------------------------------------------------
panel_map_search_dirs <- function() {
  env_dir <- Sys.getenv("JZ_FLOW_DIR", unset = "")
  if (!nzchar(env_dir)) env_dir <- Sys.getenv("FLOW_DIR", unset = "")
  dirs <- c(
    getwd(),
    script_dir,
    project_dir,
    flow_primary_data_dir,
    "E:\\R\\fuction of cell-wjz",
    "E:/R/function of cell-wjz",
    env_dir
  )
  dirs <- unique(dirs[nzchar(dirs)])
  dirs <- dirs[!vapply(dirs, jz_is_foreign_dir, logical(1))]
  dirs[dir.exists(dirs)]
}

find_panel_map_file <- function(dirs = panel_map_search_dirs()) {
  dirs <- dirs[nzchar(dirs) & dir.exists(dirs)]
  dirs <- dirs[!vapply(dirs, jz_is_foreign_dir, logical(1))]
  names <- c(
    "JZ_flow_panel_map.json",
    "JZ_flow_panel_map.json.txt",
    "JZ_flow_panel_map.txt"
  )
  exact <- unlist(lapply(dirs, function(d) file.path(d, names)), use.names = FALSE)
  hit <- exact[file.exists(exact)]
  if (length(hit) > 0) return(hit[[1]])
  fuzzy <- unlist(lapply(dirs, function(d) {
    list.files(d, pattern = "^JZ_flow_panel_map", ignore.case = TRUE, full.names = TRUE)
  }), use.names = FALSE)
  if (length(fuzzy) > 0) return(fuzzy[[1]])
  NULL
}

read_panel_map_file <- function(path) {
  parsed <- tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE), error = function(e) NULL)
  if (!is.null(parsed) && !is.null(parsed$panels)) return(parsed)
  con <- file(path, open = "r", encoding = "UTF-16LE")
  on.exit(close(con), add = TRUE)
  txt <- paste(readLines(con, warn = FALSE), collapse = "\n")
  jsonlite::fromJSON(txt, simplifyVector = FALSE)
}

load_panel_map <- function() {
  hit <- find_panel_map_file()
  if (!is.null(hit)) {
    log_msg("Panel map: ", hit)
    return(read_panel_map_file(hit))
  }
  dirs <- panel_map_search_dirs()
  listing <- paste(vapply(dirs, function(d) {
    paste0("  ", d, " -> ", paste(list.files(d), collapse = ", "))
  }, character(1)), collapse = "\n")
  stop(
    "找不到 JZ_flow_panel_map.json。\n",
    "Windows 若类型是“文本文档”、图标是记事本，真实文件名多半是 JZ_flow_panel_map.json.txt。\n",
    "在资源管理器：查看 → 去掉“隐藏已知文件类型的扩展名”，再改名为 JZ_flow_panel_map.json。\n",
    "已搜索：\n", listing
  )
}

panel_map <- load_panel_map()
if (!is.null(panel_map$qc$asinh_cofactor)) {
  asinh_cofactor <- as.numeric(panel_map$qc$asinh_cofactor)
}
if (!is.null(panel_map$groups) && length(panel_map$groups) >= 2) {
  flow_ctrl_group <- as.character(panel_map$groups[[1]])
  flow_trt_group <- as.character(panel_map$groups[[2]])
  flow_group_levels <- c(flow_ctrl_group, flow_trt_group)
}

norm_id <- function(x) {
  x <- toupper(as.character(x))
  x[is.na(x)] <- ""
  x <- gsub("\u03B1|\u0391", "A", x)
  x <- gsub("\u03B3|\u0393", "G", x)
  x <- gsub("\u03B2|\u0392", "B", x)
  x <- gsub("\u03B5|\u0395", "E", x)
  x <- gsub("GAMMA", "G", x)
  x <- gsub("ALPHA", "A", x)
  gsub("[^A-Z0-9]", "", x)
}

marker_aliases <- function(item) {
  ids <- c(item$marker, item$fluorochrome, item$aliases)
  fl <- as.character(item$fluorochrome)[1]
  if (!is.null(panel_map$fluorochrome_aliases) && nzchar(fl) && !is.null(panel_map$fluorochrome_aliases[[fl]])) {
    ids <- c(ids, unlist(panel_map$fluorochrome_aliases[[fl]], use.names = FALSE))
  }
  unique(norm_id(ids[nzchar(as.character(ids))]))
}

# Cytek/FlowJo: "BUV496-A"、"FJComp-APC-Cy7-A"、"CD45 V500"
fcs_clean_label <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("^(FJCOMP[-_ ]*|COMP[-_ ]*|UNMIX[-_ ]*)", "", x, ignore.case = TRUE)
  x <- gsub("[-_ .]([AHW])$", "", x, ignore.case = TRUE)
  x
}

fcs_tokens <- function(...) {
  raw <- paste(..., collapse = " ")
  raw <- fcs_clean_label(raw)
  raw <- toupper(raw)
  raw <- gsub("\u03B1|\u0391", "A", raw)
  raw <- gsub("\u03B3|\u0393", "G", raw)
  raw <- gsub("\u03B2|\u0392", "B", raw)
  raw <- gsub("\u03B5|\u0395", "E", raw)
  raw <- gsub("GAMMA", "G", raw)
  raw <- gsub("ALPHA", "A", raw)
  parts <- unlist(strsplit(raw, "[^A-Z0-9]+"), use.names = FALSE)
  parts <- parts[nzchar(parts) & !parts %in% c("A", "H", "W")]
  unique(c(parts, if (length(parts)) paste(parts, collapse = "") else character(0)))
}

panel_markers <- function(panel_id) {
  items <- panel_map$panels[[panel_id]]$markers
  if (is.null(items)) stop("Unknown panel: ", panel_id)
  items
}

dr_marker_names <- function(panel_id) {
  items <- panel_markers(panel_id)
  keep <- vapply(items, function(it) {
    !identical(it$compartment, "viability") && !identical(it$marker, "L/D")
  }, logical(1))
  vapply(items[keep], function(it) it$marker, character(1))
}

# -----------------------------------------------------------------------------
# 4. 通道匹配与转换
# -----------------------------------------------------------------------------
match_channels <- function(channel_names, channel_desc, panel_id) {
  items <- panel_markers(panel_id)
  nms <- as.character(channel_names)
  desc <- as.character(channel_desc)
  if (length(desc) != length(nms)) desc <- rep("", length(nms))
  nms[is.na(nms)] <- ""
  desc[is.na(desc)] <- ""
  tok_list <- lapply(seq_along(nms), function(j) fcs_tokens(nms[j], desc[j]))
  used <- rep(FALSE, length(nms))
  mapping <- data.frame(
    marker = vapply(items, function(it) it$marker, character(1)),
    fluorochrome = vapply(items, function(it) it$fluorochrome, character(1)),
    compartment = vapply(items, function(it) it$compartment, character(1)),
    channel = NA_character_,
    channel_index = NA_integer_,
    stringsAsFactors = FALSE
  )
  score_one <- function(item, idx) {
    al <- marker_aliases(item)
    toks <- tok_list[[idx]]
    hit <- al[al %in% toks]
    if (!length(hit)) return(0)
    mx <- max(nchar(hit))
    full <- max(nchar(toks))
    coverage <- mx / max(full, 1)
    base <- mx * 10 + as.integer(round(coverage * 30))
    if (any(hit == norm_id(item$marker))) {
      base <- base + 100
    } else if (any(hit == norm_id(item$fluorochrome))) {
      base <- base + 40
    }
    base
  }
  for (i in seq_along(items)) {
    best <- 0
    best_j <- NA_integer_
    for (j in seq_along(nms)) {
      if (used[j]) next
      if (grepl("^(FSC|SSC|TIME|EVENT)", nms[j], ignore.case = TRUE)) next
      sc <- score_one(items[[i]], j)
      if (sc > best) {
        best <- sc
        best_j <- j
      }
    }
    if (!is.na(best_j) && best >= 50) {
      used[best_j] <- TRUE
      mapping$channel[i] <- nms[best_j]
      mapping$channel_index[i] <- best_j
    }
  }
  mapping
}

format_channel_preview <- function(nms, desc, n = 40) {
  nms <- as.character(nms)
  desc <- as.character(desc)
  if (length(desc) != length(nms)) desc <- rep("", length(nms))
  desc[is.na(desc)] <- ""
  nms[is.na(nms)] <- ""
  nshow <- min(length(nms), n)
  lines <- sprintf("  [%d] name='%s'  desc='%s'", seq_len(nshow), nms[seq_len(nshow)], desc[seq_len(nshow)])
  extra <- if (length(nms) > nshow) paste0("  ... ", length(nms) - nshow, " more") else NULL
  paste(c(lines, extra), collapse = "\n")
}

asinh_mat <- function(mat, cofactor = asinh_cofactor) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  asinh(pmax(mat, 0) / cofactor)
}

scale_markers <- function(mat) {
  mat <- as.matrix(mat)
  mu <- colMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 2, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv < 1e-8] <- 1
  sweep(sweep(mat, 2, mu, "-"), 2, sdv, "/")
}

# 联合 tSNE/UMAP 要能看见谱系：CD4/CD8 比 CD62L/CD44 权重大，否则 naive CD4 和 naive CD8 会揉成一团
dr_lineage_marker_weights <- function(panel_id, markers) {
  markers <- as.character(markers)
  w <- setNames(rep(1, length(markers)), markers)
  core <- c("CD4", "CD8", "CD8b")
  lineage <- c("CD3", "CD19", "NKp46", "NK1.1", "CD11B", "CD11C", "His")
  if (identical(panel_id, "P2")) {
    lineage <- c("CD19", "IgD", "CD27", "IgM", "BLIMP-1")
    core <- character(0)
  }
  w[markers %in% lineage] <- 2.5
  w[markers %in% core] <- 4
  w
}

weight_dr_markers <- function(mat, panel_id) {
  mat <- as.matrix(mat)
  if (!ncol(mat)) return(mat)
  w <- dr_lineage_marker_weights(panel_id, colnames(mat))
  sweep(mat, 2, w, "*")
}

# 类内再降维：压低已经用来定大类的谱系通道，加重这一层的亚群标志（CD62L/CD44 等）
class_dr_marker_weights <- function(major, markers) {
  markers <- as.character(markers)
  w <- setNames(rep(1, length(markers)), markers)
  maj <- as.character(major)[1]
  if (is.na(maj) || !nzchar(maj)) return(w)
  w[markers %in% c("L/D", "CD45")] <- 0.2
  if (maj %in% c("CD4", "CD8")) {
    w[markers %in% c("CD3", "CD4", "CD8", "CD8b", "CD19", "NKp46", "NK1.1", "CD11B")] <- 0.35
    up <- c("CD62L", "CD44", "CD27", "CD95", "CD69", "CD25", "SCA-1",
            "LAG-3", "TIM-3", "PD-L1", "IFN-g", "TNF-a", "GZMB", "Perforin")
    w[markers %in% up] <- 3.2
  } else if (identical(maj, "NK")) {
    w[markers %in% c("CD3", "CD19", "CD4", "CD8", "CD8b")] <- 0.35
    up <- c("NKp46", "NK1.1", "CD27", "CD11B", "CD69", "GZMB", "Perforin",
            "LAG-3", "TIM-3", "PD-L1")
    w[markers %in% up] <- 3.2
  } else if (identical(maj, "NKT")) {
    w[markers %in% c("CD19", "CD11B")] <- 0.4
    up <- c("CD4", "CD8", "CD69", "GZMB", "IFN-g", "NKp46", "NK1.1", "CD44")
    w[markers %in% up] <- 3
  } else if (maj %in% c("B", "Naive_B", "Unswitched_B", "Switched_B")) {
    up <- c("IgD", "CD27", "IgM", "IgG", "BLIMP-1", "CD80", "CD86")
    w[markers %in% up] <- 3.2
  } else if (identical(maj, "Myeloid")) {
    up <- c("CD11B", "CD11C", "LY6C", "LY6G", "F4/80", "Siglec-F", "CCR3",
            "I-A/I-E", "CD103", "CD206", "FceRI")
    w[markers %in% up] <- 3
  }
  w
}

weight_class_dr_markers <- function(mat, major) {
  mat <- as.matrix(mat)
  if (!ncol(mat)) return(mat)
  w <- class_dr_marker_weights(major, colnames(mat))
  sweep(mat, 2, w, "*")
}

# CD3+ T 里按 CD4 vs CD8 象限：CD4+ CD8- / CD4- CD8+。不要比谁更大就把双阳/双阴塞进某一边
assign_t_cd4_cd8 <- function(cd4, cd8, is_t) {
  n <- length(is_t)
  out <- rep(NA_character_, n)
  is_t <- na_to_false(is_t)
  if (!any(is_t)) return(out)
  c4 <- as.numeric(cd4)[is_t]
  c8 <- as.numeric(cd8)[is_t]
  c4[!is.finite(c4)] <- -Inf
  c8[!is.finite(c8)] <- -Inf
  cut4 <- axis_pos_cut(c4)
  cut8 <- axis_pos_cut(c8)
  if (!is.finite(cut4)) cut4 <- 1.2
  if (!is.finite(cut8)) cut8 <- 1.2
  hi4 <- c4 >= cut4
  hi8 <- c8 >= cut8
  lab <- rep("other", sum(is_t))
  lab[hi4 & !hi8] <- "CD4"
  lab[!hi4 & hi8] <- "CD8"
  out[is_t] <- lab
  out
}

# -----------------------------------------------------------------------------
# 5. 读 FCS 或演示数据
# -----------------------------------------------------------------------------
read_fcs_expr <- function(path, panel_id) {
  ensure_flowcore()
  ff <- flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE, emptyValue = FALSE)
  exprs <- flowCore::exprs(ff)
  pdata <- as.data.frame(flowCore::parameters(ff)@data)
  nms <- as.character(pdata$name)
  desc <- if ("desc" %in% names(pdata)) as.character(pdata$desc) else rep("", length(nms))
  expr_nms <- colnames(exprs)
  if (length(expr_nms) == length(nms)) {
    blank <- !nzchar(nms) | is.na(nms)
    nms[blank] <- expr_nms[blank]
    # 有的 Cytek 文件 $PnN 是 BJ1-A，标志物写在 exprs 列名
    weak <- grepl("^(UV|V|B|YG|R|IR|U|G)?[0-9]+", nms, ignore.case = TRUE)
    richer <- nchar(fcs_clean_label(expr_nms)) > nchar(fcs_clean_label(nms))
    use_expr <- weak & richer & !grepl("^(FSC|SSC|TIME|EVENT)", expr_nms, ignore.case = TRUE)
    nms[use_expr] <- expr_nms[use_expr]
  }
  map <- match_channels(nms, desc, panel_id)
  n_hit <- sum(!is.na(map$channel_index))
  log_msg(basename(path), " matched ", n_hit, "/", nrow(map), " markers")
  if (n_hit == 0) {
    preview <- format_channel_preview(nms, desc)
    dump <- file.path(log_dir, paste0(panel_id, "_", tools::file_path_sans_ext(basename(path)), "_channels.csv"))
    tryCatch({
      utils::write.csv(
        data.frame(index = seq_along(nms), name = nms, desc = desc, stringsAsFactors = FALSE),
        dump, row.names = FALSE
      )
      log_msg("Wrote channel dump: ", dump)
    }, error = function(e) NULL)
    stop(
      "没有匹配到任何分析通道。不是 FCS 没读到，是通道名和 panel 地图对不上。\n",
      "Cytek unmixed 常见通道名是 BUV496-A / APC-Cy7-A（带 -A），或 desc 里才写 CD4。\n",
      basename(path), " 里的通道是：\n", preview, "\n",
      "请把这段贴回来；若只有 V1-A/B1-A 这种检测器名，需要在 SpectroFlo 把 fluorescent tag 写成标志物后再导出 unmixed。"
    )
  }
  miss <- map$marker[is.na(map$channel_index)]
  if (length(miss) > 0) {
    log_msg("Unmatched markers: ", paste(miss, collapse = ", "))
    if ("NK1.1" %in% miss) {
      log_msg("NK1.1 is AF700 on the sheet; Cytek unmixed often names that detector R718-A / APC-R700-A / Y710-A. All should match.")
      log_msg("Channel names in this file: ", format_channel_preview(nms, desc))
    }
  }
  list(exprs = exprs, names = nms, desc = desc, map = map)
}

# CD45+：有阴性岛就切在两群之间；全体都阳则只去掉左尾，不要把门切到主团上方。
qc_cd45_cut <- function(v) {
  v <- as.numeric(v)
  v <- v[is.finite(v)]
  if (!length(v)) return(0)
  if (length(v) < 20) return(as.numeric(stats::quantile(v, 0.08, names = FALSE)))
  xs <- scale(as.matrix(v))
  xs[!is.finite(xs)] <- 0
  set.seed(if (exists("seed_value")) seed_value else 42)
  km <- tryCatch(
    stats::kmeans(xs, centers = 2, nstart = 8, iter.max = 200, algorithm = "Lloyd"),
    error = function(e) NULL
  )
  if (!is.null(km) && length(unique(km$cluster)) == 2) {
    med <- tapply(v, km$cluster, median, na.rm = TRUE)
    hi <- unname(med[which.max(med)])
    lo <- unname(med[which.min(med)])
    frac_hi <- mean(km$cluster == as.integer(names(med)[which.max(med)]))
    if (is.finite(hi) && is.finite(lo) && (hi - lo) >= 0.55 && frac_hi >= 0.25 && frac_hi <= 0.97) {
      return(unname((hi + lo) / 2))
    }
  }
  as.numeric(stats::quantile(v, 0.08, names = FALSE))
}

# 清理：单细胞 →（P1 小淋巴；P2 宽单核，活化 B / 浆母更大）→ 活细胞 → CD45+
qc_filter_matrix <- function(exprs, names, map, panel_id = NA) {
  keep <- rep(TRUE, nrow(exprs))
  fsc <- grep("^FSC-A$|^FSC_A$|^FSC\\.A$", names, ignore.case = TRUE)
  ssc <- grep("^SSC-A$|^SSC_A$|^SSC\\.A$", names, ignore.case = TRUE)
  if (length(fsc) && length(ssc)) {
    fs <- exprs[, fsc[1]]
    ss <- exprs[, ssc[1]]
    keep <- keep & fs > quantile(fs, 0.01, na.rm = TRUE) & ss > quantile(ss, 0.01, na.rm = TRUE)
  }
  fsch <- grep("^FSC-H$|^FSC_H$", names, ignore.case = TRUE)
  if (length(fsc) && length(fsch)) {
    ratio <- exprs[, fsc[1]] / pmax(exprs[, fsch[1]], 1)
    keep <- keep & ratio > 0.5 & ratio < 2
  }
  fscw <- grep("^FSC-W$|^FSC_W$", names, ignore.case = TRUE)
  if (length(fscw) && any(keep)) {
    fw <- exprs[, fscw[1]]
    keep <- keep & fw <= quantile(fw[keep], 0.95, na.rm = TRUE)
  }
  if (identical(as.character(panel_id), "P1") && length(fsc) && length(ssc)) {
    fs <- exprs[, fsc[1]]
    ss <- exprs[, ssc[1]]
    keep <- keep & ss <= quantile(ss, 0.82, na.rm = TRUE)
    keep <- keep & fs >= quantile(fs, 0.05, na.rm = TRUE) & fs <= quantile(fs, 0.97, na.rm = TRUE)
  }
  if (identical(as.character(panel_id), "P2") && length(fsc) && length(ssc)) {
    fs <- exprs[, fsc[1]]
    ss <- exprs[, ssc[1]]
    keep <- keep & ss <= quantile(ss, 0.95, na.rm = TRUE)
    keep <- keep & fs >= quantile(fs, 0.02, na.rm = TRUE) & fs <= quantile(fs, 0.99, na.rm = TRUE)
  }
  ld_idx <- map$channel_index[map$marker == "L/D"]
  if (length(ld_idx) && !is.na(ld_idx[1])) {
    ld <- exprs[, ld_idx[1]]
    keep <- keep & ld <= quantile(ld, 0.95, na.rm = TRUE)
  }
  cd45_idx <- map$channel_index[map$marker == "CD45"]
  if (length(cd45_idx) && !is.na(cd45_idx[1])) {
    cd <- asinh(pmax(as.numeric(exprs[, cd45_idx[1]]), 0) / 150)
    keep <- keep & cd >= qc_cd45_cut(cd)
  }
  keep[is.na(keep)] <- FALSE
  keep
}

extract_marker_mat <- function(exprs, map, markers) {
  idx <- match(markers, map$marker)
  ch <- map$channel_index[idx]
  ok <- !is.na(ch)
  if (!any(ok)) {
    stop("没有匹配到任何分析通道（match_channels 结果全空）")
  }
  mat <- exprs[, ch[ok], drop = FALSE]
  colnames(mat) <- markers[ok]
  mat
}

# 演示数据：可区分的免疫亚群，H 改变部分亚群比例（仅用于缺 FCS 时出图）
rnorm_pop <- function(n, means, sds) {
  m <- matrix(NA_real_, n, length(means))
  colnames(m) <- names(means)
  for (nm in names(means)) {
    m[, nm] <- rnorm(n, means[[nm]], sds[[nm]])
  }
  m
}

demo_means_p1 <- function() {
  mk <- dr_marker_names("P1")
  base <- setNames(rep(0.2, length(mk)), mk)
  pop <- function(...) {
    v <- base
    alt <- list(...)
    for (nm in names(alt)) v[[nm]] <- alt[[nm]]
    v
  }
  list(
    CD4_naive = pop(CD3 = 3.2, CD4 = 3.0, CD8 = 0.1, CD8b = 0.1, CD19 = 0.1, CD62L = 3.2, CD44 = 0.3, CD27 = 2.2),
    CD4_TCM = pop(CD3 = 3.2, CD4 = 3.0, CD62L = 2.8, CD44 = 2.8, CD27 = 2.4, CD95 = 0.3),
    CD4_TSCM = pop(CD3 = 3.2, CD4 = 3.0, CD62L = 2.8, CD44 = 2.8, CD27 = 3.1, CD95 = 3.0),
    CD4_TEM = pop(CD3 = 3.2, CD4 = 3.0, CD62L = 0.3, CD44 = 3.1, CD27 = 0.3, `SCA-1` = 0.3),
    CD4_MPEC = pop(CD3 = 3.2, CD4 = 3.0, CD62L = 0.3, CD44 = 3.1, CD27 = 3.0, `SCA-1` = 0.3),
    CD4_SLEC = pop(CD3 = 3.2, CD4 = 3.0, CD62L = 0.3, CD44 = 3.1, `SCA-1` = 3.1, CD27 = 0.4),
    CD4_exh = pop(CD3 = 3.1, CD4 = 3.0, CD62L = 0.3, CD44 = 3.0, `LAG-3` = 3.1, `TIM-3` = 2.9, `PD-L1` = 2.8, CD27 = 0.3),
    Treg = pop(CD3 = 3.1, CD4 = 3.0, CD25 = 3.2, CD69 = 0.4, CD44 = 1.8),
    CD4_act = pop(CD3 = 3.2, CD4 = 3.0, CD25 = 2.8, CD69 = 3.1, CD44 = 0.3, CD62L = 3.2),
    CD4_eff = pop(CD3 = 3.2, CD4 = 3.0, CD8 = 0.1, CD62L = 0.3, CD44 = 3.0, `IFN-g` = 3.2, `TNF-a` = 3.0, CD69 = 0.4),
    CD8_naive = pop(CD3 = 3.3, CD8 = 3.1, CD8b = 2.9, CD4 = 0.1, CD62L = 3.1, CD44 = 0.3),
    CD8_act = pop(CD3 = 3.3, CD8 = 3.1, CD8b = 2.8, CD4 = 0.1, CD62L = 3.1, CD44 = 0.3, CD69 = 3.1, CD25 = 2.8),
    CD8_TCM = pop(CD3 = 3.3, CD8 = 3.1, CD8b = 2.8, CD62L = 2.7, CD44 = 2.6, CD27 = 2.2, CD95 = 0.3),
    CD8_TSCM = pop(CD3 = 3.3, CD8 = 3.1, CD8b = 2.8, CD62L = 2.7, CD44 = 2.6, CD27 = 3.1, CD95 = 3.0),
    CD8_TEM = pop(CD3 = 3.3, CD8 = 3.1, CD8b = 2.8, CD44 = 3.0, CD62L = 0.3, CD27 = 0.3, `SCA-1` = 0.3),
    CD8_MPEC = pop(CD3 = 3.3, CD8 = 3.1, CD8b = 2.8, CD44 = 3.0, CD62L = 0.3, CD27 = 3.0, `SCA-1` = 0.3),
    CD8_SLEC = pop(CD3 = 3.3, CD8 = 3.1, CD8b = 2.8, CD44 = 3.0, CD62L = 0.3, `SCA-1` = 3.1, CD27 = 0.4),
    CD8_eff = pop(CD3 = 3.3, CD8 = 3.1, CD8b = 2.8, CD44 = 2.6, GZMB = 3.0, Perforin = 2.8, `IFN-g` = 2.2),
    CD8_exh = pop(CD3 = 3.0, CD8 = 3.0, CD8b = 2.6, `LAG-3` = 3.1, `TIM-3` = 2.9, `PD-L1` = 2.8, CD44 = 2.8, CD62L = 0.3),
    NK = pop(CD3 = 0.1, `NK1.1` = 3.2, NKp46 = 3.0, NKG2D = 2.6, CD8 = 1.0, CD27 = 0.3, CD11B = 0.3),
    NK_immature = pop(CD3 = 0.1, `NK1.1` = 3.2, NKp46 = 3.0, CD27 = 3.0, CD11B = 0.3, NKG2D = 2.6, CD69 = 0.3),
    NK_DP = pop(CD3 = 0.1, `NK1.1` = 3.2, NKp46 = 3.0, CD27 = 3.0, CD11B = 2.8, NKG2D = 2.6, CD69 = 0.3),
    NK_mature = pop(CD3 = 0.1, `NK1.1` = 3.2, NKp46 = 3.0, CD27 = 0.3, CD11B = 3.0, NKG2D = 2.6, CD69 = 0.3),
    NK_act = pop(CD3 = 0.1, `NK1.1` = 3.2, NKp46 = 3.0, CD69 = 3.1, CD25 = 2.4, CD11B = 2.2, CD27 = 0.5),
    NK_eff = pop(CD3 = 0.1, `NK1.1` = 3.2, NKp46 = 3.0, GZMB = 3.1, Perforin = 2.9, CD69 = 0.4, CD11B = 2.4, CD27 = 0.4),
    NK_exh = pop(CD3 = 0.1, `NK1.1` = 3.2, NKp46 = 3.0, `LAG-3` = 3.0, `TIM-3` = 2.8, `PD-L1` = 2.8,
                CD11B = 2.6, CD27 = 0.4, CD69 = 0.3),
    NKT_CD4 = pop(CD3 = 3.0, `NK1.1` = 2.8, NKp46 = 2.6, CD4 = 3.0, CD8 = 0.2, CD44 = 2.4, CD69 = 0.3),
    NKT_DN = pop(CD3 = 3.0, `NK1.1` = 2.8, NKp46 = 2.6, CD4 = 0.3, CD8 = 0.3, CD44 = 2.4, CD69 = 0.3),
    B = pop(CD19 = 3.3, CD3 = 0.1, CD27 = 1.2),
    Myeloid = pop(CD11B = 3.2, CD3 = 0.1, CD19 = 0.1, `NK1.1` = 0.1, NKp46 = 0.1)
  )
}

demo_means_p2 <- function() {
  mk <- dr_marker_names("P2")
  base <- setNames(rep(0.2, length(mk)), mk)
  pop <- function(...) {
    v <- base
    alt <- list(...)
    for (nm in names(alt)) v[[nm]] <- alt[[nm]]
    v
  }
  list(
    Naive_B = pop(CD19 = 3.2, IgD = 3.0, IgM = 2.4, CD27 = 0.3),
    Unswitched_B = pop(CD19 = 3.1, IgD = 3.0, CD27 = 3.0, IgM = 2.6, IgG = 0.2),
    MZ_B = pop(CD19 = 3.1, IgD = 3.0, CD27 = 0.3, IgM = 3.3),
    Switched_B = pop(CD19 = 3.0, IgG = 3.1, CD27 = 2.6, IgD = 0.2, IgM = 0.3),
    Plasmablast = pop(CD19 = 2.6, IgD = 0.2, CD27 = 2.8, `BLIMP-1` = 3.2, IgG = 0.4),
    Plasma = pop(CD19 = 0.3, `BLIMP-1` = 3.3, CD27 = 3.1, IgD = 0.2),
    Activated_B = pop(CD19 = 3.0, CD80 = 2.8, CD86 = 3.0, CD40 = 2.6, CD27 = 2.4, IgD = 0.3, IgG = 0.3),
    Atypical_B = pop(CD19 = 3.0, IgD = 0.3, CD27 = 0.3, IgM = 1.6)
  )
}

demo_means_p3 <- function() {
  mk <- dr_marker_names("P3")
  base <- setNames(rep(0.2, length(mk)), mk)
  pop <- function(...) {
    v <- base
    alt <- list(...)
    for (nm in names(alt)) v[[nm]] <- alt[[nm]]
    v
  }
  list(
    Neutrophil = pop(CD11B = 3.2, LY6G = 3.3, LY6C = 1.5, CD3 = 0.1, CD19 = 0.1),
    Mono_Ly6Chi = pop(CD11B = 3.1, LY6C = 3.2, LY6G = 0.2, CD86 = 1.8, `F4/80` = 0.5),
    Mono_Ly6Clo = pop(CD11B = 3.0, LY6C = 0.3, LY6G = 0.2, `F4/80` = 0.4, CD86 = 1.2),
    Macrophage = pop(CD11B = 3.0, `F4/80` = 3.2, `I-A/I-E` = 1.8, CD86 = 1.6),
    M1_like = pop(CD11B = 3.0, `F4/80` = 3.0, iNOS = 3.2, `ARG-1` = 0.2, CD86 = 2.6, CD80 = 2.2),
    M2_like = pop(CD11B = 3.0, `F4/80` = 3.0, CD206 = 3.1, `ARG-1` = 2.8, iNOS = 0.2, `IL-10` = 2.2, `TGF-b` = 2.0),
    DC = pop(CD11C = 3.2, `I-A/I-E` = 3.3, CD80 = 2.2, CD86 = 2.4, CD11B = 0.4, CD103 = 0.2),
    cDC1 = pop(CD11C = 3.1, CD103 = 3.0, `I-A/I-E` = 3.0, CD11B = 0.4),
    cDC2 = pop(CD11C = 3.1, `I-A/I-E` = 3.1, CD11B = 0.4, CD103 = 0.2),
    Eosinophil = pop(`Siglec-F` = 3.2, CCR3 = 2.8, CD11B = 2.6, LY6G = 0.2),
    Mast = pop(FceRI = 3.1, CD200R3 = 2.8, CD11B = 2.8)
  )
}

demo_props <- function(panel_id, group) {
  ctrl <- identical(as.character(group), flow_ctrl_group)
  if (panel_id == "P1") {
    if (ctrl) {
      return(c(CD4_naive = 0.07, CD4_TCM = 0.04, CD4_TSCM = 0.03, CD4_TEM = 0.03, CD4_MPEC = 0.03,
               CD4_SLEC = 0.03, CD4_exh = 0.03, Treg = 0.03, CD4_act = 0.04, CD4_eff = 0.03,
               CD8_naive = 0.05, CD8_act = 0.03, CD8_TCM = 0.03, CD8_TSCM = 0.03, CD8_TEM = 0.03,
               CD8_MPEC = 0.03, CD8_SLEC = 0.03, CD8_eff = 0.04, CD8_exh = 0.04,
               NK = 0.03, NK_immature = 0.03, NK_DP = 0.02, NK_mature = 0.03,
               NK_act = 0.02, NK_eff = 0.02, NK_exh = 0.02, NKT_CD4 = 0.02, NKT_DN = 0.02,
               B = 0.06, Myeloid = 0.04))
    }
    return(c(CD4_naive = 0.03, CD4_TCM = 0.03, CD4_TSCM = 0.03, CD4_TEM = 0.03, CD4_MPEC = 0.03,
             CD4_SLEC = 0.04, CD4_exh = 0.04, Treg = 0.03, CD4_act = 0.05, CD4_eff = 0.04,
             CD8_naive = 0.03, CD8_act = 0.04, CD8_TCM = 0.03, CD8_TSCM = 0.03, CD8_TEM = 0.03,
             CD8_MPEC = 0.03, CD8_SLEC = 0.04, CD8_eff = 0.05, CD8_exh = 0.05,
             NK = 0.03, NK_immature = 0.03, NK_DP = 0.02, NK_mature = 0.03,
             NK_act = 0.03, NK_eff = 0.03, NK_exh = 0.03, NKT_CD4 = 0.02, NKT_DN = 0.03,
             B = 0.05, Myeloid = 0.03))
  }
  if (panel_id == "P2") {
    if (ctrl) {
      return(c(Naive_B = 0.30, Unswitched_B = 0.16, MZ_B = 0.08, Switched_B = 0.14,
               Plasmablast = 0.05, Plasma = 0.07, Activated_B = 0.08, Atypical_B = 0.12))
    }
    return(c(Naive_B = 0.16, Unswitched_B = 0.12, MZ_B = 0.08, Switched_B = 0.20,
             Plasmablast = 0.10, Plasma = 0.14, Activated_B = 0.10, Atypical_B = 0.10))
  }
  if (ctrl) {
    return(c(Neutrophil = 0.16, Mono_Ly6Chi = 0.10, Mono_Ly6Clo = 0.08, Macrophage = 0.12,
             M1_like = 0.06, M2_like = 0.06, DC = 0.08, cDC1 = 0.08, cDC2 = 0.08,
             Eosinophil = 0.10, Mast = 0.08))
  }
  c(Neutrophil = 0.11, Mono_Ly6Chi = 0.07, Mono_Ly6Clo = 0.06, Macrophage = 0.09,
    M1_like = 0.06, M2_like = 0.13, DC = 0.10, cDC1 = 0.09, cDC2 = 0.09,
    Eosinophil = 0.11, Mast = 0.09)
}

make_demo_sample <- function(panel_id, group, sample, n) {
  means_fun <- switch(panel_id, P1 = demo_means_p1, P2 = demo_means_p2, P3 = demo_means_p3)
  pops <- means_fun()
  props <- demo_props(panel_id, group)
  props <- props[names(pops)]
  props <- props / sum(props)
  counts <- pmax(1, as.integer(round(props * n)))
  extra <- n - sum(counts)
  if (extra != 0) counts[1] <- counts[1] + extra
  sds <- lapply(pops, function(v) setNames(rep(0.35, length(v)), names(v)))
  parts <- mapply(function(mu, sd, nn, lab) {
    m <- rnorm_pop(nn, as.list(mu), as.list(sd))
    data.frame(m, true_lineage = lab, stringsAsFactors = FALSE, check.names = FALSE)
  }, pops, sds, counts, names(pops), SIMPLIFY = FALSE)
  df <- do.call(rbind, parts)
  cd45 <- matrix(rnorm(nrow(df), 3.0, 0.25), ncol = 1, dimnames = list(NULL, "CD45"))
  ld <- matrix(rnorm(nrow(df), 0.2, 0.15), ncol = 1, dimnames = list(NULL, "L/D"))
  expr <- cbind(cd45, ld, as.matrix(df[, setdiff(colnames(df), "true_lineage"), drop = FALSE]))
  list(expr = expr, true_lineage = df$true_lineage, sample = sample, group = group)
}

# -----------------------------------------------------------------------------
# 6. 降维 / 聚类
# -----------------------------------------------------------------------------
run_pca <- function(mat, npcs = 20) {
  npcs <- min(npcs, ncol(mat), nrow(mat) - 1)
  pr <- prcomp(mat, center = TRUE, scale. = FALSE)
  pcs <- pr$x[, seq_len(npcs), drop = FALSE]
  list(embedding = pcs, sdev = pr$sdev)
}

run_umap <- function(mat, compact = NULL) {
  n <- nrow(as.matrix(mat))
  if (is.null(compact)) compact <- n < 2500L
  if (isTRUE(compact)) {
    nn <- max(5L, min(12L, n - 1L))
    min_dist <- 0.06
    spread <- 0.7
  } else {
    nn <- max(10L, min(20L, n - 1L))
    min_dist <- 0.18
    spread <- 1
  }
  if (has_pkg("uwot") && n >= 5L) {
    emb <- uwot::umap(
      mat, n_neighbors = nn, min_dist = min_dist, spread = spread,
      metric = "euclidean", verbose = FALSE, n_threads = 2, seed = seed_value
    )
    colnames(emb) <- c("UMAP1", "UMAP2")
    return(emb)
  }
  log_msg("uwot 未安装，UMAP 用前两主成分代替（图标题会标明 PCA fallback）")
  pca <- run_pca(mat, 2)$embedding
  colnames(pca) <- c("UMAP1", "UMAP2")
  pca
}

run_tsne <- function(mat, compact = NULL) {
  n <- nrow(as.matrix(mat))
  if (is.null(compact)) compact <- n < 2500L
  if (isTRUE(compact)) {
    perplexity <- max(5, min(12, floor((n - 1) / 6)))
  } else {
    perplexity <- max(5, min(30, floor((n - 1) / 3)))
  }
  if (has_pkg("Rtsne")) {
    set.seed(seed_value)
    ts <- Rtsne::Rtsne(
      mat, perplexity = perplexity, verbose = FALSE, check_duplicates = FALSE,
      pca = TRUE, max_iter = 1000
    )
    emb <- ts$Y
    colnames(emb) <- c("tSNE1", "tSNE2")
    return(emb)
  }
  log_msg("Rtsne 未安装，tSNE 用主成分 3–4 或 1–2 代替")
  pca <- run_pca(mat, 4)$embedding
  if (ncol(pca) >= 4) {
    emb <- pca[, 3:4, drop = FALSE]
  } else {
    emb <- pca[, 1:2, drop = FALSE]
  }
  colnames(emb) <- c("tSNE1", "tSNE2")
  emb
}

choose_k <- function(panel_id) {
  switch(panel_id, P1 = 22, P2 = 16, P3 = 20, 12)
}

cluster_cells <- function(mat, panel_id) {
  k <- choose_k(panel_id)
  if (has_pkg("FlowSOM") && has_pkg("flowCore")) {
    ff <- flowCore::flowFrame(mat)
    fs <- tryCatch(
      FlowSOM::FlowSOM(ff, colsToUse = colnames(mat), nClus = k, seed = seed_value),
      error = function(e) NULL
    )
    if (!is.null(fs)) {
      cl <- as.integer(FlowSOM::GetMetaclusters(fs))
      cl[is.na(cl)] <- 0L
      return(factor(paste0("C", cl), levels = paste0("C", sort(unique(cl)))))
    }
  }
  set.seed(seed_value)
  # Hartigan-Wong 在细胞数多时会报 Quick-TRANSfer；Lloyd 更稳，结果仍可用于分群
  km <- kmeans(mat, centers = k, nstart = 5, iter.max = 400, algorithm = "Lloyd")
  factor(paste0("C", km$cluster), levels = paste0("C", sort(unique(km$cluster))))
}

annotate_clusters <- function(med, panel_id) {
  med <- as.matrix(med)
  rawv <- function(i, m) {
    if (!m %in% colnames(med)) return(NA_real_)
    as.numeric(med[i, m])
  }
  nv <- function(i, m) {
    x <- rawv(i, m)
    if (!is.finite(x)) -Inf else x
  }
  labs <- vapply(seq_len(nrow(med)), function(i) {
    cd4 <- nv(i, "CD4")
    cd8 <- max(nv(i, "CD8"), nv(i, "CD8b"))
    cd3 <- nv(i, "CD3")
    cd19 <- nv(i, "CD19")
    cd11b <- nv(i, "CD11B")
    nk <- max(nv(i, "NK1.1"), nv(i, "NKp46"))
    t_mem <- function(prefix) {
      row <- setNames(as.list(med[i, ]), colnames(med))
      cd62 <- nv(i, "CD62L")
      cd44 <- nv(i, "CD44")
      if (!is.finite(cd62)) cd62 <- -Inf
      if (!is.finite(cd44)) cd44 <- -Inf
      mem <- if (cd62 > cd44 + 0.3) {
        paste0(prefix, "_naive")
      } else if (cd44 > cd62 + 0.3) {
        paste0(prefix, "_TEM")
      } else if (is.finite(cd62) || is.finite(cd44)) {
        paste0(prefix, "_TCM")
      } else {
        paste0(prefix, "_T")
      }
      refine_t_memory_label(row, prefix, mem)
    }
    if (panel_id == "P1") {
      nkp <- nv(i, "NKp46")
      nk11 <- nv(i, "NK1.1")
      nk <- if (is.finite(nkp) && nkp > -Inf && nkp >= 0) nkp else max(nkp, nk11)
      row <- setNames(as.list(med[i, ]), colnames(med))
      # CD45+ 内：B = CD19+；NK = CD3- NKp46+（必须在 CD11b 之前，成熟 NK 是 CD11b+）
      if (is.finite(cd19) && cd19 >= max(cd3, cd4, cd8, nk) - 0.05 && cd19 > cd3 + 0.15) return("B")
      if (is.finite(nk) && nk > cd3 + 0.3 && nk > cd19 && nk >= max(cd4, cd8) - 0.15) {
        return(label_nk_subset(row))
      }
      if (is.finite(nk) && is.finite(cd3) && nk >= 1.8 && cd3 >= 1.8 &&
          abs(cd3 - nk) < 1.2 && nk > cd19 && cd3 > cd19) {
        return(label_nkt_subset(row))
      }
      if (cd8 > cd4 + 0.15) {
        kill <- max(nv(i, "GZMB"), nv(i, "Perforin"))
        cyto <- max(nv(i, "IFN-g"), nv(i, "TNF-a"))
        cd69 <- nv(i, "CD69")
        cd62 <- nv(i, "CD62L")
        cd44 <- nv(i, "CD44")
        if (is.finite(kill) && kill > max(cd62, cd44) + 0.15) return("CD8_effector")
        if (is.finite(cyto) && cyto >= 2.5 && cyto >= max(cd62, cd44) - 0.15 && cyto > kill) {
          return("CD8_effector")
        }
        if (is.finite(cd69) && cd69 >= 2.2) return("CD8_activated")
        return(t_mem("CD8"))
      }
      if (cd4 > cd8 + 0.15) {
        cd25 <- nv(i, "CD25")
        cd69 <- nv(i, "CD69")
        tnfa <- nv(i, "TNF-a")
        ifng <- nv(i, "IFN-g")
        cd62 <- nv(i, "CD62L")
        cd44 <- nv(i, "CD44")
        if (is.finite(cd69) && cd69 >= 2.2) return("CD4_activated")
        if (is.finite(cd25) && cd25 > max(cd69, tnfa, ifng) + 0.2 &&
            cd25 >= max(cd62, cd44) - 0.6) return("Treg")
        cyto <- max(ifng, tnfa)
        if (is.finite(cyto) && cyto >= 2.5 && cyto >= max(cd62, cd44) - 0.15) return("CD4_effector")
        return(t_mem("CD4"))
      }
      if (is.finite(cd3) && cd3 > max(cd19, cd11b, nk)) return(t_mem("T"))
      if (is.finite(cd11b) && cd11b > max(cd3, cd4, cd8, cd19) + 0.2) return("Myeloid")
      return("Myeloid")
    }
    if (panel_id == "P2") {
      return(label_b_subset(setNames(vapply(colnames(med), function(nm) nv(i, nm), numeric(1)), colnames(med))))
    }
    # P3：先 CD3/CD19/NK1.1 大类，再按 CD11B 两支。嗜酸要 Siglec-F 且 CCR3；DC 在 CD11B-。
    return(label_p3_cluster(setNames(vapply(colnames(med), function(nm) nv(i, nm), numeric(1)), colnames(med))))
  }, character(1))
  data.frame(cluster = rownames(med), lineage = labs, stringsAsFactors = FALSE)
}

colv <- function(mat, name) {
  if (!name %in% colnames(mat)) return(rep(NA_real_, nrow(mat)))
  as.numeric(mat[, name])
}

# 缺通道（ICI P1 没有 CD19）时当 -Inf，不要留下 NA 把 is_t 弄成 NA 后炸在 out[is_t] <-
colv_neg <- function(mat, name) {
  x <- colv(mat, name)
  x[!is.finite(x)] <- -Inf
  x
}

na_to_false <- function(x) {
  x <- as.logical(x)
  x[is.na(x)] <- FALSE
  x
}

finite_pmax <- function(...) {
  args <- list(...)
  out <- args[[1]]
  out[!is.finite(out)] <- -Inf
  if (length(args) == 1) return(out)
  for (a in args[-1]) {
    a[!is.finite(a)] <- -Inf
    out <- pmax(out, a)
  }
  out
}

vec_get <- function(v, name) {
  if (!name %in% names(v)) return(-Inf)
  x <- as.numeric(v[[name]])
  if (!is.finite(x)) -Inf else x
}

mem_from_cd62_cd44 <- function(cd62, cd44, prefix) {
  if (cd62 > cd44 + 0.15) return(paste0(prefix, "_naive"))
  if (cd44 > cd62 + 0.15) return(paste0(prefix, "_TEM"))
  paste0(prefix, "_TCM")
}

# T_CM 里再圈 T_SCM（CD27+ CD95+）；T_EM 再分耗竭 / SLEC / MPEC / 早晚期
refine_t_memory_label <- function(v, prefix, mem) {
  if (!length(mem) || !nzchar(mem)) return(mem)
  cd27 <- vec_get(v, "CD27")
  cd95 <- vec_get(v, "CD95")
  sca <- vec_get(v, "SCA-1")
  exh <- max(vec_get(v, "LAG-3"), vec_get(v, "TIM-3"), vec_get(v, "PD-L1"))
  if (grepl("_TCM$", mem) && cd27 >= 2.0 && cd95 >= 2.0) {
    return(paste0(prefix, "_TSCM"))
  }
  if (!grepl("_TEM$", mem)) return(mem)
  if (exh >= 2.2) return(paste0(prefix, "_exhausted"))
  if (sca >= 2.2) return(paste0(prefix, "_SLEC"))
  if (cd27 >= 2.0 && sca < 2.2) return(paste0(prefix, "_MPEC"))
  if (cd27 >= 2.0) return(paste0(prefix, "_TEM_early"))
  paste0(prefix, "_TEM_late")
}

tem_family <- function(prefix) {
  paste0(prefix, c("_TEM", "_TEM_early", "_TEM_late", "_SLEC", "_MPEC", "_exhausted"))
}

tcm_family <- function(prefix) {
  paste0(prefix, c("_TCM", "_TSCM"))
}

# P1 NK 亚群（不含 NKT：NKT 是 CD3+）
nk_family <- function() {
  c("NK", "NK_activated", "NK_effector", "NK_exhausted",
    "NK_immature", "NK_DP", "NK_mature")
}

nkt_family <- function() {
  c("NKT", "NKT_CD4", "NKT_DN", "NKT_activated", "NKT_effector")
}

# 降维总图用的大类（P1 的 B/髓系、P3 的 T/B/NK 仍画在总图上，不当 dump 丢掉）
dimred_major_of <- function(panel_id, lineage, cluster_lineage = NULL) {
  s <- as.character(lineage)
  s[is.na(s)] <- ""
  n <- length(s)
  out <- rep(NA_character_, n)
  if (!is.null(cluster_lineage) && length(cluster_lineage) == n) {
    cl <- as.character(cluster_lineage)
    keep <- !is.na(cl) & nzchar(cl) & !cl %in% c("other", "Other")
    out[keep] <- cl[keep]
  }
  miss <- is.na(out) | !nzchar(out)
  if (!any(miss)) return(out)
  inf <- s
  inf[s %in% c("Target", "His_target")] <- "Target"
  inf[s %in% nkt_family() | grepl("^NKT", s)] <- "NKT"
  inf[s %in% nk_family() | grepl("^NK_", s) | s == "NK"] <- "NK"
  inf[s %in% c("Treg", "CD4_activated", "CD4_effector", "CD4") | grepl("^CD4_", s)] <- "CD4"
  inf[s %in% c("CD8") | grepl("^CD8_", s)] <- "CD8"
  inf[s %in% c("B", "Naive_B", "Atypical_B", "Memory_B", "Switched_B",
               "Unswitched_B", "IgM_memory", "MZ_B", "Plasmablast",
               "Activated_B", "Plasma")] <- "B"
  inf[s %in% c("Myeloid", "Macrophage", "M1_like_Mac", "M2_like_Mac",
               "Neutrophil", "Eosinophil", "Mono_Ly6Chi", "Mono_Ly6Clo",
               "DC", "cDC1_CD103", "cDC2", "Mast", "Basophil_mast", "Basophil")] <- "Myeloid"
  inf[s %in% c("T", "T_naive", "T_TCM", "T_TEM", "T_effector", "T_T")] <- "T"
  inf[s == "NKT"] <- "NKT"
  out[miss] <- inf[miss]
  out[is.na(out) | !nzchar(out)] <- "other"
  out
}

major_display_label <- function(major) {
  rec <- c(
    CD4 = "CD4 T", CD8 = "CD8 T", T = "T", NK = "NK", NKT = "NKT",
    B = "B cell", Myeloid = "Myeloid", Target = "His+ target",
    His_target = "His+ target", dump = "other", other = "other", Other = "other"
  )
  lab <- as.character(major)
  out <- rec[lab]
  unname(ifelse(is.na(out), lab, out))
}

major_colors <- function(levels) {
  pal <- pal_major[levels]
  miss <- levels[is.na(pal)]
  if (length(miss) > 0) {
    extra <- pal_p1_hues[((seq_along(miss) - 1L) %% length(pal_p1_hues)) + 1L]
    pal[is.na(pal)] <- extra
  }
  names(pal) <- levels
  pal
}

# P1 NK：不要用 NKG2D 当亚群（小鼠 NK 组成性表达）
label_nk_subset <- function(v) {
  cd69 <- vec_get(v, "CD69")
  kill <- max(vec_get(v, "GZMB"), vec_get(v, "Perforin"))
  exh <- max(vec_get(v, "PD-L1"), vec_get(v, "LAG-3"), vec_get(v, "TIM-3"))
  cd27 <- vec_get(v, "CD27")
  cd11b <- vec_get(v, "CD11B")
  if (cd69 >= 2.2) return("NK_activated")
  if (kill >= 2.5) return("NK_effector")
  if (exh >= 2.5) return("NK_exhausted")
  if (cd27 > 1.2 && cd11b < 1.0) return("NK_immature")
  if (cd27 > 1.2 && cd11b > 1.2) return("NK_DP")
  if (cd11b > 1.2 && cd27 < 1.2) return("NK_mature")
  "NK"
}

label_nkt_subset <- function(v) {
  cd69 <- vec_get(v, "CD69")
  kill <- max(vec_get(v, "GZMB"), vec_get(v, "IFN-g"))
  cd4 <- vec_get(v, "CD4")
  if (cd69 >= 2.2) return("NKT_activated")
  if (kill >= 2.5) return("NKT_effector")
  if (cd4 >= 1.5) return("NKT_CD4")
  "NKT_DN"
}

# 大类圈定之后，只根据这一层的定义标志给亚群命名
label_cd4_subset <- function(v) {
  cd25 <- vec_get(v, "CD25")
  cd69 <- vec_get(v, "CD69")
  tnfa <- vec_get(v, "TNF-a")
  ifng <- vec_get(v, "IFN-g")
  cd62 <- vec_get(v, "CD62L")
  cd44 <- vec_get(v, "CD44")
  # 早期活化：CD69+（常 CD25 高），可以同时是 CD62L+ CD44-
  if (cd69 >= 2.2) return("CD4_activated")
  if (cd25 > max(cd69, tnfa, ifng) + 0.15 && cd25 > max(cd62, cd44) - 0.8) return("Treg")
  cyto <- max(ifng, tnfa)
  if (cyto >= 2.5 && cyto >= max(cd62, cd44) - 0.15) return("CD4_effector")
  refine_t_memory_label(v, "CD4", mem_from_cd62_cd44(cd62, cd44, "CD4"))
}

label_cd8_subset <- function(v) {
  kill <- max(vec_get(v, "GZMB"), vec_get(v, "Perforin"))
  cyto <- max(vec_get(v, "IFN-g"), vec_get(v, "TNF-a"))
  cd69 <- vec_get(v, "CD69")
  cd62 <- vec_get(v, "CD62L")
  cd44 <- vec_get(v, "CD44")
  if (kill > max(cd62, cd44) + 0.1) return("CD8_effector")
  if (cyto >= 2.5 && cyto >= max(cd62, cd44) - 0.15) return("CD8_effector")
  if (cd69 >= 2.2) return("CD8_activated")
  refine_t_memory_label(v, "CD8", mem_from_cd62_cd44(cd62, cd44, "CD8"))
}

label_t_subset <- function(v) {
  mem_from_cd62_cd44(vec_get(v, "CD62L"), vec_get(v, "CD44"), "T")
}

label_b_subset <- function(v) {
  cd19 <- vec_get(v, "CD19")
  blimp <- vec_get(v, "BLIMP-1")
  igd <- vec_get(v, "IgD")
  igm <- vec_get(v, "IgM")
  igg <- vec_get(v, "IgG")
  cd27 <- vec_get(v, "CD27")
  act <- max(vec_get(v, "CD80"), vec_get(v, "CD86"))
  igd_hi <- igd >= 1.5
  cd27_pos <- cd27 >= 1.5
  # Plasma 只从 CD19- / dim：不要用 BLIMP 把 CD19+ 整团打成浆细胞
  if (cd19 < 1.0 && blimp >= 2.2 && cd27 >= 2.0 && !igd_hi) return("Plasma")
  if (igd_hi && !cd27_pos) {
    if (igm >= 2.8 && igm >= igd - 0.2) return("MZ_B")
    if (act >= 2.4 && act > igd) return("Activated_B")
    return("Naive_B")
  }
  if (igd_hi && cd27_pos) {
    if (igm >= 2.8 && igm >= igd - 0.2) return("MZ_B")
    if (act >= 2.4 && act > igd) return("Activated_B")
    return("Unswitched_B")
  }
  if (!igd_hi && cd27_pos) {
    if (blimp >= 2.5 && cd19 >= 1.0) return("Plasmablast")
    if (act >= 2.4 && act > max(igd, igm, igg)) return("Activated_B")
    return("Switched_B")
  }
  "Atypical_B"
}

# P2 第 1 层：IgD vs CD27 四象限。不要把 BLIMP-1 放进这一层。
label_b_major <- function(v) {
  cd19 <- vec_get(v, "CD19")
  blimp <- vec_get(v, "BLIMP-1")
  igd <- vec_get(v, "IgD")
  cd27 <- vec_get(v, "CD27")
  if (cd19 < 1.0 && blimp >= 2.2 && cd27 >= 2.0 && igd < 1.2) return("Plasma")
  igd_hi <- igd >= 1.5
  cd27_pos <- cd27 >= 1.5
  if (igd_hi && !cd27_pos) return("Naive_B")
  if (igd_hi && cd27_pos) return("Unswitched_B")
  if (!igd_hi && cd27_pos) return("Switched_B")
  "Atypical_B"
}

label_b_naive_subset <- function(v) {
  igd <- vec_get(v, "IgD")
  igm <- vec_get(v, "IgM")
  act <- max(vec_get(v, "CD80"), vec_get(v, "CD86"))
  if (igm >= 2.8 && igm >= igd - 0.2) return("MZ_B")
  if (act > igd + 0.25 && act >= 2.4) return("Activated_B")
  "Naive_B"
}

label_b_unswitched_subset <- function(v) {
  igd <- vec_get(v, "IgD")
  igm <- vec_get(v, "IgM")
  act <- max(vec_get(v, "CD80"), vec_get(v, "CD86"))
  if (igm >= 2.8 && igm >= igd - 0.2) return("MZ_B")
  if (act > igd + 0.25 && act >= 2.4) return("Activated_B")
  "Unswitched_B"
}

label_b_switched_subset <- function(v) {
  igd <- vec_get(v, "IgD")
  blimp <- vec_get(v, "BLIMP-1")
  act <- max(vec_get(v, "CD80"), vec_get(v, "CD86"))
  if (blimp >= 2.5 && blimp > igd + 0.3) return("Plasmablast")
  if (act >= 2.4 && act > igd + 0.25) return("Activated_B")
  "Switched_B"
}

# 旧 Memory 大类后备
label_b_memory_subset <- function(v) {
  label_b_switched_subset(v)
}

label_p3_cluster <- function(v) {
  cd3 <- vec_get(v, "CD3")
  cd19 <- vec_get(v, "CD19")
  nk <- vec_get(v, "NK1.1")
  lymph <- max(cd3, cd19, nk)
  siglec <- vec_get(v, "Siglec-F")
  ccr3 <- vec_get(v, "CCR3")
  ly6g <- vec_get(v, "LY6G")
  f480 <- vec_get(v, "F4/80")
  cd11c <- vec_get(v, "CD11C")
  mhc <- vec_get(v, "I-A/I-E")
  cd103 <- vec_get(v, "CD103")
  ly6c <- vec_get(v, "LY6C")
  cd11b <- vec_get(v, "CD11B")
  inos <- vec_get(v, "iNOS")
  arg1 <- vec_get(v, "ARG-1")
  cd206 <- vec_get(v, "CD206")
  fceri <- vec_get(v, "FceRI")
  cd200 <- vec_get(v, "CD200R3")
  myel_def <- max(siglec, ly6g, f480, cd11c, ly6c, inos, arg1, cd206, fceri, cd200, cd103, mhc)
  if (cd19 >= lymph && cd19 > myel_def + 0.15) return("B")
  if (cd3 >= lymph && cd3 > myel_def + 0.15) return("T")
  if (nk >= lymph && nk > myel_def + 0.15) return("NK")
  label_myeloid_subset(v)
}

# P3 髓系命名：CD11B+ 先粒/肥大/巨噬，CD11B- 才是 DC。MHCII+ 巨噬仍是巨噬。
label_myeloid_major <- function(v) {
  lab <- label_myeloid_subset(v)
  if (lab %in% c("M1_like_Mac", "M2_like_Mac")) return("Macrophage")
  if (lab %in% c("cDC1_CD103", "cDC2")) return("DC")
  lab
}

label_mac_cytokine <- function(v) {
  inos <- vec_get(v, "iNOS")
  arg1 <- vec_get(v, "ARG-1")
  cd206 <- vec_get(v, "CD206")
  if (inos >= 2.2 && arg1 < 1.5) return("M1_like_Mac")
  if (arg1 >= 2.2 && cd206 >= 2.0 && inos < 1.5) return("M2_like_Mac")
  "Macrophage"
}

label_dc_cytokine <- function(v) {
  if (vec_get(v, "CD103") >= 2.0) return("cDC1_CD103")
  "cDC2"
}

label_myeloid_subset <- function(v) {
  siglec <- vec_get(v, "Siglec-F")
  ccr3 <- vec_get(v, "CCR3")
  ly6g <- vec_get(v, "LY6G")
  f480 <- vec_get(v, "F4/80")
  cd11c <- vec_get(v, "CD11C")
  mhc <- vec_get(v, "I-A/I-E")
  cd103 <- vec_get(v, "CD103")
  ly6c <- vec_get(v, "LY6C")
  cd11b <- vec_get(v, "CD11B")
  inos <- vec_get(v, "iNOS")
  arg1 <- vec_get(v, "ARG-1")
  cd206 <- vec_get(v, "CD206")
  fceri <- vec_get(v, "FceRI")
  cd200 <- vec_get(v, "CD200R3")
  if (ly6g >= 1.2 && ly6g > siglec && ly6g > fceri) return("Neutrophil")
  if (siglec >= 1.2 && ccr3 >= 1.2 && siglec > ly6g) return("Eosinophil")
  if (fceri >= 1.2 && cd200 >= 1.2 && fceri > ly6g && fceri > siglec) return("Mast")
  if (f480 >= 2.0 && f480 > cd11c && cd11b >= 1.2) {
    if (inos >= 2.2 && arg1 < 1.5) return("M1_like_Mac")
    if (arg1 >= 2.2 && cd206 >= 2.0 && inos < 1.5) return("M2_like_Mac")
    return("Macrophage")
  }
  if (cd11c >= 1.2 && mhc >= 1.2 && f480 < cd11c && cd11b < 1.8) {
    if (cd103 >= 2.0) return("cDC1_CD103")
    return("cDC2")
  }
  if (ly6c >= 1.2 && cd11b >= 1.2 && f480 < 2.0 && ly6c > ly6g) return("Mono_Ly6Chi")
  if (cd11b >= 1.2) return("Mono_Ly6Clo")
  "Myeloid"
}

# P2 第 1 层：CD19+ 上 IgD vs CD27 四象限
#   Naive = IgD high CD27-；Unswitched = IgD high CD27+；Switched = IgD- CD27+；DN = Atypical
# 不要把 BLIMP-1 放进这一层（核背景会吞掉亚群）。CD19- 里再圈 Plasma。
p2_pos_mask <- function(v, min_sep = 0.45) {
  v <- as.numeric(v)
  n <- length(v)
  if (!n) return(logical(0))
  cut <- axis_pos_cut(v, min_sep = min_sep)
  is.finite(v) & is.finite(cut) & v >= cut
}

p2_assign_cd19neg_plasma <- function(mat, out) {
  io <- which(out == "other")
  if (length(io) < 8 || !"BLIMP-1" %in% colnames(mat)) return(out)
  bl <- colv(mat, "BLIMP-1")
  igd <- colv(mat, "IgD")
  cd27 <- colv(mat, "CD27")
  # 阈值按全体细胞切，不要只在 CD19- 内部切（整团都是浆细胞时会把门切到云团上方）
  bl_hi <- p2_pos_mask(bl, 0.4)
  cd27_hi <- p2_pos_mask(cd27, 0.4)
  igd_hi <- p2_pos_mask(igd, 0.4)
  hit <- bl_hi[io] & cd27_hi[io] & !igd_hi[io]
  if (!any(hit)) {
    hit <- bl[io] >= 2.0 & cd27[io] >= 2.0 & igd[io] < 1.2
  }
  if (any(hit)) out[io[hit]] <- "Plasma"
  out
}

gate_p2_major <- function(mat) {
  mat <- as.matrix(mat)
  n <- nrow(mat)
  out <- rep("Naive_B", n)
  if (n < 40) return(out)
  b <- rep(TRUE, n)
  if ("CD19" %in% colnames(mat)) {
    b_hi <- gate_k2_high(mat, seq_len(n), "CD19", 0.12)
    if (any(b_hi) && mean(b_hi) > 0.08 && mean(b_hi) < 0.97) {
      lo_med <- median(colv(mat, "CD19")[!b_hi], na.rm = TRUE)
      hi_med <- median(colv(mat, "CD19")[b_hi], na.rm = TRUE)
      if (is.finite(lo_med) && is.finite(hi_med) && hi_med > lo_med + 1.0 && lo_med < 0.9) {
        b <- b_hi
        out[!b] <- "other"
      }
    }
  }
  ib <- which(b)
  if (length(ib) >= 40 && all(c("IgD", "CD27") %in% colnames(mat))) {
    igd <- colv(mat, "IgD")[ib]
    cd27 <- colv(mat, "CD27")[ib]
    igd_hi <- p2_pos_mask(igd)
    cd27_pos <- p2_pos_mask(cd27)
    igd_hi[!is.finite(igd)] <- FALSE
    cd27_pos[!is.finite(cd27)] <- FALSE
    lab <- rep("Atypical_B", length(ib))
    lab[igd_hi & !cd27_pos] <- "Naive_B"
    lab[igd_hi & cd27_pos] <- "Unswitched_B"
    lab[!igd_hi & cd27_pos] <- "Switched_B"
    out[ib] <- lab
  } else if (length(ib) >= 40 && "CD27" %in% colnames(mat)) {
    mem <- p2_pos_mask(colv(mat, "CD27")[ib])
    out[ib[mem]] <- "Unswitched_B"
    out[ib[!mem]] <- "Naive_B"
  }
  p2_assign_cd19neg_plasma(mat, out)
}

# P1 NK：文档是 CD3- NKp46+；NKp46 缺失时才退回 NK1.1
p1_nk_score <- function(mat) {
  nkp <- colv(mat, "NKp46")
  nk11 <- colv(mat, "NK1.1")
  out <- nkp
  miss <- !is.finite(out)
  out[miss] <- nk11[miss]
  out[!is.finite(out)] <- -Inf
  out
}

# P3 第 1 层：CD45+ 里 CD3+ T、CD3- CD19+ B、CD3- CD19- NK1.1+ NK，其余三阴是髓系
gate_p3_major <- function(mat) {
  mat <- as.matrix(mat)
  n <- nrow(mat)
  cd3_pos <- p2_pos_mask(colv(mat, "CD3"))
  cd19_pos <- p2_pos_mask(colv(mat, "CD19"))
  nk_pos <- p2_pos_mask(colv(mat, "NK1.1"))
  out <- rep("Myeloid", n)
  out[cd3_pos] <- "T"
  out[!cd3_pos & cd19_pos] <- "B"
  out[!cd3_pos & !cd19_pos & nk_pos] <- "NK"
  out
}

# 第 1 层：只用谱系抗体圈大类（P1 T/NK，P2 B 的 IgD×CD27，P3 CD3/CD19/NK1.1）
gate_major_lineage <- function(mat, panel_id) {
  mat <- as.matrix(mat)
  n <- nrow(mat)
  if (panel_id == "P2") {
    return(gate_p2_major(mat))
  }
  if (panel_id == "P3") {
    return(gate_p3_major(mat))
  }
  cd3 <- colv_neg(mat, "CD3")
  cd4 <- colv_neg(mat, "CD4")
  cd8 <- finite_pmax(colv(mat, "CD8"), colv(mat, "CD8b"))
  cd19 <- colv_neg(mat, "CD19")
  cd11b <- colv_neg(mat, "CD11B")
  nk <- if (identical(panel_id, "P1")) p1_nk_score(mat) else finite_pmax(colv(mat, "NK1.1"), colv(mat, "NKp46"))
  nk[!is.finite(nk)] <- -Inf
  out <- rep("Myeloid", n)
  # B：CD45+ 里 CD19+。ICI P1 没有 CD19 时整列 -Inf，is_b 全 FALSE，不要留下 NA
  is_b <- na_to_false(cd19 > finite_pmax(cd3, cd4, cd8, nk) + 0.1)
  # NK：CD3- NKp46+。必须在 CD11b 之前，成熟 NK 是 CD11b+ NKp46+
  is_nk <- na_to_false(!is_b & nk > cd3 + 0.3 & nk > cd19)
  is_nkt <- na_to_false(!is_b & !is_nk & nk > 1 & cd3 > 1 & abs(cd3 - nk) < 1.2 & nk > cd19)
  # T：CD3+，再按 CD4 vs CD8 象限：CD4+ CD8- / CD4- CD8+（双阳/双阴不塞进某一边）
  is_t <- na_to_false(!is_b & !is_nk & !is_nkt & cd3 >= finite_pmax(cd19, cd11b))
  # 剩下的 CD11b+ 才是髓系 dump
  is_my <- na_to_false(!is_b & !is_nk & !is_nkt & !is_t & cd11b > finite_pmax(cd3, cd19, cd4, cd8) + 0.2)
  out[is_b] <- "B"
  out[is_nk] <- "NK"
  out[is_nkt] <- "NKT"
  t_lab <- assign_t_cd4_cd8(cd4, cd8, is_t)
  out[is_t] <- t_lab[is_t]
  out[is_t & (is.na(out) | !nzchar(out))] <- "other"
  out[is_my] <- "Myeloid"
  out
}

# 第 2 层：大类内部只用亚群标志做 kmeans，再用中位数命名（同一亚群一块颜色，不撒点）
split_within_parent <- function(mat, idx, markers, k, label_fun, fallback) {
  n <- length(idx)
  if (n == 0) return(character(0))
  markers <- intersect(markers, colnames(mat))
  if (n < 25 || length(markers) < 2) {
    med <- apply(mat[idx, , drop = FALSE], 2, median, na.rm = TRUE)
    return(rep(label_fun(med), n))
  }
  k_use <- max(2L, min(as.integer(k), max(2L, floor(n / 80L))))
  x <- mat[idx, markers, drop = FALSE]
  xs <- scale(x)
  xs[!is.finite(xs)] <- 0
  set.seed(seed_value)
  km <- tryCatch(
    kmeans(xs, centers = k_use, nstart = 8, iter.max = 250, algorithm = "Lloyd"),
    error = function(e) NULL
  )
  if (is.null(km)) return(rep(fallback, n))
  labs <- rep(fallback, n)
  for (ci in sort(unique(km$cluster))) {
    hit <- km$cluster == ci
    med <- apply(x[hit, , drop = FALSE], 2, median, na.rm = TRUE)
    labs[hit] <- label_fun(med)
  }
  labs[is.na(labs) | !nzchar(labs)] <- fallback
  labs
}

marker_score <- function(mat, idx, markers) {
  markers <- intersect(markers, colnames(mat))
  if (length(idx) == 0) return(numeric(0))
  if (length(markers) == 0) return(rep(NA_real_, length(idx)))
  x <- mat[idx, markers, drop = FALSE]
  if (ncol(x) == 1) return(as.numeric(x[, 1]))
  apply(x, 1, function(r) {
    r <- r[is.finite(r)]
    if (!length(r)) return(NA_real_)
    max(r)
  })
}

# 一维阳性阈值：优先抓住主峰右侧的尾巴（naive CD62L、CD44），
# 不要把原点一团 50/50 切开，也不要用全体 q80 在主团里切一条缝冒充阳性门。
axis_pos_cut <- function(v, min_sep = 0.55, min_pos_frac = 0.02, max_pos_frac = 0.92) {
  v <- as.numeric(v)
  v <- v[is.finite(v)]
  if (!length(v)) return(0)
  if (length(v) < 20) return(as.numeric(stats::median(v)))
  accept <- function(cut) {
    if (!is.finite(cut)) return(FALSE)
    fr <- mean(v >= cut)
    is.finite(fr) && fr >= min_pos_frac && fr <= max_pos_frac
  }
  xs <- scale(as.matrix(v))
  xs[!is.finite(xs)] <- 0
  set.seed(if (exists("seed_value")) seed_value else 42)
  km <- tryCatch(
    stats::kmeans(xs, centers = 2, nstart = 10, iter.max = 200, algorithm = "Lloyd"),
    error = function(e) NULL
  )
  if (!is.null(km) && length(unique(km$cluster)) == 2) {
    med <- tapply(v, km$cluster, median, na.rm = TRUE)
    hi <- as.integer(names(med)[which.max(med)])
    lo <- as.integer(names(med)[which.min(med)])
    frac <- mean(km$cluster == hi)
    q_hi <- as.numeric(stats::quantile(v[km$cluster == hi], 0.20, names = FALSE))
    q_lo <- as.numeric(stats::quantile(v[km$cluster == lo], 0.80, names = FALSE))
    sep_ok <- is.finite(med[as.character(hi)]) && is.finite(med[as.character(lo)]) &&
      (unname(med[as.character(hi)] - med[as.character(lo)]) >= min_sep) &&
      is.finite(q_hi) && is.finite(q_lo) && q_hi >= q_lo + 0.08
    if (sep_ok && frac >= min_pos_frac && frac <= max_pos_frac) {
      return(unname((med[as.character(hi)] + med[as.character(lo)]) / 2))
    }
  }
  d <- tryCatch(stats::density(v, n = 512, adjust = 1.15), error = function(e) NULL)
  if (!is.null(d) && length(d$y) > 10) {
    i0 <- which.max(d$y)
    y <- d$y
    x <- d$x
    peak <- y[i0]
    if (i0 < length(y) - 10) {
      i_sec <- NA_integer_
      for (i in seq.int(i0 + 6L, length(y) - 3L)) {
        if (x[i] - x[i0] < min_sep * 0.6) next
        if (y[i] >= y[i - 1L] && y[i] >= y[i + 1L] && y[i] > peak * 0.06) {
          mid_min <- min(y[i0:i])
          if (is.finite(mid_min) && mid_min < min(peak, y[i]) * 0.5) {
            i_sec <- i
            break
          }
        }
      }
      if (is.finite(i_sec)) {
        i_val <- i0 + which.min(y[i0:i_sec]) - 1L
        cut <- x[i_val]
        if (accept(cut)) return(unname(cut))
      }
      i_edge <- NA_integer_
      for (i in seq.int(i0 + 1L, length(y))) {
        if (y[i] <= peak * 0.12) {
          i_edge <- i
          break
        }
      }
      if (is.finite(i_edge)) {
        cut <- x[i_edge]
        q98 <- as.numeric(stats::quantile(v, 0.98, names = FALSE))
        if (is.finite(q98) && (q98 - cut) >= min_sep * 0.7 && accept(cut)) {
          return(unname(cut))
        }
      }
    }
  }
  q05 <- as.numeric(stats::quantile(v, 0.05, names = FALSE))
  q50 <- as.numeric(stats::quantile(v, 0.50, names = FALSE))
  q95 <- as.numeric(stats::quantile(v, 0.95, names = FALSE))
  spread <- q95 - q05
  if (is.finite(spread) && spread < min_sep * 1.2) {
    if (is.finite(q50) && q50 >= 1.2) return(unname(min(q05, q50 - min_sep)))
    return(unname(max(q95, q50 + min_sep)))
  }
  # 切不到独立阳性群时把门放在云团上方，整群阴性；不要用 q80 横切主团
  unname(max(q95, q50 + min_sep))
}

# 一层：k=2 圈出该标志高的一群；两群分不开、或高群占了多数则整层跳过
gate_k2_high <- function(mat, idx, markers, min_sep = 0.15, max_frac = 1) {
  n <- length(idx)
  out <- rep(FALSE, n)
  if (n < 40) return(out)
  sc <- marker_score(mat, idx, markers)
  ok <- is.finite(sc)
  if (sum(ok) < 40) return(out)
  xs <- scale(as.matrix(sc[ok]))
  xs[!is.finite(xs)] <- 0
  set.seed(seed_value)
  km <- tryCatch(stats::kmeans(xs, centers = 2, nstart = 10, iter.max = 250, algorithm = "Lloyd"),
                 error = function(e) NULL)
  if (is.null(km)) return(out)
  med <- tapply(sc[ok], km$cluster, median, na.rm = TRUE)
  if (length(med) < 2 || diff(range(as.numeric(med))) < min_sep) return(out)
  hi <- as.integer(names(med)[which.max(med)])
  hit <- km$cluster == hi
  if (mean(hit) > max_frac) return(out)
  out[ok] <- hit
  out
}

# CD62L vs CD44 完整象限：naive (62L+44-) / T_CM (62L+44+) / T_EM (62L-44+ 及双阴)
split_memory_3 <- function(mat, idx, prefix) {
  n <- length(idx)
  if (n == 0) return(character(0))
  cd62 <- colv(mat, "CD62L")[idx]
  cd44 <- colv(mat, "CD44")[idx]
  if (n < 8 || !any(is.finite(cd62)) || !any(is.finite(cd44))) {
    return(rep(mem_from_cd62_cd44(
      stats::median(cd62, na.rm = TRUE),
      stats::median(cd44, na.rm = TRUE),
      prefix
    ), n))
  }
  cut62 <- axis_pos_cut(cd62)
  cut44 <- axis_pos_cut(cd44)
  hi62 <- is.finite(cd62) & cd62 >= cut62
  hi44 <- is.finite(cd44) & cd44 >= cut44
  out <- rep(paste0(prefix, "_TEM"), n)
  out[hi62 & !hi44] <- paste0(prefix, "_naive")
  out[hi62 & hi44] <- paste0(prefix, "_TCM")
  out[!hi62 & hi44] <- paste0(prefix, "_TEM")
  out
}

refine_memory_subsets <- function(mat, idx, labs, prefix) {
  n <- length(labs)
  if (!n) return(labs)
  i_cm <- which(labs == paste0(prefix, "_TCM"))
  if (length(i_cm) >= 20 && all(c("CD27", "CD95") %in% colnames(mat))) {
    hi27 <- gate_k2_high(mat, idx[i_cm], "CD27", 0.15)
    hi95 <- gate_k2_high(mat, idx[i_cm], "CD95", 0.15)
    both <- hi27 & hi95
    if (any(both)) labs[i_cm[both]] <- paste0(prefix, "_TSCM")
  }
  i_em <- which(labs == paste0(prefix, "_TEM"))
  if (!length(i_em)) return(labs)
  remain <- i_em
  take_em <- function(markers, label, min_sep) {
    if (!length(remain)) return(invisible())
    hi <- gate_k2_high(mat, idx[remain], markers, min_sep)
    if (!any(hi)) return(invisible())
    labs[remain[hi]] <<- label
    remain <<- remain[!hi]
  }
  take_em(c("PD-L1", "LAG-3", "TIM-3"), paste0(prefix, "_exhausted"), 0.15)
  if ("SCA-1" %in% colnames(mat)) take_em("SCA-1", paste0(prefix, "_SLEC"), 0.15)
  mpec_lab <- if ("SCA-1" %in% colnames(mat)) paste0(prefix, "_MPEC") else paste0(prefix, "_TEM_early")
  take_em("CD27", mpec_lab, 0.15)
  if (length(remain)) labs[remain] <- paste0(prefix, "_TEM_late")
  labs
}

sequential_t_subsets <- function(mat, idx, line) {
  n <- length(idx)
  if (n == 0) return(character(0))
  labs <- rep(if (identical(line, "CD4")) "CD4_TEM" else "CD8_TEM", n)
  remain <- rep(TRUE, n)
  take_high <- function(markers, label, min_sep) {
    if (!any(remain)) return(invisible())
    hi <- gate_k2_high(mat, idx[remain], markers, min_sep)
    pos <- which(remain)
    labs[pos[hi]] <<- label
    remain[pos[hi]] <<- FALSE
  }
  if (identical(line, "CD4")) {
    take_high("CD69", "CD4_activated", 0.15)
    take_high("CD25", "Treg", 0.15)
    take_high(c("IFN-g", "TNF-a"), "CD4_effector", 0.25)
    if (any(remain)) labs[remain] <- split_memory_3(mat, idx[remain], "CD4")
    labs <- refine_memory_subsets(mat, idx, labs, "CD4")
  } else {
    take_high(c("GZMB", "Perforin"), "CD8_effector", 0.15)
    take_high(c("IFN-g", "TNF-a"), "CD8_effector", 0.25)
    take_high("CD69", "CD8_activated", 0.15)
    if (any(remain)) labs[remain] <- split_memory_3(mat, idx[remain], "CD8")
    labs <- refine_memory_subsets(mat, idx, labs, "CD8")
  }
  labs
}

# 剩余 NK：Hayakawa CD27 vs CD11b（未成熟 / 双阳 / 成熟）；双阴留 NK
split_nk_cd27_cd11b <- function(mat, idx) {
  n <- length(idx)
  if (n == 0) return(character(0))
  cd27 <- colv(mat, "CD27")[idx]
  cd11b <- colv(mat, "CD11B")[idx]
  if (n < 8 || (!any(is.finite(cd27)) && !any(is.finite(cd11b)))) {
    return(rep("NK", n))
  }
  cut27 <- if (any(is.finite(cd27))) axis_pos_cut(cd27) else Inf
  cut11 <- if (any(is.finite(cd11b))) axis_pos_cut(cd11b) else Inf
  hi27 <- is.finite(cd27) & is.finite(cut27) & cd27 >= cut27
  hi11 <- is.finite(cd11b) & is.finite(cut11) & cd11b >= cut11
  out <- rep("NK", n)
  out[hi27 & !hi11] <- "NK_immature"
  out[hi27 & hi11] <- "NK_DP"
  out[!hi27 & hi11] <- "NK_mature"
  out
}

# NK 内互斥：CD69 → 杀伤 → 耗竭 → 剩余 CD27 vs CD11b。不要用 NKG2D 当亚群。
sequential_nk_subsets <- function(mat, idx) {
  n <- length(idx)
  if (n == 0) return(character(0))
  labs <- rep("NK", n)
  remain <- rep(TRUE, n)
  take_high <- function(markers, label, min_sep) {
    if (!any(remain)) return(invisible())
    hi <- gate_k2_high(mat, idx[remain], markers, min_sep)
    pos <- which(remain)
    labs[pos[hi]] <<- label
    remain[pos[hi]] <<- FALSE
  }
  take_high("CD69", "NK_activated", 0.15)
  take_high(c("GZMB", "Perforin"), "NK_effector", 0.15)
  take_high(c("PD-L1", "LAG-3", "TIM-3"), "NK_exhausted", 0.15)
  if (any(remain)) labs[remain] <- split_nk_cd27_cd11b(mat, idx[remain])
  labs
}

# NKT 内：少数活化/效应先拿走，剩余按 CD4+ / CD4-（DN）
sequential_nkt_subsets <- function(mat, idx) {
  n <- length(idx)
  if (n == 0) return(character(0))
  labs <- rep("NKT_DN", n)
  remain <- rep(TRUE, n)
  take_high <- function(markers, label, min_sep) {
    if (!any(remain)) return(invisible())
    hi <- gate_k2_high(mat, idx[remain], markers, min_sep)
    pos <- which(remain)
    labs[pos[hi]] <<- label
    remain[pos[hi]] <<- FALSE
  }
  take_high("CD69", "NKT_activated", 0.15)
  take_high(c("GZMB", "IFN-g"), "NKT_effector", 0.25)
  if (any(remain)) {
    cd4 <- colv(mat, "CD4")[idx[remain]]
    if (sum(is.finite(cd4)) >= 8) {
      cut4 <- axis_pos_cut(cd4)
      hi4 <- is.finite(cd4) & is.finite(cut4) & cd4 >= cut4
    } else {
      hi4 <- is.finite(cd4) & cd4 >= 1.5
    }
    labs[remain] <- ifelse(hi4, "NKT_CD4", "NKT_DN")
  }
  labs
}

# P2 Naive 内：IgM 高 → MZ；少数 CD80/CD86 岛 → Activated（不要用 CD40）；剩余 Naive
sequential_b_naive <- function(mat, idx) {
  n <- length(idx)
  if (n == 0) return(character(0))
  out <- rep("Naive_B", n)
  remain <- rep(TRUE, n)
  take_high <- function(markers, label, min_sep, beat = NULL, margin = 0.25, max_frac = 1) {
    if (!any(remain)) return(invisible())
    hi <- gate_k2_high(mat, idx[remain], markers, min_sep, max_frac = max_frac)
    if (!any(hi)) return(invisible())
    pos <- which(remain)
    if (length(beat)) {
      hi_idx <- idx[remain][hi]
      sc <- median(marker_score(mat, hi_idx, markers), na.rm = TRUE)
      bt <- max(vapply(beat, function(m) median(colv(mat, m)[hi_idx], na.rm = TRUE), numeric(1)))
      if (!is.finite(sc) || !is.finite(bt) || sc <= bt + margin) return(invisible())
    }
    out[pos[hi]] <<- label
    remain[pos[hi]] <<- FALSE
  }
  take_high("IgM", "MZ_B", 0.35, max_frac = 0.55)
  take_high(c("CD80", "CD86"), "Activated_B", 0.5, beat = "IgD", margin = 0.25, max_frac = 0.35)
  if (any(remain)) out[remain] <- "Naive_B"
  out
}

# P2 Unswitched（IgD+ CD27+）内：IgM 高 → MZ；少数 CD80/CD86 → Activated
sequential_b_unswitched <- function(mat, idx) {
  n <- length(idx)
  if (n == 0) return(character(0))
  out <- rep("Unswitched_B", n)
  remain <- rep(TRUE, n)
  take_high <- function(markers, label, min_sep, beat = NULL, margin = 0.25, max_frac = 1) {
    if (!any(remain)) return(invisible())
    hi <- gate_k2_high(mat, idx[remain], markers, min_sep, max_frac = max_frac)
    if (!any(hi)) return(invisible())
    pos <- which(remain)
    if (length(beat)) {
      hi_idx <- idx[remain][hi]
      sc <- median(marker_score(mat, hi_idx, markers), na.rm = TRUE)
      bt <- max(vapply(beat, function(m) median(colv(mat, m)[hi_idx], na.rm = TRUE), numeric(1)))
      if (!is.finite(sc) || !is.finite(bt) || sc <= bt + margin) return(invisible())
    }
    out[pos[hi]] <<- label
    remain[pos[hi]] <<- FALSE
  }
  take_high("IgM", "MZ_B", 0.35, max_frac = 0.55)
  take_high(c("CD80", "CD86"), "Activated_B", 0.5, beat = "IgD", margin = 0.25, max_frac = 0.35)
  if (any(remain)) out[remain] <- "Unswitched_B"
  out
}

# P2 Switched（IgD- CD27+）内：BLIMP 少数 → Plasmablast；少数 CD80/CD86 → Activated
# IgG+ 仍是 Switched memory。CD40 不当激活。
sequential_b_switched <- function(mat, idx) {
  n <- length(idx)
  if (n == 0) return(character(0))
  out <- rep("Switched_B", n)
  remain <- rep(TRUE, n)
  take_high <- function(markers, label, min_sep, beat = NULL, margin = 0.2, max_frac = 1) {
    if (!any(remain)) return(invisible())
    hi <- gate_k2_high(mat, idx[remain], markers, min_sep, max_frac = max_frac)
    if (!any(hi)) return(invisible())
    pos <- which(remain)
    if (length(beat)) {
      hi_idx <- idx[remain][hi]
      sc <- median(marker_score(mat, hi_idx, markers), na.rm = TRUE)
      bt <- max(vapply(beat, function(m) median(colv(mat, m)[hi_idx], na.rm = TRUE), numeric(1)))
      if (!is.finite(sc) || !is.finite(bt) || sc <= bt + margin) return(invisible())
    }
    out[pos[hi]] <<- label
    remain[pos[hi]] <<- FALSE
  }
  take_high("BLIMP-1", "Plasmablast", 0.4, beat = "IgD", margin = 0.3, max_frac = 0.55)
  take_high(c("CD80", "CD86"), "Activated_B", 0.5, beat = c("IgD", "IgM", "IgG"), margin = 0.25, max_frac = 0.35)
  if (any(remain)) out[remain] <- "Switched_B"
  out
}

# 旧 Memory 大类：按 switched 路径拆；IgM 高记 Unswitched
sequential_b_memory <- function(mat, idx) {
  n <- length(idx)
  if (n == 0) return(character(0))
  out <- sequential_b_switched(mat, idx)
  remain <- out == "Switched_B"
  if (any(remain)) {
    hi <- gate_k2_high(mat, idx[remain], "IgM", 0.15, max_frac = 0.6)
    pos <- which(remain)
    if (any(hi)) {
      hi_idx <- idx[remain][hi]
      igm <- median(colv(mat, "IgM")[hi_idx], na.rm = TRUE)
      igd <- median(colv(mat, "IgD")[hi_idx], na.rm = TRUE)
      igg <- median(colv(mat, "IgG")[hi_idx], na.rm = TRUE)
      if (is.finite(igm) && igm > igd + 0.1 && igm >= igg) {
        out[pos[hi]] <- "Unswitched_B"
      }
    }
  }
  out
}

gate_k2_and <- function(mat, idx, markers, min_sep = 0.15, max_frac = 1) {
  markers <- intersect(as.character(markers), colnames(mat))
  if (!length(idx) || !length(markers)) return(rep(FALSE, length(idx)))
  hits <- lapply(markers, function(m) gate_k2_high(mat, idx, m, min_sep, max_frac))
  Reduce(`&`, hits)
}

# P3 巨噬：M1 = iNOS+ ARG-1-；M2 = ARG-1+ CD206+ iNOS-。不用 IL-10/TGF-b/TNF-a/CD86
sequential_mac <- function(mat, idx) {
  n <- length(idx)
  if (n == 0) return(character(0))
  out <- rep("Macrophage", n)
  inos_hi <- gate_k2_high(mat, idx, "iNOS", 0.2)
  arg_hi <- gate_k2_high(mat, idx, "ARG-1", 0.2)
  cd206_hi <- gate_k2_high(mat, idx, "CD206", 0.2)
  m1 <- inos_hi & !arg_hi
  if (any(m1)) out[m1] <- "M1_like_Mac"
  m2 <- !inos_hi & arg_hi & (cd206_hi | !("CD206" %in% colnames(mat)))
  if (!any(arg_hi) && any(cd206_hi)) m2 <- !inos_hi & cd206_hi
  m2 <- m2 & !m1
  if (any(m2)) out[m2] <- "M2_like_Mac"
  out
}

# CD11B+：Ly6G 中性粒 → Siglec-F 且 CCR3 嗜酸 → FceRI 且 CD200R3 肥大
# → F4/80 high 组织巨噬 → 剩余 F4/80 mid/low 里 Ly6C hi 炎症单核
sequential_cd11b_pos <- function(mat, idx) {
  n <- length(idx)
  if (n == 0) return(character(0))
  out <- rep("Mono_Ly6Clo", n)
  remain <- rep(TRUE, n)
  take_high <- function(markers, label, min_sep, beat = NULL, margin = 0.15, min_abs = 1.2, and = FALSE) {
    if (!any(remain)) return(invisible())
    hi <- if (isTRUE(and)) {
      gate_k2_and(mat, idx[remain], markers, min_sep)
    } else {
      gate_k2_high(mat, idx[remain], markers, min_sep)
    }
    if (!any(hi)) return(invisible())
    pos <- which(remain)
    hi_idx <- idx[remain][hi]
    sc <- median(marker_score(mat, hi_idx, markers), na.rm = TRUE)
    if (!is.finite(sc) || sc < min_abs) return(invisible())
    if (length(beat)) {
      bt <- max(vapply(beat, function(m) median(colv(mat, m)[hi_idx], na.rm = TRUE), numeric(1)))
      if (!is.finite(bt) || sc <= bt + margin) return(invisible())
    }
    out[pos[hi]] <<- label
    remain[pos[hi]] <<- FALSE
  }
  take_high("LY6G", "Neutrophil", 0.2, beat = c("Siglec-F", "FceRI"), margin = 0.15)
  take_high(c("Siglec-F", "CCR3"), "Eosinophil", 0.2, beat = "LY6G", margin = 0.2, and = TRUE)
  take_high(c("FceRI", "CD200R3"), "Mast", 0.2, beat = c("LY6G", "Siglec-F"), margin = 0.15, and = TRUE)
  take_high("F4/80", "Macrophage", 0.25, beat = c("LY6G", "Siglec-F"), margin = 0.1, min_abs = 1.8)
  if (any(remain)) {
    f480 <- colv(mat, "F4/80")[idx[remain]]
    # 组织巨噬 = F4/80 high；单峰高时 k=2/axis_pos 会把门切到云团上方
    hi_f480 <- is.finite(f480) & f480 >= 1.8
    if (any(hi_f480)) {
      pos <- which(remain)
      out[pos[hi_f480]] <- "Macrophage"
      remain[pos[hi_f480]] <- FALSE
    }
  }
  if (any(remain)) {
    hi <- gate_k2_high(mat, idx[remain], "LY6C", 0.15)
    pos <- which(remain)
    if (any(hi)) {
      out[pos[hi]] <- "Mono_Ly6Chi"
      remain[pos[hi]] <- FALSE
    }
    out[remain] <- "Mono_Ly6Clo"
  }
  i_mac <- which(out == "Macrophage")
  if (length(i_mac)) out[i_mac] <- sequential_mac(mat, idx[i_mac])
  out
}

# CD11B-：CD11C+ 且 MHC-II+ → DC，再 CD103+ cDC1 / CD103- cDC2
sequential_cd11b_neg <- function(mat, idx) {
  n <- length(idx)
  if (n == 0) return(character(0))
  out <- rep("other", n)
  dc <- gate_k2_and(mat, idx, c("CD11C", "I-A/I-E"), 0.15)
  if (!any(dc) && "CD11C" %in% colnames(mat)) {
    dc <- gate_k2_high(mat, idx, "CD11C", 0.2)
    if (any(dc) && "I-A/I-E" %in% colnames(mat)) {
      mhc <- median(colv(mat, "I-A/I-E")[idx[dc]], na.rm = TRUE)
      if (!is.finite(mhc) || mhc < 1.2) dc[] <- FALSE
    }
  }
  if (any(dc)) {
    out[dc] <- "cDC2"
    i_dc <- which(dc)
    hi103 <- gate_k2_high(mat, idx[i_dc], "CD103", 0.15)
    if (any(hi103)) out[i_dc[hi103]] <- "cDC1_CD103"
  }
  out
}

# P3 髓系：先按 CD11B 分两支，不要一次 kmeans 把 Ly6C/MHCII 和 Ly6G/F4/80 比大小
sequential_myeloid <- function(mat, idx) {
  n <- length(idx)
  if (n == 0) return(character(0))
  cd11b <- colv(mat, "CD11B")[idx]
  if (!"CD11B" %in% colnames(mat) || !any(is.finite(cd11b))) {
    return(sequential_cd11b_pos(mat, idx))
  }
  pos <- na_to_false(p2_pos_mask(cd11b))
  out <- rep("other", n)
  if (any(pos)) out[pos] <- sequential_cd11b_pos(mat, idx[pos])
  if (any(!pos)) out[!pos] <- sequential_cd11b_neg(mat, idx[!pos])
  out
}

# 分层圈门：P1/P2/P3 都是先圈大类，再在类内用其余抗体/细胞因子分亚群
hierarchical_gate <- function(mat, panel_id) {
  mat <- as.matrix(mat)
  major <- gate_major_lineage(mat, panel_id)
  subset <- major
  if (panel_id == "P1") {
    i4 <- which(major == "CD4")
    subset[i4] <- sequential_t_subsets(mat, i4, "CD4")
    i8 <- which(major == "CD8")
    subset[i8] <- sequential_t_subsets(mat, i8, "CD8")
    ink <- which(major == "NK")
    subset[ink] <- sequential_nk_subsets(mat, ink)
    inkt <- which(major == "NKT")
    subset[inkt] <- sequential_nkt_subsets(mat, inkt)
  } else if (panel_id == "P2") {
    i_n <- which(major == "Naive_B")
    subset[i_n] <- sequential_b_naive(mat, i_n)
    i_u <- which(major == "Unswitched_B")
    subset[i_u] <- sequential_b_unswitched(mat, i_u)
    i_s <- which(major == "Switched_B")
    subset[i_s] <- sequential_b_switched(mat, i_s)
    i_m <- which(major == "Memory_B")
    subset[i_m] <- sequential_b_memory(mat, i_m)
    i_a <- which(major == "Atypical_B")
    if (length(i_a)) {
      tmp <- sequential_b_naive(mat, i_a)
      subset[i_a] <- ifelse(tmp == "Activated_B", "Activated_B", "Atypical_B")
    }
  } else if (panel_id == "P3") {
    im <- which(major == "Myeloid")
    subset[im] <- sequential_myeloid(mat, im)
  }
  subset[is.na(subset) | !nzchar(subset)] <- major[is.na(subset) | !nzchar(subset)]
  list(major = major, subset = subset)
}

# 每个样品单独走同一套分层门，避免六个文件共用一个阈值把某管圈不完
hierarchical_gate_by_sample <- function(mat, samples, panel_id) {
  mat <- as.matrix(mat)
  n <- nrow(mat)
  samples <- as.character(samples)
  if (length(samples) != n) stop("samples length must match rows")
  major <- rep(NA_character_, n)
  subset <- rep(NA_character_, n)
  smp <- unique(samples)
  for (s in smp) {
    ii <- which(samples == s)
    if (length(ii) < 30) next
    h <- hierarchical_gate(mat[ii, , drop = FALSE], panel_id)
    major[ii] <- h$major
    subset[ii] <- h$subset
  }
  leftover <- which(is.na(major) | !nzchar(major))
  if (length(leftover)) {
    h <- hierarchical_gate(mat, panel_id)
    major[leftover] <- h$major[leftover]
    subset[leftover] <- h$subset[leftover]
  }
  list(major = major, subset = subset)
}

# -----------------------------------------------------------------------------
# 7. 出图（每个 panel 都导出 PDF + PNG）
# -----------------------------------------------------------------------------
save_gg <- function(plot, path_stub, width = 7, height = 6) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".pdf"), plot, width = width, height = height),
           error = function(e) log_msg("pdf ggsave failed: ", e$message))
  tryCatch(ggplot2::ggsave(paste0(path_stub, ".png"), plot, width = width, height = height, dpi = 300),
           error = function(e) log_msg("png ggsave failed: ", e$message))
}

theme_dr <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold", size = 13)
    )
}

pal_group <- setNames(c("#1A1A1A", "#E31A1C"), flow_group_levels)
pal_group_shape <- setNames(c(16, 15), flow_group_levels)

# 定性高对比色：同一大类内部不要共用一条色相渐变（CD8 不全绿、CD4 不全红棕、NK 不全紫、NKT 不全黄橙）
pal_p1_hues <- c(
  "#E74C3C", "#1565C0", "#2E7D32", "#F9A825", "#6A1B9A",
  "#00897B", "#EF6C00", "#AD1457", "#0277BD", "#9E9D24",
  "#5D4037", "#00ACC1", "#C62828", "#4527A0", "#558B2F",
  "#F57C00", "#37474F", "#EC407A", "#1A237E", "#FFD600"
)

pal_celltype <- c(
  "CD4 activated" = "#E74C3C",
  "CD4 T_EFF" = "#1565C0",
  "CD4 effector" = "#1565C0",
  "CD4 T_CM" = "#2E7D32",
  "CD4 TCM" = "#2E7D32",
  "CD4 T_SCM" = "#26C6DA",
  "CD4 T_EM early" = "#1DE9B6",
  "CD4 T_EM late" = "#6A1B9A",
  "CD4 T_EM" = "#8E24AA",
  "CD4 TEM" = "#8E24AA",
  "CD4 SLEC" = "#00897B",
  "CD4 MPEC" = "#7C4DFF",
  "CD4 exhausted" = "#FF00FF",
  "Treg" = "#795548",
  "CD4 naive" = "#F5B041",
  "CD4 T" = "#EF9A9A",
  "CD8 activated" = "#00E5FF",
  "CD8 T_CM" = "#2979FF",
  "CD8 TCM" = "#2979FF",
  "CD8 T_SCM" = "#18FFFF",
  "CD8 T_EM early" = "#795548",
  "CD8 T_EM late" = "#FF6D00",
  "CD8 T_EM" = "#EF6C00",
  "CD8 TEM" = "#EF6C00",
  "CD8 T_EFF" = "#6200EA",
  "CD8 effector" = "#6200EA",
  "CD8 naive" = "#76FF03",
  "CD8 SLEC" = "#1B5E20",
  "CD8 MPEC" = "#FFD600",
  "CD8 exhausted" = "#D500F9",
  "CD8 T" = "#2E7D32",
  "NK T_EFF" = "#FF6F00",
  "NK effector" = "#FF6F00",
  "NK activated" = "#C6FF00",
  "NK immature" = "#80DEEA",
  "NK DP" = "#F48FB1",
  "NK mature" = "#00695C",
  "NK exhausted" = "#7B1FA2",
  "NK" = "#1A237E",
  "NKT" = "#546E7A",
  "CD4 NKT" = "#00BFA5",
  "DN NKT" = "#FF80AB",
  "NKT activated" = "#304FFE",
  "NKT T_EFF" = "#BF360C",
  "NKT effector" = "#BF360C",
  "B cell" = "#0D47A1",
  "Macrophage" = "#0277BD",
  "T" = "#D4A017",
  "T naive" = "#F7DC6F",
  "T T_CM" = "#43A047",
  "T TCM" = "#43A047",
  "T T_EM" = "#6A1B9A",
  "T TEM" = "#6A1B9A",
  "Myeloid" = "#00897B",
  "M1-like Mac" = "#C62828",
  "M2-like Mac" = "#43A047",
  "Naive B" = "#F5B041",
  "Atypical B" = "#7E57C2",
  "IgM memory B" = "#26A69A",
  "Unswitched memory B" = "#26A69A",
  "Memory B" = "#43A047",
  "Switched B" = "#5E35B1",
  "Switched memory B" = "#5E35B1",
  "MZ B" = "#FF8A65",
  "Plasmablast" = "#C62828",
  "Activated B" = "#EC407A",
  "Plasma" = "#4A148C",
  "Eosinophil" = "#E74C3C",
  "Neutrophil" = "#FF9800",
  "Ly6C hi mono" = "#6D4C41",
  "Ly6C lo mono" = "#BCAAA4",
  "DC" = "#7E57C2",
  "CD103 DC" = "#4527A0",
  "cDC1" = "#4527A0",
  "cDC2" = "#CE93D8",
  "Basophil" = "#37474F",
  "Basophil/mast" = "#37474F",
  "Mast" = "#37474F",
  "His+ target" = "#00ACC1",
  "Target" = "#00ACC1",
  "other" = "#B0B0B0",
  "Other" = "#B0B0B0"
)

# 图1大类色：色相拉开，不要一堆近红/近绿挤在同一张总图上
pal_major <- c(
  "CD4 T" = "#D32F2F",
  "CD8 T" = "#2E7D32",
  "T" = "#B71C1C",
  "NK" = "#6A1B9A",
  "NKT" = "#F4D03F",
  "B cell" = "#1565C0",
  "Myeloid" = "#00897B",
  "His+ target" = "#00ACC1",
  "Target" = "#00ACC1",
  "other" = "#95A5A6",
  "Other" = "#95A5A6"
)

celltype_label <- function(lineage, panel_id) {
  rec <- c(
    B = "B cell",
    Naive_B = "Naive B",
    Atypical_B = "Atypical B",
    IgM_memory = "Unswitched memory B",
    Unswitched_B = "Unswitched memory B",
    Memory_B = "Memory B",
    Switched_B = "Switched memory B",
    MZ_B = "MZ B",
    Plasmablast = "Plasmablast",
    Activated_B = "Activated B",
    Plasma = "Plasma",
    CD4_naive = "CD4 naive",
    CD4_TCM = "CD4 T_CM",
    CD4_TSCM = "CD4 T_SCM",
    CD4_TEM = "CD4 T_EM",
    CD4_TEM_early = "CD4 T_EM early",
    CD4_TEM_late = "CD4 T_EM late",
    CD4_SLEC = "CD4 SLEC",
    CD4_MPEC = "CD4 MPEC",
    CD4_exhausted = "CD4 exhausted",
    CD4_activated = "CD4 activated",
    CD4_effector = "CD4 T_EFF",
    Treg = "Treg",
    CD4_T = "CD4 T",
    CD4 = "CD4 T",
    CD8 = "CD8 T",
    CD8_naive = "CD8 naive",
    CD8_TCM = "CD8 T_CM",
    CD8_TSCM = "CD8 T_SCM",
    CD8_TEM = "CD8 T_EM",
    CD8_TEM_early = "CD8 T_EM early",
    CD8_TEM_late = "CD8 T_EM late",
    CD8_SLEC = "CD8 SLEC",
    CD8_MPEC = "CD8 MPEC",
    CD8_activated = "CD8 activated",
    CD8_effector = "CD8 T_EFF",
    CD8_exhausted = "CD8 exhausted",
    CD8_T = "CD8 T",
    NK = "NK",
    NK_activated = "NK activated",
    NK_immature = "NK immature",
    NK_DP = "NK DP",
    NK_mature = "NK mature",
    NK_exhausted = "NK exhausted",
    NK_effector = "NK T_EFF",
    NKT = "NKT",
    NKT_CD4 = "CD4 NKT",
    NKT_DN = "DN NKT",
    NKT_activated = "NKT activated",
    NKT_effector = "NKT T_EFF",
    T = "T",
    T_naive = "T naive",
    T_TCM = "T T_CM",
    T_TEM = "T T_EM",
    T_effector = "T T_EFF",
    T_T = "T",
    Myeloid = "Myeloid",
    Macrophage = "Macrophage",
    M1_like_Mac = "M1-like Mac",
    M2_like_Mac = "M2-like Mac",
    Neutrophil = "Neutrophil",
    Eosinophil = "Eosinophil",
    Mono_Ly6Chi = "Ly6C hi mono",
    Mono_Ly6Clo = "Ly6C lo mono",
    DC = "DC",
    cDC1_CD103 = "cDC1",
    cDC2 = "cDC2",
    DC = "DC",
    Mast = "Mast",
    Basophil_mast = "Mast",
    Basophil = "Mast",
    Other = "other",
    Target = "His+ target",
    His_target = "His+ target"
  )
  lab <- as.character(lineage)
  out <- rec[lab]
  unname(ifelse(is.na(out), lab, out))
}

celltype_colors <- function(levels) {
  pal <- pal_celltype[levels]
  miss <- levels[is.na(pal)]
  if (length(miss) > 0) {
    used <- unique(unname(pal[!is.na(pal)]))
    extras <- setdiff(pal_p1_hues, used)
    if (length(extras) < length(miss)) extras <- c(extras, pal_p1_hues)
    pal[is.na(pal)] <- extras[seq_along(miss)]
  }
  names(pal) <- levels
  pal
}

# 同一张图里同类不要挤在一条色相上（RGB 欧氏距离；测试用）
palette_min_rgb_dist <- function(cols) {
  cols <- unname(cols)
  cols <- cols[!is.na(cols) & nzchar(cols)]
  if (length(cols) < 2L) return(Inf)
  rgb <- t(grDevices::col2rgb(cols))
  dmin <- Inf
  for (i in seq_len(nrow(rgb) - 1L)) {
    for (j in (i + 1L):nrow(rgb)) {
      d <- sqrt(sum((rgb[i, ] - rgb[j, ])^2))
      if (d < dmin) dmin <- d
    }
  }
  dmin
}

# 细胞少时点要大，否则大画布上看不见群
dr_point_size <- function(n) {
  n <- max(1, as.numeric(n)[1])
  sz <- 4.0 - 0.85 * log10(n)
  max(0.55, min(2.7, sz))
}

dr_point_alpha <- function(n) {
  n <- max(1, as.numeric(n)[1])
  if (n < 500) 1 else if (n < 2500) 0.94 else 0.8
}

embedding_axis_limits <- function(x, y, pad = 0.07) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(list(x = c(-1, 1), y = c(-1, 1)))
  xr <- range(x[ok])
  yr <- range(y[ok])
  if (diff(xr) < 1e-6) xr <- xr + c(-0.5, 0.5)
  if (diff(yr) < 1e-6) yr <- yr + c(-0.5, 0.5)
  # 正方形视野，两边一起裁空白，不要只剩稀稀拉拉的小点飘在大白纸上
  span <- max(diff(xr), diff(yr))
  mx <- mean(xr)
  my <- mean(yr)
  half <- span * (0.5 + pad)
  list(x = c(mx - half, mx + half), y = c(my - half, my + half))
}

split_dr_save_size <- function(n_keys, facet = TRUE, n_cells = Inf) {
  n_keys <- max(1L, as.integer(n_keys))
  n_cells <- max(1, as.numeric(n_cells)[1])
  ncol_leg <- if (n_keys > 18L) 3L else if (n_keys > 10L) 2L else 1L
  legend_in <- 2.65 * ncol_leg
  panel_w <- if (isTRUE(facet)) {
    if (n_cells < 400) 7.0 else if (n_cells < 1500) 8.2 else 10.0
  } else {
    if (n_cells < 400) 4.8 else if (n_cells < 1500) 5.6 else 6.6
  }
  rows <- ceiling(n_keys / ncol_leg)
  height <- if (n_cells < 400) {
    max(5.15, 4.15 + 0.38 * rows)
  } else {
    max(6.2, 4.7 + 0.42 * rows)
  }
  list(width = panel_w + legend_in, height = height, ncol_leg = ncol_leg)
}

# 图注放在右侧独立留白，不要被 coord/facet 裁掉
save_split_dr <- function(plot, path_stub, n_keys, facet = TRUE, n_cells = NULL) {
  if (is.null(n_cells) && !is.null(plot$data) && nrow(plot$data) > 0) {
    n_cells <- nrow(plot$data)
  }
  if (is.null(n_cells) || !is.finite(n_cells)) n_cells <- Inf
  sz <- split_dr_save_size(n_keys, facet, n_cells)
  plot <- plot +
    ggplot2::guides(color = ggplot2::guide_legend(
      ncol = sz$ncol_leg,
      override.aes = list(size = 4.8, alpha = 1),
      title = NULL
    )) +
    ggplot2::theme(
      legend.position = "right",
      legend.justification = c(0, 0.5),
      legend.text = ggplot2::element_text(size = 11),
      legend.key.size = grid::unit(0.5, "cm"),
      legend.spacing.y = grid::unit(0.16, "cm"),
      legend.margin = ggplot2::margin(4, 10, 4, 8),
      plot.margin = ggplot2::margin(10, 28, 14, 12)
    )
  save_gg(plot, path_stub, width = sz$width, height = sz$height)
}

theme_split_dr <- function() {
  ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 16),
      legend.title = ggplot2::element_blank(),
      legend.position = "right",
      axis.line = ggplot2::element_line(
        linewidth = 0.6,
        arrow = grid::arrow(type = "open", length = grid::unit(0.12, "inches"))
      ),
      panel.border = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )
}

plot_feature_cols <- function(df) {
  meta <- c("sample", "group", "cluster", "cluster_lineage", "lineage", "true_lineage",
            "bio_sample", "tech_rep", "dimred_major",
            "UMAP1", "UMAP2", "tSNE1", "tSNE2", "sep1", "sep2", "major_plot")
  num <- names(df)[vapply(df, is.numeric, logical(1))]
  setdiff(num, meta)
}

# color_mode="major"：全体细胞按大类（*_major_split）。
# color_mode="subset"：全体细胞按细亚群（*_lineage_split），或各大类自己的亚群图。
plot_split_lineage <- function(df, x, y, panel_id, xlab, ylab, title,
                               color_mode = c("major", "subset"), ...) {
  color_mode <- match.arg(color_mode)
  if (!x %in% names(df) || !y %in% names(df)) {
    stop("plot_split_lineage needs shared embedding columns ", x, " and ", y)
  }
  plot_df <- df
  plot_df$group <- factor(plot_df$group, levels = flow_group_levels)
  keep <- is.finite(plot_df[[x]]) & is.finite(plot_df[[y]]) & !is.na(plot_df$group)
  plot_df <- plot_df[keep, , drop = FALSE]
  if (identical(color_mode, "major")) {
    maj <- dimred_major_of(
      panel_id,
      plot_df$lineage,
      if ("cluster_lineage" %in% names(plot_df)) plot_df$cluster_lineage else NULL
    )
    plot_df$celltype <- major_display_label(maj)
    prefer <- names(pal_major)
    pal_fun <- major_colors
  } else {
    plot_df$celltype <- celltype_label(plot_df$lineage, panel_id)
    prefer <- names(pal_celltype)
    pal_fun <- celltype_colors
  }
  levs <- unique(as.character(plot_df$celltype))
  levs <- c(intersect(prefer, levs), setdiff(levs, prefer))
  plot_df$celltype <- factor(plot_df$celltype, levels = levs)
  pal <- pal_fun(levs)
  # 多数细胞先画、稀有亚群后画（叠在上面），避免耗竭/naive 被主团盖住
  freq <- table(plot_df$celltype)
  plot_df <- plot_df[order(as.integer(freq[as.character(plot_df$celltype)]),
                           decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  n <- nrow(plot_df)
  lims <- embedding_axis_limits(plot_df[[x]], plot_df[[y]])
  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = celltype)) +
    ggplot2::geom_point(size = dr_point_size(n), alpha = dr_point_alpha(n), stroke = 0.12) +
    ggplot2::facet_wrap(~group, ncol = 2, scales = "fixed") +
    ggplot2::scale_color_manual(values = pal, drop = FALSE) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(size = 3.8, alpha = 1), ncol = 1
    )) +
    ggplot2::labs(title = title, x = xlab, y = ylab, color = NULL) +
    ggplot2::coord_fixed(ratio = 1, xlim = lims$x, ylim = lims$y, expand = FALSE, clip = "off") +
    theme_split_dr() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(0.7, "lines"),
      legend.key = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 18, 10, 10)
    )
}

embed_class_cells <- function(sub) {
  if (nrow(sub) < 40L) return(sub)
  feat <- plot_feature_cols(sub)
  if (length(feat) < 2L) return(sub)
  mat <- as.matrix(sub[, feat, drop = FALSE])
  mat[!is.finite(mat)] <- 0
  sc <- tryCatch(scale(mat), error = function(e) mat)
  sc[!is.finite(sc)] <- 0
  maj <- if ("dimred_major" %in% names(sub)) {
    names(sort(table(as.character(sub$dimred_major)), decreasing = TRUE))[1]
  } else {
    NA_character_
  }
  sc <- weight_class_dr_markers(sc, maj)
  pca <- tryCatch(run_pca(sc, npcs = min(12, ncol(sc))), error = function(e) NULL)
  if (is.null(pca) || is.null(pca$embedding) || ncol(pca$embedding) < 2) return(sub)
  npcs_use <- min(ncol(pca$embedding), max(2, ncol(sc)))
  pca_use <- pca$embedding[, seq_len(npcs_use), drop = FALSE]
  um <- tryCatch(run_umap(pca_use, compact = TRUE), error = function(e) NULL)
  ts <- tryCatch(run_tsne(pca_use, compact = TRUE), error = function(e) NULL)
  if (!is.null(um) && ncol(um) >= 2 && nrow(um) == nrow(sub)) {
    sub$UMAP1 <- um[, 1]
    sub$UMAP2 <- um[, 2]
  }
  if (!is.null(ts) && ncol(ts) >= 2 && nrow(ts) == nrow(sub)) {
    sub$tSNE1 <- ts[, 1]
    sub$tSNE2 <- ts[, 2]
  }
  sub
}

export_major_subset_dimred <- function(cells, panel_id, out_dir) {
  cells$dimred_major <- dimred_major_of(
    panel_id,
    cells$lineage,
    if ("cluster_lineage" %in% names(cells)) cells$cluster_lineage else NULL
  )
  majors <- unique(as.character(cells$dimred_major))
  majors <- majors[!is.na(majors) & nzchar(majors) & !majors %in% c("other", "dump", "", "Target")]
  class_dir <- file.path(out_dir, "dimred_by_major")
  dir.create(class_dir, recursive = TRUE, showWarnings = FALSE)
  for (mj in majors) {
    hit <- !is.na(cells$dimred_major) & as.character(cells$dimred_major) == mj
    sub <- cells[hit, , drop = FALSE]
    if (nrow(sub) < 40L) {
      log_msg(panel_id, " ", mj, ": skip class dimred (n=", nrow(sub), ")")
      next
    }
    tryCatch({
      sub <- embed_class_cells(sub)
      lab <- major_display_label(mj)
      tag <- paste0(panel_id, "_", gsub("[^A-Za-z0-9_-]", "_", mj))
      ttl <- paste0(panel_id, "  ", lab, "  subsets  JZ-EVB | JZ-AB")
      n_keys <- length(unique(celltype_label(sub$lineage, panel_id)))
      save_split_dr(
        plot_split_lineage(sub, "tSNE1", "tSNE2", panel_id, "tSNE-1", "tSNE-2",
                           ttl, color_mode = "subset"),
        file.path(class_dir, paste0(tag, "_tSNE_subset_JZ_AB_vs_JZ_EVB")),
        n_keys
      )
      save_split_dr(
        plot_split_lineage(sub, "UMAP1", "UMAP2", panel_id, "UMAP-1", "UMAP-2",
                           ttl, color_mode = "subset"),
        file.path(class_dir, paste0(tag, "_UMAP_subset_JZ_AB_vs_JZ_EVB")),
        n_keys
      )
      log_msg(panel_id, " ", lab, " subset dimred n=", nrow(sub), " -> ", class_dir)
    }, error = function(e) {
      log_msg(panel_id, " skip ", mj, " class dimred: ", e$message)
    })
  }
  invisible(class_dir)
}

plot_embedding <- function(df, x, y, color_col, title, point_size = NULL) {
  plot_df <- df
  if (identical(color_col, "lineage")) {
    plot_df$lineage <- celltype_label(plot_df$lineage, NA)
  }
  nlev <- length(unique(plot_df[[color_col]]))
  if (is.null(point_size)) point_size <- dr_point_size(nrow(plot_df))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = .data[[color_col]])) +
    ggplot2::geom_point(size = point_size, alpha = dr_point_alpha(nrow(plot_df)), stroke = 0.1) +
    ggplot2::coord_equal() +
    theme_dr() +
    ggplot2::labs(title = title, color = color_col, x = x, y = y)
  if (color_col == "group") {
    p <- p + ggplot2::scale_color_manual(values = pal_group)
  } else if (identical(color_col, "lineage")) {
    p <- p + ggplot2::scale_color_manual(values = celltype_colors(unique(as.character(plot_df$lineage))))
  } else if (is.numeric(plot_df[[color_col]])) {
    p <- p + ggplot2::scale_color_gradientn(colours = c("#440154", "#21908C", "#FDE725"))
  } else {
    pal <- pal_p1_hues[((seq_len(nlev) - 1L) %% length(pal_p1_hues)) + 1L]
    p <- p + ggplot2::scale_color_manual(values = pal)
  }
  p
}

plot_marker_embedding <- function(df, x, y, marker, title) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = .data[[marker]])) +
    ggplot2::geom_point(size = dr_point_size(nrow(df)), alpha = dr_point_alpha(nrow(df)), stroke = 0.08) +
    ggplot2::scale_color_gradientn(colours = c("#0D0887", "#CC4678", "#F0F921")) +
    ggplot2::coord_equal() +
    theme_dr() +
    ggplot2::labs(title = title, color = marker, x = x, y = y)
}

plot_density_split <- function(df, x, y, title) {
  df$group <- factor(df$group, levels = flow_group_levels)
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
    ggplot2::stat_density_2d_filled(contour_var = "ndensity", bins = 12) +
    ggplot2::facet_wrap(~group, ncol = 2) +
    ggplot2::coord_equal() +
    ggplot2::guides(fill = "none") +
    theme_dr() +
    ggplot2::labs(title = title, x = x, y = y)
}

plot_cluster_heatmap <- function(med, title, outfile) {
  z <- scale(med)
  z[is.na(z)] <- 0
  df <- as.data.frame(z)
  df$cluster <- rownames(df)
  long <- tidyr_pivot(df)
  p <- ggplot2::ggplot(long, ggplot2::aes(x = marker, y = cluster, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggplot2::labs(title = title, x = NULL, y = NULL, fill = "z-score")
  save_gg(p, outfile, width = max(8, 0.35 * ncol(med) + 3), height = max(5, 0.35 * nrow(med) + 2))
  if (has_pkg("pheatmap")) {
    tryCatch({
      grDevices::pdf(paste0(outfile, "_pheatmap.pdf"), width = max(8, 0.35 * ncol(med) + 3),
                     height = max(5, 0.35 * nrow(med) + 2))
      pheatmap::pheatmap(t(z), main = title, border_color = NA)
      grDevices::dev.off()
      grDevices::png(paste0(outfile, "_pheatmap.png"), width = 1800, height = 1400, res = 180)
      pheatmap::pheatmap(t(z), main = title, border_color = NA)
      grDevices::dev.off()
    }, error = function(e) log_msg("pheatmap failed: ", e$message))
  }
}

# 按已圈好的免疫细胞类型（不是 kmeans cluster 号）算标志物中位数，做注释热图
heatmap_marker_cols <- function(cells, panel_id) {
  feat <- intersect(plot_feature_cols(cells), names(cells))
  feat <- setdiff(feat, c("L/D", "Time", "time", "FSC-A", "SSC-A", "FSC-H", "SSC-H", "FSC-W", "SSC-W"))
  pref <- tryCatch(dr_marker_names(panel_id), error = function(e) character(0))
  c(intersect(pref, feat), setdiff(feat, pref))
}

lineage_median_matrix <- function(cells, panel_id, by = c("subset", "major"), min_n = 8L) {
  by <- match.arg(by)
  feat <- heatmap_marker_cols(cells, panel_id)
  if (!length(feat) || is.null(cells) || !nrow(cells)) return(NULL)
  if (identical(by, "major")) {
    raw <- dimred_major_of(
      panel_id,
      cells$lineage,
      if ("cluster_lineage" %in% names(cells)) cells$cluster_lineage else NULL
    )
    key <- major_display_label(raw)
    class_lab <- key
  } else {
    key <- celltype_label(cells$lineage, panel_id)
    class_lab <- major_display_label(dimred_major_of(
      panel_id,
      cells$lineage,
      if ("cluster_lineage" %in% names(cells)) cells$cluster_lineage else NULL
    ))
  }
  key <- as.character(key)
  class_lab <- as.character(class_lab)
  key[is.na(key) | !nzchar(key)] <- "other"
  class_lab[is.na(class_lab) | !nzchar(class_lab)] <- "other"
  keep <- key != "other" & is.finite(rowSums(as.data.frame(cells[, feat, drop = FALSE]), na.rm = TRUE))
  if (sum(keep) < 20L) return(NULL)
  mat <- as.matrix(cells[keep, feat, drop = FALSE])
  storage.mode(mat) <- "double"
  key <- key[keep]
  class_lab <- class_lab[keep]
  ntab <- table(key)
  ok <- names(ntab)[as.integer(ntab) >= min_n]
  if (length(ok) < 2L) return(NULL)
  hit <- key %in% ok
  mat <- mat[hit, , drop = FALSE]
  key <- key[hit]
  class_lab <- class_lab[hit]
  types <- unique(key)
  med <- do.call(rbind, lapply(types, function(tp) {
    apply(mat[key == tp, , drop = FALSE], 2, median, na.rm = TRUE)
  }))
  if (is.null(dim(med))) med <- matrix(med, nrow = length(types), dimnames = list(types, feat))
  rownames(med) <- types
  colnames(med) <- feat
  cls <- vapply(types, function(tp) {
    names(sort(table(class_lab[key == tp]), decreasing = TRUE))[1]
  }, character(1))
  maj_ord <- c("CD4 T", "CD8 T", "T", "NKT", "NK", "B cell", "Myeloid", "His+ target")
  ord <- order(match(cls, maj_ord, nomatch = 99L), types)
  med <- med[ord, , drop = FALSE]
  attr(med, "cell_class") <- unname(cls[ord])
  attr(med, "cell_n") <- as.integer(ntab[rownames(med)])
  med
}

plot_annotation_heatmap <- function(med, title, outfile) {
  if (is.null(med) || !nrow(med) || !ncol(med)) return(invisible(NULL))
  z <- scale(med)
  z[is.na(z)] <- 0
  df <- as.data.frame(z)
  df$cluster <- rownames(df)
  long <- tidyr_pivot(df)
  long$cluster <- factor(long$cluster, levels = rev(rownames(med)))
  long$marker <- factor(long$marker, levels = colnames(med))
  n_type <- nrow(med)
  n_mk <- ncol(med)
  w <- max(8.5, 0.48 * n_type + 3.2)
  h <- max(5.6, 0.32 * n_mk + 2.4)
  p <- ggplot2::ggplot(long, ggplot2::aes(x = cluster, y = marker, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.15) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
      axis.text.y = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL, fill = "z-score")
  save_gg(p, outfile, width = w, height = h)
  if (has_pkg("pheatmap")) {
    ann <- NULL
    cls <- attr(med, "cell_class")
    if (!is.null(cls) && length(cls) == nrow(med)) {
      ann <- data.frame(Class = unname(cls), row.names = rownames(med), stringsAsFactors = FALSE)
    }
    pal_ann <- NULL
    if (!is.null(ann)) {
      levs <- unique(ann$Class)
      pal_ann <- list(Class = major_colors(levs))
    }
    tryCatch({
      args <- list(
        mat = t(z),
        main = title,
        border_color = NA,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        fontsize_col = 9,
        fontsize_row = 9
      )
      if (!is.null(ann)) {
        args$annotation_col <- ann
        args$annotation_colors <- pal_ann
      }
      grDevices::pdf(paste0(outfile, "_pheatmap.pdf"), width = w, height = h)
      do.call(pheatmap::pheatmap, args)
      grDevices::dev.off()
      grDevices::png(paste0(outfile, "_pheatmap.png"), width = w, height = h, units = "in", res = 180)
      do.call(pheatmap::pheatmap, args)
      grDevices::dev.off()
    }, error = function(e) log_msg("annotation pheatmap failed: ", e$message))
  }
  invisible(med)
}

export_annotation_heatmaps <- function(cells, panel_id, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  med_sub <- lineage_median_matrix(cells, panel_id, "subset")
  if (!is.null(med_sub)) {
    utils::write.csv(cbind(celltype = rownames(med_sub), cell_class = attr(med_sub, "cell_class"),
                           n = attr(med_sub, "cell_n"), as.data.frame(med_sub, check.names = FALSE)),
                     file.path(out_dir, paste0(panel_id, "_annotation_heatmap.csv")),
                     row.names = FALSE)
    plot_annotation_heatmap(
      med_sub,
      paste0(panel_id, "  immune-cell annotation (subset median, z)"),
      file.path(out_dir, paste0(panel_id, "_annotation_heatmap"))
    )
    log_msg(panel_id, " annotation heatmap: ", nrow(med_sub), " cell types x ", ncol(med_sub), " markers")
  }
  med_maj <- lineage_median_matrix(cells, panel_id, "major")
  if (!is.null(med_maj)) {
    utils::write.csv(cbind(celltype = rownames(med_maj), cell_class = attr(med_maj, "cell_class"),
                           n = attr(med_maj, "cell_n"), as.data.frame(med_maj, check.names = FALSE)),
                     file.path(out_dir, paste0(panel_id, "_annotation_heatmap_major.csv")),
                     row.names = FALSE)
    plot_annotation_heatmap(
      med_maj,
      paste0(panel_id, "  immune-cell annotation (major class median, z)"),
      file.path(out_dir, paste0(panel_id, "_annotation_heatmap_major"))
    )
  }
  invisible(list(subset = med_sub, major = med_maj))
}

tidyr_pivot <- function(df) {
  mk <- setdiff(names(df), "cluster")
  out <- do.call(rbind, lapply(mk, function(m) {
    data.frame(cluster = df$cluster, marker = m, value = df[[m]], stringsAsFactors = FALSE)
  }))
  out$cluster <- factor(out$cluster, levels = rev(unique(df$cluster)))
  out$marker <- factor(out$marker, levels = mk)
  out
}

plot_freq_box <- function(freq_df, title) {
  ylab <- if (flow_should_trim_bio()) {
    "% of cells (bio-rep; max/min dropped)"
  } else {
    "% of cells (bio-rep)"
  }
  ggplot2::ggplot(freq_df, ggplot2::aes(x = cluster, y = percent, fill = group)) +
    ggplot2::stat_summary(fun = mean, geom = "col", position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.8), size = 1.6, alpha = 0.9) +
    ggplot2::scale_fill_manual(values = pal_group) +
    theme_dr() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = title, x = NULL, y = ylab)
}

plot_freq_bar <- function(sum_df, title) {
  ggplot2::ggplot(sum_df, ggplot2::aes(x = group, y = mean_percent, fill = cluster)) +
    ggplot2::geom_col(position = "stack") +
    theme_dr() +
    ggplot2::labs(title = title, x = NULL, y = "Mean % of cells")
}

fluorochrome_display <- function(fl) {
  fl <- as.character(fl)[1]
  if (!nzchar(fl) || is.na(fl)) return(fl)
  rec <- panel_map$fluorochrome_display
  if (!is.null(rec) && !is.null(rec[[fl]]) && nzchar(as.character(rec[[fl]][1]))) {
    return(as.character(rec[[fl]][1]))
  }
  fl
}

fluorochrome_of <- function(panel_id, marker) {
  items <- panel_markers(panel_id)
  hit <- vapply(items, function(it) identical(it$marker, marker), logical(1))
  if (!any(hit)) return(NA_character_)
  as.character(items[[which(hit)[1]]]$fluorochrome)
}

axis_fl_label <- function(panel_id, marker) {
  fl <- fluorochrome_of(panel_id, marker)
  if (is.na(fl) || !nzchar(fl)) return(marker)
  paste0(marker, "-", fluorochrome_display(fl))
}

p_to_star <- function(p) {
  if (!is.finite(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

# ns 时把 P 值写出来：去极端值后 n=2，柱高看起来不同，t 检验仍常不显著
p_annot_label <- function(p) {
  star <- p_to_star(p)
  if (!is.finite(p)) return(star)
  if (identical(star, "ns")) return(sprintf("ns  P = %.2f", p))
  star
}

finite_axis_lim <- function(v, pad = 0.06) {
  v <- as.numeric(v)
  v <- v[is.finite(v)]
  if (!length(v)) return(c(0, 1))
  r <- range(v)
  if (!is.finite(diff(r)) || diff(r) < 1e-6) {
    return(c(r[1] - 0.4, r[2] + 0.4))
  }
  r + c(-1, 1) * diff(r) * pad
}

safe_kde_bw <- function(v) {
  v <- as.numeric(v)
  v <- v[is.finite(v)]
  if (length(v) < 2) return(0.25)
  h <- NA_real_
  if (requireNamespace("MASS", quietly = TRUE)) {
    h <- tryCatch(MASS::bandwidth.nrd(v), error = function(e) NA_real_)
  }
  if (!is.finite(h) || h <= 0) {
    h <- stats::sd(v)
    if (!is.finite(h) || h <= 0) h <- stats::IQR(v) / 1.34
  }
  if (!is.finite(h) || h <= 0) h <- 0.25
  h
}

# FlowJo 风格：概率等高线（圈住约 25–90% 事件）+ 最外圈以外的 outliers
flow_prob_contour <- function(x, y, n_grid = 70, bw_mult = 1.45,
                              probs = c(0.25, 0.4, 0.55, 0.7, 0.82, 0.9)) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  empty <- list(
    grid = data.frame(x = numeric(0), y = numeric(0), z = numeric(0)),
    levels = numeric(0), x = x, y = y, outlier = rep(TRUE, n)
  )
  if (n < 8 || !requireNamespace("MASS", quietly = TRUE)) return(empty)
  hx <- safe_kde_bw(x) * bw_mult
  hy <- safe_kde_bw(y) * bw_mult
  rx <- range(x)
  ry <- range(y)
  lims <- c(rx + c(-1, 1) * max(diff(rx) * 0.08, 0.05),
            ry + c(-1, 1) * max(diff(ry) * 0.08, 0.05))
  kd <- tryCatch(
    MASS::kde2d(x, y, n = n_grid, h = c(hx, hy), lims = lims),
    error = function(e) NULL
  )
  if (is.null(kd) || !is.matrix(kd$z)) return(empty)
  dx <- diff(kd$x[1:2])
  dy <- diff(kd$y[1:2])
  if (!is.finite(dx) || !is.finite(dy) || dx <= 0 || dy <= 0) return(empty)
  mass <- as.numeric(kd$z) * dx * dy
  s <- sum(mass)
  if (!is.finite(s) || s <= 0) return(empty)
  mass <- mass / s
  ord <- order(as.numeric(kd$z), decreasing = TRUE)
  cum <- cumsum(mass[ord])
  zord <- as.numeric(kd$z)[ord]
  levs <- unique(vapply(probs, function(p) {
    i <- which(cum >= p)[1]
    if (is.na(i)) zord[length(zord)] else zord[i]
  }, numeric(1)))
  levs <- sort(levs[is.finite(levs) & levs > 0])
  grid <- expand.grid(x = kd$x, y = kd$y, KEEP.OUT.ATTRS = FALSE)
  grid$z <- as.vector(kd$z)
  ix <- findInterval(x, kd$x, all.inside = TRUE)
  iy <- findInterval(y, kd$y, all.inside = TRUE)
  dens <- kd$z[cbind(ix, iy)]
  outer_lev <- if (length(levs)) min(levs) else 0
  list(grid = grid, levels = levs, x = x, y = y, outlier = !is.finite(dens) | dens < outer_lev)
}

parent_mask <- function(cells, parent) {
  n <- nrow(cells)
  if (n < 1) return(logical(0))
  cl <- if ("cluster_lineage" %in% names(cells)) as.character(cells$cluster_lineage) else rep("", n)
  lin <- as.character(cells$lineage)
  cl[is.na(cl)] <- ""
  lin[is.na(lin)] <- ""
  cd4 <- cl == "CD4" | grepl("^CD4_", lin) | lin == "Treg"
  cd8 <- cl == "CD8" | grepl("^CD8_", lin)
  nk <- cl == "NK" | lin %in% nk_family()
  nkt <- cl == "NKT" | lin %in% nkt_family()
  out <- switch(as.character(parent),
    all = rep(TRUE, n),
    CD3 = cd4 | cd8 | nkt | lin == "T" | cl == "T",
    CD4 = cd4,
    CD8 = cd8,
    NK = nk,
    NKT = nkt,
    CD4_TCM = lin %in% tcm_family("CD4"),
    CD8_TCM = lin %in% tcm_family("CD8"),
    CD4_TEM = lin %in% tem_family("CD4"),
    CD8_TEM = lin %in% tem_family("CD8"),
    Myeloid = cl == "Myeloid" | lin %in% c(
      "Neutrophil", "Eosinophil", "Mast", "Basophil_mast", "Basophil", "Macrophage",
      "M1_like_Mac", "M2_like_Mac", "DC", "cDC1_CD103", "cDC2",
      "Mono_Ly6Chi", "Mono_Ly6Clo", "Myeloid"
    ),
    CD11B = lin %in% c(
      "Neutrophil", "Eosinophil", "Mast", "Basophil_mast", "Basophil",
      "Macrophage", "M1_like_Mac", "M2_like_Mac", "Mono_Ly6Chi", "Mono_Ly6Clo"
    ),
    DC = lin %in% c("DC", "cDC1_CD103", "cDC2"),
    Macrophage = lin %in% c("Macrophage", "M1_like_Mac", "M2_like_Mac"),
    CD19 = cl %in% c("Naive_B", "Unswitched_B", "Switched_B", "Atypical_B", "Memory_B") |
      lin %in% c("Naive_B", "Unswitched_B", "Switched_B", "Atypical_B", "MZ_B",
                 "Plasmablast", "Activated_B", "IgM_memory", "Memory_B"),
    Naive_Unswitched = cl %in% c("Naive_B", "Unswitched_B") |
      lin %in% c("Naive_B", "Unswitched_B", "MZ_B", "IgM_memory"),
    Unswitched_B = cl == "Unswitched_B" | lin %in% c("Unswitched_B", "IgM_memory"),
    Switched_B = cl == "Switched_B" | lin %in% c("Switched_B", "Plasmablast"),
    Memory_B = cl %in% c("Memory_B", "Unswitched_B", "Switched_B") |
      lin %in% c("Memory_B", "IgM_memory", "Unswitched_B", "Switched_B", "Plasmablast",
                 "Plasma", "Activated_B", "MZ_B"),
    Naive_B = cl == "Naive_B" | lin %in% c("Naive_B", "MZ_B"),
    Atypical_B = cl == "Atypical_B" | lin == "Atypical_B",
    Activated_B = lin == "Activated_B" | cl == "Activated_B",
    rep(TRUE, n)
  )
  out[is.na(out)] <- FALSE
  out
}

subset_hit_mask <- function(cells, spec) {
  lin <- as.character(cells$lineage)
  lin[is.na(lin)] <- ""
  if (isTRUE(spec$use_major) && "cluster_lineage" %in% names(cells)) {
    cl <- as.character(cells$cluster_lineage)
    cl[is.na(cl)] <- ""
    return(cl == spec$lineage)
  }
  fam <- switch(as.character(spec$lineage),
    CD4_TEM = tem_family("CD4"),
    CD8_TEM = tem_family("CD8"),
    NK = nk_family(),
    NKT = nkt_family(),
    DC = c("DC", "cDC1_CD103", "cDC2"),
    Mast = c("Mast", "Basophil_mast", "Basophil"),
    spec$lineage
  )
  lin %in% fam
}

# 缺染色通道：跳过这一项，不要中断后面的亚群/大类分析
skip_if_missing_channels <- function(need, available, label) {
  need <- unique(as.character(need))
  need <- need[!is.na(need) & nzchar(need)]
  miss <- setdiff(need, available)
  if (!length(miss)) return(FALSE)
  log_msg("skip ", label, " (missing channels: ", paste(miss, collapse = ", "), ")")
  TRUE
}

subset_plot_specs <- function(panel_id) {
  mk <- function(lineage, x, y, parent, ylab, use_major = FALSE,
                 gate = "box", x_hi = NA, y_hi = NA) {
    list(lineage = lineage, x = x, y = y, parent = parent, ylab = ylab,
         use_major = use_major, gate = gate, x_hi = x_hi, y_hi = y_hi)
  }
  if (identical(panel_id, "P1")) {
    return(list(
      mk("NK", "CD3", "NKp46", "all", "NK cell (%)", TRUE, gate = "quad", x_hi = FALSE, y_hi = TRUE),
      mk("NK_activated", "CD69", "NKp46", "NK", "NK activated in NK (%)", gate = "hi_hi"),
      mk("NK_effector", "GZMB", "Perforin", "NK", "NK T_EFF in NK (%)", gate = "hi_hi"),
      mk("NK_exhausted", "LAG-3", "TIM-3", "NK", "NK exhausted in NK (%)", gate = "hi_hi"),
      mk("NK_immature", "CD27", "CD11B", "NK", "NK immature in NK (%)", gate = "quad", x_hi = TRUE, y_hi = FALSE),
      mk("NK_DP", "CD27", "CD11B", "NK", "NK DP in NK (%)", gate = "quad", x_hi = TRUE, y_hi = TRUE),
      mk("NK_mature", "CD27", "CD11B", "NK", "NK mature in NK (%)", gate = "quad", x_hi = FALSE, y_hi = TRUE),
      mk("NKT", "CD3", "NKp46", "all", "NKT (%)", TRUE, gate = "quad", x_hi = TRUE, y_hi = TRUE),
      mk("NKT_CD4", "CD8", "CD4", "NKT", "CD4 NKT in NKT (%)", gate = "quad", x_hi = FALSE, y_hi = TRUE),
      mk("NKT_DN", "CD8", "CD4", "NKT", "DN NKT in NKT (%)", gate = "quad", x_hi = FALSE, y_hi = FALSE),
      mk("NKT_activated", "CD69", "CD3", "NKT", "NKT activated in NKT (%)", gate = "hi_hi"),
      mk("NKT_effector", "GZMB", "IFN-g", "NKT", "NKT T_EFF in NKT (%)", gate = "hi_hi"),
      mk("B", "CD19", "CD3", "all", "B cell (%)", gate = "quad", x_hi = TRUE, y_hi = FALSE),
      mk("Myeloid", "CD11B", "CD3", "all", "Myeloid (%)", gate = "quad", x_hi = TRUE, y_hi = FALSE),
      mk("CD4", "CD8", "CD4", "CD3", "CD4+ T cell in CD3+ (%)", TRUE, gate = "quad", x_hi = FALSE, y_hi = TRUE),
      mk("CD8", "CD8", "CD4", "CD3", "CD8+ T cell in CD3+ (%)", TRUE, gate = "quad", x_hi = TRUE, y_hi = FALSE),
      mk("CD4_naive", "CD62L", "CD44", "CD4", "CD4 naive in CD4+ (%)", gate = "quad", x_hi = TRUE, y_hi = FALSE),
      mk("CD4_TCM", "CD62L", "CD44", "CD4", "CD4 T_CM in CD4+ (%)", gate = "quad", x_hi = TRUE, y_hi = TRUE),
      mk("CD4_TSCM", "CD27", "CD95", "CD4_TCM", "CD4 T_SCM in T_CM (%)", gate = "hi_hi"),
      mk("CD4_TEM", "CD62L", "CD44", "CD4", "CD4 T_EM in CD4+ (%)", gate = "half_x", x_hi = FALSE),
      mk("CD4_TEM_early", "CD27", "CD44", "CD4_TEM", "CD4 T_EM early in T_EM (%)", gate = "hi_x", x_hi = TRUE),
      mk("CD4_TEM_late", "CD27", "CD44", "CD4_TEM", "CD4 T_EM late in T_EM (%)", gate = "quad", x_hi = FALSE, y_hi = TRUE),
      mk("CD4_SLEC", "SCA-1", "CD44", "CD4_TEM", "CD4 SLEC in T_EM (%)", gate = "hi_hi"),
      mk("CD4_MPEC", "SCA-1", "CD27", "CD4_TEM", "CD4 MPEC in T_EM (%)", gate = "quad", x_hi = FALSE, y_hi = TRUE),
      mk("CD4_exhausted", "LAG-3", "TIM-3", "CD4_TEM", "CD4 exhausted in T_EM (%)", gate = "hi_hi"),
      mk("CD4_activated", "CD69", "CD25", "CD4", "CD4 activated in CD4+ (%)", gate = "hi_hi"),
      mk("CD4_effector", "IFN-g", "TNF-a", "CD4", "CD4 T_EFF in CD4+ (%)", gate = "hi_hi"),
      mk("Treg", "CD25", "CD4", "CD4", "Treg in CD4+ (%)", gate = "hi_x", x_hi = TRUE),
      mk("CD8_naive", "CD62L", "CD44", "CD8", "CD8 naive in CD8+ (%)", gate = "quad", x_hi = TRUE, y_hi = FALSE),
      mk("CD8_TCM", "CD62L", "CD44", "CD8", "CD8 T_CM in CD8+ (%)", gate = "quad", x_hi = TRUE, y_hi = TRUE),
      mk("CD8_TSCM", "CD27", "CD95", "CD8_TCM", "CD8 T_SCM in T_CM (%)", gate = "hi_hi"),
      mk("CD8_TEM", "CD62L", "CD44", "CD8", "CD8 T_EM in CD8+ (%)", gate = "half_x", x_hi = FALSE),
      mk("CD8_TEM_early", "CD27", "CD44", "CD8_TEM", "CD8 T_EM early in T_EM (%)", gate = "hi_x", x_hi = TRUE),
      mk("CD8_TEM_late", "CD27", "CD44", "CD8_TEM", "CD8 T_EM late in T_EM (%)", gate = "quad", x_hi = FALSE, y_hi = TRUE),
      mk("CD8_SLEC", "SCA-1", "CD44", "CD8_TEM", "CD8 SLEC in T_EM (%)", gate = "hi_hi"),
      mk("CD8_MPEC", "SCA-1", "CD27", "CD8_TEM", "CD8 MPEC in T_EM (%)", gate = "quad", x_hi = FALSE, y_hi = TRUE),
      mk("CD8_activated", "CD69", "CD25", "CD8", "CD8 activated in CD8+ (%)", gate = "hi_hi"),
      mk("CD8_effector", "GZMB", "Perforin", "CD8", "CD8 T_EFF in CD8+ (%)", gate = "hi_hi"),
      mk("CD8_exhausted", "LAG-3", "TIM-3", "CD8_TEM", "CD8 exhausted in T_EM (%)", gate = "hi_hi")
    ))
  }
  if (identical(panel_id, "P2")) {
    return(list(
      mk("Naive_B", "IgD", "CD27", "CD19", "Naive B in CD19+ (%)", gate = "quad", x_hi = TRUE, y_hi = FALSE),
      mk("Unswitched_B", "IgD", "CD27", "CD19", "Unswitched memory B in CD19+ (%)", gate = "quad", x_hi = TRUE, y_hi = TRUE),
      mk("Switched_B", "IgD", "CD27", "CD19", "Switched memory B in CD19+ (%)", gate = "quad", x_hi = FALSE, y_hi = TRUE),
      mk("Atypical_B", "IgD", "CD27", "CD19", "Atypical B in CD19+ (%)", gate = "quad", x_hi = FALSE, y_hi = FALSE),
      mk("MZ_B", "IgM", "IgD", "Naive_Unswitched", "MZ B in naive/unswitched (%)", gate = "hi_x", x_hi = TRUE),
      mk("Plasmablast", "BLIMP-1", "CD27", "Switched_B", "Plasmablast in switched (%)", gate = "hi_hi"),
      mk("Plasma", "BLIMP-1", "CD27", "all", "Plasma (%)", gate = "hi_hi"),
      mk("Activated_B", "CD86", "CD80", "CD19", "Activated B in CD19+ (%)", gate = "hi_hi")
    ))
  }
  list(
    mk("T", "CD3", "NK1.1", "all", "T cell (%)", gate = "quad", x_hi = TRUE, y_hi = FALSE),
    mk("B", "CD19", "CD3", "all", "B cell (%)", gate = "quad", x_hi = TRUE, y_hi = FALSE),
    mk("NK", "CD3", "NK1.1", "all", "NK cell (%)", gate = "quad", x_hi = FALSE, y_hi = TRUE),
    mk("Neutrophil", "LY6G", "CD11B", "CD11B", "Neutrophil in CD11B+ (%)", gate = "hi_x", x_hi = TRUE),
    mk("Eosinophil", "Siglec-F", "CCR3", "CD11B", "Eosinophil in CD11B+ (%)", gate = "hi_hi"),
    mk("Mast", "FceRI", "CD200R3", "CD11B", "Mast in CD11B+ (%)", gate = "hi_hi"),
    mk("Macrophage", "F4/80", "LY6C", "CD11B", "Tissue macrophage in CD11B+ (%)", gate = "quad", x_hi = TRUE, y_hi = FALSE),
    mk("M1_like_Mac", "iNOS", "ARG-1", "Macrophage", "M1-like Mac in macrophages (%)", gate = "quad", x_hi = TRUE, y_hi = FALSE),
    mk("M2_like_Mac", "CD206", "ARG-1", "Macrophage", "M2-like Mac in macrophages (%)", gate = "hi_hi"),
    mk("DC", "CD11C", "I-A/I-E", "Myeloid", "DC (CD11C+ MHC-II+) (%)", gate = "hi_hi"),
    mk("cDC1_CD103", "CD103", "CD11C", "DC", "cDC1 in DC (%)", gate = "hi_hi"),
    mk("cDC2", "CD103", "CD11C", "DC", "cDC2 in DC (%)", gate = "quad", x_hi = FALSE, y_hi = TRUE),
    mk("Mono_Ly6Chi", "LY6C", "F4/80", "CD11B", "Inflammatory mono in CD11B+ (%)", gate = "quad", x_hi = TRUE, y_hi = FALSE),
    mk("Mono_Ly6Clo", "LY6C", "CD11B", "CD11B", "Ly6C lo mono in CD11B+ (%)", gate = "quad", x_hi = FALSE, y_hi = TRUE)
  )
}

subset_sample_percent <- function(cells, spec) {
  par <- parent_mask(cells, spec$parent)
  hit <- subset_hit_mask(cells, spec)
  smp <- unique(as.character(cells$sample))
  rows <- lapply(smp, function(s) {
    keep_s <- as.character(cells$sample) == s
    n_par <- sum(keep_s & par)
    n_hit <- sum(keep_s & par & hit)
    grp <- as.character(cells$group[keep_s][1])
    bio <- if ("bio_sample" %in% names(cells)) as.character(cells$bio_sample[keep_s][1]) else s
    if (is.na(bio) || !nzchar(bio)) bio <- s
    tech <- if ("tech_rep" %in% names(cells)) as.character(cells$tech_rep[keep_s][1]) else NA_character_
    data.frame(
      sample = s, bio_sample = bio, tech_rep = tech, group = grp,
      n_parent = n_par, n_subset = n_hit,
      percent = if (n_par > 0) 100 * n_hit / n_par else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# 完整区域门：铺到坐标轴，不要用命中细胞的 10–90% 分位数小框。
# x_hi/y_hi: TRUE=阳性侧, FALSE=阴性侧, NA=该轴整段都要（半平面门，如 TEM=CD62L-）。
complete_gate_rect <- function(xlim, ylim, xcut, ycut, x_hi, y_hi) {
  xmin <- xlim[1]
  xmax <- xlim[2]
  ymin <- ylim[1]
  ymax <- ylim[2]
  if (isTRUE(x_hi)) xmin <- xcut
  if (identical(x_hi, FALSE)) xmax <- xcut
  if (isTRUE(y_hi)) ymin <- ycut
  if (identical(y_hi, FALSE)) ymax <- ycut
  list(
    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
    xcut = xcut, ycut = ycut
  )
}

quad_gate_rect <- function(xlim, ylim, xcut, ycut, x_hi, y_hi) {
  complete_gate_rect(xlim, ylim, xcut, ycut, x_hi, y_hi)
}

complete_gate_for <- function(x, y, spec, xlim, ylim) {
  xcut <- axis_pos_cut(x)
  ycut <- axis_pos_cut(y)
  kind <- spec$gate
  x_hi <- spec$x_hi
  y_hi <- spec$y_hi
  if (identical(kind, "half_x")) {
    y_hi <- NA
    if (length(x_hi) != 1L || is.na(x_hi)) x_hi <- FALSE
  } else if (identical(kind, "half_y")) {
    x_hi <- NA
    if (length(y_hi) != 1L || is.na(y_hi)) y_hi <- TRUE
  } else if (identical(kind, "hi_x")) {
    y_hi <- NA
    if (length(x_hi) != 1L || is.na(x_hi)) x_hi <- TRUE
  } else if (identical(kind, "hi_y")) {
    x_hi <- NA
    if (length(y_hi) != 1L || is.na(y_hi)) y_hi <- TRUE
  } else if (identical(kind, "hi_hi") || identical(kind, "box") || is.null(kind)) {
    if (length(x_hi) != 1L || is.na(x_hi)) x_hi <- TRUE
    if (length(y_hi) != 1L || is.na(y_hi)) y_hi <- TRUE
  }
  complete_gate_rect(xlim, ylim, xcut, ycut, x_hi, y_hi)
}

is_cd62_cd44_spec <- function(spec) {
  grepl("CD62L", spec$x, ignore.case = TRUE) && grepl("CD44", spec$y, ignore.case = TRUE)
}

quadrant_pct_labels <- function(x, y, xcut, ycut, xlim, ylim) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (!n || !is.finite(xcut) || !is.finite(ycut)) return(NULL)
  xx <- x[ok]
  yy <- y[ok]
  quads <- list(
    list(lab = "naive", hit = xx >= xcut & yy < ycut, px = (xcut + xlim[2]) / 2, py = (ylim[1] + ycut) / 2),
    list(lab = "T_CM", hit = xx >= xcut & yy >= ycut, px = (xcut + xlim[2]) / 2, py = (ycut + ylim[2]) / 2),
    list(lab = "T_EM", hit = xx < xcut & yy >= ycut, px = (xlim[1] + xcut) / 2, py = (ycut + ylim[2]) / 2),
    list(lab = "DN/T_EM", hit = xx < xcut & yy < ycut, px = (xlim[1] + xcut) / 2, py = (ylim[1] + ycut) / 2)
  )
  rows <- lapply(quads, function(q) {
    data.frame(
      x = q$px, y = q$py,
      label = sprintf("%s\n%.1f%%", q$lab, 100 * mean(q$hit)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

plot_subset_stat_bar <- function(sample_df, ylab, pval) {
  sample_df <- sample_df[is.finite(sample_df$percent), , drop = FALSE]
  sample_df$group <- factor(as.character(sample_df$group), levels = flow_group_levels)
  glev <- flow_group_levels
  means <- data.frame(
    group = factor(glev, levels = glev),
    mean = vapply(glev, function(g) mean(sample_df$percent[as.character(sample_df$group) == g], na.rm = TRUE), numeric(1)),
    sd = vapply(glev, function(g) {
      v <- sample_df$percent[as.character(sample_df$group) == g]
      if (length(v) < 2) 0 else stats::sd(v)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  means$mean[!is.finite(means$mean)] <- 0
  means$sd[!is.finite(means$sd)] <- 0
  y_top <- max(c(sample_df$percent, means$mean + means$sd), na.rm = TRUE)
  if (!is.finite(y_top) || y_top <= 0) y_top <- 1
  star <- p_annot_label(pval)
  ggplot2::ggplot() +
    ggplot2::geom_col(
      data = means, ggplot2::aes(x = group, y = mean, fill = group),
      width = 0.55, color = NA
    ) +
    ggplot2::geom_errorbar(
      data = means, ggplot2::aes(x = group, ymin = pmax(0, mean - sd), ymax = mean + sd),
      width = 0.16, linewidth = 0.45
    ) +
    ggplot2::geom_point(
      data = sample_df,
      ggplot2::aes(x = group, y = percent, color = group, shape = group),
      fill = "white", size = 2.7, stroke = 0.95,
      position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1)
    ) +
    ggplot2::annotate("segment", x = 1, xend = 2, y = y_top * 1.10, yend = y_top * 1.10, linewidth = 0.4) +
    ggplot2::annotate("text", x = 1.5, y = y_top * 1.18, label = star, fontface = "bold", size = 3.6) +
    ggplot2::scale_fill_manual(values = pal_group, drop = FALSE) +
    ggplot2::scale_color_manual(values = pal_group, drop = FALSE) +
    ggplot2::scale_shape_manual(values = setNames(c(21, 22), glev), drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.28))) +
    ggplot2::coord_cartesian(ylim = c(0, y_top * 1.36)) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(y = ylab)
}

plot_subset_contour <- function(df, x, y, color, xlab, ylab, pct, gate,
                                xlim = NULL, ylim = NULL, quad_labs = NULL,
                                highlight = TRUE) {
  d <- df[is.finite(df[[x]]) & is.finite(df[[y]]), c(x, y), drop = FALSE]
  names(d) <- c("x", "y")
  if (nrow(d) > 4000) {
    set.seed(if (exists("seed_value")) seed_value else 42)
    d <- d[sample.int(nrow(d), 4000), , drop = FALSE]
  }
  extra <- c(gate$xmin, gate$xmax, gate$ymin, gate$ymax, gate$xcut, gate$ycut)
  if (is.null(xlim)) xlim <- finite_axis_lim(c(d$x, extra))
  if (is.null(ylim)) ylim <- finite_axis_lim(c(d$y, extra))
  p <- ggplot2::ggplot()
  pc <- if (nrow(d) >= 40) flow_prob_contour(d$x, d$y) else NULL
  if (!is.null(pc) && length(pc$levels) >= 2 && nrow(pc$grid) > 10) {
    out <- data.frame(x = pc$x[pc$outlier], y = pc$y[pc$outlier])
    if (nrow(out) > 2500) {
      set.seed(if (exists("seed_value")) seed_value else 42)
      out <- out[sample.int(nrow(out), 2500), , drop = FALSE]
    }
    if (nrow(out)) {
      p <- p + ggplot2::geom_point(
        data = out, ggplot2::aes(x = x, y = y),
        size = 0.22, alpha = 0.45, color = color, stroke = 0, shape = 16
      )
    }
    p <- p + ggplot2::geom_contour(
      data = pc$grid, ggplot2::aes(x = x, y = y, z = z),
      breaks = pc$levels, color = color, linewidth = 0.38, lineend = "round"
    )
  } else {
    p <- p + ggplot2::geom_point(
      data = d, ggplot2::aes(x = x, y = y),
      size = 0.28, alpha = 0.35, color = color, stroke = 0
    )
  }
  if (!is.null(gate$xcut) && length(gate$xcut) && is.finite(gate$xcut[1])) {
    p <- p + ggplot2::annotate("segment", x = gate$xcut, xend = gate$xcut,
                               y = ylim[1], yend = ylim[2],
                               color = color, linewidth = 0.35, linetype = "dashed")
  }
  if (!is.null(gate$ycut) && length(gate$ycut) && is.finite(gate$ycut[1])) {
    p <- p + ggplot2::annotate("segment", x = xlim[1], xend = xlim[2],
                               y = gate$ycut, yend = gate$ycut,
                               color = color, linewidth = 0.35, linetype = "dashed")
  }
  if (isTRUE(highlight) && is.finite(gate$xmin) && is.finite(gate$xmax) &&
      is.finite(gate$ymin) && is.finite(gate$ymax)) {
    p <- p + ggplot2::annotate(
      "rect", xmin = gate$xmin, xmax = gate$xmax, ymin = gate$ymin, ymax = gate$ymax,
      color = color, fill = NA, linewidth = 0.55
    )
  }
  if (!is.null(quad_labs) && nrow(quad_labs)) {
    p <- p + ggplot2::geom_text(
      data = quad_labs, ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE, color = color, fontface = "bold", size = 2.7, lineheight = 0.92
    )
  } else if (is.finite(pct)) {
    p <- p + ggplot2::annotate(
      "text", x = if (is.finite(gate$xmax)) gate$xmax else xlim[2],
      y = if (is.finite(gate$ymax)) gate$ymax else ylim[2],
      label = sprintf("%.2f%%", pct), hjust = 1.08, vjust = -0.35,
      color = color, fontface = "bold", size = 3.4
    )
  }
  p +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE, clip = "off") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 10, 4, 4)
    ) +
    ggplot2::labs(x = xlab, y = ylab)
}

save_subset_figure <- function(bar, c_ctrl, c_trt, path_stub) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  draw_one <- function() {
    grid::grid.newpage()
    lay <- grid::grid.layout(2, 2, heights = grid::unit(c(1.05, 1.45), "null"))
    grid::pushViewport(grid::viewport(layout = lay))
    print(bar, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2), newpage = FALSE)
    print(c_ctrl, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1), newpage = FALSE)
    print(c_trt, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2), newpage = FALSE)
  }
  tryCatch({
    grDevices::pdf(paste0(path_stub, ".pdf"), width = 6.6, height = 7.2)
    draw_one()
    grDevices::dev.off()
  }, error = function(e) log_msg("subset pdf failed: ", e$message))
  tryCatch({
    grDevices::png(paste0(path_stub, ".png"), width = 6.6, height = 7.2, units = "in", res = 300)
    draw_one()
    grDevices::dev.off()
  }, error = function(e) log_msg("subset png failed: ", e$message))
}

export_subset_gate_figures <- function(cells, panel_id, out_dir) {
  specs <- subset_plot_specs(panel_id)
  sub_dir <- file.path(out_dir, "subset_stats")
  dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
  stat_rows <- list()
  for (spec in specs) {
    if (skip_if_missing_channels(
      c(spec$x, spec$y), names(cells),
      paste0(panel_id, " ", spec$lineage, " subset plot")
    )) next
    tryCatch({
      par <- parent_mask(cells, spec$parent)
      hit <- subset_hit_mask(cells, spec)
      par[is.na(par)] <- FALSE
      hit[is.na(hit)] <- FALSE
      n_par <- sum(par)
      n_hit <- sum(par & hit)
      if (!is.finite(n_par) || !is.finite(n_hit) || n_par < 20 || n_hit < 8) {
        # 细胞不够：不分析这一项，继续后面的亚群
      } else {
        samp <- subset_sample_percent(cells, spec)
        if (is.null(samp) || nrow(samp) < 2) {
          # 无样品百分比：跳过该项
        } else {
          samp_bio <- bio_percent_table(samp)
          if (is.null(samp_bio) || nrow(samp_bio) < 2) {
            # 无生物学重复表：跳过该项
          } else {
            ctrl_v <- samp_bio$percent[as.character(samp_bio$group) == flow_ctrl_group]
            trt_v <- samp_bio$percent[as.character(samp_bio$group) == flow_trt_group]
            dropped_bio <- attr(samp_bio, "dropped")
            pv <- if (length(ctrl_v) >= 2 && length(trt_v) >= 2) {
              tryCatch(stats::t.test(trt_v, ctrl_v)$p.value, error = function(e) NA_real_)
            } else {
              NA_real_
            }
            lab <- if (isTRUE(spec$use_major)) {
              celltype_label(spec$lineage, panel_id)
            } else {
              celltype_label(spec$lineage, panel_id)
            }
            ylab <- if (!is.null(spec$ylab) && nzchar(spec$ylab)) spec$ylab else paste0(lab, " (%)")
            bar <- plot_subset_stat_bar(samp_bio, ylab, pv)
            d_ctrl <- cells[par & as.character(cells$group) == flow_ctrl_group, , drop = FALSE]
            d_trt <- cells[par & as.character(cells$group) == flow_trt_group, , drop = FALSE]
            if (nrow(d_ctrl) < 8 || nrow(d_trt) < 8) {
              # EV/H 细胞不够：跳过该项
            } else {
              xlim <- finite_axis_lim(c(d_ctrl[[spec$x]], d_trt[[spec$x]]))
              ylim <- finite_axis_lim(c(d_ctrl[[spec$y]], d_trt[[spec$y]]))
              # EV、H 各自用本组细胞切完整门，不要六个样品共用一条线，也不要框命中细胞的分位数小盒
              gate_ctrl <- complete_gate_for(d_ctrl[[spec$x]], d_ctrl[[spec$y]], spec, xlim, ylim)
              gate_trt <- complete_gate_for(d_trt[[spec$x]], d_trt[[spec$y]], spec, xlim, ylim)
              xlab <- axis_fl_label(panel_id, spec$x)
              yfl <- axis_fl_label(panel_id, spec$y)
              pct_of <- function(g) {
                v <- samp_bio$percent[as.character(samp_bio$group) == g]
                if (!length(v) || all(!is.finite(v))) return(0)
                mean(v, na.rm = TRUE)
              }
              col_ctrl <- unname(pal_group[flow_ctrl_group])
              col_trt <- unname(pal_group[flow_trt_group])
              if (is.na(col_ctrl)) col_ctrl <- "#1A1A1A"
              if (is.na(col_trt)) col_trt <- "#E31A1C"
              xlim <- finite_axis_lim(c(
                d_ctrl[[spec$x]], d_trt[[spec$x]],
                gate_ctrl$xmin, gate_ctrl$xmax, gate_ctrl$xcut,
                gate_trt$xmin, gate_trt$xmax, gate_trt$xcut
              ))
              ylim <- finite_axis_lim(c(
                d_ctrl[[spec$y]], d_trt[[spec$y]],
                gate_ctrl$ymin, gate_ctrl$ymax, gate_ctrl$ycut,
                gate_trt$ymin, gate_trt$ymax, gate_trt$ycut
              ))
              labs_ctrl <- if (is_cd62_cd44_spec(spec)) {
                quadrant_pct_labels(d_ctrl[[spec$x]], d_ctrl[[spec$y]], gate_ctrl$xcut, gate_ctrl$ycut, xlim, ylim)
              } else {
                NULL
              }
              labs_trt <- if (is_cd62_cd44_spec(spec)) {
                quadrant_pct_labels(d_trt[[spec$x]], d_trt[[spec$y]], gate_trt$xcut, gate_trt$ycut, xlim, ylim)
              } else {
                NULL
              }
              c_ctrl <- plot_subset_contour(
                d_ctrl, spec$x, spec$y, col_ctrl, xlab, yfl, pct_of(flow_ctrl_group),
                gate_ctrl, xlim, ylim, labs_ctrl
              )
              c_trt <- plot_subset_contour(
                d_trt, spec$x, spec$y, col_trt, xlab, yfl, pct_of(flow_trt_group),
                gate_trt, xlim, ylim, labs_trt
              )
              stub <- paste0(panel_id, "_", gsub("[^A-Za-z0-9]+", "_", spec$lineage), "_JZ_AB_vs_JZ_EVB")
              save_subset_figure(bar, c_ctrl, c_trt, file.path(sub_dir, stub))
              utils::write.csv(samp, file.path(sub_dir, paste0(stub, "_by_sample.csv")), row.names = FALSE)
              utils::write.csv(samp_bio, file.path(sub_dir, paste0(stub, "_by_bio.csv")), row.names = FALSE)
              stat_rows[[length(stat_rows) + 1]] <- data.frame(
                panel = panel_id,
                subset = spec$lineage,
                celltype = lab,
                parent = spec$parent,
                n_EV = length(ctrl_v),
                n_H = length(trt_v),
                mean_EV = mean(ctrl_v, na.rm = TRUE),
                mean_H = mean(trt_v, na.rm = TRUE),
                sd_EV = stats::sd(ctrl_v),
                sd_H = stats::sd(trt_v),
                delta_H_minus_EV = mean(trt_v, na.rm = TRUE) - mean(ctrl_v, na.rm = TRUE),
                p_value = pv,
                dropped_EV = if (!is.null(dropped_bio) && nrow(dropped_bio)) {
                  paste(dropped_bio$dropped_bio[dropped_bio$group == flow_ctrl_group], collapse = ",")
                } else "",
                dropped_H = if (!is.null(dropped_bio) && nrow(dropped_bio)) {
                  paste(dropped_bio$dropped_bio[dropped_bio$group == flow_trt_group], collapse = ",")
                } else "",
                ylab = ylab,
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
    }, error = function(e) {
      log_msg(panel_id, " skip ", spec$lineage, " subset plot: ", e$message)
    })
  }
  if (length(stat_rows)) {
    st <- do.call(rbind, stat_rows)
    st$padj <- if (all(is.na(st$p_value))) NA_real_ else p.adjust(st$p_value, method = "BH")
    utils::write.csv(st, file.path(sub_dir, paste0(panel_id, "_subset_JZ_AB_vs_JZ_EVB_stats.csv")), row.names = FALSE)
  }
  log_msg(panel_id, " subset stat+contour figures: ", sub_dir)
  invisible(TRUE)
}

save_gating_plot <- function(plot, path_stub, width = 4.6, height = 4.4) {
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    grDevices::pdf(paste0(path_stub, ".pdf"), width = width, height = height)
    print(plot)
    grDevices::dev.off()
  }, error = function(e) log_msg("gating pdf failed: ", e$message))
  tryCatch({
    grDevices::png(paste0(path_stub, ".png"), width = width, height = height, units = "in", res = 300)
    print(plot)
    grDevices::dev.off()
  }, error = function(e) log_msg("gating png failed: ", e$message))
}

# 每个 FCS 单独出完整 2D 门：十字线铺满坐标轴，CD62L/CD44 标四个象限，不要 10–90% 小框
export_per_sample_gating_figures <- function(cells, panel_id, out_dir) {
  specs <- subset_plot_specs(panel_id)
  key <- vapply(specs, function(s) paste(s$parent, s$x, s$y, sep = "\t"), character(1))
  views <- specs[!duplicated(key)]
  views <- Filter(function(spec) {
    !skip_if_missing_channels(
      c(spec$x, spec$y), names(cells),
      paste0(panel_id, " ", spec$parent, " ", spec$x, "/", spec$y, " gating")
    )
  }, views)
  gate_dir <- file.path(out_dir, "gating")
  dir.create(gate_dir, recursive = TRUE, showWarnings = FALSE)
  smp <- unique(as.character(cells$sample))
  cut_rows <- list()
  for (s in smp) {
    keep_s <- as.character(cells$sample) == s
    grp <- as.character(cells$group[keep_s][1])
    col <- unname(pal_group[grp])
    if (is.na(col) || !nzchar(col)) col <- if (identical(grp, flow_trt_group)) "#E31A1C" else "#1A1A1A"
    s_dir <- file.path(gate_dir, s)
    dir.create(s_dir, recursive = TRUE, showWarnings = FALSE)
    for (spec in views) {
      tryCatch({
        par <- parent_mask(cells, spec$parent)
        par[is.na(par)] <- FALSE
        d <- cells[keep_s & par, , drop = FALSE]
        if (nrow(d) < 20) {
          # 该样品该门细胞不够：跳过该项
        } else {
          xlim <- finite_axis_lim(d[[spec$x]])
          ylim <- finite_axis_lim(d[[spec$y]])
          gate <- complete_gate_for(d[[spec$x]], d[[spec$y]], spec, xlim, ylim)
          xlim <- finite_axis_lim(c(d[[spec$x]], gate$xmin, gate$xmax, gate$xcut))
          ylim <- finite_axis_lim(c(d[[spec$y]], gate$ymin, gate$ymax, gate$ycut))
          gate <- complete_gate_for(d[[spec$x]], d[[spec$y]], spec, xlim, ylim)
          qlabs <- quadrant_pct_labels(d[[spec$x]], d[[spec$y]], gate$xcut, gate$ycut, xlim, ylim)
          show_quads <- is_cd62_cd44_spec(spec) || identical(spec$gate, "quad")
          xlab <- axis_fl_label(panel_id, spec$x)
          ylab <- axis_fl_label(panel_id, spec$y)
          pct <- if (isTRUE(spec$use_major) && "cluster_lineage" %in% names(d)) {
            100 * mean(as.character(d$cluster_lineage) == spec$lineage, na.rm = TRUE)
          } else if ("lineage" %in% names(d)) {
            100 * mean(as.character(d$lineage) == spec$lineage, na.rm = TRUE)
          } else {
            NA_real_
          }
          p <- plot_subset_contour(
            d, spec$x, spec$y, col, xlab, ylab, pct, gate, xlim, ylim,
            quad_labs = if (show_quads) qlabs else NULL,
            highlight = !show_quads
          )
          stub <- paste0(
            panel_id, "_", s, "_",
            gsub("[^A-Za-z0-9]+", "_", paste(spec$parent, spec$x, spec$y, sep = "_"))
          )
          save_gating_plot(p, file.path(s_dir, stub))
          cut_rows[[length(cut_rows) + 1]] <- data.frame(
            panel = panel_id, sample = s, group = grp,
            parent = spec$parent, marker_x = spec$x, marker_y = spec$y,
            xcut = gate$xcut, ycut = gate$ycut,
            xmin = gate$xmin, xmax = gate$xmax, ymin = gate$ymin, ymax = gate$ymax,
            n_parent = nrow(d),
            stringsAsFactors = FALSE
          )
        }
      }, error = function(e) {
        log_msg(panel_id, " skip ", s, " ", spec$lineage, " gating: ", e$message)
      })
    }
  }
  if (length(cut_rows)) {
    utils::write.csv(
      do.call(rbind, cut_rows),
      file.path(gate_dir, paste0(panel_id, "_per_sample_gate_cuts.csv")),
      row.names = FALSE
    )
  }
  log_msg(panel_id, " per-sample complete gates: ", gate_dir)
  invisible(TRUE)
}

# P1 活化读出：NKG2D / IFN-g / TNF-a / GZMB 在 CD4/CD8/NK 上的 MFI，不是新亚群
p1_activation_parent <- function(lin, cl = NULL) {
  lin <- as.character(lin)
  lin[is.na(lin)] <- ""
  out <- rep(NA_character_, length(lin))
  out[grepl("^CD4_", lin) | lin == "Treg"] <- "CD4"
  out[grepl("^CD8_", lin)] <- "CD8"
  out[lin %in% nk_family()] <- "NK"
  if (!is.null(cl)) {
    cl <- as.character(cl)
    cl[is.na(cl)] <- ""
    out[is.na(out) & cl == "CD4"] <- "CD4"
    out[is.na(out) & cl == "CD8"] <- "CD8"
    out[is.na(out) & cl == "NK"] <- "NK"
  }
  out
}

export_p1_activation_stats <- function(cells, out_dir) {
  if (!"lineage" %in% names(cells)) return(invisible(NULL))
  cl <- if ("cluster_lineage" %in% names(cells)) cells$cluster_lineage else NULL
  par <- p1_activation_parent(cells$lineage, cl)
  markers <- intersect(c("NKG2D", "IFN-g", "TNF-a", "GZMB"), names(cells))
  if (!length(markers) || !any(!is.na(par))) return(invisible(NULL))
  parents <- c("CD4", "CD8", "NK")
  cuts <- lapply(markers, function(mk) axis_pos_cut(as.numeric(cells[[mk]][!is.na(par)])))
  names(cuts) <- markers
  smp <- unique(as.character(cells$sample))
  rows <- list()
  for (s in smp) {
    ii <- which(as.character(cells$sample) == s)
    if (!length(ii)) next
    grp <- as.character(cells$group[ii[1]])
    bio <- if ("bio_sample" %in% names(cells)) as.character(cells$bio_sample[ii[1]]) else s
    tech <- if ("tech_rep" %in% names(cells)) as.character(cells$tech_rep[ii[1]]) else NA_character_
    for (pn in parents) {
      hit <- ii[par[ii] == pn & !is.na(par[ii])]
      if (length(hit) < 5) next
      for (mk in markers) {
        v <- as.numeric(cells[[mk]][hit])
        cut <- cuts[[mk]]
        rows[[length(rows) + 1]] <- data.frame(
          sample = s, bio_sample = bio, tech_rep = tech, group = grp,
          parent = pn, marker = mk, n_cells = length(hit),
          MFI = stats::median(v, na.rm = TRUE),
          pct_positive = 100 * mean(is.finite(v) & is.finite(cut) & v >= cut),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) return(invisible(NULL))
  tab <- do.call(rbind, rows)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(tab, file.path(out_dir, "P1_TNK_activation_by_sample.csv"), row.names = FALSE)
  mfi_df <- tab
  mfi_df$id <- paste(mfi_df$parent, mfi_df$marker, "MFI", sep = "|")
  mfi_df$percent <- mfi_df$MFI
  pct_df <- tab
  pct_df$id <- paste(pct_df$parent, pct_df$marker, "pct_pos", sep = "|")
  pct_df$percent <- pct_df$pct_positive
  stats_mfi <- compare_group_freq(mfi_df, "id")
  stats_pct <- compare_group_freq(pct_df, "id")
  stats_mfi$metric <- "MFI"
  stats_pct$metric <- "pct_positive"
  stats <- rbind(stats_mfi, stats_pct)
  utils::write.csv(stats, file.path(out_dir, "P1_TNK_activation_JZ_AB_vs_JZ_EVB_stats.csv"), row.names = FALSE)
  log_msg("P1 NKG2D/IFN-g/TNF-a/GZMB MFI written (not used as subset labels)")
  invisible(tab)
}

# P2 活化：CD40 / CD80 / CD86 是 Naive / Unswitched / Switched 上的 MFI 与阳性率，不是新亚群
export_p2_activation_stats <- function(cells, out_dir) {
  markers <- intersect(c("CD40", "CD80", "CD86"), names(cells))
  parents <- c("Naive_B", "Unswitched_B", "Switched_B")
  if (!length(markers) || !"lineage" %in% names(cells)) return(invisible(NULL))
  maj <- if ("cluster_lineage" %in% names(cells)) {
    as.character(cells$cluster_lineage)
  } else {
    as.character(cells$lineage)
  }
  maj[is.na(maj)] <- ""
  is_b <- maj %in% c(parents, "Atypical_B", "Memory_B")
  cuts <- lapply(markers, function(mk) {
    v <- as.numeric(cells[[mk]][is_b])
    axis_pos_cut(v)
  })
  names(cuts) <- markers
  smp <- unique(as.character(cells$sample))
  rows <- list()
  for (s in smp) {
    ii <- which(as.character(cells$sample) == s)
    if (!length(ii)) next
    grp <- as.character(cells$group[ii[1]])
    bio <- if ("bio_sample" %in% names(cells)) as.character(cells$bio_sample[ii[1]]) else s
    tech <- if ("tech_rep" %in% names(cells)) as.character(cells$tech_rep[ii[1]]) else NA_character_
    for (par in parents) {
      hit <- ii[maj[ii] == par]
      if (length(hit) < 5) next
      for (mk in markers) {
        v <- as.numeric(cells[[mk]][hit])
        cut <- cuts[[mk]]
        rows[[length(rows) + 1]] <- data.frame(
          sample = s, bio_sample = bio, tech_rep = tech, group = grp,
          parent = par, marker = mk, n_cells = length(hit),
          MFI = stats::median(v, na.rm = TRUE),
          pct_positive = 100 * mean(is.finite(v) & is.finite(cut) & v >= cut),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) return(invisible(NULL))
  tab <- do.call(rbind, rows)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(tab, file.path(out_dir, "P2_Bcell_activation_by_sample.csv"), row.names = FALSE)
  mfi_df <- tab
  mfi_df$id <- paste(mfi_df$parent, mfi_df$marker, "MFI", sep = "|")
  mfi_df$percent <- mfi_df$MFI
  pct_df <- tab
  pct_df$id <- paste(pct_df$parent, pct_df$marker, "pct_pos", sep = "|")
  pct_df$percent <- pct_df$pct_positive
  stats_mfi <- compare_group_freq(mfi_df, "id")
  stats_pct <- compare_group_freq(pct_df, "id")
  stats_mfi$metric <- "MFI"
  stats_pct$metric <- "pct_positive"
  stats <- rbind(stats_mfi, stats_pct)
  utils::write.csv(stats, file.path(out_dir, "P2_Bcell_activation_JZ_AB_vs_JZ_EVB_stats.csv"), row.names = FALSE)
  log_msg("P2 CD40/CD80/CD86 MFI and % positivity written (not used as subset labels)")
  invisible(tab)
}

p3_activation_parent <- function(lin) {
  lin <- as.character(lin)
  lin[is.na(lin)] <- ""
  out <- rep(NA_character_, length(lin))
  out[lin == "B"] <- "B"
  out[lin %in% c("T", "CD4", "CD8")] <- "T"
  out[lin %in% c("Macrophage", "M1_like_Mac", "M2_like_Mac")] <- "Macrophage"
  out[lin %in% c("DC", "cDC1_CD103", "cDC2")] <- "DC"
  out
}

# P3 活化：CD40/CD80/CD86 在 B/DC/巨噬上的 MFI；TNF-a 在巨噬/DC/T。都不是新亚群
export_p3_activation_stats <- function(cells, out_dir) {
  if (!"lineage" %in% names(cells)) return(invisible(NULL))
  lin <- as.character(cells$lineage)
  par <- p3_activation_parent(lin)
  markers <- intersect(c("CD40", "CD80", "CD86", "TNF-a"), names(cells))
  if (!length(markers) || !any(!is.na(par))) return(invisible(NULL))
  allowed <- list(
    B = c("CD40", "CD80", "CD86"),
    DC = c("CD40", "CD80", "CD86", "TNF-a"),
    Macrophage = c("CD40", "CD80", "CD86", "TNF-a"),
    T = c("TNF-a")
  )
  cuts <- lapply(markers, function(mk) axis_pos_cut(as.numeric(cells[[mk]][!is.na(par)])))
  names(cuts) <- markers
  smp <- unique(as.character(cells$sample))
  rows <- list()
  for (s in smp) {
    ii <- which(as.character(cells$sample) == s)
    if (!length(ii)) next
    grp <- as.character(cells$group[ii[1]])
    bio <- if ("bio_sample" %in% names(cells)) as.character(cells$bio_sample[ii[1]]) else s
    tech <- if ("tech_rep" %in% names(cells)) as.character(cells$tech_rep[ii[1]]) else NA_character_
    for (pn in names(allowed)) {
      hit <- ii[par[ii] == pn & !is.na(par[ii])]
      if (length(hit) < 5) next
      for (mk in intersect(allowed[[pn]], markers)) {
        v <- as.numeric(cells[[mk]][hit])
        cut <- cuts[[mk]]
        rows[[length(rows) + 1]] <- data.frame(
          sample = s, bio_sample = bio, tech_rep = tech, group = grp,
          parent = pn, marker = mk, n_cells = length(hit),
          MFI = stats::median(v, na.rm = TRUE),
          pct_positive = 100 * mean(is.finite(v) & is.finite(cut) & v >= cut),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) return(invisible(NULL))
  tab <- do.call(rbind, rows)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(tab, file.path(out_dir, "P3_APC_activation_by_sample.csv"), row.names = FALSE)
  mfi_df <- tab
  mfi_df$id <- paste(mfi_df$parent, mfi_df$marker, "MFI", sep = "|")
  mfi_df$percent <- mfi_df$MFI
  pct_df <- tab
  pct_df$id <- paste(pct_df$parent, pct_df$marker, "pct_pos", sep = "|")
  pct_df$percent <- pct_df$pct_positive
  stats_mfi <- compare_group_freq(mfi_df, "id")
  stats_pct <- compare_group_freq(pct_df, "id")
  stats_mfi$metric <- "MFI"
  stats_pct$metric <- "pct_positive"
  stats <- rbind(stats_mfi, stats_pct)
  utils::write.csv(stats, file.path(out_dir, "P3_APC_activation_JZ_AB_vs_JZ_EVB_stats.csv"), row.names = FALSE)
  log_msg("P3 CD40/CD80/CD86/TNF-a MFI written (not used as subset labels)")
  invisible(tab)
}

# -----------------------------------------------------------------------------
# NKT / B 亚群单独的活化、效应、耗竭（图1：上柱状图，下对照/处理等高线）
# B 的 P2 没有 PD-L1/LAG-3/TIM-3：只报 CD86/CD80/CD40，缺通道 skip，不伪造抑制门。
# -----------------------------------------------------------------------------
flow_comparison_tag <- function() {
  trt <- gsub("[^A-Za-z0-9]+", "_", as.character(flow_trt_group)[1])
  ctrl <- gsub("[^A-Za-z0-9]+", "_", as.character(flow_ctrl_group)[1])
  paste0(trt, "_vs_", ctrl)
}

func_group_short <- function(g) {
  rec <- c(
    EV = "EV", H = "H",
    "JY-EVNK" = "EVNK", "JY-NNK" = "NNK",
    "JZ-EVB" = "EVB", "JZ-AB" = "AB"
  )
  g <- as.character(g)[1]
  if (g %in% names(rec)) rec[[g]] else g
}

marker_pretty_label <- function(mk) {
  rec <- c(
    "IFN-g" = "IFN-\u03B3", "TNF-a" = "TNF-\u03B1", GZMB = "Granzyme B",
    Perforin = "Perforin", CD69 = "CD69", "PD-L1" = "PD-L1",
    "LAG-3" = "LAG-3", "TIM-3" = "TIM-3", CD80 = "CD80", CD86 = "CD86",
    CD40 = "CD40", NKG2D = "NKG2D"
  )
  mk <- as.character(mk)[1]
  if (mk %in% names(rec)) rec[[mk]] else mk
}

func_parent_pretty <- function(p) {
  rec <- c(
    NKT = "NKT cell", Naive_B = "Naive B", Unswitched_B = "Unswitched B",
    Switched_B = "Switched B", Atypical_B = "Atypical B",
    Activated_B = "Activated B"
  )
  p <- as.character(p)[1]
  if (p %in% names(rec)) rec[[p]] else p
}

functional_parent_mask <- function(cells, parent) {
  n <- nrow(cells)
  if (n < 1) return(logical(0))
  cl <- if ("cluster_lineage" %in% names(cells)) as.character(cells$cluster_lineage) else rep("", n)
  lin <- if ("lineage" %in% names(cells)) as.character(cells$lineage) else rep("", n)
  cl[is.na(cl)] <- ""
  lin[is.na(lin)] <- ""
  out <- switch(
    as.character(parent),
    NKT = cl == "NKT" | lin %in% nkt_family() | grepl("^NKT", lin),
    Naive_B = lin == "Naive_B" | (cl == "Naive_B" & !lin %in% c("MZ_B", "Activated_B")),
    Unswitched_B = lin == "Unswitched_B" | (cl == "Unswitched_B" & !lin %in% c("MZ_B", "Activated_B")),
    Switched_B = lin == "Switched_B" | (cl == "Switched_B" & !lin %in% c("Plasmablast", "Activated_B")),
    Atypical_B = lin == "Atypical_B" | cl == "Atypical_B",
    Activated_B = lin == "Activated_B" | cl == "Activated_B",
    parent_mask(cells, parent)
  )
  out[is.na(out)] <- FALSE
  out
}

func_x_marker <- function(parent, available) {
  if (identical(as.character(parent), "NKT")) {
    for (m in c("NKp46", "NK1.1", "CD3")) if (m %in% available) return(m)
    return(NA_character_)
  }
  for (m in c("CD19", "IgD", "CD27")) if (m %in% available) return(m)
  NA_character_
}

functional_state_specs <- function(panel_id) {
  if (identical(panel_id, "P1")) {
    return(list(
      list(parent = "NKT", state = "activation_effector",
           markers = c("CD69", "IFN-g", "TNF-a", "GZMB")),
      list(parent = "NKT", state = "exhaustion",
           markers = c("PD-L1", "LAG-3", "TIM-3"))
    ))
  }
  if (identical(panel_id, "P2")) {
    return(lapply(
      c("Naive_B", "Unswitched_B", "Switched_B", "Atypical_B", "Activated_B"),
      function(p) list(parent = p, state = "activation", markers = c("CD86", "CD80", "CD40"))
    ))
  }
  list()
}

functional_marker_sample_percent <- function(cells, parent, marker) {
  par <- functional_parent_mask(cells, parent)
  smp <- unique(as.character(cells$sample))
  rows <- lapply(smp, function(s) {
    keep_s <- as.character(cells$sample) == s
    hit <- keep_s & par
    v <- as.numeric(cells[[marker]][hit])
    n_par <- sum(hit)
    cut <- if (n_par >= 8) axis_pos_cut(v) else NA_real_
    n_pos <- if (is.finite(cut)) sum(is.finite(v) & v >= cut) else 0L
    grp <- as.character(cells$group[keep_s][1])
    bio <- if ("bio_sample" %in% names(cells)) as.character(cells$bio_sample[keep_s][1]) else s
    if (is.na(bio) || !nzchar(bio)) bio <- s
    tech <- if ("tech_rep" %in% names(cells)) as.character(cells$tech_rep[keep_s][1]) else NA_character_
    data.frame(
      sample = s, bio_sample = bio, tech_rep = tech, group = grp,
      parent = parent, marker = marker, n_parent = n_par, n_pos = n_pos,
      percent = if (n_par > 0 && is.finite(cut)) 100 * n_pos / n_par else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

plot_func_state_bar <- function(sample_df, ylab, pval) {
  sample_df <- sample_df[is.finite(sample_df$percent), , drop = FALSE]
  glev <- flow_group_levels
  glab <- unname(vapply(glev, func_group_short, character(1)))
  sample_df$g_lab <- factor(
    vapply(as.character(sample_df$group), func_group_short, character(1)),
    levels = glab
  )
  pal_s <- pal_group
  names(pal_s) <- glab
  means <- data.frame(
    g_lab = factor(glab, levels = glab),
    mean = vapply(glev, function(g) {
      mean(sample_df$percent[as.character(sample_df$group) == g], na.rm = TRUE)
    }, numeric(1)),
    sd = vapply(glev, function(g) {
      v <- sample_df$percent[as.character(sample_df$group) == g]
      if (length(v) < 2) 0 else stats::sd(v)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  means$mean[!is.finite(means$mean)] <- 0
  means$sd[!is.finite(means$sd)] <- 0
  y_top <- max(c(sample_df$percent, means$mean + means$sd), na.rm = TRUE)
  if (!is.finite(y_top) || y_top <= 0) y_top <- 1
  star <- p_annot_label(pval)
  ggplot2::ggplot() +
    ggplot2::geom_col(
      data = means, ggplot2::aes(x = g_lab, y = mean, color = g_lab),
      fill = "white", width = 0.52, linewidth = 0.9
    ) +
    ggplot2::geom_errorbar(
      data = means, ggplot2::aes(x = g_lab, ymin = pmax(0, mean - sd), ymax = mean + sd, color = g_lab),
      width = 0.15, linewidth = 0.45, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = sample_df,
      ggplot2::aes(x = g_lab, y = percent, color = g_lab, fill = g_lab),
      shape = 21, size = 2.5, stroke = 0.35,
      position = ggplot2::position_jitter(width = 0.07, height = 0, seed = 1)
    ) +
    ggplot2::annotate("segment", x = 1, xend = 2, y = y_top * 1.10, yend = y_top * 1.10, linewidth = 0.4) +
    ggplot2::annotate("text", x = 1.5, y = y_top * 1.20, label = star, fontface = "bold", size = 3.5) +
    ggplot2::scale_color_manual(values = pal_s, drop = FALSE) +
    ggplot2::scale_fill_manual(values = pal_s, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.30))) +
    ggplot2::coord_cartesian(ylim = c(0, y_top * 1.40)) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black"),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.4),
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 6, 2, 8)
    ) +
    ggplot2::labs(y = ylab)
}

plot_func_state_contour <- function(df, x, y, color, xlab, ylab, pct, gate,
                                    xlim, ylim, group_tag = NULL, show_x = TRUE) {
  gate$xcut <- NA_real_
  p <- plot_subset_contour(df, x, y, color, xlab, ylab, pct, gate, xlim, ylim, NULL, TRUE)
  if (!is.null(group_tag) && nzchar(group_tag)) {
    p <- p + ggplot2::labs(tag = group_tag) +
      ggplot2::theme(
        plot.tag = ggplot2::element_text(face = "bold", size = 11, color = color),
        plot.tag.position = "left"
      )
  }
  if (!isTRUE(show_x)) {
    p <- p + ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  }
  p + ggplot2::theme(
    axis.text = ggplot2::element_text(color = "black"),
    axis.title.y = ggplot2::element_text(color = "black"),
    plot.margin = ggplot2::margin(4, 8, 6, 10)
  )
}

save_func_state_figure <- function(columns, path_stub) {
  n <- length(columns)
  if (!n) return(invisible(FALSE))
  dir.create(dirname(path_stub), recursive = TRUE, showWarnings = FALSE)
  w <- max(2.9 * n + 0.4, 7.4)
  h <- 8.6
  draw_one <- function() {
    grid::grid.newpage()
    lay <- grid::grid.layout(3, n, heights = grid::unit(c(1.05, 1.28, 1.42), "null"))
    grid::pushViewport(grid::viewport(layout = lay, width = 0.98, height = 0.97))
    for (i in seq_len(n)) {
      print(columns[[i]]$bar, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i), newpage = FALSE)
      print(columns[[i]]$c_ctrl, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = i), newpage = FALSE)
      print(columns[[i]]$c_trt, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = i), newpage = FALSE)
    }
  }
  tryCatch({
    grDevices::pdf(paste0(path_stub, ".pdf"), width = w, height = h)
    draw_one()
    grDevices::dev.off()
  }, error = function(e) log_msg("functional-state pdf failed: ", e$message))
  tryCatch({
    grDevices::png(paste0(path_stub, ".png"), width = w, height = h, units = "in", res = 300)
    draw_one()
    grDevices::dev.off()
  }, error = function(e) log_msg("functional-state png failed: ", e$message))
  invisible(TRUE)
}

export_one_functional_state <- function(cells, panel_id, spec, out_dir) {
  parent <- spec$parent
  markers <- intersect(as.character(spec$markers), names(cells))
  xmk <- func_x_marker(parent, names(cells))
  label <- paste0(panel_id, " ", parent, " ", spec$state)
  if (!length(markers) || is.na(xmk) || !nzchar(xmk)) {
    skip_if_missing_channels(c(spec$markers, "lineage-x"), names(cells), label)
    return(invisible(NULL))
  }
  if (skip_if_missing_channels(c(xmk, markers), names(cells), label)) {
    return(invisible(NULL))
  }
  par <- functional_parent_mask(cells, parent)
  if (sum(par) < 20) {
    log_msg("skip ", label, " (parent n=", sum(par), ")")
    return(invisible(NULL))
  }
  d_ctrl <- cells[par & as.character(cells$group) == flow_ctrl_group, , drop = FALSE]
  d_trt <- cells[par & as.character(cells$group) == flow_trt_group, , drop = FALSE]
  if (nrow(d_ctrl) < 8 || nrow(d_trt) < 8) {
    log_msg("skip ", label, " (too few cells in one group)")
    return(invisible(NULL))
  }
  col_ctrl <- unname(pal_group[flow_ctrl_group])
  col_trt <- unname(pal_group[flow_trt_group])
  if (is.na(col_ctrl)) col_ctrl <- "#1A1A1A"
  if (is.na(col_trt)) col_trt <- "#E31A1C"
  tag_ctrl <- func_group_short(flow_ctrl_group)
  tag_trt <- func_group_short(flow_trt_group)
  parent_lab <- func_parent_pretty(parent)
  xlab <- axis_fl_label(panel_id, xmk)
  columns <- list()
  stat_rows <- list()
  hi_y_spec <- list(gate = "hi_y", x_hi = NA, y_hi = TRUE)
  for (mk in markers) {
    samp <- functional_marker_sample_percent(cells, parent, mk)
    if (is.null(samp) || nrow(samp) < 2) next
    samp_bio <- bio_percent_table(samp)
    if (is.null(samp_bio) || nrow(samp_bio) < 2) next
    ctrl_v <- samp_bio$percent[as.character(samp_bio$group) == flow_ctrl_group]
    trt_v <- samp_bio$percent[as.character(samp_bio$group) == flow_trt_group]
    pv <- if (length(ctrl_v) >= 2 && length(trt_v) >= 2) {
      tryCatch(stats::t.test(trt_v, ctrl_v)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    }
    ylab <- paste0(marker_pretty_label(mk), "\u207A ", parent_lab, " (%)")
    bar <- plot_func_state_bar(samp_bio, ylab, pv)
    yfl <- axis_fl_label(panel_id, mk)
    xlim <- finite_axis_lim(c(d_ctrl[[xmk]], d_trt[[xmk]]))
    ylim <- finite_axis_lim(c(d_ctrl[[mk]], d_trt[[mk]]))
    gate_ctrl <- complete_gate_for(d_ctrl[[xmk]], d_ctrl[[mk]], hi_y_spec, xlim, ylim)
    gate_trt <- complete_gate_for(d_trt[[xmk]], d_trt[[mk]], hi_y_spec, xlim, ylim)
    xlim <- finite_axis_lim(c(
      d_ctrl[[xmk]], d_trt[[xmk]],
      gate_ctrl$xmin, gate_ctrl$xmax, gate_trt$xmin, gate_trt$xmax
    ))
    ylim <- finite_axis_lim(c(
      d_ctrl[[mk]], d_trt[[mk]],
      gate_ctrl$ymin, gate_ctrl$ymax, gate_ctrl$ycut,
      gate_trt$ymin, gate_trt$ymax, gate_trt$ycut
    ))
    pct_of <- function(g) {
      v <- samp_bio$percent[as.character(samp_bio$group) == g]
      if (!length(v) || all(!is.finite(v))) return(0)
      mean(v, na.rm = TRUE)
    }
    c_ctrl <- plot_func_state_contour(
      d_ctrl, xmk, mk, col_ctrl, xlab, yfl, pct_of(flow_ctrl_group),
      gate_ctrl, xlim, ylim, tag_ctrl, show_x = FALSE
    )
    c_trt <- plot_func_state_contour(
      d_trt, xmk, mk, col_trt, xlab, yfl, pct_of(flow_trt_group),
      gate_trt, xlim, ylim, tag_trt, show_x = TRUE
    )
    columns[[length(columns) + 1]] <- list(bar = bar, c_ctrl = c_ctrl, c_trt = c_trt)
    dropped_bio <- attr(samp_bio, "dropped")
    stat_rows[[length(stat_rows) + 1]] <- data.frame(
      panel = panel_id, parent = parent, state = spec$state, marker = mk,
      n_ctrl = length(ctrl_v), n_trt = length(trt_v),
      mean_ctrl = mean(ctrl_v, na.rm = TRUE), mean_trt = mean(trt_v, na.rm = TRUE),
      sd_ctrl = if (length(ctrl_v) >= 2) stats::sd(ctrl_v) else NA_real_,
      sd_trt = if (length(trt_v) >= 2) stats::sd(trt_v) else NA_real_,
      p_value = pv,
      dropped_ctrl = if (!is.null(dropped_bio) && nrow(dropped_bio)) {
        paste(dropped_bio$dropped_bio[dropped_bio$group == flow_ctrl_group], collapse = ",")
      } else "",
      dropped_trt = if (!is.null(dropped_bio) && nrow(dropped_bio)) {
        paste(dropped_bio$dropped_bio[dropped_bio$group == flow_trt_group], collapse = ",")
      } else "",
      stringsAsFactors = FALSE
    )
    stub_one <- paste0(
      panel_id, "_", parent, "_", gsub("[^A-Za-z0-9]+", "_", mk), "_",
      flow_comparison_tag()
    )
    utils::write.csv(samp, file.path(out_dir, paste0(stub_one, "_by_sample.csv")), row.names = FALSE)
    utils::write.csv(samp_bio, file.path(out_dir, paste0(stub_one, "_by_bio.csv")), row.names = FALSE)
  }
  if (!length(columns)) return(invisible(NULL))
  stub <- file.path(
    out_dir,
    paste0(panel_id, "_", parent, "_", spec$state, "_", flow_comparison_tag())
  )
  save_func_state_figure(columns, stub)
  if (length(stat_rows)) {
    st <- do.call(rbind, stat_rows)
    utils::write.csv(st, paste0(stub, "_stats.csv"), row.names = FALSE)
  }
  log_msg(panel_id, " functional-state ", parent, " ", spec$state, " -> ", stub)
  invisible(stub)
}

export_functional_state_figures <- function(cells, panel_id, out_dir) {
  specs <- functional_state_specs(panel_id)
  if (!length(specs)) return(invisible(NULL))
  func_dir <- file.path(out_dir, "functional_state")
  dir.create(func_dir, recursive = TRUE, showWarnings = FALSE)
  if (identical(panel_id, "P2")) {
    log_msg("P2 B subsets: no PD-L1/LAG-3/TIM-3 on this panel; skip exhaustion, report CD86/CD80/CD40")
  }
  for (spec in specs) {
    tryCatch(
      export_one_functional_state(cells, panel_id, spec, func_dir),
      error = function(e) log_msg(panel_id, " skip ", spec$parent, " ", spec$state, ": ", e$message)
    )
  }
  log_msg(panel_id, " functional-state figures: ", func_dir)
  invisible(TRUE)
}

export_dimred_plots <- function(cells, med, annot, freq_df, panel_id, out_dir, umap_is_pca = FALSE, tsne_is_pca = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  marker_dir <- file.path(out_dir, "markers")
  dir.create(marker_dir, showWarnings = FALSE)
  umap_lab <- if (umap_is_pca) "PCA (UMAP fallback)" else "UMAP"
  tsne_lab <- if (tsne_is_pca) "PCA (tSNE fallback)" else "tSNE"
  tag <- paste0(panel_id, " JZ-AB / JZ-EVB")
  split_ttl <- paste0(panel_id, "  JZ-EVB | JZ-AB")
  log_msg(
    panel_id,
    " Figure 1 (*_major_split): all cells by major class; ",
    "*_lineage_split: same embedding, every fine subset; ",
    "dimred_by_major/: re-embed each class with high-contrast subset colors"
  )
  cells$dimred_major <- dimred_major_of(
    panel_id,
    cells$lineage,
    if ("cluster_lineage" %in% names(cells)) cells$cluster_lineage else NULL
  )
  n_major <- length(unique(major_display_label(cells$dimred_major)))
  p_major_tsne <- plot_split_lineage(
    cells, "tSNE1", "tSNE2", panel_id, "tSNE-1", "tSNE-2",
    paste0(split_ttl, "  major classes"), color_mode = "major"
  )
  p_major_umap <- plot_split_lineage(
    cells, "UMAP1", "UMAP2", panel_id, "UMAP-1", "UMAP-2",
    paste0(split_ttl, "  major classes"), color_mode = "major"
  )
  save_split_dr(p_major_tsne, file.path(out_dir, paste0(panel_id, "_JZ_AB_vs_JZ_EVB_tSNE_major_split")), n_major)
  save_split_dr(p_major_umap, file.path(out_dir, paste0(panel_id, "_JZ_AB_vs_JZ_EVB_UMAP_major_split")), n_major)

  p_subset_tsne <- plot_split_lineage(
    cells, "tSNE1", "tSNE2", panel_id, "tSNE-1", "tSNE-2",
    paste0(split_ttl, "  all subsets"), color_mode = "subset"
  )
  p_subset_umap <- plot_split_lineage(
    cells, "UMAP1", "UMAP2", panel_id, "UMAP-1", "UMAP-2",
    paste0(split_ttl, "  all subsets"), color_mode = "subset"
  )
  n_subset <- length(unique(as.character(p_subset_tsne$data$celltype)))
  save_split_dr(p_subset_tsne, file.path(out_dir, paste0(panel_id, "_JZ_AB_vs_JZ_EVB_tSNE_lineage_split")), n_subset)
  save_split_dr(p_subset_umap, file.path(out_dir, paste0(panel_id, "_JZ_AB_vs_JZ_EVB_UMAP_lineage_split")), n_subset)
  save_split_dr(p_subset_umap, file.path(out_dir, paste0(panel_id, "_JZ_AB_vs_JZ_EVB_UMAP_lineage_split_joint")), n_subset)

  tryCatch(
    export_major_subset_dimred(cells, panel_id, out_dir),
    error = function(e) log_msg(panel_id, " per-major subset dimred failed: ", e$message)
  )

  save_gg(plot_embedding(cells, "UMAP1", "UMAP2", "group", paste(tag, "-", umap_lab, "by group")),
          file.path(out_dir, paste0(panel_id, "_UMAP_by_group")))
  save_gg(plot_embedding(cells, "UMAP1", "UMAP2", "sample", paste(tag, "-", umap_lab, "by sample")),
          file.path(out_dir, paste0(panel_id, "_UMAP_by_sample")), width = 8)
  save_gg(plot_embedding(cells, "UMAP1", "UMAP2", "cluster", paste(tag, "-", umap_lab, "by cluster")),
          file.path(out_dir, paste0(panel_id, "_UMAP_by_cluster")), width = 8)
  save_gg(plot_embedding(cells, "UMAP1", "UMAP2", "lineage", paste(tag, "-", umap_lab, "by lineage")),
          file.path(out_dir, paste0(panel_id, "_UMAP_by_lineage")), width = 8)
  save_gg(plot_density_split(cells, "UMAP1", "UMAP2", paste(tag, "-", umap_lab, "density JZ-AB / JZ-EVB")),
          file.path(out_dir, paste0(panel_id, "_UMAP_density_JZ_AB_vs_JZ_EVB")), width = 10, height = 5)

  save_gg(plot_embedding(cells, "tSNE1", "tSNE2", "group", paste(tag, "-", tsne_lab, "by group")),
          file.path(out_dir, paste0(panel_id, "_tSNE_by_group")))
  save_gg(plot_embedding(cells, "tSNE1", "tSNE2", "sample", paste(tag, "-", tsne_lab, "by sample")),
          file.path(out_dir, paste0(panel_id, "_tSNE_by_sample")), width = 8)
  save_gg(plot_embedding(cells, "tSNE1", "tSNE2", "cluster", paste(tag, "-", tsne_lab, "by cluster")),
          file.path(out_dir, paste0(panel_id, "_tSNE_by_cluster")), width = 8)
  save_gg(plot_embedding(cells, "tSNE1", "tSNE2", "lineage", paste(tag, "-", tsne_lab, "by lineage")),
          file.path(out_dir, paste0(panel_id, "_tSNE_by_lineage")), width = 8)
  save_gg(plot_density_split(cells, "tSNE1", "tSNE2", paste(tag, "-", tsne_lab, "density JZ-AB / JZ-EVB")),
          file.path(out_dir, paste0(panel_id, "_tSNE_density_JZ_AB_vs_JZ_EVB")), width = 10, height = 5)

  dr_cols <- colnames(med)
  for (mk in dr_cols) {
    if (!mk %in% names(cells)) next
    save_gg(plot_marker_embedding(cells, "UMAP1", "UMAP2", mk, paste(tag, "UMAP", mk)),
            file.path(marker_dir, paste0(panel_id, "_UMAP_", gsub("[^A-Za-z0-9]+", "_", mk))),
            width = 6.5, height = 5.5)
  }

  plot_cluster_heatmap(med, paste(tag, "cluster median markers (z)"),
                       file.path(out_dir, paste0(panel_id, "_cluster_marker_heatmap")))
  tryCatch(
    export_annotation_heatmaps(cells, panel_id, out_dir),
    error = function(e) log_msg(panel_id, " annotation heatmap failed: ", e$message)
  )

  freq_df$cluster <- factor(freq_df$cluster, levels = annot$cluster)
  freq_bio <- maybe_trim_bio(aggregate_freq_by_bio(freq_df, "cluster"), "cluster")
  freq_ylab <- if (flow_should_trim_bio()) {
    paste(tag, "cluster frequency (bio-rep; max/min dropped)")
  } else {
    paste(tag, "cluster frequency (bio-rep)")
  }
  save_gg(plot_freq_box(freq_bio, freq_ylab),
          file.path(out_dir, paste0(panel_id, "_cluster_frequency_JZ_AB_vs_JZ_EVB")),
          width = max(8, 0.7 * length(unique(freq_df$cluster)) + 2), height = 5.5)

  lin_freq <- lineage_frequencies(cells)
  lin_bio <- maybe_trim_bio(aggregate_freq_by_bio(lin_freq, "lineage"), "lineage")
  names(lin_bio)[names(lin_bio) == "lineage"] <- "cluster"
  lin_ylab <- if (flow_should_trim_bio()) {
    paste(tag, "lineage frequency (bio-rep; max/min dropped)")
  } else {
    paste(tag, "lineage frequency (bio-rep)")
  }
  save_gg(plot_freq_box(lin_bio, lin_ylab),
          file.path(out_dir, paste0(panel_id, "_lineage_frequency_JZ_AB_vs_JZ_EVB")),
          width = max(7, 0.8 * length(unique(lin_bio$cluster)) + 2), height = 5.5)

  mean_df <- aggregate(percent ~ group + cluster, data = freq_bio, FUN = mean)
  names(mean_df)[names(mean_df) == "percent"] <- "mean_percent"
  save_gg(plot_freq_bar(mean_df, paste(tag, "mean cluster composition")),
          file.path(out_dir, paste0(panel_id, "_cluster_composition_stacked")),
          width = 6, height = 5.5)

  if (has_pkg("cowplot")) {
    p1 <- plot_embedding(cells, "UMAP1", "UMAP2", "group", "Group")
    p2 <- plot_embedding(cells, "UMAP1", "UMAP2", "lineage", "Lineage")
    p3 <- plot_freq_box(lin_bio, "Lineage %")
    overview <- cowplot::plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(1, 1.1, 1.2))
    title <- cowplot::ggdraw() +
      cowplot::draw_label(paste(tag, "dimensionality reduction overview"), fontface = "bold", size = 14)
    overview <- cowplot::plot_grid(title, overview, ncol = 1, rel_heights = c(0.08, 1))
    save_gg(overview, file.path(out_dir, paste0(panel_id, "_JZ_AB_vs_JZ_EVB_dimred_overview")),
            width = 16, height = 5.8)
  }
  tryCatch(
    export_subset_gate_figures(cells, panel_id, out_dir),
    error = function(e) log_msg(panel_id, " subset stat+contour figures failed: ", e$message)
  )
  tryCatch(
    export_per_sample_gating_figures(cells, panel_id, out_dir),
    error = function(e) log_msg(panel_id, " per-sample gating figures failed: ", e$message)
  )
  if (identical(panel_id, "P1")) {
    tryCatch(
      export_p1_activation_stats(cells, out_dir),
      error = function(e) log_msg(panel_id, " T/NK activation MFI failed: ", e$message)
    )
  }
  if (identical(panel_id, "P2")) {
    tryCatch(
      export_p2_activation_stats(cells, out_dir),
      error = function(e) log_msg(panel_id, " B-cell activation MFI failed: ", e$message)
    )
  }
  if (identical(panel_id, "P3")) {
    tryCatch(
      export_p3_activation_stats(cells, out_dir),
      error = function(e) log_msg(panel_id, " APC/TNF-a MFI failed: ", e$message)
    )
  }
  tryCatch(
    export_functional_state_figures(cells, panel_id, out_dir),
    error = function(e) log_msg(panel_id, " functional-state figures failed: ", e$message)
  )
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# 8. 频率统计 JZ-AB / JZ-EVB
# -----------------------------------------------------------------------------
# 旧表里的 "T" 单独成列时 read.csv 会变成 TRUE
read_embed_csv <- function(path) {
  df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  chr <- intersect(c("sample", "group", "cluster", "cluster_lineage", "lineage",
                     "bio_sample", "tech_rep"), names(df))
  for (cn in chr) {
    v <- as.character(df[[cn]])
    v[v %in% c("TRUE", "True")] <- "T"
    v[v %in% c("FALSE", "False")] <- "F"
    df[[cn]] <- v
  }
  df
}

cluster_frequencies <- function(cells) {
  tab <- as.data.frame(table(sample = cells$sample, cluster = cells$cluster), stringsAsFactors = FALSE)
  names(tab)[3] <- "n"
  meta_cols <- intersect(c("sample", "group", "bio_sample", "tech_rep"), names(cells))
  smp_group <- unique(cells[, meta_cols, drop = FALSE])
  tab <- merge(tab, smp_group, by = "sample")
  tot <- aggregate(n ~ sample, data = tab, FUN = sum)
  names(tot)[2] <- "total"
  tab <- merge(tab, tot, by = "sample")
  tab$percent <- 100 * tab$n / pmax(tab$total, 1)
  lin_col <- if ("cluster_lineage" %in% names(cells)) "cluster_lineage" else "lineage"
  if (lin_col %in% names(cells)) {
    maj <- aggregate(cells[[lin_col]], by = list(cluster = cells$cluster), FUN = function(x) {
      names(sort(table(x), decreasing = TRUE))[1]
    })
    names(maj)[2] <- "lineage"
    tab <- merge(tab, maj, by = "cluster", all.x = TRUE)
  }
  tab
}

lineage_frequencies <- function(cells) {
  tab <- as.data.frame(table(sample = cells$sample, lineage = cells$lineage), stringsAsFactors = FALSE)
  names(tab)[3] <- "n"
  meta_cols <- intersect(c("sample", "group", "bio_sample", "tech_rep"), names(cells))
  smp_group <- unique(cells[, meta_cols, drop = FALSE])
  tab <- merge(tab, smp_group, by = "sample")
  tot <- aggregate(n ~ sample, data = tab, FUN = sum)
  names(tot)[2] <- "total"
  tab <- merge(tab, tot, by = "sample")
  tab$percent <- 100 * tab$n / pmax(tab$total, 1)
  tab
}

compare_group_freq <- function(freq_df, id_col = "cluster", trim_bio = NULL) {
  freq_df <- aggregate_freq_by_bio(freq_df, id_col)
  freq_df <- maybe_trim_bio(freq_df, id_col, trim_bio)
  dropped <- attr(freq_df, "dropped")
  ids <- unique(as.character(freq_df[[id_col]]))
  g_ctrl <- flow_ctrl_group
  g_trt <- flow_trt_group
  rows <- lapply(ids, function(id) {
    sub <- freq_df[as.character(freq_df[[id_col]]) == id, ]
    ctrl_vals <- sub$percent[as.character(sub$group) == g_ctrl]
    trt_vals <- sub$percent[as.character(sub$group) == g_trt]
    if (length(ctrl_vals) < 2 || length(trt_vals) < 2) {
      pv <- NA_real_
    } else {
      pv <- tryCatch(t.test(trt_vals, ctrl_vals)$p.value, error = function(e) NA_real_)
    }
    drop_ctrl <- NA_character_
    drop_trt <- NA_character_
    if (!is.null(dropped) && nrow(dropped)) {
      dc <- dropped$dropped_bio[dropped$id == id & dropped$group == g_ctrl]
      dt <- dropped$dropped_bio[dropped$id == id & dropped$group == g_trt]
      if (length(dc)) drop_ctrl <- paste(dc, collapse = ",")
      if (length(dt)) drop_trt <- paste(dt, collapse = ",")
    }
    row <- data.frame(
      id = id,
      n_ctrl = length(ctrl_vals),
      n_trt = length(trt_vals),
      mean_ctrl = mean(ctrl_vals),
      mean_trt = mean(trt_vals),
      sd_ctrl = sd(ctrl_vals),
      sd_trt = sd(trt_vals),
      delta_trt_minus_ctrl = mean(trt_vals) - mean(ctrl_vals),
      p_value = pv,
      dropped_ctrl = drop_ctrl,
      dropped_trt = drop_trt,
      stringsAsFactors = FALSE
    )
    names(row) <- c(
      id_col,
      paste0("n_", g_ctrl), paste0("n_", g_trt),
      paste0("mean_", g_ctrl), paste0("mean_", g_trt),
      paste0("sd_", g_ctrl), paste0("sd_", g_trt),
      paste0("delta_", g_trt, "_minus_", g_ctrl),
      "p_value",
      paste0("dropped_", g_ctrl), paste0("dropped_", g_trt)
    )
    row
  })
  out <- do.call(rbind, rows)
  out$padj <- if (all(is.na(out$p_value))) NA_real_ else p.adjust(out$p_value, method = "BH")
  attr(out, "dropped") <- dropped
  out[order(out$p_value, na.last = TRUE), ]
}

# -----------------------------------------------------------------------------
# 9. 单 panel 流程
# -----------------------------------------------------------------------------
downsample_idx <- function(n, cap) {
  if (n <= cap) return(seq_len(n))
  sample.int(n, cap)
}

load_panel_cells <- function(panel_id, file_tab, use_demo, n_cap) {
  markers <- dr_marker_names(panel_id)
  all_expr <- list()
  meta_group <- character()
  meta_sample <- character()
  meta_bio <- character()
  meta_tech <- character()
  meta_true <- character()
  map_used <- NULL

  if (use_demo) {
    log_msg(panel_id, " : DEMO synthetic cells (not real FCS)")
    samples <- expand.grid(
      group = flow_group_levels, replicate = 1:3, tech = 1:2,
      stringsAsFactors = FALSE
    )
    samples$bio_sample <- paste0(samples$group, "-", samples$replicate)
    samples$sample <- paste0(samples$group, samples$replicate, "-", samples$tech)
    for (i in seq_len(nrow(samples))) {
      set.seed(seed_value + i + as.integer(factor(panel_id)) * 10)
      d <- make_demo_sample(panel_id, samples$group[i], samples$sample[i], min(n_cap, 2500))
      all_expr[[i]] <- d$expr[, c("CD45", "L/D", markers), drop = FALSE]
      meta_group <- c(meta_group, rep(samples$group[i], nrow(d$expr)))
      meta_sample <- c(meta_sample, rep(samples$sample[i], nrow(d$expr)))
      meta_bio <- c(meta_bio, rep(samples$bio_sample[i], nrow(d$expr)))
      meta_tech <- c(meta_tech, rep(as.character(samples$tech[i]), nrow(d$expr)))
      meta_true <- c(meta_true, d$true_lineage)
    }
  } else {
    sub <- file_tab[file_tab$panel == panel_id, ]
    if (nrow(sub) == 0) {
      log_msg(panel_id, " : no unmixed files (need JZ-EVB1-1_", panel_id, "_unmixed.fcs / EVB1-1_", panel_id, "_unmixed.fcs)")
      return(NULL)
    }
    for (i in seq_len(nrow(sub))) {
      log_msg("Read ", sub$file[i], "  sample=", sub$sample[i], " bio=", sub$bio_sample[i])
      rec <- read_fcs_expr(sub$path[i], panel_id)
      keep <- qc_filter_matrix(rec$exprs, rec$names, rec$map, panel_id)
      exprs <- rec$exprs[keep, , drop = FALSE]
      if (nrow(exprs) < 50) {
        log_msg("Too few events after QC: ", sub$file[i])
        next
      }
      set.seed(seed_value + i)
      take <- downsample_idx(nrow(exprs), n_cap)
      exprs <- exprs[take, , drop = FALSE]
      mat <- extract_marker_mat(exprs, rec$map, c("CD45", "L/D", markers))
      miss <- setdiff(markers, colnames(mat))
      if (length(miss) > 0) log_msg("Missing markers in ", sub$file[i], ": ", paste(miss, collapse = ", "))
      for (m in setdiff(c("CD45", "L/D", markers), colnames(mat))) {
        mat <- cbind(mat, NA_real_)
        colnames(mat)[ncol(mat)] <- m
      }
      mat <- mat[, c("CD45", "L/D", markers), drop = FALSE]
      all_expr[[length(all_expr) + 1]] <- mat
      bio <- if ("bio_sample" %in% names(sub)) as.character(sub$bio_sample[i]) else as.character(sub$sample[i])
      tech <- if ("tech_rep" %in% names(sub)) as.character(sub$tech_rep[i]) else NA_character_
      meta_group <- c(meta_group, rep(sub$group[i], nrow(mat)))
      meta_sample <- c(meta_sample, rep(sub$sample[i], nrow(mat)))
      meta_bio <- c(meta_bio, rep(bio, nrow(mat)))
      meta_tech <- c(meta_tech, rep(tech, nrow(mat)))
      if (is.null(map_used)) map_used <- rec$map
    }
  }

  if (length(all_expr) == 0) return(NULL)
  expr <- do.call(rbind, all_expr)
  # 演示数据已在 asinh 空间；真实 FCS 再 asinh(x/150)
  tr <- if (use_demo) as.matrix(expr) else asinh_mat(expr)
  list(
    transformed = tr,
    group = meta_group,
    sample = meta_sample,
    bio_sample = meta_bio,
    tech_rep = meta_tech,
    true_lineage = if (length(meta_true)) meta_true else rep(NA_character_, nrow(tr)),
    markers = markers,
    map = map_used
  )
}

analyze_one_panel <- function(panel_id, file_tab, use_demo) {
  log_msg("==== Panel ", panel_id, " : ", panel_map$panels[[panel_id]]$focus, " ====")
  if (identical(panel_id, "P1")) {
    log_msg("P1 gates: naive / T_CM / T_SCM / T_EM early-late / SLEC / MPEC / T_EFF / exhausted; CD69 activation")
  }
  if (identical(panel_id, "P2")) {
    log_msg("P2 gates: wide mononuclear FSC/SSC; CD19+ IgD vs CD27 Naive/Unswitched/Switched; MZ IgM-high; plasmablast/plasma BLIMP; CD40/CD80/CD86 as MFI not subsets")
  }
  if (identical(panel_id, "P3")) {
    log_msg("P3 gates: singlets/live/CD45+; T=CD3+ B=CD19+ NK=NK1.1+; myeloid=triple-neg then CD11B+/-; eos=Siglec-F+CCR3+; mast=FceRI+CD200R3+; F4/80 hi mac vs Ly6C hi mono; DC=CD11C+MHCII+ then cDC1/cDC2; CD80/CD86/CD40/TNF-a as MFI")
  }
  out_dir <- file.path(result_dir, panel_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dat <- load_panel_cells(panel_id, file_tab, use_demo, n_per_sample_dr)
  if (is.null(dat)) {
    writeLines("no cells", file.path(out_dir, paste0(panel_id, "_EMPTY.txt")))
    return(invisible(NULL))
  }

  feat <- dat$markers
  feat <- feat[feat %in% colnames(dat$transformed)]
  exclude_dr <- c("L/D", "CD45")
  if (!is.null(panel_map$qc$exclude_from_dr)) {
    exclude_dr <- unique(c(exclude_dr, unlist(panel_map$qc$exclude_from_dr)))
  }
  feat <- setdiff(feat, exclude_dr)
  mat_raw <- dat$transformed[, feat, drop = FALSE]
  na_frac <- colMeans(!is.finite(as.matrix(mat_raw)))
  drop_mk <- names(na_frac)[is.finite(na_frac) & na_frac > 0.2]
  if (length(drop_mk)) {
    log_msg(panel_id, " dropping unmatched/empty markers from DR (do not drop all cells): ",
            paste(drop_mk, collapse = ", "))
    feat <- setdiff(feat, drop_mk)
    mat_raw <- dat$transformed[, feat, drop = FALSE]
  }
  ok_row <- rowSums(!is.finite(as.matrix(mat_raw))) == 0
  mat_raw <- mat_raw[ok_row, , drop = FALSE]
  grp <- dat$group[ok_row]
  smp <- dat$sample[ok_row]
  bio <- if (!is.null(dat$bio_sample) && length(dat$bio_sample) >= length(ok_row)) {
    dat$bio_sample[ok_row]
  } else {
    smp
  }
  tech <- if (!is.null(dat$tech_rep) && length(dat$tech_rep) >= length(ok_row)) {
    dat$tech_rep[ok_row]
  } else {
    rep(NA_character_, length(smp))
  }
  true_lin <- dat$true_lineage[ok_row]
  mat <- scale_markers(mat_raw)
  mat <- weight_dr_markers(mat, panel_id)
  w_used <- dr_lineage_marker_weights(panel_id, colnames(mat))
  w_hi <- names(w_used)[w_used > 1]
  if (length(w_hi)) {
    log_msg(panel_id, " joint DR upweights lineage markers so CD4/CD8/NK do not collapse: ",
            paste(sprintf("%s x%s", w_hi, w_used[w_hi]), collapse = ", "))
  }

  log_msg(panel_id, " cells=", nrow(mat), " markers=", paste(feat, collapse = ","))
  log_msg(panel_id, " frequencies are from equal-n subsample used for embedding")
  pca <- run_pca(mat, npcs = min(15, ncol(mat)))
  npcs_use <- min(ncol(pca$embedding), max(5, ncol(mat)))
  pca_use <- pca$embedding[, seq_len(npcs_use), drop = FALSE]

  umap_is_pca <- !has_pkg("uwot")
  tsne_is_pca <- !has_pkg("Rtsne")
  umap <- run_umap(pca_use)
  tsne <- run_tsne(pca_use)
  cl <- cluster_cells(mat, panel_id)

  cells <- data.frame(
    sample = smp,
    bio_sample = bio,
    tech_rep = tech,
    group = factor(grp, levels = flow_group_levels),
    cluster = cl,
    UMAP1 = umap[, 1],
    UMAP2 = umap[, 2],
    tSNE1 = tsne[, 1],
    tSNE2 = tsne[, 2],
    as.data.frame(mat_raw, check.names = FALSE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  med_df <- aggregate(cells[, feat, drop = FALSE], by = list(cluster = cells$cluster), FUN = median)
  med <- as.matrix(med_df[, -1, drop = FALSE])
  rownames(med) <- as.character(med_df$cluster)

  annot <- annotate_clusters(med, panel_id)
  log_msg(panel_id, " cluster labels: ", paste(paste(annot$cluster, annot$lineage, sep = "="), collapse = ", "))
  hier <- hierarchical_gate_by_sample(as.matrix(mat_raw), smp, panel_id)
  cells$cluster_lineage <- hier$major
  cells$lineage <- hier$subset
  maj_n <- sort(table(hier$major), decreasing = TRUE)
  lin_n <- sort(table(hier$subset), decreasing = TRUE)
  log_msg(panel_id, " layer1 major: ", paste(paste(names(maj_n), as.integer(maj_n), sep = "="), collapse = ", "))
  log_msg(panel_id, " layer2 subset: ", paste(paste(names(lin_n), as.integer(lin_n), sep = "="), collapse = ", "))
  for (s in unique(as.character(smp))) {
    ss <- cells$lineage[as.character(cells$sample) == s]
    tab <- sort(table(ss), decreasing = TRUE)
    log_msg(panel_id, " ", s, " subsets: ", paste(paste(names(tab), as.integer(tab), sep = "="), collapse = ", "))
  }
  if (use_demo) cells$true_lineage <- true_lin

  freq_df <- cluster_frequencies(cells)
  stats_cl <- compare_group_freq(freq_df, "cluster")
  lin_freq <- lineage_frequencies(cells)
  stats_lin <- compare_group_freq(lin_freq, "lineage")
  if (flow_should_trim_bio()) {
    log_msg(panel_id, " JZ-AB / JZ-EVB stats: tech reps averaged, then drop 1 extreme bio-rep (max or min) per group; n=2")
  } else {
    log_msg(panel_id, " JZ-AB / JZ-EVB stats: tech reps averaged to n=3 bio (no extreme dropped)")
  }

  utils::write.csv(annot, file.path(out_dir, paste0(panel_id, "_cluster_annotation.csv")), row.names = FALSE)
  utils::write.csv(cbind(cluster = rownames(med), as.data.frame(med)),
                   file.path(out_dir, paste0(panel_id, "_cluster_median_markers.csv")), row.names = FALSE)
  utils::write.csv(freq_df, file.path(out_dir, paste0(panel_id, "_cluster_frequency_by_sample.csv")), row.names = FALSE)
  utils::write.csv(aggregate_freq_by_bio(freq_df, "cluster"),
                   file.path(out_dir, paste0(panel_id, "_cluster_frequency_by_bio.csv")), row.names = FALSE)
  if (flow_should_trim_bio()) {
    cl_trim <- trim_bio_extremes(aggregate_freq_by_bio(freq_df, "cluster"), "cluster")
    utils::write.csv(cl_trim,
                     file.path(out_dir, paste0(panel_id, "_cluster_frequency_by_bio_trimmed.csv")), row.names = FALSE)
    dropped_cl <- attr(cl_trim, "dropped")
    dropped_lin <- attr(stats_lin, "dropped")
    drop_rows <- Filter(Negate(is.null), list(
      if (!is.null(dropped_cl)) cbind(panel = panel_id, table = "cluster", dropped_cl),
      if (!is.null(dropped_lin)) cbind(panel = panel_id, table = "lineage", dropped_lin)
    ))
    if (length(drop_rows)) {
      utils::write.csv(do.call(rbind, drop_rows),
                       file.path(out_dir, paste0(panel_id, "_dropped_bio_extremes.csv")), row.names = FALSE)
    }
  }
  utils::write.csv(stats_cl, file.path(out_dir, paste0(panel_id, "_cluster_JZ_AB_vs_JZ_EVB_stats.csv")), row.names = FALSE)
  utils::write.csv(stats_lin, file.path(out_dir, paste0(panel_id, "_lineage_JZ_AB_vs_JZ_EVB_stats.csv")), row.names = FALSE)
  if (!is.null(dat$map)) {
    utils::write.csv(dat$map, file.path(out_dir, paste0(panel_id, "_channel_map.csv")), row.names = FALSE)
  }
  embed_cols <- c("sample", "bio_sample", "tech_rep", "group", "cluster", "cluster_lineage", "lineage",
                  "UMAP1", "UMAP2", "tSNE1", "tSNE2")
  extra_cols <- setdiff(names(cells), c(embed_cols, "true_lineage"))
  embed_out <- cells[, intersect(c(embed_cols, extra_cols), names(cells)), drop = FALSE]
  utils::write.csv(embed_out, file.path(out_dir, paste0(panel_id, "_cell_embeddings.csv")), row.names = FALSE)

  export_dimred_plots(cells, med, annot, freq_df, panel_id, out_dir, umap_is_pca, tsne_is_pca)
  log_msg(panel_id, " plots written to ", out_dir)
  list(annot = annot, stats_cluster = stats_cl, stats_lineage = stats_lin, n = nrow(cells))
}

# -----------------------------------------------------------------------------
# 10. 本文件只提供函数。分析入口是 JZ_Flow_dimred_pipeline.R
# -----------------------------------------------------------------------------
jz_engine_loaded <- TRUE
log_msg("JZ_flow_engine.R: function library only; run JZ_Flow_dimred_pipeline.R")
