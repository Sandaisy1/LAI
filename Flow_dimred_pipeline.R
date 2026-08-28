#!/usr/bin/env Rscript
# =============================================================================
# 流式降维：H vs EV（P1 T/NK，P2 B，P3 髓系）
# 输入：E:/R/flow J-LJY WJZ ZZX 下的 *_unmixed.fcs（不要用 raw）
# 每个 panel 单独联合 UMAP/tSNE，导出 PDF+PNG，并比较 H vs EV 细胞频率
#
# 用法：
#   setwd("E:/R/flow J-LJY WJZ ZZX")
#   source("Flow_dimred_pipeline.R")
# 三个 panel 跑完后会自动汇总全亚群频率，并按大类画轨迹：
#   source("Flow_dimred_all_subsets.R")   # 也可单独重出 results_flow/all_subsets/
#   source("Flow_dimred_trajectory.R")    # 也可单独重出 P1/P2/P3/trajectory/
# 无 FCS 时可跑演示数据（会在日志里标明 DEMO，不可当正式结果）：
#   Sys.setenv(FLOW_DEMO = "1")
#   source("Flow_dimred_pipeline.R")
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

# 读写目录：新数据夹优先；旧 E:/R/flow J 仅作回退
flow_primary_data_dir <- "E:/R/flow J-LJY WJZ ZZX"
flow_legacy_data_dir <- "E:/R/flow J"

resolve_flow_dir <- function() {
  env_dir <- Sys.getenv("FLOW_DIR", unset = "")
  preferred <- c(
    env_dir,
    flow_primary_data_dir,
    "E:\\R\\flow J-LJY WJZ ZZX"
  )
  preferred <- unique(preferred[nzchar(preferred)])
  for (d in preferred) {
    if (dir.exists(d)) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  candidates <- c(
    flow_legacy_data_dir,
    "E:\\R\\flow J",
    file.path(script_dir, "flow_J"),
    file.path(script_dir, "flow J"),
    file.path(script_dir, "flow J-LJY WJZ ZZX"),
    script_dir,
    getwd()
  )
  candidates <- unique(candidates[nzchar(candidates)])
  for (d in candidates) {
    if (!dir.exists(d)) next
    hits <- list.files(d, pattern = "(?i)_unmixed\\.fcs$", recursive = FALSE)
    if (length(hits) > 0) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  if (dir.exists(flow_legacy_data_dir)) {
    return(normalizePath(flow_legacy_data_dir, winslash = "/", mustWork = FALSE))
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
# 2. 文件名：EV1–3 与 H1–3（先匹配 EV，避免被单字母 H 误伤）
# -----------------------------------------------------------------------------
flow_ctrl_group <- "EV"
flow_trt_group <- "H"
flow_group_levels <- c("EV", "H")

parse_fcs_filename <- function(path) {
  b <- basename(path)
  # EV1_P3 / EV1-P3 / EV-1_P3 / EV1_Panel3 / EV1_P03
  m <- regexec(
    "^(EV|H)[-_ ]?([123])[-_ ]+(?:PANEL[-_ ]?)?P?0?([123])[-_ ].*(unmixed|raw)\\.fcs$",
    b, ignore.case = TRUE
  )
  hit <- regmatches(b, m)[[1]]
  if (length(hit) < 5) return(NULL)
  grp <- toupper(hit[2])
  list(
    file = b,
    path = path,
    group = grp,
    replicate = hit[3],
    sample = paste0(grp, hit[3]),
    panel = paste0("P", hit[4]),
    kind = tolower(hit[5])
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
      file = x$file, path = x$path, group = x$group, replicate = x$replicate,
      sample = x$sample, panel = x$panel, kind = x$kind, stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  df <- df[!duplicated(df$path), , drop = FALSE]
  df <- df[order(df$panel, df$group, df$replicate), ]
  rownames(df) <- NULL
  df
}

# -----------------------------------------------------------------------------
# 3. Panel 地图
# Windows 常把“另存为”弄成 flow_panel_map.json.txt（资源管理器显示“文本文档”）
# -----------------------------------------------------------------------------
panel_map_search_dirs <- function() {
  dirs <- c(
    getwd(),
    script_dir,
    project_dir,
    flow_primary_data_dir,
    "E:\\R\\flow J-LJY WJZ ZZX",
    flow_legacy_data_dir,
    "E:\\R\\flow J",
    Sys.getenv("FLOW_DIR", unset = "")
  )
  dirs <- unique(dirs[nzchar(dirs)])
  dirs[dir.exists(dirs)]
}

find_panel_map_file <- function(dirs = panel_map_search_dirs()) {
  dirs <- dirs[nzchar(dirs) & dir.exists(dirs)]
  names <- c(
    "flow_panel_map.json",
    "flow_panel_map.json.txt",
    "flow_panel_map.txt"
  )
  exact <- unlist(lapply(dirs, function(d) file.path(d, names)), use.names = FALSE)
  hit <- exact[file.exists(exact)]
  if (length(hit) > 0) return(hit[[1]])
  fuzzy <- unlist(lapply(dirs, function(d) {
    list.files(d, pattern = "^flow_panel_map", ignore.case = TRUE, full.names = TRUE)
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
    "找不到 flow_panel_map.json。\n",
    "Windows 若类型是“文本文档”、图标是记事本，真实文件名多半是 flow_panel_map.json.txt。\n",
    "在资源管理器：查看 → 去掉“隐藏已知文件类型的扩展名”，再改名为 flow_panel_map.json。\n",
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
      log_msg("NK1.1 is AF700 on the sheet; Cytek unmixed often names that detector R718-A / APC-R700-A. Both should match.")
    }
  }
  list(exprs = exprs, names = nms, desc = desc, map = map)
}

qc_filter_matrix <- function(exprs, names, map) {
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
  ld_idx <- map$channel_index[map$marker == "L/D"]
  if (length(ld_idx) && !is.na(ld_idx[1])) {
    ld <- exprs[, ld_idx[1]]
    keep <- keep & ld <= quantile(ld, 0.95, na.rm = TRUE)
  }
  cd45_idx <- map$channel_index[map$marker == "CD45"]
  if (length(cd45_idx) && !is.na(cd45_idx[1])) {
    cd <- exprs[, cd45_idx[1]]
    keep <- keep & cd >= quantile(cd, 0.2, na.rm = TRUE)
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
    CD4_TCM = pop(CD3 = 3.2, CD4 = 3.0, CD62L = 2.8, CD44 = 2.8, CD27 = 2.4),
    CD4_TEM = pop(CD3 = 3.2, CD4 = 3.0, CD62L = 0.3, CD44 = 3.1, CD27 = 0.8),
    Treg = pop(CD3 = 3.1, CD4 = 3.0, CD25 = 3.2, CD69 = 0.4, CD44 = 1.8),
    CD4_act = pop(CD3 = 3.2, CD4 = 3.0, CD25 = 1.2, CD69 = 3.1, CD44 = 2.6, CD62L = 0.4, `TNF-a` = 2.4, `IFN-g` = 1.8),
    CD8_naive = pop(CD3 = 3.3, CD8 = 3.1, CD8b = 2.9, CD4 = 0.1, CD62L = 3.1, CD44 = 0.3),
    CD8_TCM = pop(CD3 = 3.3, CD8 = 3.1, CD8b = 2.8, CD62L = 2.7, CD44 = 2.6),
    CD8_TEM = pop(CD3 = 3.3, CD8 = 3.1, CD8b = 2.8, CD44 = 3.0, CD62L = 0.3),
    CD8_eff = pop(CD3 = 3.3, CD8 = 3.1, CD8b = 2.8, CD44 = 2.6, GZMB = 3.0, Perforin = 2.8, `IFN-g` = 2.2),
    CD8_exh = pop(CD3 = 3.0, CD8 = 3.0, CD8b = 2.6, `LAG-3` = 3.1, `TIM-3` = 2.9, `PD-L1` = 2.0, CD44 = 2.8),
    NK = pop(CD3 = 0.1, `NK1.1` = 3.2, NKp46 = 3.0, NKG2D = 2.6, CD8 = 1.0, GZMB = 2.4, Perforin = 2.2),
    NKT = pop(CD3 = 3.0, `NK1.1` = 2.8, NKp46 = 2.2, CD4 = 1.2, CD44 = 2.4),
    B = pop(CD19 = 3.3, CD3 = 0.1, CD27 = 1.2),
    Myeloid = pop(CD11B = 3.2, CD3 = 0.1, CD19 = 0.1, `NK1.1` = 0.1)
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
    Naive_B = pop(CD19 = 3.2, IgD = 3.0, IgM = 2.6, CD27 = 0.3),
    IgM_memory = pop(CD19 = 3.1, CD27 = 3.0, IgM = 2.9, IgD = 0.3, IgG = 0.2),
    Memory_B = pop(CD19 = 3.1, CD27 = 3.0, IgD = 0.3, IgM = 0.3, IgG = 0.3),
    Switched_B = pop(CD19 = 3.0, IgG = 3.1, CD27 = 2.4, IgD = 0.2, IgM = 0.3),
    Activated_B = pop(CD19 = 3.0, CD80 = 2.8, CD86 = 3.0, CD40 = 2.6, CD27 = 1.5),
    Plasma = pop(CD19 = 1.2, `BLIMP-1` = 3.3, CD27 = 2.8, IgD = 0.2)
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
    Mono_Ly6Chi = pop(CD11B = 3.1, LY6C = 3.2, LY6G = 0.2, CD86 = 1.8, `F4/80` = 0.4),
    Mono_Ly6Clo = pop(CD11B = 3.0, LY6C = 0.3, LY6G = 0.2, `F4/80` = 0.4, CD86 = 1.2),
    Macrophage = pop(CD11B = 3.0, `F4/80` = 3.2, `I-A/I-E` = 1.8, CD86 = 1.6),
    M1_like = pop(CD11B = 3.0, `F4/80` = 3.0, iNOS = 3.2, CD86 = 2.6, CD80 = 2.2),
    M2_like = pop(CD11B = 3.0, `F4/80` = 3.0, CD206 = 3.1, `ARG-1` = 2.8, `IL-10` = 2.2, `TGF-b` = 2.0),
    DC = pop(CD11C = 3.2, `I-A/I-E` = 3.3, CD80 = 2.2, CD86 = 2.4, CD11B = 1.4),
    cDC1 = pop(CD11C = 3.1, CD103 = 3.0, `I-A/I-E` = 3.0, CD11B = 0.6),
    Eosinophil = pop(`Siglec-F` = 3.2, CCR3 = 2.8, CD11B = 2.4, LY6G = 0.2),
    Basophil = pop(FceRI = 3.1, CD200R3 = 2.8, CD11B = 1.6)
  )
}

demo_props <- function(panel_id, group) {
  ctrl <- identical(as.character(group), flow_ctrl_group)
  if (panel_id == "P1") {
    if (ctrl) {
      return(c(CD4_naive = 0.14, CD4_TCM = 0.07, CD4_TEM = 0.07, Treg = 0.05, CD4_act = 0.05,
               CD8_naive = 0.09, CD8_TCM = 0.06, CD8_TEM = 0.07, CD8_eff = 0.05, CD8_exh = 0.04,
               NK = 0.12, NKT = 0.04, B = 0.09, Myeloid = 0.06))
    }
    return(c(CD4_naive = 0.07, CD4_TCM = 0.05, CD4_TEM = 0.07, Treg = 0.05, CD4_act = 0.09,
             CD8_naive = 0.05, CD8_TCM = 0.05, CD8_TEM = 0.06, CD8_eff = 0.08, CD8_exh = 0.09,
             NK = 0.12, NKT = 0.05, B = 0.08, Myeloid = 0.05))
  }
  if (panel_id == "P2") {
    if (ctrl) return(c(Naive_B = 0.38, IgM_memory = 0.10, Memory_B = 0.14, Switched_B = 0.10, Activated_B = 0.16, Plasma = 0.12))
    return(c(Naive_B = 0.22, IgM_memory = 0.10, Memory_B = 0.12, Switched_B = 0.18, Activated_B = 0.20, Plasma = 0.18))
  }
  if (ctrl) {
    return(c(Neutrophil = 0.18, Mono_Ly6Chi = 0.12, Mono_Ly6Clo = 0.08, Macrophage = 0.12, M1_like = 0.06, M2_like = 0.06, DC = 0.12, cDC1 = 0.08, Eosinophil = 0.10, Basophil = 0.08))
  }
  c(Neutrophil = 0.12, Mono_Ly6Chi = 0.08, Mono_Ly6Clo = 0.07, Macrophage = 0.10, M1_like = 0.06, M2_like = 0.14, DC = 0.14, cDC1 = 0.09, Eosinophil = 0.12, Basophil = 0.08)
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

run_umap <- function(mat) {
  if (has_pkg("uwot")) {
    emb <- uwot::umap(mat, n_neighbors = 15, min_dist = 0.3, metric = "euclidean",
                      verbose = FALSE, n_threads = 2, seed = seed_value)
    colnames(emb) <- c("UMAP1", "UMAP2")
    return(emb)
  }
  log_msg("uwot 未安装，UMAP 用前两主成分代替（图标题会标明 PCA fallback）")
  pca <- run_pca(mat, 2)$embedding
  colnames(pca) <- c("UMAP1", "UMAP2")
  pca
}

run_tsne <- function(mat) {
  perplexity <- max(5, min(30, floor((nrow(mat) - 1) / 3)))
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
      cd62 <- nv(i, "CD62L")
      cd44 <- nv(i, "CD44")
      if (!is.finite(cd62)) cd62 <- -Inf
      if (!is.finite(cd44)) cd44 <- -Inf
      if (cd62 > cd44 + 0.3) return(paste0(prefix, "_naive"))
      if (cd44 > cd62 + 0.3) return(paste0(prefix, "_TEM"))
      if (is.finite(cd62) || is.finite(cd44)) return(paste0(prefix, "_TCM"))
      paste0(prefix, "_T")
    }
    if (panel_id == "P1") {
      # NK 必须明显 CD3-；不要靠 NK1.1 背景把 T 吃掉
      if (is.finite(cd19) && cd19 >= max(cd3, cd4, cd8, nk) - 0.05 && cd19 > cd3 + 0.15) return("B")
      if (is.finite(cd11b) && cd11b > max(cd3, cd4, cd8, cd19) + 0.2) return("Myeloid")
      if (is.finite(nk) && nk > cd3 + 0.3 && nk > cd19 && nk >= max(cd4, cd8) - 0.15) return("NK")
      if (is.finite(nk) && is.finite(cd3) && nk > 1 && cd3 > 1 &&
          abs(cd3 - nk) < 1.2 && nk > cd19 && cd3 > cd19 && nk >= max(cd4, cd8) - 0.4) {
        return("NKT")
      }
      if (cd8 > cd4 + 0.15) {
        exh <- max(nv(i, "LAG-3"), nv(i, "TIM-3"))
        eff <- max(nv(i, "GZMB"), nv(i, "Perforin"))
        cd62 <- nv(i, "CD62L")
        cd44 <- nv(i, "CD44")
        # IFN-g 背景不能当 effector；GZMB/Perforin 必须压过 CD44/CD62L
        if (is.finite(exh) && exh > max(eff, cd62) + 0.2 && exh >= cd44 - 0.35) return("CD8_exhausted")
        if (is.finite(eff) && eff > max(cd62, cd44) + 0.15 && eff > exh) return("CD8_effector")
        return(t_mem("CD8"))
      }
      if (cd4 > cd8 + 0.15) {
        cd25 <- nv(i, "CD25")
        cd69 <- nv(i, "CD69")
        tnfa <- nv(i, "TNF-a")
        cd62 <- nv(i, "CD62L")
        cd44 <- nv(i, "CD44")
        act <- max(cd69, tnfa)
        if (is.finite(cd25) && cd25 > max(cd69, tnfa, nv(i, "IFN-g")) + 0.2 &&
            cd25 >= max(cd62, cd44) - 0.6) return("Treg")
        # IFN-g 不能单独把 naive/TEM 标成 activated
        if (is.finite(act) && act > max(cd62, cd44) + 0.2 && act > cd25) return("CD4_activated")
        return(t_mem("CD4"))
      }
      if (is.finite(cd3) && cd3 > max(cd19, cd11b, nk)) return(t_mem("T"))
      return("Myeloid")
    }
    if (panel_id == "P2") {
      blimp <- nv(i, "BLIMP-1")
      igd <- nv(i, "IgD")
      igm <- nv(i, "IgM")
      igg <- nv(i, "IgG")
      cd27 <- nv(i, "CD27")
      act <- max(nv(i, "CD80"), nv(i, "CD86"))
      if (is.finite(blimp) && blimp > igd + 0.4 && blimp > igm + 0.25) return("Plasma")
      if (is.finite(igg) && igg > igd + 0.25 && igg > igm + 0.2) return("Switched_B")
      if (is.finite(act) && act > max(igd, igm, igg) + 0.2) return("Activated_B")
      if (is.finite(cd27) && cd27 > igd + 0.25) {
        if (is.finite(igm) && igm > igd + 0.2 && igm >= igg) return("IgM_memory")
        return("Memory_B")
      }
      if (is.finite(igd) && igd >= igm - 0.2) return("Naive_B")
      if (is.finite(igm) && igm > igd + 0.35) return("Atypical_B")
      return("Naive_B")
    }
    # P3：按定义标志先后判断。中性粒可 Ly6C 中阳，不能输给单核；巨噬可 MHCII+，不能输给 DC。
    # 嗜酸只用 Siglec-F（不用 CCR3）；M1/M2 只用 iNOS / CD206 / ARG-1（不用 IL-10/TGF-b）。
    siglec <- nv(i, "Siglec-F")
    ly6g <- nv(i, "LY6G")
    f480 <- nv(i, "F4/80")
    cd11c <- nv(i, "CD11C")
    cd103 <- nv(i, "CD103")
    ly6c <- nv(i, "LY6C")
    cd11b <- nv(i, "CD11B")
    inos <- nv(i, "iNOS")
    m2 <- max(nv(i, "CD206"), nv(i, "ARG-1"))
    baso <- max(nv(i, "FceRI"), nv(i, "CD200R3"))
    cd3 <- nv(i, "CD3")
    cd19 <- nv(i, "CD19")
    nk <- nv(i, "NK1.1")
    myel_def <- max(siglec, ly6g, f480, cd11c, ly6c, inos, m2, baso, cd103)
    lymph <- max(cd3, cd19, nk)
    if (is.finite(cd19) && cd19 >= lymph && cd19 > myel_def + 0.15) return("B")
    if (is.finite(cd3) && cd3 >= lymph && cd3 > myel_def + 0.15) return("T")
    if (is.finite(nk) && nk >= lymph && nk > myel_def + 0.15) return("NK")
    if (is.finite(ly6g) && ly6g >= 1.2 && ly6g > siglec && ly6g > baso && ly6g > cd11c) return("Neutrophil")
    if (is.finite(siglec) && siglec >= 1.2 && siglec > ly6g && siglec > baso) return("Eosinophil")
    if (is.finite(baso) && baso >= 1.2 && baso > ly6g && baso > siglec && baso > cd11c) return("Basophil_mast")
    if (is.finite(f480) && f480 >= 1.2 && f480 > cd11c && f480 > ly6g && f480 > siglec) {
      if (is.finite(inos) && inos > m2 + 0.15) return("M1_like_Mac")
      if (is.finite(m2) && m2 > inos + 0.15) return("M2_like_Mac")
      return("Macrophage")
    }
    if (is.finite(cd103) && cd103 >= 1.2 && cd103 >= cd11c - 0.15 && cd103 > f480 && cd103 > ly6g) {
      return("cDC1_CD103")
    }
    if (is.finite(cd11c) && cd11c >= 1.2 && cd11c > f480 && cd11c > ly6g) return("DC")
    if (is.finite(ly6c) && ly6c >= 1.2 && ly6c > f480 && ly6c > ly6g) return("Mono_Ly6Chi")
    if (is.finite(cd11b) && cd11b > max(cd3, cd19, nk, 0)) return("Mono_Ly6Clo")
    "Myeloid"
  }, character(1))
  data.frame(cluster = rownames(med), lineage = labs, stringsAsFactors = FALSE)
}

colv <- function(mat, name) {
  if (!name %in% colnames(mat)) return(rep(NA_real_, nrow(mat)))
  as.numeric(mat[, name])
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

# 大类圈定之后，只根据这一层的定义标志给亚群命名
label_cd4_subset <- function(v) {
  cd25 <- vec_get(v, "CD25")
  cd69 <- vec_get(v, "CD69")
  tnfa <- vec_get(v, "TNF-a")
  cd62 <- vec_get(v, "CD62L")
  cd44 <- vec_get(v, "CD44")
  act <- max(cd69, tnfa)
  if (cd25 > act + 0.15 && cd25 > max(cd62, cd44) - 0.8) return("Treg")
  if (act > max(cd62, cd44) + 0.1 && act > cd25) return("CD4_activated")
  mem_from_cd62_cd44(cd62, cd44, "CD4")
}

label_cd8_subset <- function(v) {
  exh <- max(vec_get(v, "LAG-3"), vec_get(v, "TIM-3"))
  eff <- max(vec_get(v, "GZMB"), vec_get(v, "Perforin"))
  cd62 <- vec_get(v, "CD62L")
  cd44 <- vec_get(v, "CD44")
  if (exh > eff + 0.1 && exh > cd62 + 0.1) return("CD8_exhausted")
  if (eff > max(cd62, cd44) + 0.1 && eff > exh) return("CD8_effector")
  mem_from_cd62_cd44(cd62, cd44, "CD8")
}

label_t_subset <- function(v) {
  mem_from_cd62_cd44(vec_get(v, "CD62L"), vec_get(v, "CD44"), "T")
}

label_b_subset <- function(v) {
  blimp <- vec_get(v, "BLIMP-1")
  igd <- vec_get(v, "IgD")
  igm <- vec_get(v, "IgM")
  igg <- vec_get(v, "IgG")
  cd27 <- vec_get(v, "CD27")
  act <- max(vec_get(v, "CD80"), vec_get(v, "CD86"))
  if (blimp > igd + 0.4 && blimp > igm + 0.25 && blimp >= cd27 - 0.15) return("Plasma")
  if (igg > igd + 0.2 && igg > igm + 0.15) return("Switched_B")
  if (act > max(igd, igm, igg) + 0.15) return("Activated_B")
  if (cd27 > igd + 0.15) {
    if (igm > igd + 0.1 && igm >= igg) return("IgM_memory")
    return("Memory_B")
  }
  if (igd >= igm - 0.15) return("Naive_B")
  if (igm > igd + 0.2) return("Atypical_B")
  "Naive_B"
}

# P2 第 1 层：只用表面 IgD vs CD27 圈 Naive / Memory。不要把 BLIMP-1 放进这一层，
# 核染色背景会让几乎所有团都变成 Plasma，后面再也分不出亚群。
label_b_major <- function(v) {
  igd <- vec_get(v, "IgD")
  cd27 <- vec_get(v, "CD27")
  if (cd27 > igd + 0.12) return("Memory_B")
  "Naive_B"
}

label_b_naive_subset <- function(v) {
  igd <- vec_get(v, "IgD")
  igm <- vec_get(v, "IgM")
  # CD40 在静息 B 上也有，不能当 Activated
  act <- max(vec_get(v, "CD80"), vec_get(v, "CD86"))
  if (act > igd + 0.25) return("Activated_B")
  if (igm > igd + 0.25) return("Atypical_B")
  "Naive_B"
}

label_b_memory_subset <- function(v) {
  igd <- vec_get(v, "IgD")
  igm <- vec_get(v, "IgM")
  igg <- vec_get(v, "IgG")
  act <- max(vec_get(v, "CD80"), vec_get(v, "CD86"))
  if (igg > igm + 0.15 && igg > igd) return("Switched_B")
  if (igm > igd + 0.1 && igm >= igg) return("IgM_memory")
  if (act > max(igd, igm, igg) + 0.25) return("Activated_B")
  "Memory_B"
}

# P3 第 2 层命名（后备）：先后判断，Ly6C/MHCII 不得把中性粒/巨噬抢走
label_myeloid_major <- function(v) {
  siglec <- vec_get(v, "Siglec-F")
  ly6g <- vec_get(v, "LY6G")
  f480 <- vec_get(v, "F4/80")
  cd11c <- vec_get(v, "CD11C")
  cd103 <- vec_get(v, "CD103")
  ly6c <- vec_get(v, "LY6C")
  baso <- max(vec_get(v, "FceRI"), vec_get(v, "CD200R3"))
  if (ly6g >= 1.2 && ly6g > siglec && ly6g > baso && ly6g > cd11c) return("Neutrophil")
  if (siglec >= 1.2 && siglec > ly6g && siglec > baso) return("Eosinophil")
  if (baso >= 1.2 && baso > ly6g && baso > siglec && baso > cd11c) return("Basophil_mast")
  if (f480 >= 1.2 && f480 > cd11c && f480 > ly6g && f480 > siglec) return("Macrophage")
  if (cd103 >= 1.2 && cd103 >= cd11c - 0.15 && cd103 > f480 && cd103 > ly6g) return("cDC1_CD103")
  if (cd11c >= 1.2 && cd11c > f480 && cd11c > ly6g) return("DC")
  if (ly6c >= 1.2 && ly6c > f480 && ly6c > ly6g) return("Mono_Ly6Chi")
  "Mono_Ly6Clo"
}

# P3 第 3 层：巨噬圈定后再用 iNOS / CD206 / ARG-1
label_mac_cytokine <- function(v) {
  inos <- vec_get(v, "iNOS")
  m2 <- max(vec_get(v, "CD206"), vec_get(v, "ARG-1"))
  if (inos > m2 + 0.15) return("M1_like_Mac")
  if (m2 > inos + 0.15) return("M2_like_Mac")
  "Macrophage"
}

label_dc_cytokine <- function(v) {
  cd103 <- vec_get(v, "CD103")
  cd11c <- vec_get(v, "CD11C")
  if (cd103 >= cd11c - 0.15 && cd103 > vec_get(v, "F4/80")) return("cDC1_CD103")
  "DC"
}

label_myeloid_subset <- function(v) {
  siglec <- vec_get(v, "Siglec-F")
  ly6g <- vec_get(v, "LY6G")
  f480 <- vec_get(v, "F4/80")
  cd11c <- vec_get(v, "CD11C")
  cd103 <- vec_get(v, "CD103")
  ly6c <- vec_get(v, "LY6C")
  inos <- vec_get(v, "iNOS")
  m2 <- max(vec_get(v, "CD206"), vec_get(v, "ARG-1"))
  baso <- max(vec_get(v, "FceRI"), vec_get(v, "CD200R3"))
  if (ly6g >= 1.2 && ly6g > siglec && ly6g > baso && ly6g > cd11c) return("Neutrophil")
  if (siglec >= 1.2 && siglec > ly6g && siglec > baso) return("Eosinophil")
  if (baso >= 1.2 && baso > ly6g && baso > siglec && baso > cd11c) return("Basophil_mast")
  if (f480 >= 1.2 && f480 > cd11c && f480 > ly6g && f480 > siglec) {
    if (inos > m2 + 0.15) return("M1_like_Mac")
    if (m2 > inos + 0.15) return("M2_like_Mac")
    return("Macrophage")
  }
  if (cd103 >= 1.2 && cd103 >= cd11c - 0.15 && cd103 > f480 && cd103 > ly6g) return("cDC1_CD103")
  if (cd11c >= 1.2 && cd11c > f480 && cd11c > ly6g) return("DC")
  if (ly6c >= 1.2 && ly6c > f480 && ly6c > ly6g) return("Mono_Ly6Chi")
  "Mono_Ly6Clo"
}

# P2 第 1 层：比较两个团的 CD27−IgD，不要要求同一细胞 CD27 绝对值压过 IgD
# （IgD 背景高时两个团都是 IgD>CD27，旧规则会把全部打成 Naive B）
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
  if (length(ib) < 40) return(out)
  mk <- intersect(c("IgD", "CD27"), colnames(mat))
  if (!length(mk)) return(out)
  if (length(ib) >= 40 && length(mk) >= 1) {
    x <- mat[ib, mk, drop = FALSE]
    xs <- scale(x)
    xs[!is.finite(xs)] <- 0
    set.seed(seed_value)
    km <- tryCatch(
      stats::kmeans(xs, centers = 2, nstart = 10, iter.max = 250, algorithm = "Lloyd"),
      error = function(e) NULL
    )
    if (!is.null(km) && length(unique(km$cluster)) == 2) {
      sc <- vapply(sort(unique(km$cluster)), function(ci) {
        hit <- km$cluster == ci
        cd27 <- if ("CD27" %in% mk) median(x[hit, "CD27"], na.rm = TRUE) else -Inf
        igd <- if ("IgD" %in% mk) median(x[hit, "IgD"], na.rm = TRUE) else 0
        if (!is.finite(cd27)) cd27 <- -Inf
        if (!is.finite(igd)) igd <- 0
        cd27 - igd
      }, numeric(1))
      if (diff(range(sc)) >= 0.08) {
        mem_cl <- sort(unique(km$cluster))[which.max(sc)]
        out[ib] <- ifelse(km$cluster == mem_cl, "Memory_B", "Naive_B")
        return(out)
      }
    }
  }
  mem <- gate_k2_high(mat, ib, "CD27", 0.12)
  if (any(mem)) {
    out[ib[mem]] <- "Memory_B"
    out[ib[!mem]] <- "Naive_B"
    return(out)
  }
  igd_hi <- gate_k2_high(mat, ib, "IgD", 0.12)
  if (any(igd_hi) && mean(igd_hi) < 0.95) {
    out[ib[igd_hi]] <- "Naive_B"
    out[ib[!igd_hi]] <- "Memory_B"
  }
  out
}

# 第 1 层：只用谱系抗体圈大类（P1 T/NK，P2 B 的 Naive/Memory，P3 淋巴 vs 髓系）
gate_major_lineage <- function(mat, panel_id) {
  mat <- as.matrix(mat)
  n <- nrow(mat)
  if (panel_id == "P2") {
    return(gate_p2_major(mat))
  }
  cd3 <- colv(mat, "CD3")
  cd4 <- colv(mat, "CD4")
  cd8 <- finite_pmax(colv(mat, "CD8"), colv(mat, "CD8b"))
  cd19 <- colv(mat, "CD19")
  cd11b <- colv(mat, "CD11B")
  nk <- finite_pmax(colv(mat, "NK1.1"), colv(mat, "NKp46"))
  out <- rep("Myeloid", n)
  is_b <- cd19 > finite_pmax(cd3, cd4, cd8, nk) + 0.1
  is_my <- !is_b & cd11b > finite_pmax(cd3, cd19, cd4, cd8) + 0.2
  is_nk <- !is_b & !is_my & nk > cd3 + 0.3 & nk > cd19
  is_nkt <- !is_b & !is_my & !is_nk & nk > 1 & cd3 > 1 & abs(cd3 - nk) < 1.2 & nk > cd19
  is_t <- !is_b & !is_my & !is_nk & !is_nkt & cd3 >= finite_pmax(cd19, cd11b)
  out[is_b] <- "B"
  out[is_my] <- "Myeloid"
  out[is_nk] <- "NK"
  out[is_nkt] <- "NKT"
  # T 细胞必须判到 CD4 或 CD8，不要留下 T TEM
  out[is_t] <- ifelse(cd4 >= cd8, "CD4", "CD8")
  if (panel_id == "P3") {
    out[out %in% c("CD4", "CD8", "NKT")] <- "T"
  }
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

# 一层：只用 CD62L vs CD44，k=3 按 (CD62L-CD44) 高→低标 naive / TCM / TEM
split_memory_3 <- function(mat, idx, prefix) {
  n <- length(idx)
  fallback <- paste0(prefix, "_TEM")
  if (n == 0) return(character(0))
  mk <- intersect(c("CD62L", "CD44"), colnames(mat))
  if (n < 30 || length(mk) < 2) {
    cd62 <- median(colv(mat, "CD62L")[idx], na.rm = TRUE)
    cd44 <- median(colv(mat, "CD44")[idx], na.rm = TRUE)
    return(rep(mem_from_cd62_cd44(cd62, cd44, prefix), n))
  }
  x <- mat[idx, mk, drop = FALSE]
  xs <- scale(x)
  xs[!is.finite(xs)] <- 0
  k_use <- if (n >= 90) 3L else 2L
  set.seed(seed_value)
  km <- tryCatch(stats::kmeans(xs, centers = k_use, nstart = 10, iter.max = 250, algorithm = "Lloyd"),
                 error = function(e) NULL)
  if (is.null(km)) return(rep(fallback, n))
  cl <- sort(unique(km$cluster))
  delta <- vapply(cl, function(ci) {
    hit <- km$cluster == ci
    median(x[hit, "CD62L"], na.rm = TRUE) - median(x[hit, "CD44"], na.rm = TRUE)
  }, numeric(1))
  ord <- cl[order(delta, decreasing = TRUE)]
  map <- setNames(rep(fallback, length(ord)), as.character(ord))
  map[as.character(ord[1])] <- paste0(prefix, "_naive")
  map[as.character(ord[length(ord)])] <- paste0(prefix, "_TEM")
  if (length(ord) == 3) map[as.character(ord[2])] <- paste0(prefix, "_TCM")
  labs <- unname(map[as.character(km$cluster)])
  labs[is.na(labs) | !nzchar(labs)] <- fallback
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
    take_high("CD25", "Treg", 0.15)
    take_high(c("CD69", "TNF-a"), "CD4_activated", 0.15)
    if (any(remain)) labs[remain] <- split_memory_3(mat, idx[remain], "CD4")
  } else {
    take_high(c("LAG-3", "TIM-3"), "CD8_exhausted", 0.15)
    take_high(c("GZMB", "Perforin"), "CD8_effector", 0.15)
    if (any(remain)) labs[remain] <- split_memory_3(mat, idx[remain], "CD8")
  }
  labs
}

# P2 Naive 内：只用 CD80/CD86 圈少数 Activated（不要用 CD40）；剩余 IgD vs IgM → Naive / Atypical
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
  take_high(c("CD80", "CD86"), "Activated_B", 0.5, beat = "IgD", margin = 0.25, max_frac = 0.35)
  if (any(remain)) {
    hi <- gate_k2_high(mat, idx[remain], "IgM", 0.2)
    pos <- which(remain)
    if (any(hi)) {
      hi_idx <- idx[remain][hi]
      igm <- median(colv(mat, "IgM")[hi_idx], na.rm = TRUE)
      igd <- median(colv(mat, "IgD")[hi_idx], na.rm = TRUE)
      if (is.finite(igm) && igm > igd + 0.25) {
        out[pos[hi]] <- "Atypical_B"
        remain[pos[hi]] <- FALSE
      }
    }
    out[remain] <- "Naive_B"
  }
  out
}

# P2 Memory 内：Plasma（BLIMP 少数）→ 少数 CD80/CD86 Activated 岛 → Switched → IgM memory
# 默认剩余是 Memory。CD40 不当激活；Activated 不得在分型前吞掉 Switched/IgM。
sequential_b_memory <- function(mat, idx) {
  n <- length(idx)
  if (n == 0) return(character(0))
  out <- rep("Memory_B", n)
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
  take_high("BLIMP-1", "Plasma", 0.4, beat = c("IgD", "IgM"), margin = 0.4, max_frac = 0.4)
  take_high(c("CD80", "CD86"), "Activated_B", 0.5, beat = c("IgD", "IgM", "IgG"), margin = 0.25, max_frac = 0.35)
  take_high("IgG", "Switched_B", 0.15, beat = c("IgD", "IgM"), margin = 0.05)
  if (any(remain)) {
    hi <- gate_k2_high(mat, idx[remain], "IgM", 0.15)
    pos <- which(remain)
    if (any(hi)) {
      hi_idx <- idx[remain][hi]
      igm <- median(colv(mat, "IgM")[hi_idx], na.rm = TRUE)
      igd <- median(colv(mat, "IgD")[hi_idx], na.rm = TRUE)
      igg <- median(colv(mat, "IgG")[hi_idx], na.rm = TRUE)
      if (is.finite(igm) && igm > igd + 0.1 && igm >= igg) {
        out[pos[hi]] <- "IgM_memory"
        remain[pos[hi]] <- FALSE
      }
    }
    out[remain] <- "Memory_B"
  }
  out
}

# P3 巨噬圈定后：只用 iNOS → M1，CD206/ARG-1 → M2；不用 IL-10/TGF-b/TNF-a/CD86
sequential_mac <- function(mat, idx) {
  n <- length(idx)
  if (n == 0) return(character(0))
  out <- rep("Macrophage", n)
  remain <- rep(TRUE, n)
  take_high <- function(markers, label, min_sep, beat = NULL, margin = 0.15) {
    if (!any(remain)) return(invisible())
    hi <- gate_k2_high(mat, idx[remain], markers, min_sep)
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
  take_high("iNOS", "M1_like_Mac", 0.2, beat = c("CD206", "ARG-1"), margin = 0.15)
  take_high(c("CD206", "ARG-1"), "M2_like_Mac", 0.2, beat = "iNOS", margin = 0.15)
  out
}

# P3 髓系内一层层圈：中性粒 → 嗜酸 → 嗜碱 → 巨噬 → CD103 DC → DC → 剩余单核
# 不要一次 kmeans 拿 Ly6C/MHCII 和 Ly6G/F4/80 比大小。
sequential_myeloid <- function(mat, idx) {
  n <- length(idx)
  if (n == 0) return(character(0))
  out <- rep("Mono_Ly6Clo", n)
  remain <- rep(TRUE, n)
  take_high <- function(markers, label, min_sep, beat = NULL, margin = 0.15, min_abs = 1.2) {
    if (!any(remain)) return(invisible())
    hi <- gate_k2_high(mat, idx[remain], markers, min_sep)
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
  take_high("LY6G", "Neutrophil", 0.2, beat = c("Siglec-F", "FceRI", "CD11C"), margin = 0.15)
  take_high("Siglec-F", "Eosinophil", 0.2, beat = c("LY6G", "FceRI"), margin = 0.2)
  take_high(c("FceRI", "CD200R3"), "Basophil_mast", 0.2, beat = c("LY6G", "Siglec-F", "CD11C"), margin = 0.15)
  take_high("F4/80", "Macrophage", 0.2, beat = c("CD11C", "LY6G", "Siglec-F"), margin = 0.1)
  take_high("CD103", "cDC1_CD103", 0.2, beat = c("F4/80", "LY6G", "Siglec-F"), margin = 0.15)
  take_high("CD11C", "DC", 0.2, beat = c("F4/80", "LY6G"), margin = 0.15)
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
  } else if (panel_id == "P2") {
    i_n <- which(major == "Naive_B")
    subset[i_n] <- sequential_b_naive(mat, i_n)
    i_m <- which(major == "Memory_B")
    subset[i_m] <- sequential_b_memory(mat, i_m)
  } else if (panel_id == "P3") {
    im <- which(major == "Myeloid")
    subset[im] <- sequential_myeloid(mat, im)
  }
  subset[is.na(subset) | !nzchar(subset)] <- major[is.na(subset) | !nzchar(subset)]
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

# P1 主图高对比色；P2/P3 亚群也从这套取色，不要做成全紫/全橙
pal_p1_hues <- c(
  "#E74C3C", "#7B52A5", "#5DADE2", "#E8C87A",
  "#C0392B", "#D4A017", "#27AE60", "#E67E22",
  "#922B21", "#F5CBA7", "#A9DFBF", "#1E8449",
  "#145A32", "#7D3C98", "#CA6F1E", "#85C1E9"
)

pal_celltype <- c(
  "CD4 activated" = "#E74C3C",
  "B cell" = "#7B52A5",
  "Macrophage" = "#5DADE2",
  "NKT" = "#E8C87A",
  "NK" = "#C0392B",
  "T" = "#D4A017",
  "CD8 TCM" = "#27AE60",
  "CD4 TCM" = "#E67E22",
  "Treg" = "#922B21",
  "CD4 naive" = "#F5CBA7",
  "CD8 naive" = "#A9DFBF",
  "CD8 TEM" = "#1E8449",
  "CD8 effector" = "#145A32",
  "CD8 exhausted" = "#7D3C98",
  "CD4 TEM" = "#CA6F1E",
  "CD4 T" = "#E69A3C",
  "CD8 T" = "#3D8B40",
  "T naive" = "#F7DC6F",
  "T TCM" = "#F4D03F",
  "T TEM" = "#B7950B",
  "Myeloid" = "#85C1E9",
  "M1-like Mac" = "#1A5276",
  "M2-like Mac" = "#27AE60",
  "Naive B" = "#F5CBA7",
  "Atypical B" = "#E8C87A",
  "IgM memory B" = "#A9DFBF",
  "Memory B" = "#27AE60",
  "Switched B" = "#E67E22",
  "Activated B" = "#E74C3C",
  "Plasma" = "#922B21",
  "Eosinophil" = "#E74C3C",
  "Neutrophil" = "#E67E22",
  "Ly6C hi mono" = "#CA6F1E",
  "Ly6C lo mono" = "#E8C87A",
  "DC" = "#7B52A5",
  "CD103 DC" = "#922B21",
  "Basophil" = "#A9DFBF",
  "Basophil/mast" = "#A9DFBF",
  "other" = "#B0B0B0",
  "Other" = "#B0B0B0"
)

celltype_label <- function(lineage, panel_id) {
  rec <- c(
    B = "B cell",
    Naive_B = "Naive B",
    Atypical_B = "Atypical B",
    IgM_memory = "IgM memory B",
    Memory_B = "Memory B",
    Switched_B = "Switched B",
    Activated_B = "Activated B",
    Plasma = "Plasma",
    CD4_naive = "CD4 naive",
    CD4_TCM = "CD4 TCM",
    CD4_TEM = "CD4 TEM",
    CD4_activated = "CD4 activated",
    Treg = "Treg",
    CD4_T = "CD4 T",
    CD8_naive = "CD8 naive",
    CD8_TCM = "CD8 TCM",
    CD8_TEM = "CD8 TEM",
    CD8_effector = "CD8 effector",
    CD8_exhausted = "CD8 exhausted",
    CD8_T = "CD8 T",
    NK = "NK",
    NKT = "NKT",
    T = "T",
    T_naive = "T naive",
    T_TCM = "T TCM",
    T_TEM = "T TEM",
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
    cDC1_CD103 = "CD103 DC",
    Basophil_mast = "Basophil",
    Basophil = "Basophil",
    Other = "other"
  )
  lab <- as.character(lineage)
  out <- rec[lab]
  unname(ifelse(is.na(out), lab, out))
}

celltype_colors <- function(levels) {
  pal <- pal_celltype[levels]
  miss <- levels[is.na(pal)]
  if (length(miss) > 0) {
    extra <- pal_p1_hues[((seq_along(miss) - 1L) %% length(pal_p1_hues)) + 1L]
    pal[is.na(pal)] <- extra
  }
  names(pal) <- levels
  pal
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
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )
}

# 主图：左 EV、右 H；颜色=细亚群；不要虚线
plot_split_lineage <- function(df, x, y, panel_id, xlab, ylab, title) {
  plot_df <- df
  plot_df$group <- factor(plot_df$group, levels = flow_group_levels)
  plot_df$celltype <- celltype_label(plot_df$lineage, panel_id)
  levs <- unique(plot_df$celltype)
  plot_df$celltype <- factor(plot_df$celltype, levels = levs)
  pal <- celltype_colors(levels(plot_df$celltype))
  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = celltype)) +
    ggplot2::geom_point(size = 0.45, alpha = 0.85, stroke = 0) +
    ggplot2::facet_wrap(~group, ncol = 2) +
    ggplot2::coord_equal() +
    ggplot2::scale_color_manual(values = pal, drop = FALSE) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3.2, alpha = 1))) +
    ggplot2::labs(title = title, x = xlab, y = ylab, color = NULL) +
    theme_split_dr()
}

plot_embedding <- function(df, x, y, color_col, title, point_size = 0.35) {
  plot_df <- df
  if (identical(color_col, "lineage")) {
    plot_df$lineage <- celltype_label(plot_df$lineage, NA)
  }
  nlev <- length(unique(plot_df[[color_col]]))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = .data[[color_col]])) +
    ggplot2::geom_point(size = point_size, alpha = 0.75, stroke = 0) +
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
    ggplot2::geom_point(size = 0.3, alpha = 0.8, stroke = 0) +
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
  ggplot2::ggplot(freq_df, ggplot2::aes(x = cluster, y = percent, fill = group)) +
    ggplot2::stat_summary(fun = mean, geom = "col", position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.8), size = 1.6, alpha = 0.9) +
    ggplot2::scale_fill_manual(values = pal_group) +
    theme_dr() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = title, x = NULL, y = "% of cells (per sample)")
}

plot_freq_bar <- function(sum_df, title) {
  ggplot2::ggplot(sum_df, ggplot2::aes(x = group, y = mean_percent, fill = cluster)) +
    ggplot2::geom_col(position = "stack") +
    theme_dr() +
    ggplot2::labs(title = title, x = NULL, y = "Mean % of cells")
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
  paste0(marker, "-", fl)
}

p_to_star <- function(p) {
  if (!is.finite(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

# ns 时把 P 值写出来：n=3 时柱高看起来不同，t 检验仍常不显著
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
  out <- switch(as.character(parent),
    all = rep(TRUE, n),
    CD3 = cd4 | cd8 | lin %in% c("T", "NKT") | cl %in% c("T", "NKT"),
    CD4 = cd4,
    CD8 = cd8,
    Myeloid = cl == "Myeloid" | lin %in% c(
      "Neutrophil", "Eosinophil", "Basophil_mast", "Basophil", "Macrophage",
      "M1_like_Mac", "M2_like_Mac", "DC", "cDC1_CD103", "Mono_Ly6Chi", "Mono_Ly6Clo", "Myeloid"
    ),
    Macrophage = lin %in% c("Macrophage", "M1_like_Mac", "M2_like_Mac"),
    Memory_B = cl == "Memory_B" | lin %in% c("Memory_B", "IgM_memory", "Switched_B", "Plasma", "Activated_B"),
    Naive_B = cl == "Naive_B" | lin %in% c("Naive_B", "Atypical_B"),
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
  lin == spec$lineage
}

subset_plot_specs <- function(panel_id) {
  mk <- function(lineage, x, y, parent, ylab, use_major = FALSE) {
    list(lineage = lineage, x = x, y = y, parent = parent, ylab = ylab, use_major = use_major)
  }
  if (identical(panel_id, "P1")) {
    return(list(
      mk("NK", "CD3", "NK1.1", "all", "NK cell (%)"),
      mk("NKT", "CD3", "NK1.1", "all", "NKT (%)"),
      mk("B", "CD19", "CD3", "all", "B cell (%)"),
      mk("Myeloid", "CD11B", "CD3", "all", "Myeloid (%)"),
      mk("CD4", "CD8", "CD4", "CD3", "CD4+ T cell in CD3+ (%)", TRUE),
      mk("CD8", "CD8", "CD4", "CD3", "CD8+ T cell in CD3+ (%)", TRUE),
      mk("CD4_naive", "CD62L", "CD44", "CD4", "CD4 naive in CD4+ (%)"),
      mk("CD4_TCM", "CD62L", "CD44", "CD4", "CD4 TCM in CD4+ (%)"),
      mk("CD4_TEM", "CD62L", "CD44", "CD4", "CD4 TEM in CD4+ (%)"),
      mk("CD4_activated", "CD69", "TNF-a", "CD4", "CD4 activated in CD4+ (%)"),
      mk("Treg", "CD25", "CD4", "CD4", "Treg in CD4+ (%)"),
      mk("CD8_naive", "CD62L", "CD44", "CD8", "CD8 naive in CD8+ (%)"),
      mk("CD8_TCM", "CD62L", "CD44", "CD8", "CD8 TCM in CD8+ (%)"),
      mk("CD8_TEM", "CD62L", "CD44", "CD8", "CD8 TEM in CD8+ (%)"),
      mk("CD8_effector", "GZMB", "Perforin", "CD8", "CD8 effector in CD8+ (%)"),
      mk("CD8_exhausted", "LAG-3", "TIM-3", "CD8", "CD8 exhausted in CD8+ (%)")
    ))
  }
  if (identical(panel_id, "P2")) {
    return(list(
      mk("Naive_B", "IgD", "CD27", "all", "Naive B (%)"),
      mk("Atypical_B", "IgM", "IgD", "Naive_B", "Atypical B in naive (%)"),
      mk("Memory_B", "IgD", "CD27", "all", "Memory B (%)"),
      mk("IgM_memory", "IgM", "IgD", "Memory_B", "IgM memory in memory B (%)"),
      mk("Switched_B", "IgG", "IgD", "Memory_B", "Switched B in memory B (%)"),
      mk("Activated_B", "CD80", "CD86", "all", "Activated B (%)"),
      mk("Plasma", "BLIMP-1", "IgD", "Memory_B", "Plasma in memory B (%)")
    ))
  }
  list(
    mk("T", "CD3", "NK1.1", "all", "T cell (%)"),
    mk("B", "CD19", "CD3", "all", "B cell (%)"),
    mk("NK", "CD3", "NK1.1", "all", "NK cell (%)"),
    mk("Neutrophil", "LY6G", "LY6C", "Myeloid", "Neutrophil in myeloid (%)"),
    mk("Eosinophil", "Siglec-F", "LY6G", "Myeloid", "Eosinophil in myeloid (%)"),
    mk("Basophil_mast", "FceRI", "CD200R3", "Myeloid", "Basophil in myeloid (%)"),
    mk("Macrophage", "F4/80", "CD11C", "Myeloid", "Macrophage in myeloid (%)"),
    mk("M1_like_Mac", "iNOS", "CD206", "Macrophage", "M1-like Mac in macrophages (%)"),
    mk("M2_like_Mac", "CD206", "ARG-1", "Macrophage", "M2-like Mac in macrophages (%)"),
    mk("DC", "CD11C", "F4/80", "Myeloid", "DC in myeloid (%)"),
    mk("cDC1_CD103", "CD103", "CD11C", "Myeloid", "CD103 DC in myeloid (%)"),
    mk("Mono_Ly6Chi", "LY6C", "CD11B", "Myeloid", "Ly6C hi mono in myeloid (%)"),
    mk("Mono_Ly6Clo", "LY6C", "CD11B", "Myeloid", "Ly6C lo mono in myeloid (%)")
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
    data.frame(
      sample = s, group = grp,
      n_parent = n_par, n_subset = n_hit,
      percent = if (n_par > 0) 100 * n_hit / n_par else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

gate_rect_from_xy <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 8) {
    return(list(
      xmin = if (any(ok)) min(x[ok]) else 0, xmax = if (any(ok)) max(x[ok]) else 1,
      ymin = if (any(ok)) min(y[ok]) else 0, ymax = if (any(ok)) max(y[ok]) else 1
    ))
  }
  list(
    xmin = as.numeric(stats::quantile(x[ok], 0.10, names = FALSE)),
    xmax = as.numeric(stats::quantile(x[ok], 0.90, names = FALSE)),
    ymin = as.numeric(stats::quantile(y[ok], 0.10, names = FALSE)),
    ymax = as.numeric(stats::quantile(y[ok], 0.90, names = FALSE))
  )
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
                                xlim = NULL, ylim = NULL) {
  d <- df[is.finite(df[[x]]) & is.finite(df[[y]]), c(x, y), drop = FALSE]
  names(d) <- c("x", "y")
  if (nrow(d) > 4000) {
    set.seed(if (exists("seed_value")) seed_value else 42)
    d <- d[sample.int(nrow(d), 4000), , drop = FALSE]
  }
  if (is.null(xlim)) xlim <- finite_axis_lim(c(d$x, gate$xmin, gate$xmax))
  if (is.null(ylim)) ylim <- finite_axis_lim(c(d$y, gate$ymin, gate$ymax))
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
  p +
    ggplot2::annotate(
      "rect", xmin = gate$xmin, xmax = gate$xmax, ymin = gate$ymin, ymax = gate$ymax,
      color = color, fill = NA, linewidth = 0.55
    ) +
    ggplot2::annotate(
      "text", x = gate$xmax, y = gate$ymax,
      label = sprintf("%.2f%%", pct), hjust = 1.08, vjust = -0.35,
      color = color, fontface = "bold", size = 3.4
    ) +
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
    if (!all(c(spec$x, spec$y) %in% names(cells))) next
    par <- parent_mask(cells, spec$parent)
    hit <- subset_hit_mask(cells, spec)
    par[is.na(par)] <- FALSE
    hit[is.na(hit)] <- FALSE
    n_par <- sum(par)
    n_hit <- sum(par & hit)
    if (!is.finite(n_par) || !is.finite(n_hit) || n_par < 20 || n_hit < 8) next
    samp <- subset_sample_percent(cells, spec)
    if (is.null(samp) || nrow(samp) < 2) next
    ctrl_v <- samp$percent[as.character(samp$group) == flow_ctrl_group]
    trt_v <- samp$percent[as.character(samp$group) == flow_trt_group]
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
    bar <- plot_subset_stat_bar(samp, ylab, pv)
    xy_hit <- cells[par & hit, , drop = FALSE]
    gate <- gate_rect_from_xy(xy_hit[[spec$x]], xy_hit[[spec$y]])
    xlab <- axis_fl_label(panel_id, spec$x)
    yfl <- axis_fl_label(panel_id, spec$y)
    pct_of <- function(g) {
      v <- samp$percent[as.character(samp$group) == g]
      if (!length(v) || all(!is.finite(v))) return(0)
      mean(v, na.rm = TRUE)
    }
    d_ctrl <- cells[par & as.character(cells$group) == flow_ctrl_group, , drop = FALSE]
    d_trt <- cells[par & as.character(cells$group) == flow_trt_group, , drop = FALSE]
    if (nrow(d_ctrl) < 8 || nrow(d_trt) < 8) next
    col_ctrl <- unname(pal_group[flow_ctrl_group])
    col_trt <- unname(pal_group[flow_trt_group])
    if (is.na(col_ctrl)) col_ctrl <- "#1A1A1A"
    if (is.na(col_trt)) col_trt <- "#E31A1C"
    xlim <- finite_axis_lim(c(d_ctrl[[spec$x]], d_trt[[spec$x]], gate$xmin, gate$xmax))
    ylim <- finite_axis_lim(c(d_ctrl[[spec$y]], d_trt[[spec$y]], gate$ymin, gate$ymax))
    c_ctrl <- plot_subset_contour(d_ctrl, spec$x, spec$y, col_ctrl, xlab, yfl, pct_of(flow_ctrl_group), gate, xlim, ylim)
    c_trt <- plot_subset_contour(d_trt, spec$x, spec$y, col_trt, xlab, yfl, pct_of(flow_trt_group), gate, xlim, ylim)
    stub <- paste0(panel_id, "_", gsub("[^A-Za-z0-9]+", "_", spec$lineage), "_H_vs_EV")
    save_subset_figure(bar, c_ctrl, c_trt, file.path(sub_dir, stub))
    utils::write.csv(samp, file.path(sub_dir, paste0(stub, "_by_sample.csv")), row.names = FALSE)
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
      ylab = ylab,
      stringsAsFactors = FALSE
    )
  }
  if (length(stat_rows)) {
    st <- do.call(rbind, stat_rows)
    st$padj <- if (all(is.na(st$p_value))) NA_real_ else p.adjust(st$p_value, method = "BH")
    utils::write.csv(st, file.path(sub_dir, paste0(panel_id, "_subset_H_vs_EV_stats.csv")), row.names = FALSE)
  }
  log_msg(panel_id, " subset stat+contour figures: ", sub_dir)
  invisible(TRUE)
}

export_dimred_plots <- function(cells, med, annot, freq_df, panel_id, out_dir, umap_is_pca = FALSE, tsne_is_pca = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  marker_dir <- file.path(out_dir, "markers")
  dir.create(marker_dir, showWarnings = FALSE)
  umap_lab <- if (umap_is_pca) "PCA (UMAP fallback)" else "UMAP"
  tsne_lab <- if (tsne_is_pca) "PCA (tSNE fallback)" else "tSNE"
  tag <- paste0(panel_id, " H vs EV")
  split_ttl <- paste0(panel_id, "  EV | H")

  save_gg(
    plot_split_lineage(cells, "tSNE1", "tSNE2", panel_id, "tSNE-1", "tSNE-2", split_ttl),
    file.path(out_dir, paste0(panel_id, "_H_vs_EV_tSNE_lineage_split")),
    width = 12.4, height = 5.6
  )
  save_gg(
    plot_split_lineage(cells, "UMAP1", "UMAP2", panel_id, "UMAP-1", "UMAP-2", split_ttl),
    file.path(out_dir, paste0(panel_id, "_H_vs_EV_UMAP_lineage_split")),
    width = 12.4, height = 5.6
  )

  save_gg(plot_embedding(cells, "UMAP1", "UMAP2", "group", paste(tag, "-", umap_lab, "by group")),
          file.path(out_dir, paste0(panel_id, "_UMAP_by_group")))
  save_gg(plot_embedding(cells, "UMAP1", "UMAP2", "sample", paste(tag, "-", umap_lab, "by sample")),
          file.path(out_dir, paste0(panel_id, "_UMAP_by_sample")), width = 8)
  save_gg(plot_embedding(cells, "UMAP1", "UMAP2", "cluster", paste(tag, "-", umap_lab, "by cluster")),
          file.path(out_dir, paste0(panel_id, "_UMAP_by_cluster")), width = 8)
  save_gg(plot_embedding(cells, "UMAP1", "UMAP2", "lineage", paste(tag, "-", umap_lab, "by lineage")),
          file.path(out_dir, paste0(panel_id, "_UMAP_by_lineage")), width = 8)
  save_gg(plot_density_split(cells, "UMAP1", "UMAP2", paste(tag, "-", umap_lab, "density H vs EV")),
          file.path(out_dir, paste0(panel_id, "_UMAP_density_H_vs_EV")), width = 10, height = 5)

  save_gg(plot_embedding(cells, "tSNE1", "tSNE2", "group", paste(tag, "-", tsne_lab, "by group")),
          file.path(out_dir, paste0(panel_id, "_tSNE_by_group")))
  save_gg(plot_embedding(cells, "tSNE1", "tSNE2", "sample", paste(tag, "-", tsne_lab, "by sample")),
          file.path(out_dir, paste0(panel_id, "_tSNE_by_sample")), width = 8)
  save_gg(plot_embedding(cells, "tSNE1", "tSNE2", "cluster", paste(tag, "-", tsne_lab, "by cluster")),
          file.path(out_dir, paste0(panel_id, "_tSNE_by_cluster")), width = 8)
  save_gg(plot_embedding(cells, "tSNE1", "tSNE2", "lineage", paste(tag, "-", tsne_lab, "by lineage")),
          file.path(out_dir, paste0(panel_id, "_tSNE_by_lineage")), width = 8)
  save_gg(plot_density_split(cells, "tSNE1", "tSNE2", paste(tag, "-", tsne_lab, "density H vs EV")),
          file.path(out_dir, paste0(panel_id, "_tSNE_density_H_vs_EV")), width = 10, height = 5)

  dr_cols <- colnames(med)
  for (mk in dr_cols) {
    if (!mk %in% names(cells)) next
    save_gg(plot_marker_embedding(cells, "UMAP1", "UMAP2", mk, paste(tag, "UMAP", mk)),
            file.path(marker_dir, paste0(panel_id, "_UMAP_", gsub("[^A-Za-z0-9]+", "_", mk))),
            width = 6.5, height = 5.5)
  }

  plot_cluster_heatmap(med, paste(tag, "cluster median markers (z)"),
                       file.path(out_dir, paste0(panel_id, "_cluster_marker_heatmap")))

  freq_df$cluster <- factor(freq_df$cluster, levels = annot$cluster)
  save_gg(plot_freq_box(freq_df, paste(tag, "cluster frequency")),
          file.path(out_dir, paste0(panel_id, "_cluster_frequency_H_vs_EV")),
          width = max(8, 0.7 * length(unique(freq_df$cluster)) + 2), height = 5.5)

  lin_freq <- lineage_frequencies(cells)
  names(lin_freq)[names(lin_freq) == "lineage"] <- "cluster"
  save_gg(plot_freq_box(lin_freq, paste(tag, "lineage frequency")),
          file.path(out_dir, paste0(panel_id, "_lineage_frequency_H_vs_EV")),
          width = max(7, 0.8 * length(unique(lin_freq$cluster)) + 2), height = 5.5)

  mean_df <- aggregate(percent ~ group + cluster, data = freq_df, FUN = mean)
  names(mean_df)[names(mean_df) == "percent"] <- "mean_percent"
  save_gg(plot_freq_bar(mean_df, paste(tag, "mean cluster composition")),
          file.path(out_dir, paste0(panel_id, "_cluster_composition_stacked")),
          width = 6, height = 5.5)

  if (has_pkg("cowplot")) {
    p1 <- plot_embedding(cells, "UMAP1", "UMAP2", "group", "Group")
    p2 <- plot_embedding(cells, "UMAP1", "UMAP2", "lineage", "Lineage")
    p3 <- plot_freq_box(lin_freq, "Lineage %")
    overview <- cowplot::plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(1, 1.1, 1.2))
    title <- cowplot::ggdraw() +
      cowplot::draw_label(paste(tag, "dimensionality reduction overview"), fontface = "bold", size = 14)
    overview <- cowplot::plot_grid(title, overview, ncol = 1, rel_heights = c(0.08, 1))
    save_gg(overview, file.path(out_dir, paste0(panel_id, "_H_vs_EV_dimred_overview")),
            width = 16, height = 5.8)
  }
  tryCatch(
    export_subset_gate_figures(cells, panel_id, out_dir),
    error = function(e) log_msg(panel_id, " subset stat+contour figures failed: ", e$message)
  )
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# 8. 频率统计 H vs EV
# -----------------------------------------------------------------------------
# 旧表里的 "T" 单独成列时 read.csv 会变成 TRUE
read_embed_csv <- function(path) {
  df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  chr <- intersect(c("sample", "group", "cluster", "cluster_lineage", "lineage"), names(df))
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
  smp_group <- unique(cells[, c("sample", "group")])
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
  smp_group <- unique(cells[, c("sample", "group")])
  tab <- merge(tab, smp_group, by = "sample")
  tot <- aggregate(n ~ sample, data = tab, FUN = sum)
  names(tot)[2] <- "total"
  tab <- merge(tab, tot, by = "sample")
  tab$percent <- 100 * tab$n / pmax(tab$total, 1)
  tab
}

compare_group_freq <- function(freq_df, id_col = "cluster") {
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
      stringsAsFactors = FALSE
    )
    names(row) <- c(
      id_col,
      paste0("n_", g_ctrl), paste0("n_", g_trt),
      paste0("mean_", g_ctrl), paste0("mean_", g_trt),
      paste0("sd_", g_ctrl), paste0("sd_", g_trt),
      paste0("delta_", g_trt, "_minus_", g_ctrl),
      "p_value"
    )
    row
  })
  out <- do.call(rbind, rows)
  out$padj <- if (all(is.na(out$p_value))) NA_real_ else p.adjust(out$p_value, method = "BH")
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
  meta_true <- character()
  map_used <- NULL

  if (use_demo) {
    log_msg(panel_id, " : DEMO synthetic cells (not real FCS)")
    samples <- expand.grid(group = flow_group_levels, replicate = 1:3, stringsAsFactors = FALSE)
    samples$sample <- paste0(samples$group, samples$replicate)
    for (i in seq_len(nrow(samples))) {
      set.seed(seed_value + i + as.integer(factor(panel_id)) * 10)
      d <- make_demo_sample(panel_id, samples$group[i], samples$sample[i], min(n_cap, 2500))
      all_expr[[i]] <- d$expr[, c("CD45", "L/D", markers), drop = FALSE]
      meta_group <- c(meta_group, rep(samples$group[i], nrow(d$expr)))
      meta_sample <- c(meta_sample, rep(samples$sample[i], nrow(d$expr)))
      meta_true <- c(meta_true, d$true_lineage)
    }
  } else {
    sub <- file_tab[file_tab$panel == panel_id, ]
    if (nrow(sub) == 0) {
      log_msg(panel_id, " : no unmixed files (need EV1_", panel_id, "_unmixed.fcs / H1_", panel_id, "_unmixed.fcs)")
      return(NULL)
    }
    for (i in seq_len(nrow(sub))) {
      log_msg("Read ", sub$file[i])
      rec <- read_fcs_expr(sub$path[i], panel_id)
      keep <- qc_filter_matrix(rec$exprs, rec$names, rec$map)
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
      meta_group <- c(meta_group, rep(sub$group[i], nrow(mat)))
      meta_sample <- c(meta_sample, rep(sub$sample[i], nrow(mat)))
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
    true_lineage = if (length(meta_true)) meta_true else rep(NA_character_, nrow(tr)),
    markers = markers,
    map = map_used
  )
}

analyze_one_panel <- function(panel_id, file_tab, use_demo) {
  log_msg("==== Panel ", panel_id, " : ", panel_map$panels[[panel_id]]$focus, " ====")
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
  true_lin <- dat$true_lineage[ok_row]
  mat <- scale_markers(mat_raw)

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
  hier <- hierarchical_gate(as.matrix(mat_raw), panel_id)
  cells$cluster_lineage <- hier$major
  cells$lineage <- hier$subset
  maj_n <- sort(table(hier$major), decreasing = TRUE)
  lin_n <- sort(table(hier$subset), decreasing = TRUE)
  log_msg(panel_id, " layer1 major: ", paste(paste(names(maj_n), as.integer(maj_n), sep = "="), collapse = ", "))
  log_msg(panel_id, " layer2 subset: ", paste(paste(names(lin_n), as.integer(lin_n), sep = "="), collapse = ", "))
  if (use_demo) cells$true_lineage <- true_lin

  freq_df <- cluster_frequencies(cells)
  stats_cl <- compare_group_freq(freq_df, "cluster")
  lin_freq <- lineage_frequencies(cells)
  stats_lin <- compare_group_freq(lin_freq, "lineage")

  utils::write.csv(annot, file.path(out_dir, paste0(panel_id, "_cluster_annotation.csv")), row.names = FALSE)
  utils::write.csv(cbind(cluster = rownames(med), as.data.frame(med)),
                   file.path(out_dir, paste0(panel_id, "_cluster_median_markers.csv")), row.names = FALSE)
  utils::write.csv(freq_df, file.path(out_dir, paste0(panel_id, "_cluster_frequency_by_sample.csv")), row.names = FALSE)
  utils::write.csv(stats_cl, file.path(out_dir, paste0(panel_id, "_cluster_H_vs_EV_stats.csv")), row.names = FALSE)
  utils::write.csv(stats_lin, file.path(out_dir, paste0(panel_id, "_lineage_H_vs_EV_stats.csv")), row.names = FALSE)
  if (!is.null(dat$map)) {
    utils::write.csv(dat$map, file.path(out_dir, paste0(panel_id, "_channel_map.csv")), row.names = FALSE)
  }
  embed_cols <- c("sample", "group", "cluster", "cluster_lineage", "lineage",
                  "UMAP1", "UMAP2", "tSNE1", "tSNE2")
  extra_cols <- setdiff(names(cells), c(embed_cols, "true_lineage"))
  embed_out <- cells[, intersect(c(embed_cols, extra_cols), names(cells)), drop = FALSE]
  utils::write.csv(embed_out, file.path(out_dir, paste0(panel_id, "_cell_embeddings.csv")), row.names = FALSE)

  export_dimred_plots(cells, med, annot, freq_df, panel_id, out_dir, umap_is_pca, tsne_is_pca)
  log_msg(panel_id, " plots written to ", out_dir)
  list(annot = annot, stats_cluster = stats_cl, stats_lineage = stats_lin, n = nrow(cells))
}

# -----------------------------------------------------------------------------
# 10. 主流程
# -----------------------------------------------------------------------------
if (identical(toupper(Sys.getenv("FLOW_FUNCTIONS_ONLY", "0")), "1")) {
  log_msg("FLOW_FUNCTIONS_ONLY=1: skip analysis")
} else {

log_msg("Flow dir: ", project_dir)
log_msg("Results: ", result_dir)
file_tab <- list_unmixed_files(project_dir)
use_demo <- demo_flag
if (nrow(file_tab) == 0) {
  if (demo_flag) {
    log_msg("No unmixed FCS; FLOW_DEMO=1 -> synthetic data for plot export")
    use_demo <- TRUE
  } else {
    log_msg("No *_unmixed.fcs in ", project_dir)
    log_msg("Put EV/H P1-P3 unmixed files in ", flow_primary_data_dir, ", or set FLOW_DEMO=1 to export demo plots")
    use_demo <- TRUE
    log_msg("Auto-fallback to DEMO so the script can still export figure templates")
  }
} else {
  log_msg("Found unmixed files:\n", paste(file_tab$file, collapse = "\n"))
  for (pn in c("P1", "P2", "P3")) {
    n_pn <- if (nrow(file_tab)) sum(file_tab$panel == pn) else 0L
    log_msg(pn, " files parsed: ", n_pn)
    if (n_pn == 0) {
      log_msg(pn, " missing. Put files named like EV1_", pn, "_unmixed.fcs or EV1-", pn, "_unmixed.fcs in ", project_dir)
    }
  }
  if (any(file_tab$group == "H" & grepl("^EV", file_tab$file, ignore.case = TRUE))) {
    stop("Filename parser classified an EV file as H; refusing to continue")
  }
  ensure_flowcore()
}

if (use_demo) {
  writeLines("DEMO / synthetic cells. Do not use as biological results.",
             file.path(log_dir, "DEMO_WARNING.txt"))
}

panels <- c("P1", "P2", "P3")
summaries <- list()
for (pn in panels) {
  summaries[[pn]] <- tryCatch(
    analyze_one_panel(pn, file_tab, use_demo),
    error = function(e) {
      log_msg("Panel ", pn, " failed: ", e$message)
      NULL
    }
  )
}

sum_path <- file.path(result_dir, "H_vs_EV_lineage_stats_all_panels.csv")
lin_rows <- lapply(names(summaries), function(pn) {
  x <- summaries[[pn]]
  if (is.null(x) || is.null(x$stats_lineage)) return(NULL)
  cbind(panel = pn, x$stats_lineage)
})
lin_rows <- Filter(Negate(is.null), lin_rows)
if (length(lin_rows) > 0) {
  utils::write.csv(do.call(rbind, lin_rows), sum_path, row.names = FALSE)
  log_msg("Summary table: ", sum_path)
}

extra <- file.path(script_dir, "Flow_dimred_all_subsets.R")
if (file.exists(extra)) {
  tryCatch({
    Sys.setenv(FLOW_ALL_SUBSETS_FROM_PIPELINE = "1")
    sys.source(extra, envir = .GlobalEnv)
    export_all_subsets_analysis(result_dir)
  }, error = function(e) log_msg("all-subsets summary failed: ", e$message))
  Sys.unsetenv("FLOW_ALL_SUBSETS_FROM_PIPELINE")
}

traj <- file.path(script_dir, "Flow_dimred_trajectory.R")
if (file.exists(traj)) {
  tryCatch({
    Sys.setenv(FLOW_TRAJECTORY_FROM_PIPELINE = "1")
    sys.source(traj, envir = .GlobalEnv)
    export_all_panel_trajectories(result_dir)
  }, error = function(e) log_msg("trajectory summary failed: ", e$message))
  Sys.unsetenv("FLOW_TRAJECTORY_FROM_PIPELINE")
}

log_msg("Done. Open results_flow/P1, P2, P3 for UMAP/tSNE; all_subsets for frequencies; P*/trajectory for trees.")
invisible(summaries)

}
