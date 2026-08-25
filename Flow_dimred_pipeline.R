#!/usr/bin/env Rscript
# =============================================================================
# 流式降维：T6 vs T（P1 T/NK，P2 B，P3 髓系）
# 输入：E:/R/flow J 下的 *_unmixed.fcs（不要用 raw）
# 每个 panel 单独联合 UMAP/tSNE，导出 PDF+PNG，并比较 T6 vs T 细胞频率
#
# 用法：
#   setwd("E:/R/flow J")
#   source("Flow_dimred_pipeline.R")
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

resolve_flow_dir <- function() {
  env_dir <- Sys.getenv("FLOW_DIR", unset = "")
  candidates <- c(
    env_dir,
    "E:/R/flow J",
    "E:\\R\\flow J",
    file.path(script_dir, "flow_J"),
    file.path(script_dir, "flow J"),
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
  if (dir.exists("E:/R/flow J")) {
    return(normalizePath("E:/R/flow J", winslash = "/", mustWork = FALSE))
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
# 2. 文件名：先 T6 再 T
# -----------------------------------------------------------------------------
parse_fcs_filename <- function(path) {
  b <- basename(path)
  m <- regexec("^(T6|T)-([123])_(P[123])_(unmixed|raw)\\.fcs$", b, ignore.case = TRUE)
  hit <- regmatches(b, m)[[1]]
  if (length(hit) < 5) return(NULL)
  list(
    file = b,
    path = path,
    group = toupper(hit[2]),
    replicate = hit[3],
    sample = paste0(toupper(hit[2]), "-", hit[3]),
    panel = toupper(hit[4]),
    kind = tolower(hit[5])
  )
}

list_unmixed_files <- function(root) {
  files <- list.files(root, pattern = "(?i)_unmixed\\.fcs$", full.names = TRUE)
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
    "E:/R/flow J",
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
  unique(norm_id(ids[nzchar(as.character(ids))]))
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
  desc[is.na(desc)] <- ""
  combo <- paste(nms, desc)
  combo_id <- norm_id(combo)
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
    cid <- combo_id[idx]
    if (any(al == cid) || any(vapply(al, function(a) grepl(a, cid, fixed = TRUE), logical(1)) & nchar(al) >= 3)) {
      exact <- any(al == cid)
      return(if (exact) 100 else 60)
    }
    0
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
    if (!is.na(best_j) && best >= 60) {
      used[best_j] <- TRUE
      mapping$channel[i] <- nms[best_j]
      mapping$channel_index[i] <- best_j
    }
  }
  mapping
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
  map <- match_channels(nms, desc, panel_id)
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
  if (!any(ok)) stop("没有匹配到任何分析通道")
  mat <- exprs[, ch[ok], drop = FALSE]
  colnames(mat) <- markers[ok]
  mat
}

# 演示数据：可区分的免疫亚群，T6 改变部分亚群比例（仅用于缺 FCS 时出图）
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
    CD4_naive = pop(CD3 = 3.2, CD4 = 3.0, CD8 = 0.1, CD8b = 0.1, CD19 = 0.1, CD62L = 3.0, CD44 = 0.4, CD27 = 2.2),
    CD4_act = pop(CD3 = 3.2, CD4 = 3.0, CD25 = 2.8, CD69 = 3.0, CD44 = 2.6, CD62L = 0.4, `TNF-a` = 2.2, `IFN-g` = 1.6),
    CD8_T = pop(CD3 = 3.3, CD8 = 3.1, CD8b = 2.8, CD4 = 0.1, CD44 = 2.2, GZMB = 1.8, Perforin = 1.5),
    CD8_exh = pop(CD3 = 3.0, CD8 = 3.0, CD8b = 2.6, `LAG-3` = 3.0, `TIM-3` = 2.8, `PD-L1` = 2.0, CD44 = 2.8),
    NK = pop(CD3 = 0.1, `NK1.1` = 3.2, NKp46 = 3.0, NKG2D = 2.6, CD8 = 1.0, GZMB = 2.4, Perforin = 2.2),
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
    Memory_B = pop(CD19 = 3.1, CD27 = 3.0, IgD = 0.4, IgM = 0.8),
    Switched_B = pop(CD19 = 3.0, IgG = 3.1, CD27 = 2.4, IgD = 0.2),
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
    Mono_Ly6Chi = pop(CD11B = 3.1, LY6C = 3.2, LY6G = 0.2, CD86 = 1.8),
    Macrophage = pop(CD11B = 3.0, `F4/80` = 3.2, `I-A/I-E` = 1.8, CD86 = 1.6),
    M2_like = pop(CD11B = 3.0, `F4/80` = 3.0, CD206 = 3.1, `ARG-1` = 2.8, `IL-10` = 2.2, `TGF-b` = 2.0),
    DC = pop(CD11C = 3.2, `I-A/I-E` = 3.3, CD80 = 2.2, CD86 = 2.4, CD11B = 1.4),
    cDC1 = pop(CD11C = 3.1, CD103 = 3.0, `I-A/I-E` = 3.0, CD11B = 0.6),
    Eosinophil = pop(`Siglec-F` = 3.2, CCR3 = 2.8, CD11B = 2.4),
    Basophil = pop(FceRI = 3.1, CD200R3 = 2.8, CD11B = 1.6)
  )
}

demo_props <- function(panel_id, group) {
  if (panel_id == "P1") {
    if (group == "T") return(c(CD4_naive = 0.28, CD4_act = 0.08, CD8_T = 0.22, CD8_exh = 0.06, NK = 0.16, B = 0.12, Myeloid = 0.08))
    return(c(CD4_naive = 0.16, CD4_act = 0.16, CD8_T = 0.16, CD8_exh = 0.16, NK = 0.18, B = 0.10, Myeloid = 0.08))
  }
  if (panel_id == "P2") {
    if (group == "T") return(c(Naive_B = 0.45, Memory_B = 0.18, Switched_B = 0.12, Activated_B = 0.15, Plasma = 0.10))
    return(c(Naive_B = 0.28, Memory_B = 0.16, Switched_B = 0.18, Activated_B = 0.20, Plasma = 0.18))
  }
  if (group == "T") {
    return(c(Neutrophil = 0.22, Mono_Ly6Chi = 0.16, Macrophage = 0.16, M2_like = 0.08, DC = 0.14, cDC1 = 0.08, Eosinophil = 0.10, Basophil = 0.06))
  }
  c(Neutrophil = 0.14, Mono_Ly6Chi = 0.12, Macrophage = 0.14, M2_like = 0.18, DC = 0.16, cDC1 = 0.10, Eosinophil = 0.10, Basophil = 0.06)
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
    ts <- Rtsne::Rtsne(mat, perplexity = perplexity, verbose = FALSE, check_duplicates = FALSE, pca = TRUE)
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
  switch(panel_id, P1 = 10, P2 = 6, P3 = 10, 8)
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
  km <- kmeans(mat, centers = k, nstart = 10, iter.max = 100)
  factor(paste0("C", km$cluster), levels = paste0("C", sort(unique(km$cluster))))
}

annotate_clusters <- function(med, panel_id) {
  labs <- vapply(seq_len(nrow(med)), function(i) {
    hi <- function(m, cut = 1.5) m %in% colnames(med) && med[i, m] > cut
    lo <- function(m, cut = 1.0) m %in% colnames(med) && med[i, m] < cut
    if (panel_id == "P1") {
      if (hi("CD19") && lo("CD3")) return("B")
      if ((hi("NK1.1") || hi("NKp46")) && lo("CD3")) return("NK")
      if (hi("CD11B") && lo("CD3") && lo("CD19")) return("Myeloid")
      if (hi("CD3") && (hi("NK1.1") || hi("NKp46"))) return("NKT")
      if (hi("CD3") && (hi("CD8") || hi("CD8b")) && (hi("LAG-3") || hi("TIM-3"))) return("CD8_exhausted")
      if (hi("CD3") && (hi("CD8") || hi("CD8b"))) return("CD8_T")
      if (hi("CD3") && hi("CD4") && (hi("CD69") || hi("CD25"))) return("CD4_activated")
      if (hi("CD3") && hi("CD4") && hi("CD62L") && lo("CD44")) return("CD4_naive")
      if (hi("CD3") && hi("CD4")) return("CD4_T")
      if (hi("CD3")) return("T")
      return("Other")
    }
    if (panel_id == "P2") {
      if (hi("BLIMP-1")) return("Plasma")
      if (hi("IgG")) return("Switched_B")
      if (hi("CD80") || hi("CD86")) return("Activated_B")
      if (hi("CD27") && lo("IgD")) return("Memory_B")
      if (hi("IgD") || hi("IgM")) return("Naive_B")
      if (hi("CD19")) return("B")
      return("Other")
    }
    if (hi("Siglec-F") || hi("CCR3")) return("Eosinophil")
    if (hi("FceRI") || hi("CD200R3")) return("Basophil_mast")
    if (hi("CD103") && hi("CD11C")) return("cDC1_CD103")
    if (hi("CD11C") && hi("I-A/I-E")) return("DC")
    if (hi("CD206") || hi("ARG-1")) return("M2_like_Mac")
    if (hi("iNOS") && hi("F4/80")) return("M1_like_Mac")
    if (hi("F4/80")) return("Macrophage")
    if (hi("LY6G")) return("Neutrophil")
    if (hi("LY6C") && hi("CD11B")) return("Mono_Ly6Chi")
    if (hi("CD11B")) return("Myeloid")
    "Other"
  }, character(1))
  data.frame(cluster = rownames(med), lineage = labs, stringsAsFactors = FALSE)
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

pal_group <- c(T = "#4C78A8", T6 = "#E45756")

plot_embedding <- function(df, x, y, color_col, title, point_size = 0.35) {
  nlev <- length(unique(df[[color_col]]))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = .data[[color_col]])) +
    ggplot2::geom_point(size = point_size, alpha = 0.75, stroke = 0) +
    ggplot2::coord_equal() +
    theme_dr() +
    ggplot2::labs(title = title, color = color_col, x = x, y = y)
  if (color_col == "group") {
    p <- p + ggplot2::scale_color_manual(values = pal_group)
  } else if (is.numeric(df[[color_col]])) {
    p <- p + ggplot2::scale_color_gradientn(colours = c("#440154", "#21908C", "#FDE725"))
  } else {
    pal <- grDevices::hcl.colors(max(nlev, 3), palette = "Dark 3")[seq_len(nlev)]
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

export_dimred_plots <- function(cells, med, annot, freq_df, panel_id, out_dir, umap_is_pca = FALSE, tsne_is_pca = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  marker_dir <- file.path(out_dir, "markers")
  dir.create(marker_dir, showWarnings = FALSE)
  umap_lab <- if (umap_is_pca) "PCA (UMAP fallback)" else "UMAP"
  tsne_lab <- if (tsne_is_pca) "PCA (tSNE fallback)" else "tSNE"
  tag <- paste0(panel_id, " T6 vs T")

  save_gg(plot_embedding(cells, "UMAP1", "UMAP2", "group", paste(tag, "-", umap_lab, "by group")),
          file.path(out_dir, paste0(panel_id, "_UMAP_by_group")))
  save_gg(plot_embedding(cells, "UMAP1", "UMAP2", "sample", paste(tag, "-", umap_lab, "by sample")),
          file.path(out_dir, paste0(panel_id, "_UMAP_by_sample")), width = 8)
  save_gg(plot_embedding(cells, "UMAP1", "UMAP2", "cluster", paste(tag, "-", umap_lab, "by cluster")),
          file.path(out_dir, paste0(panel_id, "_UMAP_by_cluster")), width = 8)
  save_gg(plot_embedding(cells, "UMAP1", "UMAP2", "lineage", paste(tag, "-", umap_lab, "by lineage")),
          file.path(out_dir, paste0(panel_id, "_UMAP_by_lineage")), width = 8)
  save_gg(plot_density_split(cells, "UMAP1", "UMAP2", paste(tag, "-", umap_lab, "density T vs T6")),
          file.path(out_dir, paste0(panel_id, "_UMAP_density_T_vs_T6")), width = 10, height = 5)

  save_gg(plot_embedding(cells, "tSNE1", "tSNE2", "group", paste(tag, "-", tsne_lab, "by group")),
          file.path(out_dir, paste0(panel_id, "_tSNE_by_group")))
  save_gg(plot_embedding(cells, "tSNE1", "tSNE2", "sample", paste(tag, "-", tsne_lab, "by sample")),
          file.path(out_dir, paste0(panel_id, "_tSNE_by_sample")), width = 8)
  save_gg(plot_embedding(cells, "tSNE1", "tSNE2", "cluster", paste(tag, "-", tsne_lab, "by cluster")),
          file.path(out_dir, paste0(panel_id, "_tSNE_by_cluster")), width = 8)
  save_gg(plot_embedding(cells, "tSNE1", "tSNE2", "lineage", paste(tag, "-", tsne_lab, "by lineage")),
          file.path(out_dir, paste0(panel_id, "_tSNE_by_lineage")), width = 8)
  save_gg(plot_density_split(cells, "tSNE1", "tSNE2", paste(tag, "-", tsne_lab, "density T vs T6")),
          file.path(out_dir, paste0(panel_id, "_tSNE_density_T_vs_T6")), width = 10, height = 5)

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
          file.path(out_dir, paste0(panel_id, "_cluster_frequency_T6_vs_T")),
          width = max(8, 0.7 * length(unique(freq_df$cluster)) + 2), height = 5.5)

  lin_freq <- aggregate(percent ~ sample + group + lineage, data = freq_df, FUN = sum)
  names(lin_freq)[names(lin_freq) == "lineage"] <- "cluster"
  save_gg(plot_freq_box(lin_freq, paste(tag, "lineage frequency")),
          file.path(out_dir, paste0(panel_id, "_lineage_frequency_T6_vs_T")),
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
    save_gg(overview, file.path(out_dir, paste0(panel_id, "_T6_vs_T_dimred_overview")),
            width = 16, height = 5.8)
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# 8. 频率统计 T6 vs T
# -----------------------------------------------------------------------------
cluster_frequencies <- function(cells) {
  tab <- as.data.frame(table(sample = cells$sample, cluster = cells$cluster), stringsAsFactors = FALSE)
  names(tab)[3] <- "n"
  smp_group <- unique(cells[, c("sample", "group")])
  tab <- merge(tab, smp_group, by = "sample")
  tot <- aggregate(n ~ sample, data = tab, FUN = sum)
  names(tot)[2] <- "total"
  tab <- merge(tab, tot, by = "sample")
  tab$percent <- 100 * tab$n / pmax(tab$total, 1)
  lin_map <- unique(cells[, c("cluster", "lineage")])
  tab <- merge(tab, lin_map, by = "cluster", all.x = TRUE)
  tab
}

compare_group_freq <- function(freq_df, id_col = "cluster") {
  ids <- unique(as.character(freq_df[[id_col]]))
  rows <- lapply(ids, function(id) {
    sub <- freq_df[as.character(freq_df[[id_col]]) == id, ]
    t_vals <- sub$percent[sub$group == "T"]
    t6_vals <- sub$percent[sub$group == "T6"]
    if (length(t_vals) < 2 || length(t6_vals) < 2) {
      pv <- NA_real_
    } else {
      pv <- tryCatch(t.test(t6_vals, t_vals)$p.value, error = function(e) NA_real_)
    }
    data.frame(
      id = id,
      n_T = length(t_vals),
      n_T6 = length(t6_vals),
      mean_T = mean(t_vals),
      mean_T6 = mean(t6_vals),
      sd_T = sd(t_vals),
      sd_T6 = sd(t6_vals),
      delta_T6_minus_T = mean(t6_vals) - mean(t_vals),
      p_value = pv,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$padj <- if (all(is.na(out$p_value))) NA_real_ else p.adjust(out$p_value, method = "BH")
  names(out)[1] <- id_col
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
    samples <- expand.grid(group = c("T", "T6"), replicate = 1:3, stringsAsFactors = FALSE)
    samples$sample <- paste0(samples$group, "-", samples$replicate)
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
      log_msg(panel_id, " : no unmixed files")
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
  cl <- cluster_cells(pca_use, panel_id)

  cells <- data.frame(
    sample = smp,
    group = factor(grp, levels = c("T", "T6")),
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
  cells$lineage <- annot$lineage[match(as.character(cells$cluster), annot$cluster)]
  if (use_demo) cells$true_lineage <- true_lin

  freq_df <- cluster_frequencies(cells)
  freq_df$lineage <- annot$lineage[match(as.character(freq_df$cluster), annot$cluster)]
  stats_cl <- compare_group_freq(freq_df, "cluster")
  lin_freq <- aggregate(percent ~ sample + group + lineage, data = freq_df, FUN = sum)
  stats_lin <- compare_group_freq(lin_freq, "lineage")

  utils::write.csv(annot, file.path(out_dir, paste0(panel_id, "_cluster_annotation.csv")), row.names = FALSE)
  utils::write.csv(cbind(cluster = rownames(med), as.data.frame(med)),
                   file.path(out_dir, paste0(panel_id, "_cluster_median_markers.csv")), row.names = FALSE)
  utils::write.csv(freq_df, file.path(out_dir, paste0(panel_id, "_cluster_frequency_by_sample.csv")), row.names = FALSE)
  utils::write.csv(stats_cl, file.path(out_dir, paste0(panel_id, "_cluster_T6_vs_T_stats.csv")), row.names = FALSE)
  utils::write.csv(stats_lin, file.path(out_dir, paste0(panel_id, "_lineage_T6_vs_T_stats.csv")), row.names = FALSE)
  if (!is.null(dat$map)) {
    utils::write.csv(dat$map, file.path(out_dir, paste0(panel_id, "_channel_map.csv")), row.names = FALSE)
  }
  embed_out <- cells[, c("sample", "group", "cluster", "lineage", "UMAP1", "UMAP2", "tSNE1", "tSNE2")]
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
    log_msg("Put T/T6 P1-P3 unmixed files in E:/R/flow J, or set FLOW_DEMO=1 to export demo plots")
    use_demo <- TRUE
    log_msg("Auto-fallback to DEMO so the script can still export figure templates")
  }
} else {
  log_msg("Found unmixed files:\n", paste(file_tab$file, collapse = "\n"))
  if (any(file_tab$group == "T" & grepl("^T6", file_tab$file))) {
    stop("Filename parser classified a T6 file as T; refusing to continue")
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

sum_path <- file.path(result_dir, "T6_vs_T_lineage_stats_all_panels.csv")
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

log_msg("Done. Open results_flow/P1, P2, P3 for UMAP/tSNE PDF+PNG.")
invisible(summaries)

}
