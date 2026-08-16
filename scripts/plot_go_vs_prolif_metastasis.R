#!/usr/bin/env Rscript
################################################################################
# 附加分析（不改动 analyze_go_excel_fusc.R）
# 气泡图：5 个神经/轴突 GO 通路  vs  肿瘤增殖、转移通路
#
# Y：GO:2001222, GO:1904340, GO:1902667, GO:0097374, GO:0036518
# X：细胞增殖 / 正调控增殖 / 有丝分裂周期 / 细胞迁移 / 正调控迁移 / EMT / 血管生成
# 每个通路单独 z-score 打分，再做 Spearman；禁止把通路基因合并
#
# 用法（E:/R/FUSC）：
#   1) 可先 Source 完 analyze_go_excel_fusc.R（内存里已有 expr 会直接复用）
#   2) 再 Source 本文件
#   或单独 Source 本文件（会自己读 Excel 临床 + CSV FPKM）
################################################################################

suppressPackageStartupMessages({
  for (p in c("data.table", "readxl", "ggplot2")) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
    library(p, character.only = TRUE)
  }
})

if (!exists("work_dir") || !dir.exists(as.character(work_dir))) {
  args_pm <- commandArgs(trailingOnly = TRUE)
  work_dir <- if (length(args_pm) >= 1) args_pm[[1]] else "E:/R/FUSC"
}
if (dir.exists(work_dir)) {
  setwd(work_dir)
} else {
  message("未找到 ", work_dir, " ，改用：", getwd())
  work_dir <- getwd()
}

this_script_dir <- (function() {
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) return(dirname(normalizePath(ofile)))
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[startsWith(ca, "--file=")])
  if (length(f)) return(dirname(normalizePath(f[[1]])))
  if (file.exists("/workspace/scripts")) return("/workspace/scripts")
  getwd()
})()

out_dir_pm <- if (exists("out_dir") && !exists("args_pm")) {
  file.path(out_dir, "prolif_metastasis")
} else if (exists("args_pm") && length(args_pm) >= 2) {
  args_pm[[2]]
} else {
  "results_GO_excel/prolif_metastasis"
}
dir.create(out_dir_pm, recursive = TRUE, showWarnings = FALSE)

y_go_ids <- c(
  "GO:2001222",  # regulation of neuron migration
  "GO:1904340",  # positive regulation of dopaminergic neuron differentiation
  "GO:1902667",  # regulation of axon guidance
  "GO:0097374",  # sensory neuron axon guidance
  "GO:0036518"   # chemorepulsion of dopaminergic neuron axon
)
y_go_lab <- c(
  "GO:2001222" = "GO:2001222  神经元迁移调控",
  "GO:1904340" = "GO:1904340  多巴胺能神经元分化正调控",
  "GO:1902667" = "GO:1902667  轴突导向调控",
  "GO:0097374" = "GO:0097374  感觉神经元轴突导向",
  "GO:0036518" = "GO:0036518  多巴胺能轴突化学排斥"
)

# 横坐标：肿瘤增殖、转移相关通路（每个单独打分）
x_go_ids <- c(
  "GO:0008283",  # cell proliferation
  "GO:0008284",  # positive regulation of cell proliferation
  "GO:0000278",  # mitotic cell cycle
  "GO:0016477",  # cell migration
  "GO:0030335",  # positive regulation of cell migration
  "GO:0001837",  # epithelial to mesenchymal transition
  "GO:0001525"   # angiogenesis
)
x_go_lab <- c(
  "GO:0008283" = "细胞增殖",
  "GO:0008284" = "增殖正调控",
  "GO:0000278" = "有丝分裂周期",
  "GO:0016477" = "细胞迁移",
  "GO:0030335" = "迁移正调控",
  "GO:0001837" = "EMT",
  "GO:0001525" = "血管生成"
)
x_go_group <- c(
  "GO:0008283" = "增殖",
  "GO:0008284" = "增殖",
  "GO:0000278" = "增殖",
  "GO:0016477" = "转移",
  "GO:0030335" = "转移",
  "GO:0001837" = "转移",
  "GO:0001525" = "转移"
)

# ---------------------------------------------------------------------------
# 若主脚本已读入 expr，则复用；否则按同样规则读两份数据
# ---------------------------------------------------------------------------
ensure_expr <- function() {
  if (exists("expr", inherits = TRUE)) {
    ex <- get("expr", inherits = TRUE)
    if (is.matrix(ex) && nrow(ex) > 10 && ncol(ex) > 5) {
      message("复用内存中的表达矩阵：", nrow(ex), " 基因 x ", ncol(ex), " 样本")
      return(ex)
    }
  }
  message("内存中没有 expr，按主分析相同方式读取 FUSC 数据")
  pick_col <- function(df, candidates) {
    nms <- names(df)
    hit <- candidates[candidates %in% nms]
    if (length(hit)) return(hit[[1]])
    low <- tolower(nms)
    for (cand in candidates) {
      i <- which(low == tolower(cand))
      if (length(i)) return(nms[[i[[1]]]])
    }
    NULL
  }
  normalize_id <- function(x) {
    x <- trimws(as.character(x))
    gsub("\\.", "-", gsub("\\s+", "", x))
  }
  numeric_frac <- function(x) {
    if (is.numeric(x)) return(mean(is.finite(x)))
    mean(is.finite(suppressWarnings(as.numeric(as.character(x)))))
  }
  file_ext_lower <- function(path) tolower(sub("^.*\\.", "", basename(path)))
  find_one <- function(stem, prefer, fuzzy) {
    dirs <- unique(c(getwd(), work_dir, "E:/R/FUSC"))
    cands <- character(0)
    for (d in dirs[dir.exists(dirs)]) {
      hits <- list.files(d, full.names = TRUE)
      hits <- hits[!startsWith(basename(hits), "~$")]
      bn <- basename(hits)
      sans <- tolower(sub("\\.(xlsx|xls|csv|tsv|txt)$", "", bn, ignore.case = TRUE))
      keep <- grepl("\\.(xlsx|xls|csv|tsv|txt)$", bn, ignore.case = TRUE) &
        (sans == tolower(stem) | startsWith(sans, tolower(stem)) | grepl(fuzzy, bn, ignore.case = TRUE))
      cands <- c(cands, hits[keep])
    }
    cands <- unique(cands)
    if (!length(cands)) stop("找不到 ", stem, "。请 setwd(\"E:/R/FUSC\")")
    ext <- file_ext_lower(cands)
    cands[order(match(ext, prefer, nomatch = 99L))][[1]]
  }
  clin_path <- find_one("OEP00000155_样本_Human_1786804447300", c("xlsx", "xls", "csv"), "OEP00000155")
  fpkm_path <- find_one("FUSCCTNBC_Expression_RNAseqFPKM", c("csv", "tsv", "xlsx"), "FUSCCTNBC.*FPKM")
  message("临床：", clin_path)
  message("FPKM：", fpkm_path)
  ext <- file_ext_lower(fpkm_path)
  if (ext %in% c("xlsx", "xls")) {
    fpkm <- as.data.table(read_excel(fpkm_path, .name_repair = "minimal"))
  } else {
    fpkm <- fread(fpkm_path, showProgress = FALSE)
  }
  is_num <- vapply(fpkm, function(x) numeric_frac(x) >= 0.6, logical(1))
  ann <- which(!is_num)
  num <- which(is_num)
  gene_col <- pick_col(fpkm, c("gene", "Gene", "symbol", "Symbol", "gene_name", names(fpkm)[ann]))
  if (is.null(gene_col)) gene_col <- names(fpkm)[1]
  mat <- as.matrix(data.frame(lapply(fpkm[, num, with = FALSE], function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  })))
  rownames(mat) <- make.unique(sub("\\..*$", "", as.character(fpkm[[gene_col]])))
  colnames(mat) <- normalize_id(names(fpkm)[num])
  mat <- mat[!duplicated(rownames(mat)), !duplicated(colnames(mat)), drop = FALSE]
  mx <- suppressWarnings(max(mat, na.rm = TRUE))
  if (is.finite(mx) && mx > 50) mat <- log2(mat + 1)
  mat
}

zscore_one <- function(expr_mat, genes) {
  g <- intersect(unique(genes), rownames(expr_mat))
  if (!length(g)) return(NULL)
  sub <- expr_mat[g, , drop = FALSE]
  z <- t(scale(t(sub)))
  z[!is.finite(z)] <- NA
  sc <- colMeans(z, na.rm = TRUE)
  sc[!is.finite(sc)] <- NA
  attr(sc, "n_genes") <- length(g)
  attr(sc, "genes") <- g
  sc
}

load_pm_genesets <- function(wanted) {
  dirs <- unique(c(
    this_script_dir,
    file.path(this_script_dir, "..", "data"),
    work_dir, getwd(), "data", "scripts"
  ))
  hits <- unlist(lapply(dirs, function(d) {
    if (!dir.exists(d)) return(character(0))
    list.files(d, pattern = "go_genesets_prolif_meta\\.(tsv|csv|txt)$",
               full.names = TRUE, ignore.case = TRUE)
  }), use.names = FALSE)
  extra <- unlist(lapply(dirs, function(d) {
    if (!dir.exists(d)) return(character(0))
    list.files(d, pattern = "^go_genesets\\.(tsv|csv)$", full.names = TRUE, ignore.case = TRUE)
  }), use.names = FALSE)
  paths <- unique(c(hits, extra))
  if (!length(paths)) stop("找不到 go_genesets_prolif_meta.tsv，请与本脚本一起放到 E:/R/FUSC")
  raw <- rbindlist(lapply(paths, function(p) {
    dt <- fread(p)
    go_col <- names(dt)[tolower(names(dt)) %in% c("go_id", "go", "term")][1]
    sy_col <- names(dt)[tolower(names(dt)) %in% c("symbol", "gene", "gene_name")][1]
    if (is.na(go_col) || is.na(sy_col)) return(NULL)
    data.table(GO_ID = as.character(dt[[go_col]]), symbol = as.character(dt[[sy_col]]))
  }), fill = TRUE)
  out <- lapply(wanted, function(g) unique(na.omit(raw$symbol[raw$GO_ID == g])))
  names(out) <- wanted
  out
}

spearman_two <- function(a, b) {
  common <- intersect(names(a), names(b))
  x <- as.numeric(a[common])
  y <- as.numeric(b[common])
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(list(rho = NA_real_, p = NA_real_, n = sum(ok)))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

# ---------------------------------------------------------------------------
# 打分 + 相关 + 气泡图
# ---------------------------------------------------------------------------
expr_use <- ensure_expr()
wanted <- c(y_go_ids, x_go_ids)
genesets <- load_pm_genesets(wanted)

score_list <- list()
size_rows <- list()
for (gid in wanted) {
  sc <- zscore_one(expr_use, genesets[[gid]])
  n_ann <- length(genesets[[gid]])
  n_use <- if (is.null(sc)) 0L else attr(sc, "n_genes")
  size_rows[[gid]] <- data.table(GO_ID = gid, n_annotated = n_ann, n_in_expr = n_use)
  message(gid, ": 注释 ", n_ann, " / 表达矩阵 ", n_use)
  if (!is.null(sc)) score_list[[gid]] <- sc
}
fwrite(rbindlist(size_rows), file.path(out_dir_pm, "00_geneset_sizes.csv"))

y_ok <- intersect(y_go_ids, names(score_list))
x_ok <- intersect(x_go_ids, names(score_list))
if (!length(y_ok) || !length(x_ok)) {
  stop("没有足够的通路分数。请确认 go_genesets_prolif_meta.tsv 与 FPKM 基因名为 SYMBOL。")
}

corr_rows <- list()
for (yg in y_ok) {
  for (xg in x_ok) {
    st <- spearman_two(score_list[[yg]], score_list[[xg]])
    corr_rows[[length(corr_rows) + 1L]] <- data.table(
      GO_Y = yg,
      GO_Y_lab = unname(y_go_lab[yg]),
      GO_X = xg,
      GO_X_lab = unname(x_go_lab[xg]),
      axis_group = unname(x_go_group[xg]),
      rho = st$rho,
      pvalue = st$p,
      n = st$n
    )
  }
}
corr_dt <- rbindlist(corr_rows)
corr_dt[, fdr := p.adjust(pvalue, method = "BH")]
fwrite(corr_dt, file.path(out_dir_pm, "01_fiveGO_vs_prolif_metastasis_spearman.csv"))

plot_dt <- corr_dt[is.finite(rho) & is.finite(pvalue)]
plot_dt[, neglogp := pmin(10, -log10(pmax(pvalue, 1e-12)))]
plot_dt[, GO_Y_lab := factor(GO_Y_lab, levels = rev(unname(y_go_lab[y_ok])))]
plot_dt[, GO_X_lab := factor(GO_X_lab, levels = unname(x_go_lab[x_ok]))]
lim <- max(0.3, as.numeric(stats::quantile(abs(plot_dt$rho), 0.98, na.rm = TRUE)), na.rm = TRUE)

p <- ggplot(plot_dt, aes(x = GO_X_lab, y = GO_Y_lab, size = neglogp, color = rho)) +
  geom_point(alpha = 0.92) +
  scale_size_continuous(range = c(3, 14), name = expression(-log[10](p))) +
  scale_color_gradient2(
    low = "#3C5488", mid = "grey90", high = "#E64B35",
    midpoint = 0, limits = c(-lim, lim), name = "Spearman ρ"
  ) +
  facet_grid(~ axis_group, scales = "free_x", space = "free_x") +
  labs(
    title = "五个神经/轴突 GO 通路与肿瘤增殖、转移通路的相关性",
    subtitle = "每个通路单独 z-score 打分后做 Spearman（未合并基因集）",
    caption = "横坐标：肿瘤增殖 / 转移通路活性；纵坐标：指定的五个 GO 通路。红=正相关，蓝=负相关；点越大 p 越小。",
    x = "肿瘤增殖、转移信号通路",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.caption = element_text(hjust = 0, size = 8, color = "grey30")
  )

ht <- max(5.5, 0.55 * length(y_ok) + 3.2)
ggsave(file.path(out_dir_pm, "02_fiveGO_vs_prolif_metastasis_bubble.pdf"), p, width = 11, height = ht)
ggsave(file.path(out_dir_pm, "02_fiveGO_vs_prolif_metastasis_bubble.png"), p, width = 11, height = ht, dpi = 150)

message("已保存气泡图：", normalizePath(file.path(out_dir_pm, "02_fiveGO_vs_prolif_metastasis_bubble.pdf"), mustWork = FALSE))
message("相关表：", file.path(out_dir_pm, "01_fiveGO_vs_prolif_metastasis_spearman.csv"))
