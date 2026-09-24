################################################################################
# 单独跑：下载后放到 E:/R/BRCA，RStudio 打开，从第一行 Source 本文件
#
# 不要再运行旧行：
#   ann[, distant_M := factor(meta_M, ...)]
#   ann[meta_M == "M1" | stage_simplified == "Stage IV", ...]
# 那些写法会在临床表没有 meta_M 时直接报错。
#
# 你的 TCGA-BRCA.clinical.tsv 经常只有人口学列（gender/race/primary_site），
# 没有 AJCC M/N/分期。本脚本会：
#   1) 扫描目录里所有 tsv/csv，找 GDC_phenotype / clinicalMatrix / 含 stage 的表
#   2) 没有的话尝试下载 Xena GDC phenotype
#   3) 再用条形码第 14–15 位（06/07 = 转移组织）和生存 PFI/DFI 补分组
# 转移（临床）= M1 或 Stage IV 或条形码 06/07
# 未转移（临床）= M0 或 Stage I–III（且不是 06/07）
# 若临床分期仍全空，主图改用 PFI 进展 vs 未进展
################################################################################

library(data.table)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)

if (dir.exists("E:/R/BRCA")) setwd("E:/R/BRCA")
if (!file.exists("TCGA-BRCA.clinical.tsv") && !file.exists("TCGA-BRCA.survival.tsv")) {
  stop("请把本文件放在 E:/R/BRCA（需有 clinical 或 survival tsv）")
}

# ------------------------------------------------------------------------------
# 工具
# ------------------------------------------------------------------------------
normalize_barcode <- function(x) {
  x <- toupper(gsub("\\.", "-", as.character(x)))
  x <- sub("A$", "", x)
  ifelse(nchar(x) >= 15, substr(x, 1, 15), x)
}
patient_id <- function(x) {
  x <- normalize_barcode(x)
  ifelse(nchar(x) >= 12, substr(x, 1, 12), x)
}
# TCGA 条形码第 14–15 位：01 原发，06/07 转移组织
sample_type_code <- function(x) {
  x <- normalize_barcode(x)
  ifelse(nchar(x) >= 15, substr(x, 14, 15), NA_character_)
}
first_present <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) NA_character_ else hit[1]
}
simplify_stage <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("IV|STAGE.?4", x)] <- "Stage IV"
  out[is.na(out) & grepl("III|STAGE.?3", x)] <- "Stage III"
  out[is.na(out) & grepl("II|STAGE.?2", x)] <- "Stage II"
  out[is.na(out) & grepl("I|STAGE.?1", x)] <- "Stage I"
  out[grepl("\\bX\\b|NOT AVAILABLE|UNKNOWN|NOT REPORTED|NOT APPLICABLE", x)] <- NA_character_
  out
}
classify_m <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("M1[ABC]?|\\bM1\\b", x)] <- "M1"
  out[is.na(out) & grepl("\\bM0\\b|CM0", x)] <- "M0"
  out[grepl("MX|NOT AVAILABLE|UNKNOWN|NOT REPORTED", x)] <- NA_character_
  out
}
classify_n <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("N[1-3]", x)] <- "Nplus"
  out[is.na(out) & grepl("\\bN0\\b", x)] <- "N0"
  out[grepl("NX|NOT AVAILABLE|UNKNOWN|NOT REPORTED", x)] <- NA_character_
  out
}

pick_clin_col <- function(dt, patterns, value_pat = NULL) {
  nms <- names(dt)
  skip <- grepl(
    "^(sample|patient|barcode|submitter|case_id|id|project|age_|days_|year_|uuid)",
    nms, ignore.case = TRUE
  )
  nms2 <- nms[!skip]
  for (p in patterns) {
    hit <- grep(p, nms2, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  if (is.null(value_pat)) return(NA_character_)
  for (nm in nms2) {
    v <- dt[[nm]]
    if (!(is.character(v) || is.factor(v))) next
    u <- unique(toupper(trimws(as.character(v))))
    u <- u[!(u %in% c("", "NA", "NAN", "NULL", "--", ".", "NOT AVAILABLE", "UNKNOWN"))]
    if (length(u) < 2L || length(u) > 80L) next
    if (mean(grepl(value_pat, u)) >= 0.25) return(nm)
  }
  NA_character_
}

is_huge_expr_file <- function(f) {
  grepl("fpkm|htseq|star_|tpm|count_tracking|read_group|probemap|protein",
        basename(f), ignore.case = TRUE)
}

scan_local_clin_files <- function() {
  files <- unique(c(
    list.files(".", pattern = "\\.(tsv|txt|csv)(\\.(gz|bz2))?$",
               full.names = TRUE, ignore.case = TRUE),
    list.files(".", pattern = "phenotype|clinicalMatrix|clin\\.merged|clinical",
               full.names = TRUE, ignore.case = TRUE)
  ))
  files <- files[file.exists(files) & !is_huge_expr_file(files)]
  files <- files[!grepl("results_GO_individual", files, ignore.case = TRUE)]
  if (length(files) == 0) return(data.table())

  rows <- lapply(files, function(f) {
    hdr <- tryCatch(names(fread(f, nrows = 0, header = TRUE)), error = function(e) character())
    if (length(hdr) == 0) return(NULL)
    n_hit <- sum(grepl(
      "ajcc|pathologic_[tmn]|pathologic_stage|clinical_stage|clinical_[tmn]|tumor_stage|metastas|_pm$|_pn$",
      hdr, ignore.case = TRUE
    ))
    data.table(file = f, n_hit = n_hit, n_col = length(hdr),
               hits = paste(grep("stage|ajcc|pathologic|metastas|clinical_[tmn]",
                                 hdr, ignore.case = TRUE, value = TRUE)[1:12],
                            collapse = ", "))
  })
  out <- rbindlist(Filter(Negate(is.null), rows), fill = TRUE)
  if (nrow(out) > 0) setorder(out, -n_hit, -n_col)
  out
}

try_download_phenotype <- function() {
  dests <- c("TCGA-BRCA.GDC_phenotype.tsv.gz", "TCGA-BRCA.GDC_phenotype.tsv")
  if (any(file.exists(dests))) return(dests[file.exists(dests)][1])
  urls <- c(
    "https://gdc.xenahubs.net/download/TCGA-BRCA/Xena_Matrices/TCGA-BRCA.GDC_phenotype.tsv.gz",
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.GDC_phenotype.tsv.gz",
    "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix.gz"
  )
  dest <- "TCGA-BRCA.GDC_phenotype.tsv.gz"
  for (u in urls) {
    message("尝试下载含分期的临床表：", u)
    ok <- tryCatch({
      utils::download.file(u, destfile = dest, mode = "wb", quiet = TRUE)
      file.exists(dest) && file.info(dest)$size > 10000
    }, error = function(e) FALSE)
    if (isTRUE(ok)) {
      message("已下载：", dest)
      return(dest)
    }
    if (file.exists(dest)) unlink(dest)
  }
  message("自动下载失败。请手动把 TCGA-BRCA.GDC_phenotype.tsv 放到 E:/R/BRCA")
  NA_character_
}

read_id_table <- function(path) {
  dt <- fread(path)
  idc <- first_present(names(dt), c(
    "sampleID", "sample", "sample_id", "bcr_sample_barcode",
    "bcr_patient_barcode", "submitter_id.samples", "submitter_id",
    names(dt)[1]
  ))
  dt[, sample_std := normalize_barcode(dt[[idc]])]
  dt[, patient_std := patient_id(dt[[idc]])]
  dt <- dt[!(sample_std %in% c("", "NA", "NAN")) & !is.na(sample_std)]
  dt
}

fill_columns_from <- function(ann, extra, by_x, by_y) {
  if (!by_x %in% names(ann) || !by_y %in% names(extra)) return(ann)
  extra <- extra[!is.na(extra[[by_y]]) & extra[[by_y]] != ""]
  extra <- extra[!duplicated(extra[[by_y]])]
  idx <- match(ann[[by_x]], extra[[by_y]])
  add_nms <- setdiff(names(extra), c(by_x, by_y, "sample_std", "patient_std"))
  for (nm in add_nms) {
    src <- extra[[nm]][idx]
    if (nm %in% names(ann)) {
      cur <- ann[[nm]]
      miss <- is.na(cur) | as.character(cur) %in% c("", "NA", "NAN")
      cur[miss] <- src[miss]
      ann[, (nm) := cur]
    } else {
      ann[, (nm) := src]
    }
  }
  ann
}

merge_clin_onto <- function(ann, clin, keep_extra = TRUE) {
  if (is.null(clin) || nrow(clin) == 0) return(ann)
  useful <- grep(
    "stage|ajcc|pathologic|clinical_[tmn]|metastas|_pm|_pn|_pt|tumor_stage|figo|OS$|OS\\.time|PFI|DFI|DSS",
    names(clin), ignore.case = TRUE, value = TRUE
  )
  keys <- intersect(c("sample_std", "patient_std"), names(clin))
  if (length(useful) == 0) return(ann)
  extra <- clin[, unique(c(keys, useful)), with = FALSE]
  if ("sample_std" %in% names(extra)) {
    ann <- fill_columns_from(ann, extra, "sample", "sample_std")
  }
  if ("patient_std" %in% names(extra)) {
    if (!"patient_std" %in% names(ann)) ann[, patient_std := patient_id(sample)]
    ann <- fill_columns_from(ann, extra, "patient_std", "patient_std")
  }
  ann
}

# 只用向量写列，禁止 i 里写 meta_M / stage_simplified
add_met_groups <- function(dt, sample_col = NULL) {
  dt <- as.data.table(dt)
  n <- nrow(dt)
  nms <- names(dt)
  message("======= 临床列名（全部）=======")
  message(paste(nms, collapse = " | "))
  message("======= 名字像分期/M/N 的列 =======")
  maybe <- grep("stage|ajcc|pathologic|metastas|tnm|_pm|_pn|tumor_stage", nms, ignore.case = TRUE, value = TRUE)
  message(if (length(maybe) == 0) "  （没有）" else paste("  ", maybe, collapse = "\n"))

  st <- pick_clin_col(dt, c(
    "ajcc_pathologic_tumor_stage", "ajcc_pathologic_stage", "pathologic_stage",
    "tumor_stage", "clinical_stage", "figo_stage", "ajcc_staging", "stage"
  ), value_pat = "STAGE\\s*[IVX1-4]|\\bI{1,3}V?\\b")
  mc <- pick_clin_col(dt, c(
    "ajcc_pathologic_m", "ajcc_metastasis_pathologic_pm", "pathologic_m",
    "clinical_m", "metastasis_pathologic", "ajcc_clinical_m"
  ), value_pat = "\\bM[0-1X]")
  nc <- pick_clin_col(dt, c(
    "ajcc_pathologic_n", "ajcc_nodes_pathologic_pn", "pathologic_n",
    "clinical_n", "nodes_pathologic", "ajcc_clinical_n"
  ), value_pat = "\\bN[0-3X]")
  message("选用列：stage=", st, "  M=", mc, "  N=", nc)

  st_vec <- if (!is.na(st) && st %in% names(dt)) simplify_stage(dt[[st]]) else rep(NA_character_, n)
  m_vec  <- if (!is.na(mc) && mc %in% names(dt)) classify_m(dt[[mc]]) else rep(NA_character_, n)
  n_vec  <- if (!is.na(nc) && nc %in% names(dt)) classify_n(dt[[nc]]) else rep(NA_character_, n)

  if (is.null(sample_col)) {
    sample_col <- first_present(names(dt), c("sample", "sample_std", "sampleID", "sample_id"))
  }
  type_vec <- if (!is.na(sample_col) && sample_col %in% names(dt)) {
    sample_type_code(dt[[sample_col]])
  } else {
    rep(NA_character_, n)
  }
  tissue_met <- !is.na(type_vec) & type_vec %in% c("06", "07")

  any_vec <- rep(NA_character_, n)
  any_vec[!is.na(m_vec) & m_vec == "M1"] <- "转移"
  any_vec[!is.na(st_vec) & st_vec == "Stage IV"] <- "转移"
  any_vec[tissue_met] <- "转移"
  any_vec[is.na(any_vec) & !is.na(m_vec) & m_vec == "M0"] <- "未转移"
  any_vec[is.na(any_vec) & !is.na(st_vec) & st_vec %in% c("Stage I", "Stage II", "Stage III")] <- "未转移"

  dt[, stage_simplified := st_vec]
  dt[, meta_M := m_vec]
  dt[, meta_N := n_vec]
  dt[, distant_M := factor(m_vec, levels = c("M0", "M1"))]
  dt[, node_N := factor(n_vec, levels = c("N0", "Nplus"))]
  dt[, sample_type := type_vec]
  dt[, any_met := factor(any_vec, levels = c("未转移", "转移"))]
  dt
}

add_surv_flags <- function(dt) {
  n <- nrow(dt)
  if ("PFI" %in% names(dt)) {
    pfi <- suppressWarnings(as.numeric(dt[["PFI"]]))
    prog <- rep(NA_character_, n)
    prog[!is.na(pfi) & pfi == 0] <- "未进展"
    prog[!is.na(pfi) & pfi == 1] <- "进展"
    dt[, progressed := factor(prog, levels = c("未进展", "进展"))]
  }
  if ("DFI" %in% names(dt)) {
    dfi <- suppressWarnings(as.numeric(dt[["DFI"]]))
    dfi_g <- rep(NA_character_, n)
    dfi_g[!is.na(dfi) & dfi == 0] <- "无病"
    dfi_g[!is.na(dfi) & dfi == 1] <- "复发"
    dt[, dfi_event := factor(dfi_g, levels = c("无病", "复发"))]
  }
  dt
}

# ------------------------------------------------------------------------------
# 通路分数
# ------------------------------------------------------------------------------
pick_res_dir <- function() {
  hits <- list.files(".", pattern = "^01_pathway_scores_each_GO\\.csv$",
                     recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(hits) > 0) {
    pref <- hits[grepl("results_GO_individual-1", hits, fixed = TRUE)]
    return(dirname(if (length(pref) > 0) pref[1] else hits[1]))
  }
  cands <- unique(c("results_GO_individual-1", "results_GO_individual",
                    list.files(".", pattern = "^results_GO_individual", include.dirs = TRUE)))
  cands <- cands[dir.exists(cands)]
  if (length(cands) == 0) return("results_GO_individual-1")
  has_per <- vapply(cands, function(d) {
    length(list.files(file.path(d, "per_GO"), pattern = "^pathway_score\\.csv$", recursive = TRUE)) > 0
  }, logical(1))
  if (any(has_per)) return(cands[has_per][1])
  cands[1]
}

load_or_build_score_mat <- function(res_dir) {
  candidates <- c(
    file.path(res_dir, "01_pathway_scores_each_GO.csv"),
    file.path(res_dir, "01_pathway_scores_each_GO.CSV"),
    "01_pathway_scores_each_GO.csv"
  )
  hit <- candidates[file.exists(candidates)][1]
  if (!is.na(hit) && length(hit) == 1) {
    message("读取汇总分数表：", normalizePath(hit, winslash = "/"))
    sm <- fread(hit)
    go_cols <- setdiff(names(sm), "sample")
    if (length(go_cols) == 0) stop("分数表没有 GO 列：", hit)
    mat <- as.matrix(sm[, go_cols, with = FALSE])
    storage.mode(mat) <- "double"
    rownames(mat) <- as.character(sm$sample)
    colnames(mat) <- go_cols
    return(mat)
  }

  per_dir <- file.path(res_dir, "per_GO")
  score_files <- if (dir.exists(per_dir)) {
    list.files(per_dir, pattern = "^pathway_score\\.csv$", recursive = TRUE, full.names = TRUE)
  } else character()
  if (length(score_files) == 0) {
    message("当前工作目录：", getwd())
    message("结果目录：", res_dir, " 存在=", dir.exists(res_dir))
    if (dir.exists(res_dir)) message(paste("  ", list.files(res_dir), collapse = "\n"))
    stop("找不到通路分数。请先跑完 GO_pathway_individual_analysis.R")
  }

  message("未找到汇总表，改为从 ", length(score_files), " 个 per_GO/pathway_score.csv 拼接")
  lst <- lapply(score_files, function(f) {
    dt <- fread(f)
    if (!all(c("sample", "pathway_score") %in% names(dt))) return(NULL)
    folder <- basename(dirname(f))
    go_id <- sub("^GO_", "GO:", folder)
    go_id <- sub("_.*$", "", go_id)
    if (!grepl("^GO:", go_id)) go_id <- folder
    v <- as.numeric(dt$pathway_score)
    names(v) <- as.character(dt$sample)
    attr(v, "GO") <- go_id
    v
  })
  lst <- Filter(Negate(is.null), lst)
  if (length(lst) == 0) stop("per_GO 里的 pathway_score.csv 格式不对")
  names(lst) <- vapply(lst, function(x) attr(x, "GO"), character(1))
  common <- Reduce(intersect, lapply(lst, names))
  if (length(common) < 10) common <- unique(unlist(lapply(lst, names), use.names = FALSE))
  mat <- do.call(cbind, lapply(lst, function(x) {
    out <- rep(NA_real_, length(common))
    names(out) <- common
    out[intersect(names(x), common)] <- x[intersect(names(x), common)]
    out
  }))
  colnames(mat) <- names(lst)
  rownames(mat) <- common
  out_csv <- file.path(res_dir, "01_pathway_scores_each_GO.csv")
  fwrite(data.table(sample = rownames(mat), as.data.table(mat)), out_csv)
  message("已写出拼接后的分数表：", out_csv)
  mat
}

go_name_map <- c(
  "GO:0023041" = "neuronal signal transduction",
  "GO:1904457" = "positive regulation of neuronal action potential",
  "GO:1904340" = "positive regulation of dopaminergic neuron differentiation",
  "GO:2001224" = "positive regulation of neuron migration",
  "GO:2001222" = "regulation of neuron migration",
  "GO:0019227" = "neuronal action potential propagation",
  "GO:0019228" = "neuronal action potential",
  "GO:1902847" = "regulation of neuronal signal transduction",
  "GO:0031102" = "neuron projection regeneration",
  "GO:0097492" = "sympathetic neuron axon guidance",
  "GO:0097491" = "sympathetic neuron projection guidance",
  "GO:0097374" = "sensory neuron axon guidance",
  "GO:0007158" = "neuron cell-cell adhesion",
  "GO:1902667" = "regulation of axon guidance",
  "GO:0031103" = "axon regeneration",
  "GO:0007411" = "axon guidance",
  "GO:0007409" = "axonogenesis",
  "GO:0036518" = "chemorepulsion of dopaminergic neuron axon"
)

# ------------------------------------------------------------------------------
# 读分数（当前会话已有 score_mat 就复用）
# ------------------------------------------------------------------------------
res_dir <- pick_res_dir()
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
message("使用结果目录：", normalizePath(res_dir, winslash = "/", mustWork = FALSE))

if (!(exists("score_mat", inherits = FALSE) && is.matrix(score_mat) && nrow(score_mat) > 5)) {
  score_mat <- load_or_build_score_mat(res_dir)
}
message("通路分数：", nrow(score_mat), " 样本 x ", ncol(score_mat), " 个 GO")

# ------------------------------------------------------------------------------
# 临床：扫描全部表，必要时下载 phenotype
# ------------------------------------------------------------------------------
found <- scan_local_clin_files()
if (nrow(found) > 0) {
  message("目录里可能含临床/表型的文件：")
  print(found[, .(file = basename(file), n_hit, n_col)])
}

preferred <- c(
  "TCGA-BRCA.GDC_phenotype.tsv", "TCGA-BRCA.GDC_phenotype.tsv.gz",
  "BRCA_clinicalMatrix", "BRCA_clinicalMatrix.gz",
  "TCGA-BRCA.clinicalMatrix.tsv", "TCGA-BRCA.clinical.tsv"
)
local_hits <- if (nrow(found) > 0) found$file[found$n_hit > 0] else character()
clin_files <- unique(c(preferred[file.exists(preferred)], local_hits))
if (length(clin_files) == 0 || !any(grepl("phenotype|clinicalMatrix|stage", clin_files, ignore.case = TRUE))) {
  dl <- try_download_phenotype()
  if (!is.na(dl) && file.exists(dl)) clin_files <- unique(c(dl, clin_files))
}
if (file.exists("TCGA-BRCA.clinical.tsv")) {
  clin_files <- unique(c(clin_files, "TCGA-BRCA.clinical.tsv"))
}

clin_list <- lapply(unique(clin_files), function(f) {
  message("读取临床/表型：", f)
  tryCatch(read_id_table(f), error = function(e) {
    message("  跳过（读失败）：", conditionMessage(e))
    NULL
  })
})
clin_list <- Filter(Negate(is.null), clin_list)

surv <- NULL
surv_cands <- c("TCGA-BRCA.survival.tsv", "TCGA-BRCA.survival.tsv.gz")
surv_hit <- surv_cands[file.exists(surv_cands)][1]
if (!is.na(surv_hit)) {
  surv <- read_id_table(surv_hit)
  message("读取生存表：", surv_hit, "  列=", paste(names(surv), collapse = ", "))
}

# ------------------------------------------------------------------------------
# 建 ann：已有 ann 也重新对齐分数样本，再并上所有临床/生存
# ------------------------------------------------------------------------------
ann <- data.table(sample = normalize_barcode(rownames(score_mat)))
ann[, patient_std := patient_id(sample)]
for (cl in clin_list) ann <- merge_clin_onto(ann, cl)
if (!is.null(surv)) {
  keep <- intersect(
    c("sample_std", "patient_std", "OS", "OS.time", "DSS", "DSS.time",
      "PFI", "PFI.time", "DFI", "DFI.time"),
    names(surv)
  )
  ann <- merge_clin_onto(ann, surv[, keep, with = FALSE])
}

ann <- add_met_groups(ann, sample_col = "sample")
ann <- add_surv_flags(ann)

n_met <- sum(as.character(ann$any_met) == "转移", na.rm = TRUE)
n_non <- sum(as.character(ann$any_met) == "未转移", na.rm = TRUE)
n_m1  <- sum(as.character(ann$distant_M) == "M1", na.rm = TRUE)
n_pfi1 <- if ("progressed" %in% names(ann)) sum(as.character(ann$progressed) == "进展", na.rm = TRUE) else 0L
n_pfi0 <- if ("progressed" %in% names(ann)) sum(as.character(ann$progressed) == "未进展", na.rm = TRUE) else 0L
n_06 <- sum(ann$sample_type %in% c("06", "07"), na.rm = TRUE)

message("分组人数：转移=", n_met, "  未转移=", n_non,
        "  M1=", n_m1, "  条形码06/07=", n_06,
        "  PFI进展=", n_pfi1, "  PFI未进展=", n_pfi0)

# 主图用哪一组：临床转移优先，否则 PFI
if (n_met >= 2 && n_non >= 2) {
  plot_group <- as.character(ann$any_met)
  plot_pos <- "转移"
  plot_neg <- "未转移"
  plot_title <- "神经相关 GO 在转移 vs 未转移中的表达"
  plot_sub <- "转移 = M1 或 Stage IV 或条形码 06/07；每个 GO 单独打分"
  plot_delta <- "Δ median (转移 − 未转移)"
  plot_fill <- c("未转移" = "#4DBBD5", "转移" = "#E64B35")
  plot_tag <- "any_met"
} else if (n_pfi1 >= 2 && n_pfi0 >= 2) {
  message("临床表没有可用的 M/N/分期，主图改用 PFI 进展 vs 未进展")
  plot_group <- as.character(ann$progressed)
  plot_pos <- "进展"
  plot_neg <- "未进展"
  plot_title <- "神经相关 GO 在 PFI 进展 vs 未进展中的表达"
  plot_sub <- "当前 clinical.tsv 无 AJCC M/N/分期；用生存 PFI 代替转移相关结局"
  plot_delta <- "Δ median (进展 − 未进展)"
  plot_fill <- c("未进展" = "#4DBBD5", "进展" = "#E64B35")
  plot_tag <- "PFI_progressed"
  # 让后续箱线复用 any_met 列：把 PFI 填进去，避免旧代码再找 meta_M
  ann[, any_met := factor(plot_group, levels = c(plot_neg, plot_pos))]
} else {
  sink(file.path(res_dir, "09_metastasis_grouping_log.txt"))
  cat("无法构建转移或 PFI 分组\n列名：\n")
  cat(paste(names(ann), collapse = "\n"))
  sink()
  stop(
    "无法分组：clinical.tsv 没有分期/M/N，也没有可用 PFI。",
    "请把 Xena 的 TCGA-BRCA.GDC_phenotype.tsv 放到 E:/R/BRCA 后重新 Source。",
    "列名已写入 ", file.path(res_dir, "09_metastasis_grouping_log.txt")
  )
}

fwrite(
  ann[, intersect(c(
    "sample", "any_met", "distant_M", "node_N", "progressed", "dfi_event",
    "stage_simplified", "meta_M", "meta_N", "sample_type",
    "OS", "OS.time", "PFI", "PFI.time", "DFI", "DFI.time"
  ), names(ann)), with = FALSE],
  file.path(res_dir, "09_sample_metastasis_prognosis.csv")
)
writeLines(
  c(
    paste0("res_dir=", res_dir),
    paste0("n_met=", n_met), paste0("n_non=", n_non), paste0("n_m1=", n_m1),
    paste0("n_06=", n_06), paste0("n_pfi1=", n_pfi1), paste0("n_pfi0=", n_pfi0),
    paste0("plot_tag=", plot_tag),
    paste0("clin_files=", paste(clin_files, collapse = " ; ")),
    paste0("columns=", paste(names(ann), collapse = " | "))
  ),
  file.path(res_dir, "09_metastasis_grouping_log.txt")
)

# ------------------------------------------------------------------------------
# 比较与作图
# ------------------------------------------------------------------------------
compare_one <- function(go_score_vec, group, pos, neg, grouping) {
  df <- data.frame(
    pathway_score = as.numeric(go_score_vec),
    group = as.character(group),
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$pathway_score) & df$group %in% c(pos, neg), ]
  if (sum(df$group == pos) < 2 || sum(df$group == neg) < 2) return(NULL)
  wt <- suppressWarnings(stats::wilcox.test(pathway_score ~ group, data = df))
  data.table(
    grouping = grouping, pos_level = pos, neg_level = neg,
    n_pos = sum(df$group == pos), n_neg = sum(df$group == neg),
    median_pos = stats::median(df$pathway_score[df$group == pos], na.rm = TRUE),
    median_neg = stats::median(df$pathway_score[df$group == neg], na.rm = TRUE),
    delta_median = stats::median(df$pathway_score[df$group == pos], na.rm = TRUE) -
      stats::median(df$pathway_score[df$group == neg], na.rm = TRUE),
    pvalue = wt$p.value
  )
}

stat_rows <- list()
long_rows <- list()
for (g in colnames(score_mat)) {
  gnm <- if (g %in% names(go_name_map)) unname(go_name_map[g]) else g
  sc <- as.numeric(score_mat[, g])
  names(sc) <- normalize_barcode(rownames(score_mat))
  sc <- sc[ann$sample]
  one <- rbindlist(Filter(Negate(is.null), list(
    compare_one(sc, plot_group, plot_pos, plot_neg, plot_tag),
    compare_one(sc, ann$distant_M, "M1", "M0", "distant_M"),
    compare_one(sc, ann$node_N, "Nplus", "N0", "node_N"),
    if ("progressed" %in% names(ann)) compare_one(sc, ann$progressed, "进展", "未进展", "PFI_progressed") else NULL,
    if ("dfi_event" %in% names(ann)) compare_one(sc, ann$dfi_event, "复发", "无病", "DFI_event") else NULL
  )), fill = TRUE)
  if (nrow(one) > 0) {
    one[, `:=`(GO = g, GO_name = gnm)]
    stat_rows[[g]] <- one
  }
  long_rows[[g]] <- data.table(
    sample = ann$sample, GO = g, GO_name = gnm,
    pathway_score = as.numeric(sc),
    plot_group = factor(plot_group, levels = c(plot_neg, plot_pos))
  )
}

stat_dt <- rbindlist(stat_rows, fill = TRUE)
if (nrow(stat_dt) == 0) stop("分组后没有可比较的通路（每组人数不足）")
stat_dt[, fdr := p.adjust(pvalue, method = "BH"), by = grouping]
fwrite(stat_dt, file.path(res_dir, "09_nerve_GO_score_by_metastasis.csv"))

long_any <- rbindlist(long_rows, fill = TRUE)
long_any <- long_any[!is.na(plot_group) & is.finite(pathway_score)]
long_any[, go_lab := factor(paste(GO, GO_name), levels = unique(paste(GO, GO_name)))]
p_box <- ggplot(long_any, aes(x = plot_group, y = pathway_score, fill = plot_group)) +
  geom_boxplot(outlier.size = 0.4, width = 0.65) +
  stat_compare_means(size = 2.6, label = "p.format") +
  facet_wrap(~ go_lab, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = plot_fill) +
  labs(title = plot_title, subtitle = plot_sub, x = NULL, y = "Pathway score", fill = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(size = 7),
        axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(res_dir, "09_nerve_GO_boxplot_any_met.pdf"), p_box, width = 12, height = 10)
ggsave(file.path(res_dir, "09_nerve_GO_boxplot_any_met.png"), p_box, width = 12, height = 10, dpi = 150)
print(p_box)

any_stat <- stat_dt[grouping == plot_tag]
if (nrow(any_stat) > 0) {
  any_stat[, lab := factor(paste(GO, GO_name), levels = paste(GO, GO_name)[order(delta_median)])]
  p_for <- ggplot(any_stat, aes(x = delta_median, y = lab)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    geom_point(aes(color = pvalue < 0.05, size = -log10(pmax(pvalue, 1e-12)))) +
    scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "#E64B35"), name = "p < 0.05") +
    labs(title = paste0(plot_pos, " − ", plot_neg, "：神经 GO 通路分数差"),
         x = plot_delta, y = NULL, size = expression(-log[10](p))) +
    theme_bw()
  ggsave(file.path(res_dir, "09_nerve_GO_delta_any_met.pdf"), p_for, width = 10, height = 6)
  print(p_for)
}

surv_rows <- list()
if (all(c("OS", "OS.time") %in% names(ann))) {
  pdf(file.path(res_dir, "09_nerve_GO_KM_OS_within_met.pdf"), width = 10, height = 5)
  n_km <- 0L
  for (g in colnames(score_mat)) {
    gnm <- if (g %in% names(go_name_map)) unname(go_name_map[g]) else g
    d0 <- data.frame(
      time = as.numeric(ann[["OS.time"]]),
      event = as.numeric(ann[["OS"]]),
      pathway_score = as.numeric(score_mat[match(ann$sample, normalize_barcode(rownames(score_mat))), g]),
      met = as.character(plot_group),
      stringsAsFactors = FALSE
    )
    d0 <- d0[is.finite(d0$time) & d0$time > 0 & d0$event %in% c(0, 1) &
               is.finite(d0$pathway_score) & !is.na(d0$met), ]
    for (mg in c(plot_neg, plot_pos)) {
      d <- d0[d0$met == mg, ]
      if (nrow(d) < 10 || sum(d$event) < 3) next
      d$group <- factor(
        ifelse(d$pathway_score >= stats::median(d$pathway_score, na.rm = TRUE), "High", "Low"),
        levels = c("Low", "High")
      )
      cox_g <- tryCatch(coxph(Surv(time, event) ~ group, data = d), error = function(e) NULL)
      if (!is.null(cox_g)) {
        s2 <- summary(cox_g)
        surv_rows[[paste(g, mg)]] <- data.table(
          GO = g, GO_name = gnm, met_group = mg, n = nrow(d), events = sum(d$event),
          HR = s2$conf.int[1, 1], HR_low = s2$conf.int[1, 3], HR_high = s2$conf.int[1, 4],
          pvalue = s2$coefficients[1, "Pr(>|z|)"]
        )
      }
      fit <- tryCatch(survfit(Surv(time / 30.44, event) ~ group, data = d), error = function(e) NULL)
      if (!is.null(fit)) {
        print(ggsurvplot(
          fit, data = d, pval = TRUE, risk.table = TRUE,
          legend.labs = c("Low", "High"),
          xlab = "Time (months)", ylab = "Overall survival",
          title = paste0(g, " | ", gnm, "\n", mg, " 亚组 High vs Low"),
          ggtheme = theme_bw()
        ))
        n_km <- n_km + 1L
      }
    }
  }
  invisible(dev.off())
  message("KM 页数：", n_km)
}

if (length(surv_rows) > 0) {
  surv_dt <- rbindlist(surv_rows, fill = TRUE)
  fwrite(surv_dt, file.path(res_dir, "09_nerve_GO_OS_cox_within_met.csv"))
  p_hr <- ggplot(surv_dt, aes(x = HR, y = reorder(paste(GO, met_group), HR))) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
    geom_errorbar(aes(xmin = HR_low, xmax = HR_high), orientation = "y", width = 0.2) +
    geom_point(aes(color = pvalue < 0.05), size = 2.5) +
    scale_x_log10() +
    facet_wrap(~ met_group, scales = "free_y") +
    labs(title = "OS High vs Low（分组亚组内）", x = "Hazard ratio", y = NULL, color = "p < 0.05") +
    theme_bw()
  ggsave(file.path(res_dir, "09_nerve_GO_OS_forest_within_met.pdf"), p_hr, width = 11, height = 7)
  print(p_hr)
}

message("完成。主表：", file.path(res_dir, "09_nerve_GO_score_by_metastasis.csv"))
print(stat_dt[grouping == plot_tag][order(pvalue)])
