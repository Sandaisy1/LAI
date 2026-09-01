# JY：免疫亚群降维（E:/R/fuction of cell-ljy，比较 JY-NNK / JY-EVNK）
# run: Rscript tests/test_jy_flow.R
Sys.setenv(FLOW_FUNCTIONS_ONLY = "1", JY_FUNCTIONS_ONLY = "1")
Sys.setenv(JY_FLOW_DIR = tempfile("jy_flow_"))
dir.create(Sys.getenv("JY_FLOW_DIR"), recursive = TRUE, showWarnings = FALSE)
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/test_jy_flow.R"
root <- dirname(dirname(normalizePath(this_file)))
source(file.path(root, "JY_Flow_dimred_pipeline.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
expect <- function(got, want, tag) {
  if (!identical(got, want)) {
    fail(sprintf("%s: got %s, want %s", tag, paste(got, collapse = ","), paste(want, collapse = ",")))
  }
}

# 原 flow_panel_map.json / ICI 表不得被 JY 改掉
orig_json <- jsonlite::fromJSON(file.path(root, "flow_panel_map.json"), simplifyVector = FALSE)
orig_p1 <- vapply(orig_json$panels$P1$markers, function(x) x$marker, character(1))
if (length(orig_p1) != 25L) {
  fail(sprintf("original P1 must stay 25 markers, got %s", length(orig_p1)))
}
expect(orig_json$groups[[1]], "EV", "original groups still EV")
expect(orig_json$groups[[2]], "H", "original groups still H")
expect(orig_json$data_dir, "E:/R/fuction of cell", "original data_dir unchanged")

ici_json <- jsonlite::fromJSON(file.path(root, "ICI_flow_panel_map.json"), simplifyVector = FALSE)
ici_p1 <- vapply(ici_json$panels$P1$markers, function(x) x$marker, character(1))
if (!("His" %in% ici_p1)) fail("ICI P1 must still have His")

jy_json <- jsonlite::fromJSON(file.path(root, "JY_flow_panel_map.json"), simplifyVector = FALSE)
jy_p1 <- vapply(jy_json$panels$P1$markers, function(x) x$marker, character(1))
if (length(jy_p1) != 25L) fail(sprintf("JY P1 should copy original 25 markers, got %s", length(jy_p1)))
if ("His" %in% jy_p1) fail("JY P1 must not gain His (this is not the ICI sheet)")
if (!("IgD" %in% vapply(jy_json$panels$P2$markers, function(x) x$marker, character(1)))) {
  fail("JY P2 must keep original IgD (immune-subset panel, not ICI 7-color)")
}
expect(jy_json$groups[[1]], "JY-EVNK", "JY groups JY-EVNK")
expect(jy_json$groups[[2]], "JY-NNK", "JY groups JY-NNK")
expect(jy_json$data_dir, "E:/R/fuction of cell-ljy", "JY data_dir")
expect(jy_json$cohort, "JY", "JY cohort")
expect(jy_json$comparison, "JY_NNK_vs_JY_EVNK", "JY comparison JY-NNK / JY-EVNK")

expect(flow_ctrl_group, "JY-EVNK", "loaded ctrl JY-EVNK")
expect(flow_trt_group, "JY-NNK", "loaded trt JY-NNK")
expect(as.character(flow_group_levels), c("JY-EVNK", "JY-NNK"), "group levels")
if (!isTRUE(flow_trim_bio_extremes)) {
  fail("JY must drop 1 extreme bio-rep like the immune-subset scheme (n=2)")
}

# 文件名：EVNK1-1 / NNK1-2 / JY-EVNK；生物学重复 EVNK-1；先匹配 EVNK
ev <- parse_fcs_filename("EVNK1-1_P1_unmixed.fcs")
if (is.null(ev) || !identical(ev$group, "JY-EVNK") || !identical(ev$sample, "EVNK1-1") ||
    !identical(ev$bio_sample, "EVNK-1") || !identical(ev$panel, "P1") ||
    !identical(ev$replicate, "1") || !identical(ev$tech_rep, "1")) {
  fail(sprintf("EVNK1-1_P1_unmixed.fcs -> JY-EVNK / EVNK-1 / EVNK1-1 P1, got group=%s sample=%s bio=%s tech=%s panel=%s",
               if (is.null(ev)) "NULL" else ev$group,
               if (is.null(ev)) "NULL" else ev$sample,
               if (is.null(ev)) "NULL" else ev$bio_sample,
               if (is.null(ev)) "NULL" else ev$tech_rep,
               if (is.null(ev)) "NULL" else ev$panel))
}
ev2 <- parse_fcs_filename("JY-EVNK1-2_P2_unmixed.fcs")
if (is.null(ev2) || !identical(ev2$group, "JY-EVNK") || !identical(ev2$sample, "EVNK1-2") ||
    !identical(ev2$panel, "P2") || !identical(ev2$tech_rep, "2")) {
  fail("JY-EVNK1-2_P2_unmixed.fcs should parse as JY-EVNK sample EVNK1-2 P2")
}
nnk <- parse_fcs_filename("NNK-3-1_P3_unmixed.fcs")
if (is.null(nnk) || !identical(nnk$group, "JY-NNK") || !identical(nnk$sample, "NNK3-1") ||
    !identical(nnk$bio_sample, "NNK-3") || !identical(nnk$panel, "P3")) {
  fail(sprintf("NNK-3-1_P3 must be JY-NNK / NNK-3 / NNK3-1 P3, got group=%s sample=%s bio=%s panel=%s",
               if (is.null(nnk)) "NULL" else nnk$group,
               if (is.null(nnk)) "NULL" else nnk$sample,
               if (is.null(nnk)) "NULL" else nnk$bio_sample,
               if (is.null(nnk)) "NULL" else nnk$panel))
}
jy_nnk <- parse_fcs_filename("JY-NNK-2-2_P1_unmixed.fcs")
if (is.null(jy_nnk) || !identical(jy_nnk$group, "JY-NNK") || !identical(jy_nnk$sample, "NNK2-2") ||
    !identical(jy_nnk$bio_sample, "NNK-2")) {
  fail("JY-NNK-2-2_P1_unmixed.fcs should parse as JY-NNK / NNK-2 / NNK2-2")
}
bio_only <- parse_fcs_filename("EVNK-1_P1_unmixed.fcs")
if (is.null(bio_only) || !identical(bio_only$group, "JY-EVNK") || !identical(bio_only$sample, "EVNK-1") ||
    !identical(bio_only$bio_sample, "EVNK-1") || !is.na(bio_only$tech_rep)) {
  fail("EVNK-1_P1 (no tech) should be JY-EVNK / bio EVNK-1 with no tech_rep")
}
folder <- parse_fcs_filename("/data/JY-EVNK/1-1_P1_unmixed.fcs")
if (is.null(folder) || !identical(folder$group, "JY-EVNK") || !identical(folder$sample, "EVNK1-1")) {
  fail("folder JY-EVNK/1-1_P1_unmixed.fcs should be JY-EVNK / EVNK1-1")
}

# 不要把原实验 EV/H 或 ICI EV-1 认成本方案
if (!is.null(parse_fcs_filename("EV1-1_P1_unmixed.fcs"))) {
  fail("original EV1-1_P1 must not parse as JY EVNK/NNK")
}
if (!is.null(parse_fcs_filename("H-1-1_P2_unmixed.fcs"))) {
  fail("original H-1-1_P2 must not parse as JY")
}
if (!is.null(parse_fcs_filename("ZZX-EV-1_P1_unmixed.fcs"))) {
  fail("ICI-style ZZX-EV-1_P1 must not parse as JY")
}
evnk_as_nnk <- parse_fcs_filename("EVNK1-1_P1_unmixed.fcs")
if (!identical(evnk_as_nnk$group, "JY-EVNK")) fail("EVNK must parse as JY-EVNK, not JY-NNK")
expect(jy_canon_group("NNK"), "JY-NNK", "canon NNK")
expect(jy_canon_group("JY-EVNK"), "JY-EVNK", "canon JY-EVNK")
expect(unname(pal_group["JY-EVNK"]), "#1A1A1A", "JY-EVNK black")
expect(unname(pal_group["JY-NNK"]), "#E31A1C", "JY-NNK red")

# QC 与原免疫亚群一致：CD45+ + P1 淋巴门，不是 ICI 的只去双联体/死细胞
set.seed(7)
n_leu <- 80L
n_tgt <- 80L
n_dead <- 8L
exprs_qc <- cbind(
  `FSC-A` = c(rnorm(n_leu, 80000, 4000), rnorm(n_tgt, 110000, 5000), rnorm(n_dead, 80000, 4000)),
  `SSC-A` = c(rnorm(n_leu, 20000, 2000), rnorm(n_tgt, 45000, 3000), rnorm(n_dead, 20000, 2000)),
  `FSC-H` = c(rnorm(n_leu, 80000, 4000), rnorm(n_tgt, 110000, 5000), rnorm(n_dead, 80000, 4000)),
  `L/D` = c(rnorm(n_leu + n_tgt, 200, 30), rnorm(n_dead, 8000, 200)),
  CD45 = c(rnorm(n_leu, 12000, 800), rnorm(n_tgt, 80, 20), rnorm(n_dead, 12000, 800))
)
map_qc <- data.frame(
  marker = c("L/D", "CD45"),
  channel_index = c(4L, 5L),
  stringsAsFactors = FALSE
)
chn <- colnames(exprs_qc)
keep_jy <- qc_filter_matrix(exprs_qc, chn, map_qc, "P1")
tgt_idx <- seq.int(n_leu + 1L, n_leu + n_tgt)
if (mean(keep_jy[tgt_idx]) > 0.4) {
  fail("JY QC (CD45+ and P1 lymph gate) should drop most CD45- large cells, like the original immune scheme")
}
if (mean(keep_jy[seq_len(n_leu)]) < 0.7) {
  fail("JY QC should keep most CD45+ lymphocytes")
}

# 去极端生物学重复：三个里去掉离中位数更远的那个
vals <- c(1, 2, 10)
drop <- which_extreme_bio(vals)
if (!identical(as.integer(drop), 3L)) {
  fail(sprintf("which_extreme_bio should drop the far max, got %s", paste(drop, collapse = ",")))
}

# 三套方案完全独立：JY 不得 source Flow_* 或 ICI_*
jy_bundle <- c(
  "JY_Flow_dimred_pipeline.R",
  "JY_flow_engine.R",
  "JY_flow_panel_map.json",
  "JY_Flow_dimred_all_subsets.R",
  "JY_Flow_dimred_trajectory.R"
)
for (nm in jy_bundle) {
  if (!file.exists(file.path(root, nm))) {
    fail(sprintf("JY scheme must ship %s (do not borrow Flow_* or ICI_*)", nm))
  }
}
if (exists("load_ici_engine", mode = "function")) {
  fail("JY session must not load ICI entry helpers")
}
if (exists("load_jy_engine", mode = "function")) {
  engine_body <- paste(deparse(load_jy_engine), collapse = "\n")
  if (!grepl("JY_flow_engine\\.R", engine_body)) {
    fail("load_jy_engine must load JY_flow_engine.R")
  }
  if (grepl("source\\([^)]*Flow_dimred_pipeline\\.R", engine_body) ||
      grepl("pipe <- \"Flow_dimred_pipeline\\.R\"", engine_body) ||
      grepl("ICI_flow_engine", engine_body)) {
    fail("load_jy_engine must source JY_flow_engine.R only")
  }
}
for (nm in jy_bundle[grepl("\\.R$", jy_bundle)]) {
  txt <- paste(readLines(file.path(root, nm), warn = FALSE), collapse = "\n")
  if (grepl('source\\([\'"]Flow_dimred_pipeline\\.R[\'"]', txt) ||
      grepl('source\\([\'"]ICI_flow_engine\\.R[\'"]', txt) ||
      grepl('source\\([\'"]ICI_Flow_dimred_pipeline\\.R[\'"]', txt)) {
    fail(sprintf("%s must not source Flow_dimred_pipeline.R or ICI_*", nm))
  }
}
engine_txt <- paste(readLines(file.path(root, "JY_flow_engine.R"), warn = FALSE), collapse = "\n")
if (grepl("flow_panel_map\\.json", engine_txt) &&
    !grepl("JY_flow_panel_map\\.json", engine_txt)) {
  fail("JY_flow_engine.R must read JY_flow_panel_map.json, not flow_panel_map.json")
}
if (grepl("flow_primary_data_dir\\s*<-\\s*\"E:/R/fuction of cell\"", engine_txt)) {
  fail("JY_flow_engine.R must not use E:/R/fuction of cell as the data directory")
}
if (!grepl("E:/R/fuction of cell-ljy", engine_txt, fixed = TRUE)) {
  fail("JY_flow_engine.R must point at E:/R/fuction of cell-ljy")
}
if (grepl("Internation cell immune", engine_txt) &&
    grepl("flow_primary_data_dir\\s*<-\\s*\"E:/R/Internation cell immune\"", engine_txt)) {
  fail("JY must not use the Internation data directory")
}

expect(panel_map$panels$P1$markers[[match("Perforin", jy_p1)]]$fluorochrome, "FITC", "JY P1 Perforin-FITC")

# all_subsets 会覆盖 jy_keep_cand；路径里的反斜杠必须用 fixed gsub，不能当正则
Sys.setenv(FLOW_ALL_SUBSETS_FROM_PIPELINE = "1", FLOW_FUNCTIONS_ONLY = "1")
sys.source(file.path(root, "JY_Flow_dimred_all_subsets.R"), envir = .GlobalEnv)
win_p <- "E:\\R\\fuction of cell-ljy\\JY_Flow_dimred_trajectory.R"
keep_win <- tryCatch(jy_keep_cand(win_p), error = function(e) e)
if (inherits(keep_win, "error")) {
  fail(sprintf("jy_keep_cand must accept Windows paths after all_subsets: %s", keep_win$message))
}
if (!isTRUE(keep_win)) fail("jy_keep_cand must keep E:/R/fuction of cell-ljy paths")
Sys.setenv(FLOW_TRAJECTORY_FROM_PIPELINE = "1")
sys.source(file.path(root, "JY_Flow_dimred_trajectory.R"), envir = .GlobalEnv)
keep_tr <- tryCatch(jy_keep_cand(win_p), error = function(e) e)
if (inherits(keep_tr, "error")) {
  fail(sprintf("jy_keep_cand must accept Windows paths after trajectory: %s", keep_tr$message))
}

if (!exists("export_functional_state_figures", mode = "function")) {
  fail("JY engine must ship NKT/B functional-state export (do not source Flow_*)")
}
if (!identical(flow_comparison_tag(), "JY_NNK_vs_JY_EVNK")) {
  fail("JY functional-state files must be tagged JY_NNK_vs_JY_EVNK")
}
jy_p1 <- functional_state_specs("P1")
if (!length(jy_p1) || !identical(jy_p1[[1]]$parent, "NKT")) {
  fail("JY P1 functional-state must analyze NKT")
}
if (!"Activated_B" %in% vapply(functional_state_specs("P2"), function(s) s$parent, character(1))) {
  fail("JY P2 functional-state must include Activated_B")
}

cat("PASS test_jy_flow.R\n")
