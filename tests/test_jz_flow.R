# JZ：免疫亚群降维（E:/R/fuction of cell-wjz，比较 JZ-AB / JZ-EVB）
# run: Rscript tests/test_jz_flow.R
Sys.setenv(FLOW_FUNCTIONS_ONLY = "1", JZ_FUNCTIONS_ONLY = "1")
Sys.setenv(JZ_FLOW_DIR = tempfile("jz_flow_"))
dir.create(Sys.getenv("JZ_FLOW_DIR"), recursive = TRUE, showWarnings = FALSE)
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/test_jz_flow.R"
root <- dirname(dirname(normalizePath(this_file)))
source(file.path(root, "JZ_Flow_dimred_pipeline.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
expect <- function(got, want, tag) {
  if (!identical(got, want)) {
    fail(sprintf("%s: got %s, want %s", tag, paste(got, collapse = ","), paste(want, collapse = ",")))
  }
}

# 原 flow_panel_map.json / ICI 表不得被 JZ 改掉
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

jz_json <- jsonlite::fromJSON(file.path(root, "JZ_flow_panel_map.json"), simplifyVector = FALSE)
jz_p1 <- vapply(jz_json$panels$P1$markers, function(x) x$marker, character(1))
if (length(jz_p1) != 25L) fail(sprintf("JZ P1 should copy original 25 markers, got %s", length(jz_p1)))
if ("His" %in% jz_p1) fail("JZ P1 must not gain His (this is not the ICI sheet)")
if (!("IgD" %in% vapply(jz_json$panels$P2$markers, function(x) x$marker, character(1)))) {
  fail("JZ P2 must keep original IgD (immune-subset panel, not ICI 7-color)")
}
expect(jz_json$groups[[1]], "JZ-EVB", "JZ groups JZ-EVB")
expect(jz_json$groups[[2]], "JZ-AB", "JZ groups JZ-AB")
expect(jz_json$data_dir, "E:/R/fuction of cell-wjz", "JZ data_dir")
expect(jz_json$cohort, "JZ", "JZ cohort")
expect(jz_json$comparison, "JZ_AB_vs_JZ_EVB", "JZ comparison JZ-AB / JZ-EVB")

expect(flow_ctrl_group, "JZ-EVB", "loaded ctrl JZ-EVB")
expect(flow_trt_group, "JZ-AB", "loaded trt JZ-AB")
expect(as.character(flow_group_levels), c("JZ-EVB", "JZ-AB"), "group levels")
if (!isTRUE(flow_trim_bio_extremes)) {
  fail("JZ must drop 1 extreme bio-rep like the immune-subset scheme (n=2)")
}

# 文件名：EVB1-1 / AB1-2 / JZ-EVB；生物学重复 EVB-1；先匹配 EVB
ev <- parse_fcs_filename("EVB1-1_P1_unmixed.fcs")
if (is.null(ev) || !identical(ev$group, "JZ-EVB") || !identical(ev$sample, "EVB1-1") ||
    !identical(ev$bio_sample, "EVB-1") || !identical(ev$panel, "P1") ||
    !identical(ev$replicate, "1") || !identical(ev$tech_rep, "1")) {
  fail(sprintf("EVB1-1_P1_unmixed.fcs -> JZ-EVB / EVB-1 / EVB1-1 P1, got group=%s sample=%s bio=%s tech=%s panel=%s",
               if (is.null(ev)) "NULL" else ev$group,
               if (is.null(ev)) "NULL" else ev$sample,
               if (is.null(ev)) "NULL" else ev$bio_sample,
               if (is.null(ev)) "NULL" else ev$tech_rep,
               if (is.null(ev)) "NULL" else ev$panel))
}
ev2 <- parse_fcs_filename("JZ-EVB1-2_P2_unmixed.fcs")
if (is.null(ev2) || !identical(ev2$group, "JZ-EVB") || !identical(ev2$sample, "EVB1-2") ||
    !identical(ev2$panel, "P2") || !identical(ev2$tech_rep, "2")) {
  fail("JZ-EVB1-2_P2_unmixed.fcs should parse as JZ-EVB sample EVB1-2 P2")
}
nnk <- parse_fcs_filename("AB-3-1_P3_unmixed.fcs")
if (is.null(nnk) || !identical(nnk$group, "JZ-AB") || !identical(nnk$sample, "AB3-1") ||
    !identical(nnk$bio_sample, "AB-3") || !identical(nnk$panel, "P3")) {
  fail(sprintf("AB-3-1_P3 must be JZ-AB / AB-3 / AB3-1 P3, got group=%s sample=%s bio=%s panel=%s",
               if (is.null(nnk)) "NULL" else nnk$group,
               if (is.null(nnk)) "NULL" else nnk$sample,
               if (is.null(nnk)) "NULL" else nnk$bio_sample,
               if (is.null(nnk)) "NULL" else nnk$panel))
}
jz_nnk <- parse_fcs_filename("JZ-AB-2-2_P1_unmixed.fcs")
if (is.null(jz_nnk) || !identical(jz_nnk$group, "JZ-AB") || !identical(jz_nnk$sample, "AB2-2") ||
    !identical(jz_nnk$bio_sample, "AB-2")) {
  fail("JZ-AB-2-2_P1_unmixed.fcs should parse as JZ-AB / AB-2 / AB2-2")
}
bio_only <- parse_fcs_filename("EVB-1_P1_unmixed.fcs")
if (is.null(bio_only) || !identical(bio_only$group, "JZ-EVB") || !identical(bio_only$sample, "EVB-1") ||
    !identical(bio_only$bio_sample, "EVB-1") || !is.na(bio_only$tech_rep)) {
  fail("EVB-1_P1 (no tech) should be JZ-EVB / bio EVB-1 with no tech_rep")
}
folder <- parse_fcs_filename("/data/JZ-EVB/1-1_P1_unmixed.fcs")
if (is.null(folder) || !identical(folder$group, "JZ-EVB") || !identical(folder$sample, "EVB1-1")) {
  fail("folder JZ-EVB/1-1_P1_unmixed.fcs should be JZ-EVB / EVB1-1")
}

# 不要把原实验 EV/H 或 ICI EV-1 认成本方案
if (!is.null(parse_fcs_filename("EV1-1_P1_unmixed.fcs"))) {
  fail("original EV1-1_P1 must not parse as JZ EVB/AB")
}
if (!is.null(parse_fcs_filename("H-1-1_P2_unmixed.fcs"))) {
  fail("original H-1-1_P2 must not parse as JZ")
}
if (!is.null(parse_fcs_filename("ZZX-EV-1_P1_unmixed.fcs"))) {
  fail("ICI-style ZZX-EV-1_P1 must not parse as JZ")
}
if (!is.null(parse_fcs_filename("EVNK1-1_P1_unmixed.fcs"))) {
  fail("JY EVNK1-1_P1 must not parse as JZ")
}
if (!is.null(parse_fcs_filename("JY-NNK-2-2_P1_unmixed.fcs"))) {
  fail("JY-NNK files must not parse as JZ")
}
if (!is.null(parse_fcs_filename("NNK1-1_P1_unmixed.fcs"))) {
  fail("NNK1-1 must not parse as JZ AB")
}
evb_as_ab <- parse_fcs_filename("EVB1-1_P1_unmixed.fcs")
if (!identical(evb_as_ab$group, "JZ-EVB")) fail("EVB must parse as JZ-EVB, not JZ-AB")
expect(jz_canon_group("AB"), "JZ-AB", "canon AB")
expect(jz_canon_group("JZ-EVB"), "JZ-EVB", "canon JZ-EVB")
expect(unname(pal_group["JZ-EVB"]), "#1A1A1A", "JZ-EVB black")
expect(unname(pal_group["JZ-AB"]), "#E31A1C", "JZ-AB red")

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
keep_jz <- qc_filter_matrix(exprs_qc, chn, map_qc, "P1")
tgt_idx <- seq.int(n_leu + 1L, n_leu + n_tgt)
if (mean(keep_jz[tgt_idx]) > 0.4) {
  fail("JZ QC (CD45+ and P1 lymph gate) should drop most CD45- large cells, like the original immune scheme")
}
if (mean(keep_jz[seq_len(n_leu)]) < 0.7) {
  fail("JZ QC should keep most CD45+ lymphocytes")
}

# 去极端生物学重复：三个里去掉离中位数更远的那个
vals <- c(1, 2, 10)
drop <- which_extreme_bio(vals)
if (!identical(as.integer(drop), 3L)) {
  fail(sprintf("which_extreme_bio should drop the far max, got %s", paste(drop, collapse = ",")))
}

# 四套方案完全独立：JZ 不得 source Flow_* / ICI_* / JY_*
jz_bundle <- c(
  "JZ_Flow_dimred_pipeline.R",
  "JZ_flow_engine.R",
  "JZ_flow_panel_map.json",
  "JZ_Flow_dimred_all_subsets.R",
  "JZ_Flow_dimred_trajectory.R",
  "JZ_Flow_dimred_functional_state.R"
)
for (nm in jz_bundle) {
  if (!file.exists(file.path(root, nm))) {
    fail(sprintf("JZ scheme must ship %s (do not borrow Flow_* / ICI_* / JY_*)", nm))
  }
}
if (exists("load_ici_engine", mode = "function")) {
  fail("JZ session must not load ICI entry helpers")
}
if (exists("load_jy_engine", mode = "function")) {
  fail("JZ session must not load JY entry helpers")
}
if (exists("load_jz_engine", mode = "function")) {
  engine_body <- paste(deparse(load_jz_engine), collapse = "\n")
  if (!grepl("JZ_flow_engine\\.R", engine_body)) {
    fail("load_jz_engine must load JZ_flow_engine.R")
  }
  if (grepl("source\\([^)]*Flow_dimred_pipeline\\.R", engine_body) ||
      grepl("pipe <- \"Flow_dimred_pipeline\\.R\"", engine_body) ||
      grepl("ICI_flow_engine", engine_body) ||
      grepl("JY_flow_engine", engine_body)) {
    fail("load_jz_engine must source JZ_flow_engine.R only")
  }
}
for (nm in jz_bundle[grepl("\\.R$", jz_bundle)]) {
  txt <- paste(readLines(file.path(root, nm), warn = FALSE), collapse = "\n")
  if (grepl('source\\([\'"]Flow_dimred_pipeline\\.R[\'"]', txt) ||
      grepl('source\\([\'"]ICI_flow_engine\\.R[\'"]', txt) ||
      grepl('source\\([\'"]ICI_Flow_dimred_pipeline\\.R[\'"]', txt) ||
      grepl('source\\([\'"]JY_flow_engine\\.R[\'"]', txt) ||
      grepl('source\\([\'"]JY_Flow_dimred_pipeline\\.R[\'"]', txt)) {
    fail(sprintf("%s must not source Flow_dimred_pipeline.R, ICI_*, or JY_*", nm))
  }
}
engine_txt <- paste(readLines(file.path(root, "JZ_flow_engine.R"), warn = FALSE), collapse = "\n")
if (grepl("flow_panel_map\\.json", engine_txt) &&
    !grepl("JZ_flow_panel_map\\.json", engine_txt)) {
  fail("JZ_flow_engine.R must read JZ_flow_panel_map.json, not flow_panel_map.json")
}
if (grepl("flow_primary_data_dir\\s*<-\\s*\"E:/R/fuction of cell\"", engine_txt)) {
  fail("JZ_flow_engine.R must not use E:/R/fuction of cell as the data directory")
}
if (!grepl("E:/R/fuction of cell-wjz", engine_txt, fixed = TRUE)) {
  fail("JZ_flow_engine.R must point at E:/R/fuction of cell-wjz")
}
if (grepl("Internation cell immune", engine_txt) &&
    grepl("flow_primary_data_dir\\s*<-\\s*\"E:/R/Internation cell immune\"", engine_txt)) {
  fail("JZ must not use the Internation data directory")
}

if (grepl("flow_primary_data_dir\\s*<-\\s*\"E:/R/fuction of cell-ljy\"", engine_txt)) {
  fail("JZ_flow_engine.R must not use E:/R/fuction of cell-ljy as the data directory")
}

jy_json <- jsonlite::fromJSON(file.path(root, "JY_flow_panel_map.json"), simplifyVector = FALSE)
expect(jy_json$groups[[1]], "JY-EVNK", "JY map still JY-EVNK")
expect(jy_json$groups[[2]], "JY-NNK", "JY map still JY-NNK")
expect(jy_json$data_dir, "E:/R/fuction of cell-ljy", "JY data_dir unchanged")
expect(panel_map$panels$P1$markers[[match("Perforin", jz_p1)]]$fluorochrome, "FITC", "JZ P1 Perforin-FITC")

Sys.setenv(FLOW_ALL_SUBSETS_FROM_PIPELINE = "1", FLOW_FUNCTIONS_ONLY = "1")
sys.source(file.path(root, "JZ_Flow_dimred_all_subsets.R"), envir = .GlobalEnv)
win_p <- "E:\\R\\fuction of cell-wjz\\JZ_Flow_dimred_trajectory.R"
keep_win <- tryCatch(jz_keep_cand(win_p), error = function(e) e)
if (inherits(keep_win, "error")) {
  fail(sprintf("jz_keep_cand must accept Windows paths after all_subsets: %s", keep_win$message))
}
if (!isTRUE(keep_win)) fail("jz_keep_cand must keep E:/R/fuction of cell-wjz paths")

if (!exists("export_functional_state_from_results", mode = "function")) {
  fail("JZ engine must allow rerunning functional-state from embeddings")
}
if (!identical(flow_comparison_tag(), "JZ_AB_vs_JZ_EVB")) {
  fail("JZ functional-state files must be tagged JZ_AB_vs_JZ_EVB")
}
jz_p1 <- functional_state_specs("P1")
if (!length(jz_p1) || !identical(jz_p1[[2]]$state, "exhaustion")) {
  fail("JZ P1 functional-state must include NKT exhaustion")
}

cat("PASS test_jz_flow.R\n")
