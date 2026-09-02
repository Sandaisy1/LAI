# T-T6：免疫亚群降维（E:/R/flow J，比较 T6 vs T）
# run: Rscript tests/test_t_t6_flow.R
Sys.setenv(FLOW_FUNCTIONS_ONLY = "1", T_T6_FUNCTIONS_ONLY = "1")
Sys.setenv(T_T6_FLOW_DIR = tempfile("tt6_flow_"))
dir.create(Sys.getenv("T_T6_FLOW_DIR"), recursive = TRUE, showWarnings = FALSE)
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/test_t_t6_flow.R"
root <- dirname(dirname(normalizePath(this_file)))
source(file.path(root, "T_T6_Flow_dimred_pipeline.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
expect <- function(got, want, tag) {
  if (!identical(got, want)) {
    fail(sprintf("%s: got %s, want %s", tag, paste(got, collapse = ","), paste(want, collapse = ",")))
  }
}

# 原 flow / ICI / JY / JZ 表不得被 T-T6 改掉
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
expect(jy_json$groups[[1]], "JY-EVNK", "JY groups unchanged")
expect(jy_json$data_dir, "E:/R/fuction of cell-ljy", "JY data_dir unchanged")

jz_json <- jsonlite::fromJSON(file.path(root, "JZ_flow_panel_map.json"), simplifyVector = FALSE)
expect(jz_json$groups[[1]], "JZ-EVB", "JZ groups unchanged")

tt6_json <- jsonlite::fromJSON(file.path(root, "T_T6_flow_panel_map.json"), simplifyVector = FALSE)
tt6_p1 <- vapply(tt6_json$panels$P1$markers, function(x) x$marker, character(1))
tt6_p2 <- vapply(tt6_json$panels$P2$markers, function(x) x$marker, character(1))
tt6_p3 <- vapply(tt6_json$panels$P3$markers, function(x) x$marker, character(1))
if (length(tt6_p1) != 25L) fail(sprintf("T-T6 P1 should have 25 markers, got %s", length(tt6_p1)))
if ("His" %in% c(tt6_p1, tt6_p2, tt6_p3)) fail("T-T6 panels must not gain His")
if (!("Perforin" %in% tt6_p1)) fail("T-T6 P1 must include Perforin")
if (!("IgD" %in% tt6_p2) || !("BLIMP-1" %in% tt6_p2)) {
  fail("T-T6 P2 must keep IgD and BLIMP-1 (not ICI 7-color)")
}
if (!("iNOS" %in% tt6_p3) || !("ARG-1" %in% tt6_p3)) {
  fail("T-T6 P3 must keep iNOS and ARG-1")
}
expect(tt6_json$panels$P1$markers[[match("Perforin", tt6_p1)]]$fluorochrome, "FITC", "P1 Perforin-FITC")
expect(tt6_json$panels$P2$markers[[match("IgD", tt6_p2)]]$fluorochrome, "FITC", "P2 IgD-FITC")
expect(tt6_json$panels$P2$markers[[match("CD80", tt6_p2)]]$fluorochrome, "BUV496", "P2 CD80-BUV496")
expect(tt6_json$panels$P3$markers[[match("iNOS", tt6_p3)]]$fluorochrome, "AF488", "P3 iNOS-AF488")
expect(tt6_json$groups[[1]], "T", "T-T6 groups T")
expect(tt6_json$groups[[2]], "T6", "T-T6 groups T6")
expect(tt6_json$data_dir, "E:/R/flow J", "T-T6 data_dir")
expect(tt6_json$cohort, "T-T6", "T-T6 cohort")
expect(tt6_json$comparison, "T6_vs_T", "T-T6 comparison T6 vs T")

expect(flow_ctrl_group, "T", "loaded ctrl T")
expect(flow_trt_group, "T6", "loaded trt T6")
expect(as.character(flow_group_levels), c("T", "T6"), "group levels")
if (isTRUE(flow_trim_bio_extremes)) {
  fail("T-T6 must keep all 3 bio-reps (n=3); do not drop an extreme")
}

# 文件名：T-1_P1 / T6-1_P1；先匹配 T6
t1 <- parse_fcs_filename("T-1_P1.fcs")
if (is.null(t1) || !identical(t1$group, "T") || !identical(t1$sample, "T-1") ||
    !identical(t1$bio_sample, "T-1") || !identical(t1$panel, "P1") ||
    !identical(t1$replicate, "1") || !is.na(t1$tech_rep)) {
  fail(sprintf("T-1_P1.fcs -> T / T-1 / P1, got group=%s sample=%s bio=%s tech=%s panel=%s",
               if (is.null(t1)) "NULL" else t1$group,
               if (is.null(t1)) "NULL" else t1$sample,
               if (is.null(t1)) "NULL" else t1$bio_sample,
               if (is.null(t1)) "NULL" else t1$tech_rep,
               if (is.null(t1)) "NULL" else t1$panel))
}
t6 <- parse_fcs_filename("T6-1_P1_unmixed.fcs")
if (is.null(t6) || !identical(t6$group, "T6") || !identical(t6$sample, "T6-1") ||
    !identical(t6$bio_sample, "T6-1") || !identical(t6$panel, "P1")) {
  fail(sprintf("T6-1_P1_unmixed.fcs -> T6 / T6-1 / P1, got group=%s sample=%s bio=%s panel=%s",
               if (is.null(t6)) "NULL" else t6$group,
               if (is.null(t6)) "NULL" else t6$sample,
               if (is.null(t6)) "NULL" else t6$bio_sample,
               if (is.null(t6)) "NULL" else t6$panel))
}
t6b <- parse_fcs_filename("T6-2_P3.fcs")
if (is.null(t6b) || !identical(t6b$group, "T6") || !identical(t6b$panel, "P3") ||
    !identical(t6b$bio_sample, "T6-2")) {
  fail("T6-2_P3.fcs should parse as T6 / T6-2 / P3")
}
t_un <- parse_fcs_filename("T-3_P2_unmixed.fcs")
if (is.null(t_un) || !identical(t_un$group, "T") || !identical(t_un$panel, "P2") ||
    !identical(t_un$bio_sample, "T-3")) {
  fail("T-3_P2_unmixed.fcs should parse as T / T-3 / P2")
}
tech <- parse_fcs_filename("T6-1-2_P1_unmixed.fcs")
if (is.null(tech) || !identical(tech$group, "T6") || !identical(tech$sample, "T6-1-2") ||
    !identical(tech$bio_sample, "T6-1") || !identical(tech$tech_rep, "2")) {
  fail("T6-1-2_P1_unmixed.fcs should keep group T6 with tech 2")
}
folder <- parse_fcs_filename("/data/T6/1-1_P1_unmixed.fcs")
if (is.null(folder) || !identical(folder$group, "T6") || !identical(folder$sample, "T6-1-1")) {
  fail("folder T6/1-1_P1_unmixed.fcs should be T6 / T6-1-1")
}

# 不要把其他实验认成本方案
if (!is.null(parse_fcs_filename("EV1-1_P1_unmixed.fcs"))) {
  fail("original EV1-1_P1 must not parse as T/T6")
}
if (!is.null(parse_fcs_filename("H-1-1_P2_unmixed.fcs"))) {
  fail("original H-1-1_P2 must not parse as T/T6")
}
if (!is.null(parse_fcs_filename("ZZX-EV-1_P1_unmixed.fcs"))) {
  fail("ICI-style ZZX-EV-1_P1 must not parse as T/T6")
}
if (!is.null(parse_fcs_filename("EVNK1-1_P1_unmixed.fcs"))) {
  fail("JY EVNK1-1_P1 must not parse as T/T6")
}
if (!is.null(parse_fcs_filename("EVB1-1_P1_unmixed.fcs"))) {
  fail("JZ EVB1-1_P1 must not parse as T/T6")
}
expect(tt6_canon_group("T6"), "T6", "canon T6")
expect(tt6_canon_group("T"), "T", "canon T")
expect(tt6_canon_group("T6-1"), "T6", "canon T6-1")
if (!identical(tt6_canon_group("EVNK"), NA_character_)) fail("EVNK must not canon to T/T6")
expect(unname(pal_group["T"]), "#1A1A1A", "T black")
expect(unname(pal_group["T6"]), "#E31A1C", "T6 red")

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
keep_tt6 <- qc_filter_matrix(exprs_qc, chn, map_qc, "P1")
tgt_idx <- seq.int(n_leu + 1L, n_leu + n_tgt)
if (mean(keep_tt6[tgt_idx]) > 0.4) {
  fail("T-T6 QC (CD45+ and P1 lymph gate) should drop most CD45- large cells")
}
if (mean(keep_tt6[seq_len(n_leu)]) < 0.7) {
  fail("T-T6 QC should keep most CD45+ lymphocytes")
}

# 五套方案完全独立：T-T6 不得 source Flow_* / ICI_* / JY_* / JZ_*
tt6_bundle <- c(
  "T_T6_Flow_dimred_pipeline.R",
  "T_T6_flow_engine.R",
  "T_T6_flow_panel_map.json",
  "T_T6_Flow_dimred_all_subsets.R",
  "T_T6_Flow_dimred_trajectory.R",
  "T_T6_Flow_dimred_functional_state.R"
)
for (nm in tt6_bundle) {
  if (!file.exists(file.path(root, nm))) {
    fail(sprintf("T-T6 scheme must ship %s (do not borrow Flow_* / ICI_* / JY_* / JZ_*)", nm))
  }
}
if (exists("load_ici_engine", mode = "function") ||
    exists("load_jy_engine", mode = "function") ||
    exists("load_jz_engine", mode = "function")) {
  fail("T-T6 session must not load ICI / JY / JZ entry helpers")
}
if (exists("load_tt6_engine", mode = "function")) {
  engine_body <- paste(deparse(load_tt6_engine), collapse = "\n")
  if (!grepl("T_T6_flow_engine\\.R", engine_body)) {
    fail("load_tt6_engine must load T_T6_flow_engine.R")
  }
  if (grepl("source\\([^)]*Flow_dimred_pipeline\\.R", engine_body) ||
      grepl("ICI_flow_engine", engine_body) ||
      grepl("JY_flow_engine", engine_body) ||
      grepl("JZ_flow_engine", engine_body)) {
    fail("load_tt6_engine must source T_T6_flow_engine.R only")
  }
}
for (nm in tt6_bundle[grepl("\\.R$", tt6_bundle)]) {
  txt <- paste(readLines(file.path(root, nm), warn = FALSE), collapse = "\n")
  if (grepl('source\\([\'"]Flow_dimred_pipeline\\.R[\'"]', txt) ||
      grepl('source\\([\'"]ICI_flow_engine\\.R[\'"]', txt) ||
      grepl('source\\([\'"]JY_flow_engine\\.R[\'"]', txt) ||
      grepl('source\\([\'"]JZ_flow_engine\\.R[\'"]', txt) ||
      grepl('source\\([\'"]ICI_Flow_dimred_pipeline\\.R[\'"]', txt) ||
      grepl('source\\([\'"]JY_Flow_dimred_pipeline\\.R[\'"]', txt) ||
      grepl('source\\([\'"]JZ_Flow_dimred_pipeline\\.R[\'"]', txt)) {
    fail(sprintf("%s must not source Flow_* / ICI_* / JY_* / JZ_*", nm))
  }
}
engine_txt <- paste(readLines(file.path(root, "T_T6_flow_engine.R"), warn = FALSE), collapse = "\n")
if (grepl("flow_panel_map\\.json", engine_txt) &&
    !grepl("T_T6_flow_panel_map\\.json", engine_txt)) {
  fail("T_T6_flow_engine.R must read T_T6_flow_panel_map.json, not flow_panel_map.json")
}
if (grepl("flow_primary_data_dir\\s*<-\\s*\"E:/R/fuction of cell\"", engine_txt)) {
  fail("T_T6_flow_engine.R must not use E:/R/fuction of cell as the data directory")
}
if (!grepl("E:/R/flow J", engine_txt, fixed = TRUE)) {
  fail("T_T6_flow_engine.R must point at E:/R/flow J")
}

# all_subsets 会覆盖 tt6_keep_cand；路径里的反斜杠必须用 fixed gsub
Sys.setenv(FLOW_ALL_SUBSETS_FROM_PIPELINE = "1", FLOW_FUNCTIONS_ONLY = "1")
sys.source(file.path(root, "T_T6_Flow_dimred_all_subsets.R"), envir = .GlobalEnv)
win_p <- "E:\\R\\flow J\\T_T6_Flow_dimred_trajectory.R"
keep_win <- tryCatch(tt6_keep_cand(win_p), error = function(e) e)
if (inherits(keep_win, "error")) {
  fail(sprintf("tt6_keep_cand must accept Windows paths after all_subsets: %s", keep_win$message))
}
if (!isTRUE(keep_win)) fail("tt6_keep_cand must keep E:/R/flow J paths")
if (isTRUE(tt6_keep_cand("E:\\R\\fuction of cell-ljy\\JY_flow_engine.R"))) {
  fail("tt6_keep_cand must reject JY data directory")
}
Sys.setenv(FLOW_TRAJECTORY_FROM_PIPELINE = "1")
sys.source(file.path(root, "T_T6_Flow_dimred_trajectory.R"), envir = .GlobalEnv)
keep_tr <- tryCatch(tt6_keep_cand(win_p), error = function(e) e)
if (inherits(keep_tr, "error")) {
  fail(sprintf("tt6_keep_cand must accept Windows paths after trajectory: %s", keep_tr$message))
}

if (!exists("export_functional_state_from_results", mode = "function")) {
  fail("T-T6 engine must allow rerunning functional-state from embeddings")
}
if (!identical(flow_comparison_tag(), "T6_vs_T")) {
  fail("T-T6 functional-state files must be tagged T6_vs_T")
}
tt6_p1s <- functional_state_specs("P1")
tt6_p1_parents <- unique(vapply(tt6_p1s, function(s) s$parent, character(1)))
if (!all(c("CD4", "CD8", "NK", "NKT", "Treg") %in% tt6_p1_parents)) {
  fail("T-T6 P1 functional-state must analyze CD4/CD8/NK/NKT/Treg")
}
act_mks <- unique(unlist(lapply(tt6_p1s[vapply(tt6_p1s, function(s) identical(s$state, "activation_effector"), logical(1))], function(s) s$markers)))
if (!("Perforin" %in% act_mks)) {
  fail("T-T6 P1 activation/effector must include Perforin")
}
tt6_p2_parents <- vapply(functional_state_specs("P2"), function(s) s$parent, character(1))
if (!all(c("Activated_B", "CD19", "MZ_B", "Plasma") %in% tt6_p2_parents)) {
  fail("T-T6 P2 functional-state must include CD19, Activated_B, MZ_B, Plasma")
}
tt6_p3_parents <- unique(vapply(functional_state_specs("P3"), function(s) s$parent, character(1)))
if (!all(c("Macrophage", "Neutrophil", "DC") %in% tt6_p3_parents)) {
  fail("T-T6 P3 functional-state must analyze myeloid subsets")
}
if (any(c("T", "B", "NK") %in% tt6_p3_parents)) {
  fail("T-T6 P3 must not treat dump T/B/NK as functional-state parents")
}

export_functional_state_from_results <- function(...) "stale"
if (exists("func_p3_state_parents", mode = "function")) {
  rm(list = "func_p3_state_parents", envir = .GlobalEnv)
}
tt6_engine_loaded <- TRUE
Sys.setenv(FLOW_FUNCTIONAL_STATE_FROM_PIPELINE = "1", FLOW_FUNCTIONS_ONLY = "1", T_T6_FUNCTIONS_ONLY = "1")
sys.source(file.path(root, "T_T6_Flow_dimred_functional_state.R"), envir = .GlobalEnv)
if (!exists("export_functional_state_from_results", mode = "function")) {
  fail("T-T6 functional-state script must reload the engine if the rerun helper is missing")
}
if (!exists("func_p3_state_parents", mode = "function")) {
  fail("T-T6 standalone script must reload all-subset helpers even if tt6_engine_loaded is TRUE")
}
if (!length(functional_state_specs("P3"))) {
  fail("T-T6 reloaded engine must include P3 myeloid functional-state")
}

cat("PASS test_t_t6_flow.R\n")
