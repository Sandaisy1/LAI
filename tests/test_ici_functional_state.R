# ICI His+ 亚群功能状态（run: Rscript tests/test_ici_functional_state.R）
Sys.setenv(FLOW_FUNCTIONS_ONLY = "1", ICI_FUNCTIONS_ONLY = "1")
Sys.setenv(ICI_FLOW_DIR = tempfile("ici_func_"))
dir.create(Sys.getenv("ICI_FLOW_DIR"), recursive = TRUE, showWarnings = FALSE)
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/test_ici_functional_state.R"
root <- dirname(dirname(normalizePath(this_file)))
source(file.path(root, "ICI_Flow_dimred_pipeline.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}

if (!exists("export_functional_state_from_results", mode = "function") ||
    !exists("func_p3_state_parents", mode = "function") ||
    !exists("ici_keep_his_parent_cells", mode = "function")) {
  fail("ICI engine must ship His+ functional-state helpers")
}
if (!isFALSE(flow_trim_bio_extremes)) {
  fail("ICI functional-state must keep n=3")
}

p1s <- functional_state_specs("P1")
p1_parents <- unique(vapply(p1s, function(s) s$parent, character(1)))
want_p1 <- c("CD4", "CD4_naive", "CD4_TCM", "CD4_TEM", "Treg", "CD8", "NK", "NKT", "NKT_CD4")
if (!all(want_p1 %in% p1_parents)) {
  fail(sprintf("P1 missing parents: %s", paste(setdiff(want_p1, p1_parents), collapse = ",")))
}
circ <- c("CD4_activated", "CD4_exhausted", "NK_exhausted", "NKT_activated")
if (any(circ %in% p1_parents)) fail("P1 must not use circular activated/exhausted labels as parents")
cd4_act <- Filter(function(s) identical(s$parent, "CD4") && identical(s$state, "activation_effector"), p1s)
if (!length(cd4_act) || !("CD69" %in% cd4_act[[1]]$markers)) {
  fail("P1 His+ CD4 activation must include CD69")
}

p2s <- functional_state_specs("P2")
p2_parents <- unique(vapply(p2s, function(s) s$parent, character(1)))
if (!all(c("CD19", "Naive_B", "Unswitched_B", "Switched_B", "MZ_B") %in% p2_parents)) {
  fail("P2 must list B subset parents even if CD86/CD80/CD40 are missing on this sheet")
}
if (any(vapply(p2s, function(s) identical(s$state, "exhaustion"), logical(1)))) {
  fail("P2 must not invent B-cell exhaustion")
}

p3s <- functional_state_specs("P3")
p3_parents <- unique(vapply(p3s, function(s) s$parent, character(1)))
want_p3 <- c("Neutrophil", "Macrophage", "DC", "Mono_Ly6Chi", "Mono_Ly6Clo")
if (!all(want_p3 %in% p3_parents)) {
  fail(sprintf("P3 missing myeloid parents: %s", paste(setdiff(want_p3, p3_parents), collapse = ",")))
}
if (any(c("T", "B", "NK") %in% p3_parents)) fail("P3 dump T/B/NK must not get functional-state")
if ("exhaustion" %in% vapply(p3s, function(s) s$state, character(1))) {
  fail("P3 suppression must not be named exhaustion")
}

# His- 不得进入功能统计
set.seed(4)
mk <- function(sample, group, n, his, cd69, lin = "CD4_naive") {
  data.frame(
    sample = rep(sample, n),
    group = rep(group, n),
    bio_sample = rep(sample, n),
    tech_rep = NA_character_,
    lineage = lin,
    cluster_lineage = "CD4",
    His = rnorm(n, his, 0.1),
    CD69 = rnorm(n, cd69, 0.15),
    CD4 = rnorm(n, 3.0, 0.15),
    CD3 = rnorm(n, 3.1, 0.15),
    stringsAsFactors = FALSE
  )
}
mixed <- rbind(
  mk("EV-1", "EV", 40, 3.2, 0.6),
  mk("EV-1", "EV", 40, 0.3, 3.2),
  mk("H-1", "H", 40, 3.2, 2.8),
  mk("H-1", "H", 40, 0.3, 3.2)
)
kept <- ici_keep_his_parent_cells(mixed)
if (nrow(kept) > 90 || any(kept$His < 1.2)) {
  fail(sprintf("functional-state must drop His- cells, kept n=%s", nrow(kept)))
}

# His+ CD4 有 CD69 时应出图；His- 高 CD69 不得把对照阳性率抬上去
set.seed(8)
cells_p1 <- rbind(
  mk("EV-1", "EV", 50, 3.2, 0.5),
  mk("EV-2", "EV", 50, 3.2, 0.5),
  mk("EV-3", "EV", 50, 3.2, 0.5),
  mk("H-1", "H", 50, 3.2, 3.0),
  mk("H-2", "H", 50, 3.2, 3.0),
  mk("H-3", "H", 50, 3.2, 3.0),
  mk("EV-1", "EV", 50, 0.2, 3.3)
)
td <- tempfile("ici_func_p1")
dir.create(td, recursive = TRUE)
ok <- tryCatch({
  export_functional_state_figures(cells_p1, "P1", td)
  TRUE
}, error = function(e) {
  cat("ICI P1 functional-state crashed:", e$message, "\n")
  FALSE
})
if (!isTRUE(ok)) fail("His+ P1 functional-state must not crash")
stub <- file.path(td, "functional_state", "P1_CD4_activation_effector_H_vs_EV.pdf")
if (!file.exists(stub)) fail("P1 His+ CD4 activation figure missing")
samp <- read.csv(file.path(td, "functional_state", "P1_CD4_CD69_H_vs_EV_by_sample.csv"),
                 stringsAsFactors = FALSE)
ev1 <- samp$percent[samp$sample == "EV-1" & samp$marker == "CD69"]
if (!length(ev1) || ev1 > 40) {
  fail(sprintf("His- high-CD69 cells must not inflate EV-1 percent, got %s", paste(ev1, collapse = ",")))
}
unlink(td, recursive = TRUE)

# P2 无 CD86/CD80/CD40：skip，不炸
set.seed(9)
cells_p2 <- data.frame(
  sample = rep(c("EV-1", "EV-2", "EV-3", "H-1", "H-2", "H-3"), each = 40),
  group = rep(c("EV", "EV", "EV", "H", "H", "H"), each = 40),
  bio_sample = rep(c("EV-1", "EV-2", "EV-3", "H-1", "H-2", "H-3"), each = 40),
  lineage = "Naive_B",
  cluster_lineage = "Naive_B",
  His = rnorm(240, 3.2, 0.1),
  CD19 = rnorm(240, 3.0, 0.15),
  CD27 = rnorm(240, 0.3, 0.1),
  stringsAsFactors = FALSE
)
td2 <- tempfile("ici_func_p2")
dir.create(td2, recursive = TRUE)
ok2 <- tryCatch({
  export_functional_state_figures(cells_p2, "P2", td2)
  TRUE
}, error = function(e) {
  cat("ICI P2 functional-state crashed:", e$message, "\n")
  FALSE
})
if (!isTRUE(ok2)) fail("P2 missing CD86/CD80/CD40 must skip, not crash")
if (length(list.files(file.path(td2, "functional_state"), pattern = "activation.*\\.pdf$"))) {
  fail("P2 must not invent CD86/CD80/CD40 figures on this sheet")
}
unlink(td2, recursive = TRUE)

# P3 His+ 巨噬：CD86 应出图
set.seed(10)
mk3 <- function(sample, group, n, cd86) {
  data.frame(
    sample = rep(sample, n),
    group = rep(group, n),
    bio_sample = rep(sample, n),
    lineage = "Macrophage",
    cluster_lineage = "Myeloid",
    His = rnorm(n, 3.2, 0.1),
    `F4/80` = rnorm(n, 3.1, 0.15),
    CD11B = rnorm(n, 3.0, 0.15),
    CD86 = rnorm(n, cd86, 0.15),
    CD80 = rnorm(n, cd86 - 0.2, 0.15),
    CD40 = rnorm(n, cd86 - 0.3, 0.15),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
cells_p3 <- rbind(
  mk3("EV-1", "EV", 40, 0.6), mk3("EV-2", "EV", 40, 0.6), mk3("EV-3", "EV", 40, 0.6),
  mk3("H-1", "H", 40, 2.9), mk3("H-2", "H", 40, 2.9), mk3("H-3", "H", 40, 2.9)
)
td3 <- tempfile("ici_func_p3")
dir.create(td3, recursive = TRUE)
ok3 <- tryCatch({
  export_functional_state_figures(cells_p3, "P3", td3)
  TRUE
}, error = function(e) {
  cat("ICI P3 functional-state crashed:", e$message, "\n")
  FALSE
})
if (!isTRUE(ok3)) fail("P3 His+ macrophage functional-state must not crash")
if (!file.exists(file.path(td3, "functional_state", "P3_Macrophage_activation_H_vs_EV.pdf"))) {
  fail("P3 His+ macrophage activation figure missing")
}
unlink(td3, recursive = TRUE)

txt <- paste(readLines(file.path(root, "ICI_Flow_dimred_functional_state.R"), warn = FALSE), collapse = "\n")
if (grepl('source\\([\'"]Flow_dimred_pipeline\\.R[\'"]', txt) ||
    grepl('source\\([\'"]JY_|source\\([\'"]JZ_', txt)) {
  fail("ICI_Flow_dimred_functional_state.R must not source Flow_/JY_/JZ_")
}

cat("PASS test_ici_functional_state.R\n")
