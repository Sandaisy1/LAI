# P1 T/NK subset labels (run: Rscript tests/test_flow_p1_annotate.R)
Sys.setenv(FLOW_FUNCTIONS_ONLY = "1")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/test_flow_p1_annotate.R"
source(file.path(dirname(normalizePath(this_file)), "..", "Flow_dimred_pipeline.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
expect <- function(got, want, tag) {
  if (!identical(got, want)) fail(sprintf("%s: got %s, want %s", tag, got, want))
}

mk <- function(...) {
  pops <- demo_means_p1()
  row <- pops[[1]]
  row[] <- 0.2
  extra <- list(...)
  for (nm in names(extra)) row[[nm]] <- extra[[nm]]
  as.data.frame(t(row), stringsAsFactors = FALSE)
}
lab_of <- function(row) {
  rownames(row) <- "C1"
  annotate_clusters(row, "P1")$lineage
}

expect(lab_of(mk(CD3 = 3.2, CD4 = 3.0, CD8 = 0.1, CD62L = 3.2, CD44 = 0.3)), "CD4_naive", "cd4-naive")
expect(lab_of(mk(CD3 = 3.2, CD4 = 3.0, CD8 = 0.1, CD62L = 2.8, CD44 = 2.8)), "CD4_TCM", "cd4-tcm")
expect(lab_of(mk(CD3 = 3.2, CD4 = 3.0, CD8 = 0.1, CD62L = 0.3, CD44 = 3.1)), "CD4_TEM", "cd4-tem")
expect(lab_of(mk(CD3 = 3.1, CD4 = 3.0, CD25 = 3.2, CD69 = 0.4, CD44 = 1.8)), "Treg", "treg")
expect(lab_of(mk(CD3 = 3.2, CD4 = 3.0, CD69 = 3.1, CD44 = 2.6, CD62L = 0.4, `TNF-a` = 2.4)), "CD4_activated", "cd4-act")
expect(lab_of(mk(CD3 = 3.3, CD8 = 3.1, CD8b = 2.9, CD4 = 0.1, CD62L = 3.1, CD44 = 0.3)), "CD8_naive", "cd8-naive")
expect(lab_of(mk(CD3 = 3.3, CD8 = 3.1, CD8b = 2.8, CD62L = 2.7, CD44 = 2.6)), "CD8_TCM", "cd8-tcm")
expect(lab_of(mk(CD3 = 3.3, CD8 = 3.1, CD8b = 2.8, CD44 = 3.0, CD62L = 0.3)), "CD8_TEM", "cd8-tem")
expect(lab_of(mk(CD3 = 3.3, CD8 = 3.1, CD44 = 2.6, CD62L = 0.4, GZMB = 3.0, Perforin = 2.8)), "CD8_effector", "cd8-eff")
expect(lab_of(mk(CD3 = 3.0, CD8 = 3.0, `LAG-3` = 3.1, `TIM-3` = 2.9, CD44 = 2.8)), "CD8_exhausted", "cd8-exh")
expect(lab_of(mk(CD3 = 0.1, `NK1.1` = 3.2, NKp46 = 3.0, CD8 = 1.0)), "NK", "nk")
expect(lab_of(mk(CD3 = 3.0, `NK1.1` = 2.8, NKp46 = 2.2, CD4 = 1.2, CD44 = 2.4)), "NKT", "nkt")
expect(lab_of(mk(CD19 = 3.3, CD3 = 0.1, CD4 = 0.2, CD8 = 0.2)), "B", "b")
expect(lab_of(mk(CD11B = 3.2, CD3 = 0.1, CD19 = 0.1, CD4 = 0.2)), "Myeloid", "myeloid")

# IFN-g / GZMB 背景不得把 TEM 并成 activated / effector
expect(
  lab_of(mk(CD3 = 3.2, CD4 = 3.0, CD8 = 0.1, CD62L = 0.3, CD44 = 3.1, `IFN-g` = 2.0)),
  "CD4_TEM", "ifng-not-act"
)
expect(
  lab_of(mk(CD3 = 3.3, CD8 = 3.1, CD4 = 0.1, CD62L = 0.3, CD44 = 3.0, GZMB = 0.8, `IFN-g` = 1.8)),
  "CD8_TEM", "ifng-not-eff"
)
# NK1.1 背景不得把 CD4 标成 NK
expect(
  lab_of(mk(CD3 = 3.2, CD4 = 3.0, CD8 = 0.2, `NK1.1` = 1.4, NKp46 = 1.2, CD62L = 3.0, CD44 = 0.4)),
  "CD4_naive", "nk-bg-not-nk"
)

cat("OK: P1 T/NK subset labels\n")
