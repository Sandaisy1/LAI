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

# 分层：先 CD4/CD8，再在类内分亚群；不得出现 NA
set.seed(42)
pops <- demo_means_p1()
feat <- names(pops[[1]])
bind_pop <- function(name, n) {
  m <- matrix(unlist(pops[[name]]), n, length(feat), byrow = TRUE)
  colnames(m) <- feat
  m + matrix(rnorm(n * length(feat), 0, 0.12), n)
}
mat <- rbind(
  bind_pop("CD4_naive", 80), bind_pop("CD4_TCM", 80), bind_pop("CD4_TEM", 80),
  bind_pop("Treg", 80), bind_pop("CD4_act", 80),
  bind_pop("CD8_naive", 80), bind_pop("CD8_TCM", 80), bind_pop("CD8_TEM", 80),
  bind_pop("CD8_eff", 80), bind_pop("CD8_exh", 80)
)
h <- hierarchical_gate(mat, "P1")
if (anyNA(h$subset) || any(!nzchar(h$subset))) fail("hierarchical labels contain NA")
if (any(h$major == "T") || any(grepl("^T_", h$subset))) fail("T TEM must not be used; assign CD4 or CD8")
n4 <- sum(h$major == "CD4")
n8 <- sum(h$major == "CD8")
if (n4 < 300) fail(sprintf("layer1 CD4 too few: %s", n4))
if (n8 < 300) fail(sprintf("layer1 CD8 too few: %s", n8))
if (any(h$major == "CD4" & grepl("^CD8", h$subset))) fail("CD4 parent received CD8 subset")
if (any(h$major == "CD8" & grepl("^CD4|^Treg$", h$subset))) fail("CD8 parent received CD4 subset")
cd4_labs <- unique(h$subset[h$major == "CD4"])
need4 <- c("Treg", "CD4_activated", "CD4_naive", "CD4_TEM")
miss4 <- setdiff(need4, cd4_labs)
if (length(miss4)) fail(sprintf("CD4 missing %s (got %s)", paste(miss4, collapse = ","), paste(cd4_labs, collapse = ",")))
cd8_labs <- unique(h$subset[h$major == "CD8"])
need8 <- c("CD8_effector", "CD8_naive", "CD8_TEM")
miss8 <- setdiff(need8, cd8_labs)
if (length(miss8)) fail(sprintf("CD8 missing %s (got %s)", paste(miss8, collapse = ","), paste(cd8_labs, collapse = ",")))

cat("OK: P1 T/NK subset labels\n")
