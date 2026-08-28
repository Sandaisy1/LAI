# P2 Naive / Memory split (run: Rscript tests/test_flow_p2_annotate.R)
Sys.setenv(FLOW_FUNCTIONS_ONLY = "1")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/test_flow_p2_annotate.R"
source(file.path(dirname(normalizePath(this_file)), "..", "Flow_dimred_pipeline.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}

set.seed(11)
n <- 80
naive <- data.frame(
  CD19 = rnorm(n, 3.2, 0.2),
  IgD = rnorm(n, 3.1, 0.2),
  CD27 = rnorm(n, 0.3, 0.2),
  IgM = rnorm(n, 2.4, 0.2),
  IgG = rnorm(n, 0.2, 0.15),
  CD80 = rnorm(n, 0.3, 0.15),
  CD86 = rnorm(n, 0.3, 0.15),
  CD40 = rnorm(n, 0.4, 0.15),
  `BLIMP-1` = rnorm(n, 0.2, 0.15),
  check.names = FALSE
)
memory <- data.frame(
  CD19 = rnorm(n, 3.1, 0.2),
  IgD = rnorm(n, 0.3, 0.2),
  CD27 = rnorm(n, 3.0, 0.2),
  IgM = rnorm(n, 0.4, 0.2),
  IgG = rnorm(n, 2.8, 0.2),
  CD80 = rnorm(n, 0.4, 0.15),
  CD86 = rnorm(n, 0.4, 0.15),
  CD40 = rnorm(n, 0.5, 0.15),
  `BLIMP-1` = rnorm(n, 0.3, 0.15),
  check.names = FALSE
)
mat <- as.matrix(rbind(naive, memory))
h <- hierarchical_gate(mat, "P2")
if (sum(h$major == "Naive_B") < 40) fail(sprintf("expected Naive B major, got %s", paste(table(h$major), collapse=",")))
if (sum(h$major == "Memory_B") < 40) fail(sprintf("expected Memory B major, got %s", paste(table(h$major), collapse=",")))
if (length(unique(h$subset)) < 2) fail("P2 must not collapse to a single subset")

# IgD 背景都很高时，仍应按相对 CD27 分开，不能全是 Naive
hi_igd <- mat
hi_igd[, "IgD"] <- rnorm(nrow(hi_igd), 3.0, 0.15)
h2 <- hierarchical_gate(hi_igd, "P2")
if (sum(h2$major == "Memory_B") < 20) {
  fail(sprintf("high IgD background still needs a Memory B cluster, got %s", paste(names(table(h2$major)), table(h2$major), collapse=" ")))
}

cat("OK: P2 Naive/Memory split\n")
