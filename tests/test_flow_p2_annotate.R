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

# CD40 在静息 Naive 上也高，不能把整团打成 Activated
naive_cd40 <- naive
naive_cd40$CD40 <- rnorm(n, 3.2, 0.15)
naive_cd40$CD80 <- rnorm(n, 0.35, 0.12)
naive_cd40$CD86 <- rnorm(n, 0.35, 0.12)
mat_cd40 <- as.matrix(rbind(naive_cd40, memory))
h3 <- hierarchical_gate(mat_cd40, "P2")
if (mean(h3$subset[seq_len(n)] == "Activated_B") > 0.2) {
  fail(sprintf("CD40-high naive must stay Naive, got %s", paste(names(table(h3$subset[seq_len(n)])), table(h3$subset[seq_len(n)]), collapse=" ")))
}
if (sum(h3$subset == "Naive_B") < 40) {
  fail(sprintf("CD40 must not erase Naive B, got %s", paste(names(table(h3$subset)), table(h3$subset), collapse=" ")))
}

# 真正的 CD80/CD86 少数岛可以叫 Activated，但不能占满 Naive/Memory
set.seed(12)
act_island <- data.frame(
  CD19 = rnorm(40, 3.0, 0.2),
  IgD = rnorm(40, 0.4, 0.15),
  CD27 = rnorm(40, 2.2, 0.2),
  IgM = rnorm(40, 0.3, 0.15),
  IgG = rnorm(40, 0.3, 0.15),
  CD80 = rnorm(40, 2.9, 0.15),
  CD86 = rnorm(40, 3.0, 0.15),
  CD40 = rnorm(40, 2.6, 0.15),
  `BLIMP-1` = rnorm(40, 0.25, 0.1),
  check.names = FALSE
)
switched <- memory
switched$IgG <- rnorm(n, 3.0, 0.15)
switched$IgM <- rnorm(n, 0.3, 0.12)
igm_mem <- memory
igm_mem$IgM <- rnorm(n, 2.9, 0.15)
igm_mem$IgG <- rnorm(n, 0.25, 0.1)
mat4 <- as.matrix(rbind(naive, memory, switched, igm_mem, act_island))
h4 <- hierarchical_gate(mat4, "P2")
tab4 <- table(h4$subset)
if (is.na(tab4["Naive_B"]) || tab4["Naive_B"] < 40) fail("mixed P2 lost Naive B")
if (is.na(tab4["Memory_B"]) && is.na(tab4["Switched_B"])) fail("mixed P2 lost Memory/Switched")
if (mean(h4$subset == "Activated_B") > 0.35) {
  fail(sprintf("Activated B swallowed P2: %s", paste(names(tab4), tab4, collapse=" ")))
}
if (length(unique(h4$subset)) < 4) {
  fail(sprintf("P2 should keep several subsets, got %s", paste(names(tab4), tab4, collapse=" ")))
}

cat("OK: P2 Naive/Memory split\n")
