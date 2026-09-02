# P2 IgD vs CD27 Naive / Unswitched / Switched (run: Rscript tests/test_flow_p2_annotate.R)
Sys.setenv(FLOW_FUNCTIONS_ONLY = "1")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/test_flow_p2_annotate.R"
source(file.path(dirname(normalizePath(this_file)), "..", "Flow_dimred_pipeline.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}

expect_label <- function(got, want, tag) {
  if (!identical(got, want)) fail(sprintf("%s: got %s, want %s", tag, got, want))
}

expect_label(celltype_label("Unswitched_B", "P2"), "Unswitched memory B", "label-unswitched")
expect_label(celltype_label("Switched_B", "P2"), "Switched memory B", "label-switched")
expect_label(celltype_label("MZ_B", "P2"), "MZ B", "label-mz")
expect_label(celltype_label("Plasmablast", "P2"), "Plasmablast", "label-pb")
expect_label(celltype_label("IgM_memory", "P2"), "Unswitched memory B", "label-igm-alias")

p2_specs <- subset_plot_specs("P2")
naive_sp <- Filter(function(s) identical(s$lineage, "Naive_B"), p2_specs)
if (!length(naive_sp) || !identical(naive_sp[[1]]$x, "IgD") || !identical(naive_sp[[1]]$y, "CD27")) {
  fail("Naive B contour should be IgD vs CD27")
}
unsw <- Filter(function(s) identical(s$lineage, "Unswitched_B"), p2_specs)
if (!length(unsw)) fail("P2 specs missing Unswitched_B")
if (!identical(unsw[[1]]$parent, "CD19")) fail("Unswitched parent should be CD19+")

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
switched <- data.frame(
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
mat <- as.matrix(rbind(naive, switched))
h <- hierarchical_gate(mat, "P2")
if (sum(h$major == "Naive_B") < 40) {
  fail(sprintf("expected Naive B major, got %s", paste(names(table(h$major)), table(h$major), collapse = " ")))
}
if (sum(h$major == "Switched_B") < 40) {
  fail(sprintf("expected Switched B major (IgD- CD27+), got %s", paste(names(table(h$major)), table(h$major), collapse = " ")))
}
if (any(h$major %in% c("Plasma"))) fail("CD19+ IgD/CD27 split must not call Plasma")
if (length(unique(h$subset)) < 2) fail("P2 must not collapse to a single subset")

# IgD 背景都很高时，仍应按 CD27 分开 Naive vs Unswitched，不能全是 Naive
hi_igd <- mat
hi_igd[, "IgD"] <- rnorm(nrow(hi_igd), 3.0, 0.15)
h2 <- hierarchical_gate(hi_igd, "P2")
if (sum(h2$major == "Unswitched_B") < 20) {
  fail(sprintf("high IgD background still needs Unswitched (CD27+), got %s",
               paste(names(table(h2$major)), table(h2$major), collapse = " ")))
}

# CD40 在静息 Naive 上也高，不能把整团打成 Activated
naive_cd40 <- naive
naive_cd40$CD40 <- rnorm(n, 3.2, 0.15)
naive_cd40$CD80 <- rnorm(n, 0.35, 0.12)
naive_cd40$CD86 <- rnorm(n, 0.35, 0.12)
mat_cd40 <- as.matrix(rbind(naive_cd40, switched))
h3 <- hierarchical_gate(mat_cd40, "P2")
if (mean(h3$subset[seq_len(n)] == "Activated_B") > 0.2) {
  fail(sprintf("CD40-high naive must stay Naive, got %s",
               paste(names(table(h3$subset[seq_len(n)])), table(h3$subset[seq_len(n)]), collapse = " ")))
}
if (sum(h3$subset == "Naive_B") < 40) {
  fail(sprintf("CD40 must not erase Naive B, got %s", paste(names(table(h3$subset)), table(h3$subset), collapse = " ")))
}

# CD19- BLIMP+ CD27++ IgD- → Plasma 大类；CD19+ 的 BLIMP 不得在第 1 层叫 Plasma
set.seed(13)
plasma <- data.frame(
  CD19 = rnorm(60, 0.3, 0.15),
  IgD = rnorm(60, 0.25, 0.12),
  CD27 = rnorm(60, 3.1, 0.2),
  IgM = rnorm(60, 0.3, 0.12),
  IgG = rnorm(60, 0.3, 0.12),
  CD80 = rnorm(60, 0.3, 0.12),
  CD86 = rnorm(60, 0.3, 0.12),
  CD40 = rnorm(60, 0.4, 0.12),
  `BLIMP-1` = rnorm(60, 3.2, 0.15),
  check.names = FALSE
)
h_pl <- hierarchical_gate(as.matrix(rbind(naive, switched, plasma)), "P2")
if (sum(h_pl$major == "Plasma") < 30) {
  fail(sprintf("CD19- BLIMP+ should be Plasma major, got %s",
               paste(names(table(h_pl$major)), table(h_pl$major), collapse = " ")))
}
if (any(h_pl$major[seq_len(2 * n)] == "Plasma")) {
  fail("CD19+ cells must not be Plasma at layer 1")
}

# Switched 内 BLIMP 少数岛 → Plasmablast
set.seed(14)
pb <- switched
pb$`BLIMP-1` <- rnorm(n, 3.1, 0.15)
pb$IgG <- rnorm(n, 0.4, 0.12)
h_pb <- hierarchical_gate(as.matrix(rbind(naive, switched, pb)), "P2")
if (sum(h_pb$subset == "Plasmablast") < 20) {
  fail(sprintf("switched BLIMP+ should be Plasmablast, got %s",
               paste(names(table(h_pb$subset)), table(h_pb$subset), collapse = " ")))
}

# Naive 内 IgM 高 → MZ；真正的 CD80/CD86 少数岛可以叫 Activated
set.seed(12)
mz <- naive
mz$IgM <- rnorm(40, 3.4, 0.12)
mz <- mz[seq_len(40), ]
act_island <- data.frame(
  CD19 = rnorm(40, 3.0, 0.2),
  IgD = rnorm(40, 0.4, 0.15),
  CD27 = rnorm(40, 2.6, 0.2),
  IgM = rnorm(40, 0.3, 0.15),
  IgG = rnorm(40, 0.3, 0.15),
  CD80 = rnorm(40, 2.9, 0.15),
  CD86 = rnorm(40, 3.0, 0.15),
  CD40 = rnorm(40, 2.6, 0.15),
  `BLIMP-1` = rnorm(40, 0.25, 0.1),
  check.names = FALSE
)
unsw <- data.frame(
  CD19 = rnorm(n, 3.1, 0.2),
  IgD = rnorm(n, 3.0, 0.2),
  CD27 = rnorm(n, 3.0, 0.2),
  IgM = rnorm(n, 2.6, 0.2),
  IgG = rnorm(n, 0.2, 0.12),
  CD80 = rnorm(n, 0.35, 0.12),
  CD86 = rnorm(n, 0.35, 0.12),
  CD40 = rnorm(n, 0.5, 0.15),
  `BLIMP-1` = rnorm(n, 0.25, 0.1),
  check.names = FALSE
)
mat4 <- as.matrix(rbind(naive, switched, unsw, mz, act_island))
h4 <- hierarchical_gate(mat4, "P2")
tab4 <- table(h4$subset)
if (is.na(tab4["Naive_B"]) || tab4["Naive_B"] < 40) fail("mixed P2 lost Naive B")
if (is.na(tab4["Switched_B"])) fail("mixed P2 lost Switched memory B")
if (is.na(tab4["Unswitched_B"])) fail("mixed P2 lost Unswitched memory B")
if (mean(h4$subset == "Activated_B") > 0.35) {
  fail(sprintf("Activated B swallowed P2: %s", paste(names(tab4), tab4, collapse = " ")))
}
if (length(unique(h4$subset)) < 4) {
  fail(sprintf("P2 should keep several subsets, got %s", paste(names(tab4), tab4, collapse = " ")))
}

# 活化读出是 MFI，不是新的第 1 层亚群
set.seed(15)
cells_act <- data.frame(
  sample = rep(c("EV1", "EV2", "EV3", "H1", "H2", "H3"), each = 40),
  bio_sample = rep(c("EV1", "EV2", "EV3", "H1", "H2", "H3"), each = 40),
  group = rep(c("EV", "EV", "EV", "H", "H", "H"), each = 40),
  cluster_lineage = "Naive_B",
  lineage = "Naive_B",
  CD40 = rnorm(240, 1.2, 0.2),
  CD80 = rnorm(240, 0.4, 0.15),
  CD86 = rnorm(240, 0.5, 0.15),
  stringsAsFactors = FALSE
)
td <- tempfile("p2_act")
dir.create(td)
export_p2_activation_stats(cells_act, td)
if (!file.exists(file.path(td, "P2_Bcell_activation_by_sample.csv"))) {
  fail("P2 activation MFI table missing")
}
unlink(td, recursive = TRUE)

cat("OK: P2 IgD/CD27 Naive-Unswitched-Switched split\n")
