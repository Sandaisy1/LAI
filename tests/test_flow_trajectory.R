# Per-panel major-class trajectories (run: Rscript tests/test_flow_trajectory.R)
Sys.setenv(FLOW_FUNCTIONS_ONLY = "1")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/test_flow_trajectory.R"
root <- dirname(normalizePath(this_file))
source(file.path(root, "..", "Flow_dimred_pipeline.R"))
source(file.path(root, "..", "Flow_dimred_trajectory.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}

if (!identical(infer_major_lineage("P1", "Treg"), "CD4")) fail("Treg should map to CD4")
if (!identical(infer_major_lineage("P1", "CD8_effector"), "CD8")) fail("CD8 effector should map to CD8")
if (!identical(infer_major_lineage("P1", "B"), "dump")) fail("P1 B is dump")
if (!identical(infer_major_lineage("P2", "Plasma"), "B")) fail("P2 plasma is B")
if (!identical(infer_major_lineage("P3", "Neutrophil"), "Myeloid")) fail("P3 neutrophil is myeloid")
if (!identical(infer_major_lineage("P3", "T"), "dump")) fail("P3 T is dump")
if (!identical(default_trajectory_root("P1", "CD4", c("CD4_TEM", "CD4_naive")), "CD4_naive")) {
  fail("CD4 root should be naive")
}
if (!identical(default_trajectory_root("P3", "Myeloid", c("DC", "Mono_Ly6Chi")), "Mono_Ly6Chi")) {
  fail("myeloid root should be Ly6C hi mono")
}

xy <- rbind(
  cbind(rnorm(40, 0, 0.15), rnorm(40, 0, 0.15)),
  cbind(rnorm(40, 2.4, 0.15), rnorm(40, 0.1, 0.15)),
  cbind(rnorm(40, 4.6, 0.15), rnorm(40, 1.2, 0.15))
)
ed <- mst_edges_from_points(xy[c(1, 41, 81), , drop = FALSE])
if (nrow(ed) != 2) fail(sprintf("MST should have 2 edges, got %s", nrow(ed)))

set.seed(7)
blob <- function(n, mx, my) cbind(rnorm(n, mx, 0.22), rnorm(n, my, 0.22))
xy2 <- rbind(blob(90, 0, 0), blob(90, 3, 0.15), blob(90, 5.8, 1.4))
cl <- rep(c("CD4_naive", "CD4_TCM", "CD4_TEM"), each = 90)
fit <- fit_lineage_trajectory(xy2, cl, "CD4_naive")
if (is.null(fit)) fail("trajectory fit returned NULL")
if (!identical(fit$start, "CD4_naive")) fail("start cluster not naive")
if (length(fit$curves) < 1) fail("no skeleton curves")
if (length(fit$pseudotime) != nrow(xy2)) fail("pseudotime length mismatch")
pt_n <- median(fit$pseudotime[cl == "CD4_naive"], na.rm = TRUE)
pt_e <- median(fit$pseudotime[cl == "CD4_TEM"], na.rm = TRUE)
if (!is.finite(pt_n) || !is.finite(pt_e) || pt_e <= pt_n) {
  fail(sprintf("pseudotime should increase naive -> TEM (naive=%s TEM=%s)", pt_n, pt_e))
}

td <- tempfile("flow_traj")
dir.create(file.path(td, "P1"), recursive = TRUE)
df <- data.frame(
  sample = rep(rep(c("T-1", "T-2", "T-3", "T6-1", "T6-2", "T6-3"), each = 45), 1),
  group = rep(rep(c("T", "T", "T", "T6", "T6", "T6"), each = 45), 1),
  lineage = cl,
  UMAP1 = xy2[, 1],
  UMAP2 = xy2[, 2],
  cluster_lineage = "CD4",
  stringsAsFactors = FALSE
)
# pad to 80+ per sample? 270 cells total, 45 per sample. export_one_major needs n>=80 and 2 subsets. OK.
utils::write.csv(df, file.path(td, "P1", "P1_cell_embeddings.csv"), row.names = FALSE)
export_panel_trajectories(td, "P1")
need <- c(
  "TRAJECTORY_NOTE.txt",
  "P1_CD4_trajectory.pdf",
  "P1_CD4_trajectory_T_vs_T6.pdf",
  "P1_CD4_pseudotime.csv"
)
miss <- need[!file.exists(file.path(td, "P1", "trajectory", need))]
if (length(miss)) fail(sprintf("missing outputs: %s", paste(miss, collapse = ", ")))

# dump majors should not write a tree
df3 <- df
df3$lineage <- "T"
df3$cluster_lineage <- "T"
dir.create(file.path(td, "P3"), recursive = TRUE)
utils::write.csv(df3, file.path(td, "P3", "P3_cell_embeddings.csv"), row.names = FALSE)
export_panel_trajectories(td, "P3")
if (file.exists(file.path(td, "P3", "trajectory", "P3_dump_trajectory.pdf"))) {
  fail("dump lineage should not get a trajectory plot")
}

unlink(td, recursive = TRUE)
cat("OK: per-panel major-class trajectories\n")
