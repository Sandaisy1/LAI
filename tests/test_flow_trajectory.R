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
if (!identical(infer_major_lineage("P1", "CD4_effector"), "CD4")) fail("CD4 effector should map to CD4")
if (!identical(infer_major_lineage("P1", "CD8_effector"), "CD8")) fail("CD8 effector should map to CD8")
if (!identical(infer_major_lineage("P1", "CD8_activated"), "CD8")) fail("CD8 activated should map to CD8")
if (!identical(infer_major_lineage("P1", "NK_effector"), "NK")) fail("NK effector should map to NK")
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

d_na <- matrix(c(0, NA, NA, 0), 2, 2)
ed_na <- tryCatch(mst_edges_from_points(d_na), error = function(e) e)
if (inherits(ed_na, "error")) fail(sprintf("NA distances crashed MST: %s", ed_na$message))
if (!is.matrix(ed_na) || nrow(ed_na) != 0) fail("NA distances should yield no MST edges")

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

df_plot <- data.frame(lineage = cl, stringsAsFactors = FALSE)
p_tr <- plot_trajectory_tree(df_plot, fit, "P1", "size check", facet_group = FALSE)
pt_size <- NA_real_
path_lw <- NA_real_
for (ly in p_tr$layers) {
  if (inherits(ly$geom, "GeomPoint") && !is.null(ly$aes_params$size)) pt_size <- ly$aes_params$size
  if (inherits(ly$geom, "GeomPath") && !is.null(ly$aes_params$linewidth)) path_lw <- ly$aes_params$linewidth
}
if (!is.finite(pt_size) || pt_size < 1.1) fail(sprintf("trajectory points too small: %s", pt_size))
if (!is.finite(path_lw) || path_lw > 0.45) fail(sprintf("trajectory skeleton too thick: %s", path_lw))
if (!(path_lw < 0.45 * pt_size)) fail("skeleton must stay thinner than the cell points")

td <- tempfile("flow_traj")
dir.create(file.path(td, "P1"), recursive = TRUE)
df <- data.frame(
  sample = rep(rep(c("EV1", "EV2", "EV3", "H1", "H2", "H3"), each = 45), 1),
  group = rep(rep(c("EV", "EV", "EV", "H", "H", "H"), each = 45), 1),
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
  "P1_CD4_trajectory_H_vs_EV.pdf",
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

# NK is one population: skip is expected, not a crash.
df_nk <- data.frame(
  sample = rep(c("EV1", "EV2", "EV3", "H1", "H2", "H3"), each = 90),
  group = rep(c("EV", "EV", "EV", "H", "H", "H"), each = 90),
  lineage = "NK",
  UMAP1 = rnorm(540),
  UMAP2 = rnorm(540),
  cluster_lineage = "NK",
  stringsAsFactors = FALSE
)
utils::write.csv(df_nk, file.path(td, "P1", "P1_cell_embeddings.csv"), row.names = FALSE)
skip_msg <- paste(capture.output(export_panel_trajectories(td, "P1")), collapse = "\n")
if (!grepl("skip trajectory", skip_msg)) fail("NK should log skip trajectory")
if (!grepl("only 1 subset", skip_msg)) fail("NK skip should say only 1 subset")
if (file.exists(file.path(td, "P1", "trajectory", "P1_NK_trajectory.pdf"))) {
  fail("single-subset NK should not write a trajectory tree")
}

unlink(td, recursive = TRUE)
cat("OK: per-panel major-class trajectories\n")
