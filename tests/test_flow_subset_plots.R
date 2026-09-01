# Per-subset bar + contour figures (run: Rscript tests/test_flow_subset_plots.R)
Sys.setenv(FLOW_FUNCTIONS_ONLY = "1")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/test_flow_subset_plots.R"
source(file.path(dirname(normalizePath(this_file)), "..", "Flow_dimred_pipeline.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}

if (!identical(p_to_star(0.0004), "***")) fail("p<0.001 should be ***")
if (!identical(p_to_star(0.02), "*")) fail("p<0.05 should be *")
if (!identical(p_to_star(0.2), "ns")) fail("p>=0.05 should be ns")
if (!grepl("P = 0.20", p_annot_label(0.2), fixed = TRUE)) fail("ns label should include the P value")

set.seed(5)
blob_x <- c(rnorm(300, 1.0, 0.18), runif(25, 3.8, 4.6))
blob_y <- c(rnorm(300, 1.0, 0.18), runif(25, 3.2, 4.0))
pc <- flow_prob_contour(blob_x, blob_y)
if (length(pc$levels) < 2) fail("probability contours should have several levels")
if (sum(pc$outlier) < 5) fail("far-away events should be outliers outside the last contour")

ct <- plot_subset_contour(
  data.frame(CD3 = blob_x, `NK1.1` = blob_y, check.names = FALSE),
  "CD3", "NK1.1", "#1A1A1A", "CD3", "NK1.1", 12.3,
  list(xmin = 0.7, xmax = 1.3, ymin = 0.7, ymax = 1.3)
)
has_contour <- FALSE
has_out_pts <- FALSE
for (ly in ct$layers) {
  if (inherits(ly$geom, "GeomContour")) has_contour <- TRUE
  if (inherits(ly$geom, "GeomPoint")) has_out_pts <- TRUE
}
if (!has_contour) fail("subset contour plot should draw FlowJo-style contours")
if (!has_out_pts) fail("subset contour plot should draw outlier events")

specs <- subset_plot_specs("P1")
nk <- Filter(function(s) identical(s$lineage, "NK"), specs)
if (!length(nk)) fail("P1 specs missing NK")
if (!identical(nk[[1]]$x, "CD3") || !identical(nk[[1]]$y, "NKp46")) fail("NK contour should be CD3 vs NKp46")
act <- Filter(function(s) identical(s$lineage, "CD4_activated"), specs)
if (!length(act) || !identical(act[[1]]$x, "CD69") || !identical(act[[1]]$y, "CD25")) {
  fail("CD4 activated contour should be CD69 vs CD25")
}

cells_na <- data.frame(
  lineage = c("Neutrophil", NA_character_, "T", "cDC1_CD103"),
  cluster_lineage = c("Myeloid", NA_character_, "T", "Myeloid"),
  stringsAsFactors = FALSE
)
pm_na <- parent_mask(cells_na, "Myeloid")
hm_na <- subset_hit_mask(cells_na, list(lineage = "Neutrophil", use_major = FALSE))
if (anyNA(pm_na) || anyNA(hm_na)) fail("NA lineage must not leak into parent/hit masks")
if (!identical(pm_na, c(TRUE, FALSE, FALSE, TRUE))) fail("parent_mask NA row should be FALSE")
if (!identical(hm_na, c(TRUE, FALSE, FALSE, FALSE))) fail("subset_hit_mask NA row should be FALSE")
ok_if <- tryCatch({
  par <- pm_na
  hit <- hm_na
  par[is.na(par)] <- FALSE
  hit[is.na(hit)] <- FALSE
  if (sum(par) < 20 || sum(par & hit) < 8) TRUE else TRUE
}, error = function(e) FALSE)
if (!isTRUE(ok_if)) fail("NA masks must not crash if (sum(par) < 20)")

set.seed(3)
mk <- function(sample, group, lineage, n, cd3, nk11) {
  data.frame(
    sample = rep(sample, n),
    group = rep(group, n),
    lineage = rep(lineage, n),
    cluster_lineage = rep(if (lineage == "NK") "NK" else "CD4", n),
    CD3 = rnorm(n, cd3, 0.25),
    `NK1.1` = rnorm(n, nk11, 0.25),
    NKp46 = rnorm(n, nk11, 0.25),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
cells <- rbind(
  mk("EV1", "EV", "NK", 80, 0.2, 3.1),
  mk("EV1", "EV", "CD4_naive", 120, 3.0, 0.3),
  mk("EV2", "EV", "NK", 70, 0.2, 3.1),
  mk("EV2", "EV", "CD4_naive", 130, 3.0, 0.3),
  mk("EV3", "EV", "NK", 90, 0.2, 3.1),
  mk("EV3", "EV", "CD4_naive", 110, 3.0, 0.3),
  mk("H1", "H", "NK", 30, 0.2, 3.1),
  mk("H1", "H", "CD4_naive", 170, 3.0, 0.3),
  mk("H2", "H", "NK", 35, 0.2, 3.1),
  mk("H2", "H", "CD4_naive", 165, 3.0, 0.3),
  mk("H3", "H", "NK", 28, 0.2, 3.1),
  mk("H3", "H", "CD4_naive", 172, 3.0, 0.3)
)

td <- tempfile("flow_subset_plots")
dir.create(td, recursive = TRUE)
export_subset_gate_figures(cells, "P1", td)
stub <- file.path(td, "subset_stats", "P1_NK_H_vs_EV")
if (!file.exists(paste0(stub, ".pdf")) || !file.exists(paste0(stub, ".png"))) {
  fail("NK subset figure pdf/png missing")
}
export_per_sample_gating_figures(cells, "P1", td)
ev1_dir <- file.path(td, "gating", "EV1")
if (!dir.exists(ev1_dir)) fail("per-sample gating folder missing for EV1")
ev1_png <- list.files(ev1_dir, pattern = "\\.png$")
if (!length(ev1_png)) fail("each sample should get complete 2D gating plots")
if (!file.exists(file.path(td, "gating", "P1_per_sample_gate_cuts.csv"))) {
  fail("per-sample gate cut table missing")
}
if (!file.exists(file.path(td, "subset_stats", "P1_subset_H_vs_EV_stats.csv"))) {
  fail("combined subset stats csv missing")
}
st <- utils::read.csv(file.path(td, "subset_stats", "P1_subset_H_vs_EV_stats.csv"), stringsAsFactors = FALSE)
if (!"NK" %in% st$subset) fail("NK not in subset stats")
if (st$mean_H[st$subset == "NK"] >= st$mean_EV[st$subset == "NK"]) {
  fail("H NK frequency should be lower than EV in this synthetic set")
}
bar <- plot_subset_stat_bar(data.frame(
  sample = c("EV1", "EV2", "EV3", "H1", "H2", "H3"),
  group = c("EV", "EV", "EV", "H", "H", "H"),
  percent = c(6.3, 7.0, 3.8, 4.1, 4.3, 3.0),
  stringsAsFactors = FALSE
), "CD4 T_EM in CD4+ (%)", 0.18)
pt_fill <- NA
for (ly in bar$layers) {
  if (inherits(ly$geom, "GeomPoint") && !is.null(ly$aes_params$fill)) pt_fill <- ly$aes_params$fill
}
if (!identical(pt_fill, "white")) fail("replicate points must be white-filled so they stay visible on the bars")

cells_na_plot <- cells
cells_na_plot$cluster_lineage[1] <- NA_character_
cells_na_plot$lineage[2] <- NA_character_
ok_na <- tryCatch({
  export_subset_gate_figures(cells_na_plot, "P1", file.path(tempdir(), "flow_subset_na"))
  TRUE
}, error = function(e) {
  cat("NA subset export crashed:", e$message, "\n")
  FALSE
})
if (!isTRUE(ok_na)) fail("NA cluster_lineage must not abort subset figures")

unlink(td, recursive = TRUE)
cat("OK: subset stat+contour figures\n")
