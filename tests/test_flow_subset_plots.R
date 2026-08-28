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

specs <- subset_plot_specs("P1")
nk <- Filter(function(s) identical(s$lineage, "NK"), specs)
if (!length(nk)) fail("P1 specs missing NK")
if (!identical(nk[[1]]$x, "CD3") || !identical(nk[[1]]$y, "NK1.1")) fail("NK contour should be CD3 vs NK1.1")

set.seed(3)
mk <- function(sample, group, lineage, n, cd3, nk11) {
  data.frame(
    sample = rep(sample, n),
    group = rep(group, n),
    lineage = rep(lineage, n),
    cluster_lineage = rep(if (lineage == "NK") "NK" else "CD4", n),
    CD3 = rnorm(n, cd3, 0.25),
    `NK1.1` = rnorm(n, nk11, 0.25),
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
if (!file.exists(file.path(td, "subset_stats", "P1_subset_H_vs_EV_stats.csv"))) {
  fail("combined subset stats csv missing")
}
st <- utils::read.csv(file.path(td, "subset_stats", "P1_subset_H_vs_EV_stats.csv"), stringsAsFactors = FALSE)
if (!"NK" %in% st$subset) fail("NK not in subset stats")
if (st$mean_H[st$subset == "NK"] >= st$mean_EV[st$subset == "NK"]) {
  fail("H NK frequency should be lower than EV in this synthetic set")
}

unlink(td, recursive = TRUE)
cat("OK: subset stat+contour figures\n")
