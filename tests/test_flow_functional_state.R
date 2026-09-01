# NKT / B 亚群功能状态图（run: Rscript tests/test_flow_functional_state.R）
Sys.setenv(FLOW_FUNCTIONS_ONLY = "1")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/test_flow_functional_state.R"
source(file.path(dirname(normalizePath(this_file)), "..", "Flow_dimred_pipeline.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}

if (!identical(flow_comparison_tag(), "H_vs_EV")) fail("original comparison tag should be H_vs_EV")
if (!identical(func_group_short("JY-EVNK"), "EVNK")) fail("JY-EVNK short label")
if (!identical(func_group_short("JZ-AB"), "AB")) fail("JZ-AB short label")
if (!grepl("IFN", marker_pretty_label("IFN-g"), fixed = TRUE)) fail("IFN-g pretty label")

p1s <- functional_state_specs("P1")
if (length(p1s) != 2L) fail("P1 should have NKT activation_effector and exhaustion")
if (!identical(p1s[[1]]$parent, "NKT") || !identical(p1s[[1]]$state, "activation_effector")) {
  fail("P1 first spec is NKT activation_effector")
}
if (!identical(p1s[[2]]$markers, c("PD-L1", "LAG-3", "TIM-3"))) fail("NKT exhaustion markers")
p2s <- functional_state_specs("P2")
pars <- vapply(p2s, function(s) s$parent, character(1))
want <- c("Naive_B", "Unswitched_B", "Switched_B", "Atypical_B", "Activated_B")
if (!identical(pars, want)) fail(sprintf("P2 parents got %s", paste(pars, collapse = ",")))
if (length(functional_state_specs("P3"))) fail("P3 should not get NKT/B functional-state panels")

if (!identical(func_x_marker("NKT", c("CD3", "NKp46")), "NKp46")) fail("NKT x-axis prefers NKp46")
if (!identical(func_x_marker("Naive_B", c("CD19", "IgD")), "CD19")) fail("B x-axis prefers CD19")

set.seed(11)
mk_nkt <- function(sample, group, bio, n, hi = FALSE) {
  data.frame(
    sample = rep(sample, n),
    group = rep(group, n),
    bio_sample = rep(bio, n),
    tech_rep = "1",
    lineage = "NKT_CD4",
    cluster_lineage = "NKT",
    NKp46 = rnorm(n, 3.0, 0.2),
    CD3 = rnorm(n, 3.1, 0.2),
    CD69 = rnorm(n, if (hi) 3.2 else 0.8, 0.25),
    `IFN-g` = rnorm(n, if (hi) 3.0 else 0.7, 0.25),
    `TNF-a` = rnorm(n, if (hi) 2.9 else 0.6, 0.25),
    GZMB = rnorm(n, if (hi) 3.1 else 0.7, 0.25),
    `PD-L1` = rnorm(n, if (hi) 2.8 else 0.5, 0.25),
    `LAG-3` = rnorm(n, if (hi) 2.7 else 0.5, 0.25),
    `TIM-3` = rnorm(n, if (hi) 2.6 else 0.5, 0.25),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
mk_b <- function(sample, group, bio, lineage, n, hi = FALSE) {
  data.frame(
    sample = rep(sample, n),
    group = rep(group, n),
    bio_sample = rep(bio, n),
    tech_rep = "1",
    lineage = lineage,
    cluster_lineage = lineage,
    CD19 = rnorm(n, 3.0, 0.2),
    CD86 = rnorm(n, if (hi) 3.1 else 0.7, 0.25),
    CD80 = rnorm(n, if (hi) 2.9 else 0.6, 0.25),
    CD40 = rnorm(n, if (hi) 2.8 else 1.4, 0.25),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

nkt <- rbind(
  mk_nkt("EV1", "EV", "EV-1", 80, FALSE),
  mk_nkt("EV2", "EV", "EV-2", 80, FALSE),
  mk_nkt("EV3", "EV", "EV-3", 80, FALSE),
  mk_nkt("H1", "H", "H-1", 80, TRUE),
  mk_nkt("H2", "H", "H-2", 80, TRUE),
  mk_nkt("H3", "H", "H-3", 80, TRUE)
)
pm <- functional_parent_mask(nkt, "NKT")
if (anyNA(pm) || !all(pm)) fail("all synthetic NKT rows should be in the NKT parent")
if (any(functional_parent_mask(nkt, "Naive_B"))) fail("NKT cells must not count as Naive_B")

td <- tempfile("flow_func_state")
dir.create(td, recursive = TRUE)
export_functional_state_figures(nkt, "P1", td)
act <- file.path(td, "functional_state", "P1_NKT_activation_effector_H_vs_EV")
exh <- file.path(td, "functional_state", "P1_NKT_exhaustion_H_vs_EV")
if (!file.exists(paste0(act, ".pdf")) || !file.exists(paste0(act, ".png"))) {
  fail("NKT activation_effector pdf/png missing")
}
if (!file.exists(paste0(exh, ".pdf")) || !file.exists(paste0(exh, ".png"))) {
  fail("NKT exhaustion pdf/png missing")
}
st <- utils::read.csv(paste0(act, "_stats.csv"), stringsAsFactors = FALSE)
if (!all(c("CD69", "IFN-g", "TNF-a", "GZMB") %in% st$marker)) fail("NKT activation stats missing markers")
if (!identical(as.integer(st$n_ctrl[1]), 2L) || !identical(as.integer(st$n_trt[1]), 2L)) {
  fail("functional-state bars should use n=2 after dropping 1 extreme bio-rep")
}
if (st$mean_trt[st$marker == "IFN-g"] <= st$mean_ctrl[st$marker == "IFN-g"]) {
  fail("H NKT IFN-g+ should be higher than EV in this synthetic set")
}

b <- rbind(
  mk_b("EV1", "EV", "EV-1", "Naive_B", 70, FALSE),
  mk_b("EV2", "EV", "EV-2", "Naive_B", 70, FALSE),
  mk_b("EV3", "EV", "EV-3", "Naive_B", 70, FALSE),
  mk_b("H1", "H", "H-1", "Naive_B", 70, TRUE),
  mk_b("H2", "H", "H-2", "Naive_B", 70, TRUE),
  mk_b("H3", "H", "H-3", "Naive_B", 70, TRUE),
  mk_b("EV1", "EV", "EV-1", "Atypical_B", 60, FALSE),
  mk_b("EV2", "EV", "EV-2", "Atypical_B", 60, FALSE),
  mk_b("EV3", "EV", "EV-3", "Atypical_B", 60, FALSE),
  mk_b("H1", "H", "H-1", "Atypical_B", 60, TRUE),
  mk_b("H2", "H", "H-2", "Atypical_B", 60, TRUE),
  mk_b("H3", "H", "H-3", "Atypical_B", 60, TRUE),
  mk_b("EV1", "EV", "EV-1", "Unswitched_B", 55, FALSE),
  mk_b("EV2", "EV", "EV-2", "Unswitched_B", 55, FALSE),
  mk_b("EV3", "EV", "EV-3", "Unswitched_B", 55, FALSE),
  mk_b("H1", "H", "H-1", "Unswitched_B", 55, TRUE),
  mk_b("H2", "H", "H-2", "Unswitched_B", 55, TRUE),
  mk_b("H3", "H", "H-3", "Unswitched_B", 55, TRUE),
  mk_b("EV1", "EV", "EV-1", "Switched_B", 55, FALSE),
  mk_b("EV2", "EV", "EV-2", "Switched_B", 55, FALSE),
  mk_b("EV3", "EV", "EV-3", "Switched_B", 55, FALSE),
  mk_b("H1", "H", "H-1", "Switched_B", 55, TRUE),
  mk_b("H2", "H", "H-2", "Switched_B", 55, TRUE),
  mk_b("H3", "H", "H-3", "Switched_B", 55, TRUE),
  mk_b("EV1", "EV", "EV-1", "Activated_B", 50, TRUE),
  mk_b("EV2", "EV", "EV-2", "Activated_B", 50, TRUE),
  mk_b("EV3", "EV", "EV-3", "Activated_B", 50, TRUE),
  mk_b("H1", "H", "H-1", "Activated_B", 50, TRUE),
  mk_b("H2", "H", "H-2", "Activated_B", 50, TRUE),
  mk_b("H3", "H", "H-3", "Activated_B", 50, TRUE)
)
td2 <- tempfile("flow_func_b")
dir.create(td2, recursive = TRUE)
export_functional_state_figures(b, "P2", td2)
for (p in want) {
  stub <- file.path(td2, "functional_state", paste0("P2_", p, "_activation_H_vs_EV"))
  if (!file.exists(paste0(stub, ".pdf")) || !file.exists(paste0(stub, ".png"))) {
    fail(sprintf("%s activation figure missing", p))
  }
}
if (file.exists(file.path(td2, "functional_state", "P2_Naive_B_exhaustion_H_vs_EV.pdf"))) {
  fail("P2 must not invent B-cell exhaustion figures")
}
naive_m <- functional_parent_mask(b, "Naive_B")
atyp_m <- functional_parent_mask(b, "Atypical_B")
if (any(naive_m & atyp_m)) fail("Naive_B parent must not include Atypical_B")
if (any(functional_parent_mask(b, "Naive_B") & b$lineage == "Activated_B")) {
  fail("Naive_B parent must not swallow Activated_B")
}

bar <- plot_func_state_bar(data.frame(
  sample = c("EV1", "EV2", "H1", "H2"),
  group = c("EV", "EV", "H", "H"),
  percent = c(8, 10, 22, 24),
  stringsAsFactors = FALSE
), "IFN-\u03B3\u207A NKT cell (%)", 0.04)
pt_fill <- NA
col_fill <- NA
for (ly in bar$layers) {
  if (inherits(ly$geom, "GeomCol") && !is.null(ly$aes_params$fill)) col_fill <- ly$aes_params$fill
  if (inherits(ly$geom, "GeomPoint") && !is.null(ly$aes_params$shape)) pt_fill <- ly$aes_params$shape
}
if (!identical(col_fill, "white")) fail("functional-state bars should be open (white fill)")
if (!identical(as.integer(pt_fill), 21L)) fail("replicate points should be filled circles (shape 21)")

# 总结果已出：只读 embeddings，不重跑 panel
td3 <- tempfile("flow_func_rerun")
dir.create(file.path(td3, "P1"), recursive = TRUE)
utils::write.csv(nkt, file.path(td3, "P1", "P1_cell_embeddings.csv"), row.names = FALSE)
dir.create(file.path(td3, "P2"), recursive = TRUE)
utils::write.csv(b, file.path(td3, "P2", "P2_cell_embeddings.csv"), row.names = FALSE)
export_functional_state_from_results(td3)
if (!file.exists(file.path(td3, "P1", "functional_state", "P1_NKT_activation_effector_H_vs_EV.pdf"))) {
  fail("standalone rerun should write NKT figures from embeddings")
}
if (!file.exists(file.path(td3, "P2", "functional_state", "P2_Naive_B_activation_H_vs_EV.pdf"))) {
  fail("standalone rerun should write Naive_B figures from embeddings")
}
empty <- tempfile("flow_func_empty")
dir.create(empty)
ok_stop <- tryCatch({
  export_functional_state_from_results(empty)
  FALSE
}, error = function(e) TRUE)
if (!isTRUE(ok_stop)) fail("rerun without embeddings should stop, not silently succeed")

unlink(td, recursive = TRUE)
unlink(td2, recursive = TRUE)
unlink(td3, recursive = TRUE)
unlink(empty, recursive = TRUE)
cat("OK: NKT/B functional-state figures\n")
