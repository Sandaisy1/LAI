# 全部免疫亚群功能状态图（run: Rscript tests/test_flow_functional_state.R）
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
if (!grepl("TGF", marker_pretty_label("TGF-b"), fixed = TRUE)) fail("TGF-b pretty label")

p1s <- functional_state_specs("P1")
p1_parents <- unique(vapply(p1s, function(s) s$parent, character(1)))
want_p1 <- c(
  "CD4", "CD4_naive", "CD4_TCM", "CD4_TSCM", "CD4_TEM", "CD4_TEM_early",
  "CD4_TEM_late", "CD4_SLEC", "CD4_MPEC", "Treg",
  "CD8", "CD8_naive", "CD8_TCM", "NK", "NK_immature", "NK_DP", "NK_mature",
  "NKT", "NKT_CD4", "NKT_DN"
)
if (!all(want_p1 %in% p1_parents)) {
  fail(sprintf("P1 missing parents: %s", paste(setdiff(want_p1, p1_parents), collapse = ",")))
}
circ <- c("CD4_activated", "CD4_exhausted", "CD8_effector", "NK_exhausted", "NKT_activated")
if (any(circ %in% p1_parents)) fail("P1 must not use circular activated/exhausted labels as parents")
nkt_exh <- Filter(function(s) identical(s$parent, "NKT") && identical(s$state, "exhaustion"), p1s)
if (!length(nkt_exh) || !identical(nkt_exh[[1]]$markers, c("PD-L1", "LAG-3", "TIM-3"))) {
  fail("NKT exhaustion markers")
}
p2s <- functional_state_specs("P2")
pars <- unique(vapply(p2s, function(s) s$parent, character(1)))
want <- c("Naive_B", "Unswitched_B", "Switched_B", "Atypical_B", "Activated_B")
want_p2 <- c("CD19", want, "MZ_B", "Plasmablast", "Plasma")
if (!all(want_p2 %in% pars)) {
  fail(sprintf("P2 missing parents: %s", paste(setdiff(want_p2, pars), collapse = ",")))
}
if (any(vapply(p2s, function(s) identical(s$state, "exhaustion"), logical(1)))) {
  fail("P2 must not invent B-cell exhaustion specs")
}
p3s <- functional_state_specs("P3")
if (!length(p3s)) fail("P3 myeloid subsets must get functional-state panels")
p3_parents <- unique(vapply(p3s, function(s) s$parent, character(1)))
want_p3 <- c(
  "Neutrophil", "Eosinophil", "Mast", "Macrophage", "M1_like_Mac", "M2_like_Mac",
  "DC", "cDC1_CD103", "cDC2", "Mono_Ly6Chi", "Mono_Ly6Clo"
)
if (!all(want_p3 %in% p3_parents)) {
  fail(sprintf("P3 missing parents: %s", paste(setdiff(want_p3, p3_parents), collapse = ",")))
}
if (any(c("T", "B", "NK") %in% p3_parents)) fail("P3 dump T/B/NK must not get functional-state")
p3_states <- unique(vapply(p3s, function(s) s$state, character(1)))
if (!all(c("activation", "suppression") %in% p3_states)) fail("P3 needs activation and suppression")
if ("exhaustion" %in% p3_states) fail("P3 suppression must not be named exhaustion")

if (!identical(func_x_marker("NKT", c("CD3", "NKp46")), "NKp46")) fail("NKT x-axis prefers NKp46")
if (!identical(func_x_marker("Naive_B", c("CD19", "IgD")), "CD19")) fail("B x-axis prefers CD19")
if (!identical(func_x_marker("CD4_TEM", c("CD3", "CD4")), "CD4")) fail("CD4 x-axis prefers CD4")
if (!identical(func_x_marker("CD8_naive", c("CD8", "CD3")), "CD8")) fail("CD8 x-axis prefers CD8")
if (!identical(func_x_marker("Macrophage", c("F4/80", "CD11B")), "F4/80")) fail("mac x-axis prefers F4/80")
if (!identical(func_x_marker("Neutrophil", c("LY6G", "CD11B")), "LY6G")) fail("neutrophil x-axis prefers LY6G")
if (!identical(func_x_marker("Eosinophil", c("Siglec-F", "CD11B")), "Siglec-F")) fail("eos x-axis prefers Siglec-F")
if (!identical(func_x_marker("cDC1_CD103", c("CD11C", "CD11B")), "CD11C")) fail("DC x-axis prefers CD11C")
if (!identical(func_x_marker("Mono_Ly6Chi", c("LY6C", "F4/80")), "LY6C")) fail("mono x-axis prefers LY6C")

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
if (!all(functional_parent_mask(nkt, "NKT_CD4"))) fail("NKT_CD4 lineage should count as NKT_CD4 parent")
if (any(functional_parent_mask(nkt, "CD4"))) fail("NKT_CD4 must not count as CD4 T")
if (any(functional_parent_mask(nkt, "NK"))) fail("NKT must not count as NK")

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
  mk_b("H3", "H", "H-3", "Activated_B", 50, TRUE),
  mk_b("EV1", "EV", "EV-1", "MZ_B", 40, FALSE),
  mk_b("EV2", "EV", "EV-2", "MZ_B", 40, FALSE),
  mk_b("EV3", "EV", "EV-3", "MZ_B", 40, FALSE),
  mk_b("H1", "H", "H-1", "MZ_B", 40, TRUE),
  mk_b("H2", "H", "H-2", "MZ_B", 40, TRUE),
  mk_b("H3", "H", "H-3", "MZ_B", 40, TRUE)
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
if (!file.exists(file.path(td2, "functional_state", "P2_CD19_activation_H_vs_EV.pdf"))) {
  fail("P2 CD19 (all B) activation figure missing")
}
if (!file.exists(file.path(td2, "functional_state", "P2_MZ_B_activation_H_vs_EV.pdf"))) {
  fail("P2 MZ_B activation figure missing")
}
naive_m <- functional_parent_mask(b, "Naive_B")
atyp_m <- functional_parent_mask(b, "Atypical_B")
if (any(naive_m & atyp_m)) fail("Naive_B parent must not include Atypical_B")
if (any(functional_parent_mask(b, "Naive_B") & b$lineage == "Activated_B")) {
  fail("Naive_B parent must not swallow Activated_B")
}
if (any(functional_parent_mask(b, "Naive_B") & b$lineage == "MZ_B")) {
  fail("Naive_B parent must not swallow MZ_B")
}
if (!all(functional_parent_mask(b, "CD19")[b$lineage %in% want])) {
  fail("CD19 parent should include the gated B subsets")
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

mk_t <- function(sample, group, bio, lineage, n) {
  data.frame(
    sample = rep(sample, n),
    group = rep(group, n),
    bio_sample = rep(bio, n),
    tech_rep = "1",
    lineage = lineage,
    cluster_lineage = "CD4",
    CD4 = rnorm(n, 3.0, 0.2),
    CD3 = rnorm(n, 3.1, 0.2),
    CD69 = rnorm(n, 0.8, 0.25),
    `IFN-g` = rnorm(n, 0.7, 0.25),
    `TNF-a` = rnorm(n, 0.6, 0.25),
    GZMB = rnorm(n, 0.7, 0.25),
    `PD-L1` = rnorm(n, 0.5, 0.25),
    `LAG-3` = rnorm(n, 0.5, 0.25),
    `TIM-3` = rnorm(n, 0.5, 0.25),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
cd4 <- rbind(
  mk_t("EV1", "EV", "EV-1", "CD4_naive", 40),
  mk_t("EV2", "EV", "EV-2", "CD4_naive", 40),
  mk_t("H1", "H", "H-1", "CD4_naive", 40),
  mk_t("H2", "H", "H-2", "CD4_naive", 40),
  mk_t("EV1", "EV", "EV-1", "CD4_TEM_early", 35),
  mk_t("EV2", "EV", "EV-2", "CD4_TEM_early", 35),
  mk_t("H1", "H", "H-1", "CD4_TEM_early", 35),
  mk_t("H2", "H", "H-2", "CD4_TEM_early", 35)
)
if (any(functional_parent_mask(cd4, "CD4_naive") & cd4$lineage == "CD4_TEM_early")) {
  fail("CD4_naive parent must not include T_EM early")
}
if (!all(functional_parent_mask(cd4, "CD4_TEM")[cd4$lineage == "CD4_TEM_early"])) {
  fail("CD4_TEM parent should include T_EM early")
}
if (!all(functional_parent_mask(cd4, "CD4"))) fail("CD4 parent should include naive and T_EM")

mk_mac <- function(sample, group, bio, lineage, n, hi = FALSE) {
  data.frame(
    sample = rep(sample, n),
    group = rep(group, n),
    bio_sample = rep(bio, n),
    tech_rep = "1",
    lineage = lineage,
    cluster_lineage = "Myeloid",
    `F4/80` = rnorm(n, 3.0, 0.2),
    CD11B = rnorm(n, 3.1, 0.2),
    CD86 = rnorm(n, if (hi) 3.1 else 0.7, 0.25),
    CD80 = rnorm(n, if (hi) 2.9 else 0.6, 0.25),
    CD40 = rnorm(n, if (hi) 2.8 else 1.3, 0.25),
    `TNF-a` = rnorm(n, if (hi) 2.8 else 0.5, 0.25),
    `IL-6` = rnorm(n, if (hi) 2.6 else 0.4, 0.25),
    `IL-10` = rnorm(n, if (hi) 2.7 else 0.5, 0.25),
    `TGF-b` = rnorm(n, if (hi) 2.5 else 0.4, 0.25),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
mac <- rbind(
  mk_mac("EV1", "EV", "EV-1", "Macrophage", 70, FALSE),
  mk_mac("EV2", "EV", "EV-2", "Macrophage", 70, FALSE),
  mk_mac("EV3", "EV", "EV-3", "Macrophage", 70, FALSE),
  mk_mac("H1", "H", "H-1", "Macrophage", 70, TRUE),
  mk_mac("H2", "H", "H-2", "Macrophage", 70, TRUE),
  mk_mac("H3", "H", "H-3", "Macrophage", 70, TRUE),
  mk_mac("EV1", "EV", "EV-1", "Neutrophil", 50, FALSE),
  mk_mac("H1", "H", "H-1", "Neutrophil", 50, TRUE)
)
mac$`F4/80`[mac$lineage == "Neutrophil"] <- rnorm(sum(mac$lineage == "Neutrophil"), 0.4, 0.2)
mac$LY6G <- ifelse(mac$lineage == "Neutrophil", rnorm(nrow(mac), 3.0, 0.2), rnorm(nrow(mac), 0.4, 0.2))
if (any(functional_parent_mask(mac, "Macrophage") & mac$lineage == "Neutrophil")) {
  fail("Macrophage parent must not include neutrophils")
}
tdp3 <- tempfile("flow_func_p3")
dir.create(tdp3, recursive = TRUE)
export_functional_state_figures(mac, "P3", tdp3)
mac_act <- file.path(tdp3, "functional_state", "P3_Macrophage_activation_H_vs_EV")
mac_sup <- file.path(tdp3, "functional_state", "P3_Macrophage_suppression_H_vs_EV")
if (!file.exists(paste0(mac_act, ".pdf")) || !file.exists(paste0(mac_act, ".png"))) {
  fail("P3 macrophage activation pdf/png missing")
}
if (!file.exists(paste0(mac_sup, ".pdf")) || !file.exists(paste0(mac_sup, ".png"))) {
  fail("P3 macrophage suppression pdf/png missing")
}
if (file.exists(file.path(tdp3, "functional_state", "P3_Macrophage_exhaustion_H_vs_EV.pdf"))) {
  fail("P3 must not name myeloid suppression as exhaustion")
}
st3 <- utils::read.csv(paste0(mac_act, "_stats.csv"), stringsAsFactors = FALSE)
if (!all(c("CD86", "CD80", "CD40", "TNF-a", "IL-6") %in% st3$marker)) {
  fail("P3 macrophage activation stats missing markers")
}

# 总结果已出：只读 embeddings，不重跑 panel
td3 <- tempfile("flow_func_rerun")
dir.create(file.path(td3, "P1"), recursive = TRUE)
utils::write.csv(nkt, file.path(td3, "P1", "P1_cell_embeddings.csv"), row.names = FALSE)
dir.create(file.path(td3, "P2"), recursive = TRUE)
utils::write.csv(b, file.path(td3, "P2", "P2_cell_embeddings.csv"), row.names = FALSE)
dir.create(file.path(td3, "P3"), recursive = TRUE)
utils::write.csv(mac, file.path(td3, "P3", "P3_cell_embeddings.csv"), row.names = FALSE)
export_functional_state_from_results(td3)
if (!file.exists(file.path(td3, "P1", "functional_state", "P1_NKT_activation_effector_H_vs_EV.pdf"))) {
  fail("standalone rerun should write NKT figures from embeddings")
}
if (!file.exists(file.path(td3, "P2", "functional_state", "P2_Naive_B_activation_H_vs_EV.pdf"))) {
  fail("standalone rerun should write Naive_B figures from embeddings")
}
if (!file.exists(file.path(td3, "P3", "functional_state", "P3_Macrophage_suppression_H_vs_EV.pdf"))) {
  fail("standalone rerun should write P3 macrophage figures from embeddings")
}
empty <- tempfile("flow_func_empty")
dir.create(empty)
ok_stop <- tryCatch({
  export_functional_state_from_results(empty)
  FALSE
}, error = function(e) TRUE)
if (!isTRUE(ok_stop)) fail("rerun without embeddings should stop, not silently succeed")

# 同一会话里已有旧 helper 时，独立脚本仍要从本目录覆盖最新函数
export_functional_state_from_results <- function(...) "stale"
if (exists("func_p3_state_parents", mode = "function")) {
  rm(list = "func_p3_state_parents", envir = .GlobalEnv)
}
Sys.setenv(FLOW_FUNCTIONAL_STATE_FROM_PIPELINE = "1")
source(file.path(dirname(normalizePath(this_file)), "..", "Flow_dimred_functional_state.R"))
if (!exists("func_p3_state_parents", mode = "function")) {
  fail("standalone script must reload all-subset helpers even if an old export helper exists")
}
if (!isTRUE(is.function(export_functional_state_from_results)) ||
    identical(body(export_functional_state_from_results), body(function(...) "stale"))) {
  fail("standalone script must overwrite a stale export_functional_state_from_results")
}
if (!length(functional_state_specs("P3"))) {
  fail("reloaded standalone helpers must include P3 myeloid functional-state")
}

unlink(td, recursive = TRUE)
unlink(td2, recursive = TRUE)
unlink(td3, recursive = TRUE)
unlink(tdp3, recursive = TRUE)
unlink(empty, recursive = TRUE)
cat("OK: all-subset functional-state figures\n")
