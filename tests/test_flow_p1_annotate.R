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
expect(celltype_label("CD4_TCM", "P1"), "CD4 T_CM", "label-cd4-tcm")
expect(celltype_label("CD4_TEM", "P1"), "CD4 T_EM", "label-cd4-tem")
expect(celltype_label("CD4_effector", "P1"), "CD4 T_EFF", "label-cd4-teff")
expect(celltype_label("CD8_TCM", "P1"), "CD8 T_CM", "label-cd8-tcm")
expect(celltype_label("CD8_TEM", "P1"), "CD8 T_EM", "label-cd8-tem")
expect(celltype_label("CD8_effector", "P1"), "CD8 T_EFF", "label-cd8-teff")
expect(lab_of(mk(CD3 = 3.1, CD4 = 3.0, CD25 = 3.2, CD69 = 0.4, CD44 = 1.8)), "Treg", "treg")
expect(lab_of(mk(CD3 = 3.2, CD4 = 3.0, CD69 = 3.1, CD44 = 2.6, CD62L = 0.4, `TNF-a` = 2.4)), "CD4_activated", "cd4-act")
expect(
  lab_of(mk(CD3 = 3.2, CD4 = 3.0, CD8 = 0.1, CD62L = 3.2, CD44 = 0.3, CD69 = 3.1, CD25 = 2.8)),
  "CD4_activated", "cd4-act-naive-phenotype"
)
expect(
  lab_of(mk(CD3 = 3.2, CD4 = 3.0, CD8 = 0.1, CD62L = 0.3, CD44 = 3.0, `IFN-g` = 3.2, `TNF-a` = 3.0, CD69 = 0.4)),
  "CD4_effector", "cd4-eff"
)
expect(lab_of(mk(CD3 = 3.3, CD8 = 3.1, CD8b = 2.9, CD4 = 0.1, CD62L = 3.1, CD44 = 0.3)), "CD8_naive", "cd8-naive")
expect(
  lab_of(mk(CD3 = 3.3, CD8 = 3.1, CD4 = 0.1, CD62L = 3.1, CD44 = 0.3, CD69 = 3.1, CD25 = 2.8)),
  "CD8_activated", "cd8-act-naive-phenotype"
)
expect(lab_of(mk(CD3 = 3.3, CD8 = 3.1, CD8b = 2.8, CD62L = 2.7, CD44 = 2.6)), "CD8_TCM", "cd8-tcm")
expect(lab_of(mk(CD3 = 3.3, CD8 = 3.1, CD8b = 2.8, CD44 = 3.0, CD62L = 0.3)), "CD8_TEM", "cd8-tem")
expect(lab_of(mk(CD3 = 3.3, CD8 = 3.1, CD44 = 2.6, CD62L = 0.4, GZMB = 3.0, Perforin = 2.8)), "CD8_effector", "cd8-eff")
expect(lab_of(mk(CD3 = 3.0, CD8 = 3.0, `LAG-3` = 3.1, `TIM-3` = 2.9, CD44 = 2.8)), "CD8_exhausted", "cd8-exh")
expect(lab_of(mk(CD3 = 0.1, `NK1.1` = 3.2, NKp46 = 3.0, CD8 = 1.0)), "NK", "nk")
expect(lab_of(mk(CD3 = 0.1, NKp46 = 3.1, `NK1.1` = 0.3, CD8 = 0.2)), "NK", "nk-by-nkp46")
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
  bind_pop("Treg", 80), bind_pop("CD4_act", 80), bind_pop("CD4_eff", 80),
  bind_pop("CD8_naive", 80), bind_pop("CD8_act", 80), bind_pop("CD8_TCM", 80),
  bind_pop("CD8_TEM", 80), bind_pop("CD8_eff", 80), bind_pop("CD8_exh", 80)
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
need4 <- c("Treg", "CD4_activated", "CD4_effector", "CD4_naive", "CD4_TCM", "CD4_TEM")
miss4 <- setdiff(need4, cd4_labs)
if (length(miss4)) fail(sprintf("CD4 missing %s (got %s)", paste(miss4, collapse = ","), paste(cd4_labs, collapse = ",")))
cd8_labs <- unique(h$subset[h$major == "CD8"])
need8 <- c("CD8_effector", "CD8_activated", "CD8_naive", "CD8_TCM", "CD8_TEM")
miss8 <- setdiff(need8, cd8_labs)
if (length(miss8)) fail(sprintf("CD8 missing %s (got %s)", paste(miss8, collapse = ","), paste(cd8_labs, collapse = ",")))

set.seed(4)
n <- 80
df_mix <- data.frame(
  group = rep(c("EV", "H"), each = n),
  cluster_lineage = rep(rep(c("CD4", "NK"), each = n / 2), 2),
  lineage = rep(rep(c("CD4_naive", "NK"), each = n / 2), 2),
  UMAP1 = c(rnorm(n, 0, 0.2), rnorm(n, 0.05, 0.2)),
  UMAP2 = c(rnorm(n, 0, 0.2), rnorm(n, 0.05, 0.2)),
  CD3 = c(rnorm(n / 2, 3, 0.1), rnorm(n / 2, 0.2, 0.1), rnorm(n / 2, 3, 0.1), rnorm(n / 2, 0.2, 0.1)),
  `NK1.1` = c(rnorm(n / 2, 0.2, 0.1), rnorm(n / 2, 3, 0.1), rnorm(n / 2, 0.2, 0.1), rnorm(n / 2, 3, 0.1)),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
p_sep <- plot_split_lineage(df_mix, "UMAP1", "UMAP2", "P1", "UMAP-1", "UMAP-2", "sep")
if (!("UMAP1" %in% names(p_sep$data))) fail("main split plot should use the shared embedding")
x_lab <- paste(deparse(p_sep$mapping$x), collapse = "")
if (grepl("sep1", x_lab, fixed = TRUE)) {
  fail("main tSNE/UMAP must not tile each subset into a separate box")
}
has_poly <- FALSE
has_pt <- FALSE
for (ly in p_sep$layers) {
  if (inherits(ly$geom, "GeomPolygon")) has_poly <- TRUE
  if (inherits(ly$geom, "GeomPoint")) has_pt <- TRUE
}
if (!isTRUE(has_pt)) fail("figure-1 tSNE should draw colored cells on a shared embedding")
if (isTRUE(has_poly)) fail("figure-1 tSNE should not split each subset into its own filled tile")
built <- ggplot2::ggplot_build(p_sep)
lay <- built$layout$layout
if (nrow(lay) != 2L) fail("figure-1 must be two panels (EV | H), not a per-subset grid")
if (!("group" %in% names(lay))) fail("figure-1 facets by experimental group")
if ("celltype" %in% names(lay) && length(unique(lay$celltype)) > 1L) {
  fail("figure-1 must not facet by cell type (that is figure 2)")
}
if ("lineage" %in% names(lay) && length(unique(lay$lineage)) > 1L) {
  fail("figure-1 must not facet by lineage")
}

set.seed(8)
tail_v <- c(rnorm(220, 0.45, 0.12), rnorm(50, 3.1, 0.2))
cut_t <- axis_pos_cut(tail_v)
if (!is.finite(cut_t) || cut_t < 1.0 || cut_t > 2.6) {
  fail(sprintf("CD62L-like tail cut should sit between the origin blob and the tail, got %s", cut_t))
}

set.seed(9)
n_q <- 60
mat_q <- cbind(
  CD62L = c(rnorm(n_q, 3.2, 0.15), rnorm(n_q, 2.9, 0.15), rnorm(n_q, 0.35, 0.12)),
  CD44 = c(rnorm(n_q, 0.35, 0.12), rnorm(n_q, 2.9, 0.15), rnorm(n_q, 3.1, 0.15))
)
labs_q <- split_memory_3(mat_q, seq_len(nrow(mat_q)), "CD4")
if (sum(labs_q[seq_len(n_q)] == "CD4_naive") < 40) fail("quadrant must call CD62L+ CD44- naive")
if (sum(labs_q[seq_len(n_q) + n_q] == "CD4_TCM") < 40) fail("quadrant must call CD62L+ CD44+ TCM")
if (sum(labs_q[seq_len(n_q) + 2 * n_q] == "CD4_TEM") < 40) fail("quadrant must call CD62L- CD44+ TEM")

set.seed(13)
mat_dn <- cbind(
  CD62L = rnorm(90, 0.32, 0.1),
  CD44 = rnorm(90, 0.28, 0.1)
)
labs_dn <- split_memory_3(mat_dn, seq_len(nrow(mat_dn)), "CD4")
if (mean(labs_dn == "CD4_TEM") < 0.85) {
  fail("CD62L- CD44- double-negative cells must be assigned TEM, not left ungated")
}

set.seed(21)
smear <- c(rnorm(420, 0.48, 0.16), 1.15 + stats::rexp(55, rate = 0.9))
cut_sm <- axis_pos_cut(smear)
if (!is.finite(cut_sm) || cut_sm < 0.85 || cut_sm > 2.5) {
  fail(sprintf("blob+smear cut should sit at the right edge of the origin mass, got %s", cut_sm))
}
if (mean(smear >= cut_sm) > 0.35) fail("smear cut must not slice the origin blob in half")

set.seed(10)
mat_s <- rbind(
  cbind(CD3 = rnorm(50, 3.1, 0.12), CD4 = rnorm(50, 3.0, 0.12), CD8 = rnorm(50, 0.2, 0.1),
        CD62L = rnorm(50, 3.2, 0.12), CD44 = rnorm(50, 0.3, 0.1), CD19 = 0.2, CD11B = 0.2, `NK1.1` = 0.2),
  cbind(CD3 = rnorm(50, 3.1, 0.12), CD4 = rnorm(50, 3.0, 0.12), CD8 = rnorm(50, 0.2, 0.1),
        CD62L = rnorm(50, 0.3, 0.1), CD44 = rnorm(50, 3.1, 0.12), CD19 = 0.2, CD11B = 0.2, `NK1.1` = 0.2)
)
colnames(mat_s)[colnames(mat_s) == "NK1.1"] <- "NK1.1"
hs <- hierarchical_gate_by_sample(mat_s, rep(c("EV1", "EV2"), each = 50), "P1")
if (mean(hs$subset[1:50] == "CD4_naive") < 0.6) fail("sample EV1 should be gated naive on its own")
if (mean(hs$subset[51:100] == "CD4_TEM") < 0.6) fail("sample EV2 should be gated TEM on its own")

nv <- Filter(function(s) identical(s$lineage, "CD4_naive"), subset_plot_specs("P1"))
if (!length(nv) || !identical(nv[[1]]$gate, "quad")) fail("CD4 naive subset figure should use a quadrant gate")
qr <- quad_gate_rect(c(0, 7), c(0, 6), 2, 1.2, TRUE, FALSE)
if (qr$xmin < 1.9 || qr$ymin > 0.05 || qr$ymax > 1.25) fail("naive quadrant should be CD62L-high CD44-low")
tem <- Filter(function(s) identical(s$lineage, "CD4_TEM"), subset_plot_specs("P1"))
if (!length(tem) || !identical(tem[[1]]$gate, "half_x")) fail("CD4 TEM figure should use the full CD62L- half, including DN")
tr <- complete_gate_rect(c(0, 7), c(0, 6), 2, 1.2, FALSE, NA)
if (tr$xmax > 2.05 || tr$ymin > 0.05 || tr$ymax < 5.5) {
  fail("TEM gate must cover the origin blob and the CD44+ left half, not a 10-90% hit box")
}
set.seed(14)
gx <- c(rnorm(200, 0.45, 0.12), rnorm(40, 3.2, 0.18))
gy <- c(rnorm(200, 0.40, 0.12), rnorm(40, 0.45, 0.12))
g_naive <- complete_gate_for(gx, gy, nv[[1]], c(0, 7), c(0, 6))
if ((g_naive$xmax - g_naive$xmin) < 3) fail("naive gate must span to the axis max, not a quantile box of hit cells")

cat("OK: P1 T/NK subset labels\n")
