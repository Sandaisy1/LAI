# Internation cell immune：His+ 靶细胞 + P1/P3（不改原 Flow_* 流程）
# run: Rscript tests/test_ici_flow.R
Sys.setenv(FLOW_FUNCTIONS_ONLY = "1", ICI_FUNCTIONS_ONLY = "1")
Sys.setenv(ICI_FLOW_DIR = tempfile("ici_flow_"))
dir.create(Sys.getenv("ICI_FLOW_DIR"), recursive = TRUE, showWarnings = FALSE)
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/test_ici_flow.R"
root <- dirname(dirname(normalizePath(this_file)))
source(file.path(root, "ICI_Flow_dimred_pipeline.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}
expect <- function(got, want, tag) {
  if (!identical(got, want)) {
    fail(sprintf("%s: got %s, want %s", tag, paste(got, collapse = ","), paste(want, collapse = ",")))
  }
}

# 原 flow_panel_map.json 不得被 ICI 染色表改掉（P1 仍是 25 色，FITC=Perforin）
orig_json <- jsonlite::fromJSON(file.path(root, "flow_panel_map.json"), simplifyVector = FALSE)
orig_p1 <- vapply(orig_json$panels$P1$markers, function(x) x$marker, character(1))
if (length(orig_p1) != 25L) {
  fail(sprintf("original P1 must stay 25 markers, got %s", length(orig_p1)))
}
orig_fitc <- orig_json$panels$P1$markers[[match("Perforin", orig_p1)]]$fluorochrome
expect(orig_fitc, "FITC", "original P1 Perforin remains FITC")
if ("His" %in% orig_p1) fail("original P1 must not gain His")

ici_json <- jsonlite::fromJSON(file.path(root, "ICI_flow_panel_map.json"), simplifyVector = FALSE)
ici_p1 <- vapply(ici_json$panels$P1$markers, function(x) x$marker, character(1))
ici_p3 <- vapply(ici_json$panels$P3$markers, function(x) x$marker, character(1))
if (length(ici_p1) != 14L) fail(sprintf("ICI P1 should have 14 markers, got %s", length(ici_p1)))
if (length(ici_p3) != 20L) fail(sprintf("ICI P3 should have 20 markers, got %s", length(ici_p3)))
if (!is.null(ici_json$panels$P2)) fail("ICI map must not include Panel 2")
expect(ici_json$panels$P1$markers[[match("His", ici_p1)]]$fluorochrome, "FITC", "ICI P1 His-FITC")
expect(ici_json$panels$P3$markers[[match("His", ici_p3)]]$fluorochrome, "FITC", "ICI P3 His-FITC")
expect(ici_json$panels$P3$markers[[match("CD80", ici_p3)]]$fluorochrome, "BUV496", "ICI P3 CD80-BUV496")
expect(ici_json$panels$P3$markers[[match("NK1.1", ici_p3)]]$fluorochrome, "R718", "ICI P3 NK1.1-R718")
expect(ici_json$qc$require_cd45, FALSE, "ICI QC must not require CD45")

# 当前 session 的 panel_map 应是 ICI 表
expect(panel_map$panels$P1$markers[[match("His", vapply(panel_map$panels$P1$markers, function(x) x$marker, character(1)))]]$fluorochrome,
       "FITC", "loaded ICI His-FITC")
p3_mk <- vapply(panel_map$panels$P3$markers, function(x) x$marker, character(1))
expect(panel_map$panels$P3$markers[[match("CD80", p3_mk)]]$fluorochrome, "BUV496", "loaded ICI CD80")
expect(panel_map$panels$P3$markers[[match("NK1.1", p3_mk)]]$fluorochrome, "R718", "loaded ICI NK1.1")

# 通道：His←FITC；P3 CD80←BUV496；P3 NK1.1←R718（Cytek 常把 AF700 写成 R718）
cytek_p1 <- c(
  "FSC-A", "SSC-A", "FVS450-A", "V500-A", "BUV496-A", "BUV805-A", "RB744-A",
  "RB670-A", "RY586-A", "PE-A", "AF700-A", "APC-Cy7-A", "RY703-A", "RB613-A",
  "PE-EF610-A", "FITC-A"
)
map1 <- match_channels(cytek_p1, rep("", length(cytek_p1)), "P1")
expect(map1$channel[map1$marker == "His"], "FITC-A", "P1 His<-FITC-A")
expect(map1$channel[map1$marker == "CD4"], "BUV496-A", "P1 CD4<-BUV496-A")
expect(map1$channel[map1$marker == "NK1.1"], "AF700-A", "P1 NK1.1<-AF700-A")

cytek_p3 <- c(
  "FSC-A", "SSC-A", "FVS450-A", "V500-A", "BUV805-A", "RB705-A", "RB670-A",
  "BUV496-A", "APC-Cy7-A", "BV750-A", "R718-A", "BUV737-A", "BV605-A",
  "RY703-A", "RY775-A", "RB744-A", "RB613-A", "AF647-A", "BUV661-A",
  "BV786-A", "FITC-A", "PE-EF610-A"
)
map3 <- match_channels(cytek_p3, rep("", length(cytek_p3)), "P3")
expect(map3$channel[map3$marker == "His"], "FITC-A", "P3 His<-FITC-A")
expect(map3$channel[map3$marker == "CD80"], "BUV496-A", "P3 CD80<-BUV496-A")
expect(map3$channel[map3$marker == "NK1.1"], "R718-A", "P3 NK1.1<-R718-A")
expect(map3$channel[map3$marker == "FceRI"], "PE-EF610-A", "P3 FceRI<-PE-EF610-A")

# 文件名：EV-1 / H-3 / ZZX-EV-1；不要把 P3 的 3 当成技术重复；跳过 P2
ev <- parse_fcs_filename("EV-1_P1_unmixed.fcs")
if (is.null(ev) || !identical(ev$group, "EV") || !identical(ev$sample, "EV-1") ||
    !identical(ev$bio_sample, "EV-1") || !identical(ev$panel, "P1") ||
    !identical(ev$replicate, "1")) {
  fail(sprintf("EV-1_P1_unmixed.fcs -> EV-1 P1, got group=%s sample=%s panel=%s",
               if (is.null(ev)) "NULL" else ev$group,
               if (is.null(ev)) "NULL" else ev$sample,
               if (is.null(ev)) "NULL" else ev$panel))
}
zzx <- parse_fcs_filename("ZZX-EV-1_P1_unmixed.fcs")
if (is.null(zzx) || !identical(zzx$group, "EV") || !identical(zzx$sample, "EV-1") ||
    !identical(zzx$panel, "P1")) {
  fail("ZZX-EV-1_P1_unmixed.fcs should parse as EV sample EV-1")
}
hh <- parse_fcs_filename("H-3_P3_unmixed.fcs")
if (is.null(hh) || !identical(hh$group, "H") || !identical(hh$sample, "H-3") ||
    !identical(hh$panel, "P3") || !is.na(hh$tech_rep)) {
  fail(sprintf("H-3_P3 must be bio H-3 with no tech rep, got sample=%s tech=%s panel=%s",
               if (is.null(hh)) "NULL" else hh$sample,
               if (is.null(hh)) "NULL" else hh$tech_rep,
               if (is.null(hh)) "NULL" else hh$panel))
}
h1 <- parse_fcs_filename("ZZX-H-2_P3_unmixed.fcs")
if (is.null(h1) || !identical(h1$group, "H") || !identical(h1$sample, "H-2")) {
  fail("ZZX-H-2_P3_unmixed.fcs should parse as H sample H-2")
}
if (!is.null(parse_fcs_filename("EV-1_P2_unmixed.fcs"))) {
  fail("ICI parser must skip Panel 2")
}
if (!is.null(parse_fcs_filename("ZZX_EV1-1_P1_unmixed.fcs")) &&
    identical(parse_fcs_filename("ZZX_EV1-1_P1_unmixed.fcs")$group, "H")) {
  fail("EV must not parse as H")
}

# QC：只去双联体和死细胞，不要因为 CD45- 丢掉 His+ 靶细胞；P1 不要紧淋巴门
set.seed(7)
n_leu <- 80L
n_tgt <- 80L
n_dead <- 8L
exprs_qc <- cbind(
  `FSC-A` = c(rnorm(n_leu, 80000, 4000), rnorm(n_tgt, 110000, 5000), rnorm(n_dead, 80000, 4000)),
  `SSC-A` = c(rnorm(n_leu, 20000, 2000), rnorm(n_tgt, 45000, 3000), rnorm(n_dead, 20000, 2000)),
  `FSC-H` = c(rnorm(n_leu, 80000, 4000), rnorm(n_tgt, 110000, 5000), rnorm(n_dead, 80000, 4000)),
  `L/D` = c(rnorm(n_leu + n_tgt, 200, 30), rnorm(n_dead, 8000, 200)),
  CD45 = c(rnorm(n_leu, 12000, 800), rnorm(n_tgt, 80, 20), rnorm(n_dead, 12000, 800))
)
map_qc <- data.frame(
  marker = c("L/D", "CD45"),
  channel_index = c(4L, 5L),
  stringsAsFactors = FALSE
)
chn <- colnames(exprs_qc)
keep_ici <- qc_filter_matrix(exprs_qc, chn, map_qc, "P1")
keep_orig <- qc_filter_matrix_flow(exprs_qc, chn, map_qc, "P1")
tgt_idx <- seq.int(n_leu + 1L, n_leu + n_tgt)
if (mean(keep_ici[tgt_idx]) < 0.7) {
  fail(sprintf("ICI QC must keep His+ CD45- / large-FSC targets, kept %s", mean(keep_ici[tgt_idx])))
}
if (mean(keep_orig[tgt_idx]) > 0.4) {
  fail("original QC (CD45+ and P1 lymph gate) should drop most CD45- large cells")
}

# His+ CD45- → Target；His+ CD45+ CD3+ 仍是 T，不要改成 Target
set.seed(11)
n <- 70L
mk_block <- function(his, cd45, cd3, cd4, cd8 = 0.2, nkp = 0.2) {
  cbind(
    His = rnorm(n, his, 0.12),
    CD45 = rnorm(n, cd45, 0.12),
    CD3 = rnorm(n, cd3, 0.12),
    CD4 = rnorm(n, cd4, 0.12),
    CD8 = rnorm(n, cd8, 0.08),
    CD19 = rnorm(n, 0.2, 0.08),
    CD11B = rnorm(n, 0.2, 0.08),
    NKp46 = rnorm(n, nkp, 0.08),
    `NK1.1` = rnorm(n, nkp, 0.08),
    CD62L = rnorm(n, 3.1, 0.12),
    CD44 = rnorm(n, 0.3, 0.1),
    CD25 = rnorm(n, 0.3, 0.1),
    CD69 = rnorm(n, 0.3, 0.1)
  )
}
mat <- rbind(
  mk_block(3.2, 0.3, 0.2, 0.2),
  mk_block(3.2, 3.0, 3.2, 3.0),
  mk_block(0.3, 3.0, 3.2, 3.0)
)
h <- hierarchical_gate(mat, "P1")
tgt <- seq_len(n)
his_t <- seq.int(n + 1L, 2L * n)
his_neg_t <- seq.int(2L * n + 1L, 3L * n)
if (mean(h$major[tgt] == "Target") < 0.85 || mean(h$subset[tgt] == "Target") < 0.85) {
  fail(sprintf("His+ CD45- must be Target, got major=%s subset=%s",
               paste(unique(h$major[tgt]), collapse = ","),
               paste(unique(h$subset[tgt]), collapse = ",")))
}
if (any(h$major[his_t] == "Target") || any(h$subset[his_t] == "Target")) {
  fail("His+ CD45+ T cells must keep immune labels, not Target")
}
if (mean(h$major[his_t] == "CD4") < 0.8) {
  fail(sprintf("His+ CD45+ CD3+ CD4+ should stay CD4, got %s",
               paste(unique(h$major[his_t]), collapse = ",")))
}
if (any(h$major[his_neg_t] == "Target")) {
  fail("His- CD45+ T cells must not be labeled Target")
}

# 只有 His+ CD45+ T 时也不要改成 Target
h_t_only <- hierarchical_gate(mk_block(3.2, 3.0, 3.2, 3.0), "P1")
if (any(h_t_only$major == "Target") || mean(h_t_only$major == "CD4") < 0.85) {
  fail(sprintf("His+ CD45+ only tube must stay CD4, got %s",
               paste(unique(h_t_only$major), collapse = ",")))
}
expect(celltype_label("Target", "P1"), "His+ target", "label-his-target")

cat("PASS test_ici_flow.R\n")
