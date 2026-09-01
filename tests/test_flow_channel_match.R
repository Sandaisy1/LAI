# Cytek-style FCS channel matching (run: Rscript tests/test_flow_channel_match.R)
Sys.setenv(FLOW_FUNCTIONS_ONLY = "1")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/test_flow_channel_match.R"
source(file.path(dirname(normalizePath(this_file)), "..", "Flow_dimred_pipeline.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}

expect <- function(got, want, tag) {
  if (!identical(got, want)) fail(sprintf("%s: got %s, want %s", tag, paste(got, collapse = ","), paste(want, collapse = ",")))
}

# 旧逻辑会把 BUV496-A 收成 BUV496A，从而一个都配不上
cytek_p1 <- c(
  "FSC-A", "SSC-A", "FVS450-A", "V500-A", "BUV496-A", "BUV805-A", "RB744-A",
  "BV750-A", "RB670-A", "RY586-A", "RB705-A", "BV711-A", "PE-A", "BV605-A",
  "BV650-A", "BUV661-A", "BV421-A", "BV786-A", "AF700-A", "APC-Cy7-A",
  "RY703-A", "RB613-A", "PE-EF610-A", "APC-A", "BUV737-A", "PE-Cy7-A", "FITC-A"
)
desc_blank <- rep("", length(cytek_p1))
map <- match_channels(cytek_p1, desc_blank, "P1")
if (sum(!is.na(map$channel_index)) < 20) {
  fail(sprintf("Cytek -A names matched %d/25, expected >=20", sum(!is.na(map$channel_index))))
}
expect(map$channel[map$marker == "CD4"], "BUV496-A", "CD4<-BUV496-A")
expect(map$channel[map$marker == "CD45"], "V500-A", "CD45<-V500-A")
expect(map$channel[map$marker == "CD8"], "RB744-A", "CD8<-RB744-A")
expect(map$channel[map$marker == "CD8b"], "BV711-A", "CD8b<-BV711-A")
expect(map$channel[map$marker == "NKp46"], "PE-A", "NKp46<-PE-A")
expect(map$channel[map$marker == "CD11B"], "APC-Cy7-A", "CD11B<-APC-Cy7-A")
expect(map$channel[map$marker == "L/D"], "FVS450-A", "LD<-FVS450-A")
expect(map$channel[map$marker == "NK1.1"], "AF700-A", "P1 NK1.1<-AF700-A")
expect(map$channel[map$marker == "NKG2D"], "PE-EF610-A", "P1 NKG2D<-PE-EF610-A")
expect(map$channel[map$marker == "TNF-a"], "APC-A", "P1 TNF-a<-APC-A")
expect(axis_fl_label("P1", "NKG2D"), "NKG2D-PE-eFluor 610", "P1 NKG2D axis")
expect(axis_fl_label("P1", "CD11B"), "CD11B-APC-Cy7", "P1 CD11B axis")
expect(axis_fl_label("P1", "GZMB"), "GZMB-PE-Cy7", "P1 GZMB axis")
expect(axis_fl_label("P3", "FceRI"), "FceRI-PE-eFluor 610", "P3 FceRI axis")

# 抗体写 AF700，Cytek unmixed 常把同一检测器导出成 R718-A
cytek_p1_r718 <- cytek_p1
cytek_p1_r718[cytek_p1_r718 == "AF700-A"] <- "R718-A"
map_r718 <- match_channels(cytek_p1_r718, desc_blank, "P1")
expect(map_r718$channel[map_r718$marker == "NK1.1"], "R718-A", "P1 NK1.1<-R718-A")
cytek_p1_y710 <- cytek_p1
cytek_p1_y710[cytek_p1_y710 == "AF700-A"] <- "Y710-A"
map_y710 <- match_channels(cytek_p1_y710, desc_blank, "P1")
expect(map_y710$channel[map_y710$marker == "NK1.1"], "Y710-A", "P1 NK1.1<-Y710-A")

# 表上写 PE-eFluor™ 610，Cytek 也可能导出这个全名
cytek_p1_ef <- cytek_p1
cytek_p1_ef[cytek_p1_ef == "PE-EF610-A"] <- "PE-eFluor 610-A"
map_ef <- match_channels(cytek_p1_ef, desc_blank, "P1")
expect(map_ef$channel[map_ef$marker == "NKG2D"], "PE-eFluor 610-A", "P1 NKG2D<-PE-eFluor 610-A")

# 标志物写在 name 里、荧光素写在 desc
map2 <- match_channels(
  c("FSC-A", "CD45", "CD4", "CD8", "CD8b"),
  c("", "V500", "BUV496", "RB744", "BV711"),
  "P1"
)
expect(map2$channel[map2$marker == "CD8"], "CD8", "CD8 name not stolen by CD8b")
expect(map2$channel[map2$marker == "CD8b"], "CD8b", "CD8b own channel")

# Comp 前缀 + 合在一起的 "CD45 V500"
map3 <- match_channels(
  c("FJComp-BUV496-A", "CD45 V500"),
  c("CD4", ""),
  "P1"
)
expect(map3$channel[map3$marker == "CD4"], "FJComp-BUV496-A", "FJComp CD4")
expect(map3$channel[map3$marker == "CD45"], "CD45 V500", "CD45 V500 combo")

# 检测器名对不上时必须全空，不能乱配
map4 <- match_channels(c("FSC-A", "V1-A", "B1-A", "YG1-A", "R1-A"), rep("", 5), "P1")
if (any(!is.na(map4$channel_index))) fail("detector-only names should not match P1 markers")

# 旧 bug：norm_id('BUV496-A')=='BUV496A' 不等于 'BUV496'
if (!("BUV496" %in% fcs_tokens("BUV496-A"))) fail("fcs_tokens must strip -A")
if ("CD8" %in% fcs_tokens("CD8b-A")) fail("CD8 must not be a token of CD8b-A")

# P2 新表：CD80 从 BUV496 改到 APC（与 P1 的 TNF-a 同色但不同管子）
cytek_p2 <- c(
  "FSC-A", "SSC-A", "FVS450-A", "V500-A", "BV750-A", "BV605-A", "APC-A",
  "RB670-A", "RB613-A", "RB780-A", "RB705-A", "FITC-A", "PE-A"
)
map_p2 <- match_channels(cytek_p2, rep("", length(cytek_p2)), "P2")
expect(map_p2$channel[map_p2$marker == "CD80"], "APC-A", "P2 CD80<-APC-A")
expect(map_p2$channel[map_p2$marker == "CD86"], "RB670-A", "P2 CD86<-RB670-A")
expect(map_p2$channel[map_p2$marker == "IgD"], "FITC-A", "P2 IgD<-FITC-A")
if (!is.na(map_p2$channel[map_p2$marker == "CD80"]) &&
    identical(map_p2$channel[map_p2$marker == "CD80"], "BUV496-A")) {
  fail("P2 CD80 must not still match BUV496")
}

# P3：CD80=BUV496（TNF-a 才是 APC）；NK1.1 只认 AF700
cytek_p3 <- c(
  "FSC-A", "SSC-A", "FVS450-A", "V500-A", "BUV805-A", "RB705-A", "RB670-A",
  "APC-A", "APC-Cy7-A", "BV750-A", "AF700-A", "BUV737-A", "BV605-A", "BUV496-A"
)
map_p3 <- match_channels(cytek_p3, rep("", length(cytek_p3)), "P3")
expect(map_p3$channel[map_p3$marker == "CD80"], "BUV496-A", "P3 CD80<-BUV496-A")
expect(map_p3$channel[map_p3$marker == "TNF-a"], "APC-A", "P3 TNF-a<-APC-A")
expect(map_p3$channel[map_p3$marker == "NK1.1"], "AF700-A", "P3 NK1.1<-AF700-A")
if (!is.na(map_p3$channel[map_p3$marker == "CD80"]) &&
    identical(map_p3$channel[map_p3$marker == "CD80"], "APC-A")) {
  fail("P3 CD80 must not still match APC (TNF-a is APC)")
}
cytek_p3_r718 <- cytek_p3
cytek_p3_r718[cytek_p3_r718 == "AF700-A"] <- "R718-A"
map_p3_r718 <- match_channels(cytek_p3_r718, rep("", length(cytek_p3_r718)), "P3")
expect(map_p3_r718$channel[map_p3_r718$marker == "NK1.1"], "R718-A", "P3 NK1.1<-R718-A")

# 文件名 EV1 / H-1 / ZZX_EV1-1
ev <- parse_fcs_filename("EV1_P1_unmixed.fcs")
if (is.null(ev) || !identical(ev$group, "EV") || !identical(ev$sample, "EV1")) {
  fail("EV1_P1_unmixed.fcs should parse as group EV sample EV1")
}
hh <- parse_fcs_filename("H-3_P2_unmixed.fcs")
if (is.null(hh) || !identical(hh$group, "H") || !identical(hh$sample, "H3")) {
  fail("H-3_P2_unmixed.fcs should parse as group H sample H3")
}
p3 <- parse_fcs_filename("EV1-P3_unmixed.fcs")
if (is.null(p3) || !identical(p3$panel, "P3") || !identical(p3$sample, "EV1")) {
  fail("EV1-P3_unmixed.fcs should parse as panel P3 sample EV1")
}
zzx <- parse_fcs_filename("ZZX_EV1-1_P1_unmixed.fcs")
if (is.null(zzx) || !identical(zzx$group, "EV") || !identical(zzx$sample, "EV1-1") ||
    !identical(zzx$bio_sample, "EV1") || !identical(zzx$tech_rep, "1")) {
  fail("ZZX_EV1-1_P1_unmixed.fcs should be EV tech EV1-1 / bio EV1")
}
zzxh <- parse_fcs_filename("ZZX_H2-2_P3_unmixed.fcs")
if (is.null(zzxh) || !identical(zzxh$group, "H") || !identical(zzxh$sample, "H2-2") ||
    !identical(zzxh$bio_sample, "H2")) {
  fail("ZZX_H2-2_P3_unmixed.fcs should be H tech H2-2 / bio H2")
}
if (!is.null(parse_fcs_filename("ZZX_EV1-1_P1_unmixed.fcs")) &&
    identical(parse_fcs_filename("ZZX_EV1-1_P1_unmixed.fcs")$group, "H")) {
  fail("ZZX_EV must not parse as H")
}

freq_tech <- data.frame(
  sample = c("EV1-1", "EV1-2", "EV2-1", "EV2-2", "EV3-1", "EV3-2",
             "H1-1", "H1-2", "H2-1", "H2-2", "H3-1", "H3-2"),
  bio_sample = c("EV1", "EV1", "EV2", "EV2", "EV3", "EV3",
                 "H1", "H1", "H2", "H2", "H3", "H3"),
  group = c(rep("EV", 6), rep("H", 6)),
  lineage = "NK",
  percent = c(10, 12, 11, 9, 10, 10, 4, 6, 5, 5, 5, 7),
  stringsAsFactors = FALSE
)
if (!identical(which_extreme_bio(c(1, 5, 6)), 1L)) {
  fail("1,5,6: drop the min (farther from median)")
}
if (!identical(which_extreme_bio(c(11, 10, 10)), 1L)) {
  fail("11,10,10: drop the max outlier")
}
if (!identical(which_extreme_bio(c(1, 2, 10)), 3L)) {
  fail("1,2,10: drop the max outlier")
}
if (!identical(which_extreme_bio(c(1, 5, 9)), 3L)) {
  fail("equal distance from median: drop the max")
}
if (length(which_extreme_bio(c(1, 2))) != 0L) {
  fail("n<3 must keep all bio-reps")
}

st_full <- compare_group_freq(freq_tech, "lineage", trim_bio = FALSE)
if (!identical(as.integer(st_full$n_EV[1]), 3L) || !identical(as.integer(st_full$n_H[1]), 3L)) {
  fail(sprintf("tech replicates must average to n=3 bio before trim, got n_EV=%s n_H=%s",
               st_full$n_EV[1], st_full$n_H[1]))
}

st_bio <- compare_group_freq(freq_tech, "lineage")
if (!identical(as.integer(st_bio$n_EV[1]), 2L) || !identical(as.integer(st_bio$n_H[1]), 2L)) {
  fail(sprintf("after dropping 1 extreme bio-rep, n must be 2 vs 2, got n_EV=%s n_H=%s",
               st_bio$n_EV[1], st_bio$n_H[1]))
}
if (abs(st_bio$mean_EV[1] - 10) > 1e-8) {
  fail(sprintf("EV 11/10/10 should drop 11, mean=10, got %s", st_bio$mean_EV[1]))
}
if (abs(st_bio$mean_H[1] - 5) > 1e-8) {
  fail(sprintf("H 5/5/6 should drop 6, mean=5, got %s", st_bio$mean_H[1]))
}
if (!identical(as.character(st_bio$dropped_EV[1]), "EV1")) {
  fail(sprintf("should drop EV1, got %s", st_bio$dropped_EV[1]))
}
if (!identical(as.character(st_bio$dropped_H[1]), "H3")) {
  fail(sprintf("should drop H3, got %s", st_bio$dropped_H[1]))
}

if (!grepl("fuction of cell", flow_primary_data_dir, fixed = TRUE)) {
  fail("flow data/results dir should be E:/R/fuction of cell")
}

cat("OK: Cytek channel matching\n")
