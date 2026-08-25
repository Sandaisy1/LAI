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
  "BV650-A", "BUV661-A", "BV421-A", "BV786-A", "R718-A", "APC-Cy7-A",
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

cat("OK: Cytek channel matching\n")
