# P3 cluster labels + shared P1 palette (run: Rscript tests/test_flow_p3_annotate.R)
Sys.setenv(FLOW_FUNCTIONS_ONLY = "1")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/test_flow_p3_annotate.R"
source(file.path(dirname(normalizePath(this_file)), "..", "Flow_dimred_pipeline.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}

mk <- function(...) {
  pops <- demo_means_p3()
  extra <- list(...)
  row <- pops[[1]]
  row[] <- 0.2
  for (nm in names(extra)) row[[nm]] <- extra[[nm]]
  as.data.frame(t(row), stringsAsFactors = FALSE)
}

lab_of <- function(row) {
  rownames(row) <- "C1"
  annotate_clusters(row, "P3")$lineage
}

expect <- function(got, want, tag) {
  if (!identical(got, want)) fail(sprintf("%s: got %s, want %s", tag, got, want))
}

expect(lab_of(mk(CD11B = 3.2, LY6G = 3.3, LY6C = 1.5)), "Neutrophil", "neu")
expect(lab_of(mk(CD11B = 3.1, LY6C = 3.2, LY6G = 0.2, `F4/80` = 0.4)), "Mono_Ly6Chi", "ly6c-hi")
expect(lab_of(mk(CD11B = 3.0, LY6C = 0.3, LY6G = 0.2, `F4/80` = 0.4)), "Mono_Ly6Clo", "ly6c-lo")
expect(lab_of(mk(CD11B = 3.0, `F4/80` = 3.2, `I-A/I-E` = 1.8)), "Macrophage", "mac")
expect(lab_of(mk(CD11B = 3.0, `F4/80` = 3.0, iNOS = 3.2, CD86 = 2.6)), "M1_like_Mac", "m1")
expect(lab_of(mk(CD11B = 3.0, `F4/80` = 3.0, CD206 = 3.1, `ARG-1` = 2.8)), "M2_like_Mac", "m2")
expect(lab_of(mk(CD11C = 3.2, `I-A/I-E` = 3.3, CD80 = 2.2, CD11B = 1.4)), "DC", "dc")
expect(lab_of(mk(CD11C = 3.1, CD103 = 3.0, `I-A/I-E` = 3.0, CD11B = 0.6)), "cDC1_CD103", "cdc1")
expect(lab_of(mk(`Siglec-F` = 3.2, CCR3 = 2.8, CD11B = 2.4, LY6G = 0.2)), "Eosinophil", "eos")
expect(lab_of(mk(FceRI = 3.1, CD200R3 = 2.8, CD11B = 1.6)), "Basophil_mast", "baso")
# CCR3 高但 Siglec-F 不高 → 不是嗜酸
ccr3_only <- lab_of(mk(CCR3 = 3.2, `Siglec-F` = 0.2, CD11B = 2.4, LY6C = 0.3, `F4/80` = 0.4))
if (identical(ccr3_only, "Eosinophil")) fail("CCR3-only must not be Eosinophil")
expect(lab_of(mk(CD19 = 3.3, CD3 = 0.2, CD11B = 0.3)), "B", "dump-B")
expect(lab_of(mk(CD3 = 3.2, CD19 = 0.2, CD11B = 0.3)), "T", "dump-T")
expect(lab_of(mk(`NK1.1` = 3.2, CD3 = 0.2, CD11B = 0.3)), "NK", "dump-NK")
# IL-10 / TGF-b 背景不得把巨噬打成 M2
expect(lab_of(mk(CD11B = 3.0, `F4/80` = 3.2, `IL-10` = 2.8, `TGF-b` = 2.6, CD206 = 0.2, `ARG-1` = 0.2)), "Macrophage", "il10-not-m2")

p2 <- c("Naive B", "Atypical B", "IgM memory B", "Memory B", "Switched B", "Activated B", "Plasma")
p2_cols <- unname(pal_celltype[p2])
if (length(unique(p2_cols)) < 6) fail("P2 colors collapsed to a purple ramp")
if (!identical(unname(pal_celltype[["Naive B"]]), unname(pal_celltype[["CD4 naive"]]))) {
  fail("P2 Naive B should reuse P1 CD4 naive hue")
}
if (!identical(unname(pal_celltype[["Eosinophil"]]), unname(pal_celltype[["CD4 activated"]]))) {
  fail("P3 Eosinophil should reuse P1 red, not pumpkin orange")
}
if (identical(unname(pal_celltype[["Eosinophil"]]), "#D35400")) fail("old P3 orange still in palette")

set.seed(1)
p2p <- demo_means_p2()
f2 <- names(p2p[[1]])
bp <- function(name, n) {
  m <- matrix(unlist(p2p[[name]]), n, length(f2), byrow = TRUE)
  colnames(m) <- f2
  m + matrix(rnorm(n * length(f2), 0, 0.1), n)
}
m2 <- rbind(bp("Naive_B", 100), bp("Plasma", 80), bp("Switched_B", 80), bp("Activated_B", 80), bp("IgM_memory", 80))
h2 <- hierarchical_gate(m2, "P2")
if (!all(h2$major %in% c("Naive_B", "Memory_B", "Plasma"))) {
  fail(sprintf("P2 layer1 should be Naive/Memory/Plasma, got %s", paste(unique(h2$major), collapse = ",")))
}
if (length(unique(h2$subset)) < 4) fail(sprintf("P2 subsets too few: %s", paste(unique(h2$subset), collapse = ",")))

p3p <- demo_means_p3()
f3 <- names(p3p[[1]])
bp3 <- function(name, n) {
  m <- matrix(unlist(p3p[[name]]), n, length(f3), byrow = TRUE)
  colnames(m) <- f3
  m + matrix(rnorm(n * length(f3), 0, 0.1), n)
}
m3 <- rbind(
  bp3("Neutrophil", 80), bp3("Eosinophil", 80), bp3("Macrophage", 80),
  bp3("M1_like", 70), bp3("M2_like", 70), bp3("DC", 80), bp3("Mono_Ly6Chi", 80)
)
h3 <- hierarchical_gate(m3, "P3")
if (sum(h3$major == "Myeloid") < 0.8 * nrow(m3)) fail("P3 layer1 should be mostly Myeloid")
if (length(unique(h3$subset)) < 5) fail(sprintf("P3 subsets too few: %s", paste(unique(h3$subset), collapse = ",")))

cat("OK: P3 labels and P1-shared palette\n")
