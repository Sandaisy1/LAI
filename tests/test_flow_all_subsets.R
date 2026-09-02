# Combined P1+P2+P3 subset frequencies (run: Rscript tests/test_flow_all_subsets.R)
Sys.setenv(FLOW_FUNCTIONS_ONLY = "1")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/test_flow_all_subsets.R"
root <- dirname(normalizePath(this_file))
source(file.path(root, "..", "Flow_dimred_pipeline.R"))
source(file.path(root, "..", "Flow_dimred_all_subsets.R"))

fail <- function(msg) {
  cat("FAIL:", msg, "\n")
  quit(status = 1)
}

if (!identical(subset_role("P1", "CD4_naive"), "focus")) fail("P1 T subset should be focus")
if (!identical(subset_role("P1", "NK_immature"), "focus")) fail("P1 NK subset should be focus")
if (!identical(subset_role("P1", "NKT_CD4"), "focus")) fail("P1 NKT subset should be focus")
if (!identical(subset_role("P1", "B"), "dump")) fail("P1 B should be dump")
if (!identical(subset_role("P2", "Naive_B"), "focus")) fail("P2 B subset should be focus")
if (!identical(subset_role("P3", "Neutrophil"), "focus")) fail("P3 neutrophil should be focus")
if (!identical(subset_role("P3", "T"), "dump")) fail("P3 T should be dump")

mk_cells <- function(sample, group, lineage, n) {
  data.frame(
    sample = rep(sample, n),
    group = rep(group, n),
    lineage = rep(lineage, n),
    stringsAsFactors = FALSE
  )
}

p1 <- rbind(
  mk_cells("EV1", "EV", "CD4_naive", 80),
  mk_cells("EV1", "EV", "B", 20),
  mk_cells("EV2", "EV", "CD4_naive", 70),
  mk_cells("EV2", "EV", "B", 30),
  mk_cells("EV3", "EV", "CD4_naive", 75),
  mk_cells("EV3", "EV", "B", 25),
  mk_cells("H1", "H", "CD4_naive", 40),
  mk_cells("H1", "H", "B", 60),
  mk_cells("H2", "H", "CD4_naive", 45),
  mk_cells("H2", "H", "B", 55),
  mk_cells("H3", "H", "CD4_naive", 50),
  mk_cells("H3", "H", "B", 50)
)
p2 <- rbind(
  mk_cells("EV1", "EV", "Naive_B", 90),
  mk_cells("EV1", "EV", "Plasma", 10),
  mk_cells("EV2", "EV", "Naive_B", 85),
  mk_cells("EV2", "EV", "Plasma", 15),
  mk_cells("EV3", "EV", "Naive_B", 88),
  mk_cells("EV3", "EV", "Plasma", 12),
  mk_cells("H1", "H", "Naive_B", 50),
  mk_cells("H1", "H", "Plasma", 50),
  mk_cells("H2", "H", "Naive_B", 55),
  mk_cells("H2", "H", "Plasma", 45),
  mk_cells("H3", "H", "Naive_B", 60),
  mk_cells("H3", "H", "Plasma", 40)
)
p3 <- rbind(
  mk_cells("EV1", "EV", "Neutrophil", 70),
  mk_cells("EV1", "EV", "M1_like_Mac", 30),
  mk_cells("EV2", "EV", "Neutrophil", 65),
  mk_cells("EV2", "EV", "M1_like_Mac", 35),
  mk_cells("EV3", "EV", "Neutrophil", 68),
  mk_cells("EV3", "EV", "M1_like_Mac", 32),
  mk_cells("H1", "H", "Neutrophil", 40),
  mk_cells("H1", "H", "M1_like_Mac", 60),
  mk_cells("H2", "H", "Neutrophil", 42),
  mk_cells("H2", "H", "M1_like_Mac", 58),
  mk_cells("H3", "H", "Neutrophil", 38),
  mk_cells("H3", "H", "M1_like_Mac", 62)
)

freq <- collect_all_subset_frequencies(list(P1 = p1, P2 = p2, P3 = p3))
if (is.null(freq) || nrow(freq) < 6) fail("combined frequency table too small")
if (!all(c("panel", "subset_label", "role", "percent") %in% names(freq))) {
  fail("combined frequency missing columns")
}
if (any(freq$subset_label == "CD4 naive")) fail("subset labels must keep panel prefix")
if (!any(grepl("^P1 · ", freq$subset_label))) fail("P1 labels missing")
if (!any(grepl("^P2 · ", freq$subset_label))) fail("P2 labels missing")
if (!any(grepl("^P3 · ", freq$subset_label))) fail("P3 labels missing")

sums <- aggregate(percent ~ sample + panel, data = freq, FUN = sum)
if (any(abs(sums$percent - 100) > 1e-6)) fail("within-panel percents must sum to 100")

p1b <- freq[freq$panel == "P1" & freq$lineage == "B", ]
if (!nrow(p1b) || !all(p1b$role == "dump")) fail("P1 B should be dump in combined table")
p3n <- freq[freq$panel == "P3" & freq$lineage == "Neutrophil", ]
if (!nrow(p3n) || !all(p3n$role == "focus")) fail("P3 neutrophil should be focus")

st <- all_subset_stats(freq)
if (!"delta_H_minus_EV" %in% names(st)) fail("stats missing delta_H_minus_EV")
if (length(unique(st$panel)) < 3) fail("stats should cover three panels")
if (!all(as.integer(st$n_EV) == 2L) || !all(as.integer(st$n_H) == 2L)) {
  fail("all-subset stats should drop 1 extreme bio-rep (n=2)")
}
freq_plot <- all_subset_freq_for_plots(freq)
if (is.null(freq_plot) || nrow(freq_plot) >= nrow(freq)) {
  fail("plot table should have fewer rows after dropping 1 bio-rep per group")
}

td <- tempfile("flow_all_subsets")
dir.create(file.path(td, "P1"), recursive = TRUE)
dir.create(file.path(td, "P2"), recursive = TRUE)
dir.create(file.path(td, "P3"), recursive = TRUE)
utils::write.csv(p1, file.path(td, "P1", "P1_cell_embeddings.csv"), row.names = FALSE)
utils::write.csv(p2, file.path(td, "P2", "P2_cell_embeddings.csv"), row.names = FALSE)
utils::write.csv(p3, file.path(td, "P3", "P3_cell_embeddings.csv"), row.names = FALSE)
res <- export_all_subsets_analysis(td)
if (is.null(res) || !dir.exists(file.path(td, "all_subsets"))) fail("export did not create all_subsets/")
need_files <- c(
  "ALL_SUBSETS_NOTE.txt",
  "all_subsets_frequency_by_sample.csv",
  "all_subsets_H_vs_EV_stats.csv",
  "all_subsets_frequency_H_vs_EV.pdf",
  "all_subsets_focus_delta_lollipop.pdf",
  "all_subsets_frequency_by_bio_trimmed.csv"
)
miss <- need_files[!file.exists(file.path(td, "all_subsets", need_files))]
if (length(miss)) fail(sprintf("missing outputs: %s", paste(miss, collapse = ", ")))

unlink(td, recursive = TRUE)
cat("OK: all-subset combined frequency analysis\n")
