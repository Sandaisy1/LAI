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
  mk_cells("T-1", "T", "CD4_naive", 80),
  mk_cells("T-1", "T", "B", 20),
  mk_cells("T-2", "T", "CD4_naive", 70),
  mk_cells("T-2", "T", "B", 30),
  mk_cells("T-3", "T", "CD4_naive", 75),
  mk_cells("T-3", "T", "B", 25),
  mk_cells("T6-1", "T6", "CD4_naive", 40),
  mk_cells("T6-1", "T6", "B", 60),
  mk_cells("T6-2", "T6", "CD4_naive", 45),
  mk_cells("T6-2", "T6", "B", 55),
  mk_cells("T6-3", "T6", "CD4_naive", 50),
  mk_cells("T6-3", "T6", "B", 50)
)
p2 <- rbind(
  mk_cells("T-1", "T", "Naive_B", 90),
  mk_cells("T-1", "T", "Plasma", 10),
  mk_cells("T-2", "T", "Naive_B", 85),
  mk_cells("T-2", "T", "Plasma", 15),
  mk_cells("T-3", "T", "Naive_B", 88),
  mk_cells("T-3", "T", "Plasma", 12),
  mk_cells("T6-1", "T6", "Naive_B", 50),
  mk_cells("T6-1", "T6", "Plasma", 50),
  mk_cells("T6-2", "T6", "Naive_B", 55),
  mk_cells("T6-2", "T6", "Plasma", 45),
  mk_cells("T6-3", "T6", "Naive_B", 60),
  mk_cells("T6-3", "T6", "Plasma", 40)
)
p3 <- rbind(
  mk_cells("T-1", "T", "Neutrophil", 70),
  mk_cells("T-1", "T", "M1_like_Mac", 30),
  mk_cells("T-2", "T", "Neutrophil", 65),
  mk_cells("T-2", "T", "M1_like_Mac", 35),
  mk_cells("T-3", "T", "Neutrophil", 68),
  mk_cells("T-3", "T", "M1_like_Mac", 32),
  mk_cells("T6-1", "T6", "Neutrophil", 40),
  mk_cells("T6-1", "T6", "M1_like_Mac", 60),
  mk_cells("T6-2", "T6", "Neutrophil", 42),
  mk_cells("T6-2", "T6", "M1_like_Mac", 58),
  mk_cells("T6-3", "T6", "Neutrophil", 38),
  mk_cells("T6-3", "T6", "M1_like_Mac", 62)
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
if (!"delta_T6_minus_T" %in% names(st)) fail("stats missing delta")
if (length(unique(st$panel)) < 3) fail("stats should cover three panels")

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
  "all_subsets_T6_vs_T_stats.csv",
  "all_subsets_frequency_T6_vs_T.pdf",
  "all_subsets_focus_delta_lollipop.pdf"
)
miss <- need_files[!file.exists(file.path(td, "all_subsets", need_files))]
if (length(miss)) fail(sprintf("missing outputs: %s", paste(miss, collapse = ", ")))

unlink(td, recursive = TRUE)
cat("OK: all-subset combined frequency analysis\n")
