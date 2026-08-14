#!/usr/bin/env Rscript
# One-time install for scripts/go_pathway_per_term_analysis.R

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

cran <- c("data.table", "survival")
bioc <- c("GSVA", "org.Hs.eg.db", "AnnotationDbi")

for (p in cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
for (p in bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}
message("R packages ready.")
