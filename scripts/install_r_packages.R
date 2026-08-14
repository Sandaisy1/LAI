#!/usr/bin/env Rscript
# Only CRAN packages (no GSVA, AnnotationDbi, org.Hs.eg.db, BiocManager).

cran <- c("data.table", "survival")
for (p in cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
message("R packages ready: ", paste(cran, collapse = ", "))
