#!/usr/bin/env Rscript
# Install CRAN + Bioconductor packages used by R/analyze_tg_rnaseq.R

cran_pkgs <- c(
  "ggplot2", "ggrepel", "pheatmap", "readxl", "data.table", "stringr", "scales"
)

install_if_missing <- function(pkgs, installer) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      installer(p)
    }
  }
}

install_if_missing(cran_pkgs, function(p) {
  install.packages(p, repos = "https://cloud.r-project.org")
})

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_pkgs <- c(
  "clusterProfiler", "org.Hs.eg.db", "enrichplot", "ReactomePA", "DOSE"
)
install_if_missing(bioc_pkgs, function(p) {
  BiocManager::install(p, update = FALSE, ask = FALSE)
})

message("Package check complete.")
for (p in c(cran_pkgs, bioc_pkgs)) {
  message(sprintf("  %s: %s", p, requireNamespace(p, quietly = TRUE)))
}
