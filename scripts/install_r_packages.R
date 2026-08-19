#!/usr/bin/env Rscript
# Install CRAN + Bioconductor packages used by tg_vs_ntc_deg_go.R

cran <- c("data.table", "ggplot2", "ggrepel")
bioc <- c("DESeq2", "clusterProfiler", "enrichplot", "org.Hs.eg.db")

install.packages(setdiff(cran, rownames(installed.packages())), repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE, ask = FALSE)
