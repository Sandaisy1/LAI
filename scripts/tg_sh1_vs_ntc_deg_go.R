#!/usr/bin/env Rscript
# Backward-compatible entry point. Use scripts/tg_vs_ntc_deg_go.R going forward.
this_dir <- tryCatch({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", fa[[1]])))
}, error = function(e) file.path("scripts"))
source(file.path(this_dir, "tg_vs_ntc_deg_go.R"))
