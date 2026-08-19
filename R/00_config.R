# TG knockdown RNA-seq (NTC / TG_sh1 / TG_sh5)
# Paths, group names, and the two DEG filter strategies.

CTRL_GROUP <- "NTC"
TREAT_GROUPS <- c("TG_sh1", "TG_sh5")
ALL_GROUPS <- c(CTRL_GROUP, TREAT_GROUPS)

# Linear fold-change cutoffs for upregulated genes (FC = treat / control).
# FC >= 1 means all genes with log2FC >= 0.
FC_CUTOFFS <- c(1, 1.25, 1.5, 2)

# Ranked upregulated lists (descending log2FC).
TOP_N <- c(50, 75, 100, 150, 200, 250, 300)

PSEUDOCOUNT <- 0.1
HEATMAP_MAX_GENES <- 100
ORA_MIN_GENES <- 10
GSEA_MIN_GENES <- 15
SPECIES_ORGDB <- "org.Hs.eg.db"
KEGG_ORG <- "hsa"
P_ADJUST <- "BH"
ENRICH_P_CUTOFF <- 0.05
ENRICH_Q_CUTOFF <- 0.20

EXCEL_PATTERN <- "shTG"
DIFF_FILENAME <- "gene_exp.diff"
FPKM_FILENAME <- "genes.fpkm_tracking"

default_data_dir <- function() {
  candidates <- c(
    Sys.getenv("TG_RNASEQ_DIR", unset = ""),
    "E:/R/TG_BRCA/TG",
    "E:\\R\\TG_BRCA\\TG",
    file.path(".", "data")
  )
  candidates <- candidates[nzchar(candidates)]
  for (p in candidates) {
    if (dir.exists(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  normalizePath("data", winslash = "/", mustWork = FALSE)
}

default_out_dir <- function(data_dir) {
  file.path(data_dir, "results")
}
