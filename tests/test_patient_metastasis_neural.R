#!/usr/bin/env Rscript
# 病人转移 × 神经通路：不依赖 Bioconductor 的单元测试

options(tg.patient.skip_main = TRUE, stringsAsFactors = FALSE, warn = 1)
src <- file.path("TG_RNAseq_patient_metastasis_neural.R")
if (!file.exists(src)) src <- file.path("..", "TG_RNAseq_patient_metastasis_neural.R")
if (!file.exists(src)) stop("找不到 TG_RNAseq_patient_metastasis_neural.R")
sys.source(src, envir = .GlobalEnv, keep.source = TRUE)

fail <- 0L
check <- function(ok, msg) {
  if (isTRUE(ok)) {
    cat("PASS:", msg, "\n")
  } else {
    fail <<- fail + 1L
    cat("FAIL:", msg, "\n")
  }
}

ex_dir <- file.path("examples", "patient_metastasis_neural")
if (!dir.exists(ex_dir)) ex_dir <- file.path("..", "examples", "patient_metastasis_neural")

clin <- read_table_flex(file.path(ex_dir, "patient_clinical.csv"))
prep <- prepare_clinical(clin)
check(sum(prep$metastasis_group == "Metastasis", na.rm = TRUE) == 9, "M1/yes labeled Metastasis")
check(sum(prep$metastasis_group == "NoMetastasis", na.rm = TRUE) == 10, "M0/no labeled NoMetastasis")
check(sum(is.na(prep$metastasis_group)) == 1, "MX/unknown excluded")
check(sum(prep$os_event == 1, na.rm = TRUE) >= 1, "survival events parsed")
check(sum(prep$prognosis_group == "PoorPrognosis", na.rm = TRUE) >= 1, "poor prognosis parsed")

mat <- read_expression_matrix(file.path(ex_dir, "patient_expression.csv"))
aligned <- align_patient_data(mat, clin)
check(ncol(aligned$mat) == 20, "all example samples aligned")
check(all(colnames(aligned$mat) == aligned$clinical$sample), "sample order matches")

si <- aligned$clinical
si_de <- si[!is.na(si$metastasis_group), ]
mat_de <- aligned$mat[, si_de$sample, drop = FALSE]
norm <- patient_normalize(mat_de, si_de, "fpkm")
de <- patient_de_metastasis(norm$log_mat, si_de, raw_mat = mat_de, value_type = "fpkm")
check(all(c("gene", "log2FC", "pvalue", "padj") %in% names(de)), "DE table columns")
check(any(is.finite(de$pvalue)), "p-values are estimated, not fabricated NA")

sets <- build_neural_gene_sets_symbols(rownames(norm$log_mat))
check(length(sets) >= 3, "neural gene sets map to example matrix")
scores <- pathway_mean_z(norm$log_mat, sets)
tab <- compare_scores_by_group(scores, si_de$metastasis_group)
core <- tab[tab$pathway == "NEURAL_ALL_core", ]
check(nrow(core) == 1 && is.finite(core$wilcox_p[1]), "core neural score has Wilcoxon p")
check(core$mean_Metastasis[1] > core$mean_NoMetastasis[1], "planted neural score higher in metastasis")

up <- select_up_fc_p01(de, 1)
check(all(up$pvalue < 0.01), "FC subset uses p < 0.01")
check(all(up$log2FC > 0), "only upregulated genes kept")
top <- select_up_topn_p01(de, 50)
check(nrow(top) <= 50, "topN cap")
check(all(top$pvalue < 0.01 | nrow(top) == 0), "topN also filtered by p < 0.01")

# 中文转移标签
zh <- data.frame(
  sample = c("A", "B", "C"),
  metastasis = c("有转移", "无转移", "未转移"),
  stringsAsFactors = FALSE
)
zhp <- prepare_clinical(zh)
check(identical(zhp$metastasis_group, c("Metastasis", "NoMetastasis", "NoMetastasis")),
      "Chinese metastasis labels")

files <- find_patient_files(normalizePath(ex_dir))
check(!is.null(files) && file.exists(files$expr) && file.exists(files$clinical),
      "example files discovered")

tmp <- file.path(tempdir(), "tg_patient_neural_test")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE)
out <- run_patient_metastasis_neural(project_dir = normalizePath(ex_dir), result_dir = tmp)
check(file.exists(file.path(tmp, "DE_full_metastasis_vs_nometastasis.csv")), "DE table written")
check(file.exists(file.path(tmp, "neural_pathway_Metastasis_vs_NoMetastasis.csv")), "pathway table written")
check(dir.exists(file.path(tmp, "Metastasis_vs_NoMetastasis", "FoldChange", "FC_1.5")), "FC_1.5 folder")
check(dir.exists(file.path(tmp, "Metastasis_vs_NoMetastasis", "TopRank", "top75")), "top75 folder")
check(dir.exists(file.path(tmp, "Metastasis_vs_NoMetastasis", "Focused_neural")), "Focused_neural folder")
fc15 <- file.path(tmp, "Metastasis_vs_NoMetastasis", "FoldChange", "FC_1.5")
check(length(list.files(fc15)) > 0, "FC_1.5 is not an empty placeholder")

if (fail > 0) {
  cat("FAILED checks:", fail, "\n")
  quit(status = 1)
}
cat("All patient neural / metastasis tests passed\n")
