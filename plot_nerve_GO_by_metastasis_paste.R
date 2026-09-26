## 整段复制到 RStudio。不要再跑旧的：
##   ann[, distant_M := factor(meta_M, ...)]
##   ann[meta_M == "M1" | stage_simplified == "Stage IV", ...]
## 你的 clinical.tsv 没有 meta_M / stage_simplified，那几行一定会报错。
##
## 本段会下载并 Source 最新的独立脚本（扫描 phenotype、条形码 06/07、PFI 回退）。

if (dir.exists("E:/R/BRCA")) setwd("E:/R/BRCA")

tf <- "nerve_GO_by_metastasis.R"
need_new <- !file.exists(tf) ||
  !any(grepl("try_download_phenotype|sample_type_code",
             readLines(tf, warn = FALSE, encoding = "UTF-8"), fixed = FALSE))
if (need_new) {
  url <- "https://github.com/Sandaisy1/LAI/raw/cursor/tcga-brca-go-individual-analysis-3a7a/nerve_GO_by_metastasis.R"
  message("下载最新独立脚本：", url)
  utils::download.file(url, tf, mode = "wb", quiet = FALSE)
}
source(tf, encoding = "UTF-8")
