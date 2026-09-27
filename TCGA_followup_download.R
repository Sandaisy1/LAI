################################################################################
# TCGA_followup_download.R
# 只下载「随访进展」临床表，不跑分析。
# 放到 E:/R/TCGA-BRCA-2026 后 Source；或先 Sys.setenv(TCGA_BRCA_2026_DIR=...)
#
# 1) Liu 2018 TCGA-CDR（Xena）：PFI / DFI，是后来进展，不是就诊时 M1
# 2) GDC biotab：follow_up + new tumor event（同一批 BRCA，标注更细）
################################################################################

dest <- Sys.getenv("TCGA_BRCA_2026_DIR", unset = "E:/R/TCGA-BRCA-2026")
if (!dir.exists(dest)) dir.create(dest, recursive = TRUE, showWarnings = FALSE)
message("下载目录：", dest)

files <- list(
  list(
    name = "Survival_SupplementalTable_S1_20171025_xena_sp",
    url = "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp",
    note = "Liu 2018 / Xena：PFI DFI DSS OS + new_tumor_event_type/site"
  ),
  list(
    name = "TCGA-CDR-SupplementalTableS1.xlsx",
    url = "https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81",
    note = "同一张表的官方 xlsx（GDC 出版物）"
  ),
  list(
    name = "nationwidechildrens.org_clinical_follow_up_v4.0_brca.txt",
    url = "https://api.gdc.cancer.gov/data/62d4515f-a30b-4b1a-b2dd-c8bf9476e803",
    note = "GDC follow-up v4.0：new_tumor_event_dx_indicator"
  ),
  list(
    name = "nationwidechildrens.org_clinical_nte_brca.txt",
    url = "https://api.gdc.cancer.gov/data/a88c168e-4bba-4bd2-9c0c-77934444cc1c",
    note = "GDC new tumor event：类型 + 远处复发部位"
  ),
  list(
    name = "nationwidechildrens.org_clinical_follow_up_v4.0_nte_brca.txt",
    url = "https://api.gdc.cancer.gov/data/a9baecf0-5549-4396-8805-a6d1681d11cd",
    note = "GDC follow-up v4.0 NTE 补充"
  ),
  list(
    name = "nationwidechildrens.org_clinical_follow_up_v2.1_brca.txt",
    url = "https://api.gdc.cancer.gov/data/403b5cef-8173-47c7-b56a-cc94dcfbb2e3",
    note = "GDC follow-up v2.1（较早随访表，可与 v4 合并）"
  ),
  list(
    name = "nationwidechildrens.org_clinical_patient_brca.txt",
    url = "https://api.gdc.cancer.gov/data/8162d394-8b64-4da2-9f5b-d164c54b9608",
    note = "GDC patient 临床主表"
  )
)

message("浏览器也可逐个打开：")
for (f in files) message("  ", f$name, "\n    ", f$url)

ok <- 0L
for (f in files) {
  path <- file.path(dest, f$name)
  message("下载 ", f$name, " …")
  got <- tryCatch({
    utils::download.file(f$url, destfile = path, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    message("  失败：", conditionMessage(e))
    FALSE
  })
  if (isTRUE(got) && file.exists(path) && file.info(path)$size > 1000) {
    message("  OK  ", round(file.info(path)$size / 1024), " KB  ", f$note)
    ok <- ok + 1L
  } else {
    message("  请用浏览器打开上面的链接手动保存")
  }
}

message("完成 ", ok, "/", length(files), "。浏览器用的页面：")
message("  Xena：https://xenabrowser.net/datapages/?dataset=Survival_SupplementalTable_S1_20171025_xena_sp&host=https://pancanatlas.xenahubs.net")
message("  GDC 论文：https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018")
message("  GDC 项目 Clinical 按钮：https://portal.gdc.cancer.gov/projects/TCGA-BRCA")
invisible(ok)
