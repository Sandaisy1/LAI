# 病人临床表格式

把表达矩阵和临床表放到 `E:/R/TG_BRCA/TG/`（或设置环境变量），再运行：

```r
setwd("E:/R/TG_BRCA/TG")
source("TG_RNAseq_patient_metastasis_neural.R")
```

## 必需文件

- `patient_expression.csv`：第一列基因符号，其余列为样品表达值（count / FPKM / TPM 均可）
- `patient_clinical.csv`：每行一个样品，必须能对应表达矩阵的列名

也可用环境变量指定路径：

```r
Sys.setenv(
  TG_PATIENT_EXPR = "E:/R/TG_BRCA/TG/patient_expression.csv",
  TG_PATIENT_CLINICAL = "E:/R/TG_BRCA/TG/patient_clinical.csv"
)
```

## 临床列（列名可近似匹配）

| 用途 | 可识别列名示例 |
|---|---|
| 样品 ID | `sample`, `sample_id`, `barcode`, `patient` |
| 转移 | `metastasis`, `pathologic_M`, `ajcc_pathologic_m`, `distant_metastasis`, `M` |
| 生存时间 | `os_time`, `OS.time`, `days_to_death`, `days_to_last_follow_up` |
| 生存结局 | `os_event`, `OS`, `vital_status`, `fustat`（1/Dead = 事件） |
| 预后分组（可选） | `prognosis`（`good` / `poor`） |

转移分组规则：

- `M1` / `yes` / `metastasis` → **Metastasis**
- `M0` / `no` / `non-metastasis` → **NoMetastasis**
- `MX` / `unknown` → 排除，不进入转移比较

本目录的两个 CSV 只用于演示列格式和单元测试，不是真实病人数据。
