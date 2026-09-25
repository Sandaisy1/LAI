#!/usr/bin/env Rscript
# =============================================================================
# 根据 UniProt 亚细胞定位，统计蛋白在细胞核 / 细胞质的分布并画圆饼图
#
# 输入：E:/R/Uniprot/protein.txt
#   每行一个 UniProt 登录号（如 P04637）、条目名（如 P53_HUMAN）或基因符号（如 TP53）
#   也可为带表头的 csv/tsv，默认取第一列
# 物种：基因符号默认按人（NCBI taxId 9606）检索已审核的 Swiss-Prot 条目
# 定位：只使用 UniProt 条目级 SUBCELLULAR LOCATION（不含各 isoform 的单独注释）
#
# 运行：
#   setwd("E:/R/Uniprot")
#   source("Uniprot_localization_pie.R")
# =============================================================================

options(stringsAsFactors = FALSE, timeout = 600)

cran_required <- c("ggplot2", "curl")

install_if_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) == 0) return(invisible(TRUE))
  install.packages(miss, repos = "https://cloud.r-project.org")
  still <- miss[!vapply(miss, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still) > 0) stop("缺少必需 R 包: ", paste(still, collapse = ", "))
  invisible(TRUE)
}

install_if_missing(cran_required)
suppressPackageStartupMessages(library(ggplot2))

# -----------------------------------------------------------------------------
# 路径
# -----------------------------------------------------------------------------
resolve_uniprot_dir <- function() {
  env_dir <- Sys.getenv("UNIPROT_DIR", unset = "")
  candidates <- c(
    env_dir,
    "E:/R/Uniprot",
    "E:\\R\\Uniprot",
    getwd()
  )
  candidates <- unique(candidates[nzchar(candidates)])
  for (d in candidates) {
    if (dir.exists(d) && file.exists(file.path(d, "protein.txt"))) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

uniprot_dir <- resolve_uniprot_dir()
protein_file <- file.path(uniprot_dir, "protein.txt")
out_dir <- file.path(uniprot_dir, "localization")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(protein_file)) {
  stop("找不到 protein.txt。请把文件放在 E:/R/Uniprot/protein.txt ，或先 setwd 到该目录。")
}

ORGANISM_ID <- 9606L
UNIPROT_FIELDS <- "accession,id,gene_primary,gene_names,protein_name,cc_subcellular_location"
CHUNK_SIZE <- 40L

message("读取: ", protein_file)

# -----------------------------------------------------------------------------
# 读入蛋白列表
# -----------------------------------------------------------------------------
read_protein_list <- function(path) {
  raw <- readLines(path, warn = FALSE, encoding = "UTF-8")
  raw <- trimws(raw)
  raw <- raw[!grepl("^\\s*#", raw)]
  raw <- raw[nzchar(raw)]
  if (length(raw) == 0) stop("protein.txt 是空的。")

  # 带表头或分隔符时取第一列
  split_one <- function(line) {
    parts <- unlist(strsplit(line, "[,;\\t]", perl = TRUE))
    trimws(parts[nzchar(trimws(parts))])[1]
  }
  ids <- vapply(raw, split_one, character(1), USE.NAMES = FALSE)
  ids <- ids[!is.na(ids) & nzchar(ids)]

  header_like <- grepl(
    "^(protein|gene|genesymbol|symbol|accession|uniprot|entry|id)$",
    ids[1],
    ignore.case = TRUE
  )
  if (isTRUE(header_like)) ids <- ids[-1]
  ids <- unique(ids)
  if (length(ids) == 0) stop("protein.txt 里没有可用的蛋白编号或基因名。")
  ids
}

is_accession <- function(x) {
  grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", x)
}

is_entry_name <- function(x) {
  grepl("^[A-Z0-9]{1,10}_[A-Z0-9]{1,5}$", x)
}

proteins <- read_protein_list(protein_file)
message("蛋白数: ", length(proteins))

# -----------------------------------------------------------------------------
# UniProt REST 查询（分批，避免 URL 过长）
# -----------------------------------------------------------------------------
# Windows 自带的 libcurl 常会协商 HTTP/2，随后报
# "Error in the HTTP2 framing layer"。这里强制 HTTP/1.1，失败再改用 wininet。
download_uniprot <- function(url) {
  tmp <- tempfile(fileext = ".tsv")
  last_err <- "未知错误"

  if (requireNamespace("curl", quietly = TRUE)) {
    for (attempt in 1:3) {
      handle <- curl::new_handle()
      curl::handle_setopt(
        handle,
        http_version = 2L,
        timeout = 180,
        followlocation = TRUE,
        useragent = "Uniprot_localization_pie.R"
      )
      ok <- tryCatch({
        curl::curl_download(url, tmp, handle = handle, quiet = TRUE)
        TRUE
      }, error = function(e) {
        last_err <<- conditionMessage(e)
        FALSE
      })
      if (isTRUE(ok) && file.exists(tmp) && file.info(tmp)$size > 0) return(tmp)
      Sys.sleep(1.5 * attempt)
    }
  }

  methods <- if (.Platform$OS.type == "windows") c("wininet", "libcurl") else "libcurl"
  for (method in methods) {
    ok <- tryCatch({
      utils::download.file(url, tmp, quiet = TRUE, mode = "wb", method = method)
      TRUE
    }, error = function(e) {
      last_err <<- conditionMessage(e)
      FALSE
    })
    if (isTRUE(ok) && file.exists(tmp) && file.info(tmp)$size > 0) return(tmp)
  }

  unlink(tmp)
  stop(
    "访问 UniProt 失败: ", last_err,
    "\n若仍是 HTTP2 framing layer，请换一个网络后重新 source 本脚本。"
  )
}

uniprot_search <- function(query) {
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/search?",
    "query=", utils::URLencode(query, reserved = TRUE),
    "&fields=", utils::URLencode(UNIPROT_FIELDS, reserved = TRUE),
    "&format=tsv&size=500"
  )
  tmp <- download_uniprot(url)
  on.exit(unlink(tmp), add = TRUE)
  utils::read.delim(tmp, sep = "\t", quote = "", comment.char = "", check.names = FALSE)
}

chunked_search <- function(ids, query_fun) {
  if (length(ids) == 0) return(data.frame())
  pieces <- list()
  n <- length(ids)
  starts <- seq(1L, n, by = CHUNK_SIZE)
  for (i in seq_along(starts)) {
    part <- ids[starts[i]:min(starts[i] + CHUNK_SIZE - 1L, n)]
    message(sprintf("  查询 %d/%d ...", i, length(starts)))
    pieces[[i]] <- query_fun(part)
    if (i < length(starts)) Sys.sleep(0.3)
  }
  do.call(rbind, pieces)
}

query_accessions <- function(acc) {
  q <- paste0("accession:(", paste(acc, collapse = " OR "), ")")
  uniprot_search(q)
}

query_entry_names <- function(names) {
  q <- paste0("id:(", paste(names, collapse = " OR "), ")")
  uniprot_search(q)
}

query_genes <- function(genes) {
  # 基因名里的特殊字符用引号包起来
  quoted <- paste0("gene_exact:", vapply(genes, function(g) {
    if (grepl("[^A-Za-z0-9_-]", g)) sprintf("\"%s\"", gsub("\"", "", g, fixed = TRUE)) else g
  }, character(1)))
  q <- paste0("(", paste(quoted, collapse = " OR "),
              ") AND organism_id:", ORGANISM_ID, " AND reviewed:true")
  uniprot_search(q)
}

acc_ids <- proteins[is_accession(proteins)]
entry_ids <- proteins[is_entry_name(proteins) & !is_accession(proteins)]
gene_ids <- proteins[!is_accession(proteins) & !is_entry_name(proteins)]

message("UniProt 检索: 登录号 ", length(acc_ids),
        "，条目名 ", length(entry_ids),
        "，基因符号 ", length(gene_ids))

hit_acc <- chunked_search(acc_ids, query_accessions)
hit_entry <- chunked_search(entry_ids, query_entry_names)
hit_gene <- chunked_search(gene_ids, query_genes)
hits <- rbind(hit_acc, hit_entry, hit_gene)

if (nrow(hits) == 0) {
  stop("UniProt 没有返回任何条目。请检查 protein.txt 是登录号还是基因符号，以及网络是否能访问 rest.uniprot.org。")
}

# 列名随 fields 固定，读入后按位置取，避免表头细微差别
colnames(hits)[1:6] <- c(
  "accession", "entry_name", "gene", "gene_names", "protein_name", "location_cc"
)

# UniProt 基因名在同一组里用空格分隔，不同基因用分号分隔
gene_tokens <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character(0))
  parts <- unlist(strsplit(x, "[;[:space:]]+"))
  toupper(parts[nzchar(parts)])
}

# 主基因名或别名命中都算；优先主基因名完全一致的条目
pick_gene_hit <- function(tab, gene) {
  if (nrow(tab) == 0) return(NULL)
  g <- toupper(gene)
  primary <- vapply(tab$gene, function(x) g %in% gene_tokens(x), logical(1))
  alias <- vapply(tab$gene_names, function(x) g %in% gene_tokens(x), logical(1))
  idx <- which(primary)
  if (length(idx) == 0) idx <- which(alias)
  if (length(idx) == 0) return(NULL)
  tab[idx[1], , drop = FALSE]
}

# -----------------------------------------------------------------------------
# 解析亚细胞定位，分成互斥类别（圆饼图各扇区相加为 100%）
# -----------------------------------------------------------------------------
extract_locations <- function(cc) {
  if (is.na(cc) || !nzchar(cc)) return(character(0))
  parts <- strsplit(cc, "SUBCELLULAR LOCATION:", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) return(character(0))
  block <- parts[[1]]
  block <- sub("\\s*Note=.*$", "", block)
  repeat {
    cleaned <- gsub("\\{[^{}]*\\}", "", block)
    if (identical(cleaned, block)) break
    block <- cleaned
  }
  block <- sub("^\\s*\\[[^\\]]*\\]:\\s*", "", block)
  locs <- unlist(strsplit(block, "[.;]"))
  locs <- trimws(locs)
  locs <- locs[nzchar(locs)]
  unique(locs)
}

is_nucleus_term <- function(loc) {
  grepl("^(Nucleus|Nucleoplasm|Nucleolus|Nuclear)\\b", loc, ignore.case = TRUE)
}

is_cytoplasm_term <- function(loc) {
  grepl("^(Cytoplasm|Cytosol|Cytoplasmic)\\b", loc, ignore.case = TRUE)
}

classify_locations <- function(locs) {
  in_nuc <- any(is_nucleus_term(locs))
  in_cyt <- any(is_cytoplasm_term(locs))
  if (in_nuc && in_cyt) return("细胞核和细胞质")
  if (in_nuc) return("细胞核")
  if (in_cyt) return("细胞质")
  if (length(locs) == 0) return("无定位信息")
  "其他细胞定位"
}

empty_row <- function(query) {
  data.frame(
    query = query,
    accession = NA_character_,
    entry_name = NA_character_,
    gene = NA_character_,
    protein_name = NA_character_,
    locations = NA_character_,
    in_nucleus = FALSE,
    in_cytoplasm = FALSE,
    category = "无定位信息",
    stringsAsFactors = FALSE
  )
}

rows <- vector("list", length(proteins))
for (i in seq_along(proteins)) {
  q <- proteins[i]
  if (is_accession(q)) {
    sub <- hits[hits$accession == q, , drop = FALSE]
    hit <- if (nrow(sub)) sub[1, , drop = FALSE] else NULL
  } else if (is_entry_name(q)) {
    sub <- hits[toupper(hits$entry_name) == toupper(q), , drop = FALSE]
    hit <- if (nrow(sub)) sub[1, , drop = FALSE] else NULL
  } else {
    hit <- pick_gene_hit(hits, q)
  }
  if (is.null(hit)) {
    rows[[i]] <- empty_row(q)
    next
  }
  locs <- extract_locations(hit$location_cc)
  in_nuc <- any(is_nucleus_term(locs))
  in_cyt <- any(is_cytoplasm_term(locs))
  rows[[i]] <- data.frame(
    query = q,
    accession = hit$accession,
    entry_name = hit$entry_name,
    gene = hit$gene,
    protein_name = hit$protein_name,
    locations = if (length(locs)) paste(locs, collapse = "; ") else NA_character_,
    in_nucleus = in_nuc,
    in_cytoplasm = in_cyt,
    category = classify_locations(locs),
    stringsAsFactors = FALSE
  )
}

result <- do.call(rbind, rows)
rownames(result) <- NULL

category_levels <- c("细胞核", "细胞质", "细胞核和细胞质", "其他细胞定位", "无定位信息")
result$category <- factor(result$category, levels = category_levels)

detail_file <- file.path(out_dir, "protein_localization.csv")
utils::write.csv(result, detail_file, row.names = FALSE, fileEncoding = "UTF-8")

summary_df <- as.data.frame(table(result$category), stringsAsFactors = FALSE)
colnames(summary_df) <- c("category", "n")
summary_df <- summary_df[summary_df$n > 0, , drop = FALSE]
summary_df$category <- factor(summary_df$category, levels = category_levels)
summary_df <- summary_df[order(summary_df$category), , drop = FALSE]
summary_df$percent <- summary_df$n / sum(summary_df$n) * 100
summary_df$label <- sprintf("%s\n%d (%.1f%%)", summary_df$category, summary_df$n, summary_df$percent)

n_nucleus <- sum(result$in_nucleus)
n_cytoplasm <- sum(result$in_cytoplasm)
n_total <- nrow(result)

summary_file <- file.path(out_dir, "localization_summary.csv")
utils::write.csv(summary_df[, c("category", "n", "percent")], summary_file,
                 row.names = FALSE, fileEncoding = "UTF-8")

message(sprintf("共 %d 个蛋白：注释到细胞核 %d 个，注释到细胞质 %d 个（同时定位的两边都计入）。",
                n_total, n_nucleus, n_cytoplasm))
print(summary_df[, c("category", "n", "percent")])

# -----------------------------------------------------------------------------
# 圆饼图
# -----------------------------------------------------------------------------
pie_colors <- c(
  "细胞核" = "#2F6FAD",
  "细胞质" = "#E07A3D",
  "细胞核和细胞质" = "#7B5EA7",
  "其他细胞定位" = "#6E9B6A",
  "无定位信息" = "#B0B3B8"
)

font_family <- "sans"
if (.Platform$OS.type == "windows") {
  win_fonts <- c("Microsoft YaHei", "SimHei", "SimSun")
  font_family <- win_fonts[1]
} else if (any(grepl("WenQuanYi Micro Hei", system("fc-list : family", intern = TRUE), fixed = TRUE))) {
  font_family <- "WenQuanYi Micro Hei"
}

subtitle <- sprintf(
  "共 %d 个蛋白；含细胞核注释 %d 个，含细胞质注释 %d 个",
  n_total, n_nucleus, n_cytoplasm
)

p <- ggplot(summary_df, aes(x = "", y = n, fill = category)) +
  geom_col(width = 1, color = "white", linewidth = 0.6) +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    family = font_family,
    size = 3.6,
    lineheight = 0.95
  ) +
  scale_fill_manual(values = pie_colors, drop = FALSE) +
  labs(
    title = "蛋白质亚细胞定位",
    subtitle = subtitle,
    fill = "定位"
  ) +
  theme_void(base_family = font_family) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 4)),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30", margin = margin(b = 8)),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.margin = margin(12, 12, 12, 12)
  )

pdf_file <- file.path(out_dir, "localization_pie.pdf")
png_file <- file.path(out_dir, "localization_pie.png")
ggsave(pdf_file, p, width = 8, height = 6, device = cairo_pdf)
ggsave(png_file, p, width = 8, height = 6, dpi = 180)

message("明细表: ", detail_file)
message("计数表: ", summary_file)
message("圆饼图: ", pdf_file)
message("圆饼图: ", png_file)
