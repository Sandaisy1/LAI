#!/usr/bin/env Rscript
# =============================================================================
# TG BRCA：指定人源蛋白的共同序列（de novo motif）
#   - UniProt Swiss-Prot reviewed canonical 序列（每基因一条）
#   - ZOOPS EM，目标 3–5 条非冗余 motif
#   - 显著性：序列重排 empirical p 与 E-value（不伪造 p）
#   - 每条 motif 单独出 hits 表 + bits 序列 logo（ggseqlogo chemistry）
#
# 用法：
#   # 数据与结果默认在 E:\R\Protein
#   source("TG_protein_motif_pipeline.R")
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, timeout = 600)
Sys.setenv(LANGUAGE = "en")

# -----------------------------------------------------------------------------
# 0. 依赖
# -----------------------------------------------------------------------------
cran_required <- c(
  "dplyr", "tidyr", "tibble", "stringr",
  "ggplot2", "ggseqlogo", "writexl"
)
cran_optional <- c("httr", "jsonlite")

install_if_missing <- function(pkgs, required = TRUE) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) == 0) return(invisible(TRUE))
  tryCatch(
    install.packages(miss, repos = "https://cloud.r-project.org"),
    error = function(e) message("CRAN install failed: ", e$message)
  )
  still <- miss[!vapply(miss, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still) > 0 && required) {
    stop("缺少必需 R 包: ", paste(still, collapse = ", "))
  }
  invisible(TRUE)
}

install_if_missing(cran_required, required = TRUE)
install_if_missing(cran_optional, required = FALSE)
for (p in c(cran_required, cran_optional)) {
  if (requireNamespace(p, quietly = TRUE)) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}

# -----------------------------------------------------------------------------
# 1. 路径与参数
# -----------------------------------------------------------------------------
protein_data_candidates <- function() {
  env_dir <- Sys.getenv("TG_PROTEIN_DIR", unset = "")
  out <- c(env_dir, "E:/R/Protein", "E:\\R\\Protein")
  unique(out[nzchar(out)])
}

resolve_protein_dir <- function() {
  for (d in protein_data_candidates()) {
    if (dir.exists(d)) {
      return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  for (d in c("E:/R/Protein", "E:\\R\\Protein")) {
    if (dir.exists(dirname(d)) || dir.exists("E:/") || dir.exists("E:\\")) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      if (dir.exists(d)) {
        return(normalizePath(d, winslash = "/", mustWork = FALSE))
      }
    }
  }
  fallback <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  message("未找到 E:/R/Protein，回退到: ", fallback)
  fallback
}

find_gene_list_path <- function(protein_dir) {
  candidates <- c(
    file.path(protein_dir, "TG_protein_motif_genes.txt"),
    file.path(getwd(), "TG_protein_motif_genes.txt")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0) return(hit[[1]])
  file.path(protein_dir, "TG_protein_motif_genes.txt")
}

protein_dir <- resolve_protein_dir()
result_dir  <- file.path(protein_dir, "results")
log_dir     <- file.path(result_dir, "00_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("motif_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

AA20 <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
aa_index <- setNames(seq_along(AA20), AA20)

# 宽度、motif 条数、显著性（与 Cursor 规则一致）
motif_widths      <- c(8L, 10L, 12L, 15L, 21L)
n_motifs_target   <- 5L
n_motifs_min      <- 3L
p_cutoff          <- 0.01
n_em_iter         <- 25L
n_random_starts   <- 4L
n_kmer_seeds      <- 3L
# 经验 p = (1 + n_ge) / (1 + n_shuffle)，要能出现 p < 0.01，重排次数必须 >= 100
n_shuffle         <- 200L
min_sites         <- 4L
uniprot_pause_sec <- 0.2

default_genes <- c(
  "ARHGAP35", "EFNA1", "EFNA2", "EFNA3", "FAM20B", "ITGA2", "ITGAV",
  "LAMB1", "LAMC1", "MCAM", "RBP4", "RPL3", "RPL9", "RPL18A", "RPL32",
  "RPL35", "COL6A1", "RPL10A", "RPS11", "AGRN", "AP2M1", "RPS17",
  "RPS15", "RPS19", "RPS8", "RPL4", "RPL17", "RPL6", "RPL24", "RPL11",
  "RPL18", "RPL5", "RPL7", "RPL13", "RPL21", "RPL19", "RPL35A", "RPL29",
  "RPL22", "EFNB2", "RPS13", "RPS2", "RPS27", "RPS16", "RPS14", "RPS5"
)

# -----------------------------------------------------------------------------
# 2. 基因列表与序列 I/O
# -----------------------------------------------------------------------------
read_gene_list <- function(path) {
  if (!file.exists(path)) {
    log_msg("未找到基因表，使用脚本内置列表: ", path)
    return(unique(toupper(trimws(default_genes))))
  }
  raw <- readLines(path, warn = FALSE)
  raw <- gsub("#.*$", "", raw)
  genes <- unique(toupper(trimws(raw)))
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0) stop("基因表为空: ", path)
  genes
}

write_fasta <- function(seqs, path) {
  # seqs: named character
  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  for (nm in names(seqs)) {
    writeLines(paste0(">", nm), con)
    s <- seqs[[nm]]
    if (!nzchar(s)) next
    starts <- seq(1, nchar(s), by = 60)
    writeLines(substring(s, starts, pmin(starts + 59, nchar(s))), con)
  }
}

seq_to_int <- function(seq) {
  chars <- strsplit(toupper(seq), "")[[1]]
  chars[chars == "U"] <- "C"
  as.integer(aa_index[chars])
}

int_to_seq <- function(idx) {
  paste(ifelse(is.na(idx), "X", AA20[idx]), collapse = "")
}

# -----------------------------------------------------------------------------
# 3. UniProt Swiss-Prot（reviewed, 每基因一条 canonical）
# -----------------------------------------------------------------------------
uniprot_search <- function(query, fields = "accession,gene_primary,protein_name,length,sequence,reviewed") {
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/search?",
    "query=", utils::URLencode(query, reserved = TRUE),
    "&fields=", utils::URLencode(fields, reserved = TRUE),
    "&format=tsv&size=50"
  )
  txt <- tryCatch({
    if (requireNamespace("httr", quietly = TRUE)) {
      resp <- httr::GET(url, httr::timeout(60))
      code <- httr::status_code(resp)
      if (code == 429) {
        Sys.sleep(2)
        return(uniprot_search(query, fields))
      }
      if (code >= 400) stop("HTTP ", code)
      httr::content(resp, as = "text", encoding = "UTF-8")
    } else {
      con <- url(url, open = "rb")
      on.exit(close(con), add = TRUE)
      rawToChar(readBin(con, what = "raw", n = 2e6))
    }
  }, error = function(e) e)
  if (inherits(txt, "error")) return(list(ok = FALSE, error = txt$message, table = NULL))
  if (!nzchar(trimws(txt))) return(list(ok = TRUE, error = NULL, table = NULL))
  tab <- tryCatch(
    utils::read.delim(text = txt, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = ""),
    error = function(e) NULL
  )
  list(ok = TRUE, error = NULL, table = tab)
}

pick_canonical_row <- function(tab, symbol) {
  if (is.null(tab) || nrow(tab) == 0) return(NULL)
  names(tab) <- gsub("\\s+", "_", names(tab))
  gene_col <- grep("^[Gg]ene", names(tab), value = TRUE)
  if (length(gene_col) > 0) {
    primary <- toupper(as.character(tab[[gene_col[1]]]))
    exact <- grepl(paste0("(^|[ ;])", symbol, "($|[ ;])"), primary)
    if (any(exact)) tab <- tab[exact, , drop = FALSE]
  }
  if ("Length" %in% names(tab)) {
    tab <- tab[order(-as.numeric(tab$Length)), , drop = FALSE]
  }
  tab[1, , drop = FALSE]
}

fetch_one_uniprot <- function(symbol) {
  q1 <- paste0("gene_exact:", symbol, " AND organism_id:9606 AND reviewed:true")
  res <- uniprot_search(q1)
  tab <- res$table
  query_used <- q1
  if (is.null(tab) || nrow(tab) == 0) {
    q2 <- paste0("gene:", symbol, " AND organism_id:9606 AND reviewed:true")
    res <- uniprot_search(q2)
    tab <- res$table
    query_used <- q2
  }
  row <- pick_canonical_row(tab, symbol)
  if (is.null(row)) {
    return(data.frame(
      gene = symbol,
      uniprot = NA_character_,
      protein_name = NA_character_,
      length = NA_integer_,
      sequence = NA_character_,
      status = if (!is.null(res$error)) paste0("error:", res$error) else "not_found",
      query = query_used,
      stringsAsFactors = FALSE
    ))
  }
  seq_col <- intersect(c("Sequence", "sequence"), names(row))
  acc_col <- intersect(c("Entry", "accession"), names(row))
  name_col <- intersect(c("Protein_names", "protein_name"), names(row))
  len_col <- intersect(c("Length", "length"), names(row))
  seq <- if (length(seq_col)) as.character(row[[seq_col[1]]][1]) else NA_character_
  data.frame(
    gene = symbol,
    uniprot = if (length(acc_col)) as.character(row[[acc_col[1]]][1]) else NA_character_,
    protein_name = if (length(name_col)) as.character(row[[name_col[1]]][1]) else NA_character_,
    length = if (length(len_col)) as.integer(row[[len_col[1]]][1]) else nchar(seq),
    sequence = seq,
    status = if (nzchar(seq) && !is.na(seq)) "ok" else "empty_sequence",
    query = query_used,
    stringsAsFactors = FALSE
  )
}

fetch_uniprot_sequences <- function(symbols) {
  fasta_override <- file.path(protein_dir, "TG_protein_motif_sequences.fasta")
  rows <- vector("list", length(symbols))
  for (i in seq_along(symbols)) {
    sym <- symbols[i]
    log_msg("UniProt: ", sym, " (", i, "/", length(symbols), ")")
    rows[[i]] <- fetch_one_uniprot(sym)
    Sys.sleep(uniprot_pause_sec)
  }
  fetch_tbl <- dplyr::bind_rows(rows)
  if (file.exists(fasta_override)) {
    log_msg("发现本地 FASTA 覆盖: ", fasta_override)
    local <- read_simple_fasta(fasta_override)
    for (nm in names(local)) {
      gene <- toupper(sub("\\|.*$", "", nm))
      hit <- which(fetch_tbl$gene == gene)
      if (length(hit) == 1) {
        fetch_tbl$sequence[hit] <- local[[nm]]
        fetch_tbl$length[hit] <- nchar(local[[nm]])
        fetch_tbl$status[hit] <- "ok_local_fasta"
      }
    }
  }
  fetch_tbl
}

read_simple_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  headers <- grep("^>", lines)
  if (length(headers) == 0) return(character(0))
  out <- character(length(headers))
  names_out <- character(length(headers))
  for (i in seq_along(headers)) {
    names_out[i] <- sub("^>", "", lines[headers[i]])
    end <- if (i < length(headers)) headers[i + 1] - 1 else length(lines)
    seq_lines <- lines[(headers[i] + 1):end]
    seq_lines <- seq_lines[!grepl("^>", seq_lines)]
    out[i] <- paste(gsub("\\s+", "", seq_lines), collapse = "")
  }
  names(out) <- names_out
  out
}

# -----------------------------------------------------------------------------
# 4. 背景频率、PWM、窗口打分
# -----------------------------------------------------------------------------
count_aa <- function(int_list) {
  cnt <- integer(20)
  for (x in int_list) {
    ok <- x[!is.na(x)]
    if (length(ok)) {
      tab <- tabulate(ok, nbins = 20)
      cnt <- cnt + tab
    }
  }
  cnt
}

background_freq <- function(int_list) {
  cnt <- count_aa(int_list) + 1
  as.numeric(cnt / sum(cnt))
}

pwm_from_counts <- function(counts, bg, pseudo = 1) {
  # counts: 20 x w
  den <- colSums(counts) + pseudo
  sweep(counts + (pseudo * bg), 2, den, "/")
}

log_odds_from_pwm <- function(pwm, bg) {
  log((pwm + 1e-12) / (bg + 1e-12))
}

score_windows <- function(seq_int, log_odds) {
  w <- ncol(log_odds)
  n <- length(seq_int)
  if (n < w) return(numeric(0))
  nwin <- n - w + 1
  scores <- numeric(nwin)
  for (j in seq_len(w)) {
    aa <- seq_int[j:(j + nwin - 1)]
    col <- log_odds[, j]
    add <- col[aa]
    add[is.na(add)] <- -20
    scores <- scores + add
  }
  scores
}

information_content <- function(pwm, bg) {
  sum(pwm * log2((pwm + 1e-12) / (bg + 1e-12)))
}

consensus_from_pwm <- function(pwm, min_p = 0.4) {
  paste(vapply(seq_len(ncol(pwm)), function(j) {
    a <- which.max(pwm[, j])
    if (pwm[a, j] >= min_p) AA20[a] else "x"
  }, character(1)), collapse = "")
}

# -----------------------------------------------------------------------------
# 5. ZOOPS EM（每条序列 0 或 1 个位点）
# -----------------------------------------------------------------------------
best_window_starts <- function(int_list, width) {
  # 用最常见 8-mer（或 width）做种子位置
  k <- min(8L, width)
  bucket <- vector("list", length(int_list))
  for (idx in seq_along(int_list)) {
    s <- int_list[[idx]]
    if (length(s) < k) next
    nwin <- length(s) - k + 1
    keys <- character(nwin)
    for (i in seq_len(nwin)) {
      win <- s[i:(i + k - 1)]
      keys[i] <- if (any(is.na(win))) NA_character_ else paste(win, collapse = ",")
    }
    bucket[[idx]] <- keys
  }
  kmers <- unlist(bucket, use.names = FALSE)
  kmers <- kmers[!is.na(kmers)]
  if (length(kmers) == 0) return(list())
  top <- names(sort(table(kmers), decreasing = TRUE))
  top <- head(top, n_kmer_seeds)
  seeds <- list()
  for (key in top) {
    want <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
    pos <- integer(length(int_list))
    for (i in seq_along(int_list)) {
      s <- int_list[[i]]
      pos[i] <- 0L
      if (length(s) < k) next
      for (j in seq_len(length(s) - k + 1)) {
        if (j + width - 1 > length(s)) break
        if (all(!is.na(s[j:(j + k - 1)]) & s[j:(j + k - 1)] == want)) {
          pos[i] <- as.integer(j)
          break
        }
      }
    }
    seeds[[length(seeds) + 1]] <- pos
  }
  seeds
}

random_starts <- function(int_list, width, n) {
  lapply(seq_len(n), function(dummy) {
    vapply(int_list, function(s) {
      nwin <- length(s) - width + 1
      if (nwin < 1) return(0L)
      as.integer(sample.int(nwin, 1))
    }, integer(1))
  })
}

counts_from_sites <- function(int_list, starts, width) {
  counts <- matrix(0, nrow = 20, ncol = width)
  for (i in seq_along(int_list)) {
    st <- starts[i]
    if (is.na(st) || st < 1) next
    if (st + width - 1 > length(int_list[[i]])) next
    win <- int_list[[i]][st:(st + width - 1)]
    for (j in seq_len(width)) {
      a <- win[j]
      if (!is.na(a)) counts[a, j] <- counts[a, j] + 1
    }
  }
  counts
}

em_from_start <- function(int_list, starts, width, bg, n_iter) {
  n <- length(int_list)
  lambda <- max(0.5, mean(starts > 0))
  best_llr <- -Inf
  best <- NULL
  prev_starts <- starts
  for (iter in seq_len(n_iter)) {
    counts <- counts_from_sites(int_list, starts, width)
    n_sites <- sum(starts > 0)
    if (n_sites < 2) break
    pwm <- pwm_from_counts(counts, bg)
    lod <- log_odds_from_pwm(pwm, bg)
    new_starts <- integer(n)
    post_has <- numeric(n)
    total_llr <- 0
    for (i in seq_len(n)) {
      sc <- score_windows(int_list[[i]], lod)
      if (length(sc) == 0) {
        new_starts[i] <- 0L
        next
      }
      # ZOOPS：位点先验 lambda/m，无位点 1-lambda
      m <- length(sc)
      log_site <- log(pmax(lambda, 1e-6) / m) + sc
      log_none <- log(pmax(1 - lambda, 1e-6))
      mx <- max(c(log_site, log_none))
      w_site <- exp(log_site - mx)
      w_none <- exp(log_none - mx)
      denom <- sum(w_site) + w_none
      z <- w_site / denom
      post_has[i] <- sum(z)
      j <- which.max(z)
      none_post <- w_none / denom
      new_starts[i] <- if (z[j] >= none_post && sc[j] > 0) as.integer(j) else 0L
      if (new_starts[i] > 0) total_llr <- total_llr + sc[j]
    }
    if (total_llr > best_llr) {
      best_llr <- total_llr
      best <- list(starts = new_starts, pwm = pwm, lod = lod, llr = total_llr, lambda = lambda)
    }
    if (iter > 3 && identical(prev_starts, new_starts)) break
    prev_starts <- new_starts
    starts <- new_starts
    lambda <- max(0.05, min(0.95, mean(post_has)))
  }
  best
}

extract_sites <- function(seq_df, starts, width) {
  n <- nrow(seq_df)
  out <- vector("list", n)
  k <- 0
  for (i in seq_len(n)) {
    st <- starts[i]
    if (is.na(st) || st < 1) next
    s <- seq_df$sequence[i]
    if (st + width - 1 > nchar(s)) next
    k <- k + 1
    out[[k]] <- data.frame(
      gene = seq_df$gene[i],
      uniprot = seq_df$uniprot[i],
      start = st,
      end = st + width - 1,
      site_sequence = substring(s, st, st + width - 1),
      stringsAsFactors = FALSE
    )
  }
  if (k == 0) return(NULL)
  dplyr::bind_rows(out[seq_len(k)])
}

discover_one_motif <- function(seq_df, int_list, width, bg) {
  starts_list <- c(
    best_window_starts(int_list, width),
    random_starts(int_list, width, n_random_starts)
  )
  best <- NULL
  for (st in starts_list) {
    if (length(st) != length(int_list)) next
    fit <- em_from_start(int_list, st, width, bg, n_em_iter)
    if (is.null(fit)) next
    if (is.null(best) || fit$llr > best$llr) best <- fit
  }
  if (is.null(best)) return(NULL)
  hits <- extract_sites(seq_df, best$starts, width)
  if (is.null(hits) || nrow(hits) < min_sites) return(NULL)
  lod <- best$lod
  llr <- numeric(nrow(hits))
  for (i in seq_len(nrow(hits))) {
    gi <- match(hits$gene[i], seq_df$gene)
    sc <- score_windows(int_list[[gi]], lod)
    llr[i] <- sc[hits$start[i]]
  }
  hits$llr_score <- llr
  hits <- hits[order(-hits$llr_score), , drop = FALSE]
  list(
    width = width,
    pwm = best$pwm,
    lod = lod,
    starts = best$starts,
    llr = sum(hits$llr_score),
    ic = information_content(best$pwm, bg),
    consensus = consensus_from_pwm(best$pwm),
    hits = hits,
    lambda = best$lambda
  )
}

# -----------------------------------------------------------------------------
# 6. 显著性：对序列做组成保留重排，不伪造 p
# -----------------------------------------------------------------------------
shuffle_seq <- function(seq) {
  paste(sample(strsplit(seq, "")[[1]]), collapse = "")
}

motif_score_on_seqs <- function(int_list, lod) {
  total <- 0
  n_hit <- 0
  for (s in int_list) {
    sc <- score_windows(s, lod)
    if (length(sc) == 0) next
    mx <- max(sc)
    if (mx > 0) {
      total <- total + mx
      n_hit <- n_hit + 1
    }
  }
  c(llr = total, n_hit = n_hit)
}

empirical_significance <- function(seq_df, motif, bg, n_shuffle) {
  obs <- motif_score_on_seqs(lapply(seq_df$sequence, seq_to_int), motif$lod)
  null_llr <- numeric(n_shuffle)
  for (b in seq_len(n_shuffle)) {
    shuf <- vapply(seq_df$sequence, shuffle_seq, character(1))
    shuf_int <- lapply(unname(shuf), seq_to_int)
    null_llr[b] <- motif_score_on_seqs(shuf_int, motif$lod)[["llr"]]
  }
  n_ge <- sum(null_llr >= obs[["llr"]])
  p <- (1 + n_ge) / (1 + n_shuffle)
  list(
    observed_llr = unname(obs[["llr"]]),
    n_hits = unname(obs[["n_hit"]]),
    empirical_p = p,
    null_mean = mean(null_llr),
    null_sd = stats::sd(null_llr)
  )
}

pwm_correlation <- function(a, b) {
  w <- min(ncol(a), ncol(b))
  if (w < 4) return(0)
  # 对较短 motif 在较长 motif 上滑窗，取最大相关
  long <- if (ncol(a) >= ncol(b)) a else b
  short <- if (ncol(a) >= ncol(b)) b else a
  best <- -1
  for (off in seq_len(ncol(long) - ncol(short) + 1)) {
    x <- as.vector(long[, off:(off + ncol(short) - 1)])
    y <- as.vector(short)
    r <- suppressWarnings(stats::cor(x, y))
    if (!is.na(r) && r > best) best <- r
  }
  best
}

# -----------------------------------------------------------------------------
# 7. 作图与写出
# -----------------------------------------------------------------------------
safe_ggsave <- function(plot, stub) {
  tryCatch(
    ggplot2::ggsave(paste0(stub, ".pdf"), plot, width = 9, height = 3.4),
    error = function(e) log_msg("pdf ggsave failed: ", e$message)
  )
  tryCatch(
    ggplot2::ggsave(paste0(stub, ".png"), plot, width = 9, height = 3.4, dpi = 200),
    error = function(e) log_msg("png ggsave failed: ", e$message)
  )
}

plot_seqlogo <- function(site_seqs, title) {
  ggseqlogo::ggseqlogo(site_seqs, method = "bits", col_scheme = "chemistry") +
    ggplot2::scale_y_continuous(name = "bits", limits = c(0, 4), breaks = 0:4, expand = c(0, 0.05)) +
    ggplot2::labs(title = title, x = "position") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11, hjust = 0),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6))
    )
}

write_pwm_csv <- function(pwm, path) {
  tab <- as.data.frame(t(pwm))
  names(tab) <- AA20
  tab$position <- seq_len(nrow(tab))
  tab <- tab[, c("position", AA20)]
  utils::write.csv(tab, path, row.names = FALSE)
}

mask_sites <- function(int_list, starts, width) {
  out <- lapply(int_list, identity)
  for (i in seq_along(out)) {
    st <- starts[i]
    if (is.na(st) || st < 1) next
    en <- min(length(out[[i]]), st + width - 1)
    out[[i]][st:en] <- NA_integer_
  }
  out
}

# -----------------------------------------------------------------------------
# 8. 主流程
# -----------------------------------------------------------------------------
run_protein_motif_pipeline <- function() {
  gene_file <- find_gene_list_path(protein_dir)
  genes <- read_gene_list(gene_file)
  log_msg("Protein data dir: ", protein_dir)
  log_msg("Gene list: ", gene_file)
  log_msg("Genes: ", length(genes), " -> ", paste(genes, collapse = ", "))

  fetch_tbl <- fetch_uniprot_sequences(genes)
  utils::write.csv(fetch_tbl[, setdiff(names(fetch_tbl), "sequence")],
                   file.path(result_dir, "00_sequence_fetch_log.csv"),
                   row.names = FALSE)

  ok <- fetch_tbl[fetch_tbl$status %in% c("ok", "ok_local_fasta") & !is.na(fetch_tbl$sequence), , drop = FALSE]
  ok <- ok[!duplicated(ok$gene), , drop = FALSE]
  if (nrow(ok) < min_sites) {
    stop("可用蛋白序列不足 ", min_sites, " 条（成功 ", nrow(ok), "）。见 00_sequence_fetch_log.csv")
  }
  log_msg("Sequences OK: ", nrow(ok), " / ", nrow(fetch_tbl))

  fasta_names <- paste0(ok$gene, "|", ok$uniprot)
  fasta <- ok$sequence
  names(fasta) <- fasta_names
  write_fasta(fasta, file.path(result_dir, "00_sequences.fasta"))

  int_list <- lapply(ok$sequence, seq_to_int)
  names(int_list) <- ok$gene
  bg <- background_freq(int_list)
  n_tests <- length(motif_widths) * n_motifs_target

  found <- list()
  work_int <- int_list
  for (k in seq_len(n_motifs_target)) {
    log_msg("Searching motif ", k, " / ", n_motifs_target)
    cand <- list()
    for (w in motif_widths) {
      log_msg("  width=", w)
      fit <- tryCatch(
        discover_one_motif(ok, work_int, w, bg),
        error = function(e) {
          log_msg("  EM failed width=", w, ": ", e$message)
          NULL
        }
      )
      if (!is.null(fit)) cand[[length(cand) + 1]] <- fit
    }
    if (length(cand) == 0) {
      log_msg("No more motifs at round ", k)
      break
    }
    ic_vals <- vapply(cand, function(x) x$ic, numeric(1))
    pick <- cand[[which.max(ic_vals)]]

    # 与已有 motif 去冗余
    redundant <- FALSE
    if (length(found) > 0) {
      cors <- vapply(found, function(m) pwm_correlation(m$pwm, pick$pwm), numeric(1))
      if (any(cors > 0.85)) {
        log_msg("  skip redundant motif (PWM cor=", round(max(cors), 3), ")")
        work_int <- mask_sites(work_int, pick$starts, pick$width)
        next
      }
    }

    sig <- empirical_significance(ok, pick, bg, n_shuffle)
    pick$empirical_p <- sig$empirical_p
    pick$e_value <- min(1, sig$empirical_p * n_tests)
    pick$observed_llr <- sig$observed_llr
    pick$null_mean <- sig$null_mean
    pick$significant <- is.finite(pick$empirical_p) && pick$empirical_p < p_cutoff
    log_msg(
      "  keep consensus=", pick$consensus,
      " sites=", nrow(pick$hits),
      " IC=", round(pick$ic, 2),
      " p=", signif(pick$empirical_p, 3),
      " E=", signif(pick$e_value, 3)
    )
    found[[length(found) + 1]] <- pick
    work_int <- mask_sites(work_int, pick$starts, pick$width)
  }

  if (length(found) == 0) {
    log_msg("未发现 motif。请检查序列是否拉到，或放宽宽度。")
    return(invisible(NULL))
  }

  # 只保留 p < 0.01；若不足 3 条则如实少报（不凑数、不改 p）
  keep <- found[vapply(found, function(m) isTRUE(m$significant), logical(1))]
  if (length(keep) == 0) {
    log_msg("没有任何 motif 达到 p < ", p_cutoff, "；仍写出全部候选并标明不显著")
    keep <- found
  } else if (length(keep) < n_motifs_min) {
    log_msg("显著 motif 仅 ", length(keep), " 条（目标 3–5），不凑数")
  } else if (length(keep) > n_motifs_target) {
    keep <- keep[seq_len(n_motifs_target)]
  }

  summary_rows <- list()
  all_hits <- list()
  for (i in seq_along(keep)) {
    m <- keep[[i]]
    tag <- paste0("motif", i)
    out_d <- file.path(result_dir, tag)
    dir.create(out_d, recursive = TRUE, showWarnings = FALSE)

    hits <- m$hits
    hits$motif_id <- tag
    hits$consensus <- m$consensus
    hits$motif_empirical_p <- m$empirical_p
    hits$motif_e_value <- m$e_value
    hits$genome_set_n <- nrow(ok)
    utils::write.csv(hits, file.path(out_d, paste0(tag, "_hits.csv")), row.names = FALSE)
    write_pwm_csv(m$pwm, file.path(out_d, paste0(tag, "_pwm.csv")))

    site_fa <- hits$site_sequence
    names(site_fa) <- paste0(hits$gene, "/", hits$uniprot, "/", hits$start, "-", hits$end)
    write_fasta(site_fa, file.path(out_d, paste0(tag, "_aligned_sites.fasta")))

    title <- sprintf(
      "%s  %s  n=%d/%d  IC=%.2f bits  p=%.3g  E=%.3g",
      tag, m$consensus, nrow(hits), nrow(ok), m$ic, m$empirical_p, m$e_value
    )
    logo <- tryCatch(
      plot_seqlogo(hits$site_sequence, title),
      error = function(e) {
        log_msg("seqlogo failed ", tag, ": ", e$message)
        NULL
      }
    )
    if (!is.null(logo)) {
      safe_ggsave(logo, file.path(out_d, paste0(tag, "_seqlogo")))
    }

    summary_rows[[i]] <- data.frame(
      motif_id = tag,
      width = m$width,
      consensus = m$consensus,
      n_sites = nrow(hits),
      n_sequences = nrow(ok),
      coverage = nrow(hits) / nrow(ok),
      total_ic_bits = m$ic,
      observed_llr = m$observed_llr,
      empirical_p = m$empirical_p,
      e_value = m$e_value,
      significant_p_lt_0.01 = isTRUE(m$significant),
      stringsAsFactors = FALSE
    )
    all_hits[[i]] <- hits
  }

  summary_tbl <- dplyr::bind_rows(summary_rows)
  hits_tbl <- dplyr::bind_rows(all_hits)
  utils::write.csv(summary_tbl, file.path(result_dir, "motif_significance_summary.csv"), row.names = FALSE)
  utils::write.csv(hits_tbl, file.path(result_dir, "motif_all_hits.csv"), row.names = FALSE)
  writexl::write_xlsx(
    list(summary = summary_tbl, hits = hits_tbl, fetch = fetch_tbl[, setdiff(names(fetch_tbl), "sequence")]),
    file.path(result_dir, "motif_significance_summary.xlsx")
  )

  # 总览 logo（只含显著或实际写出的 motif）
  if (length(keep) > 1 && requireNamespace("ggseqlogo", quietly = TRUE)) {
    logo_list <- lapply(seq_along(keep), function(i) {
      m <- keep[[i]]
      plot_seqlogo(
        m$hits$site_sequence,
        sprintf("motif%d  %s  p=%.3g", i, m$consensus, m$empirical_p)
      )
    })
    combo <- tryCatch(
      cowplot_or_patchwork(logo_list),
      error = function(e) NULL
    )
    if (is.null(combo) && requireNamespace("patchwork", quietly = TRUE)) {
      combo <- Reduce(`+`, logo_list) + patchwork::plot_layout(ncol = 1)
    }
    if (!is.null(combo)) {
      h <- 3.1 * length(logo_list)
      tryCatch(
        ggplot2::ggsave(file.path(result_dir, "motif_all_seqlogos.pdf"), combo, width = 9, height = h),
        error = function(e) log_msg("combo pdf failed: ", e$message)
      )
      tryCatch(
        ggplot2::ggsave(file.path(result_dir, "motif_all_seqlogos.png"), combo, width = 9, height = h, dpi = 200),
        error = function(e) log_msg("combo png failed: ", e$message)
      )
    }
  }

  log_msg("Done. Motifs written: ", nrow(summary_tbl), " -> ", result_dir)
  print(summary_tbl)
  invisible(list(summary = summary_tbl, hits = hits_tbl, fetch = fetch_tbl))
}

cowplot_or_patchwork <- function(plots) {
  if (requireNamespace("cowplot", quietly = TRUE)) {
    return(cowplot::plot_grid(plotlist = plots, ncol = 1, align = "v"))
  }
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(Reduce(`+`, plots) + patchwork::plot_layout(ncol = 1))
  }
  plots[[1]]
}

run_motif_selftest <- function() {
  set.seed(35)
  planted <- "CADCQEGGGC"
  flank <- function() paste(sample(AA20, 40, replace = TRUE), collapse = "")
  seqs <- vapply(seq_len(12), function(i) paste0(flank(), planted, flank()), character(1))
  seq_df <- data.frame(
    gene = paste0("G", seq_len(12)),
    uniprot = paste0("P", seq_len(12)),
    sequence = seqs,
    stringsAsFactors = FALSE
  )
  int_list <- lapply(seq_df$sequence, seq_to_int)
  bg <- background_freq(int_list)
  fit <- discover_one_motif(seq_df, int_list, 10L, bg)
  if (is.null(fit)) stop("selftest: planted motif not recovered")
  sig <- empirical_significance(seq_df, fit, bg, 30L)
  recovered <- grepl("C", fit$consensus) && grepl("G", fit$consensus)
  log_msg(
    "SELFTEST consensus=", fit$consensus,
    " sites=", nrow(fit$hits),
    " p=", signif(sig$empirical_p, 3),
    " recovered_CG=", recovered
  )
  if (!recovered) stop("selftest: consensus lost planted C/G")
  if (nrow(fit$hits) < 8) stop("selftest: too few planted sites")
  if (!(is.finite(sig$empirical_p) && sig$empirical_p < 0.05)) {
    stop("selftest: planted motif not significant")
  }
  log_msg("SELFTEST passed")
  invisible(TRUE)
}

if (identical(Sys.getenv("TG_MOTIF_SELFTEST"), "1")) {
  run_motif_selftest()
} else {
  run_protein_motif_pipeline()
}
