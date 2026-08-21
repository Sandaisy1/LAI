# =============================================================================
# ZZX DIA-NN：根据上清肽段丰度推断 ALDH1A1（P00352）的细胞剪切（不是胰酶）
# 在 RStudio / Cursor 的 R Console 里运行（Windows 请指定 UTF-8）：
#   source("ZZX_ALB_cleavage.R", encoding = "UTF-8")
# 或先改数据目录：
#   zzx_data_dir <- "C:/Users/Lenovo/Desktop/ZZX"
#   source("ZZX_ALB_cleavage.R", encoding = "UTF-8")
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

if (!exists("zzx_data_dir")) {
  zzx_data_dir <- "C:/Users/Lenovo/Desktop/ZZX"
}
if (!exists("zzx_n_last")) {
  zzx_n_last <- 2L
}
if (!exists("zzx_alb_auto_run")) {
  zzx_alb_auto_run <- TRUE
}

# UniProt P00352 AL1A1_HUMAN SV=2（人 ALDH1A1，501 aa；胞质酶，无信号肽）
TARGET_UNIPROT <- "P00352"
TARGET_GENE <- "ALDH1A1"
TARGET_GENE_SYNONYMS <- c("ALDH1")
TARGET_NAME <- "AL1A1_HUMAN"
TARGET_SEQUENCE <- paste0(
  "MSSSGTPDLPVLLTDLKIQYTKIFINNEWHDSVSGKKFPVFNPATEEELCQVEEGDKEDV",
  "DKAVKAARQAFQIGSPWRTMDASERGRLLYKLADLIERDRLLLATMESMNGGKLYSNAYL",
  "NDLAGCIKTLRYCAGWADKIQGRTIPIDGNFFTYTRHEPIGVCGQIIPWNFPLVMLIWKI",
  "GPALSCGNTVVVKPAEQTPLTALHVASLIKEAGFPPGVVNIVPGYGPTAGAAISSHMDID",
  "KVAFTGSTEVGKLIKEAAGKSNLKRVTLELGGKSPCIVLADADLDNAVEFAHHGVFYHQG",
  "QCCIAASRIFVEESIYDEFVRRSVERAKKYILGNPLTPGVTQGPQIDKEQYDKILDLIES",
  "GKKEGAKLECGGGPWGNKGYFVQPTVFSNVTDEMRIAKEEIFGPVQQIMKFKSLDDVIKR",
  "ANNTFYGLSAGVFTKDIDKAITISSALQAGTVWVNCYGVVSAQCPFGGFKMSGNGRELGE",
  "YGFHEYTEVKTVTVKISQKNS"
)
if (nchar(TARGET_SEQUENCE) != 501L) {
  stop("ALDH1A1 序列长度应为 501，实际 ", nchar(TARGET_SEQUENCE))
}

META_COLS <- c(
  "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
  "First.Protein.Description", "Proteotypic", "Stripped.Sequence",
  "Modified.Sequence", "Precursor.Charge", "Precursor.Id",
  "Precursor.Quantity", "Q.Value", "Protein.Q.Value", "PG.Q.Value",
  "Global.Q.Value", "Lib.Q.Value", "Proteotypic.Ids"
)
META_SUFFIXES <- c("Sequence", "Ids", "Names", "Genes", "Description", "Group")
REPORT_NAMES <- c(
  "report.pr_matrix", "report.pr_matrix.tsv", "report.pr_matrix.txt",
  "report.pr.matrix", "report.pr.matri", "report.pr.matri.tsv"
)

zzx_log <- function(logs, msg) {
  message(msg)
  logs[[length(logs) + 1L]] <- msg
  logs
}

aa_at <- function(protein, i) {
  if (i < 1L || i > nchar(protein)) {
    return("x")
  }
  substr(protein, i, i)
}

parse_num <- function(x) {
  x <- trimws(as.character(x))
  bad <- x %in% c("", "NA", "NaN", "nan", "None", "NULL", "#N/A") | is.na(x)
  x[bad] <- NA_character_
  x <- gsub(",", "", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

tokens <- function(value) {
  text <- gsub("[;|,/]", " ", as.character(value))
  text <- trimws(text)
  if (!nzchar(text) || is.na(text)) {
    return(character())
  }
  toupper(strsplit(text, "[[:space:]]+")[[1]])
}

is_target_row <- function(genes, ids, names, group) {
  g <- tokens(genes)
  id <- sub("-.*$", "", tokens(ids))
  nm <- toupper(as.character(names))
  gp <- toupper(as.character(group))
  TARGET_GENE %in% g ||
    any(TARGET_GENE_SYNONYMS %in% g) ||
    TARGET_UNIPROT %in% id ||
    grepl(TARGET_UNIPROT, gp, fixed = TRUE) ||
    grepl(TARGET_NAME, nm, fixed = TRUE)
}

find_report <- function(data_dir = NULL) {
  dirs <- c(
    data_dir,
    "C:/Users/Lenovo/Desktop/ZZX",
    "C:/User/Lenovo/Desktop/ZZX",
    file.path(getwd(), "ZZX"),
    getwd()
  )
  dirs <- unique(dirs[!is.na(dirs) & nzchar(dirs)])
  for (folder in dirs) {
    if (!dir.exists(folder)) {
      next
    }
    for (nm in REPORT_NAMES) {
      p <- file.path(folder, nm)
      if (file.exists(p)) {
        return(normalizePath(p, winslash = "/", mustWork = TRUE))
      }
    }
    hits <- list.files(folder, pattern = "^report\\.pr", full.names = TRUE)
    hits <- hits[!grepl("pg_matrix", basename(hits), ignore.case = TRUE)]
    if (length(hits) > 0L) {
      return(normalizePath(hits[[1]], winslash = "/", mustWork = TRUE))
    }
  }
  stop(
    "未找到 report.pr_matrix。请把 DIA-NN 肽段表放到 C:/Users/Lenovo/Desktop/ZZX，",
    "或先设置 zzx_data_dir 再 source。"
  )
}

sample_columns <- function(fieldnames, n_last = 2L) {
  extras <- fieldnames[!fieldnames %in% META_COLS]
  ends_meta <- vapply(extras, function(col) {
    any(vapply(META_SUFFIXES, function(sfx) endsWith(col, sfx), logical(1)))
  }, logical(1))
  numeric_like <- extras[!ends_meta]
  if (length(numeric_like) >= n_last) {
    return(tail(numeric_like, n_last))
  }
  if (length(numeric_like) > 0L) {
    return(numeric_like)
  }
  if (length(extras) > 0L) {
    return(tail(extras, n_last))
  }
  character()
}

short_sample_name <- function(col, index) {
  text <- gsub("\\", "/", as.character(col), fixed = TRUE)
  base <- basename(text)
  lower <- tolower(base)
  for (ext in c(".raw", ".mzml", ".dia", ".parquet", ".wiff2", ".wiff")) {
    if (endsWith(lower, ext)) {
      base <- substr(base, 1L, nchar(base) - nchar(ext))
      break
    }
  }
  base <- trimws(base)
  if (!nzchar(base)) {
    paste0("supernatant_", index)
  } else {
    base
  }
}

map_peptide <- function(protein, peptide) {
  pep <- toupper(trimws(peptide))
  alb <- toupper(protein)
  if (!nzchar(pep)) {
    return(data.frame(start = integer(), end = integer()))
  }
  loc <- gregexpr(pep, alb, fixed = TRUE)[[1]]
  if (loc[1] == -1L) {
    loc <- gregexpr(gsub("I", "L", pep), gsub("I", "L", alb), fixed = TRUE)[[1]]
  }
  if (loc[1] == -1L) {
    return(data.frame(start = integer(), end = integer()))
  }
  data.frame(start = as.integer(loc), end = as.integer(loc + nchar(pep) - 1L))
}

tryptic_n <- function(protein, start, block_proline = TRUE) {
  if (start <= 1L) {
    return(TRUE)
  }
  prev <- aa_at(protein, start - 1L)
  this <- aa_at(protein, start)
  if (prev %in% c("K", "R")) {
    return(!(block_proline && identical(this, "P")))
  }
  FALSE
}

tryptic_c <- function(protein, end, block_proline = TRUE) {
  if (end >= nchar(protein)) {
    return(TRUE)
  }
  aa <- aa_at(protein, end)
  nxt <- aa_at(protein, end + 1L)
  if (aa %in% c("K", "R")) {
    return(!(block_proline && identical(nxt, "P")))
  }
  FALSE
}

terminus_class <- function(n_ok, c_ok) {
  if (n_ok && c_ok) {
    return("fully_tryptic")
  }
  if (n_ok && !c_ok) {
    return("semi_tryptic_neoC")
  }
  if (!n_ok && c_ok) {
    return("semi_tryptic_neoN")
  }
  "non_tryptic"
}

motif <- function(protein, after) {
  chars <- vapply((after - 3L):(after + 4L), function(k) aa_at(protein, k), character(1))
  paste0(paste(chars[1:4], collapse = ""), "|", paste(chars[5:8], collapse = ""))
}

protease_hint <- function(protein, after) {
  n <- nchar(protein)
  if (after < 1L || after >= n) {
    return("unknown")
  }
  if (after == 1L) {
    return("initiator methionine aminopeptidase (MetAP)")
  }
  p4 <- aa_at(protein, after - 3L)
  p1 <- aa_at(protein, after)
  p1p <- aa_at(protein, after + 1L)
  if (identical(p1, "R") && identical(p4, "R")) {
    return("furin-like (RXXR)")
  }
  if (identical(p1, "R") && p4 %in% c("K", "R")) {
    return("proprotein convertase-like")
  }
  if (identical(p1, "D") && p1p %in% c("G", "A", "S")) {
    return("caspase-like (D|x)")
  }
  if (p1 %in% strsplit("ASGC", "")[[1]] && !identical(p1p, "P")) {
    return("signal-peptidase-like small P1")
  }
  if (p1 %in% c("K", "R")) {
    return("trypsin-like / K-or-R P1 (could be cellular or look tryptic after digest)")
  }
  if (p1 %in% strsplit("FLYM", "")[[1]] && !identical(p1p, "P")) {
    return("hydrophobic P1 (cathepsin/signal-like)")
  }
  "unassigned cellular protease"
}

mean_or_na <- function(values) {
  values <- values[!is.na(values)]
  if (!length(values)) {
    return(NA_real_)
  }
  mean(values)
}

region_mean <- function(profile, start, end) {
  if (end < start) {
    return(NA_real_)
  }
  mean_or_na(profile[start:end])
}

fmt_num <- function(x) {
  if (length(x) == 0L || is.na(x)) {
    return("")
  }
  format(x, digits = 6, scientific = FALSE, trim = TRUE)
}

svg_escape <- function(text) {
  text <- gsub("&", "&amp;", as.character(text), fixed = TRUE)
  text <- gsub("<", "&lt;", text, fixed = TRUE)
  text <- gsub(">", "&gt;", text, fixed = TRUE)
  text <- gsub("\"", "&quot;", text, fixed = TRUE)
  text
}

write_table <- function(path, df) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

write_svg <- function(path, lines) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = FALSE)
}

write_processing_svg <- function(path, protein, candidates, n_term_cov) {
  n <- nchar(protein)
  x_left <- function(res) {
    if (res <= 1L) 40 else 90 + 850 * (res - 2) / max(n - 1, 1)
  }
  x_right <- function(res) {
    if (res <= 1L) 90 else 90 + 850 * (res - 1) / max(n - 1, 1)
  }
  y0 <- 80
  bar_h <- 36
  lines <- c(
    '<svg xmlns="http://www.w3.org/2000/svg" width="980" height="220">',
    "<style>text{font-family:Arial,Helvetica,sans-serif;font-size:12px}</style>",
    '<text x="40" y="28" font-size="16">ALDH1A1 (P00352) cellular processing schematic</text>',
    paste0('<text x="40" y="48" fill="#555">no signal peptide; UniProt chain 2-', n, " after initiator Met removal</text>"),
    sprintf('<rect x="%.1f" y="%s" width="%.1f" height="%s" fill="#f4c7c3" stroke="#333"/>', x_left(1), y0, x_right(1) - x_left(1), bar_h),
    sprintf('<rect x="%.1f" y="%s" width="%.1f" height="%s" fill="#b7d7b0" stroke="#333"/>', x_left(2), y0, x_right(n) - x_left(2), bar_h),
    sprintf('<text x="%.1f" y="%s" font-size="11">M1</text>', x_left(1) + 4, y0 + 22),
    sprintf('<text x="%.1f" y="%s" font-size="11">CHAIN 2-%s</text>', x_left(2) + 8, y0 + 22, n)
  )
  shown <- integer()
  for (i in seq_len(nrow(candidates))) {
    after <- as.integer(candidates$after_residue[i])
    if (after %in% shown || after < 1L || after >= n) {
      next
    }
    shown <- c(shown, after)
    color <- if (identical(as.character(candidates$type[i]), "initiator_met")) "#c0392b" else "#6c3483"
    x <- x_right(after)
    lines <- c(
      lines,
      sprintf('<line x1="%.1f" y1="%s" x2="%.1f" y2="%s" stroke="%s" stroke-width="2"/>', x, y0 - 8, x, y0 + bar_h + 18, color),
      sprintf('<text x="%.1f" y="%s" font-size="10" text-anchor="middle" fill="%s">%s|%s</text>', x, y0 + bar_h + 32, color, after, after + 1L)
    )
  }
  cov <- sprintf(
    "coverage M1=%s  chain2-80=%s  rest81-%s=%s",
    fmt_num(n_term_cov$initiator), fmt_num(n_term_cov$mature_n), n, fmt_num(n_term_cov$mature_rest)
  )
  lines <- c(lines, sprintf('<text x="40" y="204" fill="#333">%s</text>', svg_escape(cov)), "</svg>")
  write_svg(path, lines)
}

write_full_length_svg <- function(path, protein, candidates, n_term_cov) {
  n <- nchar(protein)
  internals <- candidates[candidates$type == "internal" & candidates$source %in% c("neo_N", "neo_C"), , drop = FALSE]
  if (nrow(internals) > 0L) {
    internals <- internals[order(internals$after_residue), , drop = FALSE]
    internals <- internals[!duplicated(internals$after_residue) & internals$after_residue > 1 & internals$after_residue < n, , drop = FALSE]
    if (nrow(internals) > 4L) {
      internals <- internals[seq_len(4L), , drop = FALSE]
    }
  }
  n_steps <- 2L + as.integer(nrow(internals) > 0L)
  left <- 130
  bar_w <- 920
  bar_h <- 42
  zoom_w <- 70
  rest_w <- bar_w - zoom_w
  step_gap <- 152
  top <- 96
  height <- 86 + n_steps * step_gap + 56
  if (nrow(internals) > 0L) {
    height <- height + 36
  }
  x_left <- function(res) {
    if (res <= 1L) left else left + zoom_w + rest_w * (res - 2) / max(n - 1, 1)
  }
  x_right <- function(res) {
    if (res <= 1L) left + zoom_w else left + zoom_w + rest_w * (res - 1) / max(n - 1, 1)
  }
  fmt_cov <- function(v) {
    if (is.na(v)) "supernatant: not detected" else sprintf("supernatant: detected (%.3g)", v)
  }
  lines <- c(
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="%s" viewBox="0 0 1120 %s">', height, height),
    "<rect width='100%' height='100%' fill='#ffffff'/>",
    "<style>text{font-family:Microsoft YaHei,SimHei,Arial,Helvetica,sans-serif}</style>",
    '<text x="40" y="28" font-size="18" font-weight="700">Full-length ALDH1A1 cleavage  全长蛋白如何剪切</text>',
    '<text x="40" y="50" font-size="12" fill="#555">P00352 cytosolic ALDH1A1 (1-501). No signal peptide. Cellular processing only (not trypsin).</text>'
  )
  draw_domains <- function(y, keep_from, keep_to, ghost = FALSE) {
    domains <- list(
      c(1, 1, "#e74c3c", "#fadbd8", "M1", "1"),
      c(2, n, "#1e8449", "#d5f5e3", "CHAIN", paste0("2-", n))
    )
    out <- character()
    for (d in domains) {
      a <- as.integer(d[[1]])
      b <- as.integer(d[[2]])
      lo <- max(a, keep_from)
      hi <- min(b, keep_to)
      if (lo > hi) {
        next
      }
      x1 <- x_left(lo)
      x2 <- x_right(hi)
      w <- max(x2 - x1, 2)
      op <- if (ghost) 0.22 else 1
      dash <- if (ghost) ' stroke-dasharray="5 3"' else ""
      out <- c(out, sprintf(
        '<rect x="%.1f" y="%s" width="%.1f" height="%s" rx="6" fill="%s" stroke="%s" stroke-width="1.6"%s opacity="%s"/>',
        x1, y, w, bar_h, d[[4]], d[[3]], dash, op
      ))
      if (w >= 28) {
        lab <- if (lo == a && hi == b) d[[6]] else paste0(lo, "-", hi)
        out <- c(
          out,
          sprintf('<text x="%.1f" y="%s" text-anchor="middle" font-size="11" font-weight="700" fill="%s">%s</text>', x1 + w / 2, y + 17, d[[3]], d[[5]]),
          sprintf('<text x="%.1f" y="%s" text-anchor="middle" font-size="10" fill="#333">%s</text>', x1 + w / 2, y + 32, lab)
        )
      }
    }
    out
  }
  draw_scissors <- function(after, y, color, label, motif_txt) {
    x <- x_right(after)
    c(
      sprintf('<line x1="%.1f" y1="%s" x2="%.1f" y2="%s" stroke="%s" stroke-width="2.2"/>', x, y - 16, x, y + bar_h + 8, color),
      sprintf('<polygon points="%.1f,%s %.1f,%s %.1f,%s" fill="%s"/>', x - 7, y - 18, x + 7, y - 18, x, y - 4, color),
      sprintf('<text x="%.1f" y="%s" text-anchor="middle" font-size="11" font-weight="700" fill="%s">%s</text>', x, y + bar_h + 22, color, svg_escape(label)),
      sprintf('<text x="%.1f" y="%s" text-anchor="middle" font-size="10" fill="#666">%s</text>', x, y + bar_h + 36, svg_escape(motif_txt))
    )
  }
  step_badge <- function(num, y, title, subtitle) {
    cy <- y + bar_h / 2
    c(
      sprintf('<circle cx="48" cy="%.1f" r="16" fill="#2c3e50"/>', cy),
      sprintf('<text x="48" y="%.1f" text-anchor="middle" font-size="14" font-weight="700" fill="#fff">%s</text>', cy + 5, num),
      sprintf('<text x="130" y="%s" font-size="13" font-weight="700">%s</text>', y - 22, svg_escape(title)),
      sprintf('<text x="130" y="%s" font-size="11" fill="#555">%s</text>', y - 6, svg_escape(subtitle))
    )
  }
  down_arrow <- function(y_from, caption) {
    x <- left + bar_w / 2
    y1 <- y_from + bar_h + 40
    y2 <- y1 + 22
    c(
      sprintf('<line x1="%.1f" y1="%s" x2="%.1f" y2="%s" stroke="#2c3e50" stroke-width="2"/>', x, y1, x, y2),
      sprintf('<polygon points="%.1f,%s %.1f,%s %.1f,%s" fill="#2c3e50"/>', x - 6, y2, x + 6, y2, x, y2 + 10),
      sprintf('<text x="%.1f" y="%s" font-size="11" fill="#2c3e50">%s</text>', x + 12, y1 + 16, svg_escape(caption))
    )
  }
  y1 <- top
  lines <- c(
    lines,
    step_badge(1, y1, "Full-length ALDH1A1  全长蛋白", paste0("residues 1-", n, " (cytosolic; no signal peptide)")),
    draw_domains(y1, 1, n),
    draw_scissors(1, y1, "#c0392b", "cut 1|2", "xxxM|SSSG  initiator MetAP"),
    down_arrow(y1, "remove initiator Met 1")
  )
  y2 <- y1 + step_gap
  lines <- c(
    lines,
    step_badge(2, y2, "After Met removal  切掉起始 Met", paste0("released M1 (", fmt_cov(n_term_cov$initiator), "); UniProt chain 2-", n, " (", fmt_cov(n_term_cov$mature_n), ")")),
    draw_domains(y2, 1, 1, TRUE),
    draw_domains(y2, 2, n),
    sprintf('<rect x="%.1f" y="%s" width="%.1f" height="%s" rx="8" fill="none" stroke="#1e8449" stroke-width="2.4"/>', x_left(2), y2 - 4, x_right(n) - x_left(2), bar_h + 8)
  )
  if (nrow(internals) > 0L) {
    y3 <- y2 + step_gap
    cuts <- paste(sprintf("%s|%s", internals$after_residue, internals$before_residue), collapse = ", ")
    lines <- c(lines, step_badge(3, y3, "Additional cellular cuts  内部剪切", paste0("neo-N / neo-C (not trypsin): ", cuts)), draw_domains(y3, 2, n))
    for (i in seq_len(nrow(internals))) {
      after <- as.integer(internals$after_residue[i])
      lines <- c(lines, draw_scissors(after, y3, "#6c3483", paste0("cut ", after, "|", after + 1L), paste(internals$motif_P4_P4prime[i], internals$protease_hint[i])))
    }
  }
  ly <- height - 28
  lines <- c(
    lines,
    sprintf('<rect x="%s" y="%s" width="16" height="12" fill="#e74c3c" opacity="0.35" stroke="#e74c3c"/>', left, ly - 12),
    sprintf('<text x="%s" y="%s" font-size="11">initiator Met 1</text>', left + 22, ly),
    sprintf('<rect x="%s" y="%s" width="16" height="12" fill="#1e8449" opacity="0.35" stroke="#1e8449"/>', left + 200, ly - 12),
    sprintf('<text x="%s" y="%s" font-size="11">CHAIN 2-501 (UniProt)</text>', left + 222, ly),
    "</svg>"
  )
  write_svg(path, lines)
}


write_coverage_svg <- function(path, protein, pep, sample_names, profiles) {
  n <- nchar(protein)
  width <- 1000
  track_h <- 140
  pep_h <- max(80, 8 + 7 * min(nrow(pep), 40))
  height <- 70 + track_h + pep_h
  x0 <- 50
  bar_w <- 920
  x_at <- function(res) x0 + bar_w * (res - 0.5) / n
  ymax <- 1
  for (nm in names(profiles)) {
    vals <- profiles[[nm]]
    vals <- vals[!is.na(vals) & vals > 0]
    if (length(vals)) {
      ymax <- max(ymax, log10(vals))
    }
  }
  ymax <- ymax * 1.08
  lines <- c(
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%s" height="%s">', width, height),
    "<style>text{font-family:Arial,Helvetica,sans-serif;font-size:11px}</style>",
    '<text x="50" y="24" font-size="16">ALDH1A1 residue abundance (cell supernatant) and peptide map</text>',
    '<text x="50" y="42" fill="#555">grey=trypsin; orange=cellular neo-N; blue=cellular neo-C; red=both</text>',
    sprintf('<rect x="%s" y="55" width="%s" height="%s" fill="#fafafa" stroke="#ccc"/>', x0, bar_w, track_h)
  )
  for (pair in list(c(1, 1, "#f4c7c3"), c(2, n, "#eaf6e8"))) {
    xa <- x_at(as.integer(pair[[1]])) - bar_w / n / 2
    xb <- x_at(as.integer(pair[[2]])) + bar_w / n / 2
    lines <- c(lines, sprintf('<rect x="%.1f" y="55" width="%.1f" height="%s" fill="%s" opacity="0.45"/>', xa, xb - xa, track_h, pair[[3]]))
  }
  colors <- c("#1f77b4", "#d62728", "#2ca02c")
  series_names <- c(sample_names, "mean")
  for (i in seq_along(series_names)) {
    nm <- series_names[[i]]
    if (!nm %in% names(profiles)) {
      next
    }
    pts <- character()
    for (res in seq_len(n)) {
      v <- profiles[[nm]][res]
      if (is.na(v) || v <= 0) {
        next
      }
      y <- 55 + track_h - (log10(v) / ymax) * (track_h - 8)
      pts <- c(pts, sprintf("%.1f,%.1f", x_at(res), y))
    }
    if (length(pts)) {
      col <- colors[((i - 1L) %% length(colors)) + 1L]
      lines <- c(
        lines,
        sprintf('<polyline fill="none" stroke="%s" stroke-width="1.5" points="%s"/>', col, paste(pts, collapse = " ")),
        sprintf('<text x="%s" y="%s" fill="%s">%s</text>', x0 + 8, 68 + (i - 1) * 14, col, svg_escape(nm))
      )
    }
  }
  class_color <- c(
    fully_tryptic = "#bbbbbb",
    semi_tryptic_neoN = "#e67e22",
    semi_tryptic_neoC = "#2980b9",
    non_tryptic = "#c0392b"
  )
  shown <- pep[seq_len(min(nrow(pep), 80L)), , drop = FALSE]
  y_pep0 <- 55 + track_h + 16
  if (nrow(shown) > 0L) {
    for (i in seq_len(nrow(shown))) {
      y <- y_pep0 + ((i - 1L) %% 40L) * 6
      x1 <- x_at(shown$start[i])
      x2 <- x_at(shown$end[i])
      col <- class_color[[shown$terminus_class[i]]]
      if (is.null(col)) {
        col <- "#888888"
      }
      lines <- c(lines, sprintf(
        '<rect x="%.1f" y="%.1f" width="%.1f" height="4" fill="%s" opacity="0.85"/>',
        min(x1, x2), y, max(2, abs(x2 - x1)), col
      ))
    }
  }
  lines <- c(lines, sprintf('<text x="%s" y="%s" fill="#555">residue 1 ... %s (ALDH1A1 P00352)</text>', x0, height - 12, n), "</svg>")
  write_svg(path, lines)
}

cov_label <- function(v) {
  if (!is.na(v) && v > 0) "detected" else "not_detected"
}

infer_uniprot <- function(profile, neo_n, neo_c, protein) {
  n <- nchar(protein)
  met <- region_mean(profile, 1, 1)
  mat <- region_mean(profile, 2, min(80L, n))
  n2 <- neo_n[["2"]]
  c1 <- neo_c[["1"]]
  met_ok <- identical(cov_label(met), "not_detected") && identical(cov_label(mat), "detected")
  neo_ok <- (!is.null(n2) && n2$n > 0) || (!is.null(c1) && c1$n > 0)
  data.frame(
    after_residue = 1L, before_residue = 2L, type = "initiator_met", source = "uniprot_annotated",
    supported = if (met_ok || neo_ok) "yes" else "possible",
    motif_P4_P4prime = motif(protein, 1L),
    protease_hint = "initiator methionine aminopeptidase (MetAP)",
    n_neo_N_peptides = if (is.null(n2)) 0 else n2$n,
    abundance_neo_N = if (is.null(n2)) NA_real_ else n2$abundance,
    n_neo_C_peptides = if (is.null(c1)) 0 else c1$n,
    abundance_neo_C = if (is.null(c1)) NA_real_ else c1$abundance,
    abundance_left_window = met,
    abundance_right_window = mat,
    log2_step = NA_real_,
    evidence_zh = paste0(
      "UniProt 无信号肽/propeptide。注释加工是切除起始 Met，成熟链 2-", n, "。",
      "残基1 ", cov_label(met), "，残基2-80 ", cov_label(mat), "。",
      if (neo_ok) " 有 neo-N@2 或 neo-C@1，直接支持 Met 切除。" else " 无 neo 末端时，仅靠覆盖不能单独定论。",
      " 上清中检出胞质 ALDH1A1 可能来自裂解、泄漏或囊泡，内部 neo 末端才是额外细胞剪切证据。"
    ),
    note_zh = "Chain 2-501；不要把 ALDH1A1 当成有信号肽的分泌蛋白",
    stringsAsFactors = FALSE
  )
}

collect_neo <- function(pep, which) {
  bucket <- list()
  if (!nrow(pep)) {
    return(bucket)
  }
  keep <- if (identical(which, "N")) pep$cellular_neo_N else pep$cellular_neo_C
  sub <- pep[keep, , drop = FALSE]
  if (!nrow(sub)) {
    return(bucket)
  }
  for (i in seq_len(nrow(sub))) {
    pos <- if (identical(which, "N")) sub$start[i] else sub$end[i]
    key <- as.character(pos)
    if (is.null(bucket[[key]])) {
      bucket[[key]] <- list(n = 0, abundance = 0, peptides = character())
    }
    bucket[[key]]$n <- bucket[[key]]$n + 1L
    add <- sub$abundance_mean[i]
    if (!is.na(add)) {
      bucket[[key]]$abundance <- bucket[[key]]$abundance + add
    }
    bucket[[key]]$peptides <- c(bucket[[key]]$peptides, sub$Stripped.Sequence[i])
  }
  bucket
}

candidate_row <- function(after, before, type, source, supported, protein, n_neo_n, ab_n, n_neo_c, ab_c, left, right, log2_step, evidence, note) {
  data.frame(
    after_residue = as.integer(after),
    before_residue = as.integer(before),
    type = type,
    source = source,
    supported = supported,
    motif_P4_P4prime = motif(protein, as.integer(after)),
    protease_hint = protease_hint(protein, as.integer(after)),
    n_neo_N_peptides = as.integer(n_neo_n),
    abundance_neo_N = ab_n,
    n_neo_C_peptides = as.integer(n_neo_c),
    abundance_neo_C = ab_c,
    abundance_left_window = left,
    abundance_right_window = right,
    log2_step = log2_step,
    evidence_zh = evidence,
    note_zh = note,
    stringsAsFactors = FALSE
  )
}

abundance_steps <- function(profile, protein, window = 15L) {
  n <- length(profile)
  out <- list()
  if (n < (2L * window + 1L)) {
    return(out)
  }
  for (after in window:(n - window)) {
    left <- mean_or_na(profile[(after - window + 1L):after])
    right <- mean_or_na(profile[(after + 1L):(after + window)])
    if (is.na(left) || is.na(right) || (left <= 0 && right <= 0)) {
      next
    }
    lstep <- log2((right + 1) / (left + 1))
    if (abs(lstep) < 1.5) {
      next
    }
    out[[length(out) + 1L]] <- candidate_row(
      after, after + 1L, "internal", "abundance_step", "weak", protein,
      0, NA_real_, 0, NA_real_, left, right, lstep,
      sprintf("切割位点两侧 %s aa 丰度 log2 变化 %.2f。空洞也可能来自肽段过短/疏水/漏切，不能单独当作细胞剪切。", window, lstep),
      "丰度阶跃只作提示，优先看 neo-N/neo-C"
    )
  }
  if (!length(out)) {
    return(out)
  }
  steps <- do.call(rbind, out)
  steps <- steps[order(-abs(steps$log2_step)), , drop = FALSE]
  head(steps, 15L)
}

zzx_alb_run <- function(data_dir = zzx_data_dir, n_last = zzx_n_last, outdir = NULL) {
  logs <- list()
  report <- find_report(data_dir)
  logs <- zzx_log(logs, paste("[data] 读取", report))
  raw <- tryCatch(
    utils::read.delim(report, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "", fileEncoding = "UTF-8"),
    error = function(e) {
      utils::read.delim(report, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
    }
  )
  qty_cols <- sample_columns(names(raw), n_last = n_last)
  if (!length(qty_cols)) {
    stop("没有找到样品丰度列（表头 Precursor.Id 之后的上清强度列）")
  }
  sample_names <- vapply(seq_along(qty_cols), function(i) short_sample_name(qty_cols[[i]], i), character(1))
  logs <- zzx_log(logs, paste("[data] 上清丰度列:", paste(sprintf("%s <= %s", sample_names, qty_cols), collapse = ", ")))
  logs <- zzx_log(logs, "[note] 只把非胰酶末端当细胞剪切候选；K/R 胰酶切口不是细胞剪切位点。")
  logs <- zzx_log(logs, "[note] 无生物学重复时不伪造 p 值，按覆盖与丰度描述。")

  getc <- function(nm) if (nm %in% names(raw)) raw[[nm]] else rep("", nrow(raw))
  keep <- as.logical(mapply(
    is_target_row, getc("Genes"), getc("Protein.Ids"), getc("Protein.Names"), getc("Protein.Group"),
    SIMPLIFY = TRUE, USE.NAMES = FALSE
  ))
  alb <- raw[keep, , drop = FALSE]
  if (!nrow(alb)) {
    genes <- unique(getc("Genes"))
    stop("pr_matrix 里没有 ALDH1A1 / P00352。表中基因包括: ", paste(head(genes, 40), collapse = ", "))
  }
  logs <- zzx_log(logs, sprintf("[filter] ALDH1A1 前体行 %s / 总行 %s（含多电荷）", nrow(alb), nrow(raw)))

  protein <- TARGET_SEQUENCE
  pep_list <- list()
  unmapped <- 0L
  for (i in seq_len(nrow(alb))) {
    pep <- toupper(trimws(as.character(alb[["Stripped.Sequence"]][i])))
    if (!nzchar(pep) || is.na(pep)) {
      next
    }
    hits <- map_peptide(protein, pep)
    if (!nrow(hits)) {
      unmapped <- unmapped + 1L
      next
    }
    abund <- parse_num(unlist(alb[i, qty_cols, drop = TRUE], use.names = FALSE))
    if (all(is.na(abund))) {
      next
    }
    col_or_blank <- function(nm) {
      if (nm %in% names(alb)) as.character(alb[[nm]][i]) else ""
    }
    for (h in seq_len(nrow(hits))) {
      start <- hits$start[h]
      end <- hits$end[h]
      n_ok <- tryptic_n(protein, start)
      c_ok <- tryptic_c(protein, end)
      rec <- data.frame(
        Protein.Group = col_or_blank("Protein.Group"),
        Protein.Ids = col_or_blank("Protein.Ids"),
        Genes = col_or_blank("Genes"),
        Proteotypic = col_or_blank("Proteotypic"),
        Stripped.Sequence = pep,
        Modified.Sequence = col_or_blank("Modified.Sequence"),
        Precursor.Id = col_or_blank("Precursor.Id"),
        Precursor.Charge = col_or_blank("Precursor.Charge"),
        start = start,
        end = end,
        length = end - start + 1L,
        n_matches = nrow(hits),
        prev_aa = if (start > 1L) aa_at(protein, start - 1L) else "",
        next_aa = if (end < nchar(protein)) aa_at(protein, end + 1L) else "",
        tryptic_N = n_ok,
        tryptic_C = c_ok,
        cellular_neo_N = !n_ok,
        cellular_neo_C = !c_ok,
        terminus_class = terminus_class(n_ok, c_ok),
        cleavage_after_if_neoN = if (!n_ok) start - 1L else NA_integer_,
        cleavage_after_if_neoC = if (!c_ok) end else NA_integer_,
        abundance_mean = mean_or_na(abund),
        n_precursors = 1L,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      for (j in seq_along(sample_names)) {
        rec[[sample_names[[j]]]] <- abund[[j]]
      }
      pep_list[[length(pep_list) + 1L]] <- rec
    }
  }
  if (!length(pep_list)) {
    stop("ALDH1A1 肽段未能定位到 P00352 序列")
  }
  pep <- do.call(rbind, pep_list)
  logs <- zzx_log(logs, sprintf("[map] 定位到 P00352 的前体 %s 条；未映射 %s 条", nrow(pep), unmapped))

  key <- paste(pep$Stripped.Sequence, pep$start, pep$end, sep = "|")
  merged <- list()
  for (idx in split(seq_len(nrow(pep)), key)) {
    block <- pep[idx, , drop = FALSE]
    one <- block[1, , drop = FALSE]
    one$n_precursors <- nrow(block)
    one$Precursor.Charge <- "merged"
    for (nm in sample_names) {
      one[[nm]] <- if (all(is.na(block[[nm]]))) NA_real_ else sum(block[[nm]], na.rm = TRUE)
    }
    one$abundance_mean <- mean_or_na(unlist(one[1, sample_names], use.names = FALSE))
    one$terminus_class <- terminus_class(!one$cellular_neo_N, !one$cellular_neo_C)
    merged[[length(merged) + 1L]] <- one
  }
  pep <- do.call(rbind, merged)
  pep <- pep[order(pep$start, pep$end, pep$Stripped.Sequence), , drop = FALSE]
  rownames(pep) <- NULL

  n_neo <- sum(pep$cellular_neo_N | pep$cellular_neo_C)
  if (n_neo == 0L) {
    logs <- zzx_log(logs, "[warn] 没有任何 neo-N/neo-C 肽段。DIA 库很可能只搜了全胰酶肽，细胞剪切的直接末端证据会缺失；将主要依据 UniProt 加工位点 + 覆盖/丰度阶跃。")
  } else {
    logs <- zzx_log(logs, sprintf("[cleavage] 细胞剪切相关肽段（neo-N 或 neo-C）%s 条", n_neo))
  }

  n <- nchar(protein)
  profiles <- lapply(c(sample_names, "mean"), function(nm) {
    x <- rep(NA_real_, n)
    x
  })
  names(profiles) <- c(sample_names, "mean")
  counts <- integer(n)
  for (i in seq_len(nrow(pep))) {
    rng <- pep$start[i]:pep$end[i]
    counts[rng] <- counts[rng] + 1L
    for (nm in sample_names) {
      val <- pep[[nm]][i]
      if (is.na(val)) {
        next
      }
      cur <- profiles[[nm]][rng]
      profiles[[nm]][rng] <- ifelse(is.na(cur), val, cur + val)
    }
    val <- pep$abundance_mean[i]
    if (!is.na(val)) {
      cur <- profiles$mean[rng]
      profiles$mean[rng] <- ifelse(is.na(cur), val, cur + val)
    }
  }
  n_term_cov <- list(
    initiator = region_mean(profiles$mean, 1, 1),
    mature_n = region_mean(profiles$mean, 2, min(80L, n)),
    mature_rest = if (n > 80L) region_mean(profiles$mean, 81, n) else NA_real_
  )
  logs <- zzx_log(logs, sprintf(
    "[coverage] 区段平均丰度 M1=%s  链2-80=%s  其余81-%s=%s",
    fmt_num(n_term_cov$initiator), fmt_num(n_term_cov$mature_n), n, fmt_num(n_term_cov$mature_rest)
  ))

  neo_n <- collect_neo(pep, "N")
  neo_c <- collect_neo(pep, "C")
  candidates <- infer_uniprot(profiles$mean, neo_n, neo_c, protein)

  if (length(neo_n)) {
    for (key in names(neo_n)[order(-vapply(neo_n, `[[`, numeric(1), "abundance"))]) {
      pos <- as.integer(key)
      after <- pos - 1L
      if (after == 1L || after < 1L) {
        next
      }
      rec <- neo_n[[key]]
      c_rec <- neo_c[[as.character(after)]]
      candidates <- rbind(candidates, candidate_row(
        after, pos, "internal", "neo_N", "yes", protein,
        rec$n, rec$abundance,
        if (is.null(c_rec)) 0 else c_rec$n,
        if (is.null(c_rec)) NA_real_ else c_rec$abundance,
        region_mean(profiles$mean, max(1L, after - 14L), after),
        region_mean(profiles$mean, pos, min(n, pos + 14L)),
        NA_real_,
        sprintf("细胞 neo-N 肽段从残基 %s 开始（切割在 %s|%s），%s 条，Σ丰度 %s。肽段: %s", pos, after, pos, rec$n, format(rec$abundance, digits = 4), paste(head(rec$peptides, 8), collapse = ", ")),
        "非胰酶 N 端，优先视为细胞内/分泌路径蛋白酶切割"
      ))
    }
  }
  if (length(neo_c)) {
    for (key in names(neo_c)[order(-vapply(neo_c, `[[`, numeric(1), "abundance"))]) {
      after <- as.integer(key)
      if (after == 1L || after >= n) {
        next
      }
      if (any(candidates$after_residue == after & candidates$source == "neo_N")) {
        next
      }
      rec <- neo_c[[key]]
      nrec <- neo_n[[as.character(after + 1L)]]
      candidates <- rbind(candidates, candidate_row(
        after, after + 1L, "internal", "neo_C", "yes", protein,
        if (is.null(nrec)) 0 else nrec$n,
        if (is.null(nrec)) NA_real_ else nrec$abundance,
        rec$n, rec$abundance,
        region_mean(profiles$mean, max(1L, after - 14L), after),
        region_mean(profiles$mean, after + 1L, min(n, after + 15L)),
        NA_real_,
        sprintf("细胞 neo-C 肽段止于残基 %s（切割在 %s|%s），%s 条，Σ丰度 %s。肽段: %s", after, after, after + 1L, rec$n, format(rec$abundance, digits = 4), paste(head(rec$peptides, 8), collapse = ", ")),
        "非胰酶 C 端，优先视为细胞内/分泌路径蛋白酶切割"
      ))
    }
  }
  steps <- abundance_steps(profiles$mean, protein)
  if (is.data.frame(steps) && nrow(steps) > 0L) {
    candidates <- rbind(candidates, steps)
  }

  if (is.null(outdir)) {
    outdir <- file.path(dirname(report), "cleavage_ALDH1A1")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  pep_out <- pep
  write_table(file.path(outdir, "ALDH1A1_peptide_map.csv"), pep_out)
  write_table(file.path(outdir, "ALDH1A1_candidate_cleavage_sites.csv"), candidates)
  res_df <- data.frame(
    residue = seq_len(n),
    aa = strsplit(protein, "")[[1]],
    region = ifelse(seq_len(n) == 1L, "initiator_met", "chain"),
    n_peptides = counts,
    abundance_mean = profiles$mean,
    stringsAsFactors = FALSE
  )
  for (nm in sample_names) {
    res_df[[nm]] <- profiles[[nm]]
  }
  write_table(file.path(outdir, "ALDH1A1_residue_abundance.csv"), res_df)
  write_processing_svg(file.path(outdir, "ALDH1A1_processing_schematic.svg"), protein, candidates, n_term_cov)
  write_full_length_svg(file.path(outdir, "ALDH1A1_full_length_cleavage.svg"), protein, candidates, n_term_cov)
  write_coverage_svg(file.path(outdir, "ALDH1A1_peptide_coverage.svg"), protein, pep, sample_names, profiles)

  logs <- zzx_log(logs, "")
  logs <- zzx_log(logs, "======== ALDH1A1 细胞剪切推断 ========")
  met_row <- candidates[candidates$type == "initiator_met", ][1, ]
  logs <- zzx_log(logs, paste0("起始 Met 1|2 (", met_row$motif_P4_P4prime, "): ", met_row$evidence_zh))
  internals <- candidates[candidates$type == "internal" & candidates$source %in% c("neo_N", "neo_C"), , drop = FALSE]
  if (nrow(internals) > 0L) {
    logs <- zzx_log(logs, sprintf("内部细胞剪切候选 %s 个（neo 末端，已排除胰酶 K/R）:", nrow(internals)))
    for (i in seq_len(min(nrow(internals), 12L))) {
      logs <- zzx_log(logs, sprintf("  %s|%s  %s  %s", internals$after_residue[i], internals$before_residue[i], internals$motif_P4_P4prime[i], internals$protease_hint[i]))
    }
  } else {
    logs <- zzx_log(logs, "未发现内部 neo-N/neo-C。若搜库仅为胰酶特异性，内部剪切可能不可见。")
  }
  logs <- zzx_log(logs, paste("结果目录:", normalizePath(outdir, winslash = "/", mustWork = FALSE)))
  logs <- zzx_log(logs, "全长剪切示意图: ALDH1A1_full_length_cleavage.svg")
  writeLines(unlist(logs), file.path(outdir, "00_inference_log.txt"))
  invisible(outdir)
}

if (isTRUE(zzx_alb_auto_run)) {
  zzx_alb_run()
}
