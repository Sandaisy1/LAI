# =============================================================================
# ZZX DIA-NN：根据上清肽段丰度推断 ALB（P02768）的细胞剪切（不是胰酶）
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

ALB_UNIPROT <- "P02768"
ALB_GENE <- "ALB"
ALB_NAME <- "ALBU_HUMAN"
ALB_SEQUENCE <- paste0(
  "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF",
  "EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP",
  "ERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLF",
  "FAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAV",
  "ARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLK",
  "ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYAR",
  "RHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFE",
  "QLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVV",
  "LNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTL",
  "SEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLV",
  "AASQAALGL"
)
if (nchar(ALB_SEQUENCE) != 609L) {
  stop("ALB 序列长度应为 609，实际 ", nchar(ALB_SEQUENCE))
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

is_alb_row <- function(genes, ids, names, group) {
  g <- tokens(genes)
  id <- sub("-.*$", "", tokens(ids))
  nm <- toupper(as.character(names))
  gp <- toupper(as.character(group))
  ALB_GENE %in% g ||
    ALB_UNIPROT %in% id ||
    grepl(ALB_UNIPROT, gp, fixed = TRUE) ||
    grepl(ALB_NAME, nm, fixed = TRUE)
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
  if (after == 18L) {
    return("signal peptidase (ER)")
  }
  if (after == 24L) {
    return("furin-like proprotein convertase (Golgi, RXXR)")
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
  x_at <- function(res) 40 + 900 * (res - 1) / n
  y0 <- 80
  bar_h <- 36
  lines <- c(
    '<svg xmlns="http://www.w3.org/2000/svg" width="980" height="220">',
    "<style>text{font-family:Arial,Helvetica,sans-serif;font-size:12px}</style>",
    '<text x="40" y="28" font-size="16">ALB (P02768) cellular processing schematic</text>',
    paste0('<text x="40" y="48" fill="#555">signal 1-18 | propeptide 19-24 | mature 25-', n, "</text>"),
    sprintf('<rect x="%.1f" y="%s" width="%.1f" height="%s" fill="#f4c7c3" stroke="#333"/>', x_at(1), y0, x_at(19) - x_at(1), bar_h),
    sprintf('<rect x="%.1f" y="%s" width="%.1f" height="%s" fill="#ffe199" stroke="#333"/>', x_at(19), y0, x_at(25) - x_at(19), bar_h),
    sprintf('<rect x="%.1f" y="%s" width="%.1f" height="%s" fill="#b7d7b0" stroke="#333"/>', x_at(25), y0, x_at(n) + 900 / n - x_at(25), bar_h),
    sprintf('<text x="%.1f" y="%s" font-size="11">SIGNAL</text>', x_at(1) + 4, y0 + 22),
    sprintf('<text x="%.1f" y="%s" font-size="11">MATURE ALBUMIN (secreted)</text>', x_at(25) + 8, y0 + 22)
  )
  shown <- integer()
  for (i in seq_len(nrow(candidates))) {
    after <- as.integer(candidates$after_residue[i])
    if (after %in% shown || after < 1L || after >= n) {
      next
    }
    shown <- c(shown, after)
    color <- switch(as.character(candidates$type[i]),
      signal = "#c0392b",
      propeptide = "#d35400",
      "#6c3483"
    )
    x <- x_at(after + 1L)
    lines <- c(
      lines,
      sprintf('<line x1="%.1f" y1="%s" x2="%.1f" y2="%s" stroke="%s" stroke-width="2"/>', x, y0 - 8, x, y0 + bar_h + 18, color),
      sprintf('<text x="%.1f" y="%s" font-size="10" text-anchor="middle" fill="%s">%s|%s</text>', x, y0 + bar_h + 32, color, after, after + 1L)
    )
  }
  cov <- sprintf(
    "coverage 1-18=%s  19-24=%s  25-80=%s  81-%s=%s",
    fmt_num(n_term_cov$signal), fmt_num(n_term_cov$propeptide),
    fmt_num(n_term_cov$mature_n), n, fmt_num(n_term_cov$mature_rest)
  )
  lines <- c(lines, sprintf('<text x="40" y="204" fill="#333">%s</text>', svg_escape(cov)), "</svg>")
  write_svg(path, lines)
}

write_full_length_svg <- function(path, protein, candidates, n_term_cov) {
  n <- nchar(protein)
  internals <- candidates[candidates$type == "internal" & candidates$source %in% c("neo_N", "neo_C"), , drop = FALSE]
  if (nrow(internals) > 0L) {
    internals <- internals[order(internals$after_residue), , drop = FALSE]
    internals <- internals[!duplicated(internals$after_residue) & internals$after_residue > 24 & internals$after_residue < n, , drop = FALSE]
    if (nrow(internals) > 4L) {
      internals <- internals[seq_len(4L), , drop = FALSE]
    }
  }
  n_steps <- 3L + as.integer(nrow(internals) > 0L)
  left <- 130
  bar_w <- 920
  bar_h <- 42
  zoom_w <- 300
  rest_w <- bar_w - zoom_w
  step_gap <- 152
  top <- 96
  height <- 86 + n_steps * step_gap + 56
  if (nrow(internals) > 0L) {
    height <- height + 36
  }
  x_left <- function(res) {
    if (res <= 24) left + zoom_w * (res - 1) / 24 else left + zoom_w + rest_w * (res - 25) / (n - 24)
  }
  x_right <- function(res) {
    if (res <= 24) left + zoom_w * res / 24 else left + zoom_w + rest_w * (res - 24) / (n - 24)
  }
  fmt_cov <- function(v) {
    if (is.na(v)) "supernatant: not detected" else sprintf("supernatant: detected (%.3g)", v)
  }
  lines <- c(
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="%s" viewBox="0 0 1120 %s">', height, height),
    "<rect width='100%' height='100%' fill='#ffffff'/>",
    "<style>text{font-family:Microsoft YaHei,SimHei,Arial,Helvetica,sans-serif}</style>",
    '<text x="40" y="28" font-size="18" font-weight="700">Full-length ALB cleavage  全长蛋白如何剪切</text>',
    '<text x="40" y="50" font-size="12" fill="#555">P02768 prepro-albumin (1-609). Cellular processing only (not trypsin).</text>'
  )
  draw_domains <- function(y, keep_from, keep_to, ghost = FALSE) {
    domains <- list(
      c(1, 18, "#e74c3c", "#fadbd8", "SIGNAL", "1-18"),
      c(19, 24, "#e67e22", "#fdebd0", "PRO", "19-24"),
      c(25, n, "#1e8449", "#d5f5e3", "MATURE", paste0("25-", n))
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
      if (w >= 36) {
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
    lines, step_badge(1, y1, "Full-length precursor  全长前体 prepro-ALB", paste0("residues 1-", n, " (signal + propeptide + mature chain)")),
    draw_domains(y1, 1, n),
    draw_scissors(18, y1, "#c0392b", "cut 18|19", "SAYS|RGVF  signal peptidase (ER)"),
    down_arrow(y1, "ER: remove SIGNAL 1-18")
  )
  y2 <- y1 + step_gap
  lines <- c(
    lines,
    step_badge(2, y2, "After signal-peptide cleavage  切掉信号肽", paste0("released 1-18 (", fmt_cov(n_term_cov$signal), "); remaining proalbumin 19-", n)),
    draw_domains(y2, 1, 18, TRUE),
    draw_domains(y2, 19, n),
    draw_scissors(24, y2, "#d35400", "cut 24|25", "VFRR|DAHK  furin-like (Golgi, RXXR)"),
    down_arrow(y2, "Golgi: remove PROPEPTIDE 19-24")
  )
  y3 <- y2 + step_gap
  lines <- c(
    lines,
    step_badge(3, y3, "Mature albumin secreted  成熟链进入细胞上清", paste0("chain 25-", n, "  (", fmt_cov(n_term_cov$mature_n), " at 25-80)")),
    draw_domains(y3, 1, 24, TRUE),
    draw_domains(y3, 25, n),
    sprintf('<rect x="%.1f" y="%s" width="%.1f" height="%s" rx="8" fill="none" stroke="#1e8449" stroke-width="2.4"/>', x_left(25), y3 - 4, x_right(n) - x_left(25), bar_h + 8)
  )
  if (nrow(internals) > 0L) {
    y4 <- y3 + step_gap
    cuts <- paste(sprintf("%s|%s", internals$after_residue, internals$before_residue), collapse = ", ")
    lines <- c(lines, step_badge(4, y4, "Additional cellular cuts  成熟链上的内部剪切", paste0("neo-N / neo-C (not trypsin): ", cuts)), draw_domains(y4, 25, n))
    for (i in seq_len(nrow(internals))) {
      after <- as.integer(internals$after_residue[i])
      lines <- c(lines, draw_scissors(after, y4, "#6c3483", paste0("cut ", after, "|", after + 1L), paste(internals$motif_P4_P4prime[i], internals$protease_hint[i])))
    }
  }
  ly <- height - 28
  lines <- c(
    lines,
    sprintf('<rect x="%s" y="%s" width="16" height="12" fill="#e74c3c" opacity="0.35" stroke="#e74c3c"/>', left, ly - 12),
    sprintf('<text x="%s" y="%s" font-size="11">SIGNAL 1-18</text>', left + 22, ly),
    sprintf('<rect x="%s" y="%s" width="16" height="12" fill="#e67e22" opacity="0.35" stroke="#e67e22"/>', left + 160, ly - 12),
    sprintf('<text x="%s" y="%s" font-size="11">PROPEPTIDE 19-24</text>', left + 182, ly),
    sprintf('<rect x="%s" y="%s" width="16" height="12" fill="#1e8449" opacity="0.35" stroke="#1e8449"/>', left + 360, ly - 12),
    sprintf('<text x="%s" y="%s" font-size="11">MATURE 25-609 (secreted)</text>', left + 382, ly),
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
    '<text x="50" y="24" font-size="16">ALB residue abundance (cell supernatant) and peptide map</text>',
    '<text x="50" y="42" fill="#555">grey=trypsin; orange=cellular neo-N; blue=cellular neo-C; red=both</text>',
    sprintf('<rect x="%s" y="55" width="%s" height="%s" fill="#fafafa" stroke="#ccc"/>', x0, bar_w, track_h)
  )
  for (pair in list(c(1, 18, "#f4c7c3"), c(19, 24, "#ffe199"), c(25, n, "#eaf6e8"))) {
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
  lines <- c(lines, sprintf('<text x="%s" y="%s" fill="#555">residue 1 ... %s (prepro-albumin)</text>', x0, height - 12, n), "</svg>")
  write_svg(path, lines)
}

cov_label <- function(v) {
  if (!is.na(v) && v > 0) "detected" else "not_detected"
}

infer_uniprot <- function(profile, neo_n, neo_c, protein) {
  n <- nchar(protein)
  sig <- region_mean(profile, 1, 18)
  pro <- region_mean(profile, 19, 24)
  mat <- region_mean(profile, 25, min(80L, n))
  n19 <- neo_n[["19"]]
  c18 <- neo_c[["18"]]
  n25 <- neo_n[["25"]]
  c24 <- neo_c[["24"]]
  signal_ok <- identical(cov_label(sig), "not_detected") && identical(cov_label(mat), "detected")
  pro_ok <- identical(cov_label(pro), "not_detected") && identical(cov_label(mat), "detected")
  rbind(
    data.frame(
      after_residue = 18L, before_residue = 19L, type = "signal", source = "uniprot_annotated",
      supported = if (signal_ok || (!is.null(n19) && n19$n > 0) || (!is.null(c18) && c18$n > 0)) "yes" else "possible",
      motif_P4_P4prime = motif(protein, 18L),
      protease_hint = "signal peptidase (ER)",
      n_neo_N_peptides = if (is.null(n19)) 0 else n19$n,
      abundance_neo_N = if (is.null(n19)) NA_real_ else n19$abundance,
      n_neo_C_peptides = if (is.null(c18)) 0 else c18$n,
      abundance_neo_C = if (is.null(c18)) NA_real_ else c18$abundance,
      abundance_left_window = sig,
      abundance_right_window = if (is.na(pro)) mat else pro,
      log2_step = NA_real_,
      evidence_zh = paste0(
        "信号肽区 1-18 ", cov_label(sig), "，成熟 N 端 25-80 ", cov_label(mat), "。",
        if (signal_ok) "上清几乎只有成熟链覆盖，符合 ER 切掉信号肽后分泌。" else "信号肽区仍有肽段，需检查是否未加工前体或错误匹配。",
        " 仅空洞不能单独定论；neo-N@19 才是切割的直接肽段证据。"
      ),
      note_zh = "ER 信号肽酶切掉 1-18；成熟分泌蛋白不应再含这段",
      stringsAsFactors = FALSE
    ),
    data.frame(
      after_residue = 24L, before_residue = 25L, type = "propeptide", source = "uniprot_annotated",
      supported = if (pro_ok || (!is.null(n25) && n25$n > 0) || (!is.null(c24) && c24$n > 0)) "yes" else "possible",
      motif_P4_P4prime = motif(protein, 24L),
      protease_hint = "furin-like proprotein convertase (Golgi, RXXR)",
      n_neo_N_peptides = if (is.null(n25)) 0 else n25$n,
      abundance_neo_N = if (is.null(n25)) NA_real_ else n25$abundance,
      n_neo_C_peptides = if (is.null(c24)) 0 else c24$n,
      abundance_neo_C = if (is.null(c24)) NA_real_ else c24$abundance,
      abundance_left_window = pro,
      abundance_right_window = mat,
      log2_step = NA_real_,
      evidence_zh = paste0(
        "propeptide 19-24 ", cov_label(pro), "，成熟 N 端 25-80 ", cov_label(mat), "。",
        "R24|D25 在胰酶消化后是合法胰酶位点，因此成熟 N 端肽段通常不是 neo-N；",
        "要用 19-24 覆盖缺失 + 25 之后高覆盖 来支持高尔基体切除 propeptide。"
      ),
      note_zh = "高尔基体切除 RGVFRR（19-24），成熟 N 端为 D25",
      stringsAsFactors = FALSE
    )
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
    is_alb_row, getc("Genes"), getc("Protein.Ids"), getc("Protein.Names"), getc("Protein.Group"),
    SIMPLIFY = TRUE, USE.NAMES = FALSE
  ))
  alb <- raw[keep, , drop = FALSE]
  if (!nrow(alb)) {
    genes <- unique(getc("Genes"))
    stop("pr_matrix 里没有 ALB / P02768。表中基因包括: ", paste(head(genes, 40), collapse = ", "))
  }
  logs <- zzx_log(logs, sprintf("[filter] ALB 前体行 %s / 总行 %s（含多电荷）", nrow(alb), nrow(raw)))

  protein <- ALB_SEQUENCE
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
        n_matches_on_ALB = nrow(hits),
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
    stop("ALB 肽段未能定位到 P02768 序列")
  }
  pep <- do.call(rbind, pep_list)
  logs <- zzx_log(logs, sprintf("[map] 定位到 P02768 的前体 %s 条；未映射 %s 条", nrow(pep), unmapped))

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
    signal = region_mean(profiles$mean, 1, 18),
    propeptide = region_mean(profiles$mean, 19, 24),
    mature_n = region_mean(profiles$mean, 25, min(80L, n)),
    mature_rest = if (n > 80L) region_mean(profiles$mean, 81, n) else NA_real_
  )
  logs <- zzx_log(logs, sprintf(
    "[coverage] 区段平均丰度 信号肽1-18=%s  pro19-24=%s  成熟N25-80=%s  其余81-609=%s",
    fmt_num(n_term_cov$signal), fmt_num(n_term_cov$propeptide),
    fmt_num(n_term_cov$mature_n), fmt_num(n_term_cov$mature_rest)
  ))

  neo_n <- collect_neo(pep, "N")
  neo_c <- collect_neo(pep, "C")
  candidates <- infer_uniprot(profiles$mean, neo_n, neo_c, protein)

  if (length(neo_n)) {
    for (key in names(neo_n)[order(-vapply(neo_n, `[[`, numeric(1), "abundance"))]) {
      pos <- as.integer(key)
      after <- pos - 1L
      if (after %in% c(18L, 24L) || after < 1L) {
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
      if (after %in% c(18L, 24L) || after >= n) {
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
    outdir <- file.path(dirname(report), "cleavage_ALB")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  pep_out <- pep
  write_table(file.path(outdir, "ALB_peptide_map.csv"), pep_out)
  write_table(file.path(outdir, "ALB_candidate_cleavage_sites.csv"), candidates)
  res_df <- data.frame(
    residue = seq_len(n),
    aa = strsplit(protein, "")[[1]],
    region = ifelse(seq_len(n) <= 18, "signal", ifelse(seq_len(n) <= 24, "propeptide", "mature")),
    n_peptides = counts,
    abundance_mean = profiles$mean,
    stringsAsFactors = FALSE
  )
  for (nm in sample_names) {
    res_df[[nm]] <- profiles[[nm]]
  }
  write_table(file.path(outdir, "ALB_residue_abundance.csv"), res_df)
  write_processing_svg(file.path(outdir, "ALB_processing_schematic.svg"), protein, candidates, n_term_cov)
  write_full_length_svg(file.path(outdir, "ALB_full_length_cleavage.svg"), protein, candidates, n_term_cov)
  write_coverage_svg(file.path(outdir, "ALB_peptide_coverage.svg"), protein, pep, sample_names, profiles)

  logs <- zzx_log(logs, "")
  logs <- zzx_log(logs, "======== ALB 细胞剪切推断 ========")
  sig_row <- candidates[candidates$type == "signal", ][1, ]
  pro_row <- candidates[candidates$type == "propeptide", ][1, ]
  logs <- zzx_log(logs, paste0("信号肽 18|19 (", sig_row$motif_P4_P4prime, "): ", sig_row$evidence_zh))
  logs <- zzx_log(logs, paste0("Propeptide 24|25 (", pro_row$motif_P4_P4prime, "): ", pro_row$evidence_zh))
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
  logs <- zzx_log(logs, "全长剪切示意图: ALB_full_length_cleavage.svg")
  writeLines(unlist(logs), file.path(outdir, "00_inference_log.txt"))
  invisible(outdir)
}

if (isTRUE(zzx_alb_auto_run)) {
  zzx_alb_run()
}
