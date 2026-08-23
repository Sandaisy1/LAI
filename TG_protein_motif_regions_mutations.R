#!/usr/bin/env Rscript
# =============================================================================
# 蛋白 motif 区域导出 + 突变方案（不修改 TG_protein_motif_pipeline.R）
#   1) 每个蛋白的 motif 起止区域：Excel + 轨道图
#   2) 针对 motif 的氨基酸突变方案（默认焦点：RPL9 / RBP4 / ITGAV / ITGA2）
#      同一 motif 在该蛋白上的每一处命中都单独设计，不只第一处
#
# 只读已有结果：E:/R/Protein/results/ 下的 hits、PWM、FASTA
#
# 用法：
#   source("TG_protein_motif_pipeline.R")              # 先跑原分析（已跑过可跳过）
#   source("TG_protein_motif_regions_mutations.R")     # 只加区域图 / Excel / 突变方案
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

cran_required <- c("dplyr", "tidyr", "tibble", "stringr", "ggplot2", "writexl")
install_if_missing <- function(pkgs, required = TRUE) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) == 0) return(invisible(TRUE))
  tryCatch(
    install.packages(miss, repos = "https://cloud.r-project.org"),
    error = function(e) message("CRAN install failed: ", e$message)
  )
  still <- miss[!vapply(miss, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still) > 0 && required) stop("缺少必需 R 包: ", paste(still, collapse = ", "))
  invisible(TRUE)
}
install_if_missing(cran_required, required = TRUE)
for (p in cran_required) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# -----------------------------------------------------------------------------
# 路径（与原脚本同一数据根，但不 source 原脚本，避免重跑发现）
# -----------------------------------------------------------------------------
protein_data_candidates <- function() {
  env_dir <- Sys.getenv("TG_PROTEIN_DIR", unset = "")
  out <- c(env_dir, "E:/R/Protein", "E:\\R\\Protein")
  unique(out[nzchar(out)])
}

resolve_protein_dir <- function() {
  for (d in protein_data_candidates()) {
    if (dir.exists(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
  }
  for (d in c("E:/R/Protein", "E:\\R\\Protein")) {
    if (dir.exists(dirname(d)) || dir.exists("E:/") || dir.exists("E:\\")) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      if (dir.exists(d)) return(normalizePath(d, winslash = "/", mustWork = FALSE))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

protein_dir <- resolve_protein_dir()
result_dir  <- file.path(protein_dir, "results")
out_dir     <- file.path(result_dir, "regions_mutations")
fig_dir     <- file.path(out_dir, "figures")
per_dir     <- file.path(fig_dir, "per_protein")
focus_dir   <- file.path(fig_dir, "focus_mutations")
dir.create(per_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(focus_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, paste0("regions_mutations_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

AA20 <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
focus_default <- c("RPL9", "RBP4", "ITGAV", "ITGA2")
flank_n <- 10L
top_n_mut <- 3L

focus_genes_from_env <- function() {
  raw <- Sys.getenv("TG_MOTIF_FOCUS_GENES", unset = "")
  if (!nzchar(raw)) return(focus_default)
  unique(toupper(trimws(unlist(strsplit(raw, "[,;\\s]+")))))
}

# -----------------------------------------------------------------------------
# 读入原流程产物
# -----------------------------------------------------------------------------
read_simple_fasta <- function(path) {
  if (!file.exists(path)) return(character(0))
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

fasta_gene <- function(header) toupper(sub("\\|.*$", "", header))

load_motif_inputs <- function(result_dir) {
  hits_path <- file.path(result_dir, "motif_all_hits.csv")
  if (!file.exists(hits_path)) {
    stop("找不到 ", hits_path, "。请先运行 TG_protein_motif_pipeline.R")
  }
  hits <- utils::read.csv(hits_path, stringsAsFactors = FALSE)
  need <- c("gene", "start", "end", "site_sequence", "motif_id")
  miss <- setdiff(need, names(hits))
  if (length(miss) > 0) stop("hits 缺列: ", paste(miss, collapse = ", "))
  hits$gene <- toupper(as.character(hits$gene))
  hits$motif_id <- as.character(hits$motif_id)
  hits$start <- as.integer(hits$start)
  hits$end <- as.integer(hits$end)

  summary_path <- file.path(result_dir, "motif_significance_summary.csv")
  summary_tbl <- if (file.exists(summary_path)) {
    utils::read.csv(summary_path, stringsAsFactors = FALSE)
  } else {
    unique(hits[, intersect(c("motif_id", "width", "consensus", "empirical_p", "e_value"), names(hits)), drop = FALSE])
  }

  fa <- read_simple_fasta(file.path(result_dir, "00_sequences.fasta"))
  seqs <- character(0)
  if (length(fa) > 0) {
    seqs <- as.character(fa)
    names(seqs) <- vapply(names(fa), fasta_gene, character(1))
  }

  pwm <- list()
  for (mid in unique(hits$motif_id)) {
    p <- file.path(result_dir, mid, paste0(mid, "_pwm.csv"))
    if (file.exists(p)) {
      tab <- utils::read.csv(p, stringsAsFactors = FALSE)
      aa_cols <- intersect(AA20, names(tab))
      pwm[[mid]] <- as.matrix(tab[, aa_cols, drop = FALSE])
      rownames(pwm[[mid]]) <- as.character(tab$position)
    }
  }
  list(hits = hits, summary = summary_tbl, seqs = seqs, pwm = pwm)
}

protein_length <- function(gene, seqs, hits) {
  if (gene %in% names(seqs)) return(nchar(seqs[[gene]]))
  sub <- hits[hits$gene == gene, , drop = FALSE]
  if (nrow(sub) == 0) return(NA_integer_)
  max(sub$end, na.rm = TRUE)
}

flank_context <- function(seq, start, end, n) {
  if (is.na(seq) || !nzchar(seq)) return(list(left = NA_character_, right = NA_character_, window = NA_character_))
  L <- nchar(seq)
  left_i <- max(1L, start - n)
  right_i <- min(L, end + n)
  list(
    left = if (start > 1) substring(seq, left_i, start - 1) else "",
    right = if (end < L) substring(seq, end + 1, right_i) else "",
    window = substring(seq, left_i, right_i)
  )
}

# -----------------------------------------------------------------------------
# 区域表
# -----------------------------------------------------------------------------
build_region_table <- function(hits, seqs) {
  rows <- vector("list", nrow(hits))
  for (i in seq_len(nrow(hits))) {
    g <- hits$gene[i]
    seq <- if (g %in% names(seqs)) seqs[[g]] else NA_character_
    ctx <- flank_context(seq, hits$start[i], hits$end[i], flank_n)
    plen <- protein_length(g, seqs, hits)
    map <- if (is.finite(plen) && plen > 0) {
      pre <- max(0, hits$start[i] - 1)
      post <- max(0, plen - hits$end[i])
      paste0(
        if (pre > 0) paste0("[1-", pre, "]") else "",
        "{", hits$motif_id[i], ":", hits$start[i], "-", hits$end[i], ":", hits$site_sequence[i], "}",
        if (post > 0) paste0("[", hits$end[i] + 1, "-", plen, "]") else ""
      )
    } else {
      paste0(hits$motif_id[i], ":", hits$start[i], "-", hits$end[i])
    }
    rows[[i]] <- data.frame(
      gene = g,
      uniprot = if ("uniprot" %in% names(hits)) hits$uniprot[i] else NA_character_,
      protein_length = plen,
      motif_id = hits$motif_id[i],
      motif_start = hits$start[i],
      motif_end = hits$end[i],
      motif_width = hits$end[i] - hits$start[i] + 1,
      site_sequence = hits$site_sequence[i],
      left_flank = ctx$left,
      right_flank = ctx$right,
      window_seq = ctx$window,
      region_map = map,
      start_frac = if (is.finite(plen) && plen > 0) hits$start[i] / plen else NA_real_,
      end_frac = if (is.finite(plen) && plen > 0) hits$end[i] / plen else NA_real_,
      llr_score = if ("llr_score" %in% names(hits)) hits$llr_score[i] else NA_real_,
      motif_empirical_p = if ("motif_empirical_p" %in% names(hits)) hits$motif_empirical_p[i] else NA_real_,
      motif_e_value = if ("motif_e_value" %in% names(hits)) hits$motif_e_value[i] else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

protein_motif_presence <- function(regions, all_motifs) {
  genes <- unique(regions$gene)
  out <- expand.grid(gene = genes, motif_id = all_motifs, stringsAsFactors = FALSE)
  merged <- merge(out, regions[, c("gene", "motif_id", "motif_start", "motif_end", "site_sequence", "llr_score")],
                  by = c("gene", "motif_id"), all.x = TRUE)
  merged$has_motif <- !is.na(merged$motif_start)
  merged[order(merged$gene, merged$motif_id), , drop = FALSE]
}

# -----------------------------------------------------------------------------
# 突变方案：打高信息量位点，破坏 motif
# -----------------------------------------------------------------------------
col_ic <- function(p) {
  p <- as.numeric(p)
  p[!is.finite(p) | p <= 0] <- 1e-12
  p <- p / sum(p)
  log2(length(p)) + sum(p * log2(p))
}

suggest_aa <- function(wt) {
  wt <- toupper(wt)
  if (is.na(wt) || !nzchar(wt)) return("A")
  if (wt == "C") return("S")
  if (wt == "G") return("A")
  if (wt == "P") return("A")
  if (wt %in% c("K", "R")) return("E")
  if (wt %in% c("D", "E")) return("K")
  if (wt %in% c("S", "T")) return("A")
  if (wt %in% c("Y", "W")) return("F")
  if (wt == "A") return("G")
  "A"
}

mut_rationale <- function(wt, mut) {
  if (wt == "C" && mut == "S") return("C->S：去掉巯基/二硫键，体积接近")
  if (wt == "G" && mut == "A") return("G->A：降低柔性，破坏转角")
  if (wt == "P" && mut == "A") return("P->A：去掉脯氨酸折角")
  if (wt %in% c("K", "R") && mut == "E") return("碱性->酸性：电荷翻转")
  if (wt %in% c("D", "E") && mut == "K") return("酸性->碱性：电荷翻转")
  if (wt %in% c("S", "T") && mut == "A") return("去极性侧链（Ala scan）")
  if (wt %in% c("Y", "W") && mut == "F") return("去掉芳香羟基/吲哚，保留疏水环")
  if (wt == "A" && mut == "G") return("A->G：改变主链柔性")
  "Ala scan：用丙氨酸替换，尽量只去掉侧链"
}

pwm_position_stats <- function(pwm_mat) {
  if (is.null(pwm_mat) || !nrow(pwm_mat)) return(NULL)
  npos <- nrow(pwm_mat)
  data.frame(
    motif_pos = seq_len(npos),
    ic_bits = vapply(seq_len(npos), function(j) col_ic(pwm_mat[j, ]), numeric(1)),
    max_aa = apply(pwm_mat, 1, function(v) names(v)[which.max(v)]),
    max_p = apply(pwm_mat, 1, max),
    stringsAsFactors = FALSE
  )
}

build_mutation_table <- function(regions, pwm_list, seqs, hits) {
  if (nrow(regions) == 0) return(data.frame())
  out <- list()
  k <- 0
  for (i in seq_len(nrow(regions))) {
    mid <- regions$motif_id[i]
    stats <- pwm_position_stats(pwm_list[[mid]])
    site <- strsplit(as.character(regions$site_sequence[i]), "")[[1]]
    width <- length(site)
    if (!is.null(stats)) {
      stats <- stats[stats$motif_pos <= width, , drop = FALSE]
    } else {
      stats <- data.frame(
        motif_pos = seq_len(width),
        ic_bits = 0,
        max_aa = site,
        max_p = NA_real_,
        stringsAsFactors = FALSE
      )
    }
    stats <- stats[order(-stats$ic_bits, stats$motif_pos), , drop = FALSE]
    stats$ic_rank <- seq_len(nrow(stats))
    for (j in seq_len(nrow(stats))) {
      mp <- stats$motif_pos[j]
      wt <- toupper(site[mp])
      if (is.na(wt) || !nzchar(wt) || wt == "X") next
      mut <- suggest_aa(wt)
      aa_pos <- regions$motif_start[i] + mp - 1
      k <- k + 1
      out[[k]] <- data.frame(
        gene = regions$gene[i],
        uniprot = regions$uniprot[i],
        motif_id = mid,
        motif_start = regions$motif_start[i],
        motif_end = regions$motif_end[i],
        site_sequence = regions$site_sequence[i],
        motif_pos = mp,
        protein_aa_pos = aa_pos,
        wt_aa = wt,
        suggested_aa = mut,
        hgvs_p = paste0("p.", wt, aa_pos, mut),
        short_mut = paste0(wt, aa_pos, mut),
        ic_bits = stats$ic_bits[j],
        ic_rank = stats$ic_rank[j],
        pwm_max_aa = stats$max_aa[j],
        pwm_max_p = stats$max_p[j],
        priority = if (stats$ic_rank[j] <= top_n_mut) "P1_core" else "P2_support",
        rationale = mut_rationale(wt, mut),
        stringsAsFactors = FALSE
      )
    }
  }
  if (k == 0) return(data.frame())
  dplyr::bind_rows(out)
}

core_mutations_for_site <- function(mut_tbl, gene, motif_id, start, end) {
  core <- mut_tbl[
    mut_tbl$gene == gene &
      mut_tbl$motif_id == motif_id &
      mut_tbl$motif_start == start &
      mut_tbl$motif_end == end &
      mut_tbl$priority == "P1_core",
    , drop = FALSE
  ]
  if (nrow(core) == 0) return(core)
  core[order(core$ic_rank), , drop = FALSE]
}

build_combo_table <- function(mut_tbl, regions, all_motifs) {
  genes <- unique(c(regions$gene, mut_tbl$gene))
  rows <- list()
  for (g in unique(genes)) {
    for (mid in all_motifs) {
      hit <- regions[regions$gene == g & regions$motif_id == mid, , drop = FALSE]
      if (nrow(hit) == 0) {
        rows[[length(rows) + 1]] <- data.frame(
          gene = g,
          motif_id = mid,
          occurrence = NA_integer_,
          n_occurrences = 0L,
          motif_start = NA_integer_,
          motif_end = NA_integer_,
          site_sequence = NA_character_,
          occurrence_id = NA_character_,
          has_motif = FALSE,
          combo_short = NA_character_,
          combo_hgvs = NA_character_,
          n_mutations = 0L,
          target_residues = NA_character_,
          plan = paste0(g, " 未命中 ", mid, "，无需针对该 motif 突变"),
          stringsAsFactors = FALSE
        )
        next
      }
      hit <- hit[order(hit$motif_start, hit$motif_end), , drop = FALSE]
      n_occ <- nrow(hit)
      for (h in seq_len(n_occ)) {
        st <- hit$motif_start[h]
        en <- hit$motif_end[h]
        core <- core_mutations_for_site(mut_tbl, g, mid, st, en)
        combo <- paste(core$short_mut, collapse = "/")
        hgvs <- paste(core$hgvs_p, collapse = "; ")
        occ_id <- paste0(mid, "_", st, "-", en)
        rows[[length(rows) + 1]] <- data.frame(
          gene = g,
          motif_id = mid,
          occurrence = h,
          n_occurrences = n_occ,
          motif_start = st,
          motif_end = en,
          site_sequence = hit$site_sequence[h],
          occurrence_id = occ_id,
          has_motif = TRUE,
          combo_short = combo,
          combo_hgvs = hgvs,
          n_mutations = nrow(core),
          target_residues = paste(paste0(core$wt_aa, core$protein_aa_pos), collapse = ", "),
          plan = sprintf(
            "%s / %s 第%d/%d处（%d-%d，%s）：同时突变核心高IC位点 %s。建议先做单点（P1），再做三联组合。每一处命中都单独做，不要只改第一处。",
            g, mid, h, n_occ, st, en, hit$site_sequence[h], combo
          ),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  dplyr::bind_rows(rows)
}

# -----------------------------------------------------------------------------
# 作图
# -----------------------------------------------------------------------------
motif_colors <- function(ids) {
  ids <- unique(as.character(ids))
  pal <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")
  setNames(rep(pal, length.out = length(ids)), ids)
}

safe_ggsave <- function(plot, stub, width, height) {
  tryCatch(
    ggplot2::ggsave(paste0(stub, ".pdf"), plot, width = width, height = height, limitsize = FALSE),
    error = function(e) log_msg("pdf failed ", stub, ": ", e$message)
  )
  tryCatch(
    ggplot2::ggsave(paste0(stub, ".png"), plot, width = width, height = height, dpi = 160, limitsize = FALSE),
    error = function(e) log_msg("png failed ", stub, ": ", e$message)
  )
}

plot_protein_track <- function(gene, regions, seqs, hits, title_extra = "") {
  plen <- protein_length(gene, seqs, hits)
  sub <- regions[regions$gene == gene, , drop = FALSE]
  if (!is.finite(plen) || plen < 1) plen <- if (nrow(sub)) max(sub$motif_end) else 1
  backbone <- data.frame(xmin = 1, xmax = plen, ymin = 0.35, ymax = 0.65)
  cols <- motif_colors(unique(regions$motif_id))
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = backbone,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "#D0D0D0", color = NA
    )
  if (nrow(sub) > 0) {
    p <- p +
      ggplot2::geom_rect(
        data = sub,
        ggplot2::aes(xmin = motif_start, xmax = motif_end, ymin = 0.15, ymax = 0.85, fill = motif_id),
        color = "white", linewidth = 0.2, alpha = 0.9
      ) +
      ggplot2::geom_text(
        data = sub,
        ggplot2::aes(x = (motif_start + motif_end) / 2, y = 1.05, label = paste0(motif_id, "\n", motif_start, "-", motif_end)),
        size = 2.4, vjust = 0, lineheight = 0.9
      )
  }
  p +
    ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(1, max(plen, 2)), expand = c(0.02, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1.45), expand = c(0, 0)) +
    ggplot2::labs(
      title = paste0(gene, " motif regions", title_extra),
      subtitle = paste0("protein length ", plen, " aa"),
      x = "amino acid position", y = NULL, fill = "motif"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

plot_all_tracks_frac <- function(regions, seqs, hits) {
  genes <- sort(unique(regions$gene))
  lens <- vapply(genes, function(g) protein_length(g, seqs, hits), numeric(1))
  df <- regions
  df$gene <- factor(df$gene, levels = rev(genes))
  cols <- motif_colors(unique(df$motif_id))
  ggplot2::ggplot(df, ggplot2::aes(xmin = start_frac, xmax = end_frac, ymin = as.numeric(gene) - 0.35, ymax = as.numeric(gene) + 0.35, fill = motif_id)) +
    ggplot2::geom_rect(alpha = 0.9) +
    ggplot2::scale_y_continuous(breaks = seq_along(levels(df$gene)), labels = levels(df$gene)) +
    ggplot2::scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(round(100 * x), "%"), expand = c(0.01, 0)) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(title = "All proteins: motif regions (relative position)", x = "position along protein", y = NULL, fill = "motif") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

plot_presence_heatmap <- function(presence) {
  presence$gene <- factor(presence$gene, levels = rev(sort(unique(presence$gene))))
  ggplot2::ggplot(presence, ggplot2::aes(x = motif_id, y = gene, fill = has_motif)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#2C7BB6", "FALSE" = "#F0F0F0"), labels = c("TRUE" = "yes", "FALSE" = "no")) +
    ggplot2::labs(title = "Protein x motif presence", x = NULL, y = NULL, fill = "hit") +
    ggplot2::theme_minimal(base_size = 11)
}

plot_mutation_focus <- function(gene, regions, mut_tbl, seqs, hits) {
  sub_r <- regions[regions$gene == gene, , drop = FALSE]
  sub_m <- mut_tbl[mut_tbl$gene == gene & mut_tbl$priority == "P1_core", , drop = FALSE]
  track <- plot_protein_track(gene, regions, seqs, hits, " + P1 mutations")
  if (nrow(sub_m) > 0) {
    track <- track +
      ggplot2::geom_point(
        data = sub_m,
        ggplot2::aes(x = protein_aa_pos, y = 0.5),
        inherit.aes = FALSE, size = 2.2, color = "black"
      ) +
      ggplot2::geom_text(
        data = sub_m,
        ggplot2::aes(x = protein_aa_pos, y = 0.02, label = short_mut),
        inherit.aes = FALSE, size = 2.1, angle = 90, hjust = 0
      )
  }

  strips <- list()
  if (nrow(sub_r) > 0) {
    sub_r <- sub_r[order(sub_r$motif_id, sub_r$motif_start, sub_r$motif_end), , drop = FALSE]
  }
  for (i in seq_len(nrow(sub_r))) {
    hit <- sub_r[i, ]
    mid <- hit$motif_id
    occ_id <- paste0(mid, "_", hit$motif_start, "-", hit$motif_end)
    aa <- strsplit(as.character(hit$site_sequence), "")[[1]]
    mm <- core_mutations_for_site(sub_m, gene, mid, hit$motif_start, hit$motif_end)
    df <- data.frame(
      motif_pos = seq_along(aa),
      aa = aa,
      protein_pos = hit$motif_start + seq_along(aa) - 1,
      stringsAsFactors = FALSE
    )
    df$mut <- mm$suggested_aa[match(df$motif_pos, mm$motif_pos)]
    df$label <- ifelse(!is.na(df$mut), paste0(df$aa, "\n", df$protein_pos, "\n->\n", df$mut), df$aa)
    df$is_core <- !is.na(df$mut)
    same <- which(sub_r$motif_id == mid)
    n_same <- length(same)
    occ_n <- match(i, same)
    strips[[occ_id]] <- ggplot2::ggplot(df, ggplot2::aes(x = motif_pos, y = 1, fill = is_core, label = label)) +
      ggplot2::geom_tile(width = 0.92, height = 0.9, color = "white") +
      ggplot2::geom_text(size = 2.1, lineheight = 0.85) +
      ggplot2::scale_fill_manual(values = c("TRUE" = "#F4A582", "FALSE" = "#FDDBC7"), guide = "none") +
      ggplot2::scale_x_continuous(breaks = df$motif_pos) +
      ggplot2::labs(
        title = paste0(gene, "  ", mid, "  ", hit$motif_start, "-", hit$motif_end,
                       "  (", occ_n, "/", n_same, ")"),
        x = "motif position", y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
  }
  list(track = track, strips = strips)
}

# -----------------------------------------------------------------------------
# 主流程
# -----------------------------------------------------------------------------
run_regions_mutations <- function() {
  dir.create(per_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(focus_dir, recursive = TRUE, showWarnings = FALSE)
  log_msg("Protein dir: ", protein_dir)
  log_msg("Reading results: ", result_dir)
  inp <- load_motif_inputs(result_dir)
  hits <- inp$hits
  seqs <- inp$seqs
  pwm <- inp$pwm
  focus <- intersect(focus_genes_from_env(), unique(hits$gene))
  if (length(focus) == 0) {
    log_msg("焦点基因均未出现在 hits 中，改用 hits 里前 4 个基因")
    focus <- head(unique(hits$gene), 4)
  }
  log_msg("Focus genes: ", paste(focus, collapse = ", "))

  regions <- build_region_table(hits, seqs)
  all_motifs <- sort(unique(hits$motif_id))
  presence <- protein_motif_presence(regions, all_motifs)
  mut_tbl <- build_mutation_table(regions, pwm, seqs, hits)
  combo_all <- build_combo_table(mut_tbl, regions, all_motifs)
  mut_focus <- mut_tbl[mut_tbl$gene %in% focus, , drop = FALSE]
  combo_focus <- combo_all[combo_all$gene %in% focus, , drop = FALSE]
  regions_focus <- regions[regions$gene %in% focus, , drop = FALSE]

  excel_regions <- file.path(out_dir, "protein_motif_regions.xlsx")
  writexl::write_xlsx(
    list(
      regions_all = regions,
      presence = presence,
      focus_regions = regions_focus
    ),
    excel_regions
  )
  excel_mut <- file.path(out_dir, "protein_motif_mutation_plans.xlsx")
  writexl::write_xlsx(
    list(
      mutation_focus_residues = mut_focus,
      mutation_focus_combo = combo_focus,
      mutation_all_residues = mut_tbl,
      mutation_all_combo = combo_all
    ),
    excel_mut
  )
  log_msg("Excel: ", excel_regions)
  log_msg("Excel: ", excel_mut)

  ov <- plot_all_tracks_frac(regions, seqs, hits)
  safe_ggsave(ov, file.path(fig_dir, "all_proteins_motif_tracks"), width = 9, height = max(6, 0.28 * length(unique(regions$gene))))
  hm <- plot_presence_heatmap(presence)
  safe_ggsave(hm, file.path(fig_dir, "all_proteins_motif_presence"), width = 6, height = max(5, 0.22 * length(unique(presence$gene))))

  for (g in unique(regions$gene)) {
    p <- plot_protein_track(g, regions, seqs, hits)
    safe_ggsave(p, file.path(per_dir, paste0(g, "_motif_regions")), width = 9, height = 2.4)
  }

  for (g in focus) {
    parts <- plot_mutation_focus(g, regions, mut_tbl, seqs, hits)
    safe_ggsave(parts$track, file.path(focus_dir, paste0(g, "_track_mutations")), width = 9, height = 2.8)
    for (mid in names(parts$strips)) {
      safe_ggsave(parts$strips[[mid]], file.path(focus_dir, paste0(g, "_", mid, "_mutation_strip")), width = 9, height = 2.3)
    }
  }

  if (nrow(combo_focus) > 0) {
    cat("\n===== 焦点蛋白突变方案（同一 motif 的每一处命中都单独做）=====\n")
    show_cols <- intersect(
      c("gene", "motif_id", "occurrence", "n_occurrences", "motif_start", "motif_end",
        "has_motif", "combo_short", "plan"),
      names(combo_focus)
    )
    print(combo_focus[, show_cols, drop = FALSE])
  }
  log_msg("Done -> ", out_dir)
  invisible(list(regions = regions, mutations = mut_tbl, combo = combo_all, focus = focus))
}

# -----------------------------------------------------------------------------
# 自检（不访问原 pipeline、不访问 UniProt）
# -----------------------------------------------------------------------------
run_regions_selftest <- function() {
  tmp <- file.path(tempdir(), paste0("motif_regions_selftest_", as.integer(Sys.time())))
  dir.create(file.path(tmp, "results", "motif1"), recursive = TRUE, showWarnings = FALSE)
  seqs <- c(
    RPL9 = paste0(paste(rep("A", 20), collapse = ""), "CADCQEGGGC", paste(rep("A", 20), collapse = ""), "CADCQEGGGC", paste(rep("A", 10), collapse = "")),
    RBP4 = paste0(paste(rep("S", 15), collapse = ""), "CADCQEGGGC", paste(rep("S", 30), collapse = "")),
    ITGAV = paste0(paste(rep("G", 40), collapse = ""), "CADCQEGGGC", paste(rep("G", 10), collapse = "")),
    ITGA2 = paste0(paste(rep("V", 8), collapse = ""), "CADCQEGGGC", paste(rep("V", 12), collapse = ""))
  )
  fa <- file(file.path(tmp, "results", "00_sequences.fasta"), "w")
  for (nm in names(seqs)) {
    writeLines(paste0(">", nm, "|TEST"), fa)
    writeLines(seqs[[nm]], fa)
  }
  close(fa)
  hits <- data.frame(
    gene = c("RPL9", "RPL9", "RBP4", "ITGAV", "ITGA2"),
    uniprot = c("P1", "P1", "P2", "P3", "P4"),
    start = c(21, 51, 16, 41, 9),
    end = c(30, 60, 25, 50, 18),
    site_sequence = "CADCQEGGGC",
    llr_score = 10,
    motif_id = "motif1",
    motif_empirical_p = 0.005,
    motif_e_value = 0.1,
    stringsAsFactors = FALSE
  )
  utils::write.csv(hits, file.path(tmp, "results", "motif_all_hits.csv"), row.names = FALSE)
  pwm <- matrix(0.02, nrow = 10, ncol = 20, dimnames = list(NULL, AA20))
  planted <- strsplit("CADCQEGGGC", "")[[1]]
  for (j in seq_along(planted)) pwm[j, planted[j]] <- 0.64
  pwm <- sweep(pwm, 1, rowSums(pwm), "/")
  pwm_tab <- data.frame(position = seq_len(10), pwm, check.names = FALSE)
  utils::write.csv(pwm_tab, file.path(tmp, "results", "motif1", "motif1_pwm.csv"), row.names = FALSE)

  old_protein <- protein_dir
  old_result <- result_dir
  old_out <- out_dir
  old_fig <- fig_dir
  old_per <- per_dir
  old_focus <- focus_dir
  old_log <- log_file
  assign("protein_dir", tmp, envir = .GlobalEnv)
  assign("result_dir", file.path(tmp, "results"), envir = .GlobalEnv)
  assign("out_dir", file.path(tmp, "results", "regions_mutations"), envir = .GlobalEnv)
  assign("fig_dir", file.path(tmp, "results", "regions_mutations", "figures"), envir = .GlobalEnv)
  assign("per_dir", file.path(tmp, "results", "regions_mutations", "figures", "per_protein"), envir = .GlobalEnv)
  assign("focus_dir", file.path(tmp, "results", "regions_mutations", "figures", "focus_mutations"), envir = .GlobalEnv)
  dir.create(per_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(focus_dir, recursive = TRUE, showWarnings = FALSE)
  assign("log_file", file.path(out_dir, "selftest.log"), envir = .GlobalEnv)
  on.exit({
    assign("protein_dir", old_protein, envir = .GlobalEnv)
    assign("result_dir", old_result, envir = .GlobalEnv)
    assign("out_dir", old_out, envir = .GlobalEnv)
    assign("fig_dir", old_fig, envir = .GlobalEnv)
    assign("per_dir", old_per, envir = .GlobalEnv)
    assign("focus_dir", old_focus, envir = .GlobalEnv)
    assign("log_file", old_log, envir = .GlobalEnv)
  }, add = TRUE)

  res <- run_regions_mutations()
  if (!file.exists(file.path(out_dir, "protein_motif_regions.xlsx"))) stop("selftest: missing regions xlsx")
  if (!file.exists(file.path(out_dir, "protein_motif_mutation_plans.xlsx"))) stop("selftest: missing mutation xlsx")
  rpl9 <- res$mutations[res$mutations$gene == "RPL9" & res$mutations$priority == "P1_core", ]
  if (nrow(rpl9) < 1) stop("selftest: RPL9 has no P1 mutations")
  if (!any(rpl9$wt_aa == "C" & rpl9$suggested_aa == "S")) stop("selftest: expected C->S")
  rpl9_sites <- unique(paste(rpl9$motif_start, rpl9$motif_end, sep = "-"))
  if (length(rpl9_sites) < 2) {
    stop("selftest: RPL9 motif1 must have mutations at every occurrence, got ", paste(rpl9_sites, collapse = ","))
  }
  rpl9_combo <- res$combo[res$combo$gene == "RPL9" & res$combo$motif_id == "motif1" & isTRUE(res$combo$has_motif), ]
  if (nrow(rpl9_combo) < 2) {
    stop("selftest: combo table must have one row per motif1 site on RPL9")
  }
  strip_files <- list.files(focus_dir, pattern = "^RPL9_motif1_.*_mutation_strip")
  if (length(strip_files) < 2) {
    stop("selftest: expected a mutation strip for each RPL9 motif1 site")
  }
  log_msg("REGIONS SELFTEST passed")
  invisible(TRUE)
}

if (identical(Sys.getenv("TG_MOTIF_REGIONS_SELFTEST"), "1")) {
  run_regions_selftest()
} else {
  run_regions_mutations()
}
