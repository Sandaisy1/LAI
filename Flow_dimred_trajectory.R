#!/usr/bin/env Rscript
# =============================================================================
# 每个 panel、每一大类免疫细胞的轨迹（Monocle 风格骨架图）
#
# 不是把 P1/P2/P3 拼成一张矩阵。每个 panel 里按大类单独降维、单独画树：
#   P1: CD4、CD8（NK/NKT 只有一个亚群则跳过；B/髓系是 dump）
#   P2: B（Naive / Memory / Switched / Plasma / Activated 一条分化轴）
#   P3: 髓系（单核→巨噬/DC，中性粒等为旁支）
#
# 根节点按惯例：CD4/CD8 naive、Naive B、Ly6C hi 单核。轨迹表示标志物
# 状态连续体，不是发育起源的证明。
#
# 用法：
#   setwd("E:/R/flow J-LJY WJZ ZZX")
#   source("Flow_dimred_pipeline.R")       # 降维后自动调用
#   source("Flow_dimred_trajectory.R")     # 只重出轨迹
#
# 结果：results_flow/P1/trajectory/ 等
# =============================================================================

load_flow_pipeline_functions <- function() {
  if (exists("save_gg", mode = "function") && exists("celltype_label", mode = "function")) {
    return(invisible(TRUE))
  }
  pipe <- "Flow_dimred_pipeline.R"
  cands <- c(pipe, file.path(getwd(), pipe))
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    cands <- c(file.path(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))), pipe), cands)
  }
  ofile <- NULL
  if (sys.nframe() > 0) {
    for (i in sys.nframe():1) {
      ofile <- sys.frame(i)$ofile
      if (!is.null(ofile) && nzchar(ofile)) break
    }
  }
  if (!is.null(ofile) && nzchar(ofile)) {
    cands <- c(file.path(dirname(normalizePath(ofile, mustWork = FALSE)), pipe), cands)
  }
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit) || !nzchar(hit)) stop("找不到 Flow_dimred_pipeline.R")
  old <- Sys.getenv("FLOW_FUNCTIONS_ONLY", unset = NA)
  Sys.setenv(FLOW_FUNCTIONS_ONLY = "1")
  source(hit, local = FALSE)
  if (is.na(old)) Sys.unsetenv("FLOW_FUNCTIONS_ONLY") else Sys.setenv(FLOW_FUNCTIONS_ONLY = old)
  invisible(TRUE)
}

infer_major_lineage <- function(panel_id, subset) {
  s <- as.character(subset)
  n <- length(s)
  if (identical(panel_id, "P1")) {
    out <- s
    out[s %in% c("Treg", "CD4_activated") | grepl("^CD4_", s)] <- "CD4"
    out[grepl("^CD8_", s)] <- "CD8"
    out[s %in% c("B", "Myeloid")] <- "dump"
    return(out)
  }
  if (identical(panel_id, "P2")) return(rep("B", n))
  if (identical(panel_id, "P3")) {
    out <- rep("Myeloid", n)
    out[s %in% c("T", "B", "NK", "NKT", "CD4", "CD8")] <- "dump"
    return(out)
  }
  s
}

majors_for_trajectory <- function(panel_id, cells) {
  if (identical(panel_id, "P2")) return(rep("B", nrow(cells)))
  if ("cluster_lineage" %in% names(cells)) {
    maj <- as.character(cells$cluster_lineage)
    if (identical(panel_id, "P1")) maj[maj %in% c("B", "Myeloid")] <- "dump"
    if (identical(panel_id, "P3")) maj[maj %in% c("T", "B", "NK", "NKT", "CD4", "CD8")] <- "dump"
    if (any(nzchar(maj) & !is.na(maj))) return(maj)
  }
  infer_major_lineage(panel_id, cells$lineage)
}

default_trajectory_root <- function(panel_id, major, lineages) {
  lin <- unique(as.character(lineages))
  prefer <- NULL
  if (identical(panel_id, "P1") && identical(major, "CD4")) prefer <- "CD4_naive"
  if (identical(panel_id, "P1") && identical(major, "CD8")) prefer <- "CD8_naive"
  if (identical(panel_id, "P2") && identical(major, "B")) prefer <- "Naive_B"
  if (identical(panel_id, "P3") && identical(major, "Myeloid")) prefer <- "Mono_Ly6Chi"
  if (!is.null(prefer) && prefer %in% lin) return(prefer)
  names(sort(table(lineages), decreasing = TRUE))[1]
}

marker_columns <- function(df) {
  meta <- c("sample", "group", "cluster", "cluster_lineage", "lineage",
            "UMAP1", "UMAP2", "tSNE1", "tSNE2", "true_lineage",
            "Component1", "Component2", "pseudotime", "major")
  num <- names(df)[vapply(df, is.numeric, logical(1))]
  setdiff(num, meta)
}

# 类内 PCA 当 Component 1/2（接近图里的 DDRTree 坐标轴）；没有标志物则退回 UMAP
class_components <- function(df) {
  feat <- marker_columns(df)
  if (length(feat) >= 2) {
    mat <- as.matrix(df[, feat, drop = FALSE])
    mat[!is.finite(mat)] <- 0
    pr <- tryCatch(stats::prcomp(mat, center = TRUE, scale. = TRUE), error = function(e) NULL)
    if (!is.null(pr) && ncol(pr$x) >= 2) {
      xy <- pr$x[, 1:2, drop = FALSE]
      colnames(xy) <- c("Component1", "Component2")
      return(xy)
    }
  }
  if (all(c("UMAP1", "UMAP2") %in% names(df))) {
    xy <- as.matrix(df[, c("UMAP1", "UMAP2")])
    colnames(xy) <- c("Component1", "Component2")
    return(xy)
  }
  stop("need marker columns or UMAP1/UMAP2")
}

mst_edges_from_points <- function(xy) {
  xy <- as.matrix(xy)
  n <- nrow(xy)
  if (n < 2) return(matrix(integer(0), ncol = 2))
  d <- as.matrix(stats::dist(xy))
  in_tree <- rep(FALSE, n)
  in_tree[1] <- TRUE
  edges <- matrix(NA_integer_, n - 1L, 2L)
  for (e in seq_len(n - 1L)) {
    best <- Inf
    a <- b <- NA_integer_
    ii <- which(in_tree)
    jj <- which(!in_tree)
    for (i in ii) {
      for (j in jj) {
        if (d[i, j] < best) {
          best <- d[i, j]
          a <- i
          b <- j
        }
      }
    }
    in_tree[b] <- TRUE
    edges[e, ] <- c(a, b)
  }
  edges
}

smooth_edge_curve <- function(xy, idx_a, idx_b, nout = 50) {
  if (!length(idx_a) || !length(idx_b)) return(NULL)
  ca <- colMeans(xy[idx_a, , drop = FALSE])
  cb <- colMeans(xy[idx_b, , drop = FALSE])
  v <- cb - ca
  nrm <- sqrt(sum(v * v))
  if (!is.finite(nrm) || nrm < 1e-6) return(NULL)
  pts <- xy[c(idx_a, idx_b), , drop = FALSE]
  t <- as.numeric((pts[, 1] - ca[1]) * v[1] + (pts[, 2] - ca[2]) * v[2]) / (nrm * nrm)
  o <- order(t)
  t <- t[o]
  pts <- pts[o, , drop = FALSE]
  keep <- !duplicated(round(t, 4))
  t <- t[keep]
  pts <- pts[keep, , drop = FALSE]
  if (length(t) < 8) {
    ts <- seq(0, 1, length.out = nout)
    return(data.frame(x = ca[1] + ts * v[1], y = ca[2] + ts * v[2]))
  }
  df <- min(6, max(3, length(t) - 2L))
  sx <- tryCatch(stats::smooth.spline(t, pts[, 1], df = df), error = function(e) NULL)
  sy <- tryCatch(stats::smooth.spline(t, pts[, 2], df = df), error = function(e) NULL)
  if (is.null(sx) || is.null(sy)) {
    ts <- seq(0, 1, length.out = nout)
    return(data.frame(x = ca[1] + ts * v[1], y = ca[2] + ts * v[2]))
  }
  tg <- seq(min(t), max(t), length.out = nout)
  data.frame(x = as.numeric(stats::predict(sx, tg)$y), y = as.numeric(stats::predict(sy, tg)$y))
}

graph_parent <- function(edges, k, root_i) {
  adj <- vector("list", k)
  for (e in seq_len(nrow(edges))) {
    a <- edges[e, 1]
    b <- edges[e, 2]
    adj[[a]] <- c(adj[[a]], b)
    adj[[b]] <- c(adj[[b]], a)
  }
  parent <- rep(NA_integer_, k)
  dist <- rep(Inf, k)
  parent[root_i] <- 0L
  dist[root_i] <- 0
  q <- root_i
  while (length(q)) {
    u <- q[1]
    q <- q[-1]
    for (v in adj[[u]]) {
      if (is.na(parent[v])) {
        parent[v] <- u
        dist[v] <- dist[u] + 1
        q <- c(q, v)
      }
    }
  }
  list(parent = parent, dist = dist, degree = vapply(adj, length, integer(1)))
}

try_slingshot_fit <- function(xy, cluster, start) {
  if (!requireNamespace("slingshot", quietly = TRUE)) return(NULL)
  cl <- factor(cluster)
  if (!start %in% levels(cl)) return(NULL)
  sds <- tryCatch(
    slingshot::slingshot(xy, clusterLabels = cl, start.clus = start),
    error = function(e) NULL
  )
  if (is.null(sds)) return(NULL)
  crv <- tryCatch(slingshot::slingCurves(sds), error = function(e) NULL)
  if (is.null(crv) || !length(crv)) return(NULL)
  curves <- lapply(seq_along(crv), function(i) {
    s <- crv[[i]]$s
    if (is.null(s) || nrow(s) < 2) return(NULL)
    data.frame(x = s[, 1], y = s[, 2], curve = i, stringsAsFactors = FALSE)
  })
  curves <- Filter(Negate(is.null), curves)
  pt <- tryCatch(slingshot::slingPseudotime(sds), error = function(e) NULL)
  ptv <- if (is.null(pt)) rep(NA_real_, nrow(xy)) else rowMeans(as.matrix(pt), na.rm = TRUE)
  list(curves = curves, pseudotime = ptv, method = "slingshot")
}

fit_lineage_trajectory <- function(xy, cluster, start) {
  xy <- as.matrix(xy)
  if (ncol(xy) != 2) stop("xy must be n x 2")
  cl <- as.character(cluster)
  labs <- unique(cl)
  labs <- labs[nzchar(labs) & !is.na(labs)]
  if (length(labs) < 2 || nrow(xy) < 40) return(NULL)
  if (missing(start) || is.null(start) || !start %in% labs) {
    start <- names(sort(table(cl), decreasing = TRUE))[1]
  }
  sl <- try_slingshot_fit(xy, cl, start)
  cents <- t(vapply(labs, function(L) {
    colMeans(xy[cl == L, , drop = FALSE])
  }, numeric(2)))
  rownames(cents) <- labs
  edges <- mst_edges_from_points(cents)
  gp <- graph_parent(edges, nrow(cents), match(start, labs))
  curves <- list()
  if (is.null(sl)) {
    for (e in seq_len(nrow(edges))) {
      a <- edges[e, 1]
      b <- edges[e, 2]
      cr <- smooth_edge_curve(xy, which(cl == labs[a]), which(cl == labs[b]))
      if (!is.null(cr)) {
        cr$curve <- e
        curves[[length(curves) + 1]] <- cr
      }
    }
    ptv <- rep(NA_real_, nrow(xy))
    for (i in seq_along(labs)) {
      hit <- which(cl == labs[i])
      par <- gp$parent[i]
      if (is.na(par) || par == 0) {
        ptv[hit] <- 0
        next
      }
      ca <- cents[par, ]
      cb <- cents[i, ]
      v <- cb - ca
      nrm2 <- sum(v * v)
      proj <- if (nrm2 < 1e-8) 0 else {
        p <- xy[hit, , drop = FALSE]
        as.numeric((p[, 1] - ca[1]) * v[1] + (p[, 2] - ca[2]) * v[2]) / nrm2
      }
      proj <- pmin(pmax(proj, 0), 1)
      ptv[hit] <- gp$dist[par] + proj
    }
    method <- "mst"
  } else {
    curves <- sl$curves
    ptv <- sl$pseudotime
    method <- sl$method
  }
  branch <- labs[which(gp$degree >= 3)]
  list(
    xy = xy,
    cluster = cl,
    start = start,
    labels = labs,
    centroids = cents,
    edges = edges,
    curves = curves,
    pseudotime = ptv,
    branch_nodes = branch,
    degree = gp$degree,
    method = method
  )
}

arrow_rows <- function(curves) {
  if (!length(curves)) return(NULL)
  rows <- lapply(curves, function(cr) {
    if (nrow(cr) < 8) return(NULL)
    i1 <- max(2L, as.integer(round(0.55 * nrow(cr))))
    i2 <- min(nrow(cr), i1 + max(3L, as.integer(round(0.12 * nrow(cr)))))
    data.frame(x = cr$x[i1], y = cr$y[i1], xend = cr$x[i2], yend = cr$y[i2])
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

plot_trajectory_tree <- function(df, fit, panel_id, title, facet_group = FALSE) {
  plot_df <- df
  plot_df$Component1 <- fit$xy[, 1]
  plot_df$Component2 <- fit$xy[, 2]
  plot_df$celltype <- celltype_label(plot_df$lineage, panel_id)
  levs <- unique(plot_df$celltype)
  plot_df$celltype <- factor(plot_df$celltype, levels = levs)
  if ("group" %in% names(plot_df)) {
    plot_df$group <- factor(plot_df$group, levels = c("T", "T6"))
  }
  pal <- celltype_colors(levels(plot_df$celltype))
  sk <- if (length(fit$curves)) {
    do.call(rbind, lapply(seq_along(fit$curves), function(i) {
      cr <- fit$curves[[i]]
      cr$curve <- i
      cr
    }))
  } else {
    NULL
  }
  arr <- arrow_rows(fit$curves)
  lab <- NULL
  if (length(fit$branch_nodes)) {
    lab <- data.frame(
      Component1 = fit$centroids[fit$branch_nodes, 1],
      Component2 = fit$centroids[fit$branch_nodes, 2],
      label = celltype_label(fit$branch_nodes, panel_id),
      stringsAsFactors = FALSE
    )
  }
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Component1, y = Component2, color = celltype)) +
    ggplot2::geom_point(size = if (facet_group) 0.5 else 0.55, alpha = 0.8, stroke = 0)
  if (!is.null(sk) && nrow(sk)) {
    p <- p + ggplot2::geom_path(
      data = sk, ggplot2::aes(x = x, y = y, group = curve),
      color = "black", linewidth = 0.85, inherit.aes = FALSE
    )
  }
  if (!is.null(arr) && nrow(arr)) {
    p <- p + ggplot2::geom_segment(
      data = arr, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      arrow = grid::arrow(length = grid::unit(0.11, "inches"), type = "closed"),
      color = "black", linewidth = 0.45, inherit.aes = FALSE
    )
  }
  if (!is.null(lab) && nrow(lab)) {
    p <- p + ggplot2::geom_text(
      data = lab, ggplot2::aes(x = Component1, y = Component2, label = label),
      color = "black", fontface = "bold", size = 3.3, vjust = -1.1,
      inherit.aes = FALSE
    )
  }
  if (facet_group && "group" %in% names(plot_df)) {
    p <- p + ggplot2::facet_wrap(~group, ncol = 2)
  }
  p +
    ggplot2::scale_color_manual(values = pal, drop = FALSE) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 2)) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 14),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    ) +
    ggplot2::labs(title = title, x = "Component 1", y = "Component 2", color = NULL)
}

plot_pseudotime_group <- function(df, title) {
  ggplot2::ggplot(df, ggplot2::aes(x = group, y = pseudotime, fill = group)) +
    ggplot2::geom_boxplot(outlier.size = 0.4, width = 0.55) +
    ggplot2::scale_fill_manual(values = pal_group) +
    theme_dr() +
    ggplot2::labs(title = title, x = NULL, y = "Pseudotime", fill = NULL)
}

export_one_major_trajectory <- function(cells, panel_id, major, out_dir) {
  sub <- cells[as.character(cells$major) == major, , drop = FALSE]
  n_lin <- length(unique(sub$lineage))
  if (nrow(sub) < 80 || n_lin < 2) {
    log_msg(panel_id, " ", major, ": skip trajectory (n=", nrow(sub), ", subsets=", n_lin, ")")
    return(invisible(NULL))
  }
  cap <- 4000L
  if (nrow(sub) > cap) {
    set.seed(if (exists("seed_value")) seed_value else 42)
    sub <- sub[sample.int(nrow(sub), cap), , drop = FALSE]
  }
  xy <- class_components(sub)
  root <- default_trajectory_root(panel_id, major, sub$lineage)
  fit <- fit_lineage_trajectory(xy, sub$lineage, root)
  if (is.null(fit)) {
    log_msg(panel_id, " ", major, ": trajectory fit failed")
    return(invisible(NULL))
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tag <- paste0(panel_id, "_", gsub("[^A-Za-z0-9_-]", "_", major))
  sub$pseudotime <- fit$pseudotime
  sub$Component1 <- fit$xy[, 1]
  sub$Component2 <- fit$xy[, 2]
  keep <- intersect(
    c("sample", "group", "lineage", "major", "Component1", "Component2", "pseudotime", "UMAP1", "UMAP2"),
    names(sub)
  )
  utils::write.csv(sub[, keep, drop = FALSE],
                   file.path(out_dir, paste0(tag, "_pseudotime.csv")), row.names = FALSE)
  ttl <- paste0(panel_id, "  ", major, "  trajectory  (root: ",
                celltype_label(root, panel_id), ")")
  save_gg(plot_trajectory_tree(sub, fit, panel_id, ttl, facet_group = FALSE),
          file.path(out_dir, paste0(tag, "_trajectory")),
          width = 7.2, height = 6.4)
  save_gg(plot_trajectory_tree(sub, fit, panel_id, paste0(ttl, "  T | T6"), facet_group = TRUE),
          file.path(out_dir, paste0(tag, "_trajectory_T_vs_T6")),
          width = 11.6, height = 5.8)
  if (any(is.finite(sub$pseudotime))) {
    save_gg(plot_pseudotime_group(sub, paste0(panel_id, "  ", major, "  pseudotime T vs T6")),
            file.path(out_dir, paste0(tag, "_pseudotime_T_vs_T6")),
            width = 5.2, height = 4.6)
  }
  log_msg(panel_id, " ", major, " trajectory (", fit$method, ", root=", root,
          ", n=", nrow(sub), ") -> ", out_dir)
  invisible(fit)
}

read_panel_cells_for_trajectory <- function(result_dir, panel_id) {
  p <- file.path(result_dir, panel_id, paste0(panel_id, "_cell_embeddings.csv"))
  if (!file.exists(p)) return(NULL)
  df <- read_embed_csv(p)
  if (!all(c("lineage", "UMAP1", "UMAP2") %in% names(df))) return(NULL)
  df
}

export_panel_trajectories <- function(result_dir, panel_id) {
  cells <- read_panel_cells_for_trajectory(result_dir, panel_id)
  if (is.null(cells) || nrow(cells) < 80) {
    log_msg(panel_id, ": no embeddings for trajectory")
    return(invisible(NULL))
  }
  cells$major <- majors_for_trajectory(panel_id, cells)
  out_dir <- file.path(result_dir, panel_id, "trajectory")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    "Trajectory is a marker-state continuum inside one major class of one panel.",
    "P1/P2/P3 are not concatenated. Root: naive T/B or Ly6C-hi monocytes by convention.",
    "This is not proof of developmental origin. Neutrophils on the myeloid tree are a branch, not derived from monocytes.",
    paste("Method: slingshot if installed, otherwise MST of subset centroids + smooth curves.")
  ), file.path(out_dir, "TRAJECTORY_NOTE.txt"))
  majors <- setdiff(unique(cells$major), c("dump", "", NA))
  for (mj in majors) {
    tryCatch(
      export_one_major_trajectory(cells, panel_id, mj, out_dir),
      error = function(e) log_msg(panel_id, " ", mj, " trajectory failed: ", e$message)
    )
  }
  invisible(out_dir)
}

export_all_panel_trajectories <- function(result_dir) {
  if (missing(result_dir) || !nzchar(result_dir)) {
    if (exists("result_dir", envir = .GlobalEnv, inherits = FALSE)) {
      result_dir <- get("result_dir", envir = .GlobalEnv)
    } else {
      result_dir <- file.path(getwd(), "results_flow")
    }
  }
  for (pn in c("P1", "P2", "P3")) {
    export_panel_trajectories(result_dir, pn)
  }
  log_msg("Trajectories written under results_flow/P1|P2|P3/trajectory/")
  invisible(TRUE)
}

load_flow_pipeline_functions()

if (!identical(toupper(Sys.getenv("FLOW_TRAJECTORY_FROM_PIPELINE", "0")), "1") &&
    !identical(toupper(Sys.getenv("FLOW_FUNCTIONS_ONLY", "0")), "1")) {
  export_all_panel_trajectories(if (exists("result_dir")) result_dir else file.path(getwd(), "results_flow"))
}
