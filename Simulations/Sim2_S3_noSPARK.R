rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(methods)
library(Matrix)
library(RANN)
library(gaston)
library(SPARK)
library(matrixStats)
library(INSiGHT)

set.seed(111)

dat_paths <- c(
  "151507" = "./DLPFC/data/DLPFC_151507_151508_151509_151510_combined_with_ground_truth.RDS",
  "151669" = "./DLPFC/data/DLPFC_151669_151670_151671_151672_combined_with_ground_truth.RDS",
  "151671" = "./DLPFC/data/DLPFC_151669_151670_151671_151672_combined_with_ground_truth.RDS"
)
sample_ids <- c("151507", "151669", "151671")
numSignalGenes <- 5000
numGenes <- 10000
nSim <- 1
target_num_spots <- 1000
nb_size_list <- c(1.5, 3)

sample_mu <- data.frame(
  sample_id = sample_ids,
  mu_wm = c(0.25, 0.4, 0.6),
  mu_outside = c(1, 1, 1.2),
  stringsAsFactors = FALSE
)
sim_params <- data.frame(iCombo = 1, size = nb_size_list)

out_file <- "./Simulations/out/Sim2_S3Out_DLPFC_151507_151669_151671.RData"

# Simulated gene pattern plot settings: choose one mu combination and one size.
gene_to_display <- "Gene1"
gene_expression_scale <- "raw" # "raw" or "scaled"
gene_plot_iCombo <- 1
gene_plot_size <- 3
sim_gene_plot_path <- "./Simulations/results/Sim2_S3_simulated_gene_pattern.jpeg"
gene_plot_axis_text_size <- 20
gene_plot_strip_text_size <- 30
gene_plot_legend_text_size <- 20
gene_plot_legend_title_size <- 20
gene_plot_title_size <- 20
gene_plot_subtitle_size <- 20
gene_plot_region_legend_point_size <- 8

simulate_counts <- function(num_genes, num_signal_genes, wm, nb_wm, nb_outside) {
  num_spots <- length(wm)
  outside_wm <- !wm
  counts <- matrix(0, nrow = num_genes, ncol = num_spots)

  counts[, outside_wm] <- matrix(
    rnbinom(
      n = num_genes * sum(outside_wm),
      size = nb_outside[["size"]],
      mu = nb_outside[["mu"]]
    ),
    nrow = num_genes
  )

  counts[seq_len(num_signal_genes), wm] <- matrix(
    rnbinom(
      n = num_signal_genes * sum(wm),
      size = nb_wm[["size"]],
      mu = nb_wm[["mu"]]
    ),
    nrow = num_signal_genes
  )

  if (num_signal_genes < num_genes) {
    null_genes <- (num_signal_genes + 1):num_genes
    counts[null_genes, wm] <- matrix(
      rnbinom(
        n = length(null_genes) * sum(wm),
        size = nb_outside[["size"]],
        mu = nb_outside[["mu"]]
      ),
      nrow = length(null_genes)
    )
  }

  rownames(counts) <- paste0("Gene", seq_len(num_genes))
  colnames(counts) <- paste0("Spot", seq_len(num_spots))
  counts
}

plot_simulated_gene <- function(plot_designs, gene, size,
                                expression_scale = c("scaled", "raw"), out_path = NULL) {
  expression_scale <- match.arg(expression_scale)
  first_counts <- plot_designs[[1]]$counts
  gene_name <- if (is.numeric(gene)) rownames(first_counts)[gene] else gene
  if (is.na(gene_name) || !gene_name %in% rownames(first_counts)) {
    stop("Gene ", gene, " was not found in the simulated count matrix.")
  }

  raw_expression_by_design <- lapply(plot_designs, function(design) {
    counts <- design$counts
    if (!gene_name %in% rownames(counts)) {
      stop("Gene ", gene_name, " was not found in every simulated count matrix.")
    }
    as.numeric(counts[gene_name, ])
  })
  raw_expression <- unlist(raw_expression_by_design, use.names = FALSE)
  expression_range <- diff(range(raw_expression))
  plotted_expression_by_design <- if (expression_scale == "raw") {
    raw_expression_by_design
  } else if (expression_range == 0) {
    lapply(raw_expression_by_design, function(expression) rep(0, length(expression)))
  } else {
    lapply(
      raw_expression_by_design,
      function(expression) (expression - min(raw_expression)) / expression_range
    )
  }

  plot_df <- do.call(
    base::rbind,
    lapply(seq_along(plot_designs), function(iDesign) {
      design <- plot_designs[[iDesign]]
      data.frame(
        x = -design$pos[, 1],
        y = design$pos[, 2],
        expression = plotted_expression_by_design[[iDesign]],
        domain = ifelse(design$wm, "WM", "cortical"),
        panel = design$sample_id,
        stringsAsFactors = FALSE
      )
    })
  )
  plot_df$panel <- factor(plot_df$panel, levels = unique(plot_df$panel))

  p <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(aes(fill = expression, color = domain), shape = 21, size = 1.3, stroke = 0.25) +
    coord_fixed() +
    scale_fill_gradient(low = "white", high = "steelblue4") +
    scale_color_manual(values = c(WM = "firebrick3", cortical = "darkseagreen3")) +
    facet_wrap(
      ~panel,
      ncol = length(unique(vapply(plot_designs, `[[`, character(1), "sample_id"))),
      labeller = label_value
    ) +
    labs(
      x = "spatial x",
      y = "spatial y",
      fill = ifelse(expression_scale == "raw", "raw expression", "scaled expression"),
      color = "region"
    ) +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size = gene_plot_region_legend_point_size,
                                                    stroke = 1.8))) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      strip.text = element_text(size = gene_plot_strip_text_size),
      legend.position = "none",
      legend.text = element_text(size = gene_plot_legend_text_size),
      legend.title = element_text(size = gene_plot_legend_title_size, face = "bold")
    )

  if (!is.null(out_path)) {
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(p, file = out_path, width = 9, height = 4, dpi = 500)
  }
  p
}

assign_grid_cells <- function(pos, nx, ny) {
  x_bin <- pmin(floor(pos[, 1] * nx) + 1, nx)
  y_bin <- pmin(floor(pos[, 2] * ny) + 1, ny)
  paste(x_bin, y_bin, sep = "_")
}

aggregate_spots_by_grid <- function(pos, wm, target_num_spots) {
  if (target_num_spots >= nrow(pos)) {
    rownames(pos) <- paste0("Spot", seq_len(nrow(pos)))
    return(list(pos = pos, wm = wm, grid = NULL, original_num_spots = nrow(pos)))
  }

  target_side <- sqrt(target_num_spots)
  candidate_dims <- expand.grid(
    nx = seq_len(ceiling(target_side * 2)),
    ny = seq_len(ceiling(target_side * 2))
  )
  candidate_dims <- candidate_dims[candidate_dims$nx * candidate_dims$ny >= target_num_spots, ]
  candidate_dims$n_occupied <- vapply(
    seq_len(nrow(candidate_dims)),
    function(i) {
      length(unique(assign_grid_cells(pos, candidate_dims$nx[i], candidate_dims$ny[i])))
    },
    numeric(1)
  )
  candidate_dims$target_diff <- abs(candidate_dims$n_occupied - target_num_spots)
  candidate_dims$cell_diff <- abs(candidate_dims$nx * candidate_dims$ny - target_num_spots)
  best_dim <- candidate_dims[order(candidate_dims$target_diff, candidate_dims$cell_diff), ][1, ]

  grid_id <- assign_grid_cells(pos, best_dim$nx, best_dim$ny)
  grid_levels <- sort(unique(grid_id))
  aggregated_pos <- do.call(
    rbind,
    lapply(grid_levels, function(id) colMeans(pos[grid_id == id, , drop = FALSE]))
  )
  aggregated_wm <- vapply(
    grid_levels,
    function(id) any(wm[grid_id == id]),
    logical(1)
  )

  rownames(aggregated_pos) <- paste0("Spot", seq_len(nrow(aggregated_pos)))
  names(aggregated_wm) <- rownames(aggregated_pos)

  list(
    pos = aggregated_pos,
    wm = aggregated_wm,
    grid = data.frame(
      original_spot = rownames(pos),
      grid_id = grid_id,
      stringsAsFactors = FALSE
    ),
    original_num_spots = nrow(pos),
    target_num_spots = target_num_spots,
    aggregated_num_spots = nrow(aggregated_pos),
    nx = best_dim$nx,
    ny = best_dim$ny
  )
}

run_svg_package_method <- function(counts, pos, method = c("MERINGUE", "binSpect")) {
  if (!requireNamespace("SVG", quietly = TRUE)) {
    stop("Package SVG is required to run MERINGUE and binSpect.")
  }

  method <- match.arg(method)
  expr_norm <- Seurat::NormalizeData(counts)
  if (method == "MERINGUE") {
    result <- SVG::CalSVG_MERINGUE(
      expr_matrix = as.matrix(expr_norm),
      spatial_coords = as.matrix(pos),
      network_method = "delaunay",
      k = 5,
      n_threads = 1,
      adjust_method = "BY"
    )
  } else {
    result <- SVG::CalSVG_binSpect(
      expr_matrix = as.matrix(expr_norm),
      spatial_coords = as.matrix(pos),
      n_threads = 1,
      adjust_method = "BY"
    )
  }

  idx <- match(rownames(counts), result$gene)
  p_values <- result$p.value[idx]
  p_values[is.na(p_values)] <- 1
  p_values
}

check_benchmark_packages <- function() {
  required_packages <- c("scBSP", "SVG", "Seurat")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0) {
    stop(
      "Missing packages required for the full benchmark: ",
      paste(missing_packages, collapse = ", ")
    )
  }
}

combine_p_values <- function(p_list, method = c("union", "intersection", "cauchy")) {
  method <- match.arg(method)
  p_mat <- do.call(cbind, p_list)

  if (method == "union") {
    return(do.call(pmin, as.data.frame(p_mat)))
  }
  if (method == "intersection") {
    return(do.call(pmax, as.data.frame(p_mat)))
  }

  eps <- 1e-15
  p_mat <- pmin(pmax(p_mat, eps), 1 - eps)
  stat <- rowMeans(tan((0.5 - p_mat) * pi))
  p_combined <- 0.5 - atan(stat) / pi
  pmin(pmax(p_combined, 0), 1)
}

combine_method_results <- function(p_list, prefix) {
  data.frame(
    setNames(
      list(
        combine_p_values(p_list, "union"),
        combine_p_values(p_list, "intersection"),
        combine_p_values(p_list, "cauchy")
      ),
      paste0(prefix, c("_union", "_intersection", "_cauchy"))
    ),
    check.names = FALSE
  )
}

get_sample_pattern <- function(dat, sample_id, target_num_spots) {
  iSample <- match(sample_id, dat$samples.id)
  if (is.na(iSample)) {
    stop("Sample ", sample_id, " was not found in dat$samples.id.")
  }

  labels <- dat$ground_truth[[iSample]]
  annotated_spots <- !is.na(labels)
  wm <- labels[annotated_spots] == "WM"

  pos <- dat$spatial.locs[[iSample]][annotated_spots, , drop = FALSE]
  pos <- as.matrix(pos)
  pos[, 1] <- (pos[, 1] - min(pos[, 1])) / (max(pos[, 1]) - min(pos[, 1]))
  pos[, 2] <- (pos[, 2] - min(pos[, 2])) / (max(pos[, 2]) - min(pos[, 2]))
  rownames(pos) <- paste0("Spot", seq_len(nrow(pos)))

  grid_aggregation <- aggregate_spots_by_grid(
    pos = pos,
    wm = wm,
    target_num_spots = target_num_spots
  )

  list(
    sample_id = sample_id,
    pos = grid_aggregation$pos,
    wm = grid_aggregation$wm,
    grid_aggregation = grid_aggregation
  )
}

run_single_sample_methods <- function(counts, pos) {
  object_sparkx <- sparkx(counts, pos, option = "single", verbose = FALSE)
  p_values_sparkx <- object_sparkx$res_mtest$combinedPval
  rm(object_sparkx)

  p_values_meringue <- run_svg_package_method(counts, pos, method = "MERINGUE")
  p_values_binspect <- run_svg_package_method(counts, pos, method = "binSpect")

  object_scBSP <- scBSP::scBSP(Coords = pos, ExpMat_Sp = as(counts, "dgCMatrix"))
  p_values_scbsp <- object_scBSP$P_values
  rm(object_scBSP)

  list(
    SPARK_X = p_values_sparkx,
    MERINGUE = p_values_meringue,
    binSpect = p_values_binspect,
    scBSP = p_values_scbsp
  )
}

make_mu_metadata <- function(sample_mu) {
  values <- unlist(
    lapply(seq_len(nrow(sample_mu)), function(i) {
      setNames(
        list(sample_mu$mu_wm[i], sample_mu$mu_outside[i]),
        c(
          paste0("mu_wm_", sample_mu$sample_id[i]),
          paste0("mu_outside_", sample_mu$sample_id[i])
        )
      )
    }),
    recursive = FALSE
  )
  as.data.frame(values, check.names = FALSE)
}

dat_by_sample <- lapply(dat_paths, readRDS)
sample_patterns <- setNames(
  lapply(sample_ids, function(sample_id) {
    get_sample_pattern(dat_by_sample[[sample_id]], sample_id, target_num_spots)
  }),
  sample_ids
)

grid_aggregation <- lapply(sample_patterns, `[[`, "grid_aggregation")
for (sample_id in sample_ids) {
  cat(paste0(
    "Aggregated sample ", sample_id, " from ",
    grid_aggregation[[sample_id]]$original_num_spots,
    " annotated spots into ",
    grid_aggregation[[sample_id]]$aggregated_num_spots,
    " grid spots.\n"
  ))
}

gene_plot_param <- sim_params[
  sim_params$iCombo == gene_plot_iCombo & sim_params$size == gene_plot_size,
  ,
  drop = FALSE
]
if (nrow(gene_plot_param) != 1) {
  stop("gene_plot_iCombo and gene_plot_size must identify exactly one row of sim_params.")
}
rng_state_before_gene_plot <- .Random.seed
gene_plot_designs <- lapply(seq_len(nrow(sample_mu)), function(iSample) {
    sample_id <- sample_mu$sample_id[iSample]
    list(
      counts = simulate_counts(
        num_genes = numGenes,
        num_signal_genes = numSignalGenes,
        wm = sample_patterns[[sample_id]]$wm,
        nb_wm = c(size = gene_plot_param$size, mu = sample_mu$mu_wm[iSample]),
        nb_outside = c(size = gene_plot_param$size, mu = sample_mu$mu_outside[iSample])
      ),
      pos = sample_patterns[[sample_id]]$pos,
      wm = sample_patterns[[sample_id]]$wm,
      sample_id = sample_id,
      mu_wm = sample_mu$mu_wm[iSample],
      mu_outside = sample_mu$mu_outside[iSample]
    )
  })
plot_simulated_gene(
  plot_designs = gene_plot_designs,
  gene = gene_to_display,
  size = gene_plot_size,
  expression_scale = gene_expression_scale,
  out_path = sim_gene_plot_path
)
.Random.seed <- rng_state_before_gene_plot
rm(gene_plot_designs, gene_plot_param, rng_state_before_gene_plot)

nb_params <- do.call(
  base::rbind,
  lapply(seq_len(nrow(sim_params)), function(iParam) {
    param <- sim_params[iParam, , drop = FALSE]
    do.call(
      base::rbind,
      lapply(seq_len(nrow(sample_mu)), function(iSample) {
        data.frame(
          iCombo = param$iCombo,
          sample_id = sample_mu$sample_id[iSample],
          region = c("WM", "cortical"),
          size = param$size,
          mu = c(sample_mu$mu_wm[iSample], sample_mu$mu_outside[iSample]),
          stringsAsFactors = FALSE
        )
      })
    )
  })
)
print(nb_params)

method_columns <- c(
  "INSiGHT",
  paste0(rep(c("SPARK_X", "MERINGUE", "binSpect", "scBSP"), each = 3),
         c("_union", "_intersection", "_cauchy"))
)
out_null <- data.frame(matrix(nrow = 0, ncol = length(method_columns)))
colnames(out_null) <- method_columns
out_null$size <- numeric(0)
out_null$iCombo <- integer(0)
mu_template <- make_mu_metadata(sample_mu)
for (column_name in colnames(mu_template)) {
  out_null[[column_name]] <- numeric(0)
}
out_null$iSim <- integer(0)
out_power <- out_null

if (nSim > 0) {
  check_benchmark_packages()
}

for (iParam in seq_len(nrow(sim_params))) {
  param <- sim_params[iParam, , drop = FALSE]

  for (iSim in seq_len(nSim)) {
    cat(paste0(
      "Working on combo ", param$iCombo,
      ", size = ", param$size,
      ", simulation ", iSim, ".\n"
    ))

    counts_list <- setNames(
      lapply(seq_len(nrow(sample_mu)), function(iSample) {
        sample_id <- sample_mu$sample_id[iSample]
        simulate_counts(
          num_genes = numGenes,
          num_signal_genes = numSignalGenes,
          wm = sample_patterns[[sample_id]]$wm,
          nb_wm = c(size = param$size, mu = sample_mu$mu_wm[iSample]),
          nb_outside = c(size = param$size, mu = sample_mu$mu_outside[iSample])
        )
      }),
      sample_ids
    )
    spatial_locs <- lapply(sample_ids, function(sample_id) sample_patterns[[sample_id]]$pos)

    object_INSiGHT <- createINSiGHTobject(
      data.input = counts_list,
      spatial.locs = spatial_locs
    )
    object_INSiGHT <- processINSiGHT(
      object = object_INSiGHT,
      k.NN4bandwidth = 36,
      quantile.bandwidth = 0.3,
      kernel.option = "single",
      rank.tie = "average",
      NN4Kernel = FALSE,
      k.NN4Kernel = 100,
      normalization = "none"
    )
    object_INSiGHT <- testINSiGHT(object = object_INSiGHT)
    p_values_INSiGHT <- object_INSiGHT@testResults$p.values_combined
    rm(object_INSiGHT)

    p_values_by_sample <- setNames(
      lapply(sample_ids, function(sample_id) {
        run_single_sample_methods(
          counts = counts_list[[sample_id]],
          pos = sample_patterns[[sample_id]]$pos
        )
      }),
      sample_ids
    )

    combined_p <- cbind(
      INSiGHT = p_values_INSiGHT,
      combine_method_results(lapply(p_values_by_sample, `[[`, "SPARK_X"), "SPARK_X"),
      combine_method_results(lapply(p_values_by_sample, `[[`, "MERINGUE"), "MERINGUE"),
      combine_method_results(lapply(p_values_by_sample, `[[`, "binSpect"), "binSpect"),
      combine_method_results(lapply(p_values_by_sample, `[[`, "scBSP"), "scBSP")
    )

    mu_metadata <- make_mu_metadata(sample_mu)

    output_null <- data.frame(
      combined_p[(numSignalGenes + 1):numGenes, , drop = FALSE],
      size = rep(param$size, numGenes - numSignalGenes),
      iCombo = rep(param$iCombo, numGenes - numSignalGenes),
      mu_metadata[rep(1, numGenes - numSignalGenes), , drop = FALSE],
      iSim = rep(iSim, numGenes - numSignalGenes),
      check.names = FALSE
    )
    output_power <- data.frame(
      combined_p[seq_len(numSignalGenes), , drop = FALSE],
      size = rep(param$size, numSignalGenes),
      iCombo = rep(param$iCombo, numSignalGenes),
      mu_metadata[rep(1, numSignalGenes), , drop = FALSE],
      iSim = rep(iSim, numSignalGenes),
      check.names = FALSE
    )

    out_null <- base::rbind(out_null, output_null)
    out_power <- base::rbind(out_power, output_power)
    save(
      out_null,
      out_power,
      nb_params,
      grid_aggregation,
      sim_params,
      sample_mu,
      file = out_file
    )
  }
}
