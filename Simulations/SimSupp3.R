rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(methods)
library(Matrix)
library(RANN)
library(gaston)
library(matrixStats)
library(INSiGHT)

set.seed(1)

dat_paths <- c(
  "151509" = "./DLPFC/data/DLPFC_151507_151508_151509_151510_combined_with_ground_truth.RDS",
  "151669" = "./DLPFC/data/DLPFC_151669_151670_151671_151672_combined_with_ground_truth.RDS"
)
sample_ids <- c("151509", "151669")
target_num_spots <- c("151509" = 1000, "151669" = 1000)

numSignalGenes <- 5000
numGenes <- 10000
nSim <- 1

mu_wm_151509 <- 0.25
mu_outside_151509 <- 1
mu_wm_151669 <- 0.25
mu_outside_151669 <- 1
nb_size <- 1.5

random_library_sigma <- 0.5
spatial_library_sigma <- 0.1
spatial_library_factor_wm <- 0.4
spatial_library_factor_cortical <- 1.2

library_scenarios <- data.frame(
  scenario = c("none", "random", "spatial"),
  scenario_label = c("No library-size effect", "Random library-size effect", "Spatial library-size effect"),
  random_sigma = c(0, random_library_sigma, spatial_library_sigma),
  spatial_type = c("none", "none", "region"),
  spatial_factor_wm = c(1, 1, spatial_library_factor_wm),
  spatial_factor_cortical = c(1, 1, spatial_library_factor_cortical),
  stringsAsFactors = FALSE
)

sim_params <- data.frame(
  iScenario = seq_len(nrow(library_scenarios)),
  library_scenarios,
  mu_wm_151509 = mu_wm_151509,
  mu_outside_151509 = mu_outside_151509,
  mu_wm_151669 = mu_wm_151669,
  mu_outside_151669 = mu_outside_151669,
  size = nb_size,
  stringsAsFactors = FALSE
)

out_file <- "./Simulations/out/SimSupp3Out_DLPFC_151509_151669.RData"

# Simulated gene pattern plot settings: choose one library-size scenario.
gene_to_display <- "Gene1"
gene_expression_scale <- "raw" # "raw" or "scaled"
gene_plot_iScenario <- 3
sim_gene_plot_path <- "./Simulations/results/SimSupp3_simulated_gene_pattern.jpeg"
gene_plot_strip_text_size <- 24
gene_plot_region_legend_point_size <- 8

simulate_counts <- function(num_genes, num_signal_genes, wm, nb_wm, nb_outside, library_size) {
  num_spots <- length(wm)
  outside_wm <- !wm
  counts <- matrix(0, nrow = num_genes, ncol = num_spots)

  counts[, outside_wm] <- matrix(
    rnbinom(
      n = num_genes * sum(outside_wm),
      size = nb_outside[["size"]],
      mu = rep(nb_outside[["mu"]] * library_size[outside_wm], each = num_genes)
    ),
    nrow = num_genes
  )

  counts[seq_len(num_signal_genes), wm] <- matrix(
    rnbinom(
      n = num_signal_genes * sum(wm),
      size = nb_wm[["size"]],
      mu = rep(nb_wm[["mu"]] * library_size[wm], each = num_signal_genes)
    ),
    nrow = num_signal_genes
  )

  if (num_signal_genes < num_genes) {
    null_genes <- (num_signal_genes + 1):num_genes
    counts[null_genes, wm] <- matrix(
      rnbinom(
        n = length(null_genes) * sum(wm),
        size = nb_outside[["size"]],
        mu = rep(nb_outside[["mu"]] * library_size[wm], each = length(null_genes))
      ),
      nrow = length(null_genes)
    )
  }

  rownames(counts) <- paste0("Gene", seq_len(num_genes))
  colnames(counts) <- paste0("Spot", seq_len(num_spots))
  counts
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

simulate_library_size <- function(pattern, scenario_row) {
  n <- nrow(pattern$pos)
  library_size <- rep(1, n)

  if (scenario_row$random_sigma > 0) {
    sigma <- scenario_row$random_sigma
    library_size <- library_size * exp(rnorm(n, mean = -0.5 * sigma^2, sd = sigma))
  }

  if (scenario_row$spatial_type == "region") {
    spatial_factor <- ifelse(
      pattern$wm,
      scenario_row$spatial_factor_wm,
      scenario_row$spatial_factor_cortical
    )
    spatial_factor <- spatial_factor / mean(spatial_factor)
    library_size <- library_size * spatial_factor
  }

  library_size / mean(library_size)
}

normalize_counts <- function(counts, library_size, method = c("raw", "normalized")) {
  method <- match.arg(method)
  if (method == "raw") {
    return(counts)
  }

  normalized <- sweep(counts, 2, library_size, FUN = "/")

  rownames(normalized) <- rownames(counts)
  colnames(normalized) <- colnames(counts)
  normalized
}

run_insight <- function(counts_list, pos_list) {
  object <- createINSiGHTobject(
    data.input = counts_list,
    spatial.locs = pos_list,
    sample.names = sample_ids
  )
  object <- processINSiGHT(
    object = object,
    k.NN4bandwidth = 36,
    quantile.bandwidth = 0.3,
    kernel.option = "single",
    rank.tie = "average",
    rank.zero = FALSE,
    NN4Kernel = FALSE,
    k.NN4Kernel = 100,
    normalization = "none"
  )
  object <- testINSiGHT(object = object)
  object@testResults$p.values_combined
}

plot_simulated_gene <- function(plot_designs, gene, expression_scale = c("scaled", "raw"), out_path = NULL) {
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
        library_size = design$library_size,
        domain = ifelse(design$wm, "WM", "cortical"),
        panel = design$sample_id,
        stringsAsFactors = FALSE
      )
    })
  )
  plot_df$panel <- factor(plot_df$panel, levels = unique(plot_df$panel))

  p <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(aes(fill = expression, color = domain), shape = 21, size = 1.5, stroke = 0.25) +
    coord_fixed() +
    scale_fill_gradient(low = "white", high = "steelblue4") +
    scale_color_manual(values = c(WM = "firebrick3", cortical = "darkseagreen3")) +
    facet_wrap(~panel, ncol = length(plot_designs), labeller = label_value) +
    labs(
      x = "spatial x",
      y = "spatial y",
      fill = ifelse(expression_scale == "raw", "raw expression", "scaled expression"),
      color = "region"
    ) +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size = gene_plot_region_legend_point_size, stroke = 1.8))) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      strip.text = element_text(size = gene_plot_strip_text_size),
      legend.position = "none"
    )

  if (!is.null(out_path)) {
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(p, file = out_path, width = 6, height = 4, dpi = 500)
  }
  p
}

dat_by_sample <- lapply(dat_paths, readRDS)
sample_patterns <- setNames(
  lapply(sample_ids, function(sample_id) {
    get_sample_pattern(dat_by_sample[[sample_id]], sample_id, target_num_spots[[sample_id]])
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

gene_plot_param <- sim_params[sim_params$iScenario == gene_plot_iScenario, , drop = FALSE]
if (nrow(gene_plot_param) != 1) {
  stop("gene_plot_iScenario must identify exactly one row of sim_params.")
}
rng_state_before_gene_plot <- .Random.seed
gene_plot_designs <- lapply(sample_ids, function(sample_id) {
  mu_wm <- gene_plot_param[[paste0("mu_wm_", sample_id)]]
  mu_outside <- gene_plot_param[[paste0("mu_outside_", sample_id)]]
  library_size <- simulate_library_size(sample_patterns[[sample_id]], gene_plot_param)
  list(
    counts = simulate_counts(
      num_genes = numGenes,
      num_signal_genes = numSignalGenes,
      wm = sample_patterns[[sample_id]]$wm,
      nb_wm = c(size = gene_plot_param$size, mu = mu_wm),
      nb_outside = c(size = gene_plot_param$size, mu = mu_outside),
      library_size = library_size
    ),
    pos = sample_patterns[[sample_id]]$pos,
    wm = sample_patterns[[sample_id]]$wm,
    sample_id = sample_id,
    library_size = library_size
  )
})
plot_simulated_gene(
  plot_designs = gene_plot_designs,
  gene = gene_to_display,
  expression_scale = gene_expression_scale,
  out_path = sim_gene_plot_path
)
.Random.seed <- rng_state_before_gene_plot
rm(gene_plot_designs, gene_plot_param, rng_state_before_gene_plot)

nb_params <- do.call(
  base::rbind,
  lapply(seq_len(nrow(sim_params)), function(iParam) {
    param <- sim_params[iParam, , drop = FALSE]
    base::rbind(
      data.frame(
        iScenario = param$iScenario,
        scenario = param$scenario,
        sample_id = "151509",
        region = c("WM", "cortical"),
        size = param$size,
        mu = c(param$mu_wm_151509, param$mu_outside_151509),
        stringsAsFactors = FALSE
      ),
      data.frame(
        iScenario = param$iScenario,
        scenario = param$scenario,
        sample_id = "151669",
        region = c("WM", "cortical"),
        size = param$size,
        mu = c(param$mu_wm_151669, param$mu_outside_151669),
        stringsAsFactors = FALSE
      )
    )
  })
)
print(nb_params)

out_null <- data.frame(
  INSiGHT_raw = NULL,
  INSiGHT_normalized = NULL,
  size = NULL,
  iScenario = NULL,
  scenario = NULL,
  scenario_label = NULL,
  mu_wm_151509 = NULL,
  mu_outside_151509 = NULL,
  mu_wm_151669 = NULL,
  mu_outside_151669 = NULL,
  n_spots_151509 = NULL,
  n_spots_151669 = NULL,
  library_cv_151509 = NULL,
  library_cv_151669 = NULL,
  zero_prop_151509 = NULL,
  zero_prop_151669 = NULL,
  iSim = NULL
)
out_power <- out_null

for (iParam in seq_len(nrow(sim_params))) {
  param <- sim_params[iParam, , drop = FALSE]

  for (iSim in seq_len(nSim)) {
    cat(paste0(
      "Working on library-size scenario ", param$scenario,
      ", size = ", param$size,
      ", simulation ", iSim, ".\n"
    ))

    library_size_151509 <- simulate_library_size(sample_patterns[["151509"]], param)
    library_size_151669 <- simulate_library_size(sample_patterns[["151669"]], param)

    counts_151509 <- simulate_counts(
      num_genes = numGenes,
      num_signal_genes = numSignalGenes,
      wm = sample_patterns[["151509"]]$wm,
      nb_wm = c(size = param$size, mu = param$mu_wm_151509),
      nb_outside = c(size = param$size, mu = param$mu_outside_151509),
      library_size = library_size_151509
    )
    counts_151669 <- simulate_counts(
      num_genes = numGenes,
      num_signal_genes = numSignalGenes,
      wm = sample_patterns[["151669"]]$wm,
      nb_wm = c(size = param$size, mu = param$mu_wm_151669),
      nb_outside = c(size = param$size, mu = param$mu_outside_151669),
      library_size = library_size_151669
    )

    pos_list <- list(sample_patterns[["151509"]]$pos, sample_patterns[["151669"]]$pos)
    counts_list <- list(counts_151509, counts_151669)
    counts_normalized_list <- list(
      normalize_counts(counts_151509, library_size_151509, method = "normalized"),
      normalize_counts(counts_151669, library_size_151669, method = "normalized")
    )

    p_values_INSiGHT_raw <- run_insight(counts_list, pos_list)
    p_values_INSiGHT_normalized <- run_insight(counts_normalized_list, pos_list)

    combined_p <- cbind(
      INSiGHT_raw = p_values_INSiGHT_raw,
      INSiGHT_normalized = p_values_INSiGHT_normalized
    )
    zero_prop_151509 <- rowMeans(counts_151509 == 0)
    zero_prop_151669 <- rowMeans(counts_151669 == 0)

    output_null <- data.frame(
      combined_p[(numSignalGenes + 1):numGenes, , drop = FALSE],
      size = rep(param$size, numGenes - numSignalGenes),
      iScenario = rep(param$iScenario, numGenes - numSignalGenes),
      scenario = rep(param$scenario, numGenes - numSignalGenes),
      scenario_label = rep(param$scenario_label, numGenes - numSignalGenes),
      mu_wm_151509 = rep(param$mu_wm_151509, numGenes - numSignalGenes),
      mu_outside_151509 = rep(param$mu_outside_151509, numGenes - numSignalGenes),
      mu_wm_151669 = rep(param$mu_wm_151669, numGenes - numSignalGenes),
      mu_outside_151669 = rep(param$mu_outside_151669, numGenes - numSignalGenes),
      n_spots_151509 = rep(ncol(counts_151509), numGenes - numSignalGenes),
      n_spots_151669 = rep(ncol(counts_151669), numGenes - numSignalGenes),
      library_cv_151509 = rep(sd(library_size_151509) / mean(library_size_151509), numGenes - numSignalGenes),
      library_cv_151669 = rep(sd(library_size_151669) / mean(library_size_151669), numGenes - numSignalGenes),
      zero_prop_151509 = zero_prop_151509[(numSignalGenes + 1):numGenes],
      zero_prop_151669 = zero_prop_151669[(numSignalGenes + 1):numGenes],
      iSim = rep(iSim, numGenes - numSignalGenes),
      check.names = FALSE
    )
    output_power <- data.frame(
      combined_p[seq_len(numSignalGenes), , drop = FALSE],
      size = rep(param$size, numSignalGenes),
      iScenario = rep(param$iScenario, numSignalGenes),
      scenario = rep(param$scenario, numSignalGenes),
      scenario_label = rep(param$scenario_label, numSignalGenes),
      mu_wm_151509 = rep(param$mu_wm_151509, numSignalGenes),
      mu_outside_151509 = rep(param$mu_outside_151509, numSignalGenes),
      mu_wm_151669 = rep(param$mu_wm_151669, numSignalGenes),
      mu_outside_151669 = rep(param$mu_outside_151669, numSignalGenes),
      n_spots_151509 = rep(ncol(counts_151509), numSignalGenes),
      n_spots_151669 = rep(ncol(counts_151669), numSignalGenes),
      library_cv_151509 = rep(sd(library_size_151509) / mean(library_size_151509), numSignalGenes),
      library_cv_151669 = rep(sd(library_size_151669) / mean(library_size_151669), numSignalGenes),
      zero_prop_151509 = zero_prop_151509[seq_len(numSignalGenes)],
      zero_prop_151669 = zero_prop_151669[seq_len(numSignalGenes)],
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
      file = out_file
    )
  }
}
