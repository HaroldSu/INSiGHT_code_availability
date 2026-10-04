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

dat_path <- "./DLPFC/data/DLPFC_151507_151508_151509_151510_combined_with_ground_truth.RDS"
sample_id <- "151507"
numSignalGenes <- 5000
numGenes <- 10000
nSim <- 1
target_num_spots <- 1000
mu_wm <- 0.25
mu_outside <- 1
nb_size <- 3
dropout_prob_list <- c(0, 0.1, 0.3, 0.5)
sim_params <- data.frame(
  mu_wm = mu_wm,
  mu_outside = mu_outside,
  size = nb_size,
  dropout_prob = dropout_prob_list
)
# Simulated gene pattern plot settings: display all dropout probabilities.
gene_to_display <- "Gene1"
gene_expression_scale <- "raw" # "raw" or "scaled"
sim_gene_plot_path <- "./Simulations/results/SimSupp1_151507_simulated_gene_pattern.jpeg"
gene_plot_axis_text_size <- 16
gene_plot_strip_text_size <- 18
gene_plot_legend_text_size <- 25
gene_plot_legend_title_size <- 25
gene_plot_title_size <- 30
gene_plot_subtitle_size <- 30
gene_plot_region_legend_point_size <- 8

simulate_counts <- function(num_genes, num_signal_genes, wm, nb_wm, nb_outside, dropout_prob = 0) {
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

  if (dropout_prob > 0) {
    nonzero_entries <- counts != 0
    dropout_entries <- matrix(
      stats::runif(length(counts)) < dropout_prob,
      nrow = nrow(counts),
      ncol = ncol(counts)
    )
    counts[nonzero_entries & dropout_entries] <- 0
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
  required_packages <- c(
    "scBSP",
    "SVG",
    "Seurat"
  )
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

plot_simulated_gene <- function(sim_design, gene, expression_scale = c("scaled", "raw"), out_path = NULL) {
  expression_scale <- match.arg(expression_scale)
  plot_df <- do.call(
    rbind,
    lapply(sim_design, function(design) {
      counts <- design$counts
      gene_name <- gene
      if (is.numeric(gene_name)) {
        gene_name <- rownames(counts)[gene_name]
      }
      if (!gene_name %in% rownames(counts)) {
        stop("Gene ", gene_name, " was not found in the simulated count matrix.")
      }

      gene_expression <- as.numeric(counts[gene_name, ])
      expression_range <- max(gene_expression) - min(gene_expression)
      plot_expression <- if (expression_scale == "raw") {
        gene_expression
      } else if (expression_range == 0) {
        rep(0, length(gene_expression))
      } else {
        gene_expression / expression_range
      }

      mu_wm <- design$nb_params$mu[design$nb_params$region == "WM"]
      mu_cortical <- design$nb_params$mu[design$nb_params$region == "cortical"]
      gamma_size <- design$nb_params$size[design$nb_params$region == "WM"]
      dropout_prob <- design$dropout_prob
      data.frame(
        x = -design$spatial.locs[, 1],
        y = design$spatial.locs[, 2],
        expression = plot_expression,
        domain = ifelse(design$domain, "WM", "cortical"),
        combo = paste0("mu[WM]==", mu_wm, "*','~mu[cortical]==", mu_cortical),
        gamma = paste0("gamma==", gamma_size),
        dropout = paste0("dropout==", dropout_prob),
        stringsAsFactors = FALSE
      )
    })
  )

  p <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(
      aes(fill = expression, color = domain),
      shape = 21,
      size = 1.5,
      stroke = 0.5
    ) +
    coord_fixed() +
    scale_fill_gradient(low = "white", high = "steelblue4") +
    scale_color_manual(
      values = c(WM = "firebrick3", cortical = "darkseagreen3")
    ) +
    facet_grid(
      combo + dropout ~ gamma,
      labeller = labeller(combo = label_parsed, gamma = label_parsed, dropout = label_parsed)
    ) +
    labs(
      title = "Example pattern: Supplementary Simulation 1",
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
      axis.text = element_text(size = gene_plot_axis_text_size),
      axis.title = element_text(size = gene_plot_axis_text_size),
      strip.text = element_text(size = gene_plot_strip_text_size),
      legend.text = element_text(size = gene_plot_legend_text_size),
      legend.title = element_text(size = gene_plot_legend_title_size, face = "bold"),
      plot.title = element_text(size = gene_plot_title_size),
      plot.subtitle = element_text(size = gene_plot_subtitle_size)
    )

  if (!is.null(out_path)) {
    ggsave(p, file = out_path, width = 10, height = 16, dpi = 300)
  }
  p
}

dat <- readRDS(dat_path)
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
pos <- grid_aggregation$pos
wm <- grid_aggregation$wm
cat(paste0(
  "Aggregated ", grid_aggregation$original_num_spots,
  " annotated spots into ", nrow(pos),
  " grid spots using target_num_spots = ", target_num_spots,
  ".\n"
))

nb_params <- do.call(
  base::rbind,
  lapply(seq_len(nrow(sim_params)), function(iParam) {
    mu_outside <- sim_params$mu_outside[iParam]
    mu_wm <- sim_params$mu_wm[iParam]
    nb_size <- sim_params$size[iParam]
    dropout_prob <- sim_params$dropout_prob[iParam]
    data.frame(
      sample_id = sample_id,
      region = c("WM", "cortical"),
      size = nb_size,
      mu = c(mu_wm, mu_outside),
      dropout_prob = dropout_prob,
      stringsAsFactors = FALSE
    )
  })
)
print(nb_params)

sim_design <- lapply(seq_len(nrow(sim_params)), function(iParam) {
  mu_outside <- sim_params$mu_outside[iParam]
  mu_wm <- sim_params$mu_wm[iParam]
  nb_size <- sim_params$size[iParam]
  dropout_prob <- sim_params$dropout_prob[iParam]
  nb_wm <- c(size = nb_size, mu = mu_wm)
  nb_outside <- c(size = nb_size, mu = mu_outside)
  counts <- simulate_counts(
    num_genes = numGenes,
    num_signal_genes = numSignalGenes,
    wm = wm,
    nb_wm = nb_wm,
    nb_outside = nb_outside,
    dropout_prob = dropout_prob
  )
  list(
    counts = counts,
    spatial.locs = pos,
    domain = wm,
    sample_id = sample_id,
    dropout_prob = dropout_prob,
    grid_aggregation = grid_aggregation,
    nb_params = data.frame(
      sample_id = sample_id,
      region = c("WM", "cortical"),
      size = nb_size,
      mu = c(mu_wm, mu_outside),
      dropout_prob = dropout_prob,
      stringsAsFactors = FALSE
    ),
    signal_genes = rownames(counts)[seq_len(numSignalGenes)],
    null_genes = rownames(counts)[(numSignalGenes + 1):numGenes]
  )
})
names(sim_design) <- paste0(
  "mu_wm_", sim_params$mu_wm,
  "_mu_nonwm_", sim_params$mu_outside,
  "_size_", sim_params$size,
  "_dropout_", sim_params$dropout_prob
)
# saveRDS(sim_design, file = "./Simulations/out/SimSupp1_DLPFC_WM_151507_design.RDS")

plot_simulated_gene(
  sim_design = sim_design,
  gene = gene_to_display,
  expression_scale = gene_expression_scale,
  out_path = sim_gene_plot_path
)

out_null <- data.frame(INSiGHT_rank_zero_FALSE = NULL,
                       INSiGHT_rank_zero_TRUE = NULL,
                       SPARK_X = NULL,
                       MERINGUE = NULL,
                       binSpect = NULL,
                       scBSP = NULL,
                       sample_id = NULL,
                       size = NULL,
                       mu_wm = NULL,
                       mu_outside = NULL,
                       dropout_prob = NULL,
                       zero_prop = NULL,
                       iSim = NULL)
out_power <- data.frame(INSiGHT_rank_zero_FALSE = NULL,
                        INSiGHT_rank_zero_TRUE = NULL,
                        SPARK_X = NULL,
                        MERINGUE = NULL,
                        binSpect = NULL,
                        scBSP = NULL,
                        sample_id = NULL,
                        size = NULL,
                        mu_wm = NULL,
                        mu_outside = NULL,
                        dropout_prob = NULL,
                        zero_prop = NULL,
                        iSim = NULL)

if (nSim > 0) {
  check_benchmark_packages()
}

for (iParam in seq_len(nrow(sim_params))) {
  mu_outside <- sim_params$mu_outside[iParam]
  mu_wm <- sim_params$mu_wm[iParam]
  nb_size <- sim_params$size[iParam]
  dropout_prob <- sim_params$dropout_prob[iParam]
  nb_wm <- c(size = nb_size, mu = mu_wm)
  nb_outside <- c(size = nb_size, mu = mu_outside)

  for (iSim in seq_len(nSim)) {
    cat(paste0(
      "Working on mu_wm = ", mu_wm,
      ", mu_outside = ", mu_outside,
      ", size = ", nb_size,
      ", dropout_prob = ", dropout_prob,
      ", simulation ", iSim, ".\n"
    ))

    counts <- simulate_counts(
      num_genes = numGenes,
      num_signal_genes = numSignalGenes,
      wm = wm,
      nb_wm = nb_wm,
      nb_outside = nb_outside,
      dropout_prob = dropout_prob
    )
    zero_prop <- rowMeans(counts == 0)

    object_INSiGHT_rank_zero_FALSE <- createINSiGHTobject(data.input = list(counts), spatial.locs = list(pos))
    object_INSiGHT_rank_zero_FALSE <- processINSiGHT(object = object_INSiGHT_rank_zero_FALSE,
                                                     k.NN4bandwidth = 36, quantile.bandwidth = 0.3,
                                                     kernel.option = "single", rank.tie = "average",
                                                     rank.zero = FALSE,
                                                     NN4Kernel = FALSE, k.NN4Kernel = 100,
                                                     normalization = "none"
                                                     )
    object_INSiGHT_rank_zero_FALSE <- testINSiGHT(object = object_INSiGHT_rank_zero_FALSE)
    p_values_INSiGHT_rank_zero_FALSE <- object_INSiGHT_rank_zero_FALSE@testResults$p.values_combined
    rm(object_INSiGHT_rank_zero_FALSE)

    object_INSiGHT_rank_zero_TRUE <- createINSiGHTobject(data.input = list(counts), spatial.locs = list(pos))
    object_INSiGHT_rank_zero_TRUE <- processINSiGHT(object = object_INSiGHT_rank_zero_TRUE,
                                                    k.NN4bandwidth = 36, quantile.bandwidth = 0.3,
                                                    kernel.option = "single", rank.tie = "average",
                                                    rank.zero = TRUE,
                                                    NN4Kernel = FALSE, k.NN4Kernel = 100,
                                                    normalization = "none"
                                                    )
    object_INSiGHT_rank_zero_TRUE <- testINSiGHT(object = object_INSiGHT_rank_zero_TRUE)
    p_values_INSiGHT_rank_zero_TRUE <- object_INSiGHT_rank_zero_TRUE@testResults$p.values_combined
    rm(object_INSiGHT_rank_zero_TRUE)

    object_sparkx <- sparkx(counts, pos, option = "single", verbose = F)
    p_values_sparkx <- object_sparkx$res_mtest$combinedPval
    rm(object_sparkx)

    p_values_meringue <- run_svg_package_method(counts, pos, method = "MERINGUE")
    p_values_binspect <- run_svg_package_method(counts, pos, method = "binSpect")

    object_scBSP <- scBSP::scBSP(Coords = pos, ExpMat_Sp = as(counts, "dgCMatrix"))
    p_values_scBSP <- object_scBSP$P_values
    rm(object_scBSP)

    output_null <- data.frame(INSiGHT_rank_zero_FALSE = p_values_INSiGHT_rank_zero_FALSE[(numSignalGenes + 1):numGenes],
                              INSiGHT_rank_zero_TRUE = p_values_INSiGHT_rank_zero_TRUE[(numSignalGenes + 1):numGenes],
                              SPARK_X = p_values_sparkx[(numSignalGenes + 1):numGenes],
                              MERINGUE = p_values_meringue[(numSignalGenes + 1):numGenes],
                              binSpect = p_values_binspect[(numSignalGenes + 1):numGenes],
                              scBSP = p_values_scBSP[(numSignalGenes + 1):numGenes],
                              sample_id = rep(sample_id, numGenes - numSignalGenes),
                              size = rep(nb_size, numGenes - numSignalGenes),
                              mu_wm = rep(mu_wm, numGenes - numSignalGenes),
                              mu_outside = rep(mu_outside, numGenes - numSignalGenes),
                              dropout_prob = rep(dropout_prob, numGenes - numSignalGenes),
                              zero_prop = zero_prop[(numSignalGenes + 1):numGenes],
                              iSim = rep(iSim, numGenes - numSignalGenes))
    output_power <- data.frame(INSiGHT_rank_zero_FALSE = p_values_INSiGHT_rank_zero_FALSE[seq_len(numSignalGenes)],
                               INSiGHT_rank_zero_TRUE = p_values_INSiGHT_rank_zero_TRUE[seq_len(numSignalGenes)],
                               SPARK_X = p_values_sparkx[seq_len(numSignalGenes)],
                               MERINGUE = p_values_meringue[seq_len(numSignalGenes)],
                               binSpect = p_values_binspect[seq_len(numSignalGenes)],
                               scBSP = p_values_scBSP[seq_len(numSignalGenes)],
                               sample_id = rep(sample_id, numSignalGenes),
                               size = rep(nb_size, numSignalGenes),
                               mu_wm = rep(mu_wm, numSignalGenes),
                               mu_outside = rep(mu_outside, numSignalGenes),
                               dropout_prob = rep(dropout_prob, numSignalGenes),
                               zero_prop = zero_prop[seq_len(numSignalGenes)],
                               iSim = rep(iSim, numSignalGenes))

    out_null <- base::rbind(out_null, output_null)
    out_power <- base::rbind(out_power, output_power)
    save(
      out_null,
      out_power,
      nb_params,
      grid_aggregation,
      file = "./Simulations/out/SimSupp1Out_DLPFC_151507.RData"
    )
  }
}
