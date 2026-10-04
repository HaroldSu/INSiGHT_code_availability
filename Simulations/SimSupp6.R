rm(list = ls())

library(methods)
library(Matrix)
library(RANN)
library(gaston)
library(INSiGHT)

set.seed(111)

dat_path <- "./DLPFC/data/DLPFC_151507_151508_151509_151510_combined_with_ground_truth.RDS"
sample_id <- "151507"
numSignalGenes <- 5000
numGenes <- 10000
nSim <- 1
target_num_spots <- 1000
mu_wm_list <- c(0.25, 8)
mu_outside_list <- c(1, 4)
nb_size_list <- c(1.5, 3)
k_s_list <- c(6, 18, 36, 100)
out_file <- "./Simulations/out/SimSupp6Out_DLPFC_151507_k_s_sensitivity.RData"

sim_params <- merge(
  data.frame(
    iCombo = seq_along(mu_wm_list),
    mu_wm = mu_wm_list,
    mu_outside = mu_outside_list
  ),
  data.frame(size = nb_size_list)
)
sim_params <- sim_params[order(sim_params$iCombo, sim_params$size), ]

method_columns <- c(
  paste0("INSiGHT_NN4Kernel_TRUE_k_", k_s_list),
  "INSiGHT_NN4Kernel_FALSE"
)

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

assign_grid_cells <- function(pos, nx, ny) {
  x_bin <- pmin(floor(pos[, 1] * nx) + 1, nx)
  y_bin <- pmin(floor(pos[, 2] * ny) + 1, ny)
  paste(x_bin, y_bin, sep = "_")
}

aggregate_spots_by_grid <- function(pos, wm, target_num_spots) {
  if (target_num_spots >= nrow(pos)) {
    rownames(pos) <- paste0("Spot", seq_len(nrow(pos)))
    names(wm) <- rownames(pos)
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

run_insight <- function(counts, pos, NN4Kernel, k_s = NULL) {
  object <- createINSiGHTobject(data.input = list(counts), spatial.locs = list(pos))
  object <- processINSiGHT(
    object = object,
    k.NN4bandwidth = 36,
    quantile.bandwidth = 0.3,
    kernel.option = "single",
    rank.tie = "average", rank.zero = FALSE,
    NN4Kernel = NN4Kernel,
    k.NN4Kernel = ifelse(is.null(k_s), 100, k_s),
    normalization = "none"
  )
  object <- testINSiGHT(object = object)
  object@testResults$p.values_combined
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
    data.frame(
      sample_id = sample_id,
      region = c("WM", "cortical"),
      size = nb_size,
      mu = c(mu_wm, mu_outside),
      stringsAsFactors = FALSE
    )
  })
)
print(nb_params)

out_null <- data.frame(
  matrix(numeric(0), ncol = length(method_columns), dimnames = list(NULL, method_columns)),
  sample_id = character(),
  size = numeric(),
  mu_wm = numeric(),
  mu_outside = numeric(),
  iSim = integer(),
  check.names = FALSE
)
out_power <- out_null

for (iParam in seq_len(nrow(sim_params))) {
  mu_outside <- sim_params$mu_outside[iParam]
  mu_wm <- sim_params$mu_wm[iParam]
  nb_size <- sim_params$size[iParam]
  nb_wm <- c(size = nb_size, mu = mu_wm)
  nb_outside <- c(size = nb_size, mu = mu_outside)

  for (iSim in seq_len(nSim)) {
    cat(paste0(
      "Working on mu_wm = ", mu_wm,
      ", mu_outside = ", mu_outside,
      ", size = ", nb_size,
      ", simulation ", iSim, ".\n"
    ))

    counts <- simulate_counts(
      num_genes = numGenes,
      num_signal_genes = numSignalGenes,
      wm = wm,
      nb_wm = nb_wm,
      nb_outside = nb_outside
    )

    p_values_by_method <- lapply(k_s_list, function(k_s) {
      cat(paste0("  Running INSiGHT with NN4Kernel = TRUE, k_s = ", k_s, ".\n"))
      run_insight(
        counts = counts,
        pos = pos,
        NN4Kernel = TRUE,
        k_s = k_s
      )
    })
    names(p_values_by_method) <- paste0("INSiGHT_NN4Kernel_TRUE_k_", k_s_list)

    cat("  Running INSiGHT with NN4Kernel = FALSE.\n")
    p_values_by_method[["INSiGHT_NN4Kernel_FALSE"]] <- run_insight(
      counts = counts,
      pos = pos,
      NN4Kernel = FALSE
    )

    combined_p <- as.data.frame(p_values_by_method, check.names = FALSE)

    output_null <- data.frame(
      combined_p[(numSignalGenes + 1):numGenes, , drop = FALSE],
      sample_id = rep(sample_id, numGenes - numSignalGenes),
      size = rep(nb_size, numGenes - numSignalGenes),
      mu_wm = rep(mu_wm, numGenes - numSignalGenes),
      mu_outside = rep(mu_outside, numGenes - numSignalGenes),
      iSim = rep(iSim, numGenes - numSignalGenes),
      check.names = FALSE
    )
    output_power <- data.frame(
      combined_p[seq_len(numSignalGenes), , drop = FALSE],
      sample_id = rep(sample_id, numSignalGenes),
      size = rep(nb_size, numSignalGenes),
      mu_wm = rep(mu_wm, numSignalGenes),
      mu_outside = rep(mu_outside, numSignalGenes),
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
      k_s_list,
      method_columns,
      file = out_file
    )
  }
}
