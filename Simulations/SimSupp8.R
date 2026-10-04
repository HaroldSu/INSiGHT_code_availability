# Cross-section dependence sensitivity analysis. Run from the workspace root.
# No counts are shared across genes or spatial locations. Only the three
# sections at a matching gene/location are coupled by the Gaussian copula.

# ---- User settings ----
settings <- list(
  num_samples = 3L,
  num_spots = 1000L,
  num_genes = 10000L,
  num_signal_genes = 5000L,
  n_replicates = 5L,
  rho_values = c(0, 0.2, 0.5, 0.8),
  mu_wm = 0.25,
  mu_cortical = 1,
  mu_null = 1,
  nb_size = 1.5,
  base_seed = 111L,
  genes_per_chunk = 250L,
  uniform_epsilon = 1e-10,
  # Same domain in all sections; use the curved WM boundary from SimSupp5.
  wm_boundary_height = 0.35,
  wm_boundary_amplitude = 0.08,
  process_args = list(
    spot.threshold = 10,
    gene.threshold = 0.05,
    normalization = "none",
    rank.zero = FALSE,
    rank.tie = "average",
    kernel.option = "single",
    NN4Kernel = TRUE,
    k.NN4Kernel = 100,
    # The supplementary manuscript specifies k_b = 100 and quantile = 0.1.
    # Several older simulation scripts instead use 36 and 0.3.
    k.NN4bandwidth = 36,
    quantile.bandwidth = 0.3,
    eigen.less = FALSE,
    eigen.prop = 0.1
  ),
  test_args = list(QuadForm.method = "davies", p.adjust.method = "BY")
)
result_file <- "./Simulations/out/SimSupp8_cross_section_dependence.rds"
# An existing checkpoint is resumed automatically. Choose another filename
# to start a separate analysis. Increasing n_replicates permits extension.
# -----------------------

validate_settings <- function(cfg) {
  integer_fields <- c("num_samples", "num_spots", "num_genes",
                      "num_signal_genes", "n_replicates", "base_seed", "genes_per_chunk")
  for (field in integer_fields) {
    value <- cfg[[field]]
    if (length(value) != 1L || !is.finite(value) || value < 1 || value != floor(value)) {
      stop(field, " must be a positive integer.")
    }
  }
  if (cfg$num_samples != 3L) stop("This experiment and its summaries require three samples.")
  if (cfg$num_signal_genes >= cfg$num_genes) stop("Both null genes and SVGs are required.")
  if (!length(cfg$rho_values) || any(!is.finite(cfg$rho_values)) ||
      any(cfg$rho_values < 0 | cfg$rho_values >= 1) || anyDuplicated(cfg$rho_values)) {
    stop("rho_values must contain unique values in [0, 1).")
  }
  if (cfg$n_replicates >= 100000 ||
      cfg$base_seed + length(cfg$rho_values) * 100000 >= .Machine$integer.max) {
    stop("Replicate count or base_seed exceeds the supported seed range.")
  }
  for (field in c("mu_wm", "mu_cortical", "mu_null", "nb_size")) {
    if (length(cfg[[field]]) != 1L || !is.finite(cfg[[field]]) || cfg[[field]] <= 0) {
      stop(field, " must be a positive finite number.")
    }
  }
  if (length(cfg$uniform_epsilon) != 1L || !is.finite(cfg$uniform_epsilon) ||
      cfg$uniform_epsilon <= 0 || cfg$uniform_epsilon >= 0.5) {
    stop("uniform_epsilon must lie in (0, 0.5).")
  }
  if (cfg$process_args$k.NN4bandwidth >= cfg$num_spots ||
      (isTRUE(cfg$process_args$NN4Kernel) && cfg$process_args$k.NN4Kernel >= cfg$num_spots)) {
    stop("Neighbor counts must be smaller than num_spots.")
  }
}

construct_spatial_domain <- function(cfg) {
  nx <- ceiling(sqrt(cfg$num_spots))
  ny <- ceiling(cfg$num_spots / nx)
  pos <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = nx), y = seq(0, 1, length.out = ny)
  ))[seq_len(cfg$num_spots), , drop = FALSE]
  rownames(pos) <- paste0("Spot", seq_len(nrow(pos)))
  boundary <- cfg$wm_boundary_height + cfg$wm_boundary_amplitude * sin(2 * pi * pos[, 1])
  wm <- pos[, 2] <= boundary
  if (!any(wm) || all(wm)) stop("The spatial design must contain WM and cortical spots.")
  list(pos = pos, wm = wm, nx = nx, ny = ny)
}

generate_correlated_nb <- function(cfg, design, rho) {
  sample_names <- paste0("section", seq_len(cfg$num_samples))
  gene_names <- paste0("Gene", seq_len(cfg$num_genes))
  counts <- setNames(lapply(sample_names, function(x) {
    matrix(0L, cfg$num_genes, nrow(design$pos),
           dimnames = list(gene_names, rownames(design$pos)))
  }), sample_names)
  sigma <- rho^abs(outer(seq_len(cfg$num_samples), seq_len(cfg$num_samples), "-"))
  chol_sigma <- chol(sigma)

  # Chunks avoid allocating a genes x spots x sections latent array at once.
  for (first in seq.int(1L, cfg$num_genes, by = cfg$genes_per_chunk)) {
    idx <- seq.int(first, min(first + cfg$genes_per_chunk - 1L, cfg$num_genes))
    mu <- matrix(cfg$mu_null, length(idx), nrow(design$pos))
    signal <- idx <= cfg$num_signal_genes
    mu[signal, design$wm] <- cfg$mu_wm
    mu[signal, !design$wm] <- cfg$mu_cortical
    latent <- matrix(rnorm(length(mu) * cfg$num_samples), ncol = cfg$num_samples) %*% chol_sigma
    for (m in seq_len(cfg$num_samples)) {
      uniforms <- pmin(1 - cfg$uniform_epsilon,
                       pmax(cfg$uniform_epsilon, pnorm(latent[, m])))
      values <- qnbinom(uniforms, mu = as.vector(mu), size = cfg$nb_size)
      if (any(!is.finite(values)) || any(values > .Machine$integer.max)) {
        stop("Invalid or oversized NB counts; check the generation parameters.")
      }
      counts[[m]][idx, ] <- matrix(as.integer(values), nrow = length(idx))
    }
  }
  counts
}

simulate_dataset <- function(cfg, design, rho, seed) {
  set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  counts <- generate_correlated_nb(cfg, design, rho)
  list(
    counts = counts,
    positions = setNames(rep(list(design$pos), cfg$num_samples), names(counts)),
    is_svg = setNames(seq_len(cfg$num_genes) <= cfg$num_signal_genes, rownames(counts[[1]]))
  )
}

measure_correlations <- function(dataset, chunk_size) {
  pairs <- list(c(1L, 2L), c(2L, 3L), c(1L, 3L))
  names(pairs) <- c("cor12", "cor23", "cor13")
  values <- matrix(NA_real_, length(dataset$is_svg), length(pairs),
                   dimnames = list(names(dataset$is_svg), names(pairs)))
  for (first in seq.int(1L, nrow(values), by = chunk_size)) {
    idx <- seq.int(first, min(first + chunk_size - 1L, nrow(values)))
    ranks <- lapply(dataset$counts, function(x) {
      r <- matrixStats::rowRanks(x[idx, , drop = FALSE],
                                ties.method = "average", preserveShape = TRUE)
      sweep(r, 1, rowMeans(r), "-")
    })
    for (j in seq_along(pairs)) {
      a <- ranks[[pairs[[j]][1]]]
      b <- ranks[[pairs[[j]][2]]]
      denominator <- sqrt(rowSums(a^2) * rowSums(b^2))
      values[idx, j] <- ifelse(denominator > 0, rowSums(a * b) / denominator, NA_real_)
    }
  }
  # Null genes isolate stochastic dependence: SVGs share a spatial mean pattern
  # even at rho = 0. Average gene-wise Spearman correlations, not pooled counts.
  groups <- list(null = !dataset$is_svg, svg = dataset$is_svg, all = rep(TRUE, nrow(values)))
  output <- list()
  for (group in names(groups)) {
    for (pair in names(pairs)) {
      v <- values[groups[[group]], pair]
      valid <- is.finite(v)
      suffix <- if (group == "null") "" else paste0("_", group)
      output[[paste0(pair, suffix)]] <- if (any(valid)) mean(v[valid]) else NA_real_
      output[[paste0("n_", pair, suffix)]] <- sum(valid)
    }
  }
  output
}

run_insight <- function(dataset, cfg) {
  object <- INSiGHT::createINSiGHTobject(
    data.input = dataset$counts, spatial.locs = dataset$positions,
    sample.names = names(dataset$counts)
  )
  object <- do.call(INSiGHT::processINSiGHT, c(list(object = object), cfg$process_args))
  expected_genes <- names(dataset$is_svg)
  if (!setequal(object@genes, expected_genes) ||
      !all(vapply(object@spaCoord, nrow, integer(1)) == cfg$num_spots)) {
    stop("INSiGHT QC removed genes or spots. Check thresholds before changing assessment denominators.")
  }
  object <- do.call(INSiGHT::testINSiGHT, c(list(object = object), cfg$test_args))
  result <- object@testResults
  if (!all(c("ID", "p.values_combined") %in% names(result)) ||
      anyDuplicated(result$ID) || !setequal(as.character(result$ID), expected_genes)) {
    stop("INSiGHT testResults do not contain uniquely identified results for every simulated gene.")
  }
  p <- result$p.values_combined[match(expected_genes, as.character(result$ID))]
  if (any(!is.finite(p)) || any(p < 0 | p > 1)) {
    stop("INSiGHT returned invalid joint p-values; inspect the package output before proceeding.")
  }
  setNames(p, expected_genes)
}

summarize_results <- function(p_values, dataset, correlations, rho, replicate, seed) {
  null_p <- p_values[!dataset$is_svg]
  svg_p <- p_values[dataset$is_svg]
  metrics <- data.frame(rho = rho, replicate = replicate, seed = seed)
  for (threshold in c(0.01, 0.05, 0.10)) {
    suffix <- sprintf("%03d", as.integer(round(100 * threshold)))
    metrics[[paste0("type1_", suffix)]] <- mean(null_p < threshold)
    metrics[[paste0("power_", suffix)]] <- mean(svg_p < threshold)
  }
  for (name in names(correlations)) metrics[[name]] <- correlations[[name]]
  metrics$cor_adjacent <- mean(c(metrics$cor12, metrics$cor23))
  metrics$n_null <- length(null_p)
  metrics$n_svg <- length(svg_p)
  for (m in seq_along(dataset$counts)) {
    x <- dataset$counts[[m]][!dataset$is_svg, , drop = FALSE]
    metrics[[paste0("null_mean_section", m)]] <- mean(x)
    metrics[[paste0("null_zero_prop_section", m)]] <- mean(x == 0)
  }
  metrics
}

save_checkpoint <- function(state, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".SimSupp8_", tmpdir = dirname(path), fileext = ".rds")
  saveRDS(state, temporary)
  # Same-directory rename atomically replaces the checkpoint on macOS/Linux.
  if (!file.rename(temporary, path)) {
    stop("Cannot replace checkpoint. Previous checkpoint retained; new results are in ", temporary)
  }
}

implementation_signature <- function() {
  function_names <- c("createINSiGHTobject", "processINSiGHT", "testINSiGHT",
                      "constructGaussianW", "QualityControl", "ComputeBandwidth", "CenterKernel", "HWH_matvec")
  package_names <- c("INSiGHT", "Matrix", "matrixStats", "RANN", "RSpectra", "CompQuadForm")
  list(
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    packages = setNames(vapply(package_names, function(p) {
      as.character(utils::packageVersion(p))
    }, character(1)), package_names),
    functions = setNames(lapply(function_names, function(name) {
      f <- get(name, envir = asNamespace("INSiGHT"), inherits = FALSE)
      list(formals = formals(f), body = deparse(body(f)))
    }), function_names)
  )
}

run_simulation <- function(cfg, path) {
  validate_settings(cfg)
  for (package in c("INSiGHT", "matrixStats")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Required package is missing: ", package)
  }
  signature <- implementation_signature()
  design <- construct_spatial_domain(cfg)
  if (file.exists(path)) {
    state <- readRDS(path)
    if (!identical(state$schema_version, 1L)) stop("Unrecognized SimSupp8 checkpoint format.")
    old_config <- state$settings
    new_config <- cfg
    old_config$n_replicates <- new_config$n_replicates <- NULL
    if (!identical(old_config, new_config) || !identical(state$implementation, signature)) {
      stop("Checkpoint settings or INSiGHT implementation differ. Choose a new result_file for this analysis.")
    }
    if (nrow(state$replicate_results) && max(state$replicate_results$replicate) > cfg$n_replicates) {
      stop("n_replicates is smaller than the existing completed replicate numbers.")
    }
    state$settings <- cfg
    message("Resuming ", length(state$runs), " completed replicate(s).")
  } else {
    state <- list(
      schema_version = 1L, settings = cfg, implementation = signature,
      design = design, runs = list(), replicate_results = data.frame(),
      failures = data.frame(), session_info = utils::sessionInfo(),
      correlation_definition = paste(
        "cor12/cor23/cor13: mean gene-wise Spearman correlation across matched",
        "locations for null genes. _svg and _all fields are additional diagnostics.")
    )
    save_checkpoint(state, path)
  }

  for (i_rho in seq_along(cfg$rho_values)) {
    rho <- cfg$rho_values[[i_rho]]
    for (i_rep in seq_len(cfg$n_replicates)) {
      key <- sprintf("rho%02d_rep%03d", i_rho, i_rep)
      if (!is.null(state$runs[[key]])) next
      seed <- as.integer(cfg$base_seed + (i_rho - 1L) * 100000L + i_rep)
      message("SimSupp8: rho = ", rho, ", replicate ", i_rep, "/", cfg$n_replicates)
      started <- proc.time()[["elapsed"]]
      warnings_seen <- character()
      one_run <- tryCatch(withCallingHandlers({
        dataset <- simulate_dataset(cfg, design, rho, seed)
        correlations <- measure_correlations(dataset, cfg$genes_per_chunk)
        p <- run_insight(dataset, cfg)
        metrics <- summarize_results(p, dataset, correlations, rho, i_rep, seed)
        metrics$elapsed_seconds <- unname(proc.time()[["elapsed"]] - started)
        list(metrics = metrics, null_p_values = p[!dataset$is_svg])
      }, warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
      }), error = function(e) e)

      if (inherits(one_run, "error")) {
        state$failures <- rbind(state$failures, data.frame(
          rho = rho, replicate = i_rep, seed = seed,
          error = conditionMessage(one_run), stringsAsFactors = FALSE
        ))
        save_checkpoint(state, path)
        stop("Replicate failed; completed results are saved. Resume retries this replicate.\n",
             conditionMessage(one_run), call. = FALSE)
      }
      one_run$warnings <- unique(warnings_seen)
      state$runs[[key]] <- one_run
      state$replicate_results <- do.call(rbind, lapply(state$runs, `[[`, "metrics"))
      rownames(state$replicate_results) <- NULL
      save_checkpoint(state, path)
      message("Saved: Type I error (0.05) = ", round(one_run$metrics$type1_005, 4),
              "; power (0.05) = ", round(one_run$metrics$power_005, 4))
      rm(dataset, correlations, p, metrics, one_run)
      invisible(gc())
    }
  }
  message("Simulation complete: ", path)
  invisible(state)
}

simulation_results <- run_simulation(settings, result_file)
