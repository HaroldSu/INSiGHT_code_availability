rm(list = ls())

library(methods)
library(Matrix)
library(RANN)
library(gaston)
library(SPARK)
library(matrixStats)
library(INSiGHT)

set.seed(111)

sample_id <- "synthetic"
numGenes <- 10000
numSignalGenes <- 1000
nSim <- 5
target_num_spots_list <- c(1000, 2000, 3000, 5000)
mu_wm <- 0.25
mu_outside <- 1
nb_size <- 3
out_file <- "./Simulations/out/SimSupp4Out_synthetic_runtime.RData"

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

generate_spatial_design <- function(num_spots) {
  nx <- ceiling(sqrt(num_spots))
  ny <- ceiling(num_spots / nx)
  pos <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = nx),
    y = seq(0, 1, length.out = ny)
  ))
  pos <- pos[seq_len(num_spots), , drop = FALSE]
  colnames(pos) <- c("x", "y")
  rownames(pos) <- paste0("Spot", seq_len(num_spots))

  wm_boundary <- 0.35 + 0.08 * sin(2 * pi * pos[, 1])
  wm <- pos[, 2] <= wm_boundary
  names(wm) <- rownames(pos)

  list(
    pos = pos,
    wm = wm,
    original_num_spots = num_spots,
    target_num_spots = num_spots,
    aggregated_num_spots = num_spots,
    nx = nx,
    ny = ny
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

elapsed_time <- function(expr) {
  start_time <- proc.time()
  value <- force(expr)
  elapsed <- as.numeric((proc.time() - start_time)[["elapsed"]])
  list(value = value, elapsed = elapsed)
}

timed_method <- function(method, counts, pos) {
  status <- "completed"
  error_message <- NA_character_
  timing <- c(create = NA_real_, preprocess = NA_real_, test = NA_real_, total = NA_real_)

  result <- tryCatch({
    gc()
    if (method %in% c("INSiGHT", "INSiGHT_fast")) {
      total_start <- proc.time()
      use_fast_eigen <- identical(method, "INSiGHT_fast")

      create_result <- elapsed_time(
        createINSiGHTobject(data.input = list(counts), spatial.locs = list(pos))
      )
      object <- create_result$value

      process_args <- list(
        object = object,
        k.NN4bandwidth = 36,
        quantile.bandwidth = 0.3,
        kernel.option = "single",
        NN4Kernel = TRUE,
        k.NN4Kernel = 100,
        eigen.less = use_fast_eigen,
        rank.zero = F,
        rank.tie = "average",
        normalization = "none"
      )
      if (use_fast_eigen) process_args$eigen.prop <- 0.1

      process_result <- elapsed_time(
        do.call(processINSiGHT, process_args)
      )
      object <- process_result$value

      test_result <- elapsed_time(testINSiGHT(object = object))
      object <- test_result$value
      p_values <- object@testResults$p.values_combined

      timing[["create"]] <- create_result$elapsed
      timing[["preprocess"]] <- process_result$elapsed
      timing[["test"]] <- test_result$elapsed
      timing[["total"]] <- as.numeric((proc.time() - total_start)[["elapsed"]])
      rm(object)
    } else if (method == "SPARK_X") {
      total_result <- elapsed_time(sparkx(counts, pos, option = "mixture", verbose = FALSE))
      p_values <- total_result$value$res_mtest$combinedPval
      timing[["total"]] <- total_result$elapsed
      rm(total_result)
    } else if (method == "MERINGUE") {
      total_result <- elapsed_time(run_svg_package_method(counts, pos, method = "MERINGUE"))
      p_values <- total_result$value
      timing[["total"]] <- total_result$elapsed
      rm(total_result)
    } else if (method == "binSpect") {
      total_result <- elapsed_time(run_svg_package_method(counts, pos, method = "binSpect"))
      p_values <- total_result$value
      timing[["total"]] <- total_result$elapsed
      rm(total_result)
    } else if (method == "scBSP") {
      total_result <- elapsed_time(scBSP::scBSP(Coords = pos, ExpMat_Sp = as(counts, "dgCMatrix")))
      p_values <- total_result$value$P_values
      timing[["total"]] <- total_result$elapsed
      rm(total_result)
    } else {
      stop("Unknown method: ", method)
    }

    list(timing = timing, num_p_values = length(p_values))
  }, error = function(e) {
    list(
      timing = timing,
      num_p_values = NA_integer_,
      status = "error",
      error_message = conditionMessage(e)
    )
  })

  if (!is.null(result$status)) {
    status <- result$status
    error_message <- result$error_message
  }

  data.frame(
    method = method,
    create_time = result$timing[["create"]],
    preprocess_time = result$timing[["preprocess"]],
    test_time = result$timing[["test"]],
    total_time = result$timing[["total"]],
    num_p_values = result$num_p_values,
    status = status,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

check_benchmark_packages()

nb_wm <- c(size = nb_size, mu = mu_wm)
nb_outside <- c(size = nb_size, mu = mu_outside)
nb_params <- data.frame(
  sample_id = sample_id,
  region = c("WM", "cortical"),
  size = nb_size,
  mu = c(mu_wm, mu_outside),
  stringsAsFactors = FALSE
)

methods_to_run <- c(
  "INSiGHT", "INSiGHT_fast", "SPARK_X", "MERINGUE", "binSpect", "scBSP"
)
runtime_results <- data.frame()
grid_aggregations <- list()

for (target_num_spots in target_num_spots_list) {
  grid_aggregation <- generate_spatial_design(target_num_spots)
  pos <- grid_aggregation$pos
  wm <- grid_aggregation$wm
  grid_aggregations[[as.character(target_num_spots)]] <- grid_aggregation

  cat(paste0(
    "Target spots = ", target_num_spots,
    "; actual spots = ", nrow(pos), ".\n"
  ))

  for (iSim in seq_len(nSim)) {
    cat(paste0(
      "Running single-sample runtime benchmark: target spots = ", target_num_spots,
      ", actual spots = ", nrow(pos),
      ", simulation ", iSim, ".\n"
    ))

    counts <- simulate_counts(
      num_genes = numGenes,
      num_signal_genes = numSignalGenes,
      wm = wm,
      nb_wm = nb_wm,
      nb_outside = nb_outside
    )

    for (method in methods_to_run) {
      cat(paste0("  Method: ", method, "\n"))
      method_result <- timed_method(method, counts, pos)
      method_result$sample_id <- sample_id
      method_result$target_num_spots <- target_num_spots
      method_result$num_spots <- nrow(pos)
      method_result$num_genes <- numGenes
      method_result$num_signal_genes <- numSignalGenes
      method_result$size <- nb_size
      method_result$mu_wm <- mu_wm
      method_result$mu_outside <- mu_outside
      method_result$iSim <- iSim

      runtime_results <- base::rbind(runtime_results, method_result)
      save(
        runtime_results,
        nb_params,
        grid_aggregations,
        target_num_spots_list,
        methods_to_run,
        file = out_file
      )
    }
  }
}
