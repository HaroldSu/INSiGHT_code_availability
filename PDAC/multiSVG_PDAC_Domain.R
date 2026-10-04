suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SpatialPCA)
})

setwd('/mnt/home/suhaohao/Research/INSiGHT/PDAC')
DATA_PATH <- '/mnt/home/suhaohao/Research/INSiGHT/PDAC/data'
PROCESSED_PATH <- '/mnt/home/suhaohao/Research/INSiGHT/PDAC/out'
source('/mnt/home/suhaohao/Research/INSiGHT/PDAC/domainDetection.R')

args <- commandArgs(trailingOnly = TRUE)
features_file <- if (length(args) >= 1) args[[1]] else file.path(PROCESSED_PATH, 'svgList_INSiGHT_3.RDS')
feature_set <- if (length(args) >= 2) args[[2]] else 'INSiGHT'
clusternum <- if (length(args) >= 3) as.integer(strsplit(args[[3]], ',', fixed = TRUE)[[1]]) else 8L
output_suffix <- if (length(args) >= 4) args[[4]] else 'INSiGHT_3sample_k8'
bandwidth_common <- if (length(args) >= 5) as.numeric(args[[5]]) else 0.1
min_features <- if (length(args) >= 6) as.integer(args[[6]]) else 20L
num_cores_spark <- if (length(args) >= 7) as.integer(args[[7]]) else 10L

if (length(clusternum) == 0 || anyNA(clusternum) || any(clusternum < 2)) {
  stop('Region numbers must be integers greater than 1.')
}
if (!is.finite(bandwidth_common) || bandwidth_common <= 0) {
  stop('bandwidth_common must be a positive numeric value.')
}
if (is.na(min_features) || min_features < 1 || is.na(num_cores_spark) || num_cores_spark < 1) {
  stop('min_features and num_cores_spark must be positive integers.')
}

if (!grepl('^/', features_file)) {
  features_file <- file.path(getwd(), features_file)
}
if (!file.exists(features_file)) {
  stop('Cannot find SVG-list RDS: ', features_file)
}
features_list <- readRDS(features_file)
if (!is.list(features_list)) {
  features_list <- list(input_genes = as.character(features_list))
}
if (is.null(names(features_list)) || any(names(features_list) == '')) {
  names(features_list) <- paste0('gene_set_', seq_along(features_list))
}
if (feature_set != 'all') {
  if (!feature_set %in% names(features_list)) {
    stop('Cannot find gene set ', feature_set, ' in ', features_file)
  }
  features_list <- features_list[feature_set]
}

SAMPLES <- c('S1A_NAA16', 'S1B_NAA12', 'S1C_NAA11')
message('Using SVG-list RDS: ', features_file)
message('Gene sets: ', paste(names(features_list), collapse = ', '))
message('Samples: ', paste(SAMPLES, collapse = ', '))
message('Cluster number(s): ', paste(clusternum, collapse = ', '))
message('bandwidth_common: ', bandwidth_common, '; min.features: ', min_features,
        '; numCores_spark: ', num_cores_spark)

seu.obj <- readRDS(file.path(DATA_PATH, 'ST8.rds'))
missing_images <- setdiff(SAMPLES, names(seu.obj@images))
if (length(missing_images) > 0) {
  stop('Missing images in ST8.rds: ', paste(missing_images, collapse = ', '))
}

sample_cells <- unlist(lapply(SAMPLES, function(sample) {
  rownames(seu.obj@images[[sample]]@coordinates)
}), use.names = FALSE)
seu.obj <- subset(seu.obj, cells = sample_cells)
seu.obj@images <- seu.obj@images[SAMPLES]

counts_all <- GetAssayData(seu.obj, assay = 'Spatial', layer = 'counts')
count_list <- list()
location_list <- list()
for (sample in SAMPLES) {
  locations <- seu.obj@images[[sample]]@coordinates[, c('row', 'col'), drop = FALSE]
  cells <- intersect(rownames(locations), colnames(counts_all))
  if (length(cells) == 0) {
    stop('No expression cells were found for image ', sample)
  }
  count_list[[sample]] <- counts_all[, cells, drop = FALSE]
  location_list[[sample]] <- as.matrix(locations[cells, , drop = FALSE])
}

for (name in names(features_list)) {
  genes <- intersect(as.character(features_list[[name]]), rownames(counts_all))
  if (length(genes) == 0) {
    stop('No genes for ', name, ' are present in the Spatial assay.')
  }
  message('Running SpatialPCA for ', name, ' with ', length(genes), ' genes.')
  seu.obj <- runSpatialPCA_Multi(
    seu.obj = seu.obj,
    count_list = count_list,
    location_list = location_list,
    genes = genes,
    names = name,
    clusternum = clusternum,
    numCores_spark = num_cores_spark,
    min.features = min_features,
    bandwidth_common = bandwidth_common
  )
}

out_file <- file.path(
  PROCESSED_PATH,
  paste0(Sys.Date(), '_PDAC_S1A_S1B_S1C_spatialPCA_', output_suffix, '.rds')
)
saveRDS(seu.obj, out_file)
message('Saved ', out_file)
