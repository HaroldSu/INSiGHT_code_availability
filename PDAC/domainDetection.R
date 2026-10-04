suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SpatialPCA)
})

# Run multi-sample SpatialPCA using a supplied SVG vector and add domains to the Seurat object.
runSpatialPCA_Multi <- function(seu.obj, count_list, location_list, genes, names,
                                clusternum = 3, numCores_spark = 8,
                                min.locations = 20, min.features = 3,
                                bandwidth_common = 0.1) {
  if (length(count_list) != length(location_list)) {
    stop("The number of count matrices must match the number of location matrices.")
  }
  if (length(clusternum) == 0 || any(!is.finite(clusternum)) || any(clusternum < 2)) {
    stop("clusternum must contain integers greater than 1.")
  }

  spa_multi_result <- SpatialPCA_Multiple_Sample(
    count_list = count_list,
    location_list = location_list,
    gene.type = "spatial",
    sparkversion = "spark",
    numCores_spark = numCores_spark,
    gene.number = length(genes),
    customGenelist = genes,
    min.loctions = min.locations,
    min.features = min.features,
    bandwidth_common = bandwidth_common
  )

  spatial_pcs <- Reduce(cbind, spa_multi_result[[2]])
  output_cells <- sub("^Sample[0-9]+_", "", colnames(spatial_pcs))
  if (!all(output_cells %in% colnames(seu.obj))) {
    stop("SpatialPCA output cells cannot be aligned to cells in the PDAC Seurat object.")
  }
  knearest <- round(sqrt(ncol(spatial_pcs)))

  for (cl in clusternum) {
    clusterlabel <- walktrap_clustering(
      clusternum = cl,
      latent_dat = spatial_pcs,
      knearest = knearest
    )
    metadata_df <- data.frame(x = clusterlabel, row.names = output_cells)
    colnames(metadata_df) <- paste0(names, "_spatialPCA_", cl)
    seu.obj <- AddMetaData(seu.obj, metadata = metadata_df)
  }

  seu.obj
}
