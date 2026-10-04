library(ggplot2)
library(Seurat)
library(SPARK)
library(stringr)
library(mclust)
library(SpatialPCA)

paperColors =c("gold", "skyblue", "orchid4", "#868603", "royalblue2","darkorange",  "violet",'cyan','#6c0404',
              "#90bca6","#bbca81","#f1f9c7","#efc8c2","#6EA0B8","#5c92b6","#e2c987")

SpatialPCA_Multiple_Sample_tune = function(count_list,location_list,gene.type="spatial",sparkversion="spark",numCores_spark=5,gene.number=3000, customGenelist=NULL,min.loctions = 20, min.features=20,bandwidth_common=0.1,SpatialPCnum = 20){

    require(Seurat)
    if(length(count_list) != length(location_list)) {
            stop("The number of count matrix should match with the number of location matrix. ")
    }# end fi

    # create seurat object for each dataset
    seurat_list = list()
    for(count in 1:length(count_list)){
        colnames(count_list[[count]]) = paste0("Sample",count,"_",colnames(count_list[[count]]))
        rownames(location_list[[count]]) = colnames(count_list[[count]])
        seurat_list[[count]] <- CreateSeuratObject(counts = count_list[[count]], project = "multiple")
    }

    # normalize and identify variable features for each dataset independently
    seurat_list <- lapply(X = seurat_list, FUN = function(x) {
        print(paste0("Normalizing dataset"))
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 10000)
    })

    # select features that are repeatedly variable across datasets for integration
    features <- Seurat::SelectIntegrationFeatures(object.list = seurat_list,nfeatures = 10000)
    # integration, find anchors
    MultipleSample.anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
    # this command creates an 'integrated' data assay
    MultipleSample.combined <- Seurat::IntegrateData(anchorset = MultipleSample.anchors)
    DefaultAssay(MultipleSample.combined) <- "integrated"
    # Run the standard workflow for visualization and clustering
    MultipleSample.combined <- Seurat::ScaleData(MultipleSample.combined, verbose = FALSE)
    # obtain integrated and normalized data
    integrated_data = MultipleSample.combined@assays$integrated@scale.data
    # create each SpatialPCA object, this step is essentially for spatial gene selection
    spatialpca_list = list()
    for(count in 1:length(seurat_list)){
        print(paste0("Creating SpatialPCA object for dataset ",count))
        spatialpca_list[[count]] <- CreateSpatialPCAObject(counts=count_list[[count]], location=location_list[[count]], project = "SpatialPCA",gene.type=gene.type,sparkversion=sparkversion,numCores_spark=numCores_spark,gene.number=gene.number, customGenelist=customGenelist,min.loctions = min.loctions, min.features=min.features)
    }

    matched_data = list()
    matched_location = list()
    Kernal_mat = list()
    for(count in 1:length(count_list)){
        print(paste0("Creating SpatialPCA kernel for dataset ",count))
        match_gene = na.omit(match(rownames( spatialpca_list[[count]]@normalized_expr), rownames(integrated_data)))
        match_spot = grep(paste0("1_",count),colnames(integrated_data))
        matched_data[[count]] = integrated_data[match_gene,match_spot]
        matched_location[[count]] = spatialpca_list[[count]]@location
        spatialpca_list[[count]]@normalized_expr = t(scale(t(matched_data[[count]])))
        spatialpca_list[[count]]@location = scale(matched_location[[count]])
        spatialpca_list[[count]] = SpatialPCA_buildKernel(spatialpca_list[[count]], kerneltype="gaussian", bandwidth.set.by.user=bandwidth_common)
        Kernal_mat[[count]] = spatialpca_list[[count]]@kernelmat
    }

    common_genes = unique(unlist(lapply(matched_data, function(x){ row.names(x)})))
    MultipleSample_merge = spatialpca_list[[1]] # initialize using the first data, will replace with integrated data later
    MultipleSample_merge@normalized_expr = t(scale(t(integrated_data[match(common_genes, rownames(integrated_data)),])))
    MultipleSample_merge@kernelmat = Matrix::bdiag(Kernal_mat)
    MultipleSample_merge@kernelmat = as.matrix(MultipleSample_merge@kernelmat)
    # MultipleSample_merge@kernelmat = as(MultipleSample_merge@kernelmat, "sparseMatrix")
    # MultipleSample_merge@sparseKernel = TRUE
    MultipleSample_merge@params$expr = MultipleSample_merge@normalized_expr
    MultipleSample_merge@location = do.call("rbind",lapply(spatialpca_list, function(x){ x@location}))

    # perform SpatialPCA on integrated expression and integrated block diagonal kernel
    MultipleSample_merge = SpatialPCA_EstimateLoading(MultipleSample_merge,fast=TRUE,SpatialPCnum=SpatialPCnum)
    MultipleSample_merge = SpatialPCA_SpatialPCs(MultipleSample_merge, fast=TRUE)
    colnames(MultipleSample_merge@SpatialPCs) = rownames(MultipleSample_merge@location)
    
    id_list = list()
    SpatialPC_list = list()
    Location_spatialpc_list = list()
    for(sample in 1:length(count_list)){
	    id_list[[sample]] = grep(paste0("Sample",sample), colnames(MultipleSample_merge@SpatialPCs))
	    SpatialPC_list[[sample]] = MultipleSample_merge@SpatialPCs[,id_list[[sample]]]
            Location_spatialpc_list[[sample]] = spatialpca_list[[sample]]@location
    }

    names(SpatialPC_list) = paste0("Sample",1:length(SpatialPC_list))
    names(Location_spatialpc_list) = paste0("Sample",1:length(Location_spatialpc_list))
    return(list("MultipleSample_SpatialPCA_object"=MultipleSample_merge, 
		"SpatialPC_list" = SpatialPC_list, 
		"Location_spatialpc_list"=Location_spatialpc_list))
}

#' @title Run SPARK SVGs
#' @param seu.obj Spatial Seurat object
#' @param gene_pct_thres  QC threshold, to filter out genes express less than gene_pct_thres percent of spots
#' @param spot_num_thres  QC threshold, to filter out spots with less than spot_num_thres number of genes
#' @return return a vector of SVGs 
#' @import SPARK
run_SPARK = function(seu.obj,gene_pct_thres = 0.1, spot_num_thres = 10, genes=NULL,spots = NULL ){
    counts = LayerData(seu.obj, assays = 'Spatial', layer = 'counts')
    if(class(seu.obj@images[[1]]) == 'VisiumV2'){
      print('VisiumV2')
      locations = GetTissueCoordinates(seu.obj)[,c('x','y')]
      colnames(locations) = c('row','col')
      # locations = Reduce(rbind, lapply(names(seu.obj@images),function(sample){
      #     coord_mtx = seu.obj@images[[sample]]@ boundaries$centroids@coords                        
      #     coord_mtx 
      # }))
      # colnames(locations) = c('row','col')
      # print(dim(locations))
    }else{
      print('VisiumV1')
      locations = seu.obj@images[[1]]@coordinates[,c('row','col')]
    }
    

    if (!is.null(genes)){
      counts = counts[na.omit(match(genes,rownames(counts))),]
      # counts = counts[genes,]
      gene_pct_thres = 0 
    }
    if (!is.null(spots)){
      counts = counts[,na.omit(match(spots,colnames(counts)))]
      locations = locations[na.omit(match(spots,row.names(locations))),]
      # counts = counts[,spots]
      # locations = locations[spots, ]
      spot_num_thres = 0
    }
    spark <- CreateSPARKObject(counts=counts,
                            location=locations,
                            percentage = gene_pct_thres , ## to filter gene 
                            min_total_counts = spot_num_thres) ## to filter spots
    spark@lib_size <- apply(spark@counts, 2, sum)
    spark <- spark.vc(spark,
                    covariates = NULL,
                    lib_size = spark@lib_size,
                    fit.model = 'gaussian',
                    num_core = 16,
                    verbose = F)
    spark <- spark.test(spark,
                        check_positive = T,
                        verbose = F)
    Spark_SVG_df = spark@res_mtest %>% arrange(combined_pvalue) 
    Spark_SVG_df
    # SVGs_spark = rownames(Spark_SVG_df ) [Spark_SVG_df $adjusted_pvalue < 0.05]
    # SVGs_spark
}

## For SeuratPCA
#' @title Run Seurat PCA domain detection use selected features
#' @param seu.obj Spatial Seurat object
#' @param genes  Genes used as input features
#' @param names Name of the setting, used as domain name 
#' @return return a seurat object with SeuratPCA domains, named as  paste0(names,'_seuratPCA)
runSeuratPCA = function(seu.obj, genes, names, clusternum = NULL,npcs = 30, verbose = FALSE,algorithm = 1){

    if(class(seu.obj@images[[1]]) == 'VisiumV2'){
      locations = GetTissueCoordinates(seu.st)[,c('x','y')]
      colnames(locations) = c('row','col')
    }else{
      locations = seu.obj@images[[1]]@coordinates[,c('row','col')]
    }

    if (is.null(clusternum)){
      if(!'ground_truth' %in% colnames(seu.obj@meta.data)) {
          stop('please provide ground truth label in ground_truth column')
          }
      clusternum = length(unique(seu.obj$ground_truth[!is.na(seu.obj$ground_truth)]))
      domain_name = paste0(names,'_seuratPCA')
    }
    
    if (names == 'hvg2000'){
        seu.obj <- SCTransform(seu.obj, assay = "Spatial", verbose =  verbose, variable.features.n = 2000)
    }else{ ## mkgenes
        ## return.only.var.genes: so scale.data return all genes, for cluster DE heatmap 
        seu.obj <- SCTransform(seu.obj, assay = "Spatial", verbose =  verbose, residual.features = NULL, return.only.var.genes = FALSE ) 
    }
    for(cl in clusternum){
      domain_name = paste0(names,'_seuratPCA_', cl)
      seu.obj <- RunPCA(seu.obj, assay = "SCT", verbose = FALSE, features = intersect(genes , rownames(seu.obj)),npcs = npcs)
      seu.obj <- FindNeighbors(seu.obj,dims = 1:npcs,verbose = verbose)
      best_res <- res_search(seu.obj, target_k = cl, algorithm = algorithm ,verbose = verbose)
      seu.obj <- FindClusters(seu.obj,resolution = best_res, algorithm = algorithm, verbose = verbose)
      seu.obj@meta.data[,domain_name]  = seu.obj$seurat_clusters
      # seu.obj@meta.data[,paste0(domain_name,'_refined')] = refine_cluster_10x(clusterlabels=seu.obj@meta.data[,domain_name],
      #               location=as.matrix(GetTissueCoordinates(seu.obj)[,c('x','y')]),shape="hexagon")
      seu.obj@meta.data[,paste0(domain_name,'_refined')] = refine_cluster_10x(clusterlabels=seu.obj@meta.data[,domain_name],
                    location=as.matrix(locations),shape="hexagon")
    }
    seu.obj
    
}
#' @title Run SpatialPCA using user-defined genes as input features
#' @param seu.obj Spatial Seurat object
#' @param genes  Genes used as input features
#' @param names Name of the setting, used as domain name 
#' @return return a seurat object with spatialPCA domains, named as  paste0(names,'_spatialPCA')
#' @import SpatialPCA
runSpatialPCA = function(seu.obj, genes,names, clusternum = NULL, min.loctions = 20, min.features = 20,PCnum = 20 ){
  if (is.null(clusternum)){
    if(!'ground_truth' %in% colnames(seu.obj@meta.data)) {
        stop('please provide ground truth label in ground_truth column')
        }
    clusternum = length(unique(seu.obj$ground_truth[!is.na(seu.obj$ground_truth)]))
  }
  ## extract location
  if(class(seu.obj@images[[1]]) == 'VisiumV2'){
      locations = GetTissueCoordinates(seu.st)[,c('x','y')]
      colnames(locations) = c('row','col')
    }else{
      locations = seu.obj@images[[1]]@coordinates[,c('row','col')]
  }

  spca.obj = CreateSpatialPCAObject(counts = LayerData(seu.obj, assays = 'Spatial', layer = 'counts'), 
               location= as.matrix(locations), 
               customGenelist=genes,min.loctions = min.loctions, min.features = min.features)
  spca.obj = SpatialPCA_buildKernel(spca.obj, kerneltype="gaussian", bandwidthtype="SJ",bandwidth.set.by.user=NULL)
  spca.obj = SpatialPCA_EstimateLoading(spca.obj,fast=FALSE,SpatialPCnum = PCnum ) 
  spca.obj = SpatialPCA_SpatialPCs(spca.obj, fast=FALSE)
  knearest = round(sqrt(dim(spca.obj@location)[1]))
  for (cl in clusternum){
    domain_name = paste0(names,'_spatialPCA_', cl)
    clusterlabel = walktrap_clustering(clusternum=cl,latent_dat=spca.obj@SpatialPCs,knearest=knearest) 
    clusterlabel_refine = refine_cluster_10x(clusterlabels=clusterlabel,location=spca.obj@location,shape="hexagon")
    metadata_df = data.frame(x = clusterlabel_refine,row.names = colnames(spca.obj@SpatialPCs) )
    colnames(metadata_df) = c(domain_name ) 
    # metadata_df = data.frame(x = clusterlabel_refine, y = clusterlabel,row.names = colnames(spca.obj@SpatialPCs) )
    # colnames(metadata_df) = c(domain_name , paste0(domain_name,'_origin') ) 
    seu.obj = AddMetaData(seu.obj, metadata = metadata_df)
  }
  seu.obj
}

#' @title Run SpatialPCA_Multiple_Sample using user-defined genes as input features
#' @param seu.obj Spatial Seurat object
#' @param count_list A list of count matrix of all samples
#' @param location_list A list of location matrix of all samples
#' @param genes  Genes used as input features
#' @param names Name of the setting, used as domain name 
#' @return return a seurat object with spatialPCA domains, named as  paste0(names,'_spatialPCA_',cl)
#' @import SpatialPCA
runSpatialPCA_Multi <- function(seu.obj,count_list,location_list, genes, names, clusternum){
    if (is.null(clusternum)){
      if(!'ground_truth' %in% colnames(seu.obj@meta.data)) {
          stop('please provide ground truth label in ground_truth column')
          }
      clusternum = length(unique(seu.obj$ground_truth[!is.na(seu.obj$ground_truth)]))
      domain_name = paste0(names,'_spatialPCA')
    }

    spa_multi_result <- SpatialPCA_Multiple_Sample(count_list = count_list  ,location_list = location_list,gene.type="spatial",
                                      sparkversion="spark",numCores_spark=5,gene.number=3000, 
                                      customGenelist= genes,min.loctions = 20, min.features=3,bandwidth_common=0.1)
    print('here')                                  
    SpatialPCs <- Reduce(cbind, spa_multi_result[[2]])
    knearest = round(sqrt(ncol(SpatialPCs)))
    for (cl in clusternum){
        clusterlabel = walktrap_clustering(clusternum=cl,latent_dat=SpatialPCs,knearest=knearest) 
        metadata_df = data.frame(x = clusterlabel,row.names = colnames(seu.obj) )
        domain_name <- paste0(names,'_spatialPCA_', cl)
        colnames(metadata_df) = c(domain_name ) 
        seu.obj = AddMetaData(seu.obj, metadata = metadata_df)
    }
    seu.obj
}


#' @title Run SpaceFlow domain detection use selected features
#' @param seu.obj Spatial Seurat object
#' @param genes  Genes used as input features
#' @param names Name of the setting, used as domain name 
#' @return return a seurat object with SpaceFlow domains, named as paste0(names,'_SpaceFlow)
runSpaceFlow = function(seu.obj,  genes, names, clusternum = NULL){ 
  if (is.null(clusternum)){
    if(!'ground_truth' %in% colnames(seu.obj@meta.data)) {
        stop('please provide ground truth label in ground_truth column')
        }
    clusternum = length(unique(seu.obj$ground_truth[!is.na(seu.obj$ground_truth)]))
  }
  domain_result = SpaceFlow_py(count_mtx = t(as.matrix(seu.obj[['Spatial']]$count)), genes= rownames(seu.obj[['Spatial']]$count),
            spots = colnames(seu.obj[['Spatial']]$count), loc_df = seu.obj@images[[1]]@ coordinates[,c('imagecol','imagerow')],
              SVGs = genes, target_cl = clusternum )
  seu.obj@meta.data[,paste0(names, '_', colnames(domain_result))] = domain_result #'SpaceFlow_'
  seu.obj
}


#' @title Run StLearn domain detection use selected features
#' @param seu.obj Spatial Seurat object
#' @param genes  Genes used as input features
#' @param names Name of the setting, used as domain name 
#' @return return a seurat object with SpaceFlow domains, named as paste0(names,'_StLearn)
runStLearn = function(seu.obj, adata_path,  genes, names, clusternum = NULL){ 
    if(is.null(clusternum)){
      clusternum = length(unique(seu.obj$ground_truth[!is.na(seu.obj$ground_truth)]))
    }
    domain_result = StLearn_py(adata_path = adata_path, genes= rownames(seu.obj[['Spatial']]$count),
             spots = colnames(seu.obj[['Spatial']]$count),
                SVGs = genes, target_cl = clusternum )
    # domain_result
    seu.obj@meta.data[,paste0(names,'_StLearn')] = domain_result
    seu.obj
}









### For BayesSpace
#' @title Run BayesSpace spatialPreprocess() function, but replace hvg by SVGs
#' @param SVGs A vector of SVGs, used as input features to get PCA embedding
spatialPreprocess_svg <- function(sce, platform=c("Visium", "ST"),
                              n.PCs=15, SVGs, skip.PCA=FALSE,
                              log.normalize=TRUE, assay.type="logcounts") {
    
    ## Set BayesSpace metadata
    metadata(sce)$BayesSpace.data <- list()
    metadata(sce)$BayesSpace.data$platform <- match.arg(platform)
    metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
    # metadata(sce)$BayesSpace.data$use_dimred <- use.dimred
    # metadata(sce)$BayesSpace.data$d <- n.PCs

    ## Run PCA on HVGs, log-normalizing if necessary
    if (!skip.PCA) {
        if (log.normalize)
            sce <- logNormCounts(sce)
   
        dec <- modelGeneVar(sce, assay.type=assay.type)
        # top <- getTopHVGs(dec, n=n.HVGs)
        sce <- runPCA(sce, subset_row=SVGs, ncomponents=n.PCs, exprs_values=assay.type)
        rowData(sce)[["is.HVG"]] <- (rownames(sce) %in% SVGs)
    }

    sce
}



#' @title Run BayesSpace domain detection use selected features
#' @param seu.obj Spatial Seurat object
#' @param genes  Genes used as input features
#' @param names Name of the setting, used as domain name 
#' @return return a seurat object with BayesSpace domains, named as paste0(names,'_BayesSpace)
runBayesSpace = function(seu.obj,  genes, names, clusternum = NULL){ 
    if(is.null(clusternum)){
      clusternum = length(unique(seu.obj$ground_truth[!is.na(seu.obj$ground_truth)]))
    }
    counts = LayerData(seu.obj, assays = 'Spatial', layer = 'counts')
    sce <- SingleCellExperiment(assays=list(counts=counts),colData = seu.obj@images[[1]]@ coordinates[,c('col','row')])

    set.seed(102)
    ## the default setting of BayesSpace selects is.Hvg by HVGs, but here I replace them by SVGs 
    sce <- spatialPreprocess_svg(sce, platform="ST", 
                                n.PCs=30, SVGs = genes, log.normalize = TRUE)
    sce <- spatialCluster(sce, q=clusternum, platform="ST", d=30,
                            init.method="mclust", model="t", gamma=2,
                            nrep=5000, burn.in=100,
                            save.chain=TRUE)
    seu.obj@meta.data[,paste0(names,'_BayesSpace')] = colData(sce)$spatial.cluster
    seu.obj
}


# runBayesSpace = function(seu.obj,  genes, names, clusternum = NULL){ 
#     if(is.null(clusternum)){
#       clusternum = length(unique(seu.obj$ground_truth[!is.na(seu.obj$ground_truth)]))
#     }
#     counts = LayerData(seu.obj, assays = 'Spatial', layer = 'counts')
#     sce <- SingleCellExperiment(assays=list(counts=counts),colData = seu.obj@images[[1]]@ coordinates[,c('col','row')])

#     set.seed(102)
#     ## the default setting of BayesSpace selects is.Hvg by HVGs, but here I replace them by SVGs 
#     sce <- spatialPreprocess_svg(sce, platform="ST", 
#                                 n.PCs=50, SVGs = genes, log.normalize = TRUE)
#     sce <- spatialCluster(sce, q=clusternum, platform="ST", d=30,
#                             init.method="mclust", model="t", gamma=2,
#                             nrep=1000, burn.in=100,
#                             save.chain=TRUE)
#     seu.obj@meta.data[,paste0(names,'_BayesSpace')] = colData(sce)$spatial.cluster
#     seu.obj
# }






res_search <- function(seu.obj, target_k = 7, algorithm = 4, res_start = 0.5, res_step = 0.05, res_epochs = 20, random_state = 0, verbose = TRUE) {
  if (verbose) {
    cat("searching resolution to k =", target_k, "\n")
  }
  res <- res_start
  seu.obj <- suppressWarnings(FindClusters(seu.obj, verbose = FALSE, resolution = res, random.seed = random_state, algorithm = algorithm))
  old_k <- length(unique(seu.obj$seurat_clusters))
  if (verbose) {
    cat(sprintf("Res = %f, Num of clusters = %d\n", res, old_k))
  }

  run <- 0
  while (old_k != target_k) {
    old_sign <- ifelse(old_k < target_k, 1, -1)
    seu.obj <- suppressWarnings(FindClusters(seu.obj, resolution = res + res_step * old_sign, random.seed = random_state, verbose = FALSE, algorithm = algorithm))
    new_k <- length(unique(seu.obj$seurat_clusters))

    if (verbose) {
      cat(sprintf("Res = %f, Num of clusters = %d\n", res + res_step * old_sign, new_k))
    }

    if (new_k == target_k) {
      res <- res + res_step * old_sign
      if (verbose) {
        cat(sprintf("recommended res = %f\n", res))
      }
      break
    }

    new_sign <- ifelse(new_k < target_k, 1, -1)
    if (new_sign == old_sign) {
      res <- res + res_step * old_sign
      if (verbose) {
        cat(sprintf("Res changed to %f\n", res))
      }
      old_k <- new_k
    } else {
      res_step <- res_step / 2
      if (verbose) {
        cat(sprintf("Res step halved to %f\n", res_step))
      }
    }

    if (run > res_epochs) {
      if (verbose) {
        cat("Exact resolution not found\n")
        cat(sprintf("Recommended res = %f\n", res))
      }
      break
    }

    run <- run + 1
  }

  if (verbose) {
    cat(sprintf("Final recommended res = %f\n", res))
  }
  
  return(res)
}
