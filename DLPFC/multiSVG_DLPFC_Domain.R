library(Seurat)
library(stringr)
library(SPARK)
library(dplyr)
library(tibble)
library(ggvenn)
library(SpatialPCA)
library(patchwork)

setwd('/mnt/home/suhaohao/Research/INSiGHT/DLPFC')
data_path = '/mnt/ufs18/rs-037/Cui_ST_group/Datasets/DLPFC'
PROCESSED_PATH = '/mnt/home/suhaohao/Research/INSiGHT/DLPFC/'
source('/mnt/home/suhaohao/Research/INSiGHT/DLPFC/domainDetection.R')

args <- commandArgs(trailingOnly = TRUE)
features_file <- if (length(args) >= 1) args[[1]] else 'out/svgList_1.RDS'
method <- if (length(args) >= 2) args[[2]] else 'all'
strategy_arg <- if (length(args) >= 3) args[[3]] else 'all'
valid_strategy_args <- c('all', 'a', 'u', 'union', 'i', 'intersection', 'c', 'cauchy')
if (!tolower(strategy_arg) %in% valid_strategy_args && length(args) < 4) {
  output_suffix <- strategy_arg
  strategy_arg <- 'all'
} else {
  output_suffix <- if (length(args) >= 4) args[[4]] else paste(method, strategy_arg, sep = '_')
}
if (!grepl('^/', features_file)) {
  features_file <- file.path(PROCESSED_PATH, features_file)
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

normalize_method <- function(x) {
  tolower(gsub('[^A-Za-z0-9]+', '', x))
}
normalize_strategy <- function(x) {
  x <- tolower(x)
  if (x %in% c('all', 'a')) return('all')
  if (x %in% c('u', 'union')) return('U')
  if (x %in% c('i', 'intersection')) return('I')
  if (x %in% c('c', 'cauchy')) return('C')
  stop('Strategy must be one of: all, U, I, C, union, intersection, cauchy.')
}
available_methods <- unique(sub(' \\([UIC]\\)$', '', names(features_list)))
strategy <- normalize_strategy(strategy_arg)
if (tolower(method) == 'all') {
  if (strategy != 'all') stop('Method all requires strategy all.')
  # Run every named gene vector, including INSiGHT parameter variants.
} else {
  method_match <- available_methods[normalize_method(available_methods) == normalize_method(method)]
  if (length(method_match) != 1) {
    stop('Cannot match method "', method, '". Available methods: ',
         paste(available_methods, collapse = ', '), '. Use all to run every list.')
  }
  method <- method_match
  strategy_codes <- if (strategy == 'all') c('U', 'I', 'C') else strategy
  strategy_names <- paste0(method, ' (', strategy_codes, ')')
  if (all(strategy_names %in% names(features_list))) {
    features_list <- features_list[strategy_names]
  } else if (strategy == 'all' && method %in% names(features_list)) {
    features_list <- features_list[method]
  } else {
    stop('Cannot find requested SVG strategy/list for method ', method,
         '. Expected: ', paste(strategy_names, collapse = ', '))
  }
}
if (!length(features_list) || anyDuplicated(names(features_list))) {
  stop('SVG lists must have unique names and contain at least one gene set.')
}
valid_genes <- vapply(features_list, function(x) {
  is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(x))
}, logical(1))
if (!all(valid_genes)) {
  stop('Expected nonempty character gene vectors: ',
       paste(names(features_list)[!valid_genes], collapse = ', '))
}
features_list <- lapply(features_list, unique)
names(features_list) <- gsub('[^A-Za-z0-9]+', '_', names(features_list))
names(features_list) <- gsub('(^_+|_+$)', '', names(features_list))
output_suffix <- gsub('[^A-Za-z0-9_.-]+', '_', output_suffix)

message('Using SVG-list RDS: ', features_file)
message('Method: ', method)
message('Strategy: ', strategy_arg)
message('Gene sets: ', paste(names(features_list), collapse = ', '))

SAMPLES = c('151507','151508','151509','151510')

## load DLPFC and add ground truth
seu.obj_list <- list()
for(sample in SAMPLES){
  seu.st <- Load10X_Spatial(data.dir = file.path(data_path,sample),
                            filename="filtered_feature_bc_matrix.h5")
  meta_df <- read.csv(file.path(data_path,sample,'Human_Brain_Maynard_02082021_Visium_anno.csv'))
  meta_df$spot_name = lapply(meta_df$V1, function(x){
    tail(str_split(x,'_')[[1]],1)
  })
  meta_df = meta_df %>% dplyr::rename("ground_truth" = "V2") %>%
        column_to_rownames('spot_name') %>% select('ground_truth')

  seu.st <- AddMetaData(seu.st, metadata = meta_df)
  seu.st$Sample_ID = sample
  seu.st <- RenameCells(seu.st, new.names = paste0(sample,':',Cells(seu.st)))
  names(seu.st@images) <- sample
  seu.obj_list[[sample]] = seu.st
}
seu.obj <- Reduce(merge,seu.obj_list)
seu.obj <- JoinLayers(seu.obj)

## extract count & spatial
SAMPLES = names(seu.obj_list)
count_list <- lapply(SAMPLES,function(sample){
    seu.obj <- seu.obj_list[[sample]]
    seu.obj[['Spatial']]$counts
})
names(count_list) <- SAMPLES
location_list <- lapply(SAMPLES,function(sample){
    seu.obj <- seu.obj_list[[sample]]
    if(class(seu.obj@images[[1]]) == 'VisiumV2'){
      locations = GetTissueCoordinates(seu.obj)[,c('x','y')]
      colnames(locations) = c('row','col')
    }else{
      locations = seu.obj@images[[1]]@coordinates[,c('row','col')]
    }
    as.matrix(locations)
})
names(location_list) <- SAMPLES

str(features_list)
for (name in names(features_list)){
    set.seed(1234)
    seu.obj <- runSpatialPCA_Multi(seu.obj = seu.obj,count_list = count_list,location_list =location_list,
                               genes = features_list[[name]], names = name, clusternum = NULL)
}
saveRDS(seu.obj, file = file.path(PROCESSED_PATH,paste0(Sys.Date(),'_DLPFC_151507_151508_151509_151510_spatialPCA_', output_suffix, '.rds')))
