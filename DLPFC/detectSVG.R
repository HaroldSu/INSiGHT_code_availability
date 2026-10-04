rm(list = ls())

# library(Seurat)
library(dplyr)
library(tidyr)
library(methods)
library(SPARK)
library(scBSP)
library(INSiGHT)
library(SVG)
library(ACAT)

dat <- readRDS("./DLPFC/data/DLPFC_151507_151508_151509_151510_combined.RDS")

# ==============================================================================
# INSiGHT
## -----------------------------------------------------------------------------
object_INSiGHT <- createINSiGHTobject(data.input = dat$counts,
                                      spatial.locs = dat$spatial.locs,
                                      sample.names = dat$samples.id)
object_INSiGHT <- processINSiGHT(object = object_INSiGHT, spot.threshold = 10, gene.threshold = 0.01,
                                 NN4Kernel = TRUE, k.NN4Kernel = 100,
                                 rank.tie = "average",rank.zero = FALSE,
                                 normalization = "library-size")
object_INSiGHT <- testINSiGHT(object = object_INSiGHT, QuadForm.method = "davies",
                              p.adjust.method = "BY")

p.values.adj <- object_INSiGHT@testResults$p.values_adjusted
SVG.INSiGHT <- object_INSiGHT@genes[p.values.adj < 0.05]
## -----------------------------------------------------------------------------
## Output the SVG list for INSiGHT
## -----------------------------------------------------------------------------
SVG.INSiGHT.list <- list()
SVG.INSiGHT.list[[1]] <- SVG.INSiGHT
names(SVG.INSiGHT.list) <- "SVG.INSiGHT"
saveRDS(SVG.INSiGHT.list, file = "./DLPFC/out/geneList/svgList_INSiGHT.RDS")
saveRDS(object = object_INSiGHT@genes, file = "./DLPFC/out/geneList_QC.RDS")
# ==============================================================================
# SPARK-X
## -----------------------------------------------------------------------------
genes.afterQC <- readRDS("./DLPFC/out/geneList_QC.RDS")
n.genes.afterQC <- length(genes.afterQC)
p.values.allSample <- list()

for(iSample in 1:4){
  counts <- dat$counts[[iSample]][genes.afterQC,]
  pos <- dat$spatial.locs[[iSample]]
  object_sparkx <- sparkx(counts, pos,
                          option = "mixture")
  p.values.allSample[[iSample]] <- object_sparkx$res_mtest$combinedPval
  rm(counts, pos)
}

p.values.Cauchy <- apply(matrix(unlist(p.values.allSample), nrow = 4, byrow = T), 2, ACAT::ACAT)
p.values.Cauchy.adj <- p.adjust(p.values.Cauchy, method = "BY")

p.values.allSample.adj <- p.adjust(unlist(p.values.allSample), method = "BY")

SVG_SPARKX.allSample <- list()
for(iSample in 1:4){
  p.values.adj.perSample <- p.values.allSample.adj[((iSample - 1) * n.genes.afterQC + 1) : (iSample * n.genes.afterQC)]
  SVG_SPARKX.allSample[[iSample]] <- genes.afterQC[p.values.adj.perSample < 0.05]
}

SVG.SPARKX_union <- c()
SVG.SPARKX_intersection <- genes.afterQC
for(iSample in 1:4){
  SVG.SPARKX_union <- c(SVG.SPARKX_union, SVG_SPARKX.allSample[[iSample]])
  SVG.SPARKX_intersection <- intersect(SVG.SPARKX_intersection, SVG_SPARKX.allSample[[iSample]])
}
SVG.SPARKX_union <- unique(SVG.SPARKX_union)

SVG.SPARKX_cauchy <- genes.afterQC[p.values.Cauchy.adj < 0.05]
## -----------------------------------------------------------------------------
## Output the SVG list for SPARK-X
## -----------------------------------------------------------------------------
SVG.SPARKX.list <- list()
SVG.SPARKX.list[[1]] <- SVG.SPARKX_union
SVG.SPARKX.list[[2]] <- SVG.SPARKX_intersection
SVG.SPARKX.list[[3]] <- SVG.SPARKX_cauchy
names(SVG.SPARKX.list) <- c("SPARKX_union",
                               "SPARKX_intersection",
                               "SPARKX_Cauchy")
saveRDS(SVG.SPARKX.list, file = "./DLPFC/out/svgList_SPARKX.RDS")
# ==============================================================================
# scBSP
## -----------------------------------------------------------------------------
genes.afterQC <- readRDS("./DLPFC/out/geneList_QC.RDS")
n.genes.afterQC <- length(genes.afterQC)
p.values.allSample <- list()

for(iSample in 1:4){
  counts <- dat$counts[[iSample]][genes.afterQC,]
  pos <- dat$spatial.locs[[iSample]]
  object_scBSP <- scscBSP(Coords = pos, ExpMat_Sp = as(counts, "dgCMatrix"))
  p.values.allSample[[iSample]] <- object_scBSP$P_values
}

p.values.Cauchy <- apply(matrix(unlist(p.values.allSample), nrow = 4, byrow = T), 2, ACAT::ACAT)
p.values.Cauchy.adj <- p.adjust(p.values.Cauchy, method = "BY")
p.values.allSample.adj <- p.adjust(unlist(p.values.allSample), method = "BY")

SVG_scBSP.allSample <- list()
for(iSample in 1:4){
  p.values.adj.perSample <- p.values.allSample.adj[((iSample - 1) * n.genes.afterQC + 1) : (iSample * n.genes.afterQC)]
  SVG_scBSP.allSample[[iSample]] <- genes.afterQC[p.values.adj.perSample < 0.05]
}

SVG.scBSP_union <- c()
SVG.scBSP_intersection <- genes.afterQC
for(iSample in 1:4){
  SVG.scBSP_union <- c(SVG.scBSP_union, SVG_scBSP.allSample[[iSample]])
  SVG.scBSP_intersection <- intersect(SVG.scBSP_intersection, SVG_scBSP.allSample[[iSample]])
}
SVG.scBSP_union <- unique(SVG.scBSP_union)
SVG.scBSP_cauchy <- genes.afterQC[p.values.Cauchy.adj < 0.05]
## -----------------------------------------------------------------------------
## Output the SVG list for scBSP
## -----------------------------------------------------------------------------
SVG.scBSP.list <- list()
SVG.scBSP.list[[1]] <- SVG.scBSP_union
SVG.scBSP.list[[2]] <- SVG.scBSP_intersection
SVG.scBSP.list[[3]] <- SVG.scBSP_cauchy
names(SVG.scBSP.list) <- c("scBSP_union",
                         "scBSP_intersection",
                         "scBSP_Cauchy")
saveRDS(SVG.scBSP.list, file = "./DLPFC/out/svgList_scBSP.RDS")
# ==============================================================================
# Moran's I (MERINGUE)
## -----------------------------------------------------------------------------
genes.afterQC <- readRDS("./DLPFC/out/geneList_QC.RDS")
n.genes.afterQC <- length(genes.afterQC)
p.values.allSample <- list()

for(iSample in 1:4){
  cat(paste0("Working on sample ", dat$samples.id[iSample], "...\n"))
  # 1. Extract counts and spatial coordinates
  counts <- dat$counts[[iSample]][genes.afterQC,]
  pos <- as.matrix(dat$spatial.locs[[iSample]])

  # 2. Normalize the expression matrix (MERINGUE expects normalized data)
  # Standard log-normalization equivalent to Seurat's LogNormalize
  expr_norm <- normalizeCounts(counts)

  # 3. Compute Moran's I utilizing the SVG package's C++ acceleration
  # The 'knn' network method naturally implements a k-NN truncated adjacency matrix
  res_meringue <- CalSVG_MERINGUE(
    expr_matrix = as.matrix(expr_norm),
    spatial_coords = pos,
    network_method = "delaunay",
    k = 5,                 # Adjust to match your intended environmental bandwidth
    n_threads = 1,          # Multithreading for speed
    adjust_method = "BY"  # P-value adjustment is done downstream with BY
  )

  idx <- match(rownames(counts), res_meringue$gene)
  # 4. Extract unadjusted p-values (ensure the exact column name matches your SVG version, e.g., 'p_value' or 'pval')
  p_vals <- res_meringue$p.value[idx]

  # 5. Handle potential NAs (e.g., genes with 0 spatial variance) to prevent ACAT failure
  p_vals[is.na(p_vals)] <- 1.0

  p.values.allSample[[iSample]] <- p_vals

  # Clean up environment to manage memory across loops
  rm(counts, pos, lib_sizes, expr_norm, res_meringue, p_vals)
}

p.values.Cauchy <- apply(matrix(unlist(p.values.allSample), nrow = 4, byrow = T), 2, ACAT::ACAT)
p.values.Cauchy.adj <- p.adjust(p.values.Cauchy, method = "BY")

p.values.allSample.adj <- p.adjust(unlist(p.values.allSample), method = "BY")

SVG_MERINGUE.allSample <- list()
for(iSample in 1:4){
  p.values.adj.perSample <- p.values.allSample.adj[((iSample - 1) * n.genes.afterQC + 1) : (iSample * n.genes.afterQC)]
  SVG_MERINGUE.allSample[[iSample]] <- genes.afterQC[p.values.adj.perSample < 0.05]
}

SVG.MERINGUE_union <- c()
SVG.MERINGUE_intersection <- genes.afterQC
for(iSample in 1:4){
  SVG.MERINGUE_union <- c(SVG.MERINGUE_union, SVG_MERINGUE.allSample[[iSample]])
  SVG.MERINGUE_intersection <- intersect(SVG.MERINGUE_intersection, SVG_MERINGUE.allSample[[iSample]])
}
SVG.MERINGUE_union <- unique(SVG.MERINGUE_union)

SVG.MERINGUE_cauchy <- genes.afterQC[p.values.Cauchy.adj < 0.05]
## -----------------------------------------------------------------------------
## Output the SVG list for MERINGUE
## -----------------------------------------------------------------------------
SVG.MERINGUE.list <- list()
SVG.MERINGUE.list[[1]] <- SVG.MERINGUE_union
SVG.MERINGUE.list[[2]] <- SVG.MERINGUE_intersection
SVG.MERINGUE.list[[3]] <- SVG.MERINGUE_cauchy
names(SVG.MERINGUE.list) <- c("MERINGUE_union",
                              "MERINGUE_intersection",
                              "MERINGUE_Cauchy")
saveRDS(SVG.MERINGUE.list, file = "./DLPFC/out/svgList_MERINGUE.RDS")
# ==============================================================================
# binSpect
## -----------------------------------------------------------------------------
genes.afterQC <- readRDS("./DLPFC/out/geneList_QC.RDS")
n.genes.afterQC <- length(genes.afterQC)
p.values.allSample <- list()

for(iSample in 1:4){
  cat(paste0("Working on sample ", dat$samples.id[iSample], "...\n"))
  # 1. Extract counts and spatial coordinates
  counts <- dat$counts[[iSample]][genes.afterQC,]
  pos <- as.matrix(dat$spatial.locs[[iSample]])

  # 2. Normalize the expression matrix (MERINGUE expects normalized data)
  # Standard log-normalization equivalent to Seurat's LogNormalize
  expr_norm <- normalizeCounts(counts)

  # 3. Compute Moran's I utilizing the SVG package's C++ acceleration
  res_binSpect <- CalSVG_binSpect(
    expr_matrix = as.matrix(expr_norm),
    spatial_coords = pos,
    n_threads = 1,
    adjust_method = "BY"
  )

  idx <- match(rownames(counts), res_binSpect$gene)
  # 4. Extract unadjusted p-values (ensure the exact column name matches your SVG version, e.g., 'p_value' or 'pval')
  p_vals <- res_binSpect$p.value[idx]

  # 5. Handle potential NAs (e.g., genes with 0 spatial variance) to prevent ACAT failure
  p_vals[is.na(p_vals)] <- 1.0

  p.values.allSample[[iSample]] <- p_vals

  # Clean up environment to manage memory across loops
  rm(counts, pos, lib_sizes, expr_norm, res_binSpect, p_vals, idx)
}

p.values.allSample.adj <- p.adjust(unlist(p.values.allSample), method = "BY")

p.values.allSample.matrix <- matrix(unlist(p.values.allSample), nrow = 4, byrow = T)
p.values.allSample.matrix[p.values.allSample.matrix <= 0] <- 1e-100
p.values.allSample.matrix[p.values.allSample.matrix >= 1] <- 1
p.values.Cauchy <- apply(p.values.allSample.matrix, 2, ACAT::ACAT)
p.values.Cauchy.adj <- p.adjust(p.values.Cauchy, method = "BY")

SVG_binSpect.allSample <- list()
for(iSample in 1:4){
  p.values.adj.perSample <- p.values.allSample.adj[((iSample - 1) * n.genes.afterQC + 1) : (iSample * n.genes.afterQC)]
  SVG_binSpect.allSample[[iSample]] <- genes.afterQC[p.values.adj.perSample < 0.05]
}

SVG.binSpect_union <- c()
SVG.binSpect_intersection <- genes.afterQC
for(iSample in 1:4){
  SVG.binSpect_union <- c(SVG.binSpect_union, SVG_binSpect.allSample[[iSample]])
  SVG.binSpect_intersection <- intersect(SVG.binSpect_intersection, SVG_binSpect.allSample[[iSample]])
}
SVG.binSpect_union <- unique(SVG.binSpect_union)

SVG.binSpect_cauchy <- genes.afterQC[p.values.Cauchy.adj < 0.05]

## -----------------------------------------------------------------------------
## Output the SVG list for binSpect
## -----------------------------------------------------------------------------
SVG.binSpect.list <- list()
SVG.binSpect.list[[1]] <- SVG.binSpect_union
SVG.binSpect.list[[2]] <- SVG.binSpect_intersection
SVG.binSpect.list[[3]] <- SVG.binSpect_cauchy
names(SVG.binSpect.list) <- c("binSpect_union",
                            "binSpect_intersection",
                            "binSpect_Cauchy")
saveRDS(SVG.binSpect.list, file = "./DLPFC/out/svgList_binSpect.RDS")
# ==============================================================================
# Combine data
## -----------------------------------------------------------------------------
SVG.INSiGHT.list <- readRDS(file = "./DLPFC/out/geneList/svgList_INSiGHT.RDS")
SVG.SPARKX.list <- readRDS(file = "./DLPFC/out/geneList/svgList_SPARKX.RDS")
SVG.SPARK.list <- readRDS(file = "./DLPFC/out/geneList/svgList_SPARK.RDS")
SVG.scBSP.list <- readRDS(file = "./DLPFC/out/geneList/svgList_scBSP.RDS")
SVG.MERINGUE.list <- readRDS(file = "./DLPFC/out/geneList/svgList_MERINGUE.RDS")
SVG.binSpect.list <- readRDS(file = "./DLPFC/out/geneList/svgList_binSpect.RDS")
SVG.SOMDE.list <- readRDS(file = "./DLPFC/out/geneList/svgList_SOMDE.RDS")


SVG.list <- list()
SVG.list <- append(SVG.list, SVG.INSiGHT.list)
SVG.list <- append(SVG.list, SVG.SPARKX.list)
SVG.list <- append(SVG.list, SVG.SPARK.list)
SVG.list <- append(SVG.list, SVG.scBSP.list)
SVG.list <- append(SVG.list, SVG.MERINGUE.list)
SVG.list <- append(SVG.list, SVG.binSpect.list)
SVG.list <- append(SVG.list, SVG.SOMDE.list)

METHODS <- c(
  "INSiGHT",
  "SPARK-X (U)", "SPARK-X (I)", "SPARK-X (C)",
  "SPARK (U)", "SPARK (I)", "SPARK (C)",
  "scBSP (U)", "scBSP (I)", "scBSP (C)",
  "MERINGUE (U)", "MERINGUE (I)", "MERINGUE (C)",
  "binSpect (U)", "binSpect (I)", "binSpect (C)",
  "SOMDE (U)", "SOMDE (I)", "SOMDE (C)"
)
names(SVG.list) <- METHODS

saveRDS(SVG.list, file = "./DLPFC/out/geneList/svgList.RDS")

