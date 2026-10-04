## =============================================================================
## Gene set analysis
## =============================================================================
rm(list = ls())

library(Seurat)
library(dplyr)
library(tidyr)
library(methods)
library(ggplot2)
library(INSiGHT)
library(clusterProfiler)
library(enrichplot)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

SVG.list <- readRDS("./AD/out/svgList/svgList_AD.RDS")
SVG.INSiGHT.AD <- SVG.list$INSiGHT
SVG.list <- readRDS("./AD/out/svgList/svgList_earlyAD.RDS")
SVG.INSiGHT.earlyAD <- SVG.list$INSiGHT
SVG.list <- readRDS("./AD/out/svgList/svgList_Control.RDS")
SVG.INSiGHT.Control <- SVG.list$INSiGHT

SVG.INSiGHT.list <- list(
  "Control" = SVG.INSiGHT.Control,
  "Early-stage" = SVG.INSiGHT.earlyAD,
  "Late-stage" = SVG.INSiGHT.AD
)
rm(SVG.INSiGHT.Control, SVG.INSiGHT.earlyAD, SVG.INSiGHT.AD)

SVG.INSiGHT.unique <- list(
  "Control" = setdiff(SVG.INSiGHT.list[["Control"]],
                    union(SVG.INSiGHT.list[["Early-stage"]], SVG.INSiGHT.list[["Late-stage"]])),
  "Early-stage" = setdiff(SVG.INSiGHT.list[["Early-stage"]],
                    union(SVG.INSiGHT.list[["Control"]], SVG.INSiGHT.list[["Late-stage"]])),
  "Late-stage" = setdiff(SVG.INSiGHT.list[["Late-stage"]],
               union(SVG.INSiGHT.list[["Control"]], SVG.INSiGHT.list[["Early-stage"]]))
)

background.gene.list <- list(
  "Control" = readRDS("./AD/out/geneList_QC_Control.RDS"),
  "Early-stage" = readRDS("./AD/out/geneList_QC_earlyAD.RDS"),
  "Late-stage" = readRDS("./AD/out/geneList_QC_AD.RDS")
)

# ------------------------------------------------------------------------------
target.list <- SVG.INSiGHT.list
# ------------------------------------------------------------------------------
## Step 1. Run GO enrichment for each
ego.list <- lapply(names(target.list), function(m){
  ego <- enrichGO(gene         = target.list[[m]],
                  OrgDb        = org.Hs.eg.db,
                  keyType      = "SYMBOL",   # adjust if IDs are ENSEMBL/ENTREZ
                  ont          = "BP",
                  pAdjustMethod= "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  universe = background.gene.list[[m]])
  ego@result$Stage <- m   # tag the method
  return(ego@result)
})
saveRDS(ego.list, file = "./AD/out/GSEA_GO_output")
ego.list <- readRDS("./AD/out/GSEA_GO_output")
# ------------------------------------------------------------------------------
## Step 2. Combine results
ego.df <- bind_rows(ego.list)
# ------------------------------------------------------------------------------
## Step 3. Make a dotplot across methods
# Take top 5 enriched terms per method
ego.top <- ego.df %>%
  group_by(Stage) %>%
  filter(p.adjust < 0.05) %>%
  slice_min(p.adjust, n = 10) %>%
  ungroup()

ego.top$logPadjust <- -log10(ego.top$p.adjust)

# Set order of x-axis
ego.top$Stage <- factor(ego.top$Stage, levels = c("Control", "Early-stage", "Late-stage"))

ego.top <- ego.top %>%
  arrange(Stage, p.adjust)
# ------------------------------------------------------------------------------
plot_gsea <- ggplot(ego.top, aes(x = Stage, y = Description)) +
  geom_point(aes(size = Count, color = logPadjust)) +
  scale_color_continuous(low = "deepskyblue2", high = "firebrick1", name = "-log10(Adj. p-value)") +
  scale_size_continuous(name = "Gene count", range = c(13,23)) +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(size = 60, face = "bold", angle = 25, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 60, face = "bold"),
        legend.title = element_text(size = 60, face = "bold"),
        legend.text = element_text(size = 60)) +
  labs(x = "", y = "", title = "") +
  guides(color = guide_colourbar(barwidth = 3, barheight = 18))

ggsave(plot_gsea, file = "./AD/results/GSEA_plot_INSiGHT_allSVG.pdf",
       width = 49, height = 20, dpi = 200)
# ------------------------------------------------------------------------------
ego.top$Description[ego.top$Stage == "Control"]
ego.top$Description[ego.top$Stage == "Early-stage"]
ego.top$Description[ego.top$Stage == "Late-stage"]
