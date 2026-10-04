rm(list = ls())

library(Seurat)
library(dplyr)
library(tidyr)
library(methods)
library(ggplot2)
library(INSiGHT)
SVG.list <- readRDS("./PDAC/out/svgList.RDS")

names(SVG.list)
names(SVG.list) <- c(
  "INSiGHT",
  "SPARK-X (U)", "SPARK-X (I)", "SPARK-X (C)",
  "SPARK (U)", "SPARK (I)","SPARK (C)",
  "scBSP (U)", "scBSP (I)", "scBSP (C)",
  "MERINGUE (U)", "MERINGUE (I)", "MERINGUE (C)",
  "binSpect (U)", "binSpect (I)","binSpect (C)",
  "SOMDE (U)", "SOMDE (I)", "SOMDE (C)"
)

## =============================================================================
## GO enrichment analysis
## =============================================================================
library(clusterProfiler)
library(enrichplot)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
background_genes <- readRDS("./PDAC/out/geneList_QC.RDS")

## Step 1. Run GO enrichment for each
ego.list <- lapply(names(SVG.list), function(m){
  ego <- enrichGO(gene         = SVG.list[[m]],
                  OrgDb        = org.Hs.eg.db,
                  keyType      = "SYMBOL",   # adjust if IDs are ENSEMBL/ENTREZ
                  ont          = "BP",
                  pAdjustMethod= "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  universe = background_genes)
  ego@result$Method <- m   # tag the method
  return(ego@result)
})
# saveRDS(ego.list, file = "./PDAC/out/GSEA_GO_output.rds")
# ego.list <- readRDS("./PDAC/out/GSEA_GO_output.rds")

## Step 2. Combine results
ego.df <- bind_rows(ego.list)

## Step 3. Make a dotplot across methods
# Take top 5 enriched terms per method
ego.top <- ego.df %>%
  group_by(Method) %>%
  filter(p.adjust < 0.05) %>% 
  slice_min(p.adjust, n = 10) %>%
  ungroup()

ego.top$logPadjust <- -log10(ego.top$p.adjust)
ego.top$Method <- factor(ego.top$Method, levels = c(
  "INSiGHT",
  "SPARK-X (U)", "SPARK-X (I)", "SPARK-X (C)",
  "SPARK (U)", "SPARK (I)","SPARK (C)",
  "scBSP (U)", "scBSP (I)", "scBSP (C)",
  "MERINGUE (U)", "MERINGUE (I)", "MERINGUE (C)",
  "binSpect (U)", "binSpect (I)","binSpect (C)",
  "SOMDE (U)", "SOMDE (I)", "SOMDE (C)"
))

plot_gsea <- ggplot(ego.top, aes(x = Method, y = Description)) +
  geom_point(aes(size = Count, color = logPadjust)) +
  scale_color_continuous(low = "deepskyblue2", high = "firebrick1", name = "-log10(adj. p-value)") +
  scale_size_continuous(name = "Gene count", range = c(20,30)) +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(size = 70, face = "bold", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 80, face = "bold"),
        legend.title = element_text(size = 60, face = "bold"),
        legend.text = element_text(size = 60)) +
  labs(x = "", y = "",
       title = "") + 
  guides(
    color = guide_colourbar(barwidth = 3, barheight = 18)
  )

ggsave(plot_gsea, file = "./PDAC/results/GSEA_allSVG_plot.pdf",
       width = 70, height = 50, dpi = 200, limitsize = F)

ego.top$Description[ego.top$Method == "INSiGHT"]
ego.top$Description[ego.top$Method == "scBSP (U)"]
ego.top$Description[ego.top$Method == "scBSP (I)"]
ego.top$Description[ego.top$Method == "scBSP (C)"]
ego.top$Description[ego.top$Method == "SPARK (U)"]
ego.top$Description[ego.top$Method == "SPARK (I)"]
ego.top$Description[ego.top$Method == "SPARK (C)"]
ego.top$Description[ego.top$Method == "SPARK-X (U)"]
ego.top$Description[ego.top$Method == "SPARK-X (I)"]
ego.top$Description[ego.top$Method == "SPARK-X (C)"]
ego.top$Description[ego.top$Method == "MERINGUE (U)"]
ego.top$Description[ego.top$Method == "MERINGUE (I)"]
ego.top$Description[ego.top$Method == "MERINGUE (C)"]
ego.top$Description[ego.top$Method == "binSpect (U)"]
ego.top$Description[ego.top$Method == "binSpect (I)"]
ego.top$Description[ego.top$Method == "binSpect (C)"]
ego.top$Description[ego.top$Method == "SOMDE (U)"]
ego.top$Description[ego.top$Method == "SOMDE (I)"]
ego.top$Description[ego.top$Method == "SOMDE (C)"]
