  # Load packages
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

  # Join all meta data with the object
  # Read the original Seurat object
obj <- readRDS("obj.RDS")
load("joined_meta.RData")
obj@meta.data <- left_join(obj@meta.data, joined_meta)
rownames(obj@meta.data) <- colnames(obj@assays$RNA$counts)

  # Set the pericyte marker genes of interest
peri_markers <- c("DLC1", "DCN", "SNAI2", "TWIST1")

  # Plot the marker gene expression for each time point individually
for (day in unique(obj$Day_fixed)) {
  obj_peri <- subset(obj, Day_fixed == day)
    
    # Each plot shows the expression per cell type
  dotplot <- DotPlot(obj_peri,
                     features = peri_markers,
                     group.by = "twoPassAnnotation_clean",
                     cols = c("slateblue1", "hotpink")) +
    RotatedAxis() +
    ggtitle(paste0("Marker gene expression on day ", day)) +
    scale_color_gradientn(colors=c("#1600A3", "slateblue1", "white", "hotpink", "#D10069"),
                          rescaler = ~ scales::rescale_mid(.x, mid = 0))
  
  pdf(file=paste0("Gene_expression_analyses/Pericyte_marker_expression_d", day, ".pdf"),
                  paper = "legal")
  plot(dotplot)
  dev.off()
}



