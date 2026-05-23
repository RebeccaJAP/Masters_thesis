  # Load packages
library(GSEABase)
library(GSVA)
library(dplyr)
library(ggplot2)
library(tidyr)

  # Load the filtered pseudobulk object
load("Pseudobulk_data/pseudobulk_obj_filtered.RData")

  # Load the gene set files
gmt_file <- "Gene_expression_analyses/Gene_set_files/c2.cp.v2025.1.Hs.symbols.gmt"
gene_sets_c2 <- getGmt(gmt_file)

gmt_file2 <- "Gene_expression_analyses/Gene_set_files/c5.go.bp.v2025.1.Hs.symbols.gmt"
gene_sets_c5 <- getGmt(gmt_file2)

hallmark_file <- "Gene_expression_analyses/Gene_set_files/h.all.v2025.1.Hs.symbols.gmt"
gene_sets_h <- getGmt(hallmark_file)


  # Get the genes of the pathways of interest
wnt_pathway <- geneIds(gene_sets_c2[["KEGG_MEDICUS_REFERENCE_WNT_SIGNALING_PATHWAY"]])                 
bmp_pathway <- geneIds(gene_sets_c2[["KEGG_MEDICUS_REFERENCE_BMP_SIGNALING_PATHWAY"]])   
tgfb_pathway <- geneIds(gene_sets_c2[["KEGG_MEDICUS_REFERENCE_TGF_BETA_SIGNALING_PATHWAY"]])  
fgf_pathway <- geneIds(gene_sets_c2[["PID_FGF_PATHWAY"]])
shh_pathway <-  geneIds(gene_sets_c2[["BIOCARTA_SHH_PATHWAY"]])
egf_pathway <-  geneIds(gene_sets_c2[["BIOCARTA_EGF_PATHWAY"]])
prc2_pathway <- geneIds(gene_sets_c2[["BIOCARTA_PRC2_PATHWAY"]])
ra_pathway <-  geneIds(gene_sets_c2[["PID_RETINOIC_ACID_PATHWAY"]])
nodal_pathway <- geneIds(gene_sets_c2[["KEGG_MEDICUS_REFERENCE_NODAL_SIGNALING_PATHWAY"]])
emt_hallmark <- geneIds(gene_sets_h[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]])
cellTocell_adhesion <-  geneIds(gene_sets_c5[["GOBP_CELL_CELL_ADHESION"]])

  # List the pathways names
pathwaysToSelect <- c("wnt_pathway",
                      "bmp_pathway",
                      "tgfb_pathway",
                      "fgf_pathway",
                      "shh_pathway",
                      "egf_pathway",
                      "prc2_pathway",
                      "ra_pathway",
                      "nodal_pathway",
                      "emt_hallmark", 
                      "cellTocell_adhesion")

  # Combine the pathway genes into one object
gs <- sapply(pathwaysToSelect, function(x) get(x))

  # Create the object containing the GSVA parameters
gsvaPar <- gsvaParam(pseudobulk_obj@assays$RNA["data"], gs)

  # Run the GSVA
gsva.es <- gsva(gsvaPar, verbose=FALSE)

  # Combine the GSVA values for each pathway into one dataframe
fullMatrix_results <- t(gsva.es[1:dim(gsva.es)[1],1:dim(gsva.es)[2]])
fullMatrix_results <- as.data.frame(fullMatrix_results)

  # Get get the SampleID and donor from the rownames
fullMatrix_results$new_SampleID <- sapply(strsplit(as.character(rownames(fullMatrix_results)),"_"), function(x) x[1])
fullMatrix_results$new_donor <- sapply(strsplit(as.character(rownames(fullMatrix_results)),"_"), function(x) x[2])
fullMatrix_results$new_SampleID <- gsub("g", "", fullMatrix_results$new_SampleID)
rownames(fullMatrix_results) <- NULL

  # Add metadata from the pseudobulk object
fullMatrix_results <- left_join(fullMatrix_results,
                                distinct(subset(pseudobulk_obj@meta.data,
                                                select = c(new_SampleID, SampleID, Day_fixed, pool, donor, new_donor))))
fullMatrix_results$donor_short <- sub("^[^i]*?i-", "", fullMatrix_results$donor)

  # Change into a long dataframe
gsva_long <- as.data.frame(fullMatrix_results %>%
                             pivot_longer(-c("SampleID", "new_SampleID", "new_donor", "donor", "donor_short", "Day_fixed", "pool"),
                                          names_to = "Pathways",
                                          values_to = "GSVA_scores"))

  # Adjust the pathway names
gsva_long <- gsva_long %>%
  mutate(
    Pathways = case_when(
      Pathways %in% c("bmp_pathway", "fgf_pathway", "shh_pathway", "egf_pathway", "prc2_pathway", "ra_pathway") ~ paste(toupper(mapply(`[`, strsplit(Pathways, "_", fixed = TRUE), 1)), "signaling pathway", sep = " "),
      Pathways == "wnt_pathway" ~ "Wnt signaling pathway",
      Pathways == "cellTocell_adhesion" ~ "Cell-cell adhesion",
      Pathways == "tgfb_pathway" ~ "TGF\u03b2 signaling pathway",
      Pathways == "nodal_pathway" ~ "Nodal signaling pathway",
      Pathways == "emt_hallmark" ~ "EMT defining genes"))


  # Plot the GSVA values across time points and pools for each pathway
all_relev_paths_plot <- lapply(unique(gsva_long$Pathways), function(path) {
  data <- subset(gsva_long, Pathways == path)
  allrelevPaths_gsva <- ggplot(data, aes(y = donor_short, x = Day_fixed, fill = GSVA_scores))+
    geom_tile(col = "black") +
    facet_wrap(~pool, scales = "free_y") +
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      
      strip.text=element_text(size = 8, margin = margin(0,0.07,0,0.07, "cm")),
      
      legend.position="bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.margin = margin(0,0,0,0, "cm"),
      
      panel.spacing.x = unit(4, 'points'),
      panel.spacing.y = unit(1, 'points')
    ) +
    xlab("Day")+
    scale_fill_gradientn("GSVA score",
                         colours = c("#1600A3", "slateblue1", "white", "hotpink", "#D10069"),
                         values = NULL,
                         rescaler = ~ scales::rescale_mid(.x, mid = 0)) +
    ylab(NULL) +
    ggtitle(path)
  
  pdf(file=paste0("Gene_expression_analyses/gsvaScores_", path, "all_days.pdf"), width = 9.0, height = 5.8)
  plot(allrelevPaths_gsva)
  dev.off()
  
})
                                       

    ## Check if the GSVA values correlate with the pericyte proportion or R number

  # Combine the GSVA results with rest of the metadata
gsva_long <- left_join(gsva_long, distinct(subset(pseudobulk_obj@meta.data, select = c(SampleID, donor, proportion_Per, R))))

  # Plot the relationships
cor_plot <- ggplot(gsva_long, aes(x = GSVA_scores, y = proportion_Per, color = as.factor(Day_fixed))) +
  geom_point(alpha = 0.7) + 
  labs(colour = "Day") +
  xlab("GSVA score") +
  ylab("Pericyte proportion") +
  facet_wrap(~Pathways) +
  scale_color_manual(values = c("slateblue1", "hotpink",  "olivedrab2", "goldenrod1", "chocolate1"))

  # Save the plot
pdf(file=paste0("Gene_expression_analyses/gsva_cor_plot_all_pathways.pdf"), paper = "a5r")
plot(cor_plot)
dev.off()
