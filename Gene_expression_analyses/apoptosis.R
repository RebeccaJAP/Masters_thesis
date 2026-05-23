  # Load packages
library(GSEABase)
library(GSVA)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Seurat)


    ## Apoptosis and GSVA

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
kegg_apop_pathway <- geneIds(gene_sets_c2[["KEGG_APOPTOSIS"]])   
hm_apop_pathway <- geneIds(gene_sets_h[["HALLMARK_APOPTOSIS"]])
gobp_apop_sign_pathway <- geneIds(gene_sets_c5[["GOBP_APOPTOTIC_SIGNALING_PATHWAY"]])
gobp_apop_devl_pathway <- geneIds(gene_sets_c5[["GOBP_APOPTOTIC_PROCESS_INVOLVED_IN_DEVELOPMENT"]])
gobp_apop_pathway <- geneIds(gene_sets_c5[["GOBP_APOPTOTIC_PROCESS"]])

  # List the pathways names
pathwaysToSelect <- c("kegg_apop_pathway",
                      "hm_apop_pathway", 
                      "gobp_apop_sign_pathway",
                      "gobp_apop_devl_pathway",
                      "gobp_apop_pathway")

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

  # Add the pathway gene counts
gsva_long <- gsva_long %>% mutate(pathway_gene_count = case_when(
  Pathways == "gobp_apop_devl_pathway" ~ 41,
  Pathways == "gobp_apop_pathway" ~ 1714,
  Pathways == "gobp_apop_sign_pathway" ~ 572,
  Pathways == "kegg_apop_pathway" ~ 79,
  Pathways == "hm_apop_pathway" ~ 151),
  pathway_gene_count = paste0("Pathway gene count: ", pathway_gene_count))


    # Check if the GSVA values correlate with the pericyte proportion

  # Combine the GSVA results with rest of the metadata
gsva_long <- left_join(gsva_long, distinct(subset(pseudobulk_obj@meta.data, select = c(SampleID, donor, proportion_Per))))

  # Plot the relationships
cor_plot <- ggplot(gsva_long, aes(x = GSVA_scores, y = proportion_Per, color = as.factor(Day_fixed))) +
  geom_point(alpha = 0.7) + 
  labs(colour = "Day") +
  xlab("GSVA score") +
  ylab("Pericyte proportion") +
  facet_wrap(vars(Pathways, pathway_gene_count),
             labeller = labeller(Pathways = 
                                   c(gobp_apop_devl_pathway = "GOBD apoptotic process involved in development",
                                     gobp_apop_pathway = "GOBD apoptotic process",
                                     gobp_apop_sign_pathway = "GOBD apoptotic signaling",
                                     hm_apop_pathway = "HM apoptotic process",
                                     kegg_apop_pathway = "KEGG apoptotic process"))) +
  scale_color_manual(values = c("slateblue1", "hotpink",  "olivedrab2", "goldenrod1", "chocolate1")) +
   stat_cor(inherit.aes = F,
    data = gsva_long,
     method = "kendall",
     cor.coef.name = "tau",
     aes(x=GSVA_scores, y=proportion_Per)) 

  # Save the plot
pdf(file=paste0("Gene_expression_analyses/gsva_cor_plot_apoptosis.pdf"), paper = "a5r")
plot(cor_plot)
dev.off()


    ## Pericytes and apoptosis module scores

  # Join all meta data with the object
obj <- readRDS("otPools_all_integrated.RDS")
load("joined_meta.RData")
obj@meta.data <- joined_meta


  # Add the module scores for the hallmark apoptosis pathway
obj <- AddModuleScore(
  object = obj,
  features = list(hm_apop_pathway),
  name = 'HM_apop_module_scores')

hm_apoptosis_violin <- ggplot(data = obj@meta.data, aes(x = as.factor(Day_fixed), y = HM_apop_module_scores1, fill = as.factor(Day_fixed))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  facet_wrap(~cell_type_group,
             labeller = labeller(cell_type_group = 
                                   c(nonPericyte = "Non-pericyte cell types",
                                     pericyte = "Pericytes"))) +
  scale_fill_manual(values = c("slateblue1", "hotpink", "limegreen", "mediumorchid1", "chocolate1")) +
  ylab("HM apoptosis module score") +
  xlab("Day") +
  theme(legend.position="none")

  # Save the plot
pdf(file=paste0("Gene_expression_analyses/apoptosis_violin.pdf"), paper = "a5r")
plot(hm_apoptosis_violin)
dev.off()
