  # Load packages
library(GSEABase)
library(GSVA)
library(dplyr)
library(Seurat)


  # Join all meta data with the object
obj <- readRDS("otPools_all_integrated.RDS")
load("joined_meta.RData")
obj@meta.data <- joined_meta

  # Load the gene set files
gmt_file <- "Gene_expression_analyses/Gene_set_files/c2.cp.v2025.1.Hs.symbols.gmt"
gene_sets_c2 <- getGmt(gmt_file)

gmt_file2 <- "Gene_expression_analyses/Gene_set_files/c5.go.bp.v2025.1.Hs.symbols.gmt"
gene_sets_c5 <- getGmt(gmt_file2)

hallmark_file <- "Gene_expression_analyses/Gene_set_files/h.all.v2025.1.Hs.symbols.gmt"
gene_sets_h <- getGmt(hallmark_file)


    ## Pericyte typical
  
  # Get the genes defining the EMT
emt_hallmark <- geneIds(gene_sets_h[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]])

  # Add the EMT module scores
obj <- AddModuleScore(
  object = obj,
  features = list(emt_hallmark),
  name = 'HM_EMT_module_scores')


  ## Apoptosis

  # Get the genes from two different apoptosis gene sets
hm_apop_pathway <- geneIds(gene_sets_h[["HALLMARK_APOPTOSIS"]])
kegg_apop_pathway <- geneIds(gene_sets_c2[["KEGG_APOPTOSIS"]])   

  # Add the module scores for each gene set
obj <- AddModuleScore(
  object = obj,
  features = list(hm_apop_pathway),
  name = 'HM_apop_module_scores')

obj <- AddModuleScore(
  object = obj,
  features = list(kegg_apop_pathway),
  name = 'KEGG_apop_module_scores')


    ## Cell-to-cell adhesion

GOBP_cellTocell_adhesion <-  geneIds(gene_sets_c5[["GOBP_CELL_CELL_ADHESION"]])  
kegg_cam_pathway <- geneIds(gene_sets_c2[["KEGG_CELL_ADHESION_MOLECULES_CAMS"]])  

obj <- AddModuleScore(
  object = obj,
  features = list(GOBP_cellTocell_adhesion),
  name = 'GOBP_cellTocell_adhesion_module_scores')

obj <- AddModuleScore(
  object = obj,
  features = list(kegg_cam_pathway),
  name = 'KEGG_cellTocell_adhesion_module_scores')


    ## Neural specification

GOBP_anterior_posterior_pathway <- geneIds(gene_sets_c5[["GOBP_ANTERIOR_POSTERIOR_PATTERN_SPECIFICATION"]])  
GOBP_neuronal_commitment_pathway <- geneIds(gene_sets_c5[["GOBP_COMMITMENT_OF_NEURONAL_CELL_TO_SPECIFIC_NEURON_TYPE_IN_FOREBRAIN"]]) 
GOBP_negreg_npc_prolif_pathway <- geneIds(gene_sets_c5[["GOBP_NEGATIVE_REGULATION_OF_NEURAL_PRECURSOR_CELL_PROLIFERATION"]])  

obj <- AddModuleScore(
  object = obj,
  features = list(GOBP_anterior_posterior_pathway),
  name = 'GOBP_anterior_posterior_module_scores')

obj <- AddModuleScore(
  object = obj,
  features = list(GOBP_neuronal_commitment_pathway),
  name = 'GOBP_neuronal_commitment_module_scores')

obj <- AddModuleScore(
  object = obj,
  features = list(GOBP_negreg_npc_prolif_pathway),
  name = 'GOBP_negreg_npc_prolif_module_scores')


  ## Extracellular matrix

GOBP_posreg_ecm_organization_pathway <- geneIds(gene_sets_c5[["GOBP_POSITIVE_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION"]])  
GOBP_reg_ecm_organization_pathway <- geneIds(gene_sets_c5[["GOBP_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION"]])  

obj <- AddModuleScore(
  object = obj,
  features = list(GOBP_posreg_ecm_organization_pathway),
  name = 'GOBP_posreg_ecm_organization_module_scores')

obj <- AddModuleScore(
  object = obj,
  features = list(GOBP_reg_ecm_organization_pathway),
  name = 'GOBP_reg_ecm_organization_module_scores')


  # Extract and save the metadata with added module scores
meta <- obj@meta.data                           

save(meta, file = "Gene_expression_analyses/meta_module_scores.RData")



