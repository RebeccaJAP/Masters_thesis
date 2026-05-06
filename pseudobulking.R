library(Seurat)
library(dplyr)
library(tidyr)


  # Create the pseudobulk object and add metadata
# Choose type = 1, 2, 3, or 4:
# 1 pseudobulked by donor & SampleID, not filtered
# 2 pseudobulked by donor & SampleID, filtered
# 3 pseudobulked by donor, SampleID & cell type category (pericyte/non-pericyte), filtered
# 4 pseudobulked by donor, SampleID & cell type category (pericyte/non-pericyte), filtered, and split by cell type category

pseudobulking <- function(type = 1) {
  
    # Load the original object
  obj <- readRDS("otPools_all_integrated.RDS")
  
    # Load the whole meta data
  load("joined_meta.RData")
  
    # Ensure the original rownames are not lost while joining variables to the original object
  rows <- rownames(obj@meta.data)
  
    # Select the variables necessary for pseudobulking & later joining
  add_meta <- joined_meta %>%
    select(robustID, new_donor, new_SampleID, freq_donor, cell_type_group)
  
  obj@meta.data <- left_join(obj@meta.data, add_meta)
  
    # Compute the number of cells per group
  cell_counts <- obj@meta.data %>%
    group_by(SampleID, donor, cell_type_group) %>%
    mutate(n_cell_type_group=n())
  
  obj@meta.data <- left_join(obj@meta.data, cell_counts)
  rownames(obj@meta.data) <- rows
  
    # Filter the data
  if(type == 2) {
    obj <- obj[, which(obj$freq_donor >= 50)]
  }
  else if(type >= 3) {
    obj <- obj[, which(obj$n_cell_type_group >= 50)]
  }
  
    # Pseudobulk the data
  if(type <= 2) {
    pseudobulk_obj <- AggregateExpression(obj,
                                          group.by = c("new_SampleID", "new_donor"),
                                          return.seurat = T)
  }
  else {
    pseudobulk_obj <- AggregateExpression(obj,
                                          group.by = c("new_SampleID", "new_donor", "cell_type_group"),
                                          return.seurat = T)
  }
  
    # Remove the extra letter added to some SampleIDs by AggregateExpression
  pseudobulk_obj@meta.data <- pseudobulk_obj@meta.data %>%
    mutate(new_SampleID = gsub("g", "", new_SampleID))
  
    # Update the variable with rownames after filtering
  rows <- rownames(pseudobulk_obj@meta.data)
  
  
      # Select relevant variables:
      
    # All donor-level metadata fits the pseudobulk objects
  donor_meta <- read.table("hipsci.qc1_sample_info.20170927.tsv", sep="\t", header=T)
  colnames(donor_meta)[which(names(donor_meta) == "donor")] <- "donor_short"
  colnames(donor_meta)[which(names(donor_meta) == "name")] <- "donor"
  cols <- intersect(colnames(joined_meta), colnames(donor_meta))
  
    # Select certain variables from the original metadata and from the ones computed later
  cols <- append(cols, c("new_SampleID", "new_donor", "SampleID", "Day_fixed", "donorSample", 
            "pool", "geneMutatedOrCorrected", "proportion_d0",
            "proportion_day_fixed", "first_day", "donors_per_pool",
            "prop_patient_orig", "prop_patient", "R",
            "sendai", "sendai_approx", "freq_donor",
            "diff_method", "z_R", "z_pluri_novelty"))
  
  if(type <= 2) {
      # Add cell type proportions to the variables of interest when all cell types are pseudobulked together
    cols <- append(cols, c("proportion_oRG", "proportion_vRG",
                           "proportion_panRG_O", "proportion_Per", 
                           "proportion_PgS", "proportion_PgG2M",
                           "z_proportion_Per", "z_proportion_vRG",
                           "z_proportion_panRG_O", "z_proportion_PgS",
                           "z_proportion_oRG", "z_proportion_PgG2M"))
  }
  
    # Extract the relevant meta data
  pseudobulk_meta <- joined_meta %>%
    select(all_of(cols)) %>%
    distinct()
  
  pseudobulk_obj@meta.data <- left_join(pseudobulk_obj@meta.data, pseudobulk_meta)
  
  rownames(pseudobulk_obj@meta.data) <- rows
  
  if(type == 4) {
      # Split the object to separate pericytes and non-pericytes
    pseudobulk_obj <- SplitObject(pseudobulk_obj, split.by = "cell_type_group")
  }
  
  return(pseudobulk_obj)
}


  # Add cell cycle scores
    # Input pericytes and other cell types of the split object separately 
pseudobulk_cell_cycle_scoring <- function(pseudobulk_obj) {
  
    # Load the cell cycle phase marker genes
  s.genes <- cc.genes$s.genes 
  g2m.genes <- cc.genes$g2m.genes 
  
  pseudobulk_counts <- pseudobulk_obj@assays$RNA$counts 
  
    # Normalize the counts
  norm_pseudobulk_counts <- pseudobulk_counts %*% diag (mean(colSums(pseudobulk_counts))/colSums(pseudobulk_counts))  
  
    # Compute the cell cycle scores for synthesis and G2M
  d_score <- cbind.data.frame(   
    orig.ident = colnames(pseudobulk_counts),  
    s = colSums(log10(norm_pseudobulk_counts+1)[s.genes[s.genes %in% rownames(norm_pseudobulk_counts)],]),  
    g2m = colSums(log10(norm_pseudobulk_counts+1)[g2m.genes[g2m.genes %in% rownames(norm_pseudobulk_counts)],]))
  
    # Add the cell cycle scores to the pseudobulk object
  pseudobulk_obj@meta.data <- left_join(pseudobulk_obj@meta.data, d_score, by = "orig.ident")
  
    # Compute the z scores for the cell cycle scores
  pseudobulk_obj$z_s <- (d_score$s - mean(d_score$s) ) / sd(d_score$s)  
  pseudobulk_obj$z_g2m <- (d_score$g2m - mean(d_score$g2m) ) / sd(d_score$g2m)
  
    # Ensure the rownames are retained
  row.names(pseudobulk_obj@meta.data) <- pseudobulk_obj$orig.ident
  
  return(pseudobulk_obj)
  
}


# Create and save all necessary pseudobulk objects

pseudobulk_obj <- pseudobulking(type = 1)
pseudobulk_obj <- pseudobulk_cell_cycle_scoring(pseudobulk_obj)
save(pseudobulk_obj, file = "Pseudobulk/pseudobulk_obj_unfiltered.RData")


pseudobulk_obj <- pseudobulking(type = 2)
pseudobulk_obj <- pseudobulk_cell_cycle_scoring(pseudobulk_obj)
save(pseudobulk_obj, file = "Pseudobulk/pseudobulk_obj_filtered.RData")


pseudobulk_obj <- pseudobulking(type = 3)
pseudobulk_obj <- pseudobulk_cell_cycle_scoring(pseudobulk_obj)
save(pseudobulk_obj, file = "Pseudobulk/pseudobulk_obj_pericytes_others_filtered.RData")


pseudobulk_obj <- pseudobulking(type = 4)
pseudobulk_obj <- pseudobulk_cell_cycle_scoring(pseudobulk_obj$nonPericyte)
save(pseudobulk_obj, file = "Pseudobulk/pseudobulk_obj_others_filtered.RData")

pseudobulk_obj <- pseudobulk_cell_cycle_scoring(pseudobulk_obj$pericyte)
save(pseudobulk_obj, file = "Pseudobulk/pseudobulk_obj_pericytes_filtered.RData")



