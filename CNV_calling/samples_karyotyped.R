  # Load packages
library(stringi)
library(Seurat)
library(dplyr)


  # Choose the cell type annotation used
annot <- "twopass"
# annot <- "marker"

  # Read the original Seurat object
obj <- readRDS("obj.RDS")

  # Use the updated metadata to get the marker-based cell types:
load("joined_meta.RData")

if (annot == "twopass") {
  joined_meta$cell_type_annot <- joined_meta$twoPassAnnotation_clean
} else {
  joined_meta$cell_type_annot <- joined_meta$cell_type_pred_y
}

obj@meta.data <- joined_meta

  # Only keep samples with the relevant donors and more than 50 cells

karyotyped_donors <- c("HPSI1016i-riwg_2", "HPSI0314i-hoik_1", "HPSI0316i-ierp_4", "HEL_318.3")

karyotyped_meta <- obj@meta.data %>%
  filter(donor %in% karyotyped_donors &
           Day_fixed == 20)  %>%
  group_by(new_SampleID, new_donor, cell_type_annot) %>%
  summarize(n=n()) %>%
  filter(n>50) %>%
  arrange(desc(n))

for (sample in unique(karyotyped_meta$new_SampleID)) {
    # Filter the metadata so that only samples of interest remain
  samples_meta <- karyotyped_meta %>%
    filter(new_SampleID == sample)
  
  for (donor in unique(samples_meta$new_donor)) {
      # Filter the metadata so that only donors of interest remain
    sample_donor_meta <- samples_meta %>%
      filter(new_donor == donor)
    
    for (cell_type in unique(sample_donor_meta$cell_type_annot)) {
      # Get the entire data based on the filtered metadata
      # and filter so that each subset consists of one cell type
      
      donor_subset_sample <- obj[,which(obj$new_SampleID == sample
                                        & obj$cell_type_annot == cell_type)]
      
        # Mark references and observations
      donor_subset_sample@meta.data$status <- NA
      donor_subset_sample@meta.data$status[which(donor_subset_sample@meta.data$new_donor != donor)] <- "reference"
      donor_subset_sample@meta.data$status[which(donor_subset_sample@meta.data$new_donor == donor)] <- "observation"
      
        # Ensure the pool name and day are separated by "_" for easier processing later
      sample <- stri_replace_last(sample, "_", regex = "-")
        # Ensure there are no spaces in the name of the cell type
      cell_type <- stri_replace(cell_type, "-", regex = " ")
      
      save(donor_subset_sample,
           file = paste0("CNV_calling/Data/", annot, "/", annot, "_", sample, "_", donor, "_", cell_type, ".RData"))
    }
  }
}