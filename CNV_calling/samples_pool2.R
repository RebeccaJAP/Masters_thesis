  # Load packages
library(stringi)
library(Seurat)
library(dplyr)


  # Choose the cell type annotation used
annot <- "twopass"
# annot <- "marker"

  # Read the original Seurat object
obj <- readRDS("otPools_all_integrated.RDS")

  # Use the updated metadata to get the marker-based cell types:
load("joined_meta.RData")

if (annot == "twopass") {
  joined_meta$cell_type_annot <- joined_meta$twoPassAnnotation_clean
} else {
  joined_meta$cell_type_annot <- joined_meta$cell_type_pred_y
}

obj@meta.data <- joined_meta

  # Only keep the pool and timepoint of interest
  # and only keep donor and cell type if the sample has more than 50 cells
sample_meta <- obj@meta.data %>%
  filter(SampleID == "OT_NI23-4_POOL2_d20")  %>%
  group_by(new_SampleID, new_donor, cell_type_annot) %>%
  summarize(n=n()) %>%
  filter(n>50) %>%
  arrange(desc(n))

for (donor in unique(sample_meta$new_donor)) {
    # Filter the metadata so that only donor of interest remain
  sample_donor_meta <- samples_meta %>%
    filter(new_donor == donor)
    
  for (cell_type in unique(sample_donor_meta$cell_type_annot)) {
      # Get the entire data based on the filtered metadata
    donor_subset_sample <- obj[,which(obj$SampleID == "OT_NI23-4_POOL2_d20" &
                                          obj$new_donor == donor)]
      
      # Mark references and observations
    donor_subset_sample@meta.data$status <- NA
    donor_subset_sample@meta.data$status[which(donor_subset_sample@meta.data$cell_type_annot != cell_type)] <- "reference"
    donor_subset_sample@meta.data$status[which(donor_subset_sample@meta.data$cell_type_annot == cell_type)] <- "observation"
      
      # Ensure the pool name and day are separated by "_" for easier processing later
    sample <- stri_replace_last(sample, "_", regex = "-")
      # Ensure there are no spaces in the name of the cell type
    cell_type <- stri_replace(cell_type, "-", regex = " ")
    
    save(donor_subset_sample,
           file = paste0("CNV_calling/Data/Pool2/", annot, "_", sample, "_", donor, "_", cell_type, ".RData"))
    }
  }
