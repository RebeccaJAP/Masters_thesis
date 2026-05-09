  # Load packages
library(SCEVAN)
library(Seurat)

set.seed(97)

  # Set the annotation type used:
annot <- "twopass"
# annot <- "marker"

  # Set reference type
ref_type <- "pseudobulk" # reference cells pseudobulked into one
#ref_type <- "sce_database" # no additional reference input
#ref_type <- "individual_cells" # reference cells used without pseudobulking; SCEVAN matches to closest reference

  # This file is intended to be run as a batch job
dataname <- commandArgs(trailingOnly = TRUE)

  # Load data with the correct cell type annotations
if (annot == "twopass") {
  load(paste0("CNV_calling/Data/twopass/",dataname))
} else {
  load(paste0("CNV_calling/Data/marker/",dataname))
}

  # Separate other donors' cells from the observations
reference <- donor_subset_sample[,which(donor_subset_sample$status=="reference")]
observation <- donor_subset_sample[,which(donor_subset_sample$status=="observation")]

  # Remove the ".RData" ending from the folder name
folder_name <- sub(".RData", "", dataname)


if (ref_type == "pseudobulk") {
    # Set version number
  V <- "V1"
  
    # Create a reference by pseudobulking all cells of status "reference"
  ref_cell <- as(rowMeans(as.matrix(reference@assays$RNA$counts)), "sparseMatrix")
  dimnames(ref_cell) <- list(rownames(donor_subset_sample@assays$RNA$counts), c("reference_cell"))
  
    # Combine the matrices of counts of the reference and the observations.
  donor_subset_sample_mat <- cbind(observation@assays$RNA$counts, ref_cell)
  
    # Run SCEVAN
  sce_donor_subset_sample <- SCEVAN::pipelineCNA(donor_subset_sample_mat,
                                                 ClonalCN = TRUE,
                                                 par_cores = 6,
                                                 norm_cell = "reference_cell",
                                                 organism = "human",
                                                 sample = folder_name,
                                                 SCEVANsignatures = FALSE)
  
} else if(ref_type == "individual_cells") {
    # Set version number
  V <- "V2"
  
    # Get the names of reference cells
  reference_names <- colnames(reference)
  
    # Run SCEVAN
  sce_donor_subset_sample <- SCEVAN::pipelineCNA(donor_subset_sample@assays$RNA$counts,
                                                 ClonalCN = TRUE,
                                                 par_cores = 6,
                                                 norm_cell = reference_names,
                                                 organism = "human",
                                                 sample = folder_name,
                                                 SCEVANsignatures = FALSE)
  
} else if(ref_type == "sce_database") {
    # Set version number
  V <- "V3"
  
    # Run SCEVAN
  sce_donor_subset_sample <- SCEVAN::pipelineCNA(observation@assays$RNA$counts,
                                                 ClonalCN = TRUE,
                                                 par_cores = 6,
                                                 organism = "human",
                                                 sample = folder_name,
                                                 SCEVANsignatures = TRUE) # Reference from SCEVAN database used
}

save(sce_donor_subset_sample, file=paste0("CNV_calling/SCEVAN/Results/", V, "_", dataname))
