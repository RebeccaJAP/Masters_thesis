  # Load packages
library(SCEVAN)
library(Seurat)

set.seed(97)

  # Set the annotation type used:
annot <- "twopass"
# annot <- "marker"

  # This file is intended to be run as a batch job
dataname <- commandArgs(trailingOnly = TRUE)

  # Load data with the correct cell type annotations
if (annot == "twopass") {
  load(paste0("CNV_calling/Data/twopass/",dataname))
} else {
  load(paste0("CNV_calling/Data/marker/",dataname))
}


#data <- read.table("ref_annotated_cell_types.txt", sep = "\n", header = FALSE)
#data <- t(data)


  # Separate other donors' cells from the observations
reference <- donor_subset_sample[,which(donor_subset_sample$status=="reference")]
observation <- donor_subset_sample[,which(donor_subset_sample$status=="observation")]

  # Create a reference by pseudobulking all cells of status "reference"
ref_cell <- as(rowMeans(reference@assays$RNA$counts), "sparseMatrix")
dimnames(ref_cell) <- list(rownames(donor_subset_sample@assays$RNA$counts), c("reference_cell"))

  # Combine the matrices of counts of the reference and the observations.
donor_subset_sample_mat <- cbind(observation@assays$RNA$counts, ref_cell)

  # Remove the ".RData" ending from the folder name
folder_name <- sub(".RData", "", dataname)

  # Run SCEVAN
sce_donor_subset_sample <- SCEVAN::pipelineCNA(donor_subset_sample_mat,
                                               ClonalCN = TRUE,
                                               par_cores = 6,
                                               norm_cell = "reference_cell",
                                               organism = "human",
                                               sample = folder_name,
                                               SCEVANsignatures = FALSE)

save(sce_donor_subset_sample, file=paste0("CNV_calling/SCEVAN/Results/", folder_name))
