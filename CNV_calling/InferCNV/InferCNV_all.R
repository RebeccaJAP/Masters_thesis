# Load packages
library(infercnv)
library(Seurat)

set.seed(97)

  # Set the annotation type used:
annot <- "twopass"
# annot <- "marker"

  # Should the reference cells be of the same cell type of different donors or
  # different cell types of the same donor (only used in pool OT_NI23-4_POOL2_d20)?
ref_variable <- "donor"
# ref_variable <- "cell_type"

  # Set version number (update to V2, V3...)
V <- V1

  # This file is intended to be run as a batch job
dataname <- commandArgs(trailingOnly = TRUE)

  # Load data with the correct type of reference
if (ref_variable == "cell_type") {
  load(paste0("CNV_calling/Data/Pool2/",dataname))
} else {
  load(paste0("CNV_calling/Data/",dataname))
}

# Both files contain a Seurat object called "donor_subset_sample".
# Each cells of interest have value "observation" in column "status", while the others have value "reference".

  # Create the annotations file for inferCNV:
annotations <- subset(donor_subset_sample@meta.data, select=c(status))

  # Add the cell type (or donor name) to the annotation to make reading the results easier
if(ref_variable == "donor") {
  constant_var <- donor_subset_sample$cell_type_annot[1]
} else if (ref_variable == "cell_type") {
  constant_var <- donor_subset_sample$donor[1]
}
annotations$status <- paste(annotations$status, constant_var, sep="_")

  # Create the inferCNV object
inf_donor_subset_sample = CreateInfercnvObject(raw_counts_matrix = donor_subset_sample@assays$RNA$counts,
                                               annotations_file = annotations,
                                               delim="\t",
                                               gene_order_file = "hg38_gencode_v27.txt",
                                               ref_group_names = c(paste0("reference_", constant_var))) 


  # Remove the ".RData" ending from the folder name
folder_name <- sub(".RData", "", dataname)
folder_name <- paste(V, folder_name, sep="_")

  # Ensure the results go to the correct folder
if (ref_variable == "cell_type") {
  folder_name <- paste0("Pool2/", folder_name)
}

  # Run inferCNV
inf_donor_subset_sample = infercnv::run(inf_donor_subset_sample,
                                        cutoff=0.1,                          # Recommended value for 10x data
                                        out_dir= paste0("CNV_calling/InferCNV/Results/", folder_name),
                                        # cluster_by_groups = F,               # Not relevant when only one observation group is run per time
                                        denoise=T,
                                        # BayesMaxPNormal = 0.5,               # Smaller value filters out more of the uncertain CNVs
                                        HMM=T,
                                        # window_length = 31,                  # Needs to be an odd integer. With smaller value smaller CNVs should be detected. Default 101.
                                        # HMM_type = "i6",                     
                                        ref_subtract_use_mean_bounds = F,    # Use the average of all reference cells as the reference. With = T the average of subclusters would be used.
                                        # scale_data = T,
                                        # sd_amplifier = 1.5,
                                        # smooth_method = "coordinates",
                                        analysis_mode = "samples",            # "cells" should be better for sc but did not lead to better results here (and is slower); "subclusters" works best when running all donors at once
                                        # leiden_resolution = 0.01,            # Smaller value leads to less subclusters
                                        # tumor_subcluster_pval = 0.05,
                                        resume_mode = F)