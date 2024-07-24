.libPaths(c("/projappl/project_2010414/project_rpackages_4.3.2", .libPaths()))
libpath <- .libPaths()[1]

library(SCEVAN)
library(Seurat)


# SCEVAN paper warns the results might be very bad if the initial choice of normal cells is wrong. Thus, more than one seed should be tested.

#set.seed(34)
set.seed(97)


dataname <- commandArgs(trailingOnly = TRUE)


# Load either clustered or original data subsetted by sample and cluster / cell type.
# Both files contain a Seurat object called "donor_subset_sample".
# The object has all of the sample's donors but the "status" column has value "observation" for the main donor and "reference" for the others.

load(paste0("/scratch/project_2010414/ot/saved/Method_tests/Karyotyped_donors/ref_annotated_data_clust/",dataname))
# load(paste0("/scratch/project_2010414/ot/saved/Method_tests/Karyotyped_donors/ref_annotated_data_d20/",dataname))


# Separate other donors' cells from the observations

reference <- donor_subset_sample[,which(donor_subset_sample$status=="reference")]
observation <- donor_subset_sample[,which(donor_subset_sample$status=="observation")]


# Create a dummy cell to use as the reference

dummy_cell <- as(rowMeans(reference@assays$RNA$counts), "sparseMatrix")
dimnames(dummy_cell) <- list(rownames(donor_subset_sample@assays$RNA$counts), c("reference_cell"))


# Combine the matrices of counts of the dummy cell and the observations.

donor_subset_sample_mat <- cbind(observation@assays$RNA$counts, dummy_cell)


sce_donor_subset_sample <- SCEVAN::pipelineCNA(donor_subset_sample_mat,
                                               ClonalCN = TRUE,
                                               par_cores = 6,
                                               norm_cell = "reference_cell",
                                            #  SUBCLONES = T,               # This makes sure the CNVs are analysed in addition to separating normal and non-normal cells from each other.
                                            #  beta_vega = 0.5,             # Higher value results in coarser segmentation.
                                            #  AdditionalGeneSets =         # This could be used to add signatures of normal cells but I have not used it or found an example of using it.
                                               organism = "human",
                                               sample = paste0("V4_", dataname),
                                               SCEVANsignatures = FALSE)    # This stops SCEVAN from using references other than the dummy cell.

save(sce_donor_subset_sample, file=paste0("sce_results/V4_", dataname))

