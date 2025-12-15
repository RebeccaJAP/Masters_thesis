  # Load packages
library(SCEVAN)
library(Seurat)

set.seed(97)

  # Set the annotation type used:
# annot <- "twopass"
annot <- "marker"



dataname <- commandArgs(trailingOnly = TRUE)


  # Load data subsetted by SampleID and cell type

if (annot == "twopass") {
   load(paste0("Karyotyped_donors/ref_annotated_data_d20/",dataname))}
else (annot == "marker") {
  load(paste0("Karyotyped_donors/ref_annotated_data_new_cell_types/",dataname))
  }

  # Both files contain a Seurat object called "donor_subset_sample".
  # The object has all of the sample's donors but the "status" column has value "observation" for the main donor and "reference" for the others.


  # Separate other cell lines' cells from the observations

reference <- donor_subset_sample[,which(donor_subset_sample$status=="reference")]
observation <- donor_subset_sample[,which(donor_subset_sample$status=="observation")]


  # Pseudobulk the reference cells

dummy_cell <- as(rowMeans(reference@assays$RNA$counts), "sparseMatrix")
dimnames(dummy_cell) <- list(rownames(donor_subset_sample@assays$RNA$counts), c("reference_cell"))


  # Combine the matrices of counts of the pseudobulk reference and the observation cells.

donor_subset_sample_mat <- cbind(observation@assays$RNA$counts, dummy_cell)

  # Run SCEVAN

sce_donor_subset_sample <- SCEVAN::pipelineCNA(donor_subset_sample_mat,
                                               ClonalCN = T,
                                               par_cores = 6,
                                               norm_cell = "reference_cell",
                                               SUBCLONES = T,               # Make sure the CNVs are analysed in addition to separating normal and non-normal cells from each other.
                                               beta_vega = 0.5,             # Higher value results in coarser segmentation.
                                               organism = "human",
                                               sample = paste(dataname, annot, sep = "_"),
                                               SCEVANsignatures = F)    # Stop SCEVAN from using references other than the pseudobulk reference input.

  # Save the results
save(sce_donor_subset_sample, file=paste0("sce_results/", dataname, "_", annot))
