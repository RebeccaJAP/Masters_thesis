.libPaths(c("/projappl/project_2010414/project_rpackages_4.4.0", .libPaths()))

library(Seurat)


	# Load the Seurat object with new clusters

load("/scratch/project_2010414/ot/saved/Reclustering/reclustered_obj.RData")


	# Read the new cell types

cell_types <- read.csv("/scratch/project_2010414/ot/saved/desired_columns.csv")


	# Matching is not necessary here: the order of cell types is the same as the order of the object

# cell_types_matched <- cell_types[match(cell_types$robustID, obj$robustID),]


	# Add the cell type column to the Seurat object meta data

obj@meta.data <- merge(obj@meta.data, cell_types)


save(obj, file=paste0("/scratch/project_2010414/ot/saved/Reclustering/reclustered_obj_cell_types.RData"))
