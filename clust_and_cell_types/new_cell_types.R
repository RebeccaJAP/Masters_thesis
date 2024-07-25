.libPaths(c("/projappl/project_2010414/project_rpackages_4.4.0", .libPaths()))

library(Seurat)


	# Load the reclustered Seurat object

load("/scratch/project_2010414/ot/saved/Reclustering/V4_obj.RData")


	# Read the new cell types

cell_types <- read.csv("/scratch/project_2010414/ot/saved/desired_columns.csv")


	# Matching is not needed: the order of cell types is the same as the order of the object

# cell_types_matched <- cell_types[match(cell_types$robustID, new_obj_V4$robustID),]


	# Add the cell type column to the Seurat object meta data

new_obj_V4@meta.data <- merge(new_obj_V4@meta.data, cell_types)


save(new_obj_V4, file=paste0("/scratch/project_2010414/ot/saved/Reclustering/V4_obj_cell_types.RData"))
