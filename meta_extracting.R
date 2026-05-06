	# Extract metadata from the original Seurat object to minimize the computational
	# resources required when only metadata is handled

obj <- readRDS(”otPools_all_integrated.RDS")

meta <- obj@meta.data

save(meta, file = "meta.RData")
