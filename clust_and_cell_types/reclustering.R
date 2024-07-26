.libPaths(c("/projappl/project_2010414/project_rpackages_4.4.0", .libPaths()))

library(Seurat)
library(dplyr)

set.seed(27)

obj <- readRDS("/scratch/project_2010414/ot/saved/otPools_all_integrated_v3_annotNour.RDS")

obj <- FindNeighbors(obj, dims = 1:14)

obj <- FindClusters(
	obj,
	resolution = 0.2,
	method = "matrix",
	algorithm = 1,
	n.start = 10,
	n.iter = 10,
	group.singletons = TRUE,
	temp.file.location = "/scratch/project_2010414",
	verbose = TRUE)

obj <- RunUMAP(obj, dims = 0:12)

save(obj, file=paste0("/scratch/project_2010414/ot/saved/Reclustering/reclustered_obj.RData"))
