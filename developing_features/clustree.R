library(clustree)
# install.packages("remotes")
# remotes::install_github("lazappi/clustree@develop", dependencies = TRUE,
#                        build_vignettes = TRUE)

plaqviewobj <- readRDS(file = "data/Alencar_2020/source_files/UNPROCESSED.rds")

#### STEP2: SEURAT PROCESS ####
# Run the standard workflow for visualization and clustering
plaqviewobj <- NormalizeData(plaqviewobj)
plaqviewobj <- FindVariableFeatures(plaqviewobj, verbose = T, nfeatures = 2000)
plaqviewobj <- ScaleData(plaqviewobj, verbose = T)

plaqviewobj <- RunPCA(plaqviewobj, npcs = 30, verbose = FALSE)
plaqviewobj <- RunUMAP(plaqviewobj, reduction = "pca", dims = 1:20)
plaqviewobj <- RunTSNE(plaqviewobj, reduction = "pca", dims = 1:20)
plaqviewobj <- FindNeighbors(plaqviewobj, reduction = "pca", dims = 1:20)

plaqviewobj <- FindClusters(plaqviewobj, resolution = 0.1)
plaqviewobj <- FindClusters(plaqviewobj, resolution = 0.2)
plaqviewobj <- FindClusters(plaqviewobj, resolution = 0.3)
plaqviewobj <- FindClusters(plaqviewobj, resolution = 0.4)

Seurat::Idents(plaqviewobj) <- plaqviewobj$Author_Provided

clustree(plaqviewobj)
clustree_overlay(plaqviewobj, x_value = "UMAP_1", y_value = "UMAP_2")

