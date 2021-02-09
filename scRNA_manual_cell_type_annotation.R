library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)




#  Helper function to plot genes expression in UMAP dims
plotgenesUMAP = function(object, genes, treatmentID=NULL) {
  if (is.null(treatmentID)) {
    FeaturePlot(object = object, features = genes, reduction = "umap", 
                sort.cell = TRUE, label = TRUE, label.size = 3, pt.size = 1)
  } else {
    FeaturePlot(object = object, features = genes, reduction = "umap", split.by = treatmentID, 
                sort.cell = TRUE, label = TRUE, label.size = 3, pt.size = 1)
  }
}

###########################
# Read Wei's seurat object
mod_smc_obj = readRDS("final_stanford_labeled.rds")
head(mod_smc_obj@meta.data)

# Plot clusters to get an idea of how the clusters are labeled
rna_clusters = DimPlot(mod_smc_obj, reduction = "umap", group.by = "seurat_clusters",
                       pt.size = 0.1,label.size = 5, label = TRUE)


# Find cluster markers
Idents(mod_smc_obj) = mod_smc_obj@meta.data$seurat_clusters
cluster_markers = FindAllMarkers(mod_smc_obj)

# Explore marker genes in each cluster
c0_genes = cluster_markers[cluster_markers$cluster=="0", ] # These look like Endothelial cells
head(c0_genes, n=20)

c1_genes = cluster_markers[cluster_markers$cluster=="1", ] # These look like fibroblasts
head(c1_genes, n=20)

c2_genes = cluster_markers[cluster_markers$cluster=="2", ] # These are definitely Macrophages
head(c2_genes, n=20)

c3_genes = cluster_markers[cluster_markers$cluster=="3", ] # These look like T cells
head(c3_genes, n=20)

c4_genes = cluster_markers[cluster_markers$cluster=="4", ] # These look like Fibromyocytes
head(c4_genes, n=30)

c5_genes = cluster_markers[cluster_markers$cluster=="5", ] # These look like Pericytes1
head(c5_genes, n=30)

c6_genes = cluster_markers[cluster_markers$cluster=="6", ] # These look like SMC
head(c6_genes, n=30)

c7_genes = cluster_markers[cluster_markers$cluster=="7", ] # These look like B cells
head(c7_genes, n=30)

c8_genes = cluster_markers[cluster_markers$cluster=="8", ] # These look like Plasma cells1
head(c8_genes, n=30)

c9_genes = cluster_markers[cluster_markers$cluster=="9", ] # These look like Pericytes2
head(c9_genes, n=30)

c10_genes = cluster_markers[cluster_markers$cluster=="10", ] # Unknown according to supplementary material (cluster11): These look like monocytes imo
head(c10_genes, n=30)

c11_genes = cluster_markers[cluster_markers$cluster=="11", ] # Unknown according to supplmentary material (cluster18), set it as unknown for now
head(c11_genes, n=30)

c12_genes = cluster_markers[cluster_markers$cluster=="12", ] # These look like Neurons
head(c12_genes, n=30)

c13_genes = cluster_markers[cluster_markers$cluster=="13", ] # These look like Mast cells
head(c13_genes, n=30)


# Annotate clusters
manually_annotated_clusters_ids = c("Endothelial", "Fibroblast", "Macrophage", "T cell", "Fibromyocyte", "Pericyte1", "SMC",
                                    "B cell", "Plasma cell", "Pericyte2", "Monocyte", "Unknown", "Neuron", "Mast cell")
names(manually_annotated_clusters_ids) = levels(mod_smc_obj)
mod_smc_obj[["manually_annotated_labels"]] = manually_annotated_clusters_ids[match(mod_smc_obj@meta.data$seurat_clusters,
                                                                                   names(manually_annotated_clusters_ids))]

DimPlot(mod_smc_obj, group.by = "manually_annotated_labels" , repel=TRUE,
        label.size=4, pt.size = 0.1,
        label = TRUE) + ggtitle("scRNA Athero lesion dataset manually assigned cell types")


