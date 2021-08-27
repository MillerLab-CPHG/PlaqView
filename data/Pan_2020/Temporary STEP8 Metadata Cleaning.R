#### Library and data loading ----
library(Seurat)
library(patchwork)
library(readr)
library(scCATCH)
library(SingleR)
library(tidyverse)
library(monocle3)
library(SeuratData)
library(magrittr)
library(ggrepel)
library(dyno)
library(readxl)
library(SeuratDisk)

original_color_list <-
  {c("rosybrown2",
     "cadetblue1",
     "lemonchiffon3",
     "darkseagreen",
     "skyblue3",
     "thistle3",
     "cadetblue3",
     "darkseagreen1",
     "palevioletred3",
     "palevioletred1",
     "darkseagreen2",
     "rosybrown3",
     "thistle2",
     "lightsteelblue3",
     "salmon1",
     "palevioletred4",
     "lemonchiffon4",
     "cadetblue2"
  )}

color_function <- colorRampPalette(original_color_list)
manual_color_list <- color_function(40) # change this if clusters >40




#### Clean-Up Metadata ####

plaqviewobj <- readRDS("Pan_2020.rds")

# show all metadata columns
names(plaqviewobj@meta.data)

# last chance to rename any 
plaqviewobj@meta.data$Seurat_Clusters <- plaqviewobj@meta.data$seurat_clusters
plaqviewobj@meta.data$Author_Provided <- plaqviewobj@meta.data$manually_annotated_labels
plaqviewobj@meta.data$Seurat_with_Tabula_Ref <- plaqviewobj@meta.data$predicted.id_tabulus.sapien

# choose which ones to keep for display
plaqviewobj@meta.data <- 
  plaqviewobj@meta.data[, which(colnames(plaqviewobj@meta.data)
                                %in% c(
                                  "Seurat_Clusters",
                                  "scCATCH_Blood",
                                  "scCATCH_BV",
                                  "scCATCH_Heart",
                                  "Author_Provided",
                                  "SingleR.calls",
                                  "Seurat_with_Tabula_Ref"  
                                ))]
plaqviewobj <- DietSeurat(plaqviewobj, counts = T, data = T, dimreducs = c('umap'))
summary(table(plaqviewobj$Author_Provided))


#### output ####
saveRDS(plaqviewobj, file = "renamed_metadata.RDS")
