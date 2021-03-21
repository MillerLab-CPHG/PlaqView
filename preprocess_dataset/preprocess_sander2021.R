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
library(colourpicker)
library(rDGIdb)
library(ggrepel)


## seurat doesnt support GEO files, but can do it with a matrix
# note that normally geoGEO is for the matrix, but for NGS they are 
# placed in the supplementary files
# also note to UNZIP the files (it can come as .gz)

gset <- read.delim("GSE131778_DEcompressed_human_coronary_scRNAseq_wirka_et_al_GEO.txt",
                   row.names=1)

# construct the seurat object
stanford <- CreateSeuratObject(
  gset,
  project = "GSE131778_source",
  assay = "RNA",
  min.cells = 0,
  min.features = 0
)

# here you can load the file using the directory shown in the output in the previous line

#### SEURAT: QC and reduction---- 

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
stanford[["percent.mt"]] <- PercentageFeatureSet(stanford, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(stanford@meta.data, 5)
#Low-quality / dying cells often exhibit extensive mitochondrial contamination

# Visualize QC metrics as a violin plot
VlnPlot(stanford, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(stanford, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(stanford, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

stanford <- subset(stanford, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
stanford <- NormalizeData(stanford, normalization.method = "LogNormalize", scale.factor = 10000)
stanford <- FindVariableFeatures(stanford, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top20 <- head(VariableFeatures(stanford), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(stanford)
plot1 
LabelPoints(plot = plot1, points = top20, repel = TRUE)

all.genes <- rownames(stanford)
# scaling the data enables the subsequent PCAs
stanford <- ScaleData(stanford, features = all.genes)

# PCA; finding the number of PCA axis for DS analysis
# the following are three separate ways to look at PCA axis
stanford <- RunPCA(stanford, features = VariableFeatures(object = stanford))
print(stanford[["pca"]], dims = 1:5, nfeatures = 5)

# dim shows the PCA#, you can increase it to show more if needed
# this is helpful to show which genes are involved in the respsective PCA
VizDimLoadings(stanford, dims = 1:2, reduction = "pca")


DimPlot(stanford, reduction = "pca")

# this can be helpful to determine which PCA to include DS
DimHeatmap(stanford, dims = 1:15, cells = 500, balanced = TRUE)

# Jacksaw procedure to determine PCA relevance
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
stanford <- JackStraw(stanford, num.replicate = 100)
stanford <- ScoreJackStraw(stanford, dims = 1:30)
JackStrawPlot(stanford, dims = 1:15)

# alternatively you can run ELbow plot
ElbowPlot(stanford) # this shows you can use up to 20

# now for the clustering
stanford <- FindNeighbors(stanford, dims = 1:10)
stanford <- FindClusters(stanford, resolution = 0.5)

# this looks at the cluster ID for the first 5 cells
# not very important
head(Idents(stanford), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)
# reticulate::py_install(packages = 'umap-learn')
stanford <- RunUMAP(stanford, dims = 1:30)
stanford <- RunTSNE(stanford, dims = 1:30)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(stanford, reduction = "umap", label = T)

saveRDS(stanford, file = "stanford_presingleR.rds")
#### SINGLE-R LABELING (post Seurat automated cell cluster annotation alternative) ----
# BiocManager::install("SingleR")
# here we are using Human Primary Cell Atlas design for blood
# https://bioconductor.org/packages/3.12/data/experiment/vignettes/celldex/inst/doc/userguide.html#2_General-purpose_references
hpca.se <- celldex::HumanPrimaryCellAtlasData() # build the reference
hpca.se

# now run the prediction using the reference
# singleR requires that it be in a 'singlecellexperiment' format
# they are workout agnostic

for_singleR_input <- GetAssayData(stanford)
pred.stanford <- SingleR(test = for_singleR_input, 
                         ref = hpca.se, 
                         label = hpca.se$label.main) # reference cell types
pred.stanford
# summarize distribution
table(pred.stanford$labels)

# to show annotation confidence map
plotScoreHeatmap(pred.stanford)

# to show # that are pruned due to low score
summary(is.na(pred.stanford$pruned.labels))

### to place the singleR predictions into Seurat as a sep unit ###
# seurat.obj[["SingleR.labels"]] <- singler.results$labels
stanford[["SingleR.labels"]] <- pred.stanford$labels # this nest under metadata

# Copy over the labels and pruned.labels (Note: any other column of the results could be used as well)
stanford$SingleR.pruned.calls <- pred.stanford$pruned.labels
stanford$SingleR.calls <- pred.stanford$labels

#### RECODE SINGLE-R NAMES ----
stanford@meta.data[["SingleR.calls"]] <- recode(stanford@meta.data[["SingleR.calls"]], Smooth_muscle_cells = "SMC")
stanford@meta.data[["SingleR.calls"]] <- recode(stanford@meta.data[["SingleR.calls"]], Endothelial_cells = "EC")
stanford@meta.data[["SingleR.calls"]] <- recode(stanford@meta.data[["SingleR.calls"]], NK_cell = "NK")
stanford@meta.data[["SingleR.calls"]] <- recode(stanford@meta.data[["SingleR.calls"]], Chondrocytes = "CH")
stanford@meta.data[["SingleR.calls"]] <- recode(stanford@meta.data[["SingleR.calls"]], Fibroblasts = "FB")
stanford@meta.data[["SingleR.calls"]] <- recode(stanford@meta.data[["SingleR.calls"]], Monocyte = "Mono")
stanford@meta.data[["SingleR.calls"]] <- recode(stanford@meta.data[["SingleR.calls"]], B_cell = "B_Cells")
stanford@meta.data[["SingleR.calls"]] <- recode(stanford@meta.data[["SingleR.calls"]], Macrophage = "Mø")
stanford@meta.data[["SingleR.calls"]] <- recode(stanford@meta.data[["SingleR.calls"]], Tissue_stem_cells = "SC")
stanford@meta.data[["SingleR.calls"]] <- recode(stanford@meta.data[["SingleR.calls"]], T_cells = "T_Cells")
stanford@meta.data[["SingleR.calls"]] <- recode(stanford@meta.data[["SingleR.calls"]], 'Pre-B_cell_CD34-' = "PreB_CD34-")
stanford@meta.data[["SingleR.calls"]] <- recode(stanford@meta.data[["SingleR.calls"]], 'Pro-B_cell_CD34+' = "ProB_CD34+")


stanford@meta.data[["SingleR.calls"]]


#### COLOR SCHEME (for repr oducible external plots) ----
# the following function creates 'color_list' that shows the default color scheme
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
color_list <- ggplotColours(n=18)
# or you can try to define your own palette then pass the argument 'cols =' in Dimplot
# install.packages("colourpicker") # run this package via 'tools'addins' in rstudio-
manual_color_list <-
  c("rosybrown2",
    "palevioletred1",
    "lemonchiffon3",
    "darkseagreen",
    "cadetblue1",
    
    "skyblue3",
    "cadetblue3",
    "lemonchiffon4",
    "darkseagreen2",
    "rosybrown3",
    "thistle2",
    "darkseagreen1",
    "palevioletred3",
    
    "lightsteelblue3",
    "cadetblue2",
    "salmon1",
    
    "palevioletred4",
    "thistle3"
  )

#### SEURAT DIM PLOTS #### 
# now you can call the dimplot

pdf("Figure_images/umap_automatic_annotation.pdf", width=8, height=6)
plot <- DimPlot(
  stanford,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  repel = T, # repel labels
  pt.size = 1,
  cols = manual_color_list,
  group.by = "SingleR.calls") +# group.by is important, use this to call metadata separation
  ggtitle(NULL)
plot
dev.off()


pdf("Figure_images/umap_wirka_annotation.pdf", width=8, height=6)
plot <- DimPlot(
  stanford,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  repel = T, # repel labels
  pt.size = 1,
  cols = manual_color_list,
  group.by = "manually_annotated_labels") +# group.by is important, use this to call metadata separation
  ggtitle(NULL)
plot
dev.off()

# # if you want to subset some cells
# # then un-comment these codes
# CellSelector(plot)
#### POPULATION STATISTICS ####
cellcounts <- (stanford@meta.data[["SingleR.calls"]]) # extract cell names
cellcounts <- as.factor(cellcounts) # convert to factors

cellcounts.summary <- as.data.frame(summary(cellcounts)) # for levels laters

cellcounts <- as.data.frame((cellcounts)) # summarize factors into table
cellcounts <- dplyr::rename(cellcounts, Celltype = '(cellcounts)' ) # just to rename the columns

# this is to create levels so i can flip the graph to the way i want
cellcounts.summary <- rownames_to_column(cellcounts.summary)
cellcounts.summary <- reorder(cellcounts.summary$rowname, cellcounts.summary$`summary(cellcounts)`)
sortedlevel <- levels(cellcounts.summary)

pdf("Figure_images/umap_populationpercentage.pdf", width=6, height=4)
ggplot(data = cellcounts, aes(y = Celltype)) +
  geom_bar(fill = manual_color_list) +
  xlab("Counts") +
  ylab("Cell Types") +
  xlim(c(0,1650)) +
  theme_light() +
  scale_y_discrete(limits=sortedlevel) +
  stat_count(geom = "text", # this stat_count function gives percentages
           aes(label = paste(round((..count..)/sum(..count..)*100,2),"%")))
dev.off()


#### SINGLE-R DIAGNOSTICS diagnostic plots ----
plotScoreHeatmap(pred.stanford) #inspect the confidence of the predicted labels 
plotScoreDistribution(pred.stanford, ncol = 3)


#### MONOCLE3 TRAJECTORY INFERENCE ----
# in previous versions we tried the seurat wrapper it just didnt work
# below we manually wrap the data ourselves

# convert to monocle cds object 
# Extract data, phenotype data, and feature data from the SeuratObject
expressiondata <- stanford@assays[["RNA"]]@data

cellmd <- stanford@meta.data

genemd <- data.frame(gene_short_name = row.names(expressiondata), 
                     row.names = row.names(expressiondata))

# Construct monocle cds
stanford.cds <- new_cell_data_set(expression_data = expressiondata,
                              cell_metadata = cellmd,
                              gene_metadata = genemd)
stanford.cds <- preprocess_cds(stanford.cds, num_dim = 30) # we used 30 in earlier seurat scripts

# 
# run clustering again (didnt transfer from seurat)
stanford.cds <- reduce_dimension(stanford.cds, reduction_method = "UMAP")
stanford.cds <- cluster_cells(stanford.cds, reduction_method = "UMAP")


#### TRANSFER SEURAT EMBEDDINGS #####
# Note that these may be calculated on the Integrated object, not the counts
#   and thus will involve fewer genes
temp.cds <- ProjectDim(stanford, reduction = "pca") # this will be removed
reducedDim(stanford.cds, type = "PCA") <- temp.cds@reductions$pca@cell.embeddings
stanford.cds@preprocess_aux$prop_var_expl <- temp.cds@reductions$pca@stdev
plot_pc_variance_explained(stanford.cds)

# Transfer Seurat UMAP embeddings
stanford.cds@int_colData@listData$reducedDims$UMAP <- temp.cds@reductions$umap@cell.embeddings

## transfer singleR labels to moncle3 object
colData(stanford.cds)$assigned_cell_type <- stanford@meta.data[["SingleR.calls"]] # call this by opening the object

#### MONOCLE3 CONT. ----
# now learn the PATH (trajectory)
stanford.cds <- learn_graph(stanford.cds)

# this calls up a shiny app, choose the ROOT NODE
stanford.cds <- order_cells(stanford.cds, reduction_method = "UMAP")

# finally, you can visualize the learned path
pdf("Figure_images/monocle3_RNAvelocity_seuratpartition.pdf", width=6, height=6)
plot_cells(stanford.cds,
           color_cells_by = "assigned_cell_type",
           label_groups_by_cluster=F,
           show_trajectory_graph = T,
           trajectory_graph_segment_size = 1,
           label_leaves=F, # this gives a little node label (outcome)
           label_roots = T,
           label_branch_points = F,
           graph_label_size = 1, # size of # in circle
           group_label_size = 3,
           cell_size = 1,
           alpha = 0.7,
           scale_to_range = T) +
  scale_color_manual(values = manual_color_list) # sync color scheme
dev.off()

# now you can show pseudotime
pdf("Figure_images/monocle3_pseudotime_seuratpartition.pdf", width=7, height=6)
plot_cells(stanford.cds,
           color_cells_by = "pseudotime",
           show_trajectory_graph = F,
           trajectory_graph_segment_size = 1,
           label_leaves=F, # this gives a little node label (outcome)
           label_roots = T,
           label_branch_points = F,
           graph_label_size = 1, # size of # in circle
           group_label_size = 3,
           cell_size = 1,
           alpha = 0.7,
           scale_to_range = T) 
dev.off()


#### Subset Trajectory & analysis of SMC----
stanford.cds_subset <- choose_cells(stanford.cds) # calls up shiny app

plot_cells(stanford.cds_subset,
           color_cells_by = "pseudotime",
           show_trajectory_graph = T,
           trajectory_graph_segment_size = 1,
           label_leaves=F, # this gives a little node label (outcome)
           label_roots = T,
           label_branch_points = F,
           graph_label_size = 1, # size of # in circle
           group_label_size = 3,
           cell_size = 1,
           alpha = 0.7,
           scale_to_range = T) 

#### MORAN's I Test of Autocorrelation ####
# now we can extrapolate genes that are differentially expressed in this region
# Moran’s I is a measure of multi-directional and multi-dimensional spatial autocorrelation. 
# the statistic tells you whether cells at nearby positions on a 
# trajectory will have similar (or dissimilar) +
# expression levels for the gene being tested.
## first lets do the whole dataset
# a special gene module score heatmap (for the whole dataset)
# pr_graph_test_res <- graph_test(stanford.cds, neighbor_graph="principal_graph", cores=2)
write.csv(pr_graph_test_res, file = "moransI_all_clusters.csv")
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.00000001)) # you can adjust the p-value here
head(pr_deg_ids)
gene_module_df <- find_gene_modules(stanford.cds[pr_deg_ids,], resolution=1e-3)
cell_group_df <- tibble::tibble(cell=row.names(colData(stanford.cds)), 
                                cell_group=colData(stanford.cds)$assigned_cell_type)
agg_mat <- aggregate_gene_expression(stanford.cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

# which then can be visualized like so;
# this can show you the different gene modules that can are responsible for changes over pseudotime
plot_cells(stanford.cds,
           genes=gene_module_df %>% filter(module %in% c(2,3,7)), # specify the module you want to examine
           label_cell_groups=T,
           show_trajectory_graph=F)

subset(gene_module_df, module == 2)

## now lets do the subsets
# pr_graph_test_res.sub <- graph_test(stanford.cds_subset, neighbor_graph="principal_graph", cores=2)
pr_deg_ids.sub <- row.names(subset(pr_graph_test_res.sub, q_value < 0.00000001))
write.csv(pr_graph_test_res.sub, file = "moransI_subset_cluster.csv")
head(pr_deg_ids.sub)

# collect the trajectory-variable genes into modules
gene_module_df.sub <- find_gene_modules(stanford.cds_subset[pr_deg_ids.sub,], resolution=1e-3)
# visualize these genes
# here I am just pulling out genes that have high moran's i and might be helpful in the paper
# SELECTED FOR PUBLICATIONS
pdf("Figure_images/monocle3_genesoverpseudotime_seuratpartition_extended.pdf", width=7, height=6)
plot_cells(stanford.cds_subset, 
           genes=c("MYH11", 'IGFBP2',"PPP1R14A","CNN1", "TNFRSF11B",
                   "C7", "C3",
                   "SERPINF1",  "FBLN1", 
                   "CXCL12", "MMP2", 
                   "FN1"), # this is faceting by the genes that are DE
           show_trajectory_graph=FALSE, 
           label_cell_groups=F, cell_size = 1)

dev.off()

# recluster at higher definition
stanford.cds_subset = cluster_cells(stanford.cds_subset, resolution=1e-2)

pdf("Figure_images/monocle3_RNAvelocitySUBSET_seuratpartition.pdf", width=6, height=6)
plot_cells(stanford.cds_subset, 
           color_cells_by="cluster",
           label_groups_by_cluster=F,
           show_trajectory_graph = T,
           trajectory_graph_segment_size = 1,
           label_leaves=F, # this gives a little node label (outcome)
           label_roots = F,
           label_branch_points = F,
           graph_label_size = 1, # size of # in circle
           group_label_size = 4,
           cell_size = 1,
           alpha = 0.5,
           scale_to_range = T)
dev.off()

#### SPECIAL LOOK UP ----
# use this and compare to the annotated map
# this you can use to more accurately differentiate the CELLS expression
# this is how to use singleR imported cell identity
diffexSMC_chondro_BYCELL <- FindMarkers(stanford, group.by = "SingleR.calls",
                                 ident.1 = "Smooth_muscle_cells",
                                 ident.2 = "Fibroblasts")

write.table(diffexSMC_chondro_BYCELL, file = "SMCvsChondro_BYCELL.csv", sep = ",", 
            col.names = NA) # col.names = NA sets the first cell to be empty so its offset correctly

# generate plots for genes of interests
features <- c("GATA6", "CTAGE1", "PLCE1", "DUSP13", "SAMD8")
RidgePlot(stanford, features = features, ncol = 2, group.by = "SingleR.calls", log = TRUE, y.max = 20)
VlnPlot(stanford, features = features, ncol = 2, group.by = "SingleR.calls", log = TRUE)
FeaturePlot(stanford, features = features)
DotPlot(stanford, features = features, group.by = "SingleR.calls")



#### MORE LOOK UP ----
# first set active idents
Idents(stanford) <- stanford@meta.data[["SingleR.calls"]]
markers_SMC <- FindMarkers(stanford, ident.1 = "SMC") # since we didnt specific ident.2, this shows relative to all other cell types

write.csv(markers_SMC, file = "SMC_specific_DEGs.csv")

features <- c("MYH11", "FN1",  "COL6A1", 'COL6A2', "PPP1R14A",
              "TNFRSF11B", #TOP HEAVY
              "FBLN1","LUM", "TCF21", "C7", "C6", # BOTTOM HEAVY
              "SERPINF1" )

RidgePlot(stanford, features = features, ncol = 1, group.by = "SingleR.calls", log = TRUE, y.max = 20)
VlnPlot(stanford, features = features, ncol = 1, group.by = "SingleR.calls", log = TRUE)
FeaturePlot(stanford, features = features)
DotPlot(stanford, features = features, group.by = "SingleR.calls")

pdf("Figure_images/monocle3_genesoverpseudotime_seuratpartition_extended_v2.pdf", width=7, height=6)
plot_cells(stanford.cds_subset, 
           genes=features, # this is faceting by the genes that are DE
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)

dev.off()

# 03192021
CAC_genes_for_Wei <- read_csv("~/Desktop/CAC_genes_for_Wei.csv", 
                             col_names = FALSE)

features <- CAC_genes_for_Wei$X1
RidgePlot(stanford, features = features, ncol = 4, group.by = "SingleR.calls", log = TRUE, y.max = 20)
VlnPlot(stanford, features = features, ncol = 1, group.by = "SingleR.calls", log = TRUE)
FeaturePlot(stanford, features = features, ncol = 4)
DotPlot(stanford, features = features, group.by = "SingleR.calls") +
  theme(axis.text.x = element_text(angle = 90))

Idents(stanford) <- stanford@meta.data[["manually_annotated_labels"]]
DotPlot(stanford, features = features) +
  theme(axis.text.x = element_text(angle = 90))

#### STACKED POPULATION PLOT ####
Idents(stanford) <- stanford@meta.data[["SingleR.calls"]]
pop1 <- as.data.frame(prop.table(table(Idents(stanford))))
pop1$Method <- "SingleR"

Idents(stanford) <- stanford@meta.data[["manually_annotated_labels"]]
pop2 <- as.data.frame(prop.table(table(Idents(stanford))))
pop2$Method <- "ManualCluster"

plot1 <- ggplot(pop1, aes(y = Freq, x =Method,)) + 
  geom_bar(position="stack", stat="identity",
           fill = manual_color_list,
           width = 0.1) + 
  xlab("") +
  ylab("Cell Porportions") +
  geom_label_repel(aes(label = Var1), colour = "black",
            position = position_stack(vjust = 0.5),
            max.overlaps = 15 ,
            force = 10,
            force_pull = 3,
            max.iter = 999999,
            max.time = 1,
            min.segment.length = 0.01,
            xlim = c(1, 3)) +
  theme(axis.text.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 
  
ggsave(plot1, file = "verticle_cellpop_plot_singleR.pdf",
       width = 5, height = 7)


plot2 <- ggplot(pop2, aes(y = Freq, x =Method,)) + 
  geom_bar(position="stack", stat="identity",
           fill = manual_color_list[1:14],
           width = 0.1) + 
  xlab("") +
  ylab("Cell Porportions") +
  geom_label_repel(aes(label = Var1), colour = "black",
                   position = position_stack(vjust = 0.5),
                   max.overlaps = 15 ,
                   force = 10,
                   force_pull = 3,
                   max.iter = 999999,
                   max.time = 1,
                   min.segment.length = 0.001,
                   xlim = c(1, 3)) +
  theme(axis.text.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 


ggsave(plot2, file = "verticle_cellpop_plot_manual.pdf",
       width = 5, height = 7)



#### OUTPUT RDS FILES ####
saveRDS(stanford, file = "final_stanford_labeled.rds")
stanford <- readRDS(file = "final_stanford_labeled.rds")

 saveRDS(stanford.cds, file = "final_stanford_labeled_CDS.rds")

#### FINAL OUTPUTS ####
# this extended label rds has the wirka labels too
stanford <- readRDS(file = "final_stanford_extendedlabels_02082021.rds")
stanford.cds <- readRDS(file = "final_stanford_labeled_CDS.rds")

