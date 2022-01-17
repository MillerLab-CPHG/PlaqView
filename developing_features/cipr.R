library(Seurat)
library(CIPR)


allmarkersfromcsv <- read.csv(file = "../DataProcessing/data/Litvinukova_2020/diff_by_author.csv")

# Plot summarizing top scoring references per cluster (logFC comparison)
CIPR(input_dat = allmarkersfromcsv,
     comp_method = "logfc_dot_product", 
     reference = "immgen", 
     plot_ind = F,
     plot_top = T)











