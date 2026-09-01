library(Seurat)
library(Signac)
library(scCustomize)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(tidydr)
library(Matrix)
library(ggplot2)
setwd('/data/bio_data/')
####### read and process public data (Chengxiang Qiu et al, Systematic reconstruction of cellular trajectories across mouse embryogenesis) #######
data1 <- readRDS('embryo/rds/seurat_object_E13.5.rds')
data2 <- readRDS('embryo/rds/seurat_object_E12.5.rds')
data3 <- readRDS('embryo/rds/seurat_object_E11.5.rds')
data4 <- readRDS('embryo/rds/seurat_object_E10.5.rds')
data5 <- readRDS('embryo/rds/seurat_object_E9.5.rds')
data6 <- readRDS('embryo/rds/seurat_object_E8.5a.rds')
data7 <- readRDS('embryo/rds/seurat_object_E8.5b.rds')
data = merge(x= data1,y=c(data2,data3,data4,data5,data6,data7))
# modify gene names
library(org.Mm.eg.db)
# Ensembl ID
ensembl_ids <- rownames(data)
# Ensembl -> Symbol
gene_names <- mapIds(org.Mm.eg.db,keys = ensembl_ids,column = "SYMBOL",keytype = "ENSEMBL",multiVals = "first")
# 
id_map <- data.frame(
  ensembl = names(gene_names),
  symbol = as.character(gene_names),
  stringsAsFactors = FALSE)
id_map$final_name <- ifelse(is.na(id_map$symbol), id_map$ensembl, id_map$symbol)
id_map$final_name <- make.unique(id_map$final_name)
# extract Counts, Data  
counts <- GetAssayData(data, slot = "counts")
rownames(counts) <- id_map$final_name
# new seurat：
new_assay <- CreateAssayObject(counts = counts)
data[["RNA"]] <- new_assay
data <- NormalizeData(data)
data <- FindVariableFeatures(data,nfeatures =3000)
data <- ScaleData(data)
data <- RunPCA(data)
data <- FindNeighbors(data, dims = 1:30)
data <- FindClusters(data)
data <- RunUMAP(data, dims = 1:30)
var_int <- "sample"
data <- RunHarmony(data, group.by.vars = var_int, assay.use="RNA", reduction.use = "pca", plot_convergence = FALSE)
red_use <- "harmony"
dims_use <- 1:30
nneigh <- 30
data <- RunUMAP(object = data, reduction = red_use, dims = dims_use,n.neighbors = nneigh)
data <- FindNeighbors(data, verbose = T, reduction = red_use, dims = dims_use)
data <- FindClusters(data, verbose = T, algorithm = 1, resolution = c(0.3,0.5,0.8))

################# label transfer ###################
# read our datasets
embryo <- readRDS('embryo/Membryo_rna_atac.rds')
DefaultAssay(embryo) <- "RNA"
# transfer public annotation to our datasets
anchor <- FindTransferAnchors(
  reference = data,query = embryo,
  reduction = "rpca",dims = 1:30,k.anchor=30)
predicted.labels <- TransferData(anchorset = anchor,refdata = data$cell_type,
  weight.reduction = "pcaproject", # rpca.ref, rpca.query
  #dims = 1:30
)
embryo <- AddMetaData(object = embryo, metadata = predicted.labels)
# compare public celltype annotation to ours
df <- as.data.frame(table(embryo$predicted.id,embryo$L2_anno))
colnames(df) <- c("Predicted", "L2_Anno", "Count")
# reorder
col_order <- df %>% group_by(Predicted) %>%
  slice_max(order_by = Count, n = 1, with_ties = FALSE) %>% ungroup() %>%
  arrange(Predicted) %>% pull(L2_Anno) %>% as.character()
all_cols <- levels(df$L2_Anno)
remaining_cols <- setdiff(all_cols, col_order)
final_col_levels <- unique(c(col_order, remaining_cols))
df$L2_Anno <- factor(df$L2_Anno, levels = final_col_levels)
# plot
df_plot <- df %>%
  group_by(Predicted) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()
options(repr.plot.width=12, repr.plot.height=15)
p <- ggplot(df_plot, aes(x = L2_Anno, y = Predicted, fill = Proportion)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#e31a1c", name = "Proportion") + theme_minimal() +
  labs(title = "Mapping Confusion Matrix (Ordered Columns)", x = "L2 Annotation", y = "Predicted ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
  geom_text(aes(label = ifelse(Count > 0, Count, "")), size = 3, 
            color = ifelse(df_plot$Proportion > 0.5, "white", "black"))
