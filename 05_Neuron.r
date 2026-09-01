library(Seurat)
library(Signac)
library(scCustomize)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(Matrix)
library(tidydr)
setwd('/data/bio_data/')
set.seed(123)
embryo <- readRDS('embryo/Membryo_rna_atac.rds')
embryo_neuron <- embryo[,embryo$L1_anno == 'Neuroectoderm']
Process_RNA <- function(proj=proj, nfeatures=2000, pc.nums=50, dims=1:50, pca.name = 'pca', umap.name= 'umap.rna', res = 0.8){
  DefaultAssay(proj) <- "RNA"
  #
  proj <- NormalizeData(proj)
  proj <- FindVariableFeatures(proj, nfeatures = nfeatures)
  # VariableFeatures(proj) <- VariableFeatures(proj)[!VariableFeatures(proj) %in% blacklist.genes]
  proj <- ScaleData(proj)
  proj <- RunPCA(proj, npcs = pc.nums, reduction.name = pca.name)
  #
  return(proj)
}
embryo_neuron <- Process_RNA(embryo_neuron,nfeatures=3000,pc.nums=50, dims=1:50, pca.name = 'pca', umap.name= 'umap.rna', res = c(0.3,0.5,0.8))
var_int <- "sample"
embryo_neuron <- RunHarmony(embryo_neuron, group.by.vars = var_int, assay.use="RNA", reduction.use = "pca", plot_convergence = FALSE)
red_use <- "harmony"
dims_use <- 1:50
nneigh <- 50
embryo_neuron <- RunUMAP(object = embryo_neuron, reduction = red_use, dims = dims_use,n.neighbors = nneigh) # ,min.dist=0.2
embryo_neuron <- FindNeighbors(embryo_neuron, verbose = T, reduction = red_use, dims = dims_use)
embryo_neuron <- FindClusters(embryo_neuron, verbose = T, algorithm = 1, resolution = c(0.5,1.0,1.2,1.5))

# track plot
DefaultAssay(embryo_neuron) <- "ATAC"
Idents(embryo_neuron) <- 'Sample'
DefaultAssay(embryo_neuron) <- 'ATAC'
# cluster_cols = c('#6778AE','#F7F398','#91D0BE','#B2DF8A','#58A4C3',   '#aa8282','#666666','#476D87','#55A1B1','#E39A35','#DCC1DD','#E5D2DD','#A6761D','#CCE0F5','#BEAED4')
options(repr.plot.width=14, repr.plot.height=14)
cov_plot <- CoveragePlot(object = embryo_neuron, region = 'Dcx', features = "Dcx", annotation = TRUE, peaks = TRUE, links = FALSE)  & scale_fill_manual(values = sample_cols) 
cov_plot