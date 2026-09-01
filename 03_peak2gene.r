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
library(TFBSTools)
library(presto)
library(motifmatchr)
library(qlcMatrix)
library(RColorBrewer)
library(ggrepel)
library(JASPAR2020)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(chromVAR)
setwd('/data/bio_data/')
options(future.globals.maxSize = 1000 * 1024^3)
set.seed(1234)
embryo <- readRDS('embryo/Membryo_rna_atac.rds')
# SCT
DefaultAssay(embryo) <- "RNA"
embryo <- SCTransform(embryo, vars.to.regress = "percent.mt", verbose = FALSE)
# chromVAR
DefaultAssay(embryo) <- "ATAC"
pfm <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
embryo <- AddMotifs(object = embryo, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)
library(BiocParallel)
# register(SerialParam())
register(MulticoreParam(6))
embryo <- RunChromVAR(object = embryo, genome = BSgenome.Mmusculus.UCSC.mm10)

######## peak to gene ###############
# downsample to 2000 cells / celltype
sample <- unlist(lapply(unique(embryo$stage),function(x){ 
        cell <- colnames(embryo)[embryo$stage==x]
        if(length(cell)>2000){sample(cell,2000)} else {cell}
}))
embryo <- embryo[,sample]
# Find DEGs bw clusters ####
DefaultAssay(embryo) <- "RNA"
Idents(embryo) <- embryo$stage
de_genes <- FindAllMarkers(embryo,only.pos = T,logfc.threshold = 0.1,max.cells.per.ident = 1000)
deg <- unique(de_genes[de_genes$p_val_adj<0.01 & de_genes$avg_log2FC>0.25,'gene'])
# Linking peaks to genes ####
DefaultAssay(embryo) <- "ATAC"
embryo <- RegionStats(embryo, genome = BSgenome.Mmusculus.UCSC.mm10)
deg <- unique(de_genes[de_genes$p_val_adj<0.05 & de_genes$avg_log2FC>0.1,'gene'])
embryo <- LinkPeaks(embryo,peak.assay = "ATAC",expression.assay = "SCT",genes.use = deg) # 非常耗时 投递处理
link <- Links(embryo)
link <- data.frame(link)
link$adj_pval <- p.adjust(link$pvalue,method = 'BH')
link$gene_cluster <- de_genes[match(link$gene,de_genes$gene),'cluster']
lvl <- levels(embryo)
link <- arrange(link,factor(gene_cluster,levels=lvl))
link <- link[link$adj_pval<0.05 & link$score>0,]
t <- link %>% group_by(gene) %>% summarise(n=n())
summary(t$n)
link_processed <- link[link$zscore>0.1,]
link_processed <- arrange(link_processed,factor(gene_cluster,levels=lvl))
# figure 2E
## gene expression and chromatin accessibility for significantly linked CRE-gene pairs
library(zoo)
DefaultAssay(embryo) <- "RNA"
gene.use <- unique(link$gene)
# caculate cellcycle score
cc.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
embryo <- CellCycleScoring(embryo, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) 
embryo <- ScaleData(embryo,vars.to.regress = c("S.Score", "G2M.Score"),features=gene.use)
rna <- embryo@assays$RNA['scale.data'][gene.use,]
asplit_cells <- split(colnames(embryo), embryo$Sample)
n <- 50
means <- do.call(cbind, lapply(lvl, function(x){
        df <- rna[,asplit_cells[[x]]]
        t(apply(df,1,function(x){rollapply(x,n,mean,by=n)}))
}))
celltype_major <- unlist(lapply(lvl, function(x) rep(x, length(asplit_cells[[x]])%/%n)))
## plot RNA 
anno_col <- data.frame(celltype_major)
rownames(anno_col) <- colnames(means) <- paste0(seq(1:ncol(means)),colnames(means))
anno_col <- arrange(anno_col,factor(celltype_major,levels=lvl))
anno_col$celltype_major <- factor(anno_col$celltype_major,levels=lvl)
means <- means[link$gene,rownames(anno_col)]
genes_mark <- list(E8  = c('Myh6','Myocd','Ttn','Ankrd1','T','Cdx4','Fst','Sfrp5','Slc2a3','Slc16a3','Car4','Prtg','Emb','Itgb6'),
  E9  = c('Vwf','Apob','Apoa1','Cubn','Crabp2','Scd1','Msmo1','Grik1'),
  E10 = c('Afp','Ttr','Apoa1','Hes1','Ccn1','Hbb-y'),
  E11 = c('Hbb-bs','Hbb-bt','Hemgn','Robo3','Slc17a6','Slc6a5','Itgam','Il1r2','Ucp1','Tktl1','Ctse','Ampd1'),
  E12 = c('Neurod2','Gad2','Gria2','Myt1l','Dlx6os1','Pcdh9','Ncam2','Nxph1','Slc17a6','Fstl5','Acan','Lncppara'))
genes_mark <- unique(unlist(genes_mark))
genes_mark_df <- as.data.frame(genes_mark)
# 
first_occurrence <- sapply(genes_mark, function(x) which(rownames(means) == x)[1])
first_occurrence
means[link$gene,rownames(anno_col)]
options(repr.plot.width = 7, repr.plot.height = 6)
# yellow2blue = colorRampPalette(colors = c("blue","black","yellow"))(100)
library(circlize)
BluYl = colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B"))
sample_cols <- c("#2C7BB6","#ABD9E9","#FFFFBF","#FDAE61","#D7191C")
label_colors <- list(celltype_major = setNames(sample_cols, levels(embryo$Sample)))
p1 = pheatmap(means, cluster_rows = F, cluster_cols = F, scale = "row", show_rownames = F,
         breaks = seq(-2, 2, length = 101), col = BluYl,
         annotation_col = anno_col,show_colnames = F,use_raster = TRUE, annotation_colors = label_colors, 
         main='Gene expression') + rowAnnotation(link = anno_mark(at = first_occurrence, labels = genes_mark, labels_gp = gpar(fontsize = 8))) 
print(p1)
dev.off()
## plot ATAC
peak.use <- unique(link$peak)
DefaultAssay(embryo) <- "ATAC"
Idents(embryo) <- 'Sample'
atac <- RunTFIDF(embryo@assays$ATAC@counts[peak.use,])
asplit_cells <- split(colnames(embryo), embryo$Sample)
means_peak <- do.call(cbind, lapply(lvl, function(x){
        df <- atac[,asplit_cells[[x]]]
        t(apply(df,1,function(x){rollapply(x,n,mean,by=n)}))
}))
celltype_major <- unlist(lapply(lvl, function(x) rep(x, length(asplit_cells[[x]])%/%n)))
anno_col2 <- data.frame(celltype_major)
rownames(anno_col2) <- colnames(means_peak) <- paste0(seq(1:ncol(means_peak)),colnames(means_peak))
anno_col2$celltype_major <- factor(anno_col2$celltype_major,levels=lvl)
anno_col2 <- arrange(anno_col2,factor(celltype_major,levels=lvl))
options(repr.plot.width = 6, repr.plot.height = 6)
p2 = pheatmap(means_peak[link_processed$peak,rownames(anno_col2)],cluster_rows = F, cluster_cols = F, scale = "row", show_rownames = F,
               breaks = seq(-2, 2, length = 101), color = colorRamp2(c(-2, 0, 2), c("#440154", "#21908C", "#FDE725")), # colorRamp2(c(0, 1, 3), c("#ECEFF1", "#FB8C00", "#B71C1C"))
               annotation_col = anno_col, show_colnames = F, use_raster=TRUE, annotation_colors = label_colors, 
               main='Chromatin accessibility') 
print(p2)
dev.off()

