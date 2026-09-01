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
# common cell intersect
embryo_rna <- readRDS("/data/bio_data/embryo/rawdata/embro_RNA_soupx_doubletfinder_filtered.rds")
embryo_atac <- readRDS("/data/bio_data/embryo/rawdata/embro_ATAC_filtered.rds")
commoncell <- intersect(colnames(embryo_atac),colnames(embryo_rna))
embryo_atac <- embryo_atac[,commoncell]
embryo_rna <- embryo_rna[,commoncell]
# merge into one seurat object
embryo <- embryo_rna
embryo[['aRNA']] <- embryo_atac[['aRNA']]
embryo[['ATAC']] <- embryo_atac[['ATAC']]
embryo@meta.data[,c('nCount_ATAC','nFeature_ATAC','fragments','FRiP','TSS.enrichment','TSS.percentile','nucleosome_signal','ATAC_snn_res.0.1','ATAC_snn_res.0.3','ATAC_snn_res.0.5','ATAC_snn_res.1','nCount_aRNA','nFeature_aRNA','ATAC_L2_anno')] <- embryo_atac@meta.data[commoncell,c('nCount_ATAC','nFeature_ATAC','fragments','FRiP','TSS.enrichment','TSS.percentile','nucleosome_signal','ATAC_snn_res.0.1','ATAC_snn_res.0.3','ATAC_snn_res.0.5','ATAC_snn_res.1','nCount_aRNA','nFeature_aRNA','ATAC_L2_anno')]
embryo[['lsi']]=embryo_atac[['lsi']]
embryo[['harmony_lsi']]=embryo_atac[['harmony_lsi']]
embryo[['umap.atac']]=embryo_atac[['umap.atac']]
saveRDS(embryo,'embryo/Membryo_rna_atac.rds')
# wnn
embryo <- FindMultiModalNeighbors(embryo, reduction.list = list("harmony", "harmony_lsi"), dims.list = list(1:40, 2:40))
embryo <- RunUMAP(embryo, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
embryo <- FindClusters(embryo, graph.name = "wsnn", algorithm = 3, resolution = c(0.1, 0.3, 0.5), verbose = FALSE)
options(repr.plot.width=13, repr.plot.height=8)
#fig1B
# cluster_cols = c('#6778AE','#F7F398','#91D0BE','#B2DF8A','#58A4C3',   '#aa8282','#666666','#476D87','#55A1B1','#E39A35','#DCC1DD','#E5D2DD','#A6761D','#CCE0F5','#BEAED4')
p1 <- DimPlot(embryo, group.by = "L2_anno", label = TRUE, repel = TRUE,pt.size=0.8, label.size = 4,cols = mycolors,reduction = "wnn.umap",) + 
        theme_dr(xlength = 0.2, ylength = 0.2, arrow = arrow(length = unit(0.2, "inches"),type = "closed")) + 
        theme(panel.grid = element_blank(), axis.title = element_text(face = 2,hjust = 0.03),plot.title = element_blank())

## ATAC anno vs RNA anno
# 
df <- as.data.frame(table(embryo$ATAC_L2_anno, embryo$L2_anno))
colnames(df) <- c("ATAC_L2_anno", "L2_Anno", "Count")

# 
col_order <- df %>% 
  group_by(ATAC_L2_anno) %>%
  # 
  slice_max(order_by = Count, n = 1, with_ties = FALSE) %>% 
  ungroup() %>%
  # 
  arrange(ATAC_L2_anno) %>% 
  pull(L2_Anno) %>% 
  as.character()

# 
all_cols <- levels(df$L2_Anno)
remaining_cols <- setdiff(all_cols, col_order)
final_col_levels <- unique(c(col_order, remaining_cols))

# 
df$L2_Anno <- factor(df$L2_Anno, levels = final_col_levels)

# 
df_plot <- df %>%
  group_by(ATAC_L2_anno) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

# 
options(repr.plot.width = 16, repr.plot.height = 12)
p <- ggplot(df_plot, aes(x = L2_Anno, y = ATAC_L2_anno, fill = Proportion)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#e31a1c", name = "Proportion") + 
  theme_minimal() +
  labs(title = "Mapping Confusion Matrix (Ordered Columns)", 
       x = "L2 Annotation", 
       y = "ATAC L2 Annotation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
  geom_text(aes(label = ifelse(Count > 0, Count, ""), 
                color = Proportion > 0.5), 
            size = 3) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none")
