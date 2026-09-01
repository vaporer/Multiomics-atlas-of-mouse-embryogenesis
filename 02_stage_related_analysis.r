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
set.seed(123)
library(future)
library(future.apply)
plan("multisession", workers = 4)
options(future.globals.maxSize = 40 * 1024^3) # 10GB
embryo <- readRDS('/data/bio_data/zhangzhongjin/embryo/Membryo_rna_atac.rds')
embryo
# figure 1C
## proportion of developmental stages across each cell type
sample_cols <- c(
  "E8"  = "#2C7BB6", 
  "E9"  = "#ABD9E9",
  "E10" = "#FFFFBF",
  "E11" = "#FDAE61",
  "E12" = "#D7191C")
df <- as.data.frame(prop.table(table(embryo$L2_anno, embryo$Sample), margin = 1))
colnames(df) <- c("CellType", "Sample", "Proportion")
options(repr.plot.width=8, repr.plot.height=8)
p1 <- ggplot(df, aes(x = CellType, y = Proportion, fill = Sample)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +  
  theme_bw() +
  scale_fill_manual(values = sample_cols) + #
  labs(title = "Cell Type Composition across Samples",
       x = "Cell Type (L2)",
       y = "Proportion (0-1)") +
  theme(axis.text.y = element_text(size = 8)) 
# figure 2A
## Streamgraph showing dynamic proportion changes of cell lineages from E8.5 to E12.5
sankey_df <- embryo@meta.data %>% group_by(timepoint,L2_anno) %>%
  summarise(n = n(), .groups = 'drop')
library(ggalluvial)
library(scales) 
sankey_prop <- sankey_df %>%
  group_by(timepoint) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()
first_timepoint <- levels(factor(sankey_prop$timepoint))[1]
options(repr.plot.width=8, repr.plot.height=9)
p_sankey <- ggplot(sankey_prop, aes(x = timepoint, stratum = L2_anno, alluvium = L2_anno,
                                   y = proportion, fill = L2_anno, label = L2_anno)) +
  geom_flow(width = 0.4, alpha = 0.4) +
  geom_stratum(width = 0.4, color = 'grey50', size = 0.3) +
  geom_text(stat = 'stratum', 
            aes(label = ifelse(timepoint == first_timepoint, as.character(L2_anno), "")), 
            nudge_x = -0.25, 
            hjust = 1,       
            size = 4, 
           # min.y = 0.01
           ) + 
  scale_fill_manual(values = palette) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = 'Cell Type Composition Across Embryonic Development (E8-E12)',
       x = 'Timepoint', y = 'Cell Proportion', fill = 'Cell Type') +
  scale_x_discrete(expand = expansion(mult = c(0.4, 0.05))) + 
  coord_cartesian(clip = "off") + 
  theme(legend.position = 'none', 
    plot.title = element_text(hjust = 0.5, face = 'bold', size = 12),
    plot.margin = margin(l = 30, r = 10, t = 10, b = 10),
    panel.grid = element_blank(),
  )

# differentially expressed genes (DEGs) for each later stage relative to the preceding stage. 
DefaultAssay(embryo) <- "RNA"
Idents(embryo) <- 'Sample'
target_groups <- c("E8", "E9", "E10", "E11",'E12')
combinations <- combn(target_groups, 2)
num_comparisons <- ncol(combinations)
all_deg_results_list <- future_lapply(seq_len(num_comparisons), function(i) {
  group1 <- combinations[2, i]
  group2 <- combinations[1, i]
  comparison_name <- paste0(group1, "_vs_", group2)
  # message(sprintf("正在处理: %s", comparison_name)) 
  result <- NULL
  tryCatch({
    deg <- FindMarkers(
      object = embryo,
      ident.1 = group1,
      ident.2 = group2,
      min.pct = 0.1,
      logfc.threshold = 0.1,
      only.pos = FALSE,
      verbose = FALSE
    )
    if (nrow(deg) > 0) {
      deg$gene <- rownames(deg)
      deg$comparison <- comparison_name
      deg$group1 <- group1
      deg$group2 <- group2
      result <- deg 
    }
  }, error = function(e) {
    return(NULL) 
  })
  return(result)
}, future.seed = TRUE) 
all_deg_results <- do.call(rbind,all_deg_results_list)
# figure 2B
## summarizing differentially expressed genes (DEGs) for each later stage relative to the preceding stage. 
plot_data <- age_degs %>% mutate(regulation = ifelse(avg_log2FC > 0, "Upregulated", "Downregulated")) %>%
  group_by(comparison, regulation) %>% summarise(gene_count = n(), .groups = 'drop')
plot_data$count_group <- factor(plot_data$comparison, levels = c('E9_vs_E8', 'E10_vs_E9', 'E11_vs_E10', 'E12_vs_E11'))
plot_data$regulation <- factor(plot_data$regulation, levels = c("Upregulated", "Downregulated"))
options(repr.plot.width=5, repr.plot.height=6)
p1 <- ggplot(plot_data, aes(x = comparison, y = gene_count, fill = regulation)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) + 
  scale_fill_manual(values = c('Upregulated' = '#b2182b', 'Downregulated' = '#2166ac')) +
  scale_y_continuous(name = "Number of DEGs (x1000)", 
                     expand = expansion(mult = c(0, 0.1)),
                     labels = function(x) format(x/1000)) +
  labs(title = "DEG Distribution", x = NULL) +
  theme_minimal(base_size = 12) + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey95"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_blank(),
    legend.position = "top")

# differentially expressed regins (DERs) for each later stage relative to the preceding stage. 
DefaultAssay(embryo) <- "ATAC"
Idents(embryo) <- 'Sample'
target_groups <- c("E8", "E9", "E10", "E11",'E12')
combinations <- combn(target_groups, 2)
num_comparisons <- ncol(combinations)
all_der_results_list <- future_lapply(seq_len(num_comparisons), function(i) {
  group1 <- combinations[2, i]
  group2 <- combinations[1, i]
  comparison_name <- paste0(group1, "_vs_", group2)
  result <- NULL
  tryCatch({
    deg <- FindMarkers(
      object = embryo,
      ident.1 = group1,
      ident.2 = group2,
      min.pct = 0.1,
      logfc.threshold = 0.1,
      only.pos = FALSE,
      verbose = FALSE
    )
    if (nrow(deg) > 0) {
      deg$gene <- rownames(deg)
      deg$comparison <- comparison_name
      deg$group1 <- group1
      deg$group2 <- group2
      result <- deg 
    }
  }, error = function(e) {
    # 捕获错误但不中断整个并行过程
    return(NULL) 
  })
  return(result)
}, future.seed = TRUE) 
all_der_results <- do.call(rbind,all_der_results_list)
