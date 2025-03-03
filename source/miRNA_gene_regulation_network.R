library(Matrix)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(tibble)
library(rtracklayer)
library(tidyr)
library(dplyr)
library(data.table)
library(tidyverse)
library(gprofiler2)
library(igraph)
library(ggraph)
library(VennDiagram)
library(grid)
library(extrafont)
library(eulerr)
library(RColorBrewer)
library(circlize)
library(scales)

library(paletteer)


clean_miRNA <- function(x) {
  x <- trimws(x)          # Remove leading/trailing whitespace
  x <- gsub("\\.", "-", x) # Replace periods with hyphens
}



plot_network_circle_blue <- function(network_data, path_name, gene_list, result_folder, file_name){
  
  network_data_sub <- network_data %>% filter(gene %in% gene_list)
  
  g <- graph_from_data_frame(network_data_sub, directed = FALSE)
  E(g)$weight <- network_data_sub$correlation
  E(g)$color <- "blue"
  node_types <- data.frame(name = V(g)$name,
                           type = ifelse(V(g)$name %in% gene_list, "mRNA", "miRNA"))
  node_types$log2FC <- ifelse(
    node_types$type == "miRNA",
    network_data_sub$miRNA_log2FC[match(node_types$name, network_data_sub$miRNA)],
    network_data_sub$gene_log2FC[match(node_types$name, network_data_sub$gene)]
  )

  
  # Assign `type` and `in_gene_22q` attributes
  V(g)$type <- node_types$type
  V(g)$logFC <- node_types$log2FC
  
  min_log2FC <- min(node_types$log2FC, na.rm = TRUE)
  max_log2FC <- max(node_types$log2FC, na.rm = TRUE)
  
  # Assign `type` and `in_gene_22q` attributes
  V(g)$type <- node_types$type
  V(g)$log2FC <- node_types$log2FC
  
  E(g)$color <- "blue"
  
  max_correlation <- max(network_data_sub$correlation)
  V(g)$name <- gsub("mir", "miR", V(g)$name)
  
  # Modify the scale_edge_width to show only integer values in the legend
  p <- ggraph(g, 'linear', circular = TRUE) + 
    geom_edge_link(aes(width = weight, color = I(E(g)$color)), alpha = 0.7) + 
    geom_node_point(aes(shape = type, color = log2FC), size = 5) + 
    geom_node_text(aes(label = name), repel = TRUE, size = 5) + 
    scale_shape_manual(values = c("mRNA" = 15, "miRNA" = 16), 
                       guide = guide_legend(title = "Type", order = 1)) + 
    scale_color_gradientn(colors = c("blue", "white", "red"),
                          values = scales::rescale(c(min_log2FC, 0, max_log2FC)), 
                          guide = guide_colorbar(title = "log2FC", order = 2)) + 
    scale_edge_width(breaks = c(1, 2, 3, 4, 5, 6), range = c(0.5, 3), 
                     guide = guide_legend(title = "", order = 3)) + 
    labs(title = paste(path_name)) + 
    expand_limits(x = c(-1.5, 1.5)) +  # Adds space on left and right
    theme_void() + 
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 24, margin = margin(b = 20)), # Adds space below title
      plot.margin = margin(t = 30, r = 10, b = 10, l = 10)  # Adds extra margin to push everything down
    )
  
  print(p)
  
  ggsave(file.path(result_folder, paste0(file_name, ".png")), p, width = 11, height = 8, dpi = 300)
  ggsave(file.path(result_folder, paste0(file_name, ".pdf")), p, width = 11, height = 8)
}
  