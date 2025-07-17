library(Seurat)
library(ggplot2)
library(patchwork)

plot_sc_feature <- function(seurat_obj, res = NULL, features = "Test", 
                            group_by = "seurat_clusters",
                            group_by_label = "Resolution: 0.2",
                            reduction_label = "umap_harmony",
                            save = TRUE, output_prefix = "sc_feature_plot") {
  
  # UMAP plot by clusters
  p1 <- DimPlot(seurat_obj, 
                reduction = reduction_label, 
                group.by = group_by, 
                label = TRUE) + 
    ggtitle(group_by_label)
  
  # Feature expression plot
  p2 <- FeaturePlot(seurat_obj, 
                    features = features, 
                    pt.size = 0.1, 
                    reduction = reduction_label) + 
    scale_color_gradientn(colors = c("lightgrey", "blue")) + 
    theme(legend.position = "right") + 
    labs(title = features)
  
  # Violin plot
  p3 <- RidgePlot(seurat_obj, features = features, group.by = group_by)
  
  # Ridge plot
  p4 <- VlnPlot(seurat_obj, features = features, group.by = group_by, pt.size = 0.001)
  
  
  # Combine plots
  combined_plot <- (p1 | p2) / (p3 | p4)
  
  if (save) {
    # ggsave(paste0(output_prefix, ".pdf"), combined_plot, width = 12, height = 10)
    ggsave(paste0(output_prefix, ".png"), combined_plot, width = 12, height = 10, dpi = 300)
    message("Plot saved as PDF and PNG")
  } else {
    message("Plot was not saved")
  }
  
  return(combined_plot)
}

