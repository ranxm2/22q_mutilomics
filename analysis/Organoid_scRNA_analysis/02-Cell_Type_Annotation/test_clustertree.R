
library(Seurat)
library(clustree)
library(dplyr)
library(tidyverse)
library(tidyr)

seurat_obj <- readRDS(here::here("data", "synaptosomes_scRNA", "seurat-cluster_2025-05-15.rds"))
seurat_obj$source <- seurat_obj$orig.ident
seurat_obj <- subset(seurat_obj, subset = source != "C3" & source != "Q3")
seurat_obj$group <- case_when(
  seurat_obj$source %in% c("C1", "C2") ~ "CTRL",
  seurat_obj$source %in% c("Q1", "Q2") ~ "22q11DS",
  TRUE                                      ~ seurat_obj$source
)

# 2. Compute nearest‐neighbor graph (if you haven’t already)
#    Here we use the first 20 PCs; tweak dims as appropriate.
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

# 3. Run clustering across a sequence of resolutions
#    e.g. from 0.1 to 1.0 in steps of 0.1
resolutions <- seq(0.1, 1.0, by = 0.1)
seurat_obj <- FindClusters(seurat_obj, resolution = resolutions)

# 4. Inspect the column names in metadata to find your cluster labels
head(colnames(seurat_obj@meta.data), 20)
# You should see names like "integrated_snn_res.0.1", "integrated_snn_res.0.2", …

# 5. Plot the cluster tree
#    clustree() will automatically look for meta.data columns
#    that start with the given prefix.
p<-clustree(seurat_obj, prefix = "integrated_snn_res.")

# save the plot
ggsave("result/04-Cluster-Marker-4sample/clustertree_plot.png", plot = p, width = 10, height = 8)

# plot the UMAP for each resolution
for (res in resolutions) {
  p <- DimPlot(seurat_obj, reduction = "umap", group.by = paste0("integrated_snn_res.", res), label=T) +
    ggtitle(paste("UMAP for resolution", res))
  
  # save the plot
  ggsave(paste0("result/04-Cluster-Marker-4sample/umap_res_", res, ".png"), plot = p, width = 10, height = 8)
}



markers_df <- read.csv(
  here::here("data", "synaptosomes_scRNA", "ct_mannual_1126.csv"),
  row.names = 1
)






markers <- markers_df %>%
  mutate(marker = str_split(marker, ",")) %>%
  unnest(marker)

DefaultAssay(seurat_obj) <- "RNA"





library(Seurat)
library(ggplot2)
library(patchwork)

plot_sc_feature <- function(seurat_obj, res = NULL, features = "Test", 
                            group_by = "seurat_clusters",
                            group_by_label = "Resolution: 0.3",
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




res=0.3
Idents(seurat_obj) <- seurat_obj$integrated_snn_res.0.3
seurat_obj$seurat_clusters <- seurat_obj$integrated_snn_res.0.3




for (i in  1:nrow(markers)) {
  gene <- markers[i, "marker"] %>% pull()
  cell_type <- markers[i, "cell_type"]%>% pull()
  print(gene)
  print(cell_type)

  # plot the gene
  file_path <- sprintf("result/04-Cluster-Marker-4sample/res_%s_v2/Gene_%s", res, gene)
  plot_sc_feature(seurat_obj=seurat_obj, res = res, features = gene, save = TRUE,
                  group_by = "seurat_clusters",
                  group_by_label = sprintf("Res_%s", res),
                  reduction_label = "umap",
                  output_prefix = file_path)
  print(sprintf("plot %s for %s", gene, cell_type))
}



markers_df <- read.csv(
  here::here("data", "synaptosomes_scRNA", "ct_mannual_0530.csv"),
  row.names = 1
)

markers <- markers_df %>%
  mutate(marker = str_split(marker, ",")) %>%
  unnest(marker)



for (i in  2:nrow(markers)) {
  gene <- markers[i, "marker"] %>% pull()
  cell_type <- markers[i, "name"]%>% pull()
  print(gene)
  print(cell_type)
  
  # plot the gene
  file_path <- sprintf("result/04-Cluster-Marker-4sample/res_%s_v2/Gene_%s", res, gene)
  plot_sc_feature(seurat_obj=seurat_obj, res = res, features = gene, save = TRUE,
                  group_by = "seurat_clusters",
                  group_by_label = sprintf("Res_%s", res),
                  reduction_label = "umap",
                  output_prefix = file_path)
  print(sprintf("plot %s for %s", gene, cell_type))
}




library(openxlsx)
markers_db <- read.xlsx(here::here("data","ref", "Cell_marker_Seq.xlsx"),
                        check.names = F, sheet = 1)
markers_db %>% filter(
  tissue_type == "Brain",
  species == "Human",
  cancer_type == "Normal") -> markers_db
markers_db %>% group_by(cell_name) %>%
  summarise(total = n()) -> cell_type_marker_count


# Find marker gene for each cluster
cluster_marker <- FindAllMarkers(
  seurat_obj, assay = "RNA", only.pos = T, densify = T) 


seurat_obj_marker <- JoinLayers(seurat_obj, assay = "RNA")

FindAllMarkers(seurat_obj_marker, assay = "RNA", only.pos = T, densify = T) -> cluster.markers

cluster_marker<-cluster.markers
# save the cluster marker
write.csv(
  cluster_marker,
  file = "result/04-Cluster-Marker-4sample/res_0.3_v2/cluster_marker_resolution_0.3.csv",
  row.names = F
)


cluster_marker %>%
  filter(p_val_adj < 0.05,
         avg_log2FC > 0) %>%
  group_by(cluster)-> sig_markers
  # group_by(cluster) %>%
  # slice_max(avg_log2FC, n = 2000) -> sig_markers
table(sig_markers$cluster)

lapply(unique(cluster_marker$cluster), function(cl){
  # cl <- 0
  print(sprintf("processing resolution %s cluster %s", res, cl))
  sig_markers %>%
    filter(cluster == cl) %>%
    arrange(desc(avg_log2FC)) -> cl_marker_genes

  markers_db %>%
    filter(Symbol %in% cl_marker_genes$gene) %>%
    left_join(cl_marker_genes, by = c("Symbol" = "gene")) %>%
    group_by(cell_name) %>%
    summarise(n = n(),
              markers = paste(unique(Symbol), collapse = ",")) %>%
    left_join(cell_type_marker_count) %>%
    mutate(proportion = n / total) %>%
    arrange(desc(proportion)) -> res_anno
  res_anno$cluster <- cl

  # add gene logfc to the annotation res
  message("adding log2FC to the annotation res....")

  sapply(res_anno$markers, function(genes){
    genes <- genes %>% str_split(",", simplify = T) %>% .[1,]
    cluster_marker %>%
      mutate(avg_log2FC = round(avg_log2FC, 4)) %>%
      filter(cluster == cl,
             gene %in% genes) %>%
      arrange(match(gene, genes)) %>%
      dplyr::select(avg_log2FC) %>%
      unlist() %>% paste0(., collapse = ",")
  }) -> res_anno$log2FC

  # res_anno$log2FC %>% sapply(function(x){
  #   str_split(x ,",", simplify = T) %>% .[1,] %>%
  #     as.numeric() %>% sum()
  # }) -> res_anno$log2FC_sum

  # create a new column for avg_log2FC
  res_anno$log2FC %>% sapply(function(x){
    str_split(x ,",", simplify = T) %>% .[1,] %>%
      as.numeric() %>% mean()
  }) -> res_anno$avg_log2FC

  # add gene avg expr to the annotation res
  message("adding avg expr to the annotation res....")

  res_anno %>%
    arrange(desc(proportion))
}) -> cell.type.annotation

names(cell.type.annotation) <- paste0("cluster ", unique(cluster_marker$cluster))

openxlsx::write.xlsx(
  cell.type.annotation,
  file = sprintf("./result/04-Cluster-Marker-4sample/res_%s_v2/cell_type_annotation_resolution_%s.xlsx", res, res)
)

