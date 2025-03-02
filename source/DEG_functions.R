library(Seurat)
library(Signac)
library(Matrix)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(apeglm)
library(EnhancedVolcano)
library(pheatmap)
library(tibble)
library(rtracklayer)
library(tidyr)
library(ggfortify)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(gprofiler2)
library(biomaRt)
library(GenomicRanges)
library(GenomicFeatures)
library(edgeR)
library(parallel)
library(BiocParallel)
library(GSVA)
library(decoupleR)
library(ggrepel)
library(patchwork)
library(ggsignif) 
library(extrafont)



Result_folder_structure <- function(result_folder) {
  dir.create(result_folder, showWarnings = FALSE)
  dir.create(paste0(result_folder, "/01-Sample_info"), showWarnings = FALSE)
  dir.create(paste0(result_folder, "/02-DEG"), showWarnings = FALSE)
  dir.create(paste0(result_folder, "/03-Enrichment"), showWarnings = FALSE)
  dir.create(paste0(result_folder, "/04-GSEA"), showWarnings = FALSE)
  dir.create(paste0(result_folder, "/05-22q_Gene"), showWarnings = FALSE)
}

# Test for the function
# result_folder = "results/02-DEG-Vehicle"
# Result_folder_structure(result_folder)



DEAnalysis <- function(counts =NULL, 
                       reference_group = NULL, 
                       compare_group = NULL, 
                       condition_list = NULL, 
                       target_gene =NULL, 
                       result_folder = NULL) {
  # Filter the count data 
  sample_info <- condition_list %>%
    mutate(sample = rownames(.)) %>%
    filter(group %in% c(reference_group, compare_group))
  count_input <- counts[, sample_info$sample]
  
  # filter gene withall 0
  count_input <- count_input[rowSums(count_input) > 0, ]
  
  # Create a DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = count_input,
                                colData = sample_info,
                                design = ~ group)
  dds$group <- relevel(dds$group, ref = reference_group)
  

  ##############################
  ### Run DESeq analysis    ####
  ###------------------------###
  
  dds <- DESeq(dds)
  res <- results(dds)
  resOrdered <- res[order(res$padj), ]
  
  # omit the NA values
  resOrdered <- resOrdered[!is.na(resOrdered$padj),]
  dds <- dds[rownames(resOrdered),]
  write.csv(resOrdered, file.path(result_folder,"02-DEG", "01_all_gene_results.csv"))
  
  # Filter the DEG with padj < 0.05 and log2FoldChange > 1
  dem_1 <- resOrdered %>% as.data.frame() %>% rownames_to_column(var = "gene") %>% 
    filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% arrange(padj)
  dem_1 <- dem_1[!is.na(dem_1$padj),]
  write.csv(dem_1, file.path(result_folder,"02-DEG","02_DEG_log2fc_1.csv"), row.names = FALSE)
  
  # 
  dem_1_5 <- resOrdered %>% as.data.frame() %>% rownames_to_column(var = "gene") %>% 
    filter(padj < 0.05 & abs(log2FoldChange) > 1.5) %>% arrange(padj)
  dem_1_5 <- dem_1_5[!is.na(dem_1_5$padj),]
  write.csv(dem_1_5, file.path(result_folder,"02-DEG","03_DEG_log2fc_1_5.csv"), row.names = FALSE)
  print("DEG analysis is done")

  return(dds)
}



plot_sample_heatmap <- function(dds_obj,figure_folder, file_name) {
   # Now apply variance stabilizing transformation
   vsd.obj <- varianceStabilizingTransformation(dds_obj, blind = TRUE)
   sampleDists <- dist(t(assay(vsd.obj)))
   sampleDistMatrix <- as.matrix( sampleDists )
   rownames(sampleDistMatrix) <- paste( vsd.obj$group )
   colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
   p <- pheatmap::pheatmap(sampleDistMatrix,
                           clustering_distance_rows = sampleDists,
                           clustering_distance_cols = sampleDists,
                           col = colors) 
   ggsave(file.path(figure_folder, paste0(file_name, ".png")), 
          p, width = 8, height = 6, units = "in", dpi = 300)
  ggsave(file.path(figure_folder, paste0(file_name, ".pdf")), 
         p, width = 8, height = 6, units = "in", dpi = 300)
   print("Sample distance heatmap is done")
}

# # Test for the function
# plot_sample_heatmap(dds_obj, figure_folder = 
#                       file.path(result_folder,"01-Sample_info"), 
#                     file_name = "01_sample_distance_heatmap")




plot_sample_PCA_plot <- function(dds_obj, figure_folder, file_name, 
                               reference_group, compare_group){

  vsd.obj <- varianceStabilizingTransformation(dds_obj, blind = TRUE)
  pcaData <- plotPCA(vsd.obj,  intgroup = c("group"), returnData = T)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  # Define unique groups and their corresponding colors
  unique_groups <- c(reference_group, compare_group)
  color_palette <- c("#10d5da", "#fe867f")  # Extend this list if you have more groups
  assigned_colors <- setNames(color_palette[seq_along(unique_groups)], unique_groups)
  
  p <-ggplot(pcaData, aes(PC1, PC2, color=group)) +
    geom_point(size=3) +
    labs(x = paste0("PC1: ",percentVar[1],"% variance"),
         y = paste0("PC2: ",percentVar[2],"% variance"),
    ) +
    stat_ellipse(level = 0.95)+
    theme_bw() +
    # theme_classic()+
    theme(text = element_text(family = "Arial", colour = "black")) +
    scale_color_manual(values = assigned_colors) +
    ggrepel::geom_text_repel(aes(label = name), color = "black")
  ggsave(file.path(figure_folder, paste0(file_name, ".png")), 
         p,width = 6, height = 4, units = "in", dpi = 300)
  ggsave(file.path(figure_folder, paste0(file_name, ".pdf")), 
         p,width = 6, height = 4, units = "in")
  print("PCA plot is done")
}

# # test for the function
# plot_sample_PCA_plot(dds_obj,figure_folder = file.path(result_folder,"01-Sample_info"),
#                    file_name = "02_sample_PCA_plot", 
#                    reference_group, compare_group)


plot_volcanoplot <- function(result_df, figure_folder, file_name, label_gene = NULL) {
  
  deg <- results(dds)
  deg <- as.data.frame(deg)
  result_df <- result_df  %>% filter(!is.na(padj))
  max_p  <- max(-log10(deg$padj), na.rm = TRUE)
  max_fc <- max(abs(deg$log2FoldChange), na.rm = TRUE)
  
  deg$diffexpressed <- "NO"
  deg$diffexpressed[deg$log2FoldChange > 1 & deg$padj < 0.05] <- "UP"
  deg$diffexpressed[deg$log2FoldChange < -1 & deg$padj < 0.05] <- "DOWN"
  num_up   <- sum(deg$diffexpressed == "UP")
  num_down <- sum(deg$diffexpressed == "DOWN")
  p <-ggplot(deg, aes(log2FoldChange, -log10(padj), color = diffexpressed)) +
    geom_point(show.legend = F) +
    scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
    scale_x_continuous(limits = c(- max_fc , max_fc )) +
    scale_y_continuous(limits = c(0,max_p  )) +
    # geom_label_repel(max.overlaps = 20, show.legend = F, nudge_y = 2, nudge_x = 2) +
    geom_vline(xintercept = c(-1, 1)) +
    geom_hline(yintercept = -log10(0.05)) +
    theme_classic() +
    theme(legend.position = "right") +
    annotate("text", x = 0.2 * max_fc + 1, y = 0.9 * max_p, 
             label = paste("UP:", num_up), color = "red",size=6 ) +
    annotate("text", x = -0.2 * max_fc - 1, y = 0.9 * max_p, 
             label = paste("DOWN:", num_down), color = "blue",size=6) +
    theme( text = element_text(family = "Arial", colour = "black"),
           axis.title = element_text(color = "black", family = "Arial", size = 16),  # Increase axis title size
           axis.text = element_text(color = "black", family = "Arial", size = 14),   # Increase axis text size
           legend.position = "right"
    )
  
  ggsave(file.path(output_dir,"01-DEG" , "05_volcano_plot_all_gene_log2fc_1.png"), p, width = 8, height = 8, units = "in", dpi = 300)
  ggsave(file.path(output_dir, "01-DEG" ,"05_volcano_plot_all_gene_log2fc_1.pdf"), p, width = 8, height = 8, units = "in", dpi = 300)
  
  }





