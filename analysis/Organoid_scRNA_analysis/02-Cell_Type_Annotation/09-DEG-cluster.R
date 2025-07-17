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

# # Import fonts (this step may take some time and may need admin privileges)
# font_import(prompt = FALSE)
# 
# # Load fonts for PDF outputs
# loadfonts(device = "pdf")

setwd("~/local_research/Wen/Weibo/02-22q/2-analysis/3-pipeline/scRNA")
set.seed(2024)



# load the data
merged <- readRDS("process_data/04_merged_with_clustering_4_mannual_correct.rds")



####################### cluster DEG
merged$orig.group <- ifelse(merged$orig.ident %in% c("C1", "C2"), "CTRL", "22q")
merged$deg_group <- paste0(merged$orig.group, "_", merged$seurat_clusters_merge)
table(merged$deg_group)

DefaultAssay(merged) <- "RNA"
Idents(merged) <- "deg_group"
lapply(unique(merged$seurat_clusters_merge), function(cl){
  print(cl)
  g1 <- paste0("CTRL_", cl)
  g2 <- paste0("22q_", cl)
  deg <- FindMarkers(merged, ident.1 = g1, ident.2 = g2)
  deg$cluster <- cl
  deg$gene <- rownames(deg)
  deg
}) -> deg_list
names(deg_list) <- paste0(unique(merged$seurat_clusters_merge))

# make. avg_log2FC the -avg_log2FC
deg_list <- lapply(deg_list, function(deg){
  deg$avg_log2FC <- -deg$avg_log2FC
  deg
})




library(stringr)
# assume deg_list is a named list, one element per cluster
old_names <- names(deg_list)
# truncate every name to at most 31 characters (no “…”)
new_names <- str_trunc(old_names, width = 31, side = "right", ellipsis = "")
names(deg_list) <- new_names

# now write
openxlsx::write.xlsx(
  deg_list,
  "result/07-DEG-cluster/cluster_DEGs.xlsx"
)



library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# your vector of logFC cutoffs
logFCs <- c(0.1, 0.25, 0.5, 1)

for (th in logFCs) {
  # 1. tally up/down per cluster at this threshold
  counts_df <- imap_dfr(deg_list, ~{
    df2 <- .x %>%
      filter(p_val_adj < 0.05, abs(avg_log2FC) > th)
    tibble(
      cluster = .y,
      up      = sum(df2$avg_log2FC >  0, na.rm = TRUE),
      down    = sum(df2$avg_log2FC <  0, na.rm = TRUE)
    )
  })
  
  # 2. pivot longer and make 'down' negative
  counts_long <- counts_df %>%
    pivot_longer(c(up, down), names_to = "direction", values_to = "n") %>%
    mutate(
      n         = ifelse(direction == "down", -n, n),
      direction = factor(direction, levels = c("down","up"))
    )
  
  # figure out padding
  pad <- max(abs(counts_long$n)) * 0.1
  
  # 3. build the plot
  p <- ggplot(counts_long, aes(x = n, y = cluster, fill = direction)) +
    geom_col(width = 0.7) +
    geom_text(
      aes(
        x     = n,
        label = abs(n),
        hjust = ifelse(direction == "up", -0.1, 1.1)
      ),
      size = 3
    ) +
    scale_x_continuous(
      labels = abs,
      # add extra room on both sides (you may need to tweak these multipliers)
      expand = expansion(add = c(max(abs(counts_long$n)) * 0.2,
                                 max(abs(counts_long$n)) * 0.2))
    ) +
    scale_fill_manual(
      values = c(down = "steelblue", up = "firebrick"),
      labels = c("Down-regulated","Up-regulated")
    ) +
    labs(
      title = paste0("DEG counts (|log₂FC| > ", th, ")"),
      x     = "Number of DE genes",
      y     = "Cell type (cluster)",
      fill  = "Direction"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position     = "top",
      axis.text.y         = element_text(size = 10),
      # bump out the right margin (t, r, b, l)
      plot.margin         = margin(5, 30, 5, 5, "pt")
    ) +
    # allow drawing of labels beyond the panel
    coord_cartesian(clip = "off")
  
  # 4. save it
  ggsave(
    filename = paste0("result/07-DEG-cluster/diverging_bar_logFC_", th, ".png"),
    plot     = p,
    width    = 8, 
    height   = 6,
    dpi      = 300
  )
}






lapply(deg_list, function(deg){
  deg %>%
    filter(p_val_adj < 0.05) %>%
    filter(abs(avg_log2FC) > 1)
}) -> deg_list_sign_1


str(deg_list_sign)


lapply(deg_list_sign, function(deg){
  print(deg$cluster[1])
  up.enrichment <- gost(deg %>% filter(avg_log2FC > 0) %>% rownames(),
                        organism = "hsapiens", correction_method = "fdr", evcodes = T)
  down.enrichment <- gost(deg %>% filter(avg_log2FC < 0) %>% rownames(),
                        organism = "hsapiens", correction_method = "fdr", evcodes = T)
  write.csv(subset(up.enrichment$result, select = -c(parents)),
            file = sprintf("result/07-DEG-cluster/cluster_DEG_enrichment_logFC_1/%s_up_enrichment.csv", deg$cluster[1]),
            row.names = F)
  write.csv(subset(down.enrichment$result, select = -c(parents)),
            file = sprintf("result/07-DEG-cluster/cluster_DEG_enrichment_logFC_1/%s_down_enrichment.csv", deg$cluster[1]),
            row.names = F)
})



library(dplyr)
library(purrr)
library(gprofiler2)

# your vector of log₂FC thresholds:
logFCs <- c(1, 0.5, 0.25)

for (th in logFCs) {
  message("Running enrichment for |log2FC| > ", th)
  # make a clean directory for this threshold
  outdir <- file.path("result/07-DEG-cluster",
                      paste0("cluster_DEG_enrichment_logFC_", th))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # filter each cluster once
  deg_list_sign <- map(deg_list, ~
                         .x %>%
                         filter(p_val_adj < 0.05, abs(avg_log2FC) > th)
  )
  
  # run enrichment per cluster
  walk(deg_list_sign, function(deg) {
    # get a safe name
    cluster_name <- unique(deg$cluster)
    safe_name    <- gsub("\\s+", "_", cluster_name)
    
    # gene sets
    up_genes   <- rownames(filter(deg, avg_log2FC >  0))
    down_genes <- rownames(filter(deg, avg_log2FC <  0))
    
    # run gost
    up_enr   <- gost(up_genes,   organism = "hsapiens",
                     correction_method = "fdr", evcodes = TRUE)
    
    # write results
    write.csv(
      dplyr::select(up_enr$result, -parents),
      file = file.path(outdir,
                       paste0(safe_name, "_up_enrichment.csv")),
      row.names = FALSE
    )
    
    
    down_enr <- gost(down_genes, organism = "hsapiens",
                     correction_method = "fdr", evcodes = TRUE)
    
    write.csv(
      dplyr::select(down_enr$result, -parents),
      file = file.path(outdir,
                       paste0(safe_name, "_down_enrichment.csv")),
      row.names = FALSE
    )
  })
  message("Done enrichment at |log2FC| > ", th)
}
