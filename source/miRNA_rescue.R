library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(here)
library(dplyr)
library(tidyr)
library(ggsci)

# assume these objects are already defined in your workspace:
#   target_gene_214      (the data.frame you read & renamed)
#   counts                (your raw count matrix)
#   sample_list_raw       (your sample info tibble)
#   comparisons_214       (your list of pairs)
#   result_root           ("./results/DEG-comparisons")
#   Result_folder_structure()  (your helper to mkdir -p)
#   DEAnalysis(), varianceStabilizingTransformation()


make_miRNA_target_heatmap <- function(pair, 
                                      target_gene_df, 
                                      counts, 
                                      suffix = "",
                                      filter_mean = TRUE,
                                      sample_list_raw,
                                      delete_patterns = character(),
                                      result_folder = character(),
                                      width    = 8, height = 16) {
  
  compare_group   <- pair[1]
  reference_group <- pair[2]
  pair_name       <- paste0(compare_group, "_vs_", reference_group)

  # 1. filter samples & run DE + VST
  filter_info   <- sample_list_raw %>% filter(group %in% c(reference_group, compare_group))
  rownames(filter_info) <- filter_info$"Sample.ID.."
  filter_counts <- counts[, filter_info$'Sample.ID..', drop = FALSE]

  
  # filter gene withall 0
  filter_counts  <- filter_counts [rowSums(filter_counts ) > 0, ]
  
  # Create a DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = filter_counts,
                                colData = filter_info,
                                design = ~ group)
  dds$group <- relevel(dds$group, ref = reference_group)
  dds <- DESeq(dds)
  vsd  <- varianceStabilizingTransformation(dds, blind = TRUE)
  
  # 2. pull out the assay & center
  anno_col <- as.data.frame(colData(vsd)[, "group", drop = FALSE])
  stab     <- assay(vsd) - rowMeans(assay(vsd), na.rm = TRUE)

  
  # 3. join with your target‐gene table
  tg <- target_gene_df %>%
    mutate(term_name = gsub("Neuro", "neuro", term_name))
  
  intersect_genes <- intersect(tg$gene, rownames(stab))
  
  # restrict to only genes we care about
  tg_stab <- stab[intersect_genes, rownames(anno_col), drop = FALSE ]
  tg_stab <- as.data.frame(tg_stab ) %>%
    rownames_to_column(var = "gene")
  
  tg_df   <- left_join(
    tg ,
    tg_stab,
    by ="gene")
  
  # 4. compute group means & filter by direction
  ref_samps  <- rownames(anno_col)[ anno_col$group == reference_group ]
  comp_samps <- rownames(anno_col)[ anno_col$group == compare_group ]
  
  tg_df <- tg_df %>%
    rowwise() %>%
    mutate(
      ref_mean  = mean(c_across(all_of(ref_samps))),
      comp_mean = mean(c_across(all_of(comp_samps)))
    ) %>%
    ungroup() 
  
  
  if (filter_mean) {
    tg_df <-   tg_df %>% filter(if (grepl("OE", compare_group)) ref_mean > comp_mean else ref_mean < comp_mean)
  }
  
  # 5. drop any unwanted gene‐modules
  if (length(delete_patterns)) {
    for (pat in delete_patterns) {
      tg_df <- tg_df %>% filter(!grepl(pat, gene_module))
    }
  }
  
  # 6. save the filtered table
  write.csv(
    tg_df,
    file = file.path(result_folder, paste0(pair_name, suffix,"_miRNA_target_gene_module.csv")),
    row.names = FALSE
  )
  
  # 7. build heatmap matrix & annotations
  mat <- tg_df %>%
    dplyr::select(all_of(rownames(anno_col))) %>%
    as.matrix()
  rownames(mat) <- tg_df$gene_module
  
  annotation_row <- data.frame(Module = tg_df$term_name, row.names = tg_df$gene_module)
  module_cols    <- setNames(
    brewer.pal(min(n_distinct(tg_df$term_name),12),"Set3"),
    unique(tg_df$term_name)
  )
  anno_colors <- list(Module = module_cols, group = setNames(
    c("#10d5da","#fe867f"),
    c(reference_group, compare_group)
  ))
  
  # 8. color ramp
  mr <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(256)
  
  # 9. draw & save transparent PNG + PDF
  for (cluster_opt in c(FALSE, TRUE)) {
      suffix_label <- if (cluster_opt) paste0("_clustered",suffix) else suffix
      
      for (fmt in c( "pdf","png")) {

      out_file <- file.path(
        result_folder,
        paste0(pair_name, "_miRNA_target_heatmap", suffix_label , ".", fmt)
      )
      
      # open the appropriate device
      if (fmt == "png") {
        png(
          filename = out_file,
          width    = width, height = height, units = "in", res = 300,
          bg       = "transparent"
        )
      } else {
        pdf(
          file = out_file,
          width = width, height = height,
          bg    = "transparent"
        )
      }
      
      # draw the heatmap
      pheatmap(
        mat,
        scale            = "row",
        color            = mr,
        annotation_row   = annotation_row,
        annotation_colors= anno_colors,
        annotation_col   = anno_col,
        cluster_rows     = cluster_opt,
        cluster_cols     = FALSE,
        fontsize_row     = 10,
        fontsize_col     = 10,
        treeheight_row   = 0,
        show_rownames    = TRUE,
        border_color     = NA
      )
      # close the device
      
      dev.off()
      print(paste("Saved heatmap to", out_file))
      # if (fmt == "png") {
      #   p<-pheatmap(
      #     mat,
      #     scale            = "row",
      #     color            = mr,
      #     annotation_row   = annotation_row,
      #     annotation_colors= anno_colors,
      #     annotation_col   = anno_col,
      #     cluster_rows     = cluster_opt,
      #     cluster_cols     = FALSE,
      #     fontsize_row     = 10,
      #     fontsize_col     = 10,
      #     treeheight_row   = 0,
      #     show_rownames    = TRUE,
      #     border_color     = NA
      #   )
      #   print(p)
      # }
      }
  }
  
  
  
  
}



make_miRNA_target_heatmap_2_way <- function(pair, 
                                      target_gene_df, 
                                      counts, 
                                      suffix = "",
                                      filter_mean = TRUE,
                                      sample_list_raw,
                                      delete_patterns = character(),
                                      result_folder = character(),
                                      width    = 8, height = 16) {
  png(
    filename = file.path(result_folder,"clean.plot.png"),
    width    = 8, height = 16, units = "in", res = 300,
    bg       = "transparent"
  )
  dev.off()
  
  compare_group   <- pair[1]
  reference_group <- pair[2]
  pair_name       <- paste0(compare_group, "_vs_", reference_group)
  
  # 1. filter samples & run DE + VST
  filter_info   <- sample_list_raw %>% filter(group %in% c(reference_group, compare_group))
  rownames(filter_info) <- filter_info$"Sample.ID.."
  filter_counts <- counts[, filter_info$'Sample.ID..', drop = FALSE]
  
  
  # filter gene withall 0
  filter_counts  <- filter_counts [rowSums(filter_counts ) > 0, ]
  
  # Create a DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = filter_counts,
                                colData = filter_info,
                                design = ~ group)
  dds$group <- relevel(dds$group, ref = reference_group)
  dds <- DESeq(dds)
  vsd  <- varianceStabilizingTransformation(dds, blind = TRUE)
  
  # 2. pull out the assay & center
  anno_col <- as.data.frame(colData(vsd)[, "group", drop = FALSE])
  stab     <- assay(vsd) - rowMeans(assay(vsd), na.rm = TRUE)
  
  
  # 3. join with your target‐gene table
  tg <- target_gene_df %>%
    mutate(term_name = gsub("Neuro", "neuro", term_name))
  
  intersect_genes <- intersect(tg$gene, rownames(stab))
  
  # restrict to only genes we care about
  tg_stab <- stab[intersect_genes, rownames(anno_col), drop = FALSE ]
  tg_stab <- as.data.frame(tg_stab ) %>%
    rownames_to_column(var = "gene")
  
  tg_df   <- left_join(
    tg ,
    tg_stab,
    by ="gene")
  
  # 4. compute group means & filter by direction
  ref_samps  <- rownames(anno_col)[ anno_col$group == reference_group ]
  comp_samps <- rownames(anno_col)[ anno_col$group == compare_group ]
  
  tg_df <- tg_df %>%
    rowwise() %>%
    mutate(
      ref_mean  = mean(c_across(all_of(ref_samps))),
      comp_mean = mean(c_across(all_of(comp_samps)))
    ) %>%
    ungroup() 
  

  
  
  # 5. drop any unwanted gene‐modules
  if (length(delete_patterns)) {
    for (pat in delete_patterns) {
      tg_df <- tg_df %>% filter(!grepl(pat, gene_module))
    }
  }
  
  # 6. save the filtered table
  write.csv(
    tg_df,
    file = file.path(result_folder, paste0(pair_name, suffix,"_miRNA_target_gene_module.csv")),
    row.names = FALSE
  )
  
  
  
  for (direction in c("up", "down")) {
    if (direction == "up") {
      tg_df_sub <- tg_df %>% 
        filter(ref_mean  < comp_mean)
    } else {
      tg_df_sub <- tg_df %>% 
        filter(ref_mean  > comp_mean)
    }
    
    # save the filtered table
    write.csv(
      tg_df_sub,
      file = file.path(result_folder, paste0(pair_name, suffix,"_miRNA_target_gene_module_", direction, ".csv")),
      row.names = FALSE
    )
    
    # 7. build heatmap matrix & annotations
    mat <- tg_df_sub %>%
      dplyr::select(all_of(rownames(anno_col))) %>%
      as.matrix()
    rownames(mat) <- tg_df_sub$gene_module
    
    annotation_row <- data.frame(Module = tg_df_sub$term_name, row.names = tg_df_sub$gene_module)
    module_cols    <- setNames(
      brewer.pal(min(n_distinct(tg_df_sub$term_name),12),"Set3"),
      unique(tg_df_sub$term_name)
    )
    anno_colors <- list(Module = module_cols, group = setNames(
      c("#10d5da","#fe867f"),
      c(reference_group, compare_group)
    ))
    
    # 8. color ramp
    mr <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(256)
    
for (fmt in c("pdf","png")) {
  for (cluster_opt in c(FALSE, TRUE)) {
    suffix_label <- if (cluster_opt) paste0("_clustered", suffix) else suffix
    out_file <- file.path(
      result_folder,
      paste0(pair_name, "_miRNA_target_heatmap", suffix_label, ".", fmt)
    )
    # 1) open device first:
    if (fmt == "png") {
      png(out_file,
          width    = width, height = height,
          units    = "in", res = 300,
          bg       = "transparent")
    } else {
      pdf(out_file,
          width    = width, height = height,
          bg       = "transparent")
    }

    # 2) immediately draw the heatmap to that device
    pheatmap(
      mat,
      scale             = "row",
      color             = mr,
      annotation_row    = annotation_row,
      annotation_col    = anno_col,
      annotation_colors = annotation_colors,
      annotation_legend = FALSE,
      cluster_rows      = cluster_opt,
      cluster_cols      = FALSE,
      fontsize_row      = 10,
      fontsize_col      = 10,
      treeheight_row    = 0,
      show_rownames     = TRUE,
      border_color      = NA,
      angle_col         = 45
    )

    # 3) close it
    dev.off()
    message("Saved heatmap to ", out_file)
  }
}

    
    
}
}





make_miRNA_target_heatmap_label <- function(pair, 
                                      target_gene_df, 
                                      counts, 
                                      suffix = "",
                                      filter_mean = TRUE,
                                      sample_list_raw,
                                      delete_patterns = character(),
                                      result_folder = character(),
                                      width    = 8, height = 16) {
  
  compare_group   <- pair[1]
  reference_group <- pair[2]
  pair_name       <- paste0(compare_group, "_vs_", reference_group)
  
  # 1. filter samples & run DE + VST
  filter_info   <- sample_list_raw %>% filter(group %in% c(reference_group, compare_group))
  rownames(filter_info) <- filter_info$"Sample.ID.."
  filter_counts <- counts[, filter_info$'Sample.ID..', drop = FALSE]
  
  # filter gene withall 0
  filter_counts  <- filter_counts [rowSums(filter_counts ) > 0, ]
  
  # Create a DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = filter_counts,
                                colData = filter_info,
                                design = ~ group)
  dds$group <- relevel(dds$group, ref = reference_group)
  dds  <- DESeq(dds)
  vsd  <- varianceStabilizingTransformation(dds, blind = TRUE)
  
  # 2. pull out the assay & center
  anno_col <- as.data.frame(colData(vsd)[, "group", drop = FALSE])
  stab     <- assay(vsd) - rowMeans(assay(vsd), na.rm = TRUE)
  
  # 3. join with your target‐gene table
  tg <- target_gene_df %>%
    mutate(term_name = gsub("Neuro", "neuro", term_name))
  
  clean_df <- target_gene_df %>%
    # 1. remove the gene_module column
    dplyr::select(-gene_module) %>% 
    # 2. keep only one row per gene‐term combination (in case of duplicates)
    distinct(gene, term_name, .keep_all = TRUE) %>%  
    # 3. create a 'present' flag for pivoting
    mutate(present = 1) %>%               
    # 4. pivot term_name into one‐hot columns, filling missing with 0
    pivot_wider(
      id_cols      = gene,
      names_from   = term_name,
      values_from  = present,
      values_fill  = list(present = 0)
    )
  
  # left join the tg and clean_df
  tg <- left_join(tg, clean_df, by = "gene")
  
  intersect_genes <- intersect(tg$gene, rownames(stab))
  
  # restrict to only genes we care about
  tg_stab <- stab[intersect_genes, rownames(anno_col), drop = FALSE ]
  tg_stab <- as.data.frame(tg_stab ) %>%
    rownames_to_column(var = "gene")
  
  tg_df   <- left_join(
    tg ,
    tg_stab,
    by ="gene")
  
  # 4. compute group means & filter by direction
  ref_samps  <- rownames(anno_col)[ anno_col$group == reference_group ]
  comp_samps <- rownames(anno_col)[ anno_col$group == compare_group ]
  
  tg_df <- tg_df %>%
    rowwise() %>%
    mutate(
      ref_mean  = mean(c_across(all_of(ref_samps))),
      comp_mean = mean(c_across(all_of(comp_samps)))
    ) %>%
    ungroup() 
  
  
  if (filter_mean) {
    tg_df <-   tg_df %>% filter(if (grepl("OE", compare_group)) ref_mean > comp_mean else ref_mean < comp_mean)
  }
  
  # 5. drop any unwanted gene‐modules
  if (length(delete_patterns)) {
    for (pat in delete_patterns) {
      tg_df <- tg_df %>% filter(!grepl(pat, gene_module))
    }
  }
  
  # 6. save the filtered table
  write.csv(
    tg_df,
    file = file.path(result_folder, paste0(pair_name, suffix,"_miRNA_target_gene_module.csv")),
    row.names = FALSE
  )
  
  # 7. build heatmap matrix & annotations
  mat <- tg_df %>%
    dplyr::select(all_of(rownames(anno_col))) %>%
    as.matrix()
  rownames(mat) <- tg_df$gene
  
  # only keep the unique gene names for the heatmap
  mat <- mat[!duplicated(tg_df$gene), , drop = FALSE]
  
  # how many distinct pathways?
  n_terms <- n_distinct(tg_df$term_name)
  # get exactly that many colors from the D3 “category10” palette
  d3_cols <- pal_d3("category10")(n_terms)
  
  # assign them to your term names
  module_cols <- setNames(
    d3_cols,
    unique(substr(tg_df$term_name, 1, 2)) # use first two letters of term_name
  )
  
  annotation_row <- clean_df %>%
    filter(gene %in% rownames(mat)) %>%
    column_to_rownames("gene") 
  colnames(annotation_row) <- substr(colnames(annotation_row), 1, 2)
  # reverse the column order to match the original term_name order
  annotation_row <- annotation_row[, rev(colnames(annotation_row))]
  
  annotation_row  <- annotation_row %>% dplyr::select(all_of(names(module_cols))) %>%
    mutate(across(everything(),
                  ~ factor(as.character(.), levels = c("0","1"))))
  
  # ---- 2. rebuild annotation_colors as a flat list ----
  annotation_colors <- list()
  
  # for each module (M1, M2, …), map 0->white, 1->its color
  for(term in names(rev(module_cols))) {
    annotation_colors[[term]] <- c(
      "0" = "white",
      "1" = module_cols[[term]]
    )
  }
  
  # then add your 'group' mapping, making sure the names match your factor levels
  grp_levels <- levels(anno_col$group)
  annotation_colors[["group"]] <- setNames(
    c("#10d5da", "#fe867f"),   # your two colors
    grp_levels                 # e.g. c("22q002","22q002-128-OE")
  )
  
  

  # 8. color ramp
  mr <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(256)
  annotation_row <- annotation_row[, rev(colnames(annotation_row))]
  # 9. draw & save transparent PNG + PDF

  

  for (fmt in c("pdf","png")) {
    for (cluster_opt in c(FALSE, TRUE)) {
      suffix_label <- if (cluster_opt) paste0("_clustered", suffix) else suffix
      out_file <- file.path(
        result_folder,
        paste0(pair_name, "_miRNA_target_heatmap", suffix_label, ".", fmt)
      )
      # 1) open device first:
      if (fmt == "png") {
        png(out_file,
            width    = width, height = height,
            units    = "in", res = 300,
            bg       = "transparent")
      } else {
        pdf(out_file,
            width    = width, height = height,
            bg       = "transparent")
      }
      
      # 2) immediately draw the heatmap to that device
      pheatmap(
        mat,
        scale             = "row",
        color             = mr,
        annotation_row    = annotation_row,
        annotation_col    = anno_col,
        annotation_colors = annotation_colors,
        annotation_legend = FALSE,
        cluster_rows      = cluster_opt,
        cluster_cols      = FALSE,
        fontsize_row      = 10,
        fontsize_col      = 10,
        treeheight_row    = 0,
        show_rownames     = TRUE,
        border_color      = NA,
        angle_col         = 45
      )
      
      # 3) close it
      dev.off()
      message("Saved heatmap to ", out_file)
    }
  }
  

  
  
  
  
}

pair <- comparisons_128 [[2]]
make_miRNA_target_heatmap_label(pair= pair, 
                                target_gene_df = target_gene_df, 
                                counts = counts, 
                                sample_list_raw = sample_list_raw,
                                delete_patterns = c("^TEAD4"),
                                result_folder = './results/rescue_target_gene_module_label',
                                width  = 8, height = 8
) 
