library(Seurat)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)


scHeatmap <- function(seu, 
                      features, 
                      group.by = c("cell_type"), 
                      colors = c("#4A90E2", "#D9D9D9", "#D0021B"),
                      col_range = c(-2, 2),
                      anno_colors = NULL,        
                      custom_set_colors = NULL,  
                      column_fontsize = 10,      
                      row_fontsize = 10,
                      legend_fontsize = 10,
                      title_fontsize = 10,
                      name = 'Z-score',
                      save_path = NULL,
                      width = 6, height = 5, units = "in", dpi = 300,
                      split_first = FALSE
) {
  
  meta_df <- seu@meta.data
  
  # --- 1. Handle Input Features (List vs Vector) ---
  if (is.list(features)) {
    gene_sets <- features
    flat_features <- unlist(gene_sets)
    use_gene_sets <- TRUE
  } else {
    flat_features <- features
    gene_sets <- list("Genes" = features)
    use_gene_sets <- FALSE
  }

  # --- 2. Enforce Factor Levels in Metadata ---
  for (col in group.by) {
    if (!is.factor(meta_df[[col]])) {
      meta_df[[col]] <- factor(meta_df[[col]], levels = sort(unique(meta_df[[col]])))
    }
  }

  # --- 3. Run Seurat DotPlot to Get Scaled Data ---
  if (length(group.by) == 1) {
    seu$combined_group <- as.character(meta_df[, group.by])
  } else {
    seu$combined_group <- apply(meta_df[, group.by, drop = FALSE], 1, paste, collapse = " - ")
  }
  
  p <- Seurat::DotPlot(seu, features = flat_features, group.by = "combined_group")
  plot_data <- p$data

  # --- 4. Define Strict Group Ordering ---
  level_list <- lapply(rev(group.by), function(col) levels(meta_df[[col]]))
  group_combos <- expand.grid(level_list, stringsAsFactors = FALSE)
  ordered_group_names <- do.call(paste, c(rev(group_combos), sep = " - "))
  valid_group_names <- ordered_group_names[ordered_group_names %in% unique(plot_data$id)]
  # --- 5. Construct and Order the Matrix ---
  mat <- dcast(plot_data, features.plot ~ id, value.var = "avg.exp.scaled")
  rownames(mat) <- mat$features.plot
  mat <- mat[, -1] 
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0 
  
  # Strict Reordering
  mat <- mat[flat_features, , drop = FALSE]
  mat <- mat[, valid_group_names, drop = FALSE]

  # ============================================================
  # Heatmap Rendering Logic
  # ============================================================

  # --- 6. Reconstruct Column Annotation DataFrame ---
  group_labels <- colnames(mat)
  
  if (length(group.by) == 1) {
    anno_df <- data.frame(v1 = group_labels, stringsAsFactors = FALSE)
    colnames(anno_df) <- group.by
    rownames(anno_df) <- group_labels
  } else {
    group_split <- do.call(rbind, strsplit(group_labels, " - "))
    anno_df <- as.data.frame(group_split)
    colnames(anno_df) <- group.by
    rownames(anno_df) <- group_labels
  }
  
  for (col in group.by) {
    anno_df[[col]] <- factor(anno_df[[col]], levels = levels(meta_df[[col]]))
  }

  # --- 7. Generate Annotation Colors ---
  if (is.null(anno_colors)) {
    anno_colors <- list()
    palette_all <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", 
                     "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
                     "#CCEBC5", "#FFED6F")
    ptr <- 1
    
    for (col in group.by) {
      lvls <- levels(meta_df[[col]])
      n_levels <- length(lvls)
      pal_indices <- ((ptr - 1) + seq_len(n_levels) - 1) %% length(palette_all) + 1
      current_pal <- palette_all[pal_indices]
      names(current_pal) <- lvls
      anno_colors[[col]] <- current_pal
      ptr <- ptr + n_levels 
    }
  }

  # --- 8. Define Color Mapping (Using Parameter) ---
  # NEW: Uses col_range[1] and col_range[2]
  col_fun <- colorRamp2(c(col_range[1], 0, col_range[2]), colors)
  
  # Top Annotation (Group Info)
  top_anno <- HeatmapAnnotation(
    df = anno_df,
    col = anno_colors,
    annotation_name_gp = gpar(fontsize = legend_fontsize),
    show_annotation_name = FALSE, # Hides top annotation names
    gap = unit(2, "mm")
  )
  
  # Row Annotation (Gene Sets)
  row_anno <- NULL
  row_split_var <- NULL
  
  if (use_gene_sets) {
    row_split_var <- rep(names(gene_sets), times = sapply(gene_sets, length))
    row_split_var <- factor(row_split_var, levels = names(gene_sets), ordered = TRUE)
    
    if (is.null(custom_set_colors)) {
      set_pal <- RColorBrewer::brewer.pal(max(3, length(gene_sets)), "Set2")[1:length(gene_sets)]
      custom_set_colors <- setNames(set_pal, names(gene_sets))
    }
    
    row_anno <- rowAnnotation(
      "Gene Set" = row_split_var,
      col = list("Gene Set" = custom_set_colors),
      show_legend = TRUE,
      show_annotation_name = FALSE # <--- NEW: Hides left annotation names
    )
  }

  if (length(group.by) == 1) {
    split_first <- FALSE
  }


  if (split_first) {
    n_splits <- length(unique(anno_df[[group.by[1]]]))
    column_split <- anno_df[[group.by[1]]]
    column_title <- rep("", n_splits)
  }
  else {
    column_split <- NULL
    column_title <- NULL
  }
  # --- 9. Draw Heatmap ---
  ht_params <- list(
    matrix = mat,
    name = name,
    col = col_fun,
    top_annotation = top_anno,
    cluster_rows = FALSE,      
    cluster_columns = FALSE,   
    show_row_names = TRUE,
    show_column_names = FALSE,
    column_split = column_split,
    column_gap = unit(2, "mm"),
    row_names_gp = gpar(fontsize = row_fontsize),
    column_names_gp = gpar(fontsize = column_fontsize),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = title_fontsize, fontface = "bold"),
      labels_gp = gpar(fontsize = legend_fontsize)
    ),
    column_title = column_title
  )
  
  if (use_gene_sets) {
    ht_params$left_annotation <- row_anno
    ht_params$row_split <- row_split_var
    ht_params$row_gap <- unit(2, "mm")
    ht_params$row_title <- rep("", length(unique(row_split_var))) 
  }
  
  ht <- do.call(Heatmap, ht_params)
  
draw_func <- function() {
    draw(ht, heatmap_legend_side = "right", 
         annotation_legend_side = "right", 
         padding = unit(c(2, 2, 2, 2), "mm"))
}

  # Saving options
  if (!is.null(save_path)) {
    ext <- tools::file_ext(save_path)
    if (ext == "pdf") pdf(save_path, width = width, height = height)
    else if (ext == "png") png(save_path, width = width, height = height, units = units, res = dpi)
    else if (ext == "tiff") tiff(save_path, width = width, height = height, units = units, res = dpi)
    else stop("Unsupported file extension: use .pdf, .png, or .tiff")
    
    draw_func()
    dev.off()
  }

  draw_func()
  
  # return(list(mat = mat, anno_df=anno_df))
}