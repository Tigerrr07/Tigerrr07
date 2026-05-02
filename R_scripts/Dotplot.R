library(Seurat)
library(ggplot2)
# library(viridis)


scDotPlot <- function(seu, features, group.by, my_colors = c("#4A90E2", "#D9D9D9", "#D0021B")) {
  
  meta_df <- seu@meta.data
  for (col in group.by) {
    if (!is.factor(meta_df[[col]])) {
      meta_df[[col]] <- factor(meta_df[[col]], levels = sort(unique(meta_df[[col]])))
    }
  }
  if (length(group.by) == 1) {
    seu$combined_group <- as.character(meta_df[, group.by])
  } else {
    seu$combined_group <- apply(meta_df[, group.by, drop = FALSE], 1, paste, collapse = " - ")
  }
  
  p <- Seurat::DotPlot(
    seu,
    features = features, 
    group.by = "combined_group"
  ) + 
    scale_color_gradientn(colours = my_colors)
  plot_data <- p$data
  
  level_list <- lapply(rev(group.by), function(col) levels(meta_df[[col]]))
  group_combos <- expand.grid(level_list, stringsAsFactors = FALSE)
  ordered_group_names <- do.call(paste, c(rev(group_combos), sep = " - "))
  valid_group_names <- ordered_group_names[ordered_group_names %in% unique(plot_data$id)]
  plot_data$id <- factor(plot_data$id, levels = rev(valid_group_names))
  
  pct_range <- range(plot_data$pct.exp, na.rm = TRUE)
  if (pct_range[1] == pct_range[2]) {
    size_limits <- c(0, pct_range[2])
    if (size_limits[2] == 0) size_limits[2] <- 1
  } else {
    size_limits <- pct_range
  }
  size_breaks <- pretty(size_limits, n = 4)
  size_breaks <- size_breaks[size_breaks >= size_limits[1] & size_breaks <= size_limits[2]]
  if (length(size_breaks) == 0) size_breaks <- size_limits
  p <- ggplot(plot_data, aes(x = features.plot, y = id)) +
    geom_point(
      aes(size = pct.exp, fill = avg.exp.scaled),
      shape = 21,
      color = "black",
      stroke = 0.5
    ) +
    scale_fill_gradientn(colours = my_colors) +
    scale_size(range = c(1, 6), limits = size_limits, breaks = size_breaks) +
    guides(
      fill = guide_colorbar(order = 1),
      size = guide_legend(order = 2)
    ) +
    theme_classic() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      plot.title = element_text(hjust = 0.5, color = "black"),
      legend.text = element_text(color = "black", size = 14),
      legend.title = element_text(color = "black", size = 14),
      axis.title = element_text(color = "black", size = 14),
      axis.text.x = element_text(color = "black", size = 14,  vjust = 0.5, hjust = 1, angle = 90),
      axis.text.y = element_text(color = "black", size = 14),
      axis.line.x = element_line(color = "black", linewidth = 0.7, lineend = "square"),
      axis.line.y = element_line(color = "black", linewidth = 0.7, lineend = "square"),
      axis.ticks.x = element_line(color = "black", linewidth = 0.7),
      axis.ticks.y = element_line(color = "black", linewidth = 0.7),
      axis.ticks.length = unit(4, "pt")
    ) + 
    labs(x = NULL, y = NULL, size = "Percent Expressed", fill = "Average Expression")
  return(p)
}