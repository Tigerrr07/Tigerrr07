library(Seurat)
library(ggplot2)
# library(viridis)


scDotPlot <- function(seu, features, group.by, my_colors = c("#4A90E2", "#D9D9D9", "#D0021B")) {
  p <- Seurat::DotPlot(
    seu,
    features = features, 
    group.by = group.by
  ) + 
    scale_color_gradientn(colours = my_colors)
  
  p <- ggplot(p$data, aes(x = features.plot, y = id)) +
    geom_point(
      aes(size = pct.exp, fill = avg.exp.scaled),
      shape = 21,
      color = "black",
      stroke = 0.5
    ) +
    scale_fill_gradientn(colours = my_colors) +
    scale_size(range = c(1, 6)) +
    theme_classic() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, color = "black"),
      axis.title = element_text(color = "black", size = 12),
      axis.text.y = element_text(color = "black", size = 12),
      axis.text.x = element_text(color = "black", size = 12,  vjust = 0.5, hjust = 1, angle = 90),
      legend.text = element_text(color = "black", size = 12),
      legend.title = element_text(color = "black", size = 12),
      panel.grid = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 1, lineend = "square"),
      axis.line.y = element_line(color = "black", linewidth = 1, lineend = "square"),
      axis.ticks.x = element_line(color = "black", linewidth = 1),
      axis.ticks.y = element_line(color = "black", linewidth = 1),
      axis.ticks.length = unit(4, "pt")
    ) + 
    labs(x = NULL, y = NULL, size = "Percent Expressed", fill = "Average Expression")
  return(p)
}