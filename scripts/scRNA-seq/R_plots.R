setwd("~/Tigerrr07/scripts/scRNA-seq")
getwd()


source("../R_scripts/seurat_anndata_transformation.R")
source("../R_scripts/Dotplot.R")
source("../R_scripts/Heatmap.R")

seu <- qread("write/seurat.qs")

marker_genes <- list(
  `CD4 T`            = c("IL7R"),
  `CD8 T`            = c("CD8A", "CD8B"),
  B                 = c("CD79A", "MS4A1"),
  NK                = c("GNLY", "NKG7"),
  `FCGR3A+ Monocytes`= c("FCGR3A"),
  `CD14+ Monocytes`  = c("CD14", "LGALS3", "S100A8"),
  `Dendritic cell`   = c("FCER1A", "CST3"),
  Megakaryocyte      = c("PPBP")
)
marker_genes_vec <- unique(unlist(marker_genes, use.names = FALSE))
seu$cell_type <- factor(
  seu$cell_type,
  levels = c("CD4 T", "CD8 T", "B", "NK", "FCGR3A+ Monocytes",
             "CD14+ Monocytes", "Dendritic cell", "Megakaryocyte")
)

# Set lineage
new_labels <- c(
  "CD4 T"             = "Lymphoid",
  "CD8 T"             = "Lymphoid",
  "B"                 = "Lymphoid",
  "NK"                = "Lymphoid",
  "FCGR3A+ Monocytes" = "Myeloid",
  "CD14+ Monocytes"   = "Myeloid",
  "Dendritic cell"    = "Myeloid",
  "Megakaryocyte"     = "Myeloid"
)

seu$lineage <- unname(new_labels[as.character(seu$cell_type)])
seu$lineage <- factor(seu$lineage, levels = c("Lymphoid", "Myeloid"))


p <- scDotPlot(seu, features = marker_genes_vec, group.by = c("lineage", "cell_type"))
ggsave("plots/scDotPlot.pdf", p, width = 8, height = 5)   

pdf("plots/sc_lineage_ct_heatmap.pdf", width = 7, height = 5)
scHeatmap(seu, features = marker_genes, group.by = c("lineage", "cell_type"))
dev.off()