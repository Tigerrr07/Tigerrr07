# For detailed documentation, see: Scanpy2Seurat_transformation.md

library(Seurat)
library(Matrix)
library(qs)

load_exported_anndata <- function(export_path = "export", save_path = NULL) {
  # Read raw counts matrix
  raw_counts <- readMM(paste0(export_path, "/raw_counts.mtx"))
  # Read normalized expression matrix
  norm_expr <- readMM(paste0(export_path, "/norm_expr.mtx"))
  
  # Read gene and cell names
  genes <- read.csv(paste0(export_path, "/genes.csv"), header = FALSE, stringsAsFactors = FALSE)$V1
  cells <- read.csv(paste0(export_path, "/cells.csv"), header = FALSE, stringsAsFactors = FALSE)$V1
  
  # Read metadata
  metadata <- read.csv(paste0(export_path, "/metadata.csv"), row.names = 1)
  
  # Transpose raw counts and set row/col names
  rownames(raw_counts) <- genes
  colnames(raw_counts) <- cells
  rownames(norm_expr) <- genes
  colnames(norm_expr) <- cells
  
  # Create Seurat object
  seu <- CreateSeuratObject(counts = raw_counts, meta.data = metadata)
  # Set normalized data using SetAssayData
  seu <- SetAssayData(seu, layer = "data", new.data = norm_expr)
  
  # Load UMAP coordinates if exists
  umap_file <- paste0(export_path, "/umap_coordinates.csv")
  if (file.exists(umap_file)) {
    umap_coords <- read.csv(umap_file, row.names = 1)
    seu@reductions$umap <- Seurat::CreateDimReducObject(
      embeddings = as.matrix(umap_coords),
      key = "UMAP_",
      assay = DefaultAssay(seu)
    )
  }
  
  tsne_file <- paste0(export_path, "/tsne_coordinates.csv")
  if (file.exists(tsne_file)) {
    tsne_coords <- read.csv(tsne_file, row.names = 1)
    seu@reductions$tsne <- Seurat::CreateDimReducObject(
      embeddings = as.matrix(tsne_coords),
      key = "TSNE_",
      assay = DefaultAssay(seu)
    )
  }
  
  # Save with qs if requested
  if (!is.null(save_path)) {
    qs_path <- paste0(save_path, ".qs")
    qsave(seu, qs_path)
    cat("Seurat object saved to", qs_path, "\n")
  }

  return(seu)
}

