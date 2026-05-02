# For detailed documentation, see: Scanpy2Seurat_transformation.md

library(Seurat)
library(Matrix)
library(qs)

load_exported_anndata <- function(export_path = "export_anndata", save_path = NULL) {
  # Read raw counts matrix
  raw_path <- paste0(export_path, "/raw_counts.mtx")
  raw_gz <- paste0(raw_path, ".gz")
  raw_counts <- if (file.exists(raw_gz)) {
    readMM(gzfile(raw_gz, "rt"))
  } else {
    readMM(raw_path)
  }
  # Read normalized expression matrix
  norm_path <- paste0(export_path, "/norm_expr.mtx")
  norm_gz <- paste0(norm_path, ".gz")
  norm_expr <- if (file.exists(norm_gz)) {
    readMM(gzfile(norm_gz, "rt"))
  } else {
    readMM(norm_path)
  }

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
  seu <- CreateSeuratObject(counts = raw_counts, meta.data = metadata, assay = "RNA")
  DefaultAssay(seu) <- "RNA"

  counts_mat <- Seurat::GetAssayData(seu, layer = "counts")

  # Seurat may sanitize feature names in counts (e.g. '_' -> '-').
  # Force exact dimnames from counts to prevent feature union on SetAssayData.
  rownames(norm_expr) <- rownames(counts_mat)
  colnames(norm_expr) <- colnames(counts_mat)

  # Set normalized data using SetAssayData on the RNA assay
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

# Export a Seurat object to files that can be loaded into Scanpy/AnnData
export_seurat <- function(
  seu,
  export_path = "export_seurat",
  assay = "RNA",
  counts_layer = "counts",
  data_layer = "data",
  save_umap = TRUE,
  save_tsne = TRUE,
  coerce_sparse = TRUE,
  compress_mtx = TRUE
) {
  dir.create(export_path, showWarnings = FALSE, recursive = TRUE)

  gzip_file <- function(src, dest) {
    in_con <- file(src, "rb")
    out_con <- gzfile(dest, "wb")
    on.exit({
      try(close(in_con), silent = TRUE)
      try(close(out_con), silent = TRUE)
    }, add = TRUE)
    repeat {
      buf <- readBin(in_con, what = "raw", n = 1e6)
      if (!length(buf)) break
      writeBin(buf, out_con)
    }
  }

  # Helper to support Seurat v4 (slot) and v5 (layer)
  get_assay_matrix <- function(obj, assay_name, layer_name, slot_name) {
    if ("layer" %in% names(formals(Seurat::GetAssayData))) {
      return(Seurat::GetAssayData(obj, assay = assay_name, layer = layer_name))
    }
    return(Seurat::GetAssayData(obj, assay = assay_name, layer = slot_name))
  }

  counts_mat <- get_assay_matrix(seu, assay, counts_layer, "counts")
  data_mat <- get_assay_matrix(seu, assay, data_layer, "data")

  if (coerce_sparse) {
    counts_mat <- Matrix::drop0(as(counts_mat, "CsparseMatrix"))
    data_mat <- Matrix::drop0(as(data_mat, "CsparseMatrix"))
  }

  # Write matrices (features x cells), optionally gzipped
  counts_path <- paste0(export_path, "/raw_counts.mtx")
  data_path <- paste0(export_path, "/norm_expr.mtx")
  if (compress_mtx) {
    tmp_counts <- tempfile(pattern = "raw_counts_", fileext = ".mtx")
    tmp_data <- tempfile(pattern = "norm_expr_", fileext = ".mtx")
    writeMM(counts_mat, tmp_counts)
    writeMM(data_mat, tmp_data)
    gzip_file(tmp_counts, paste0(counts_path, ".gz"))
    gzip_file(tmp_data, paste0(data_path, ".gz"))
    unlink(c(tmp_counts, tmp_data))
  } else {
    writeMM(counts_mat, counts_path)
    writeMM(data_mat, data_path)
  }

  # Genes and cells
  genes <- rownames(counts_mat)
  cells <- colnames(counts_mat)
  write.table(
    genes,
    paste0(export_path, "/genes.csv"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = ","
  )
  write.table(
    cells,
    paste0(export_path, "/cells.csv"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = ","
  )

  # Metadata
  write.csv(seu@meta.data, paste0(export_path, "/metadata.csv"))

  # UMAP coordinates if present
  if (save_umap && "umap" %in% names(seu@reductions)) {
    umap_coords <- Embeddings(seu, reduction = "umap")
    write.csv(umap_coords, paste0(export_path, "/umap_coordinates.csv"))
  }

  # TSNE coordinates if present
  if (save_tsne && "tsne" %in% names(seu@reductions)) {
    tsne_coords <- Embeddings(seu, reduction = "tsne")
    write.csv(tsne_coords, paste0(export_path, "/tsne_coordinates.csv"))
  }

  cat("All files exported to", export_path, "\n")
}
