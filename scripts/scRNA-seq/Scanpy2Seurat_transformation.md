## Scanpy to Seurat transformation


### Python script:
Required:
* `adata`
* `adata.layers['count']`, `adata.X`

``` python
from scipy.io import mmwrite
import scanpy as sc
import os
 
def export_anndata(adata, export_path="export"):
    os.makedirs(export_path, exist_ok=True)
    
    # store raw counts in layers
    mmwrite(f'{export_path}/raw_counts.mtx', adata.layers['counts'].T)
    mmwrite(f'{export_path}/norm_expr.mtx', adata.X.T)
    
    pd.Series(adata.var.index).to_csv(f'{export_path}/genes.csv', header=False, index=False)
    pd.Series(adata.obs.index).to_csv(f'{export_path}/cells.csv', header=False, index=False)
    adata.obs.to_csv(f'{export_path}/metadata.csv')
    
    if 'X_umap' in adata.obsm:
        umap_df = pd.DataFrame(
            adata.obsm['X_umap'], 
            columns=['UMAP_1', 'UMAP_2'], 
            index=adata.obs.index
        )
        umap_df.to_csv(f'{export_path}/umap_coordinates.csv')
        print(f"UMAP coordinates saved to {export_path}/umap_coordinates.csv")
    else:
        print("Warning: No UMAP coordinates found in adata.obsm")

    if 'X_tsne' in adata.obsm:
        umap_df = pd.DataFrame(
            adata.obsm['X_tsne'], 
            columns=['TSNE_1', 'TSNE_2'], 
            index=adata.obs.index
        )
        umap_df.to_csv(f'{export_path}/tsne_coordinates.csv')
        print(f"T-STNE coordinates saved to {export_path}/tsne_coordinates.csv")
    else:
        print("Warning: No TSNE coordinates found in adata.obsm")
        
    print(f"All files exported to {export_path}")

export_anndata(adata, export_path="export") # Save into current folder
```

The folder will look like this:
```
export
├── cells.csv
├── genes.csv
├── metadata.csv
├── norm_expr.mtx
├── raw_counts.mtx
├── tsne_coordinates.csv
└── umap_coordinates.csv
```

# R script

``` R
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
  seurat_obj <- CreateSeuratObject(counts = raw_counts, meta.data = metadata)
  # Set normalized data using SetAssayData
  seurat_obj <- SetAssayData(seurat_obj, layer = "data", new.data = norm_expr)
  
  # Load UMAP coordinates if exists
  umap_file <- paste0(export_path, "/umap_coordinates.csv")
  if (file.exists(umap_file)) {
    umap_coords <- read.csv(umap_file, row.names = 1)
    seurat_obj@reductions$umap <- Seurat::CreateDimReducObject(
      embeddings = as.matrix(umap_coords),
      key = "UMAP_",
      assay = DefaultAssay(seurat_obj)
    )
  }
  
  tsne_file <- paste0(export_path, "/tsne_coordinates.csv")
  if (file.exists(tsne_file)) {
    tsne_coords <- read.csv(tsne_file, row.names = 1)
    seurat_obj@reductions$tsne <- Seurat::CreateDimReducObject(
      embeddings = as.matrix(tsne_coords),
      key = "TSNE_",
      assay = DefaultAssay(seurat_obj)
    )
  }
  
  # Save with qs if requested
  if (!is.null(save_path)) {
    qs_path <- paste0(save_path, ".qs")
    qsave(seurat_obj, qs_path)
    cat("Seurat object saved to", qs_path, "\n")
  }

  return(seurat_obj)
}


# Set python data path: data.path/export
data.path <- ""
export_path <- file.path(data.path, "export")

# Seurat object will be saved as : data.path/seurat.qs
save_path <- file.path(data.path, "seurat")

load_exported_anndata(export_path, save_path)
```

