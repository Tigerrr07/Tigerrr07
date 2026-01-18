## Scanpy to Seurat transformation


### Python script:
Required:
* `adata`
* `adata.layers['count']`, `adata.X`

``` python
from export_anndata import export_anndata

export_anndata(adata, export_path='export')
```

The export folder will look like this:
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

### R script

``` R
source("load_exported_anndata.R")
# Set python data path: data.path/export
data.path <- ""
export_path <- file.path(data.path, "export")

# Seurat object will be saved as : data.path/seurat.qs
save_path <- file.path(data.path, "seurat")

load_exported_anndata(export_path, save_path)
```