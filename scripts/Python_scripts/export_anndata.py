# For detailed documentation, see: Scanpy2Seurat_transformation.md

from scipy.io import mmwrite
import scanpy as sc
import pandas as pd
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
        tsne_df = pd.DataFrame(
            adata.obsm['X_tsne'], 
            columns=['TSNE_1', 'TSNE_2'], 
            index=adata.obs.index
        )
        tsne_df.to_csv(f'{export_path}/tsne_coordinates.csv')
        print(f"T-SNE coordinates saved to {export_path}/tsne_coordinates.csv")
    else:
        print("Warning: No TSNE coordinates found in adata.obsm")
        
    print(f"All files exported to {export_path}")

