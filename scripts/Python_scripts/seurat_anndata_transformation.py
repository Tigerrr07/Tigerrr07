# For detailed documentation, see: Scanpy2Seurat_transformation.md

from scipy.io import mmwrite, mmread
import scanpy as sc
import pandas as pd
import os
import gzip
import numpy as np
 
def export_anndata(adata, export_path="export", compress_mtx=True):
    os.makedirs(export_path, exist_ok=True)
    
    # store raw counts in layers
    raw_path = f'{export_path}/raw_counts.mtx'
    norm_path = f'{export_path}/norm_expr.mtx'
    if compress_mtx:
        with gzip.open(raw_path + ".gz", "wb") as fh:
            mmwrite(fh, adata.layers['counts'].T)
        with gzip.open(norm_path + ".gz", "wb") as fh:
            mmwrite(fh, adata.X.T)
    else:
        mmwrite(raw_path, adata.layers['counts'].T)
        mmwrite(norm_path, adata.X.T)
    
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


def _read_mtx(path_base):
    if os.path.exists(path_base + ".gz"):
        with gzip.open(path_base + ".gz", "rb") as fh:
            return mmread(fh)
    return mmread(path_base)


def load_anndata_from_seurat_export(export_path="export"):
    counts_path = os.path.join(export_path, "raw_counts.mtx")
    data_path = os.path.join(export_path, "norm_expr.mtx")

    counts = _read_mtx(counts_path).T.tocsr()
    print(counts.shape)
    data = _read_mtx(data_path).T.tocsr()
    print(data.shape)

    genes = pd.read_csv(os.path.join(export_path, "genes.csv"), header=None)[0].astype(str).to_numpy()
    cells = pd.read_csv(os.path.join(export_path, "cells.csv"), header=None)[0].astype(str).to_numpy()

    obs = pd.read_csv(os.path.join(export_path, "metadata.csv"), index_col=0)
    obs.index = obs.index.astype(str)

    adata = sc.AnnData(X=data, obs=obs)
    adata.var_names = genes
    adata.obs_names = cells
    adata.layers["counts"] = counts

    umap_file = os.path.join(export_path, "umap_coordinates.csv")
    if os.path.exists(umap_file):
        umap_df = pd.read_csv(umap_file, index_col=0)
        adata.obsm["X_umap"] = umap_df.loc[adata.obs_names].to_numpy()

    tsne_file = os.path.join(export_path, "tsne_coordinates.csv")
    if os.path.exists(tsne_file):
        tsne_df = pd.read_csv(tsne_file, index_col=0)
        adata.obsm["X_tsne"] = tsne_df.loc[adata.obs_names].to_numpy()

    return adata
