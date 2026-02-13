#!/usr/bin/env python3

import argparse
import os
import sys
import numpy as np
sys.path.append("/users/PAS2148/tiger/Tigerrr07/scripts/Python_scripts")
from seurat_anndata_transformation import load_anndata_from_seurat_export


def main():
    parser = argparse.ArgumentParser(
        description="Test loader for Seurat->AnnData export"
    )
    parser.add_argument(
        "--export-path",
        default="export",
        help="Path containing raw_counts.mtx.gz, norm_expr.mtx.gz, genes.csv, cells.csv, metadata.csv",
    )
    args = parser.parse_args()

    export_path = args.export_path
    if not os.path.isdir(export_path):
        print(f"ERROR: export path not found: {export_path}")
        return 2

    adata = load_anndata_from_seurat_export(export_path)

    if adata.n_obs == 0 or adata.n_vars == 0:
        print("ERROR: AnnData has zero cells or genes")
        return 3

    if "counts" not in adata.layers:
        print("ERROR: missing counts layer")
        return 4

    if adata.layers["counts"].shape != adata.shape:
        print("ERROR: counts layer shape mismatch")
        return 5

    if not np.all(adata.obs_names == adata.obs.index):
        print("ERROR: obs_names and obs index mismatch")
        return 6

    if adata.obsm.get("X_umap") is not None and adata.obsm["X_umap"].shape[0] != adata.n_obs:
        print("ERROR: UMAP rows do not match number of cells")
        return 7

    if adata.obsm.get("X_tsne") is not None and adata.obsm["X_tsne"].shape[0] != adata.n_obs:
        print("ERROR: TSNE rows do not match number of cells")
        return 8

    print("OK")
    print(f"cells={adata.n_obs} genes={adata.n_vars}")
    print(f"X shape={adata.X.shape} counts shape={adata.layers['counts'].shape}")
    print(f"UMAP={'X_umap' in adata.obsm} TSNE={'X_tsne' in adata.obsm}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
