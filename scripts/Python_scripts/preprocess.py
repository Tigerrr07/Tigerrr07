import os
os.environ["RENV_CONFIG_AUTOLOADER_ENABLED"] = "FALSE"
import logging
# R-Python interop imports
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.packages import importr
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

R_HOME = os.path.expanduser("/users/PAS2148/tiger/miniconda3/envs/RNA_pipline/lib/R")
os.environ["R_HOME"] = R_HOME

print("R_HOME:", os.environ["R_HOME"])

# Activate rpy2 converters
rcb.logger.setLevel(logging.ERROR)  # Suppress Rpy2 warnings
ro.pandas2ri.activate()  # Enable Pandas↔R conversion
anndata2ri.activate()  # Enable AnnData↔R conversion


# Print R version to verify connection
print(ro.r("version"))

# Activate automatic conversion between R and Python objects
numpy2ri.activate()
pandas2ri.activate()

# Import required R packages
sce_pkg = importr("SingleCellExperiment")
scdblfinder = importr("scDblFinder")
BiocParallel = importr("BiocParallel")


def run_scdblfinder(adata, seed):
    data = adata.X  # (cells x genes)
    genes = adata.var_names.tolist()
    cells = adata.obs_names.tolist()


    # Transpose to genes x cells for R compatibility
    data = data.T

    # Pass variables to R environment
    ro.globalenv["data_mat"] = data
    ro.globalenv["genes"] = StrVector(genes)
    ro.globalenv["cells"] = StrVector(cells)
    ro.r(f"set.seed({seed})")

    # Run scDblFinder in R
    ro.r('''
    library(SingleCellExperiment)
    library(scDblFinder)

    rownames(data_mat) <- genes
    colnames(data_mat) <- cells

    sce <- SingleCellExperiment(list(counts = data_mat))
    sce <- scDblFinder(sce)

    doublet_score <- sce$scDblFinder.score
    doublet_class <- sce$scDblFinder.class
    ''')

    # Retrieve output
    score = np.array(ro.r("doublet_score"))
    label = np.array(ro.r("as.character(doublet_class)"))

    result = pd.DataFrame({
        "scDblFinder_score": score,
        "scDblFinder_class": label
    }, index=cells)

    return result


def pre_process(adata, lo_lim_nfeature=500, up_lim_nfeature=6000, up_lim_mt=10, filter_doublets=False, seed=0):
    adata.var_names_make_unique()
    print(f"Before Pre-processing: {adata.n_obs}")
    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.upper().str.startswith(("HBA", "HBB"))

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], percent_top=[20], log1p=True, inplace=True)

    with plt.rc_context({'font.size': 20}):
        sc.pl.violin(
        adata,
        ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_hb"],
        jitter=0.4,
        multi_panel=True,
        size=2
        )



    with plt.rc_context({'font.size': 12}):
        fig, axs = plt.subplots(1, 2, figsize=(6.5, 3), layout="constrained")
        sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", show=False, ax=axs[0])
        sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", show=False, ax=axs[1])
        plt.show()

    if lo_lim_nfeature is not None:
        adata = adata[adata.obs['n_genes_by_counts'] > lo_lim_nfeature].copy()
        
    if up_lim_nfeature is not None:
        adata = adata[adata.obs['n_genes_by_counts'] < up_lim_nfeature].copy()

    if up_lim_mt is not None:
        adata = adata[adata.obs['pct_counts_mt'] < up_lim_mt].copy()

    print(f"After filtering of low quality cells: {adata.n_obs}")


    with plt.rc_context({'font.size': 20}):
        sc.pl.violin(
        adata,
        ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_hb"],
        jitter=0.4,
        multi_panel=True,
        size=2
        )
    
    with plt.rc_context({'font.size': 12}):
        fig, axs = plt.subplots(1, 2, figsize=(6.5, 3), layout="constrained")
        sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", show=False, ax=axs[0])
        sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", show=False, ax=axs[1])
        plt.show()

    df_doublets = run_scdblfinder(adata, seed=seed)
    adata.obs["scDblFinder_score"] = df_doublets["scDblFinder_score"]
    adata.obs["scDblFinder_class"] = df_doublets["scDblFinder_class"]
    print(adata.obs["scDblFinder_class"].value_counts())

    if filter_doublets:
        adata = adata[adata.obs['scDblFinder_class'] == "singlet"].copy()
    
        print(f"After filtering dobulet: {adata.n_obs}")
        with plt.rc_context({'font.size': 20}):
            sc.pl.violin(
            adata,
            ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_hb"],
            jitter=0.4,
            multi_panel=True,
            size=2
            )
        
    with plt.rc_context({'font.size': 12}):
        fig, axs = plt.subplots(1, 2, figsize=(6.5, 3), layout="constrained")
        sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", show=False, ax=axs[0])
        sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", show=False, ax=axs[1])
        plt.show()        

    print(f"After Pre-processing: {adata.n_obs}")

    return adata