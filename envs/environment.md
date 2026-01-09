## bulk RNA-seq analysis, Python and R, conda and renv environment

Conda command:
``` bash
conda create -n rnaseq_env -c bioconda kallisto fastp
conda activate rnaseq_env
```

R command:
``` R
renv::init(bare = TRUE) # Empty project library
getOption("repos")
options(repos = c(CRAN = "https://cran.r-project.org"))
install.packages("BiocManager", type = "source")
options(repos = BiocManager::repositories())

renv::install("tximport")
renv::install("DESeq2")
renv::install("GenomicFeatures")
renv::install("txdbmaker")

renv::install("ComplexHeatmap")
renv::install("rhdf5")
renv::install("extrafont")

renv::install("showtext")
renv::install("edgeR")
renv::install("org.Mm.eg.db")
renv::install("clusterProfiler")

renv::snapshot()
```

## 
## scRNA-seq & ST data analysis, Python and R, conda envrionment
Manually installation:
``` bash
conda create -n scRNA_pipline python=3.10

# Python packages
pip install ipykernel # For jupyter notebook
pip install ipywidgets
python -m ipykernel install --user --name scRNA_pipline

pip install 'scanpy[leiden]'
pip install harmonypy # batch correction
pip install decoupler
pip install openpyxl # For excel
pip install gseapy
pip install squidpy
pip install bbknn # optional
pip install scvi-tools # optional


# R packages
conda install -c conda-forge r-base=4.3.3 -y
conda install anndata2ri -y
conda install -c conda-forge r-seurat==5.0.3 r-seuratobject=5.0.2 r-sp==2.1_3 r-matrix==1.6_5 -y
conda install bioconductor-scdblfinder=1.16.0 -y

```

Install from enviroment file:
```bash
conda env create -f scRNA_analysis.yml -n bio_analysis_env # Replace with any name
```

