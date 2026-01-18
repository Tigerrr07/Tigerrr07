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

## scRNA-seq & ST data analysis, Python and R, conda abd renv envrionment
Manually installation, conda environment:
``` bash
conda create -n scRNA_pipline python=3.10 -y
conda activate scRNA_pipline
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


renv environment:
``` R
renv::init(bare = TRUE) # Empty project library
getOption("repos")
options(repos = c(CRAN = "https://cran.r-project.org"))
install.packages("BiocManager", type = "source")
options(repos = BiocManager::repositories())

renv::install("Seurat@5.4.0")
renv::install("ragg")
renv::install("qs@0.27.3") # Save object
renv::install("rmarkdown")
renv::install("markdown")
renv::install("ComplexHeatmap")


renv::snapshot()
```

Bash command in OSC environment:
``` bash
moudle load gcc/12.3.0
module load R/4.4.0
```


## Prolem when install Seurat
### Problem 1: igraph need glbk
**Problem**: igraph need glbk: https://r.igraph.org/articles/installation-troubleshooting.html#libglpk-so-40-cannot-open-shared-object-file-no-such-file-or-directory
``` bash
The following required system packages are not installed:
- glpk-devel  [required by igraph]
The R packages depending on these system packages may fail to install.

An administrator can install these packages with:
- sudo dnf install glpk-devel
```

**Solution**: download glbk manually，add into environment
``` bash
wget https://ftp.gnu.org/gnu/glpk/glpk-5.0.tar.gz
tar -xzf glpk-5.0.tar.gz
cd glpk-5.0
./configure --prefix=$HOME/.local
make -j8
make install
echo 'export LD_LIBRARY_PATH=$HOME/.local/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
```

### Problem 2: igraph dynamic library error (libicui18n)
**Problem**: `igraph` installation failed due to a `libicui18n.so.73` missing error. This was caused by an environment conflict where R was detecting incompatible libraries from `miniconda` in the `PATH` instead of the correct system libraries.

```R
Error: package or namespace load failed for ‘igraph’ in dyn.load(file, DLLpath = DLLpath, ...):
 unable to load shared object '.../igraph.so':
  libicui18n.so.73: cannot open shared object file: No such file or directory
```
**Solution**: Force clean the `PATH` environment variable within the R session to remove `miniconda` references, then reinstall.

```R

# 1. Remove all 'miniconda' paths from the environment
Sys.setenv(PATH = paste(grep("miniconda", strsplit(Sys.getenv("PATH"), ":")[[1]], 
                             value = TRUE, invert = TRUE), collapse = ":"))

# 2. Install
renv::install("Seurat@5.0.3")
```