
## Nextflow installation
* [Tutorial link](https://nf-co.re/docs/get_started/environment_setup/nextflow#install-nextflow:~:text=handles%20dependencies%20automatically.-,Self%2Dinstalling%20package,-To%20install%20Nextflow)
* Make sure Java is installed as well (usally HPC servers has the Java installed)


## Create R environment

Select a working path first, then:
``` bash
# For DESeq2
renv::init()

install.packages("BiocManager", type = "source")
options(repos = BiocManager::repositories())
renv::install("DESeq2")
```

## Motif analysis using HOMER
* [HOMER installation](http://homer.ucsd.edu/homer/introduction/install.html#:~:text=Installing%20the%20basic%20HOMER%20software)
* Make sure Perl is installed as well (usally HPC servers has the Perl installed)

## Running workflow

Run the ATAC-seq analysis in the following order:

1. Run `run_atac.sh` to execute the Nextflow ATAC-seq pipeline and generate peak/count outputs.
2. Run `deseq2_peaks.R` to perform KO vs WT differential peak analysis on consensus peaks.
3. Run `load_deseq2_results.R` to summarize significant peaks and export up/down gene lists.
4. Run `run_motif.sh` to perform HOMER motif enrichment on the DESeq2 significant peak sets.

