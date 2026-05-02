#!/usr/bin/bash
#SBATCH --account XXXX
#SBATCH --job-name=atac
#SBATCH --time=20:00:00
#SBATCH --mem=128G
#SBATCH --output=slurm-%j.out

cd /path/to/your/project

mkdir -p results

echo "Started at: $(date)"

# For mouse, read length is 150
nextflow run nf-core/atacseq \
    --input ./samplesheet.csv \
    --outdir ./results \
    --genome mm10 \
    --read_length 150 \
    -profile singularity \
    -resume

echo "Finished at: $(date)"

