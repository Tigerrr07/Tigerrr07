#!/usr/bin/bash
#SBATCH --job-name bulk_RNAseq_job
#SBATCH --account PAS2148
#SBATCH --time=12:00:00
#SBATCH --nodes=1 
#SBATCH --mem=256GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chen.13589@osu.edu

source ~/.bashrc

conda activate rnaseq_env # for kallisto, fastp

# Set working directory
cd /path/to/workdir

data_dir="/path/to/data"  # change to your data path
out_dir="/path/to/out" # change to your path, for saving data

# Download transcripts file from GENCODE
ref_fa="/path/to/genome/gencode.vM37.transcripts.fa" # change to your path, transcripts files
ref_index="/path/to/genome/gencode.vM37.transcripts.index" # change to your path, best at the same directory as ref_fa

if [[ -d $out_dir ]]
then
echo "Output directory exist."
else
    mkdir -p $out_dir
fi

# build index
if [[ -e $ref_index ]]
then
    echo "Ref directory exist."
else
    kallisto index -i $ref_index $ref_fa
fi


input_file="metadata.txt"

tail -n +2 "$input_file" | while IFS=$'\t' read -r sample_name condition data_name; do
    echo "Processing ${sample_name}..."
    fastq1=${data_dir}/${data_name}_R1_001.fastq.gz
    fastq2=${data_dir}/${data_name}_R2_001.fastq.gz
    echo "Processing $fastq1 and $fastq2"

    if [[ -d ${out_dir}/${sample_name} ]]; then
        echo "Directory exist."
    else
        echo "Not exist."
        mkdir ${out_dir}/${sample_name}
    fi

    # Quality control and adapter trimming
    fastp -w 16 -i $fastq1 -I $fastq2 -o ${out_dir}/${sample_name}/${sample_name}_R1.fastq.gz -O ${out_dir}/${sample_name}/${sample_name}_R2.fastq.gz -h ${out_dir}/${sample_name}/$sample_name.html -j ${out_dir}/${sample_name}/$sample_name.fastp.json

    # Quantification
    kallisto quant -i $ref_index -o ${out_dir}/${sample_name} --bias ${out_dir}/${sample_name}/${sample_name}_R1.fastq.gz ${out_dir}/${sample_name}/${sample_name}_R2.fastq.gz

    rm ${out_dir}/${sample_name}/*.fastq.gz
done

