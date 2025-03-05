#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=5G
#SBATCH --time=48:00:00
#SBATCH --partition=preemptable
#SBATCH --job-name=cellranger_job
#SBATCH --error=/home/xran2/cellranger_out/logs/job.%J_error.txt
#SBATCH --output=/home/xran2/cellranger_out/logs/job.%J_out.txt
#SBATCH --mail-user=ximing.ran@emory.edu
#SBATCH --mail-type=END,FAIL

# Load necessary modules
# module load some-module if needed

# Define variables
SAMPLE=$1  # Sample name passed as an argument
OUTPUT_DIR=/home/xran2/cellranger_out/${SAMPLE}  # Output directory for the sample
FASTQ_DIR=/panfs/compbio/users/xran2/wen/22q/scRNA/raw/ZWE12492-121402  # FASTQ file directory
TRANSCRIPTOME=/projects/compbio/users/xran2/refdata-gex-GRCh38-2020-A  # Transcriptome reference

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p /home/xran2/cellranger_out/logs

# Run Cell Ranger
/projects/compbio/users/xran2/tool/cellranger/cellranger-8.0.1/cellranger count \
    --id=${SAMPLE} \
    --sample=${SAMPLE} \
    --transcriptome=${TRANSCRIPTOME} \
    --fastqs=${FASTQ_DIR} \
    --localcores=32 \
    --localmem=160 \
    --nosecondary \
    --disable-ui


