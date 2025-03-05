#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=5G
#SBATCH --time=48:00:00
#SBATCH --partition=preemptable
#SBATCH --job-name=cellranger_${SAMPLE}
#SBATCH --error=/panfs/compbio/users/xran2/wen/22q/scRNA/logs/job.${SAMPLE}_%J_error.txt
#SBATCH --output=/panfs/compbio/users/xran2/wen/22q/scRNA/logs/job.${SAMPLE}_%J_out.txt
#SBATCH --mail-user=ximing.ran@emory.edu
#SBATCH --mail-type=END,FAIL

# Define variables for a single sample
SAMPLE="${SAMPLE}"  # Sample to test
OUTPUT_DIR="/panfs/compbio/users/xran2/wen/22q/scRNA/cellranger_out/${SAMPLE}"  # Output directory for the sample
FASTQ_DIR="/panfs/compbio/users/xran2/wen/22q/scRNA/raw/ZWE12492-121402"  # FASTQ file directory
TRANSCRIPTOME="/projects/compbio/users/xran2/tool/refdata-gex-GRCh38-2024-A" # Transcriptome reference

# Create output and logs directories if they don't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p /panfs/compbio/users/xran2/wen/22q/scRNA/logs

# Run Cell Ranger for the single sample
/projects/compbio/users/xran2/tool/cellranger-8.0.1/bin/cellranger count \
    --id=${SAMPLE} \
    --sample=${SAMPLE} \
    --transcriptome=${TRANSCRIPTOME} \
    --fastqs=${FASTQ_DIR} \
    --localcores=32 \
    --localmem=160 \
    --nosecondary \
    --disable-ui \
    --create-bam false

