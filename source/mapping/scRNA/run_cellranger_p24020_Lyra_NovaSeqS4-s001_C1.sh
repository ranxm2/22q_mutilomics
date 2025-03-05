#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=5G
#SBATCH --time=48:00:00
#SBATCH --partition=preemptable
#SBATCH --job-name=cellranger_p24020_Lyra_NovaSeqS4-s001_C1
#SBATCH --error=/panfs/compbio/users/xran2/wen/22q/scRNA/logs/job.p24020_Lyra_NovaSeqS4-s001_C1_%J_error.txt
#SBATCH --output=/panfs/compbio/users/xran2/wen/22q/scRNA/logs/job.p24020_Lyra_NovaSeqS4-s001_C1_%J_out.txt
#SBATCH --mail-user=ximing.ran@emory.edu
#SBATCH --mail-type=END,FAIL

# Define variables for this sample
SAMPLE="p24020_Lyra_NovaSeqS4-s001_C1"
OUTPUT_DIR="/panfs/compbio/users/xran2/wen/22q/scRNA/cellranger_out/p24020_Lyra_NovaSeqS4-s001_C1"
FASTQ_DIR="/panfs/compbio/users/xran2/wen/22q/scRNA/raw/ZWE12492-121402"
TRANSCRIPTOME="/projects/compbio/users/xran2/tool/refdata-gex-GRCh38-2024-A"

# Create output and logs directories
mkdir -p ${OUTPUT_DIR}
mkdir -p /panfs/compbio/users/xran2/wen/22q/scRNA/logs

# Run Cell Ranger for the sample
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
