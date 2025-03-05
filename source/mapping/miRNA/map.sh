#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=10G
#SBATCH --time=48:00:00
#SBATCH --partition=week-long-cpu,month-long-cpu,preemptable
#SBATCH --job-name=miRNA
#SBATCH --error=/projects/compbio/users/xran2/wen/22q/smallRNA/logs/job.%J_error.txt
#SBATCH --output=/projects/compbio/users/xran2/wen/22q/smallRNA/logs/job.%J_out.txt
#SBATCH --mail-user=ximing.ran@emory.edu
#SBATCH --mail-type=END,FAIL

echo "Job ID: $SLURM_JOB_ID"
echo "Node list: $SLURM_JOB_NODELIST"
echo "Submit directory: $SLURM_SUBMIT_DIR"
echo "Submit host: $SLURM_SUBMIT_HOST"
echo "Queue: $SLURM_JOB_PARTITION"

cd /projects/compbio/users/xran2/wen/22q/smallRNA

conda activate env_nf

echo "Running nextflow..."
module load singularity

nextflow run nf-core/smrnaseq \
  -r 2.3.0 \
  -profile singularity \
  --input fastq_samples.csv \
  --genome 'GRCh38' \
  --mirtrace_species 'hsa' \
  --protocol 'illumina' \
  --outdir result
