#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=MA20-20
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=trm76056@uga.edu

cd $SLURM_SUBMIT_DIR

ml Bowtie2/2.4.1-GCC-8.3.0

bowtie2 -x AllReads -1 files/trimMA20-20_1.fastq.gz -2 files/trimMA20-20_2.fastq.gz -S MA20-20.sam
