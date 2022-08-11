#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=vcf
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G

cd $SLURM_SUBMIT_DIR

module load R/3.6.2-foss-2019b

Rscript vcf.R
