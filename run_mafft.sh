#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=mafft
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=trm76056@uga.edu
#SBATCH --time=100:00:00
#SBATCH --mem=100G

cd /scratch/trm76056/Geuk21/metabolic_genes

ml MAFFT/7.470-GCC-8.3.0-with-extensions

#mafft --thread 4 --auto --nuc all_lysin.fa > all_lysin_mafft.fa

for x in noMA*
do
    echo $x
    mafft --thread 4 --auto --nuc $x > $x.mafft.fa
done
