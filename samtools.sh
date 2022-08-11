#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=samtools
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=trm76056@uga.edu
#SBATCH --time=100:00:00

cd $SLURM_SUBMIT_DIR

ml SAMtools/1.10-GCC-8.3.0
ml BCFtools/1.13-GCC-8.3.0

for x in bam_files/*.sorted.bam;
do
     $echo 'doing samtools sort'
     samtools sort -n -o $x.namesort.bam $x
     $echo 'doing samtools fixmate'
     samtools fixmate -rm -O bam --threads 4 $x.namesort.bam $x.fix.bam
     $echo 'doing samtools sort again'
     samtools sort -o $x.positionsort.bam $x.fix.bam
     $echo 'doing samtools markdup'
     samtools markdup -r --threads 4 $x.positionsort.bam $x.mark.bam
     $echo 'doing samtools mpileup'
     samtools mpileup -q 20 -Q 20 -C 50 -f AllReads.ORP.fasta.fai | bcftools call -mv > $x.raw.vcf
     $echo 'doing bcftools filter'
     bcftools filter -s LowQual -e '%QUAL<20 || DP>100' $x.raw.vcf > $x.flt.vcf
done
