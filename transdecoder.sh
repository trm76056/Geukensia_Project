#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=MA14-30
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=150:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=trm76056@uga.edu

cd /scratch/trm76056/Geuk21/sample_assemblies

module load TransDecoder/5.5.0-foss-2019b-Perl-5.30.0
module load BLAST+/2.11.0-gompi-2019b
module load HMMER/3.2.1-GCC-8.3.0

species="MA14-30"

TransDecoder.LongOrfs -t $species.ORP.fasta

blastp -query $species.ORP.fasta.transdecoder_dir/longest_orfs.pep \
     -db uniprot_sprot.fasta -max_target_seqs 1 \
     -outfmt 6 -evalue 1e-5 -num_threads 10 > $species.blastp.outfmt6

hmmpress -f /scratch/trm76056/Geuk21/Pfam-A.hmm

hmmscan --cpu 10 --domtblout $species.pfam.domtblout /scratch/trm76056/Geuk21/Pfam-A.hmm $species.ORP.fasta.transdecoder_dir/longest_orfs.pep

TransDecoder.Predict -t $species.ORP.fasta --retain_pfam_hits $species.pfam.domtblout --retain_blastp_hits $species.blastp.outfmt6