#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=salmon
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=100:00:00

cd $SLURM_SUBMIT_DIR

ml Salmon/1.1.0-gompi-2019b
#this indexes the transcriptome
#salmon index -t AllReads.ORP.fasta -i AllReads.salmon_index

#this will go through each RNA seq read and produce quant files
#declare -a my_list=("GA1-30","GA3-30","GA4-30","GA5-30","GA6-20","GA7-20","GA8-20","GA9-20","GA10-20","GAC1","GAC2","GAC3","GAC4","GAC5","MA12-30","MA13-30","MA15-30","MA16-20","MA17-20","MA18-20","MA20-20","MAC1","MAC2","MAC3","MAC5")

#declare -a my_list=("GA1-30" "GA3-30" "GA4-30" "GA5-30" "GA6-20" "GA7-20" "GA8-20" "GA9-20" "GA10-20" "GAC1" "GAC2" "GAC3" "GAC4" "GAC5" "MA12-30" "MA13-30" "MA15-30" "MA16-20" "MA17-20" "MA18-20" "MA20-20" "MAC1" "MAC2" "MAC3" "MAC5") 

#for samp in "${my_list[@]}"
#do
#echo "Processing $samp"
#salmon quant -i AllReads.salmon_index -l A \
#    -1 ./files/trim${samp}_1.fastq.gz \
#    -2 ./files/trim${samp}_2.fastq.gz \
#    -p 8 --validateMappings -o quants/${samp}_quant
#done

salmon quant -i AllReads.salmon_index -l A \
    -1 ./files/trimMA14-30_1.fastq.gz \
    -2 ./files/trimMA14-30_2.fastq.gz \
    -p 8 --validateMappings -o quants/MA14-30_quant    
