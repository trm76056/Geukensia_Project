#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=blast
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=trm76056@uga.edu
#SBATCH --time=160:00:00

cd $SLURM_SUBMIT_DIR

module load BLAST+/2.11.0-gompi-2020b

#species="MAC1"
#list=("GAC1","GAC2","GAC3","GAC4","GAC5","MAC1","MAC2","MAC3","MAC5","GA1-30","GA3-30","GA4-30","GA5-40","GA6-20","GA7-20","GA8-20","GA9-20","GA10-20","MA12-30","MA13-30","MA14-30","MA15-30","MA16-20","MA17-20","MA18-20","MA20-20")

#list=("GAC1" "GAC2" "GAC3" "GAC4" "GAC5" "MAC1" "MAC2" "MAC3" "MAC5" "GA1-30" "GA3-30" "GA4-30" "GA5-40" "GA6-20" "GA7-20" "GA8-20" "GA9-20" "GA10-20" "MA12-30" "MA13-30" "MA14-30" "MA15-30" "MA16-20" "MA17-20" "MA18-20" "MA20-20")

#for x in ${list[@]}
#do
#    echo $x
#    tblastn -db /scratch/trm76056/Geuk21/ORP_assemblies/$x.ORP.fasta -query lysin.fa -out $x.lysin.out -outfmt '6 sseqid evalue'
#done

makeblastdb -in AllReads_edited.ORP.fasta -parse_seqids -blastdb_version 5 -dbtype nucl

cd metabolic_genes

for x in *.fa
do
    tblastx -db /scratch/trm76056/Geuk21/AllReads_edited.ORP.fasta -query $x -out AllReads.$x.out -outfmt '6 qseqid sseqid sseq' -evalue 1e-5 -num_threads 4 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -max_target_seqs 1
done

cd $SLURM_SUBMIT_DIR

#list=("new_OG0016514.fa", "new_OG0016518.fa", "new_OG0016528.fa", "new_OG0016537.fa", "new_OG0016540.fa", "new_OG0016547.fa", "new_OG0016554.fa", "new_OG0016572.fa", "new_OG0016573.fa", "new_OG0016575.fa", "new_OG0016582.fa", "new_OG0016585.fa", "new_OG0016588.fa", "new_OG0016603.fa", "new_OG0016618.fa", "new_OG0016619.fa", "new_OG0016630.fa", "new_OG0016638.fa", "new_OG0016668.fa", "new_OG0016699.fa", "new_OG0016708.fa", "new_OG0016718.fa", "new_OG0016722.fa", "new_OG0016741.fa", "new_OG0016745.fa", "new_OG0016753.fa", "new_OG0016758.fa", "new_OG0016759.fa", "new_OG0016761.fa", "new_OG0016772.fa", "new_OG0016773.fa", "new_OG0016788.fa", "new_OG0016789.fa", "new_OG0016801.fa", "new_OG0016822.fa", "new_OG0016825.fa", "new_OG0016827.fa", "new_OG0016831.fa", "new_OG0016844.fa", "new_OG0016845.fa", "new_OG0016850.fa", "new_OG0016851.fa", "new_OG0016910.fa", "new_OG0016922.fa", "new_OG0016934.fa", "new_OG0016939.fa", "new_OG0016956.fa", "new_OG0016972.fa", "new_OG0016981.fa", "new_OG0017000.fa", "new_OG0017010.fa", "new_OG0017014.fa", "new_OG0017027.fa", "new_OG0017028.fa", "new_OG0017062.fa", "new_OG0017122.fa", "new_OG0017138.fa", "new_OG0017172.fa", "new_OG0017175.fa", "new_OG0017218.fa", "new_OG0017220.fa", "new_OG0017255.fa", "new_OG0017257.fa", "new_OG0017262.fa", "new_OG0017264.fa", "new_OG0017281.fa", "new_OG0017298.fa", "new_OG0017320.fa", "new_OG0017322.fa", "new_OG0017326.fa", "new_OG0017346.fa", "new_OG0017356.fa", "new_OG0017364.fa", "new_OG0017370.fa", "new_OG0017373.fa", "new_OG0017376.fa", "new_OG0017377.fa", "new_OG0017380.fa", "new_OG0017382.fa", "new_OG0017394.fa", "new_OG0017404.fa", "new_OG0017406.fa", "new_OG0017420.fa", "new_OG0017426.fa", "new_OG0017449.fa", "new_OG0017477.fa", "new_OG0017480.fa", "new_OG0017484.fa", "new_OG0017526.fa", "new_OG0017533.fa", "new_OG0017557.fa", "new_OG0017573.fa", "new_OG0017590.fa", "new_OG0017640.fa", "new_OG0017652.fa", "new_OG0017658.fa", "new_OG0017660.fa", "new_OG0017665.fa", "new_OG0017668.fa", "new_OG0017682.fa", "new_OG0017704.fa", "new_OG0017715.fa", "new_OG0017723.fa", "new_OG0017724.fa", "new_OG0017749.fa", "new_OG0017817.fa", "new_OG0017823.fa", "new_OG0017826.fa", "new_OG0017827.fa", "new_OG0017884.fa", "new_OG0017890.fa", "new_OG0017913.fa", "new_OG0017928.fa")

#for str in ${list[@]}:
#do
#    str=${str%?}
#    echo $str
#    blastx -db refseq_protein -query coding_files/OrthoFinder/Results_Nov29/Single_Copy_Orthologue_Sequences/$str -out $str.out -outfmt '6 qseqid sseqid evalue' -num_threads 4
#done
