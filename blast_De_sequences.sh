#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=blast
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=100:00:00
#SBATCH --mem=5G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=trm76056@uga.edu

cd $SLURM_SUBMIT_DIR

ml SAMtools/1.14-GCC-8.3.0
ml BLAST+/2.12.0-gompi-2020b

list=('NODE_79258_length_1109_cov_4931.262089_g45399_i0' 'NODE_146163_length_638_cov_1890.644940_g94544_i0' 'NODE_132818_length_709_cov_2763.802752_g31329_i2' 'TRINITY_DN8214_c0_g1_i8' 'NODE_107352_length_810_cov_2761.454422_g1131_i2' 'NODE_73781_length_1246_cov_532.526448_g39157_i2' 'TRINITY_DN1444_c0_g2_i3' 'R30245533' 'NODE_125206_length_755_cov_8521.335714_g76585_i0' 'TRINITY_DN14868_c0_g1_i5' 'TRINITY_DN4186_c0_g1_i3' 'NODE_63356_length_1407_cov_415.693047_g33282_i0' 'NODE_353088_length_278_cov_0.748768_g307517_i0' 'NODE_196290_length_456_cov_3684.783042_g142171_i0' 'NODE_457423_length_217_cov_1920.962963_g403074_i0' 'TRINITY_DN35973_c0_g1_i1' 'NODE_39014_length_1943_cov_5443.280191_g17692_i3' 'TRINITY_DN1168_c0_g1_i6' 'NODE_68470_length_1323_cov_3323.331230_g36309_i0' 'NODE_98459_length_964_cov_72595.322332_g56003_i0' 'TRINITY_DN2735_c0_g1_i1' 'NODE_32126_length_2164_cov_17374.526316_g15578_i2' 'NODE_114830_length_827_cov_2292.042746_g68229_i0' 'R30275356' 'NODE_231317_length_383_cov_2554.801829_g176997_i0' 'NODE_54182_length_1530_cov_356.033677_g29542_i0' 'NODE_72005_length_1208_cov_12316.528685_g38054_i2')

for item in ${list[@]}
do
    echo $item
    samtools faidx AllReads.ORP.fasta $item > $item.samout
    echo "blast "$item
    blastx -db nr -query $item.samout -out $item.blastout -outfmt '6 qseqid sseqid evalue' -remote
done
