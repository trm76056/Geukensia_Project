#this script will get the nucleotide sequence from teh transcriptomes for the metabolic loci

file='creatine_kinase.all.out'

while read line
do
    stringarray=($line)
    samtools faidx ORP_assemblies/${stringarray[0]}.ORP.fasta ${stringarray[1]} 
done < $file > $file.fasta

