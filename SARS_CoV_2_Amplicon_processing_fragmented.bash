#!/bin/bash

#Unique name given to each library. This will be provided in a list file for each library.  
file_base=$1

#Creates a log for each processing step for each library 
log=${file_base}.pipeline.log

{

echo "***********************************" 
echo "begin consensus generation for sample: $file_base" 
echo "***********************************" 

#Hard coded paths to the reference genome you are aligning to, the .bed file for primer postions, .fasta file containing all primer sequences, and .tsv file containg primer pair info 
#NOTE must provide this path for your own system
#Easy way to determine path is to enter directory where each file is located, type "pwd". Copy path below and add /your_file_name for each 
reference_genome=SARS-CoV-2_Reference_Genome.fa
primer_bed=nCoV_2019_PRIMER_bedfile.bed
#primer_fasta=/home/jrf69/project/SARS-CoV-2_Seq/References/nCoV_2019_V3_primers.fa
#pair_information=/home/jrf69/project/SARS-CoV-2_Seq/References/nCoV_2019_V3_pair_info.tsv

#Change into directory with symlinks to paired read data for specific sample
cd ${file_base}/Unaligned/ 

#variables set of each Read 1 and Read 2 files with same file base. Must be done for paired end sequencing
f1=${file_base}_*_R1_001.fastq.gz
f2=${file_base}_*_R2_001.fastq.gz

#Align paired end reads to the reference covid virus genome
bwa mem -t 16 $reference_genome $f1 $f2 | samtools view -b -F 4 -F 2048 | samtools sort -o ${file_base}_aln.bam

#Perform QC and soft clip primer sequences 
ivar trim -e -i ${file_base}_aln.bam -b $primer_bed -p ${file_base}_aln_trimmed.bam 

#sort output of ivar 
samtools sort ${file_base}_aln_trimmed.bam -o ${file_base}_aln_trimmed_sorted.bam 

#index output of ivar 
samtools index ${file_base}_aln_trimmed_sorted.bam  

#Use samtools to identify variants in reads compared to references, call consensus sequence with ivar 
samtools mpileup -A -d 1000 -Q 0 ${file_base}_aln_trimmed_sorted.bam | ivar consensus -t 0.60 -m 20 -p ${file_base}_consensus

#Use samtools to identify variants in reads compared to reference, output variants in tsv
samtools mpileup -A -d 0 -Q 0 ${file_base}_aln_trimmed_sorted.bam | ivar variants -p ${file_base}_variants -q 20 -t 0.02 -r $reference_genome 

awk -F '\t' '(0.2 < $11 && $11 < 0.80) && ($14=="TRUE")' ${file_base}_variants.tsv > ${file_base}_frequency.tsv

wc -l ${file_base}_frequency.tsv > ${file_base}_frequency_count.txt

#remove intermediate files not neccesary for analysis 
rm ${file_base}_aln_trimmed.bam
rm ${file_base}_aln_trimmed_sorted.bam.bai
rm ${file_base}_aln.bam.bai

echo "***********************************" 
echo "DONE WITH CONSENSUS GENERATION"
echo "***********************************" 

#finish log file for pipeline 
} 2>&1  | tee -a $log

