#!/usr/bin/env bash

i=$1 #file base
reference_genome_gff=$2 #Reference genome

picard MarkDuplicates --REMOVE_DUPLICATES true --INPUT ${i}/Unaligned/${i}_aln_trimmed_sorted.bam --OUTPUT ${i}/Unaligned/${i}.marked.bam --VALIDATION_STRINGENCY LENIENT --METRICS_FILE ${i}/Unaligned/${i}.metrics.txt > ${i}/Unaligned/${i}.MarkDuplicates.log; 

samtools sort -o ${i}/Unaligned/${i}.marked_sort.bam ${i}/Unaligned/${i}.marked.bam; samtools index ${i}/Unaligned/${i}.marked_sort.bam; 

echo -e ${i}"\t"`samtools view -F 0x40 ${i}/Unaligned/${i}.marked_sort.bam | cut -f 1 | sort | uniq | wc -l` > ${i}/Unaligned/${i}.reads.txt; 

bedtools genomecov -bga -ibam ${i}/Unaligned/${i}.marked_sort.bam > ${i}/Unaligned/${i}.depth.txt; 
bedtools genomecov -bga -ibam ${i}/Unaligned/${i}.marked_sort.bam | awk '$4 < 10' > ${i}/Unaligned/low_coverage_sites.bed

rm -rf ${i}/Unaligned/${i}.marked.bam; 

module load BCFtools/1.12-GCCcore-10.2.0 

bcftools mpileup --max-depth 1000000 --max-idepth 1000000 --per-sample-mF --count-orphans --no-BAQ -Ov --min-BQ 0 --annotate AD,ADF,ADR,DP,SP,AD,ADF,ADR -f ancestral_pers.fasta ${i}/Unaligned/${i}.marked_sort.bam | bcftools +fill-tags - | bcftools call -mv -Ou --prior-freqs AN,AC -A --variants-only --keep-alts --keep-masked-ref -o ${i}/Unaligned/${i}.tmp.vcf; bcftools +fill-tags -o ${i}/Unaligned/${i}.tmp1.vcf ${i}/Unaligned/${i}.tmp.vcf -- -t FORMAT/VAF; bcftools filter -i 'FORMAT/DP>=50' ${i}/Unaligned/${i}.tmp1.vcf > ${i}/Unaligned/${i}.tmp5.vcf;  bcftools norm --multiallelics -any -Ov -o ${i}/Unaligned/${i}.vcf --check-ref s --fasta-ref ancestral_pers.fasta ${i}/Unaligned/${i}.tmp5.vcf;  

rm -rf ${i}/Unaligned/${i}.marked_sort.bam; 

bcftools view -Oz -o ${i}/Unaligned/${i}.vcf.gz ${i}/Unaligned/${i}.vcf; bcftools index  ${i}/Unaligned/${i}.vcf.gz; bcftools stats ${i}/Unaligned/${i}.vcf.gz > ${i}/Unaligned/${i}.stats.ALL.txt; cat ${i}/Unaligned/${i}.stats.ALL.txt | grep -A 12 "\[3\]type" | perl -p -e "s/\#\W+//g" | perl -p -e "s/\[[0-9]+\]//g" > ${i}/Unaligned/${i}.stats.mut.ALL.txt; cat ${i}/Unaligned/${i}.stats.ALL.txt | grep -A 1 "\[5\]ts/tv" | perl -p -e "s/\#\W+//g" | perl -p -e "s/\[[0-9]+\]//g" | tr " " "_" | sed "s:(::g" | sed "s:)::g" | sed "s:\/::g" > ${i}/Unaligned/${i}.stats.tstv.ALL.txt; 

cat ancestral_pers.bed | while read j; do bcftools view -Oz -o ${i}/Unaligned/${i}.vcf.gz ${i}/Unaligned/${i}.vcf; bcftools index -f ${i}/Unaligned/${i}.vcf.gz; refname=`echo $j | awk '{print $1}'`; start=`echo $j | awk '{print $2}'`; end=`echo $j | awk '{print $3}'`; gene=`echo $j | awk '{print $4}'`; bcftools stats -r "${refname}:${start}-${end}"   ${i}/Unaligned/${i}.vcf.gz > ${i}/Unaligned/${i}.${gene}.gene.stats.txt;  cat ${i}/Unaligned/${i}.${gene}.gene.stats.txt | grep -A 12 "\[3\]type" | perl -p -e "s/\#\W+//g" | perl -p -e "s/\[[0-9]+\]//g" > ${i}/Unaligned/${i}.${gene}.gene.stats.mut.txt; cat ${i}/Unaligned/${i}.${gene}.gene.stats.txt | grep -A 1 "\[5\]ts/tv" | perl -p -e "s/\#\W+//g" | perl -p -e "s/\[[0-9]+\]//g" | tr " " "_" | sed "s:(::g" | sed "s:)::g" | sed "s:\/::g" > ${i}/Unaligned/${i}.${gene}.gene.stats.tstv.txt; done

source ~/project/miniconda3/bin/activate vcf
vcf-annotator --output ${i}/Unaligned/${i}_annot.vcf ${i}/Unaligned/${i}.vcf $reference_genome_gff
source ~/project/miniconda3/bin/activate base

bcftools view -Oz -o ${i}/Unaligned/${i}_annot.vcf.gz ${i}/Unaligned/${i}_annot.vcf; bcftools index -f ${i}/Unaligned/${i}_annot.vcf.gz; 

