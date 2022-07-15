#!/usr/bin/env bash

file=variant_frequencies.txt; 

echo -e 'POS\tREF\tALT\tGene\tVariantType\tFeatureType\tIsPseudo\tIsGenic\tIsTransition\tIsSynonymous\tAminoAcidChange\tSNPCodonPosition\tCodonPosition\tAltAminoAcid\tRefAminoAcid\tAltCodon\tRefCodon\tVAF\n' > ${file}; 

bcftools merge `ls */Unaligned/*_annot.vcf.gz | grep -v tmp` | cut -f2,4,5,6,10- | sed "s:QUAL:GENE:g" | grep -v "##" | grep "ALT" | sed "s:\#::g" | cut -f5- | tr "/\t" "\n\n" | grep -v -E "(Unaligned|marked_sort)" | tr "\n" "\t" >> ${file}; cat ${file} | grep -v ^$ | tr "\n" "\t" > tmp.010101; echo -e "\n" >> tmp.010101; mv tmp.010101 ${file}; bcftools merge `ls */Unaligned/*_annot.vcf.gz | grep -v tmp` | bcftools norm -f ancestral_pers.fasta -m -any | bcftools filter -i 'FORMAT/VAF>0.01' | bcftools query -f '%POS\t%REF\t%ALT\t%Gene\t%VariantType\t%FeatureType\t%IsPseudo\t%IsGenic\t%IsTransition\t%IsSynonymous\t%AminoAcidChange\t%SNPCodonPosition\t%CodonPosition\t%AltAminoAcid\t%RefAminoAcid\t%AltCodon\t%RefCodon[\t%VAF]\n' | perl -p -e "s/\t\./\t0/ig" >> ${file}; 

cat ${file} | grep -v ^$ | grep -v -E "(Insertion|Deletion)" | sed "s:VAF\t::g" > tmp.010101; 

mv tmp.010101 ${file}
