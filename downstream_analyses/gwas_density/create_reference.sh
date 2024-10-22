#All annotations 
zgrep -v "^#" gencode.v27.primary_assembly.annotation.gtf.gz | cut -d";" -f1 | awk -F"\t" -v OFS="\t" '{ print $1,$4-1,$5,$9}' | sort --parallel=4 -k1,1 -k2n,2 | uniq > annotation_for_intergenic_bases.bed   
cat <(zgrep -v "^#" gencode.v27.primary_assembly.annotation.gtf.gz) CLS3i.intergenic.annotation.v27.gtf | cut -d";" -f1 | awk -F"\t" -v OFS="\t" '{ print $1,$4-1,$5,$9}' | sort --parallel=4 -k1,1 -k2n,2 | uniq > enhanced_for_intergenic_bases.bed 

#protein_coding 
awk ' $3 ~ /gene/' <(zcat gencode.v27.primary_assembly.annotation.gtf.gz) | awk -F";" '$2 ~ /protein_coding/' | cut -d";" -f1 | awk -F"\t" -v OFS="\t" '$1 ~ /chr/ { print $1,$4-1,$5,$9}' | grep -v chrM | sort --parallel=4 -k1,1 -k2n,2 | uniq > annotation_proteincoding.bed 
awk ' $3 ~ /exon/' <(zcat gencode.v27.primary_assembly.annotation.gtf.gz) | awk -F";" '$3 ~ /protein_coding/' | cut -d";" -f1 | awk -F"\t" -v OFS="\t" '$1 ~ /chr/ { print $1,$4-1,$5,$9}' | grep -v chrM | sort --parallel=4 -k1,1 -k2n,2 | uniq > annotation_proteincoding_exon.bed    

#lncRNA 
awk ' $3 ~ /gene/' <(zcat gencode.v27.long_noncoding_RNAs.gtf.gz) | cut -d";" -f1 | awk -F"\t" -v OFS="\t" '$1 ~ /chr/ { print $1,$4-1,$5,$9}' | grep -v chrM | sort --parallel=4 -k1,1 -k2n,2 | uniq > annotation_lncRNA.bed 
grep -Ff <(bedtools intersect -a annotation_lncRNA.bed -b annotation_proteincoding.bed -wao -f 0.1 | awk -F"\t" '$9 == 0' | cut -f4) annotation_lncRNA.bed > annotation_lncRNA.nooverlap.bed && mv annotation_lncRNA.nooverlap.bed annotation_lncRNA.bed
awk ' $3 ~ /exon/' <(zcat gencode.v27.long_noncoding_RNAs.gtf.gz) | cut -d";" -f1 | awk -F"\t" -v OFS="\t" '$1 ~ /chr/ { print $1,$4-1,$5,$9}' | grep -Ff <( cut -f4 annotation_lncRNA.bed) | grep -v chrM | sort --parallel=4 -k1,1 -k2n,2 | uniq > annotation_lncRNA_exon.bed 

#Decoy models
awk ' $3 ~ /exon/' ../decoy_models/random_replicates_locirelocation.gtf | cut -d";" -f1 | awk -F"\t" -v OFS="\t" '{ print $1,$4-1,$5,$9 }' | sort --parallel=4 -k1,1 -k2,2n > annotation_decoy.v27.exon.bed 
awk ' ! $3 ~ /exon/' ../decoy_models/random_replicates_locirelocation.gtf | cut -d";" -f1 | awk -F"\t" -v OFS="\t" '{ print $1,$4-1,$5,$9 }' | sort --parallel=4 -k1,1 -k2,2n > annotation_decoy.v27.transcript.bed 
 
for i in *bed; do sed -i 's/gene_id //g' $i; done 
for i in *bed; do sed -i 's/"//g' $i; done 
for i in *bed; do bedtools merge -i $i -c 4 -o distinct > ${i/.bed/.merged.bed}; rm $i; done 

#Intergenic space 
grep -Ff <( cut -f1 enhanced_for_intergenic_bases.bed | sort -u) hg38.chrom.sizes.bed > assembly.chr.txt 
grep -Ff <( cut -f1 enhanced_for_intergenic_bases.bed | sort -u) complex_regions.merged.bed > complex_regions.merged.chr.bed

bedtools subtract -a assembly.chr.txt -b annotation_for_intergenic_bases.merged.bed > intergenic.gencode_v27.bed 
bedtools subtract -a intergenic.gencode_v27.bed -b complex_regions.merged.chr.bed > intergenic.gencode_v27_nocomplexregions.bed 

bedtools subtract -a assembly.chr.txt -b enhanced_for_intergenic_bases.merged.bed > intergenic.enhanced_v27.bed 
bedtools subtract -a intergenic.enhanced_v27.bed -b complex_regions.merged.chr.bed > intergenic.enhanced_v27_nocomplexregions.bed 

#Extract length
for i in annotation*merged.bed; do echo $i; awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $i; done   
#Total length of intergenic region 
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' intergenic.enhanced_v27_nocomplexregions.bed
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' intergenic.gencode_v27_nocomplexregions.bed
#Put results in file: covered_area.tsv

#Extract Number of hits
#Whole catalog
for i in bedtools_gwas_*.bed; do IFS="_" read -r x y biotype <<< "${i/.bed/}"; echo -e "$biotype\t$(awk ' $11 != 0 ' $i | cut -f5,7 | sort -u | wc -l)"; done 
#Pruned catalog 
for i in bedtools_gwas_*.bed; do IFS="_" read -r x y biotype <<< "${i/.bed/}"; echo -e "$biotype\t$(awk ' $11 != 0 ' $i | cut -f5,10 | sort -u | wc -l)"; done
#Put results in file: hits.tsv

 

 
