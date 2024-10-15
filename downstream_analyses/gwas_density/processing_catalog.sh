## Processing the GWAS catalog

# Fix entries with special characters;
#1. (1956 ";" separated SNP -> identifying haplotypes) 
awk -F"\t" ' { for (i = 1; i <= NF; ++i) {if($i ~ ";") print i}} ' <(cat gwas_catalog_v1.0.2-associations_e100_r2023-10-11.tsv | awk -F"\t" '$22~";"') | sort | uniq -c 

#2. Fields with comma separated values: 12, 13, 15, 21, 22, 24, 25, 30 
cat gwas_catalog_v1.0.2-associations_e100_r2023-10-11.tsv | awk 'BEGIN{OFS=FS="\t"; colsI="12,13,15,21,22,25"; split(colsI, cols, ",")} {n=split($22, a, ";"); for(col in cols){d[cols[col]]=$(cols[col])}; for(i=1;i<=n;i++){for(col in cols){split(d[cols[col]], b, ";"); gsub(/ /,"", b[i]); $(cols[col])=b[i]} print }}' > gwas.catalog.extHap.tsv 

#3. Replace NAs
sed 's/\t\t/\tNA\t/g' gwas.catalog.extHap.tsv | sed 's/\t\t/\tNA\t/g' > gwas.catalog.extHap.NAs.tsv 

#4. Remove interaction terms
cat gwas.catalog.extHap.NAs.tsv | awk -F"\t" '$22!~/rs.+xrs.+/' > gwas.catalog.extHap.NAs.noInt.tsv

#5. Format 
cat  gwas.catalog.extHap.NAs.noInt.tsv | awk 'BEGIN{OFS=FS="\t"} NR>1 && $13!="NA" {if($23==1 && $24!="NA"){print "chr"$12, $13-1, $13, "rs"$24, "0", "+", $8, $(NF-2)} else {print "chr"$12, $13-1, $13, $22, "0", "+", $8, $(NF-2)}}' | sort -Vk1,1 | uniq > gwas.catalog.efo.bed7 

#6. Assign EFO mapping. Needed to reduce number of terms during semantic clustering in next step. 
Rscript efo_mapping.R 

#7. Remove temrs in parenthesis and simplify the semantic
sed 's/ (.*)//g' gwas.catalog.mappedefo.bed7 | sed 's/\[.*\]//g' | tr -d '[:punct:]' | sort -u > gwas.catalog.strip.bed7 

Rscript format_gwas_catalog.R
