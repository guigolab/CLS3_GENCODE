## Objective: 
# Test the enrichment for GWAS signal in CLS detected rather than background, and compare with annotation v27. 

Requirements:

1. ENCODE Blacklist region
```
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz 
gunzip ENCFF356LFX.bed.gz 
```

2. hg38 chromosome sizes
```
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
sed 's/chr[0-9]*_//g' hg38.chrom.sizes | sed 's/chr[a-Z]*_//g' | sed 's/_alt//g' | sed 's/_random//g' | sed 's/v/./g' > hg38.chrom.sizes.adapted 
awk -v OFS="\t" '{ print $1, 0, $2}' hg38.chrom.sizes.adapted > hg38.chrom.sizes.bed
```

3. hg38 centromeres
```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz
zcat centromeres.txt.gz | cut -f2,3,4 > centromeres.bed
```

4. hg38 gaps
```
http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=2307526118_6zkBZXXowMh9lw5r9xhw4vnmKW[%E2%80%A6]871&hgta_outputType=bed&hgta_outFileName=hg38.gaps.bed 
>> Get output: hg38.gaps.original.bed << 
sed 's/chr[0-9]*_//g' hg38.gaps.original.bed | sed 's/chr[a-Z]*_//g' | sed 's/_alt//g' | sed 's/_random//g' | sed 's/v/./g' > hg38.gaps.bed 
```

5. Gencode v27 primary assembly and lncRNA annotation separatedly.
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.primary_assembly.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.long_noncoding_RNAs.gtf.gz
```

6. GWAS Catalogue
```
wget --output-document gwas_catalog_v1.0.2-associations_e100_r2023-10-11.tsv https://www.ebi.ac.uk/gwas/api/search/downloads/alternative
wget http://www.ebi.ac.uk/efo/efo.obo
wget https://www.ebi.ac.uk/gwas/api/search/downloads/trait_mappings

./processing_catalog.sh

#Associate each SNP to its proxy in 5kb window.
format_snps_ld.R -i gwas_coordinates_resolvedName_catalog.bed 
awk -F"\t" -v OFS="\t" '{ print $1, $2-1, $2, $3, $4, $5}' gwas_coordinates_resolvedName_catalog.noLD.5k.bed > gwas_pruned_coordinates.bed
```

7. CLS - GENCODEv47 mapping
```
wget https://public-docs.crg.es/rguigo/Data/gkaur/CLS3_finalFiles/v47-CLS3mapping_status.txt
```

8. Extended CLS IDs
```
wget https://public-docs.crg.es/rguigo/Data/gkaur/CLS3_finalFiles/gencodev47Files/gencode.v47.primary_assembly.annotation.enhanced.gtf
grep "CLS3i" /users/rg/tperteghella/Gencode/InputData/Final/enhanced_annotation_v47.refined.gtf | grep -w exon | cut -d"\"" -f10 | tr "," "\n" | sort -u > ic_extending47.ids 
grep "CLS3i" /users/rg/tperteghella/Gencode/InputData/Final/enhanced_annotation_v47.refined.gtf | grep -w exon | cut -d"\"" -f2,10 | tr "\"" "\t" | awk -v OFS="\t" '{ split($2, ics, ","); for (i in ics) { print $1, ics[i] } }' | sort -u > mapping_ic_ID.tsv
```

9. CLS Data
```
cat Hv3_splicedmasterTable_refined.gtf Hv3_unsplicedmasterTable_refined.gtf > Hv3_intronchains_refined.gtf
awk '$9 ~ /Intergenic/' /users/rg/tperteghella/Gencode/InputData/Final/v47-CLS3mapping_status.txt | cut -f4 | tr "," "\n" | sort -u | grep -v OLD | grep -v UNMAPPED > ic_intergenic_in47.ids 
grep -Ff ic_intergenic_in47.ids <(cut -f1,4 ~/Gencode/InputData/Final/v47-CLS3mapping_status.txt | awk -v OFS="\t" '{ split($2, ics, ","); for (i in ics) { print $1, ics[i] } }' ) >> mapping_ic_ID.tsv 

grep -Ff <(cat ic_intergenic_in47.ids ic_extending47.ids) Hv3_intronchains_refined.gtf | awk ' $0 ~ /chr/ ' > CLS3i.intergenicICv27.annotation.gtf

grep -w exon CLS3i.intergenicICv27.annotation.gtf | awk -F"\t" -v OFS="\t" '{split($9, tags, "\""); print $1,$4,$5,tags[2] }' | sort -k1,1 -k2,2n > CLS3i.intergenic.loci.v27.exon.bed 
grep -v exon CLS3i.intergenicICv27.annotation.gtf | awk -F"\t" -v OFS="\t" '{split($9, tags, "\""); print $1,$4,$5,tags[2] }' | sort -k1,1 -k2,2n > CLS3i.intergenic.loci.v27.transcript.bed 

bedtools merge -i CLS3i.intergenic.loci.v27.exon.bed -c 4 -o distinct > CLS3i.intergenic.loci.v27.exon.merged.bed 
bedtools merge -i CLS3i.intergenic.loci.v27.transcript.bed -c 4 -o distinct > CLS3i.intergenic.loci.v27.transcript.merged.bed

awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2+1 }END{print SUM}' CLS3i.intergenic.loci.v27.transcript.merged.bed  #260.006.639 bp 
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2+1 }END{print SUM}' CLS3i.intergenic.loci.v27.exon.merged.bed #22.816.468 bp 
#+1 because were taken from GTF without transforming in 0-based as BED requires. 
```
   
10. Formatted references: GENCODE v27 annotation, Intergenic space, and Decoy Models
```
./create_references.sh
```

11. Bedtools 2.29.2
12. bedGraphToBigWig
13. Deeptools2 3.5.2

# Intersect references with the GWAS catalog
```
# GENCODE v27 annotation and decoy models
for i in annotation_*merged.bed; do name=$(basename $i); name=${name/merged./}; bedtools intersect -a gwas_pruned_coordinates.bed -b $i -wo -f 1.0 > ${name/annotation/bedtools_gwas}; done

#CLS3 models
for i in CLS3i*.merged.bed; do name=$(basename $i); name=${name/merged./}; bedtools intersect -a gwas_pruned_coordinates.bed -b $i -wo -f 1.0 > ${name/CLS3i/bedtools_gwas}; done

# Intergenic space before and after enxtending with CLS
for i in *nocomplexregions.bed; do name=$(basename $i); name=${name/_nocomplexregions/}; bedtools intersect -a gwas_pruned_coordinates.bed -b $i -wo -f 1.0 > bedtools_gwas.$name; done
```
Statistics are collected at: [GWAS_Stats]()

# Compare the density of GWAS across main GENCODE v27 annotation, CLS and decoy models.
```
gwas_pvaldensity.R -t bedtools_gwas_CLS3i.intergenic.loci.v27.transcript.bed -r bedtools_gwas_decoy.transcript.bed -c bedtools_gwas_lncRNA.bed -p bedtools_gwas_proteincoding.bed 
```

# Visualise density by Targeted Catalog
```
#Prepare annotation file:
grep -Ff <( cut -f4 CLS3i.intergenic.loci.v27.transcript.bed ) Hv3_intronchains_refined.gtf | grep -v exon | cut -d"\"" -f2,6,18 | tr "\"" "\t" | \
awk -F"\t" '{ split($2,db,","); if ( length(db) > 1 ) for(t in db) print $1"\t"db[t]"\t"$3; else print $0}'| \
awk -F"\t" '{ if ($3 ~ /SIDH....C/ && $3 ~ /SIDH....P/) print $1"\t"$2"\tCommon\t"$3; else if ($3 ~ /SIDH....C/) print $1"\t"$2"\tPostCapture\t"$3; else print $1"\t"$2"\tPreCapture\t"$3 }' | \
sort | uniq > CLS3i.intergenicICv27.annotation.tsv

plot_annotation_biotype.R -i bedtools_gwas_CLS3i.intergenic.loci.v27.exon.bed -a CLS3i.intergenicICv27.annotation.tsv -k covered_areas.tsv -w hits.tsv -e T -b CLS_exonic 
```

# Visualise density along gene body and surrounding
```
#Get gene bodies to test; groups

awk -F"\t" -v OFS="\t" ' $3 ~ /transcript/ {split($9, tags, "\""); print $1,$4-1,$5,$7,tags[4] }' random_replicates_locirelocation.gtf > decoy.genebody.bed
grep -v exon CLS3i.intergenicICv27.annotation.gtf | awk -F"\t" -v OFS="\t" '{split($9, tags, "\""); print $1,$4-1,$5,$7,tags[2] }' > CLS.genebody.bed 

zgrep -v exon gencode.v27.primary_assembly.annotation.gtf.gz | grep -v "^#" | awk -F"\t" -v OFS="\t" ' $3 ~ /gene/ {split( $9, tags, "\""); print $1, $4, $5, $7, tags[2]}' > gencodev27.genebody.bed 
grep -Ff <( cut -f4 ../bedfiles/annotation_proteincoding.bed ) gencodev27.genebody.bed > proteincoding.genebody.bed
grep -Ff <( cut -f4 ../bedfiles/annotation_lncRNA.bed ) gencodev27.genebody.bed > lncrna.genebody.bed

# Get the bigwig of the GWAS catalog 
awk -F"\t" -v OFS="\t" '{ print $1,$6-1,$6,1 }' gwas_pruned_coordinates.bed | sort -k1,1 -k2,2n -u > gwas_prunedcoordinates.bedGraph 
bedGraphToBigWig gwas_prunedcoordinates.bedGraph hg38.chrom.sizes gwas_prunedcoordinates.bw 

awk -F"\t" -v OFS="\t" '{ print $1,$2,$3,1 }' gwas_pruned_coordinates.bed | sort -k1,1 -k2,2n -u > gwas_coordinates.bedGraph  
bedGraphToBigWig gwas_coordinates.bedGraph hg38.chrom.sizes gwas_coordinates.bw 

#Compute Matrix
computeMatrix scale-regions -R  CLS.genebody.bed lncrna.genebody.bed proteincoding.genebody.bed decoy.genebody.bed -S gwas_coordinates.bw -o gwas.tsv.gz --upstream 20000 --downstream 20000 --sortRegions ascend --missingDataAsZero -p 8 --smartLabels --binSize 200 --regionBodyLength 5000 --averageTypeBins sum

#Compute matrix per catalog. Needed for clustering heatmap later.
computeMatrix scale-regions -R CLS.genebody.bed -S gwas_coordinates.bw -o CLS.tsv.gz --upstream 15000 --downstream 15000 --sortRegions ascend --missingDataAsZero -p 8 --smartLabels --binSize 500 --regionBodyLength 5000 --skipZeros
computeMatrix scale-regions -R lncrna.genebody.bed -S gwas_coordinates.bw -o lncrna.tsv.gz --upstream 15000 --downstream 15000 --sortRegions ascend --missingDataAsZero -p 8 --smartLabels --binSize 500 --regionBodyLength 5000 --skipZeros
computeMatrix scale-regions -R proteincoding.genebody.bed -S gwas_coordinates.bw -o proteincoding.tsv.gz --upstream 15000 --downstream 15000 --sortRegions ascend --missingDataAsZero -p 8 --smartLabels --binSize 500 --regionBodyLength 5000 --skipZeros
computeMatrix scale-regions -R decoy.genebody.bed -S gwas_coordinates.bw -o decoy.tsv.gz --upstream 15000 --downstream 15000 --sortRegions ascend --missingDataAsZero -p 8 --smartLabels --binSize 500 --regionBodyLength 5000 --skipZeros

#Plot
plotProfile -m gwas.tsv.gz -o gwas.profile.svg --perGroup --plotType heatmap --regionsLabel CLS lncRNA proteinCoding Decoy --yMin 0 --yMax 0.0001 --endLabel TTS --plotFileFormat svg
```
