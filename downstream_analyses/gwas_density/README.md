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

5. Gencode primary assembly and v27 lncRNA annotation.
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v47.primary_assembly.annotation.gtf.gz
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

7. [Disjoin Gene IDs](https://github.com/guigolab/CLS3_GENCODE/tree/main/complementary_data/gencode_biotypes_dataset).
8. Bedtools 2.29.2
9. bedGraphToBigWig
10. Deeptools2 3.5.2

# Intersect references with the GWAS catalog
```
for i in exon transcript; do 
    zgrep -wFf protein_coding.transcripts.v47.disjoint.ids gencode.v47.primary_assembly.annotation.gtf.gz | awk -F"\t" -v OFS="\t" -v b=$i '$3 == b { split($9, tags, "\""); print $1, $4-1, $5, tags[4]}' > protein_coding.transcripts.v47.disjoint.$i.bed 
   zgrep -wFf lncRNA.transcripts.v27.disjoint.ids gencode.v27.long_noncoding_RNAs.gtf.gz | awk -F"\t" -v OFS="\t" -v b=$i '$3 == b { split($9, tags, "\""); print $1, $4-1, $5, tags[4]}' > lncRNA.transcripts.v27.disjoint.$i.bed 
   zgrep -wFf cls.transcripts.disjoint.ids gencode.v47.primary_assembly.annotation.gtf.gz | awk -F"\t" -v OFS="\t" -v b=$i '$3 == b { split($9, tags, "\""); print $1, $4-1, $5, tags[4]}' > cls.transcripts.disjoint.$i.bed 
   grep -wFf decoys.transcripts.disjoint.ids random_replicates.gtf | awk -F"\t" -v OFS="\t" -v b=$i '$3 == b { split($9, tags, "\""); print $1, $4-1, $5, tags[4]}' > decoys.transcripts.disjoint.$i.bed
done
```

# Extending to proximal regions +2kb on both sides 
```
grep -v chrM protein_coding.transcripts.v47.disjoint.transcript.bed | bedtools slop -i - -g hg38.chrom.sizes -b 2000 > protein_coding.transcripts.v47.disjoint.transcript.proximal.bed 
bedtools slop -i lncRNA.transcripts.v27.disjoint.transcript.bed -g hg38.chrom.sizes -b 2000 > lncRNA.transcripts.v27.disjoint.transcript.proximal.bed 
bedtools slop -i decoys.transcripts.disjoint.transcript.bed -g hg38.chrom.sizes -b 2000 > decoys.transcripts.disjoint.transcript.proximal.bed 
bedtools slop -i cls.transcripts.disjoint.transcript.bed -g hg38.chrom.sizes -b 2000 > cls.transcripts.disjoint.transcript.proximal.bed 
```

# Computing the non-overlapping span and intersecting the catalog 
```
for i in *bed; do echo $i; sort --parallel=4 -k1,1 -k2,2n $i | bedtools merge -i - -c 4 -o distinct | awk '{ sum += $3-$2 } END { print sum }'; done 
for i in *.bed; do newname=$(basename $i); echo $newname; bedtools intersect -a $i -b gwas_pruned_coordinates.bed -wao > ${newname/disjoint/gwas}; done  
for i in *gwas.bed; do echo $i; awk '$11 != 0' $i | cut -f5,10 | sort --parallel=4 -u | wc -l; done 
```


# Visualise density along gene body and surrounding
```
zgrep -wFf protein_coding.loci.v47.disjoint.ids gencode.v47.primary_assembly.annotation.gtf.gz | awk -F"\t" -v OFS="\t" '$3 == "gene" { split($9, tags, "\""); print $1, $4-1, $5, tags[2], "0", $7 }' > protein_coding.v47.genebody.bed 
zgrep -wFflncRNA.loci.v27.disjoint.ids gencode.v27.long_noncoding_RNAs.gtf.gz | awk -F"\t" -v OFS="\t" '$3 == "gene" { split($9, tags, "\""); print $1, $4-1, $5, tags[2], "0", $7 }' > lncRNA.v27.genebody.bed 
zgrep -wFf cls.loci.disjoint.ids gencode.v47.primary_assembly.annotation.gtf.gz | awk -F"\t" -v OFS="\t" '$3 == "gene" { split($9, tags, "\""); print $1, $4-1, $5, tags[2], "0", $7 }' > cls.genebody.bed 
grep -wFf decoys.loci.disjoint.ids random_replicates.gtf | awk -F"\t" -v OFS="\t" '$3 == "transcript" { split($9, tags, "\""); print $1, $4-1, $5, tags[2], "0", $7 }' | sort --parallel=4 -k1,1 -k2,2n | bedtools merge -i - -c 5 -o distinct > decoys.genebody.bed 

# Get the bigwig of the GWAS catalog 
awk -F"\t" -v OFS="\t" '{ print $1,$6-1,$6,1 }' gwas_pruned_coordinates.bed | sort -k1,1 -k2,2n -u > gwas_prunedcoordinates.bedGraph 
bedGraphToBigWig gwas_prunedcoordinates.bedGraph hg38.chrom.sizes gwas_prunedcoordinates.bw 

#Compute Matrix
computeMatrix scale-regions -R  CLS.genebody.bed lncrna.genebody.bed proteincoding.genebody.bed decoy.genebody.bed -S gwas_prunedcoordinates.bw -o gwas_pruned.tsv.gz --upstream 20000 --downstream 20000 --sortRegions ascend --missingDataAsZero -p 8 --smartLabels --binSize 200 --regionBodyLength 5000

#Plot
plotProfile -m gwas_pruned.tsv.gz -o gwas.profile.svg --perGroup --plotType heatmap --regionsLabel "CLS transcripts" "lncRNAs (v27)" "protein-coding (v27)" "decoy models" --yMin 0 --yMax 0.0001 --endLabel TTS --plotFileFormat svg
```
