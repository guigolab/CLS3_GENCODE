## Objective: 
# Generate a random set of transcripts by relocating in the intergenic space (according to GENCODEv27 extended with CLS models) the loci obtained through CLS. 

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

5. Gencode v27 primary assembly
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.primary_assembly.annotation.gtf.gz
awk -F"\t" -v OFS="\t" '{ split($9, tags, "\""); print $1,$4,$5,tags[2],tags[4]}' <(zcat gencode.v27.primary_assembly.annotation.gtf.gz) | sort -k1,1 -k2,2n > gencode_v27.sorted.gtf 
```

6. CLS models
```
zgrep -Ff <(zgrep "spliced \"1\"" Hv3+withinTmerge_gencodev27.loci.gtf.gz | cut -d"\"" -f2 | sort -u) Hv3+withinTmerge_gencodev27.loci.gtf.gz > spliced_cls.loci.gtf
awk -F"\t" -v OFS="\t" '{ split($9, tags, "\""); print $1,$4,$5,tags[2],tags[4]}' spliced_cls.loci.gtf | sort -k1,1 -k2,2n > cls.sorted.gtf
```

7. Bedtools 2.29.2



# Data preprocessing

Merge the CLS loci with the annotation v27 and extend boundaries 10kb
```
cat <( cut -f1,2,3,4 gencode_v27.sorted.gtf) <(cut -f1,2,3,4 cls.sorted.gtf ) | sort -k1,1 -k2,2n > gencode_v27+cls.sorted.bed  
bedtools merge -i gencode_v27+cls.sorted.bed -c 4 -o distinct > gencode_v27+cls.merged.bed
bedtools slop -i gencode_v27+cls.merged.bed -g hg38.chrom.sizes.adapted -b 10000 > gencode_v27+cls.merged.extended10kb.bed
bedtools merge -i gencode_v27+cls.merged.extended10kb.bed -c 4 -o distinct > gencode_v27+cls.merged.final.bed  
``` 

# Obtain the intergenic space
```
bedtools subtract -a hg38.chrom.sizes.bed -b gencode_v27+cls.merged.final.bed | sort -k1,1 -k2,2n > intergenic.gencode_v27+cls.bed
```

Exclude the complex regions; ENCODE blacklist regions, centromeres and gaps.
```
cat ENCFF356LFX.bed centromeres.bed hg38.gaps.bed | sort -k1,1 -k2,2n > complex_regions.bed 
bedtools merge -i complex_regions.bed > complex_regions.merged.bed
bedtools subtract -a intergenic.gencode_v27+cls.bed -b complex_regions.merged.bed | sort -k1,1 -k2,2n > intergenic.gencode_v27+cls.no_complex_regions.bed

#Exclude regions not of interest for relocation (alts, patches, and chrM)
grep "chr" intergenic.gencode_v27+cls.no_complex_regions.bed | grep -v "chrM" > intergenic.gencode_v27+cls.no_complex_regions.canonical_assembly.bed

#Check size of intergenic space
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' intergenic.gencode_v27+cls.no_complex_regions.canonical_assembly.bed
#354.642.528bp 
```
 
# Relocate all loci at random  
```
python3 relocate_loci.py intergenic.gencode_v27+cls.no_complex_regions.canonical_assembly.bed spliced_cls.loci.gtf
```

# Output
Find the set of randomly relocated, strand-less transcripts [here](https://zenodo.org/api/records/13946596/draft/files/random_replicates_locirelocation.gtf.gz/content).
