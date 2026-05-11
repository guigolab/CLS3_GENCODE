## Download data on constrained regions for 29-way alignment and Merge together regions that are separated by just one base pair
```
wget http://www.broadinstitute.org/ftp/pub/assemblies/mammals/29mammals/29way_omega_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.gz
zcat 29way_omega_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.gz | sort -k1,1 -k2,2n | bedtools merge -i - -d 1 | gzip > 29way_omega_lods_elements_12mers.chr_specific.fdr_0.1.bed.gz
```
  
## Download data on constrained regions for 241-way alignment and merge together regions that are separated by just one base pair with an FDR < 5% (equivalent to a phyloP score > 2.27)
```
wget https://hgdownload.gi.ucsc.edu/goldenPath/hg38/cactus241way/cactus241way.phyloP.bw 
bigWigToWig cactus241way.phyloP.bw cactus241way.phyloP.wig
wig2bed --zero-indexed --do-not-sort < cactus241way.phyloP.wig | gzip > cactus241way.phyloP.bed.gz
zcat cactus241way.phyloP.bed.gz | awk -F'\t' '$5 > 2.27 {print $0}' | sort -k1,1 -k2,2n | bedtools merge -i - -d 1 | gzip > cactus241way.phyloP_fdr_0.05.bed.gz
```
