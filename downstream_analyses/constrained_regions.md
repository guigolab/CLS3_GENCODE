# The constrained genome as a proxy for functionality
Constrained regions in the genome are those with a significant lack of genetic variation, something that can be linked to their functional relevance. We aimed to explore how genes derived from the CLS data are associated with this highly-conserved regions, as a potential explanation their functionality. We have leveraged two sets of constrained regions obtained from the alignment of 29[^1] and 240[^2] mammals, respectively.

## 29-way alignment
We downloaded the set of reported constrained regions, computed at 12-bp resolution. Regions separated by only one base pair were merged together.
```
wget http://www.broadinstitute.org/ftp/pub/assemblies/mammals/29mammals/29way_omega_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.gz
zcat 29way_omega_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.gz | sort -k1,1 -k2,2n | bedtools merge -i - -d 1 | gzip > 29way_omega_lods_elements_12mers.chr_specific.fdr_0.1.bed.gz
```
## 240-way alignment
We downloaded the single-base phyloP p-values associated with this alignment and subsetted those regions with a score >2.27 (equivalent to a threshold of >5% FDR). Regions separated by only one base pair were merged together.
```
wget https://hgdownload.gi.ucsc.edu/goldenPath/hg38/cactus241way/cactus241way.phyloP.bw
bigWigToWig cactus241way.phyloP.bw cactus241way.phyloP.wig
wig2bed --zero-indexed --do-not-sort < cactus241way.phyloP.wig | gzip > cactus241way.phyloP.bed.gz
zcat cactus241way.phyloP.bed.gz | awk -F'\t' '$5 > 2.27 {print $0}' | sort -k1,1 -k2,2n | bedtools merge -i - -d 1 | gzip > cactus241way.phyloP_fdr_0.05.bed.gz
```
[^1]: Lindblad-Toh, K., Garber, M., Zuk, O. _et al_. A high-resolution map of human evolutionary constraint using 29 mammals. _Nature_ **478**, 476–482 (2011). https://doi.org/10.1038/nature10530
[^2]: Matthew J. Christmas _et al_. Evolutionary constraint and innovation across hundreds of placental mammals. _Science_ **380**, eabn3943 (2023). DOI:10.1126/science.abn3943 https://doi.org/10.1126/science.abn3943
