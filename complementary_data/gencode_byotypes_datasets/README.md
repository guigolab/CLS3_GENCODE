# Datasets used for analysis
Accurate gene annotations are fundamental to functionally interpreting the activity (transcription, chromatin modifications, folding, protein binding) and the sequence variation of genomes. We wanted to show how the incorporation of the CLS data into GENCODE greatly enhances the functional interpretability of the human genome. 

## 1. Novel CLS models
In general, we have considered the novel CLS models with respect to GENCODE v27, although some analyses have been restricted to the novel human lncRNA loci in v47, complemented with intergenic spliced CLS models. 

Original files available [here](https://github.com/guigolab/CLS3_GENCODE/blob/main/data_release/README.md).

### Novel CLS models
```
egrep "refCompare \"Extends|refCompare \"Intronic|refCompare \"revIntronic|refCompare \"Antisense|refCompare \"Intergenic" Hv3_unsplicedmasterTable_refined.gtf Hv3_splicedmasterTable_refined.gtf | cut -d"\"" -f2 > novel_cls_ids.txt
```
The creation of TSS is documented [here](https://github.com/guigolab/CLS3_GENCODE/tree/main/complementary_data/tss_representatives).

### novel human lncRNA loci now in v47
```
pending
```

### intergenic spliced CLS models
We decided to complement the analyses with all those models not included in v47 solely because they did not reach the minimal recount threshold (9,772 genes, 22,211 transcripts)
```
pending
```

## 2. GENCODE v27 reference
GENCODE [v27](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.primary_assembly.annotation.gtf.gz) was used as reference given that was the geneset available at the time of the design. We have compared them to previously annotated lncRNAs (8,922 genes, 15,922 transcripts) and protein-coding genes (19,823 genes, 146,877 transcripts), so defined. 

### lncRNAs
In order to reduce the signal coming from lncRNAs overlapping protein-coding genes, we decided to discard all those lncRNAs which overlapped > 10% of their lenght with an annotated protein-coding gene.

```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.long_noncoding_RNAs.gtf.gz
awk ' $3 ~ /gene/' <(zcat gencode.v27.long_noncoding_RNAs.gtf.gz) | cut -d";" -f1 | awk -F"\t" -v OFS="\t" '$1 ~ /chr/ { print $1,$4-1,$5,$9}' | grep -v chrM | sort --parallel=4 -k1,1 -k2n,2 | uniq > annotation_lncRNA.bed 
grep -Ff <(bedtools intersect -a annotation_lncRNA.bed -b annotation_proteincoding.bed -wao -f 0.1 | awk -F"\t" '$9 == 0' | cut -f4) annotation_lncRNA.bed > annotation_lncRNA.nooverlap.bed && mv annotation_lncRNA.nooverlap.bed annotation_lncRNA.bed
awk ' $3 ~ /exon/' <(zcat gencode.v27.long_noncoding_RNAs.gtf.gz) | cut -d";" -f1 | awk -F"\t" -v OFS="\t" '$1 ~ /chr/ { print $1,$4-1,$5,$9}' | grep -Ff <( cut -f4 annotation_lncRNA.bed) | grep -v chrM | sort --parallel=4 -k1,1 -k2n,2 | uniq > annotation_lncRNA_exon.bed 
```

### protein-coding genes
```
awk ' $3 ~ /gene/' <(zcat gencode.v27.primary_assembly.annotation.gtf.gz) | awk -F";" '$2 ~ /protein_coding/' | cut -d";" -f1 | awk -F"\t" -v OFS="\t" '$1 ~ /chr/ { print $1,$4-1,$5,$9}' | grep -v chrM | sort --parallel=4 -k1,1 -k2n,2 | uniq > annotation_proteincoding.bed 
awk ' $3 ~ /exon/' <(zcat gencode.v27.primary_assembly.annotation.gtf.gz) | awk -F";" '$3 ~ /protein_coding/' | cut -d";" -f1 | awk -F"\t" -v OFS="\t" '$1 ~ /chr/ { print $1,$4-1,$5,$9}' | grep -v chrM | sort --parallel=4 -k1,1 -k2n,2 | uniq > annotation_proteincoding_exon.bed    
```

## 3. Background
For some analyses, we have also employed a set of decoy models (17,223 genes, 85,283 transcripts) that attempt to mimic the background (non-genic) behavior of the genome.
Find them [here](https://github.com/guigolab/CLS3_GENCODE/tree/main/complementary_data/decoy_models).



