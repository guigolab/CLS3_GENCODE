# Datasets used for analysis
Accurate gene annotations are fundamental to functionally interpreting the activity (transcription, chromatin modifications, folding, protein binding) and the sequence variation of genomes. We wanted to show how the incorporation of the CLS data into GENCODE greatly enhances the functional interpretability of the human genome. 

## 1. Novel CLS models
In general, we have considered the novel CLS models with respect to GENCODE v27, although some analyses have been restricted to the novel human lncRNA CLS loci in v47. 
Original files available [here](https://github.com/guigolab/CLS3_GENCODE/blob/main/data_release/README.md).

### CLS transcripts: transcripts - genes

The lncRNA transcripts and genes now in v47 because CLS data, as mapped [here](https://zenodo.org/records/15004659/files/v47-CLS3mapping_status.txt?download=1), have been assigned a novelty category as detailed in [GENCODE-CLS3 Mappings](https://github.com/guigolab/CLS3_GENCODE/tree/main/data_release#gencode-cls3-mappings). Results are available in the folder for [transcripts](https://github.com/guigolab/CLS3_GENCODE/blob/main/complementary_data/gencode_byotypes_datasets/v47_CLS3_mapping.transcripts) and [genes](https://github.com/guigolab/CLS3_GENCODE/blob/main/complementary_data/gencode_byotypes_datasets/v47_CLS3_mapping.genes). 

```
pending
```

## 2. GENCODE reference
The universe of transcripts and genes for comparison has been defined as follow.

### protein-coding: 89,832 transcripts - 19,744 genes
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz
zcat gencode.v47.primary_assembly.annotation.gtf.gz | awk -F'\t' '$1 ~ /^chr([0-9]+|[XYM])$/' | awk -F'\t' '$3=="transcript" && $9 ~ /gene_type "protein_coding"/ && $9 ~ /transcript_type "protein_coding"/' | cut -d"\"" -f4 | sort -u > protein_coding.transcripts.v47.ids
zcat gencode.v47.primary_assembly.annotation.gtf.gz | awk -F'\t' '$1 ~ /^chr([0-9]+|[XYM])$/' | awk -F'\t' '$3=="transcript" && $9 ~ /gene_type "protein_coding"/ && $9 ~ /transcript_type "protein_coding"/' | cut -d"\"" -f2 | sort -u > protein_coding.loci.v47.ids
```

### lncRNAs: 26,709 transcripts - 14,680 genes
GENCODE v27 was used as reference given that was the geneset available at the time of the design.

```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.long_noncoding_RNAs.gtf.gz
zcat gencode.v27.long_noncoding_RNAs.gtf.gz | awk -F'\t' '$1 ~ /^chr([0-9]+|[XY])$/' | awk -F'\t' '$3=="transcript" && !($9 ~ /TEC/)' | cut -d"\"" -f4 | sort -u > lncRNA.transcripts.v27.ids 
zcat gencode.v27.long_noncoding_RNAs.gtf.gz | awk -F'\t' '$1 ~ /^chr([0-9]+|[XY])$/' | awk -F'\t' '$3=="transcript" && !($9 ~ /TEC/)' | cut -d"\"" -f2 | sort -u > lncRNA.loci.v27.ids 
```

## 3. Background
For some analyses, we have also employed a set of decoy models (17,223 genes, 85,283 transcripts) that attempt to mimic the background (non-genic) behavior of the genome.
Find them [here](https://github.com/guigolab/CLS3_GENCODE/tree/main/complementary_data/decoy_models).

```
cut -d"\"" -f4 random_replicates_locirelocation.original.gtf | sort -u > decoys.transcripts.ids
cut -d"\"" -f2 random_replicates_locirelocation.original.gtf | sort -u > decoys.loci.ids
```

The creation of the list of TSS list from these files is documented [here](https://github.com/guigolab/CLS3_GENCODE/tree/main/complementary_data/tss_representatives).

# Disjoint sets
From these sets, non-overlapping subset have been extracted for analyses that would have otherwise been confounded.

### Gather transcripts annotation in BED format
```
zgrep -wFf protein_coding.transcripts.v47.ids gencode.v47.primary_assembly.annotation.gtf.g | awk -F"\t" -v OFS="\t" '$3 == "transcript" { split($9, tags, "\""); print $1,$4,$5,tags[2],tags[4] }' | sort --parallel=4 -k1,1 -k4,4n > proteincodingv47.bed
zgrep -wFf lncRNA.transcripts.v27.ids gencode.v27.primary_assembly.annotation.gtf.gz | awk -F"\t" -v OFS="\t" '$3 == "transcript" { split($9, tags, "\""); print $1,$4,$5,tags[2],tags[4] }' | sort --parallel=4 -k1,1 -k2,2n > lncRNAv27.bed
grep -wFf decoys.transcripts.ids random_replicates_locirelocation.original.gtf | awk -F"\t" -v OFS="\t" '$3 == "transcript" { split($9, tags, "\""); print $1,$4,$5,tags[2],tags[4] }' | sort --parallel=4 -k1,1 -k4,4n > decoy.bed
zgrep -wFf cls.transcripts.ids gencode.v47.primary_assembly.annotation.gtf.gz | awk -F"\t" -v OFS="\t" '$3 == "transcript" { split($9, tags, "\""); print $1,$4,$5,tags[2],tags[4] }' | sort --parallel=4 -k1,1 -k4,4n > intergenicCLS.bed 
```

### protein-coding: 67,119 transcripts - 13,883 genes
```
bedtools intersect -a proteincodingv47.bed -b lncRNAv27.bed decoy.bed intergenicCLS.bed -v | cut -f5 > protein_coding.transcripts.v47.disjoint.ids
grep -vFf <(bedtools intersect -a proteincodingv47.bed -b lncRNAv27.bed decoy.bed intergenicCLS.bed -wa | cut -f4) protein_coding.v47.loci.ids > protein_coding.loci.v47.disjoint.ids
```

### lncRNAs: 14,908 transcripts - 7,774 genes
```
bedtools intersect -a lncRNAv27.bed -b proteincodingv47.bed decoy.bed intergenicCLS.bed -v | cut -f5 > lncRNA.transcripts.v27.disjoint.id
grep -vFf <(bedtools intersect -a lncRNAv27.bed -b proteincodingv47.bed decoy.bed intergenicCLS.bed -wa | cut -f4) lncRNA.v27.loci.ids > lncRNA.loci.v27.disjoint.ids
```

## cls: XXX transcripts - 8,706 genes
```
bedtools intersect -a intergenicCLS.bed -b proteincodingv47.bed decoy.bed lncRNAv27.bed -v | cut -f5 > cls.transcripts.disjoint.ids
awk '$9 == "Intergenic"' v47-CLS3_mapping_status | cut -f1 | sort -u > cls.loci.ids
```

## decoys: 84,063 transcripts - 17,005 genes
As they were designed in the intergenic space of v27, a minimal number of decoy genes are discarded in this stage.
```
bedtools intersect -a decoy.bed -b intergenicCLS.bed proteincodingv47.bed lncRNAv27.bed -v | cut -f5 > decoy.transcripts.disjoint.ids
grep -vFf <(bedtools intersect -a decoy.bed -b proteincodingv47.bed lncRNAv27.bed intergenicCLS.bed -wa | cut -f4) decoy.loci.ids > decoy.loci.disjoint.ids
```
