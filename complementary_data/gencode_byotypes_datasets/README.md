# Datasets used for analysis
Accurate gene annotations are fundamental to functionally interpreting the activity (transcription, chromatin modifications, folding, protein binding) and the sequence variation of genomes. We wanted to show how the incorporation of the CLS data into GENCODE greatly enhances the functional interpretability of the human genome. 

## 1. Novel CLS models
In GENCODE v47, a total of 151,618 transcripts (listed [here](https://zenodo.org/records/15004659/files/v47-CLS3mapping_status.txt?download=1)) have been either created (140,268 lncRNAs and 293 others) or modified (10,852 lncRNAs and 205 others) thanks to CLS data. Each of these transcripts and genes have been assigned a novelty category as detailed in [GENCODE-CLS3 Mappings](https://github.com/guigolab/CLS3_GENCODE/tree/main/data_release#gencode-cls3-mappings).

### Novel transcripts: 149,010 transcripts - 17,777 genes
In these analyses, novel transcripts are restricted to lncRNAs transcripts, including those extending previously annotated ones, while novel loci are instead only those genes introduced as a consequence of CLS data. We also restricted the set to loci on main chromosomes only.

```
wget https://zenodo.org/records/15004659/files/CLS3_transcripts_in_v47.all_biotypes.chr.gencode_versions.genes.txt?download=1 && mv CLS3_transcripts_in_v47.all_biotypes.chr.gencode_versions.genes.txt\?download\=1 v47-CLS3_extended_mappings
awk -F "\t" '$5 ~ /lncRNA/ && $6 ~ /^([0-9]+|[XYM])$/' v47-CLS3_extended_mappings | cut -f1 | sort -u > novel.transcripts.ids
awk -F "\t" '$9 ~ /created_gene/ && $2 ~ /CLS3_created/ && $5 ~ /lncRNA/ && $6 ~ /^([0-9]+|[XYM])$/' v47-CLS3_extended_mappings | cut -f8 | sort -u > novel.loci.ids
```

### CLS loci: 19,473 transcripts - 8,576 genes
CLS transcripts and loci are a subset of the total novel genes, corresponding to transcripts and loci intergenic with respect to GENCODE v27 (the reference at the time of the design).

```
wget https://zenodo.org/records/13946596/files/v47-CLS3mapping_status.txt?download=1 && mv v47-CLS3mapping_status.txt\?download\=1 v47-CLS3_status
zgrep -Ff <(awk '$9 == "Intergenic"' v47-CLS3_status | cut -f1 | sort -u) gencode.v47.primary_assembly.annotation.gtf.gz | awk -F'\t' '$3 == "transcript" && $1 ~ /^chr([0-9]+|[XYM])$/' | cut -d"\"" -f4 | sort -u > cls.transcripts.ids
zgrep -Ff <(awk '$9 == "Intergenic"' v47-CLS3_status | cut -f1 | sort -u) gencode.v47.primary_assembly.annotation.gtf.gz | awk -F'\t' '$3 == "transcript" && $1 ~ /^chr([0-9]+|[XYM])$/' | cut -d"\"" -f2 | sort -u > cls.loci.ids
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
GENCODE v27 was used as reference being the geneset available at the time of the design.

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

The creation of the list of TSS list from these files is documented [here](https://github.com/guigolab/CLS3_GENCODE/tree/main/complementary_data/tss_representatives). For novel TSSs, we rescued the novel transcripts belonging to biotypes other than lncRNAs (151,618 total transcripts). 

# Disjoint sets
From these sets, non-overlapping subset have been extracted for analyses that would have otherwise been confounded.

### Gather transcripts annotation in BED format
```
zgrep -wFf protein_coding.transcripts.v47.ids gencode.v47.primary_assembly.annotation.gtf.gz | awk -F"\t" -v OFS="\t" '$3 == "transcript" { split($9, tags, "\""); print $1,$4,$5,tags[2],tags[4] }' | sort --parallel=4 -k1,1 -k4,4n > proteincodingv47.bed
zgrep -wFf lncRNA.transcripts.v27.ids gencode.v27.long_noncoding_RNAs.gtf.gz | awk -F"\t" -v OFS="\t" '$3 == "transcript" { split($9, tags, "\""); print $1,$4,$5,tags[2],tags[4] }' | sort --parallel=4 -k1,1 -k2,2n > lncRNAv27.bed
grep -wFf decoys.transcripts.ids random_replicates_locirelocation.original.gtf | awk -F"\t" -v OFS="\t" '$3 == "transcript" { split($9, tags, "\""); print $1,$4,$5,tags[2],tags[4] }' | sort --parallel=4 -k1,1 -k4,4n > decoy.bed
zgrep -wFf cls.transcripts.ids gencode.v47.primary_assembly.annotation.gtf.gz | awk -F"\t" -v OFS="\t" '$3 == "transcript" { split($9, tags, "\""); print $1,$4,$5,tags[2],tags[4] }' | sort --parallel=4 -k1,1 -k4,4n > intergenicCLS.bed 
```

### protein-coding: 67,119 transcripts - 13,883 genes
```
bedtools intersect -a proteincodingv47.bed -b lncRNAv27.bed decoy.bed intergenicCLS.bed -v | cut -f5 > protein_coding.transcripts.v47.disjoint.ids
grep -vFf <(bedtools intersect -a proteincodingv47.bed -b lncRNAv27.bed decoy.bed intergenicCLS.bed -wa | cut -f4) protein_coding.loci.v47.ids > protein_coding.loci.v47.disjoint.ids
```

### lncRNAs: 14,908 transcripts - 7,774 genes
```
bedtools intersect -a lncRNAv27.bed -b proteincodingv47.bed decoy.bed intergenicCLS.bed -v | cut -f5 > lncRNA.transcripts.v27.disjoint.ids
grep -vFf <(bedtools intersect -a lncRNAv27.bed -b proteincodingv47.bed decoy.bed intergenicCLS.bed -wa | cut -f4) lncRNA.loci.v27.ids > lncRNA.loci.v27.disjoint.ids
```

### cls: 19,140 transcripts - 8,392 genes
```
bedtools intersect -a intergenicCLS.bed -b proteincodingv47.bed decoy.bed lncRNAv27.bed -v | cut -f5 > cls.transcripts.disjoint.ids
grep -vFf <(bedtools intersect -a intergenicCLS.bed -b proteincodingv47.bed lncRNAv27.bed decoy.bed -wa | cut -f4 | cut -d"." -f1) cls.loci.ids > cls.loci.disjoint.ids
```

### decoys: 84,063 transcripts - 17,005 genes
As they were designed in the intergenic space of v27, a minimal number of decoy genes are discarded in this stage.
```
bedtools intersect -a decoy.bed -b intergenicCLS.bed proteincodingv47.bed lncRNAv27.bed -v | cut -f5 > decoys.transcripts.disjoint.ids
grep -vFf <(bedtools intersect -a decoy.bed -b proteincodingv47.bed lncRNAv27.bed intergenicCLS.bed -wa | cut -f4) decoys.loci.ids > decoys.loci.disjoint.ids
```
