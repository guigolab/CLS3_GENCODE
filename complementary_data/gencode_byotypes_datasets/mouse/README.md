# Datasets used for analysis
The same criteria and filtering designed for Human were applied in Mouse.

## 1. Novel CLS models
In GENCODE vM16, a total of 142,428 transcripts (listed [here](https://zenodo.org/records/15004659/files/vM36-CLS3mapping_status.txt?download=1)) have been either created (136,169 lncRNAs and 137 others) or modified (6,095 lncRNAs and 27 others) thanks to CLS data. Each of these transcripts and genes have been assigned a novelty category as detailed in [GENCODE-CLS3 Mappings](https://github.com/guigolab/CLS3_GENCODE/tree/main/data_release#gencode-cls3-mappings).

### Novel transcripts: 142,217 transcripts - 22,779 genes
### CLS loci: 25,955 transcripts - 10,143 genes
CLS transcripts and loci are a subset of the total novel genes, corresponding to transcripts and loci intergenic with respect to GENCODE vM16 (the reference at the time of the design).

```
wget https://zenodo.org/records/13946596/files/vM36-CLS3mapping_status.txt?download=1 && mv vM36-CLS3mapping_status.txt\?download\=1 vM36-CLS3_status
zgrep -Ff <(awk '$9 == "Intergenic"' vM36-CLS3_status | cut -f1 | sort -u) gencode.vM36.primary_assembly.annotation.gtf.gz | awk -F'\t' '$3 == "transcript" && $1 ~ /^chr([0-9]+|[XYM])$/' | cut -d"\"" -f4 | sort -u > cls.transcripts.ids
zgrep -Ff <(awk '$9 == "Intergenic"' vM36-CLS3_status | cut -f1 | sort -u) gencode.vM36.primary_assembly.annotation.gtf.gz | awk -F'\t' '$3 == "transcript" && $1 ~ /^chr([0-9]+|[XYM])$/' | cut -d"\"" -f2 | sort -u > cls.loci.ids
```

## 2. GENCODE reference
The universe of transcripts and genes for comparison has been defined as follow.

### protein-coding: 58,547 transcripts - 21,526 genes
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_M36/gencode.vM36.primary_assembly.annotation.gtf.gz
zcat gencode.vM36.primary_assembly.annotation.gtf.gz | awk -F'\t' '$1 ~ /^chr([0-9]+|[XYM])$/' | awk -F'\t' '$3=="transcript" && $9 ~ /gene_type "protein_coding"/ && $9 ~ /transcript_type "protein_coding"/' | cut -d"\"" -f4 | sort -u > protein_coding.transcripts.vM36.ids
zcat gencode.vM36.primary_assembly.annotation.gtf.gz | awk -F'\t' '$1 ~ /^chr([0-9]+|[XYM])$/' | awk -F'\t' '$3=="transcript" && $9 ~ /gene_type "protein_coding"/ && $9 ~ /transcript_type "protein_coding"/' | cut -d"\"" -f2 | sort -u > protein_coding.loci.vM36.ids
```

### lncRNAs: 14,094 transcripts - 9,279 genes
GENCODE vM16 was used as reference being the geneset available at the time of the design.

```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_M16/gencode.vM16.long_noncoding_RNAs.gtf.gz
zcat gencode.vM16.long_noncoding_RNAs.gtf.gz | awk -F'\t' '$1 ~ /^chr([0-9]+|[XY])$/' | awk -F'\t' '$3=="transcript" && !($9 ~ /TEC/)' | cut -d"\"" -f4 | sort -u > lncRNA.transcripts.vM16.ids 
zcat gencode.vM16.long_noncoding_RNAs.gtf.gz | awk -F'\t' '$1 ~ /^chr([0-9]+|[XY])$/' | awk -F'\t' '$3=="transcript" && !($9 ~ /TEC/)' | cut -d"\"" -f2 | sort -u > lncRNA.loci.vM16.ids 
```

