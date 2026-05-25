# Datasets used for generating the representative TSSs:
The criteria used for selection of the representative TSSs is decribed [here]() _pending_.

## Representative TSS sets all transcripts: 
Please note that the gene numbers correspond to the transcripts, and are not equilant to the gene datasets used for some analyses.

### 1. protein coding TSS: 48,477 TSSs - 89,832 transcripts - 19,744 genes
TSS of [protein coding transcripts](https://github.com/guigolab/CLS3_GENCODE/blob/main/complementary_data/gencode_byotypes_datasets/README.md#protein-coding-89832-transcripts---19744-genes) from GENCODE v47 reference.

### 2. lncRNA TSS: 20,345 TSSs - 26,709 transcripts - 14,680 genes
TSS of [lncRNA transcripts](https://github.com/guigolab/CLS3_GENCODE/blob/main/complementary_data/gencode_byotypes_datasets/README.md#lncrnas-26709-transcripts---14680-genes) from GENCODE v27 reference.

### 3. decoy TSS: 60,143 TSSs - 85,283(*2) transcripts - 17,223(*2) genes
TSS of the [decoy transcripts]().
Decoy transcripts were generated in a strand-agnostic manner and therefore duplicated on both strands prior to representative TSS set generation, resulting in a doubling of the total decoy transcript count.

### 4. novel TSS: 52,355 TSSs - 149,392 transcripts - 26,946 genes
novel TSS of the entirity of transcripts added or extended in gencode v47 thanks to CLS.
These include the below described novel protein-coding (4.a) and novel lncRNA transcripts (4.b).

### 4a. novel protein coding TSS: 198 TSSs - 371 transcripts - 141 genes
novel TSS associated to protein-coding genes.

### 4b. novel lncRNA transcripts TSS: 52,146 TSSs - 149,010 transcripts - 26,794 genes
novel TSS associated to lncRNA genes. 
These transcripts are further subdivided into novel lncRNA (4.b.i) and CLS transcripts (4.b.ii).

### 4b.i. novel lncRNA TSS: 40,470 TSSs - 129,537 transcripts - 18,218 genes
TSS of the novel lncRNA transcripts (4.b), removing the intergenic CLS transcripts (4.b.ii). 

### 4b.ii. CLS transcripts TSS: 11,676 TSSs - 19,473 transcripts - 8,576 genes
novel TSS associated with CLS transcripts and loci intergenic w.r.t. GENCODE v27 (the reference at the time of the design)
The TSSs from CLS transcripts are subdivided into embryonic. adult, common and placenta derived TSS sets.


## Representative TSS sets disjoint transcripts:
### 1. protein coding TSS: 37,007 TSSs - 67,119 transcripts - 16,686 genes
TSS of [disjoint protein coding transcripts](https://github.com/guigolab/CLS3_GENCODE/blob/main/complementary_data/gencode_byotypes_datasets/README.md#protein-coding-67119-transcripts---13883-genes) from GENCODE v47 reference.

### 2.lncRNAs: 11,519 TSSs - 14,908 transcripts - 8,282 genes
TSS of [disjoint lncRNA transcripts](https://github.com/guigolab/CLS3_GENCODE/blob/main/complementary_data/gencode_byotypes_datasets/README.md#lncrnas-14908-transcripts---7774-genes) from GENCODE v27 reference.

### 3. decoys: 59416 TSSs - 84,063(*2) transcripts - 17,005(*2) genes
TSS of the [disjoint decoy transcripts](https://github.com/guigolab/CLS3_GENCODE/blob/main/complementary_data/gencode_byotypes_datasets/README.md#decoys-84063-transcripts---17005-genes).
Decoy transcripts were generated in a strand-agnostic manner and therefore duplicated on both strands prior to representative TSS set generation, resulting in a doubling of the total decoy transcript count.

### 4b.i. novel lncRNA: 11,644 TSSs - 34,835 transcripts - 5,710 genes
TSS of the [disjoint novel lncRNA transcripts](https://github.com/guigolab/CLS3_GENCODE/blob/main/complementary_data/gencode_byotypes_datasets/README.md#cls-19140-transcripts---8392-genes) (4.b), without the intergenic CLS transcripts (4.b.ii). 
Therefore these are the novel TSS w.r.t. lncRNAs v27.

### 4.b.ii. cls: 11,462 TSSs - 19,140 transcripts - 8,406 genes
novel TSS associated with [disjoint CLS transcripts](https://github.com/guigolab/CLS3_GENCODE/blob/main/complementary_data/gencode_byotypes_datasets/README.md#cls-19140-transcripts---8392-genes) and loci intergenic w.r.t. GENCODE v27 (the reference at the time of the design)
The TSSs from CLS transcripts are subdivided into embryonic. adult, common and placenta derived TSS sets.





