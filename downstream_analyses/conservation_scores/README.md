# Sequence conservation across mammals. 
We used mammalian conservation of human lncRNAs as a metric to compare the GENCODE novel CLS annotations with previously existing lncRNAs and protein-coding genes. We obtained PhyloP scores for the Zoonomia 241-way mammalian alignment [^1,^2] from the UCSC Genome Browser[^3]. Since we aimed to compare the nature of GENCODE annotations rather than examine precise conservation levels, we evaluate conservation per-transcript, for which we computed mean exon and splice-junction PhyloP scores.

Transcripts lacking PhyloP scores for more than 10% of the exon bases or any of its splice-junction were excluded. We set the neutrally evolving range based on the PhyloP score distribution from -1.0 to 1.0. To avoid confounding conservation signals from protein-coding genes, we used the [disjoint set of annotations described at _complementary_data_](https://github.com/guigolab/CLS3_GENCODE/tree/main/complementary_data/gencode_byotypes_datasets#disjoint-sets). The sets of conserved transcripts (either by sequence or splice-junctions conservation) are released here and detailed below. A gene is considered conserved if one of its transcript is conserved.

### protein-coding: 86,952 transcripts
Of which xxx show conservation signal for their exonic sequence, while xx have conserved splice-junctions. These correspond to 18,691 genes (94.7% of the total genes present in the disjoint set).

### lncRNAs: 5,944 transcripts
Of which xxx show conservation signal for their exonic sequence, while xx have conserved splice-junctions. These correspond to 3,144 genes (40.4% of the total genes present in the disjoint set).

### cls: 2,210 transcripts
Of which xxx show conservation signal for their exonic sequence, while xx have conserved splice-junctions. These correspond to 1,230 genes (14.6% of the total genes present in the disjoint set).

### decoys: 2,113 transcripts
Of which xxx show conservation signal for their exonic sequence, while xx have conserved splice-junctions. These correspond to 1,015 genes (5.9% of the total genes present in the disjoint set).

[^1] https://doi.org/10.1126/science.abn3943
[^2] https://doi.org/10.1038/s41586-020-2876-6
[^3]: Raney, B. J. et al. The UCSC Genome Browser database: 2024 update. Nucleic Acids Res. 52, D1082–D1088 (2024).


# data sets:

These id file list the conserved transcript or loci from:
  complementary_data/gencode_byotypes_datasets/*.ids

Files for different metrics for conservation are available:
  - both - splice junctions and exons pass conservation metric
  - either - splice junctions or exons pass conservation metric
  - sj - splice junctions pass conservation metric, exons can be either
  - exon - exons pass conservation metric, splice junctions can be either
