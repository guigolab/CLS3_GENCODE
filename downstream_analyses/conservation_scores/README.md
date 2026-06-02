# Sequence conservation across mammals. 
We used mammalian conservation of human lncRNAs as a metric to compare the GENCODE novel CLS annotations with previously existing lncRNAs and protein-coding genes. We obtained PhyloP scores for the Zoonomia 241-way mammalian alignment [^1],[^2] from the UCSC Genome Browser[^3]. Since we aimed to compare the nature of GENCODE annotations rather than examine precise conservation levels, we evaluate conservation per-transcript, for which we computed mean exon and splice-junction PhyloP scores.

Transcripts lacking PhyloP scores for more than 10% of the exon bases or any of its splice-junction were excluded. We set the neutrally evolving range based on the PhyloP score distribution from -1.0 to 1.0. To avoid confounding conservation signals from protein-coding genes, we used the [disjoint set of annotations described at _complementary_data_](https://github.com/guigolab/CLS3_GENCODE/tree/main/complementary_data/gencode_byotypes_datasets#disjoint-sets). The sets of conserved _disjoint_ transcripts (either by sequence or splice-junctions conservation) are released here and detailed below. A gene is considered conserved if one of its transcript is conserved.

### protein-coding: 54,643 transcripts (96%)
Of which 46,701 show conservation signal for their exonic sequence, while 53,705 have conserved splice-junctions. These correspond to 12,989 genes (94% of the total genes present in the disjoint set).

### lncRNAs: 1,924 transcripts (14.7%)
Of which 819 show conservation signal for their exonic sequence, while 1,621 have conserved splice-junctions. These correspond to 1,010 genes (12% of the total genes present in the disjoint set).

### CLS: 2,146 transcripts (11.2%)
Of which 1,126 show conservation signal for their exonic sequence, while 1,522 have conserved splice-junctions. These correspond to 1,196 genes (14% of the total genes present in the disjoint set).

### Decoys: 2,074 transcripts (2.5%)
Of which 769 show conservation signal for their exonic sequence, while 1,622 have conserved splice-junctions. These correspond to 1,001 genes (5.9% of the total genes present in the disjoint set).

# Data sets:
These id file list the conserved transcript or loci from _complementary_data/gencode_byotypes_datasets/*.ids_

Files for different metrics for conservation are available:
<li> both - splice junctions and exons pass conservation metric </li>
<li> either - splice junctions or exons pass conservation metric </li>
<li> sj - splice junctions pass conservation metric, exons can be either </li>
<li> exon - exons pass conservation metric, splice junctions can be either </li>

[^1]: https://doi.org/10.1126/science.abn3943
[^2]: https://doi.org/10.1038/s41586-020-2876-6
[^3]: Raney, B. J. et al. The UCSC Genome Browser database: 2024 update. Nucleic Acids Res. 52, D1082–D1088 (2024).

