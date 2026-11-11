# Single-cell data analyses

## Prerequisites
Data has been gathered as described at [Human Common Cell Atlas](https://github.com/cellatlas/human)[^1] and processed using the set of kallisto, bustools and kb-python algorithms[^2][^3] with GENCODE v47 as reference (see Snakefile_cls3_GENCODEv47).

In this folder, the scripts used to process the matrices and filter out low quality cells (_process_matrices.R_), evaluate intronic-read fraction (_check_intronic.R_), and flag lower-quality/outlier clusters (_flag_clusters.R_).
Analyses on clustering and cell mapping across different solutions (_cluster_jaccard.R_), clustering agreement evaluation (_tau_clustering.R_), and marker discovery (_get_markers.R_) are also available.

[^1]: Ángel Galvez-Merchán, A. Sina Booeshaghi, Lior Pachter, A human commons cell atlas reveals cell type specificity for OAS1 isoforms, bioRxiv 2024.03.23.586412; doi: https://doi.org/10.1101/2024.03.23.586412
[^2]: Sullivan DK, Min KHJ, Hjörleifsson KE, Luebbert L, Holley G, Moses L, Gustafsson J, Bray NL, Pimentel H, Booeshaghi AS, Melsted P, Pachter L. kallisto, bustools and kb-python for quantifying bulk, single-cell and single-nucleus RNA-seq. Nat Protoc. 2025 Mar;20(3):587-607. doi: 10.1038/s41596-024-01057-0. Epub 2024 Oct 10. PMID: 39390263.
[^3]: Delaney K Sullivan, Kristján Eldjárn Hjörleifsson, Nikhila P Swarna, Conrad Oakes, Guillaume Holley, Páll Melsted, Lior Pachter, Accurate quantification of nascent and mature RNAs from single-cell and single-nucleus RNA-seq, Nucleic Acids Research, Volume 53, Issue 1, 13 January 2025, gkae1137, https://doi.org/10.1093/nar/gkae1137
