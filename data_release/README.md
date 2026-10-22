# CLS3 data release
## Files description
Summary of files generated from the CLS3 data:
    <li><a href="#CLS-anchored-models">CLS anchored models</a></li>
    <li><a href="#CLS-transcripts">CLS transcripts</a></li>
    <li><a href="#GENCODE-CLS3-Mappings">GENCODE-CLS Mappings</a></li>
    <li><a href="#Target-files">Target files</a></li>
    
<img width="794" height="595" alt="image" src="https://github.com/user-attachments/assets/d3f081a0-5b96-482a-9053-477dcfb45f1c" />
 
## CLS anchored models

The CLS3 long-read data is processed through LyRic to obtaing high-confidence transcript models per sample. Further, to collect a comprehensive set of transcripts from the different tissues and technologies, the models were "anchored" according to their support at 5' and 3' ends ([anchorTranscriptsEnds.pl](https://github.com/guigolab/LyRic/blob/master/utils/anchorTranscriptsEnds.pl)), and merged together ([tmerge](https://github.com/guigolab/tmerge)). In this way we reduce the redundancy in the dataset while also preserving potential alternative start and end sites. 

The final set of transcript models, anchor-merged across samples, for human and mouse can be downloaded at:
  - [Human](https://zenodo.org/records/15004659/files/Hv3_masterTable_refined.gtf.gz?download=1)
  - [Mouse](https://zenodo.org/records/15004659/files/Mv2_masterTable_refined.gtf.gz?download=1)
  
The attribute tags description can be found [here](gtf_tags_explained.md).

 ## CLS transcripts
To further reduce the redundancy of the aforementioned transcripts collection, models sharing the same structure were merged together, ignoring eventual variation at the terminal exons, and therefore neglecting end-support information. 
 
 Different strategies were followed for spliced and monoexonic transcripts.
 1. For the _spliced transcripts_, those with identical intron chains were merged into a single one, assigning the furthest start and end sites as TSS and TTS.
 2. For the _monoexonic transcripts_, when sharing more that 50% overlap with each other, they were merged together.
   
 The chain GTF for human and mouse can be downloaded using the following links:
   - [Human - spliced](https://zenodo.org/records/15004659/files/Hv3_splicedmasterTable_refined.gtf.gz?download=1)
   - [Human - monoexonic](https://zenodo.org/records/15004659/files/Hv3_unsplicedmasterTable_refined.gtf.gz?download=1)
   - [Mouse - spliced](https://zenodo.org/records/15004659/files/Mv2_splicedmasterTable_refined.gtf.gz?download=1)
   - [Mouse - monoexonic](https://zenodo.org/records/15004659/files/Mv2_unsplicedmasterTable_refined.gtf.gz?download=1)


## GENCODE-CLS3 Mappings
The latest human and mouse annotations have reported a huge increase in the overall number of lncRNA genes and transcripts, together with a great improvement in the annotation of already exising lncRNA genes, thanks to the incorporation of CLS3 models. Such transcripts have been processed by the HAVANA team of manual annotators at EBI, and the refinement and augmentation of GENCODE v47 and vM36 annotation attributable to CLS data have been documented here:

  - [v47-CLS3_extended_mappings](https://zenodo.org/records/15004659/files/CLS3_transcripts_in_v47.all_biotypes.chr.gencode_versions.genes.txt?download=1)
  - [vM36-CLS3_extended_mappings](https://zenodo.org/records/15004659/files/CLS3_transcripts_in_vM36.all_biotypes.chr.gencode_versions.genes.txt?download=1)

The original numbers released in the publication, can be obtained by querying the files as follow;
```
#Total novel lncRNA genes (17,931):
awk -F "\t" '$9 ~ /created_gene/ && $2 ~ /CLS3_created/ && $5 ~ /lncRNA/' v47-CLS3_extended_mappings | cut -f8 | sort -u | wc -l

#Total novel lncRNA transcripts (140,268):
awk -F "\t" '$2 ~ /CLS3_created/ && $5 ~ /lncRNA/' v47-CLS3_extended_mappings | cut -f1 | sort -u | wc -l
```

## Target files
Probes were designed to target poorly annotated non-coding genomic regions, including lncRNAs, small non-coding RNAs, enhancers, conserved RNA structures, GWAS-associated variants, regions with predicted coding potential, and evolutionarily conserved sequences. Using mainly GENCODE v27 (hg38) and multiple genomic databases, target regions were collected, generally extended by 200 bp, and merged. Regions overlapping protein-coding genes, repeats, or mitochondrial/alternative chromosomes were excluded. Highly expressed regions (RPKM ≥ 50) were also removed using GTEx RNA-seq data from blood, brain, heart, and liver. A similar approach was used for mouse, using GENCODE vM16 (mm10) and leveraging ENCODE data for expression testing.

The final targets used for probe design:
- [Human targeted regions](https://zenodo.org/records/15004659/files/hs.allNonPcgTargetsMerged.targets.gtf.gz?download=1)
- [Mouse targeted regions](https://zenodo.org/records/15004659/files/mm.allNonPcgTargetsMerged.targets.gtf.gz?download=1)

For mouse, 8 catalogs (CMfinderCRSs, GWAScatalog, UCE, VISTAenhancers, fantomCat, fantomEnhancers, bigTranscriptome, miTranscriptome) were liftedOver from human. Mappings for all such targets can be found here:
  - [LiftedOver targets](https://zenodo.org/records/15004659/files/final.liftedOverTargets.mapping.txt?download=1)
