# CLS3 data release
## Files description
Summary of files generated from the CLS3 data:
    <li><a href="#CLS-anchored-models">CLS anchored models</a></li>
    <li><a href="#CLS-transcripts">CLS transcripts</a></li>
    <li><a href="#CLS-loci">CLS loci</a>
    <li><a href="#GENCODE-CLS3-Mappings">GENCODE-CLS3 Mappings</a>
    <li><a href="#Extended-GENCODE-v47">Extended GENCODE v47</a>
    <li><a href="#Target-files">Target files</a>
    
![CLS workflow](https://github.com/user-attachments/assets/47c87440-8aca-44fc-a3bc-0bf27d5527f8)
 
## CLS anchored models

The CLS3 long-read data is processed through LyRic to obtaing high-confidence transcript models per sample through starting from pacBioSII- and ONT-sequenced long reads. Further, to collect a comprehensive set of transcripts
from the different tissues and technologies, the models were "anchored" according to their support at 5' and 3' ends ([anchorTranscriptsEnds.pl](https://github.com/guigolab/LyRic/blob/master/utils/anchorTranscriptsEnds.pl)), and merged together ([tmerge](https://github.com/guigolab/tmerge)). 
In this way we reduce the redundancy in the dataset while also preserving potential alternative start and end sites. 

The final set of transcript models, anchor-merged across samples, for human and mouse can be downloaded at:
  - [Human](https://zenodo.org/records/13946596/files/Hv3_masterTable_refined.gtf.gz?download=1)
  - [Mouse](https://zenodo.org/records/13946596/files/Mv2_masterTable_refined.gtf.gz?download=1)
  
The attribute tags description can be found [here](gtf_tags_explained.md).

 ## CLS transcripts
 To further reduce the redundancy of the aforementioned transcripts collection, models sharing the same structure were merged together, ignoring eventual variation at the terminal exons, and therefore neglecting end-support information. 
 
 Different strategies were followed for spliced and monoexonic transcripts.
 1. For the _spliced transcripts_, those with identical intron chains were merged into a single one, assigning the furthest start and end sites as TSS and TTS.
 2. For the _monoexonic transcripts_, they were merged together in case they share more that 50% overlap with each other.
   
 The chain GTF for human and mouse can be downloaded using the following links:
   - [Human - spliced](https://zenodo.org/records/13946596/files/Hv3_splicedmasterTable_refined.gtf.gz?download=1)
   - [Human - monoexonic](https://zenodo.org/records/13946596/files/Hv3_unsplicedmasterTable_refined.gtf.gz?download=1)
   - [Mouse - spliced](https://zenodo.org/records/13946596/files/Mv2_splicedmasterTable_refined.gtf.gz?download=1)
   - [Mouse - monoexonic](https://zenodo.org/records/13946596/files/Mv2_unsplicedmasterTable_refined.gtf.gz?download=1)

Get an [overview](https://github.com/guigolab/CLS3_GENCODE/tree/main/data_release/overview) of the transcripts distribution across stages, tissues, and technology.

## CLS loci
With the intent of grouping together different models in uniquely identifiable loci, we clustered CLS transcripts into regions of continuous transcription.
Eventually, transcripts sharing any overlap on the same strand have been brought together into a single locus, preserving their structure. 

The loci GTFs can be downloaded at: 
  - [Human - gencode v27 tagged](https://zenodo.org/records/13946596/files/Hv3_masterTable_refined_+withinTmerge_gencodev27_tagged.loci.gtf.gz?download=1)
  - [Mouse - gencode vM16 tagged](https://zenodo.org/records/13946596/files/Mv2_masterTable_refined_+withinTmerge_gencodevM16_tagged.loci.gtf.gz?download=1)

## GENCODE-CLS3 Mappings
The latest human and mouse annotations have reported a huge increase in the overall number of lncRNA genes and transcripts, together with a great improvement in the annotation of already exising lncRNA genes, thanks to the incorporation of CLS3 models. Such transcripts have been processed by the HAVANA team of manual annotators at EBI, and the refinement and augmentation of GENCODE v47 and vM36 annotation attributable to CLS data have been documented here:

  - [v47-CLS3_extended_mappings] (coming soon on Zenodo)
  - [vM36-CLS3_extended_mappings] (coming soon on Zenodo)

The original numbers released in the publication, can be obtained by quierying the files as follow;
```
#Total novel lncRNA genes (17,931):
awk -F "\t" '$9 ~ /created_gene/ && $2 ~ /CLS3_created/ && $5 ~ /lncRNA/' v47-CLS3_extended_mappings | cut -f8 | sort -u | wc -l

#Total novel lncRNA transcripts (140,268):
awk -F "\t" '$2 ~ /CLS3_created/ && $5 ~ /lncRNA/' v47-CLS3_extended_mappings | cut -f1 | sort -u | wc -l
```

The novelty brought by these transcripts, with respect to version v27 (human) and vM16 (mouse) are reported the following table, and assigned asdefined in the schema below.

  - [v47-CLS3 mapping](https://zenodo.org/records/13946596/files/v47-CLS3mapping_status.txt?download=1)
  - [vM36-CLS3 mapping] (coming soon on Zenodo)

The mapping across v47 ENSTs and the CLS3 anchICs they were extended/created from, with added details like novelty at the transcript as well as gene level. 

For each transcript (ENST) created/extended due to CLS3 (anchICs), the file lists: <br />
||||
|-|-|-|
**geneID_v47** | v47 gene (ENSG) that the transcript belongs to <br />
**transcriptID_v47** | v47 transcript ID (ENST) <br />
**created/extended** | tag specifying whether the transcript was created or extended using TAGENE/manually <br />
**CLS3_anchIC** | CLS3 anchIC(s) that led to the addition of the transcript to v47 <br />
**CLS3_anchIC_gffComparev27** | gffcompare classification for the anchIC(s) w.r.t. v27 (reference annotation) <br />
**v47-CLS3_mappingTag** | states the mapping strategy used; details in the above section <br />
**v47_biotype** | v47 biotype <br />
**transcriptClassification** | transcript (ENST) novelty status taking into account the different gffcompare classifications from all the underlying anchICs <br />
**geneClassification** | gene (ENSG) novelty status taking into account the different gffcompare classifications from all the underlying transcripts <br />
**CLS3_anchTM** | CLS3 anchTM(s) that led to the addition of the transcript to v47. Mapped through the anchICs. <br />

The gffcompare novelty status definitions w.r.t. v27 for the anchICs, ENSTs and ENSGs.
![Mapping](https://github.com/user-attachments/assets/7bbfea20-27d5-4bf5-8a55-494e5991943b)
<br />

## Extended GENCODE v47
For the sole purpose of analysing those CLS transcripts that were not yet incorporated into the annotation, for some analyses we leveraged an extended version of the GENCODE v47.
This refers to a tailored-made GTF file enhanced by adding CLS3 loci built from non artefactual spliced transcripts.

This can be downloaded here:
  - [Human v47 extended annotation](https://zenodo.org/records/13946596/files/enhanced_annotation_v47.refined.gtf.gz?download=1)
    
## Target files
Targets used for probe design:
- [Human targeted regions](https://zenodo.org/records/13946596/files/hs.allNonPcgTargetsMerged.targets.gtf.gz?download=1)
- [Mouse targeted regions](https://zenodo.org/records/13946596/files/mm.allNonPcgTargetsMerged.targets.gtf.gz?download=1)

For the target design, 8 catalogs (CMfinderCRSs, GWAScatalog, UCE, VISTAenhancers, fantomCat, fantomEnhancers, bigTranscriptome, miTranscriptome) were liftedOver from human to mouse. Mappings for all such targets can be found here:
  - [LiftedOver targets](https://zenodo.org/records/13946596/files/final.liftedOverTargets.mapping.txt?download=1)
