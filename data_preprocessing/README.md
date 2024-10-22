# Data Preprocessing and Quality Assessment

## Preprocessing

### ONT
Basecalling for ONT was performed using Guppy v6 SUP. NanoPlot was utilized to generate metrics, including read length distributions, quality scores (Q-scores), and the total number of reads, providing insights into the overall performance of the sequencing run. Additionally, the split_on_adapter utility from the Duplex Tools suite was employed to evaluate the presence of concatamers in the data. See the [resolving_concatamers](https://github.com/guigolab/CLS3_GENCODE/blob/main/data_preprocessing/resolving_concatamers.md) documentation for a more detailed description. 

### PacBio
To make the PacBio data compatible with the downstream LyRic processing, PacBio FASTQ files containing CCS (Circular Consensus Sequencing) reads were generated using the [pb_gen](https://github.com/guigolab/pb_gen) workflow.

Raw PacBio and ONT reads are available through ArrayExpress accession E-MTAB-14562.

## Quality Assessment

### Read length distribution
Read length per sample was addressed using the script:
```
./get_read_length.sh
```

### SQANTI Evaluation
The accuracy and quality of the models have been assessed using SQANTI3[^27], in comparison to GENCODE v27 and vM16, for [human](https://guigolab.github.io/CLS3_GENCODE/SQANTI_reports/Human_CLStranscripts_v27.html) and [mouse](https://guigolab.github.io/CLS3_GENCODE/SQANTI_reports/Mouse_CLStranscripts_vM16.html), respectively. Run parameters as follow:
```
sqanti command line
```

[^27]: Pardo-Palacios, F. J. et al. SQANTI3: curation of long-read transcriptomes for accurate identification of known and novel isoforms. Nat. Methods (2024) doi:10.1038/s41592-024-02229-2.
