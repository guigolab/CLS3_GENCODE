#!/bin/bash
set -e
set -o pipefail

#----remove old data if any
rm -rf data/ plots/ stats/
#----download required files
mkdir -p data/ data/masterTable/ data/targets/ data/metadata/ stats/ plots/
echo -e "#### created needed directories ####"

for spec in Hv3 Mv2; do
       wget -q https://github.com/guigolab/gencode-cls-master-table/releases/latest/download/${spec}_masterTable_refined.gtf.gz -O data/masterTable/${spec}_masterTable_refined.gtf.gz

        #download CLS transcript (IC) files
      wget https://github.com/guigolab/gencode-cls-master-table/releases/download/GencodeCLS_v3.0/${spec}_splicedmasterTable_refined.gtf.gz -O data/masterTable/${spec}_splicedmasterTable_refined.gtf.gz
       wget https://github.com/guigolab/gencode-cls-master-table/releases/download/GencodeCLS_v3.0/${spec}_unsplicedmasterTable_refined.gtf.gz -O data/masterTable/${spec}_unsplicedmasterTable_refined.gtf.gz
       gunzip data/masterTable/${spec}*

       #download target catalog gtfs
	wget https://github.com/guigolab/gencode-cls-master-table/releases/download/Supplementary/${spec}_CLS3_targetDesign_mergedRegions.gtf.gz -O data/targets/${spec}_CLS3_targetDesign_mergedRegions.gtf.gz
	gunzip data/targets/${spec}_CLS3_targetDesign_mergedRegions.gtf.gz

       #download metadata
       wget https://github.com/guigolab/gencode-cls-master-table/releases/download/GencodeCLS_v1.0/${spec}_metadata.tsv.gz -O data/metadata/${spec}_metadata.tsv.gz
       zcat data/metadata/${spec}_metadata.tsv.gz | cut -f1 | cut -b1-7 | sort | uniq > data/metadata/${spec}_samples
       echo -e "downloaded the needed files $spec"
done

#----data generation
#source activate forPlotting
bash rcd_catalogDetectionStats_IClevel.sh

# After the job finishes, continue with the following
#----collate into a final stats file
for spec in Hv3 Mv2; do
	cat stats/${spec}.MergedTargetRegions.proportionDetected.* | grep -v sncRNA >> stats/${spec}.MergedTargetRegions.proportionDetected_AllTissues
done

Rscript rcd_catalogDetectionStats_ICs.R
