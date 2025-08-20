########################
###### Description
########################
# Download GENCODE - CLS Master Table data. Used to produce Figure S17A & S17C

# Handle errors
set -e          # exit on any non-0 exit status
set -o pipefail # exit on any non-0 exit status in pipe

########################
###### Download
########################

# masterTable - full
## Human
wget "https://github.com/guigolab/gencode-cls-master-table/releases/download/GencodeCLS_v3.0/Hv3_masterTable_refined.gtf.gz"
## Mouse
wget "https://github.com/guigolab/gencode-cls-master-table/releases/download/GencodeCLS_v3.0/Mv2_masterTable_refined.gtf.gz"

# masterTable - spliced
## Human
wget "https://github.com/guigolab/gencode-cls-master-table/releases/download/GencodeCLS_v3.0/Hv3_splicedmasterTable_refined.gtf.gz"
## Mouse
wget "https://github.com/guigolab/gencode-cls-master-table/releases/download/GencodeCLS_v3.0/Mv2_splicedmasterTable_refined.gtf.gz"

# samplesMetadata
## Human
wget "https://github.com/guigolab/gencode-cls-master-table/releases/download/GencodeCLS_v1.0/Hv3_metadata.tsv.gz"
## Mouse
wget "https://github.com/guigolab/gencode-cls-master-table/releases/download/GencodeCLS_v1.0/Mv2_metadata.tsv.gz"

# targetDesign
## Human
wget "https://github.com/guigolab/gencode-cls-master-table/releases/download/Supplementary/Hv3_CLS3_targetDesign.gtf.gz"
## Mouse
wget "https://github.com/guigolab/gencode-cls-master-table/releases/download/Supplementary/Mv2_CLS3_targetDesign.gtf.gz"

