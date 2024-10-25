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

# masterTable
wget "https://github.com/guigolab/gencode-cls-master-table/releases/download/GencodeCLS_v3.0/Hv3_masterTable_refined.gtf.gz"

# targetDesign
wget "https://github.com/guigolab/gencode-cls-master-table/releases/download/Supplementary/Hv3_CLS3_targetDesign.gtf.gz"
