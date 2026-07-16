########################
###### Description
########################
# Download GENCODE - CLS Master Table data. Used to produce catalog specificity and intron chain analysis plots

# Handle errors
set -e          # exit on any non-0 exit status
set -o pipefail # exit on any non-0 exit status in pipe

########################
###### Download
########################

# masterTable - full
## Human
wget -O "Hv3_masterTable_refined.gtf.gz" "https://zenodo.org/records/15004659/files/Hv3_masterTable_refined.gtf.gz?download=1"
## Mouse
wget -O "Mv2_masterTable_refined.gtf.gz" "https://zenodo.org/records/15004659/files/Mv2_masterTable_refined.gtf.gz?download=1"

# masterTable - spliced
## Human
wget -O "Hv3_splicedmasterTable_refined.gtf.gz" "https://zenodo.org/records/15004659/files/Hv3_splicedmasterTable_refined.gtf.gz?download=1"
## Mouse
wget -O "Mv2_splicedmasterTable_refined.gtf.gz" "https://zenodo.org/records/15004659/files/Mv2_splicedmasterTable_refined.gtf.gz?download=1"
