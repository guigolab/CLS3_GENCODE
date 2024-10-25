########################
###### Description
########################
# Additional files. Too big to upload. Used to produce:
# 1) Figure 3B (polyA signals)

########################
###### Download and process
########################

# Handle errors
set -e          # exit on any non-0 exit status
set -o pipefail # exit on any non-0 exit status in pipe

############ hg38.polyAsignals.bed ############
# Download:
wget "https://public-docs.crg.es/rguigo/Papers/2017_lagarde-uszczynska_CLS/data/polyA/signals/hg38.polyAsignals.bed.gz"

# Extract:
gzip -d hg38.polyAsignals.bed.gz

############ recount3_final.pass1.gff3.gz ############
# Download:


# Extract:
