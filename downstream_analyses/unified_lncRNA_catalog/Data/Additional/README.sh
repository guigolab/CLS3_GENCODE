#############################
###### Description
#############################
# Additional files. Too big to upload.

#############################
###### Download and process
#############################

# Handle errors
set -e          # exit on any non-0 exit status
set -o pipefail # exit on any non-0 exit status in pipe

############ hg38.polyAsignals.bed ############
# Download:
wget "https://public-docs.crg.es/rguigo/Papers/2017_lagarde-uszczynska_CLS/data/polyA/signals/hg38.polyAsignals.bed.gz"

# Extract:
gzip -d hg38.polyAsignals.bed.gz

############ mm10.polyAsignals.bed ############
# Download:
wget "https://public-docs.crg.es/rguigo/Papers/2017_lagarde-uszczynska_CLS/data/polyA/signals/mm10.polyAsignals.bed.gz"

# Extract:
gzip -d mm10.polyAsignals.bed.gz

############ hg38 CAGE peaks ############
wget "https://github.com/guigolab/cls3-final-files-description/releases/download/data.files.v1/hg38_fair+new_CAGE_peaks_phase1and2.bed"

############ mm10 CAGE peaks ############


############ Human recount3 data (recount3_final.pass1.gff3.gz) ############
# Download:
wget "https://github.com/guigolab/cls3-final-files-description/releases/download/data.files.v1/recount3_final.pass1.gff3.gz"

############ Mouse recount3 data (m38_recount3.bed.gz) ############
# Download 
wget "https://zenodo.org/records/13946596/files/m38_recount3.bed.gz?download=1" -O m38_recount3.bed.gz

############ Gencode v47 ############
wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.chr_patch_hapl_scaff.annotation.gtf.gz"

############ Gencode vM36 ############
wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.chr_patch_hapl_scaff.annotation.gtf.gz"

############ Human chromosome sizes ############
wget "https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes"

############ Mouse chromosome sizes ############
wget "https://hgdownload.soe.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes"
