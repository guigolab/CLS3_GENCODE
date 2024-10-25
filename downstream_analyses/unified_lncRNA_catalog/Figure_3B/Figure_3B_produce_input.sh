#!/usr/bin/env bash

# +++++++++++++++++++++++++++++++++++++++++++++
# REQUIREMENTS
# +++++++++++++++++++++++++++++++++++++++++++++
# 1) Catalogs in gtf.gz file format in ./Data/Source/Catalogs/ directory.
# 2) Additional files in ./Data/Source/Additional/ directory.
# 3) Installed bedtools v2.31.0
# 4) Utils in ./Utils/ directory:
#     - bed12togff
#     - buildLoci.pl
#     - gffToHash.pm
#     - hashToGff.pm
#     - extractTranscriptEndsFromBed12.pl
#     - gff2bed_full.pl
#     - join.py
#     - sortbed

# +++++++++++++++++++++++++++++++++++++++++++++
# DESCRIPTION
# +++++++++++++++++++++++++++++++++++++++++++++
# Takes data in gtf format from /Source/Raw_gtfs and calculates Figure 3B plot input


##########################################################################################
# Clean-up generated data
##########################################################################################
echo '++++++ Cleaning up data from previous run ++++++'
rm -r ./Data/.temp/Figure_3B/
echo '>> Done <<'

##########################################################################################
# Create directory structure
##########################################################################################
echo '++++++ Creating necessary directory structure ++++++'
# Dir for .temp data
mkdir -p ./Data/.temp/Figure_3B/
mkdir ./Data/.temp/Figure_3B/Cage/
mkdir ./Data/.temp/Figure_3B/PolyA/
mkdir ./Data/.temp/Figure_3B/Fl/
mkdir ./Data/.temp/Figure_3B/Decompressed/
mkdir ./Data/.temp/Figure_3B/Recount_support/
echo '>> Done <<'

# Handle errors
set -e          # exit on any non-0 exit status
set -o pipefail # exit on any non-0 exit status in pipe

##########################################################################################
##########################################################################################
# Generate data for FIGURE 3B
##########################################################################################
##########################################################################################

##########################################################################################
# % Support
##########################################################################################

echo 'Get fl-tx'

while read file || [ -n "$file" ]
    do
        echo $file
        zcat ./Data/Source/Catalogs/$file.gtf.gz | ./Utils/gff2bed_full.pl - | awk 'BEGIN{FS=OFS="\t"}$6!="."' | awk '$1 ~ /^chr[0-9XYM]{1,2}$/ {print $0}' | awk '$10 > 1 {print $0}' | gzip > ./Data/.temp/Figure_3B/$file.bed12.gz

        while read end dist
        do
            echo $end
            zcat ./Data/.temp/Figure_3B/$file.bed12.gz | ./Utils/extractTranscriptEndsFromBed12.pl $end | ./Utils/sortbed | gzip > ./Data/.temp/Figure_3B/$file.$end.bed.gz
        done < ./Data/Configs/Figure_3B_ends.config

        ## polyA
        echo '++++++ Looking for polyA supported TM ++++++'
        while read end dist
        do
            zcat ./Data/.temp/Figure_3B/$file.$end.bed.gz | ./Utils/sortbed | bedtools slop -s -l 50 -r 10 -i stdin -g ./Data/Source/Additional/chromInfo_hg38.txt | bedtools intersect -u -s -a stdin -b ./Data/Source/Additional/hg38.polyAsignals.bed | gzip > ./Data/.temp/Figure_3B/$file.$end.bed.vspolyAsignals.bedtsv.gz
        done < ./Data/Configs/Figure_3B_ends.config
        echo '>> Done <<'

        ## CAGE FANTOM phase 1+2
        echo '++++++ Looking for CAGE supported TM ++++++'
        while read end dist
        do
            zcat ./Data/.temp/Figure_3B/$file.$end.bed.gz | ./Utils/sortbed | bedtools slop -s -l 50 -r 50 -i stdin -g ./Data/Source/Additional/chromInfo_hg38.txt | bedtools intersect -u -s -a stdin -b ./Data/Source/Additional/hg38_fair+new_CAGE_peaks_phase1and2.bed | gzip > ./Data/.temp/Figure_3B/$file.$end.bed.vsCage.fantom.bedtsv.gz
        done < ./Data/Configs/Figure_3B_ends.config
        echo '>> Done <<'

        ## Extract Cage and polyA supported transcripts
        echo '++++++ Extracting CAGE and polyA supported transcripts ++++++'
        zcat ./Data/.temp/Figure_3B/$file.5.bed.vsCage.fantom.bedtsv.gz | cut -f4 | sed 's/,/\n/g' | sort | uniq > ./Data/.temp/Figure_3B/$file.tx.cage
        zcat ./Data/.temp/Figure_3B/$file.3.bed.vspolyAsignals.bedtsv.gz | cut -f4 | sed 's/,/\n/g' | sort | uniq > ./Data/.temp/Figure_3B/$file.tx.polyA
        echo '>> Done <<'

        ## Get fl-tx for Unified
        ./Utils/join.py -a ./Data/.temp/Figure_3B/$file.tx.cage -b ./Data/.temp/Figure_3B/$file.tx.polyA -x 1 -y 1 > ./Data/.temp/Figure_3B/Fl/$file.fl.tx.tsv

        echo '++++++ Aquiring geneIds and txIds ++++++'
        # --------------------------------------------------------------------
        # GTF from hg38.bed12
        # --------------------------------------------------------------------
        zcat ./Data/.temp/Figure_3B/$file.bed12.gz | ./Utils/bed12togff - > ./Data/.temp/Figure_3B/$file.lncRNAs.hg38.gtf

        # --------------------------------------------------------------------
        # build the gene ids
        # --------------------------------------------------------------------
        echo $full_name
        bedtools intersect -s -wao -a ./Data/.temp/Figure_3B/$file.lncRNAs.hg38.gtf -b ./Data/.temp/Figure_3B/$file.lncRNAs.hg38.gtf | ./Utils/buildLoci.pl - > ./Data/.temp/Figure_3B/$file.lncRNAs.hg38.gene.tmp.gtf
        ./Utils/gff2gff.pl ./Data/.temp/Figure_3B/$file.lncRNAs.hg38.gene.tmp.gtf > ./Data/.temp/Figure_3B/$file.lncRNAs.hg38.gene.gtf

        # --------------------------------------------------------------------
        # get the gene ids from the hg38.gtfs
        # --------------------------------------------------------------------
        cat ./Data/.temp/Figure_3B/$file.lncRNAs.hg38.gene.gtf | awk '{print $10}' | sed 's/;//g' | sed 's/"//g' | sort | uniq > ./Data/.temp/Figure_3B/$file.buildLoci.geneIds

        # --------------------------------------------------------------------
        # Get the txs Ids
        # --------------------------------------------------------------------
        cat ./Data/.temp/Figure_3B/$file.lncRNAs.hg38.gene.gtf | awk '$3=="exon"{print $12}' | sed 's/;//g'| sed 's/"//g'| sort | uniq > ./Data/.temp/Figure_3B/$file.buildLoci.txIds

        echo '>> Done <<'
        ###################################
        echo 'Generating plot input data'
        # --------------------------------------------------------------------
        fl=`cat ./Data/.temp/Figure_3B/Fl/$file.fl.tx.tsv | wc -l`
        total=`zcat ./Data/.temp/Figure_3B/$file.bed12.gz | cut -f4 | sort | uniq | wc -l`
        propFl=`echo $fl / $total | bc -l`
        gene=`cat ./Data/.temp/Figure_3B/$file.buildLoci.geneIds | sort | uniq | wc -l`
        tx=`cat ./Data/.temp/Figure_3B/$file.buildLoci.txIds | sort | uniq | wc -l`
        nbTx=`echo $tx / $gene | bc -l`
        echo -e "$file\t$propFl\t$gene\t$nbTx" >> ./Figure_3B/Figure_3B_input.tsv

    done < ./Data/Configs/Figure_3B_catalogs.config


##########################################################################################
# Recount Support
##########################################################################################
# Copy necessary files and transform them inplace
## Recount data
zcat ./Data/Source/Additional/recount3_final.pass1.gff3.gz | awk '{print $1"\t"$4"\t"$5"\t"$6"\t"$7}' > ./Data/.temp/Figure_3B/Recount_support/recount3_final.pass1.pseudo.bed

# Run interecting: catalogs with recount data
echo "Running intersection"
# Create results.tsv with a header
echo -e "catalogue\ttotal\tsupported\tunsupported\tsupported_percentage" > ./Data/.temp/Figure_3B/Recount_support/results.tsv
while read file || [ -n "$file" ]
  do
      echo "Name: $file"
      # Convert to gtf for compatibility
      zcat ./Data/.temp/Figure_3B/$file.bed12.gz | ./Utils/bed12togff - > ./Data/.temp/Figure_3B/Recount_support/$file.gtf

      echo "~~~ Preparing introns ~~~"
      # For introns to be compatible with recount data, they need to be modified: start+1 stop-1
      # I included this in modified version on make_introns.awk script
      sort -k12,12 -k4,4n -k5,5n ./Data/.temp/Figure_3B/Recount_support/$file.gtf | awk -v fldgn=10 -v fldtr=12 -f ./Utils/make_introns_for_recount.awk | awk '{print $1"\t"$4"\t"$5"\t"$7"\t"$12}' | sed 's/;//g' | sed 's/"//g' > ./Data/.temp/Figure_3B/Recount_support/$file.introns
      echo "+++ Done +++"
      echo "~~~ Right-outer joining recount data with ./Data/.temp/Figure_3B/Recount_support/$file introns  ~~~"
      # Right-outer join seems easier for join.py script than left-outer join which is I guess not supported
      # Results are as expected, but be sure that u submit recount as "file a"
      ./Utils/join.py -a ./Data/.temp/Figure_3B/Recount_support/recount3_final.pass1.pseudo.bed -b ./Data/.temp/Figure_3B/Recount_support/$file.introns -x 1,2,3,5 -y 1,2,3,4 -u > ./Data/.temp/Figure_3B/Recount_support/${file}_with_recount.tsv
      # Calculate supported and unsupported tx
      total_tx=`cat ./Data/.temp/Figure_3B/Recount_support/$file.introns | awk '{print $5}' | sort | uniq | wc -l`
      # Fortunately as < 50 awk prints also "-" that were assigned to unmatched introns
      not_supported_tx=`cat ./Data/.temp/Figure_3B/Recount_support/${file}_with_recount.tsv | awk '$6 < 50 {print $5}' | sort | uniq | wc -l`
      supported_tx=$((total_tx-not_supported_tx))
      supported_tx_percentage=`echo "scale=2 ; $supported_tx / $total_tx" | bc`
      # Save results to tsv
      echo -e "$file\t$total_tx\t$supported_tx\t$not_supported_tx\t$supported_tx_percentage" >> ./Data/.temp/Figure_3B/Recount_support/results.tsv
      echo "+++ Done +++"
  done < ./Data/Configs/Figure_3B_catalogs.config

# Copy results to /Data/Processed/Recount
cat ./Data/.temp/Figure_3B/Recount_support/results.tsv > ./Figure_3B/Figure_3B_recount.tsv

# Return
echo '################## COMPLETED ##################'
