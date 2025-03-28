#!/usr/bin/bash

set -e
set -o pipefail

working=data/bam/

echo -e "capDesign\tpre_post\ttissue\ttotalMappedReads\ttotalMappedReadsnoSIRV\ttotalMappedReadsnoERCC\ttotalMappedReadsnoSIRVnoERCC\ttotalonlyERCCReads\tonTargetNonERCCReads\tonTargetwithERCCReads\tonTargetonlyERCCReads\tOTRnonERCC\tOTRwithERCC\tOTRonlyERCC" > stats/header_OTR

for file in $working/pacBio*bam $working/ont*bam;
do
  dataset=`basename $file`
  echo $dataset

  if [[ $dataset == *"Mv2"* || $dataset == *"Mpre"* ]]; then
  	targetsFile="/users/rg/gkaur/targetFiles/mm.allNonPcgTargetsMerged.targets.gtf"
  else
  	targetsFile="/users/rg/gkaur/targetFiles/hs.allNonPcgTargetsMerged.targets.gtf"
  fi

nonT_ERCC='/users/project/gencode_006070/gkaur/ERCCfiles'
capDesign=`echo $dataset | awk -F "_" '{print $1}'`
pre_post=`echo $dataset | awk -F "_" '{print $2}'`
tissue=`echo $dataset | awk -F "_" '{print $4}'`

totalMappedReads=`samtools view $file |cut -f1|sort --parallel=4|uniq|wc -l`
totalMappedReadsnoSIRV=`samtools view $file | cut -f1,3 | fgrep -v SIRV | cut -f1 | sort --parallel=4 | uniq | wc -l`
totalMappedReadsnoERCC=`samtools view $file | cut -f1,3 | fgrep -v ERCC | cut -f1 | sort --parallel=4 | uniq | wc -l`
totalMappedReadsnoSIRVnoERCC=`samtools view $file | cut -f1,3 | fgrep -v SIRV | fgrep -v ERCC | cut -f1 | sort --parallel=4 | uniq | wc -l`
totalonlyERCCReads=`samtools view $file | cut -f1,3 | fgrep ERCC | cut -f1 | sort --parallel=4 | uniq | wc -l`

ml load BEDTools/2.29.2-GCC-9.3.0;
####
cat $targetsFile | bedtools intersect -split -u -bed -abam $file -b stdin |cut -f 1-6 > `basename $file .bam`.vs.`basename $targetsFile .gtf`.bed

onTargetNonERCCReads=`cat \`basename $file .bam\`.vs.\`basename $targetsFile .gtf\`.bed | fgrep -v ERCC| cut -f4 |perl -ne '$_=~s/\/\d{1}$//; print' |sort --parallel=4|uniq|wc -l`
	#let onTargetNonERCCReadsRate=$onTargetNonERCCReads/$totalMappedReadsnoSIRVnoERCC

onTargetwithERCCReads=`cat \`basename $file .bam\`.vs.\`basename $targetsFile .gtf\`.bed | cut -f4 |perl -ne '$_=~s/\/\d{1}$//; print' |sort --parallel=4|uniq|wc -l`
	#let onTargetwithERCCReadsRate=$onTargetwithERCCReads/$totalMappedReadsnoSIRV

onTargetonlyERCCReads=`cat \`basename $file .bam\`.vs.\`basename $targetsFile .gtf\`.bed | fgrep ERCC| cut -f4 |perl -ne '$_=~s/\/\d{1}$//; print' |sort --parallel=4|uniq|wc -l`
	#let onTargetonlyERCCReadsRate=$onTargetonlyERCCReads/$totalonlyERCCReads


	#cat $nonT_ERCC | bedtools intersect -split -u -bed -abam $file -b stdin |cut -f 1-6 > `basename $file .bam`.vs.`basename $nonT_ERCC .gtf`.bed
	#nT_ERCC='ERCCfiles/nonTarget_ERCC'

	#offTargetERCCReads=`cat \`basename $file .bam\`.vs.\`basename $nonT_ERCC .gtf\`.bed | fgrep ERCC| cut -f4 |sort --parallel=4|uniq|wc -l`
	#let offTargetReads=$totalMappedReads_withERCC-$onTargetNonERCCReads-$onTargetERCCReads
	#let offTargetnonERCCReads=$totalMappedReads-$onTargetNonERCCReads
	#let onTargetwithERCC=$onTargetERCCReads+$onTargetNonERCCReads

echo -e "$capDesign\t$pre_post\t$tissue\t$totalMappedReads\t$totalMappedReadsnoSIRV\t$totalMappedReadsnoERCC\t$totalMappedReadsnoSIRVnoERCC\t$totalonlyERCCReads\t$onTargetNonERCCReads\t$onTargetwithERCCReads\t$onTargetonlyERCCReads"|awk '{print $0"\t"$9/$7"\t"$10/$5"\t"$11/$8}' > stats/${dataset}_OTR
done
