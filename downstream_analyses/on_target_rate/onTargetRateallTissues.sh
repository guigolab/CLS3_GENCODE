#!/usr/bin/bash
#$ -N OTRcalc
#$ -t 1-8
#$ -tc 8
#$ -q rg-el7,long-centos79
#$ -cwd
#$ -l virtual_free=32G,h_rt=12:00:00
#$ -pe smp 4
#$ -P prj006070
#$ -V
#$ -e errOTRcalc2
#$ -o outOTRcalc2

working=/users/project/gencode_006070/gkaur/CapTrap-CLSfigs/OTRresults
echo -e "capDesign\tpre_post\ttissue\ttotalMappedReads\ttotalMappedReadsnoSIRV\ttotalMappedReadsnoERCC\ttotalMappedReadsnoSIRVnoERCC\ttotalonlyERCCReads\tonTargetNonERCCReads\tonTargetwithERCCReads\tonTargetonlyERCCReads\tOTRnonERCC\tOTRwithERCC\tOTRonlyERCC" > header_OTR

file=$(ls -1 ${working}/*allT*bam | sed -n ${SGE_TASK_ID}p)

dataset=`basename $file`
echo $dataset

if [[ $dataset == *"Mv2"* || $dataset == *"Mpre"* ]]; then
        targetsFile="/users/rg/gkaur/targetFiles/mm.allNonPcgTargetsMerged.targets.gtf"
else
        targetsFile="/users/rg/gkaur/targetFiles/hs.allNonPcgTargetsMerged.targets.gtf"
fi

nonT_ERCC='/users/project/gencode_006070/gkaur/ERCCfiles'
#ALLregions="hs.gencode24+SIRVS.gff"
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
onTargetwithERCCReads=`cat \`basename $file .bam\`.vs.\`basename $targetsFile .gtf\`.bed | cut -f4 |perl -ne '$_=~s/\/\d{1}$//; print' |sort --parallel=4|uniq|wc -l`
onTargetonlyERCCReads=`cat \`basename $file .bam\`.vs.\`basename $targetsFile .gtf\`.bed | fgrep ERCC| cut -f4 |perl -ne '$_=~s/\/\d{1}$//; print' |sort --parallel=4|uniq|wc -l`

echo -e "$capDesign\t$pre_post\t$tissue\t$totalMappedReads\t$totalMappedReadsnoSIRV\t$totalMappedReadsnoERCC\t$totalMappedReadsnoSIRVnoERCC\t$totalonlyERCCReads\t$onTargetNonERCCReads\t$onTargetwithERCCReads\t$onTargetonlyERCCReads"|awk '{print $0"\t"$9/$7"\t"$10/$5"\t"$11/$8}' > ${dataset}_OTR                
