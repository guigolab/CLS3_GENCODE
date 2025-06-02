#!/bin/bash
set -e
set -o pipefail

#----remove old data if any
rm -rf data/ plots/ stats/
#----download required files
mkdir -p data/ data/masterTable/ data/recount/ stats/ stats/filtered_ICs/ stats/filtered_ICtables/ plots/
echo -e "#### created needed directories ####"

for spec in Hv3 Mv2; do
        #----download CLS transcript (IC) files
        wget https://github.com/guigolab/gencode-cls-master-table/releases/download/GencodeCLS_v3.0/${spec}_splicedmasterTable_refined.gtf.gz -O data/masterTable/${spec}_splicedmasterTable_refined.gtf.gz 
        gunzip data/masterTable/${spec}_splicedmasterTable_refined.gtf.gz
	echo "downloaded IC table"
	#----download non-redundant spliced transcripts
	wget https://github.com/guigolab/cls3-final-files-description/releases/download/data.files.v1/Trmap${spec}_Pairs_Spliced.txt.gz -O data/Trmap${spec}_Pairs_Spliced.txt.gz
	gunzip data/Trmap${spec}_Pairs_Spliced.txt.gz
	echo "downloaded Trmap pairs"
	#--download recount files
	if [[ $spec == "Hv3" ]]; then
        	wget https://github.com/guigolab/cls3-final-files-description/releases/download/data.files.v1/recount3_final.pass1.gff3.gz -O data/recount/recount3_Hv3.gff3.gz
                gunzip data/recount/recount3_Hv3.gff3.gz
	else
                #recountFile=$currnoHiss/recount_support/m39_recount3_final.pass1_formatted.gff3
                #--recount file for liftOver that I did from m39 to m38--# did not have the original from refSeq earlier.
                #recountFile=$currnoHiss/recount_support/mm10_recount3_final.pass1_bedcoords_highestScore.gff3   #mm39 coords converted to mm10; highest score for overlapping coords selected
                #---original mouse recount data is in mm10/m38.. use that instead to be more precise 
         	wget https://github.com/guigolab/cls3-final-files-description/releases/download/data.files.v1/m38_recount3.gff3.gz -O data/recount/recount3_Mv2.gff3.gz
                gunzip data/recount/recount3_Mv2.gff3.gz
	fi
		echo "downloaded all"
done
##-----recount support; supported models curve
for filter in {1..500}; 
do	#creates recount support intron files for 1 to 500 read support filter. for the support curve
	for spec in Hv3 Mv2;
	do
		recountFile=data/recount/recount3_${spec}.gff3
		echo "downloaded recount ${spec}"
#-------get filtered recount junctions-----#
		awk -v f="$filter" -F"\t" '{if($6 > f) print $0}' ${recountFile} > stats/filtered_ICs/${spec}.recount3_final_filteredMT${filter}reads.pass1.gff3
		echo "recount junctions filtered $filter reads"
#-------tag the recount support in the IC table based on the filters (# of read support in recount)
		perl getIntrons_tagRecountSupport.pl data/masterTable/${spec}_splicedmasterTable_refined.gtf data/Trmap${spec}_Pairs_Spliced.txt stats/filtered_ICs/${spec}.recount3_final_filteredMT${filter}reads.pass1.gff3 > stats/filtered_ICtables/${spec}_splicedmasterTable_introns_recountICtable_filterMT${filter}reads.gtf
	
#-------calculate the stats for each filter for the cumulative curve-----#
		currentFile=stats/filtered_ICtables/${spec}_splicedmasterTable_introns_recountICtable_filterMT${filter}reads.gtf
	
#----get stats for all/novel/intergenic recount3 support--------------------------------################
		for type in all novel intergenic;
		do	if [[ $type == "all" ]]; then
				get(){
				cat
				}
			elif [[ $type == *novel* ]]; then
				get(){
				egrep "refCompare \"Intergenic\"|refCompare \"Extends\"|refCompare \"Intronic\"|refCompare \"Antisense\""
				}
			elif [[ $type == *intergenic* ]]; then 
				get(){
				grep "refCompare \"Intergenic\""
				}
			fi

		Nsupported=`cat $currentFile | grep -w transcript | get |grep -c "transcriptRecountSupport \"supported\""`
#		pacBioS=`cat $currentFile | grep -w transcript | get |grep "transcriptRecountSupport \"supported\"" | grep -c "tech \"pacBioOnly\""`
#		ontS=`cat $currentFile | grep -w transcript | get |grep "transcriptRecountSupport \"supported\"" |grep -c "tech \"ontOnly\""`
#		bothTechS=`cat $currentFile | grep -w transcript | get |grep "transcriptRecountSupport \"supported\"" |grep -c "tech \"pacBio+ont\""`
#	echo "$type ##### $Nsupported\t$pacBioS\t$ontS\t$bothTechS"
		Nunsupported=`cat $currentFile | grep -w transcript | get |grep -c "transcriptRecountSupport \"unsupported\""`
#		pacBioU=`cat $currentFile | grep -w transcript | get |grep "transcriptRecountSupport \"unsupported\"" | grep -c "tech \"pacBioOnly\""`
#		ontU=`cat $currentFile | grep -w transcript | get |grep "transcriptRecountSupport \"unsupported\"" |grep -c "tech \"ontOnly\""`
#		bothTechU=`cat $currentFile | grep -w transcript | get |grep "transcriptRecountSupport \"unsupported\"" |grep -c "tech \"pacBio+ont\""`

		(( Total=$Nsupported+$Nunsupported))
#		(( pacBioT=$pacBioU+$pacBioS ))
#		(( ontT=$ontU+$ontS ))
#		(( bothTechT=$bothTechU+$bothTechS ))
	#echo $a $b|awk '{print $1/$2*100"%"}'
		allSperc=`echo $Nsupported $Total | awk '{print $1/$2*100}'` # | awk '{print $1/$2*100"%"}'` #$(( Nsupported/Total ))`
#		pacBioSperc=`echo $pacBioS $pacBioT | awk '{print $1/$2*100}'`
#		ontSperc=`echo $ontS $ontT | awk '{print $1/$2*100}'` #$ontS/$ontT 
#		bothTechSperc=`echo $bothTechS $bothTechT | awk '{print $1/$2*100}'` #$bothTechU/$bothTechT 
	#	echo -e "spec\tfilterThreshold\tsupportedICs\ttech\tsupport\n" > stats/filterThreshold_vs_supportedTMs.txt
		echo -e "${spec}\t$filter\t$Nsupported\t$allSperc\t${type}_ICs\tsupported >> stats/${spec}_filterThreshold_${filter}_supportedTMs_3classes"
#${spec}\t$filter\t$pacBioS\t$pacBioSperc\t${type}_pacBioOnly_ICs\tsupported
#${spec}\t$filter\t$ontS\t$ontSperc\t${type}_ONTonly_ICs\tsupported
#${spec}\t$filter\t$bothTechS\t$bothTechSperc\t${type}_pacBio+ONT_ICs\tsupported" >> stats/${spec}_filterThreshold_${filter}_supportedTMs_3classes
		echo "got stats for $filter reads filter"
		done
	done
#echo -e "$NsupportedNovel\t$pacBioSNovel\t$ontSNovel\t$bothTechSNovel" > stats/filterThreshold_${filter}_NovelsupportedTMs
done


#=======
###echo -e "spec\tfilterThreshold\tsupportedICs\tpercsupportedICs\ttech\tsupport" > stats/${spec}_filterThreshold_vs_supportedTMs_3classes.txt
###cat stats/${spec}_filterThreshold_*_supportedTMs_3classes >> stats/${spec}_filterThreshold_vs_supportedTMs_3classes.txt
###rm stats/${spec}_filterThreshold_*_supportedTMs_3classes
