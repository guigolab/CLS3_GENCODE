#rm stats/*.MergedTargetRegions.proportionDetected.Novel_perTissue_*

#!/usr/bin/bash
#$ -N targetRegionsperTissue
#$ -t 1-12
#$ -tc 12
#$ -P prj006070
#$ -q rg-el7,long-centos79,short-centos79
#$ -cwd
#$ -l virtual_free=32G,h_rt=01:00:00
#$ -pe smp 4
#$ -V
#$ -e targetRegionsperTlogs/$TASK_ID.err
#$ -o targetRegionsperTlogs/$TASK_ID.out

#cut -f1 $currnoHiss/tmp1.Hv3_sampleCoding | cut -b1-7 | sortU > Hv3_samples
#cut -f1 $currnoHiss/tmp1.Mv2_sampleCoding | cut -b1-7 | sortU > Mv2_samples

for spec in Hv3 Mv2;
do
	tissue=$(cat ${spec}_samples | sed -n ${SGE_TASK_ID}p)

#mkdir $currnoHiss/plots/TargetAnalaysis/stats
echo $currnoHiss
target=/users/rg/gkaur/LyRic_gencode_copy/targets/
out=$currnoHiss/plots/TargetAnalaysis/stats/

echo -e "spec\tclass\tregionType\tdetectedRegions\tallRegions\tdetectedTMs\tdetectedICs\tdetectedLoci\ttissue\tpropDetectedRegions" > $out/Hv3.refined_MergedTargetRegions.proportionDetected.Novel_perTissue
echo -e "spec\tclass\tregionType\tdetectedRegions\tallRegions\tdetectedTMs\tdetectedICs\tdetectedLoci\ttissue\tpropDetectedRegions" > $out/Mv2.refined_MergedTargetRegions.proportionDetected.Novel_perTissue

echo $tissue 

	for class in NONCODE miTranscriptome fantomCat gencodeLncRna refSeq bigTranscriptome CMfinderCRSs phyloCSF GWAScatalog UCE sncRNA fantomEnhancers VISTAenhancers allCatalogs uncaptured;
	do
			if [[ $spec == "Hv3" ]]; then
  	            		sp="hs"
				refv="v27"
        		else
                		sp="mm"
				refv="vM16"
        		fi
	
			if [[ $class == "sncRNA" ]]; then
				get="egrep \"miRNA|snoRNA|snRNA\""
			elif [[ $class == "allCatalogs" || $class == "uncaptured" ]]; then
				get="cat"
			else
				get="grep $class"
			fi
	
			if [[ $class == "uncaptured" ]]; then
				targets="grep 'target \"\"'"
			else
				targets="grep -v 'target \"\"'"
			fi
		
			cTs=$currnoHiss/IntronChainMT/${spec}_splicedmasterTable_refined.gtf
		        cTu=$currnoHiss/IntronChainMT/${spec}_unsplicedmasterTable_refined.gtf

		        noveltySpecify() {
		                awk '{if ($20 ~ /Intergenic/ || $20 ~ /Intronic/ || $20 ~ /Extends/ || $20 ~ /Antisense/ || $20 ~ /revIntronic/ || $20 ~ /runOn/) print $0}'
		        }

		currentDetected=`grep $tissue $currnoHiss/$spec\_masterTable_refined.gtf | grep -w transcript | eval "$targets" | awk -F"\"" '{print $6}' | grep . | perl -pi -e "s/,/\n/g" | eval "$get" | sort --parallel=4 -T /tmp/ | uniq | grep -Ff - $target/$sp\.allNonPcgTargetsMerged.targets.gtf | wc -l`
		currentDetectedICs=`grep $tissue <(cat ${cTs} ${cTu}) | grep -w transcript | eval "$targets" | eval "$get" | awk -F"\"" '{print $2}' | wc -l`
		currentDetectedNovelICs=`cat ${cTs} ${cTu} | eval "$get" | grep -w transcript | eval "$targets" | noveltySpecify | grep $tissue | wc -l`
		currentDetectedNovel_regions=`cat ${cTs} ${cTu} | noveltySpecify | grep $tissue | grep -w transcript | eval "$targets" | awk -F"\"" '{print $6}' | grep . | perl -pi -e "s/,/\n/g" | eval "$get" | grep -c -Ff - ${target}/${sp}.allNonPcgTargetsMerged.targets.gtf`

		#--------was a bit broad.. brought in from other tissues as well--## 
		#currentDetectedNovel_regions=`grep $tissue $currnoHiss/plots/${spec}_CLSlociToAnchsToSIDs${refv} | grep none | cut -f2 | perl -pi -e "s/,/\n/g" | grep -Ff - <(grep -w transcript $currnoHiss/${spec}_masterTable.gtf) | awk -F"\"" '{print $6}' | perl -pi -e "s/,/\n/g" | eval "$get" | sort --parallel=4 -T /tmp/ | uniq | grep -c -Ff - ${target}/${sp}.allNonPcgTargetsMerged.targets.gtf`
		echo $tissue $currentDetectedNovelICs $spec $class step1

		############
		
		currentDetectedNonNovel_regions=$(($currentDetected-$currentDetectedNovel_regions))
		currentDetectedNonNovelICs=$(($currentDetectedICs-$currentDetectedNovelICs))
		currentTotalRegions=`cat $target/$sp\.allNonPcgTargetsMerged.targets.gtf | eval "$get" | wc -l`
		echo step2

		echo -e "$spec\t$class\tdetectedNonNovel\t$currentDetectedNonNovel_regions\t$currentTotalRegions\t0\t$currentDetectedNonNovelICs\t0\t$tissue"|awk '{print $0"\t"($4/$5)*100}'>>$out/$spec\.MergedTargetRegions.proportionDetected.Novel_perTissue_${tissue}
		echo -e "$spec\t$class\tdetectedNovel\t$currentDetectedNovel_regions\t$currentTotalRegions\t0\t$currentDetectedNovelICs\t0\t$tissue"|awk '{print $0"\t"($4/$5)*100}'>>$out/$spec\.MergedTargetRegions.proportionDetected.Novel_perTissue_${tissue}
		echo -e "$spec\t$class\tdetectedALL\t$currentDetected\t$currentTotalRegions\t0\t$currentDetectedICs\t0\t$tissue"|awk '{print $0"\t"($4/$5)*100}'>>$out/$spec\.MergedTargetRegions.proportionDetected.Novel_perTissue_${tissue}
	done
done
