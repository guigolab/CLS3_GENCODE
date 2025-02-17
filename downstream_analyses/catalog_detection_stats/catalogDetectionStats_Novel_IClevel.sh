#!/usr/bin/bash
#$ -N targetRegions
#$ -t 1-15
#$ -tc 15
#$ -P prj006070
#$ -q rg-el7,long-centos79,short-centos79
#$ -cwd
#$ -l virtual_free=32G,h_rt=01:00:00
#$ -pe smp 4
#$ -V
#$ -e targetRegionslogs/$TASK_ID.err
#$ -o targetRegionslogs/$TASK_ID.out

class=$(echo -e "NONCODE\nmiTranscriptome\nfantomCat\ngencodeLncRna\nrefSeq\nbigTranscriptome\nCMfinderCRSs\nphyloCSF\nGWAScatalog\nUCE\nsncRNA\nfantomEnhancers\nVISTAenhancers\nallCatalogs\nuncaptured" | sed -n ${SGE_TASK_ID}p)

mkdir $currnoHiss/plots/TargetAnalaysis/stats/
echo $currnoHiss
target=/users/rg/gkaur/LyRic_gencode_copy/targets/
out=$currnoHiss/plots/TargetAnalaysis/stats/

#############

echo -e "spec\tclass\tregionType\tdetectedRegions\tallRegions\tdetectedTMs\tdetectedICs\tdetectedLoci\ttissue\tpropDetectedRegions" > $out/Hv3.refined_MergedTargetRegions.proportionDetected.Novel_AllTissues
echo -e "spec\tclass\tregionType\tdetectedRegions\tallRegions\tdetectedTMs\tdetectedICs\tdetectedLoci\ttissue\tpropDetectedRegions" > $out/Mv2.refined_MergedTargetRegions.proportionDetected.Novel_AllTissues

	for spec in Hv3 Mv2; do
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

	currentDetected=`cat $currnoHiss/$spec\_masterTable_refined.gtf | grep -w transcript | eval "$targets" | awk -F"\"" '{print $6}' | grep . | perl -pi -e "s/,/\n/g" | eval "$get" | sort --parallel=4 -T /tmp/ | uniq | grep -Ff - $target/$sp\.allNonPcgTargetsMerged.targets.gtf | wc -l`
	currentDetectedICs=`cat ${cTs} ${cTu} | grep -w transcript | eval "$targets" | eval "$get" | awk -F"\"" '{print $2}' | wc -l`

	#-----Novel include intergenic, intronic, revIntronic, antisense, runOn, extends
	currentDetectedNovelICs=`cat ${cTs} ${cTu} | eval "$get" | grep -w transcript | eval "$targets" | noveltySpecify | wc -l`
	currentDetectedNovel_regions=`cat ${cTs} ${cTu} | noveltySpecify | eval "$targets" | awk -F"\"" '{print $6}' | grep . | perl -pi -e "s/,/\n/g" | eval "$get" | sort --parallel=4 -T /tmp/ | uniq | grep -c -Ff - ${target}/${sp}.allNonPcgTargetsMerged.targets.gtf`
	
	echo step1	

##############
	
	currentDetectedNonNovel_regions=$(($currentDetected-$currentDetectedNovel_regions))
	currentDetectedNonNovelICs=$(($currentDetectedICs-$currentDetectedNovelICs))

	currentTotalRegions=`cat $target/$sp\.allNonPcgTargetsMerged.targets.gtf | eval "$get" | wc -l`
	echo step2

	echo -e "$spec\t$class\tdetectedNonNovel\t$currentDetectedNonNovel_regions\t$currentTotalRegions\t0\t$currentDetectedNonNovelICs\t0\tallTissues"|awk '{print $0"\t"($4/$5)*100}'>>$out/$spec\.MergedTargetRegions.proportionDetected.Novel.${class}
	echo -e "$spec\t$class\tdetectedNovel\t$currentDetectedNovel_regions\t$currentTotalRegions\t0\t$currentDetectedNovelICs\t0\tallTissues"|awk '{print $0"\t"($4/$5)*100}'>>$out/$spec\.MergedTargetRegions.proportionDetected.Novel.${class}  #log-> "\t"log($6)/log(2)
	echo -e "$spec\t$class\tdetectedALL\t$currentDetected\t$currentTotalRegions\t0\t$currentDetectedICs\t0\tallTissues"|awk '{print $0"\t"($4/$5)*100}'>>$out/$spec\.MergedTargetRegions.proportionDetected.Novel.${class}

	done
echo "XXdoneXX" >&2
