#-----header
for spec in Hv3 Mv2; do
        echo -e "spec\tclass\tregionType\tdetectedRegions\tallRegions\tdetectedTMs\tdetectedICs\tdetectedLoci\ttissue\tpropDetectedRegions" > stats/${spec}.MergedTargetRegions.proportionDetected_AllTissues

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
	cTs=data/masterTable/${spec}_splicedmasterTable_refined.gtf
        cTu=data/masterTable/${spec}_unsplicedmasterTable_refined.gtf
	
	noveltySpecify() {
        	awk '{if ($20 ~ /Intergenic/ || $20 ~ /Intronic/ || $20 ~ /Extends/ || $20 ~ /Antisense/ || $20 ~ /revIntronic/ || $20 ~ /runOn/) print $0}'
        }

	currentDetected=`cat data/masterTable/${spec}_masterTable_refined.gtf | grep -w transcript | eval "$targets" | awk -F"\"" '{print $6}' | grep . | perl -pi -e "s/,/\n/g" | eval "$get" | sort --parallel=4 -T /tmp/ | uniq | grep -Ff - data/targets/${spec}_CLS3_targetDesign_mergedRegions.gtf | wc -l`
	currentDetectedICs=`cat ${cTs} ${cTu} | grep -w transcript | eval "$targets" | eval "$get" | awk -F"\"" '{print $2}' | wc -l`

	#-----Novel include intergenic, intronic, revIntronic, antisense, runOn, extends
	currentDetectedNovelICs=`cat ${cTs} ${cTu} | eval "$get" | grep -w transcript | eval "$targets" | noveltySpecify | wc -l`
	currentDetectedNovel_regions=`cat ${cTs} ${cTu} | noveltySpecify | eval "$targets" | awk -F"\"" '{print $6}' | grep . | perl -pi -e "s/,/\n/g" | eval "$get" | sort --parallel=4 -T /tmp/ | uniq | grep -c -Ff - data/targets/${spec}_CLS3_targetDesign_mergedRegions.gtf`
	
	echo step1	

##############
	
	currentDetectedNonNovel_regions=$(($currentDetected-$currentDetectedNovel_regions))
	currentDetectedNonNovelICs=$(($currentDetectedICs-$currentDetectedNovelICs))

	currentTotalRegions=`cat data/targets/${spec}_CLS3_targetDesign_mergedRegions.gtf | eval "$get" | wc -l`
	echo step2

	echo -e "$spec\t$class\tdetectedNonNovel\t$currentDetectedNonNovel_regions\t$currentTotalRegions\t0\t$currentDetectedNonNovelICs\t0\tallTissues"|awk '{print $0"\t"($4/$5)*100}'>> stats/${spec}.MergedTargetRegions.proportionDetected.${class}
	echo -e "$spec\t$class\tdetectedNovel\t$currentDetectedNovel_regions\t$currentTotalRegions\t0\t$currentDetectedNovelICs\t0\tallTissues"|awk '{print $0"\t"($4/$5)*100}'>> stats/${spec}.MergedTargetRegions.proportionDetected.${class}  #log-> "\t"log($6)/log(2)
	echo -e "$spec\t$class\tdetectedALL\t$currentDetected\t$currentTotalRegions\t0\t$currentDetectedICs\t0\tallTissues"|awk '{print $0"\t"($4/$5)*100}'>> stats/${spec}.MergedTargetRegions.proportionDetected.${class}

	done
done
echo "XXdoneXX"
