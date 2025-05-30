	echo -e "spec\tclass\tregions\tallTargs\tprop" > stats/MergedTargetRegions.proportionAllTogetherACTUAL #stats/$spec\.MergedTargetRegions.proportionAllTogetherACTUAL
for spec in Hv3 Mv2; do

        if [[ $spec == "Hv3" ]]; then
                sp="hs"
        else
                sp="mm"
        fi

	#---- # of individual target elements detected-----#
        cat data/masterTable/${spec}_masterTable_refined.gtf | grep -w transcript | grep SID.....C.... | grep -v SID.....P.... | grep -v "target \"\"" | awk -F"\"" '{print $6}'| perl -pi -e "s/,/\n/g" | sort | uniq > stats/$spec\_postCap_targets
        cat data/masterTable/${spec}_masterTable_refined.gtf | grep -w transcript | grep SID.....C.... | grep SID.....P.... | grep -v "target \"\"" | awk -F"\"" '{print $6}'| perl -pi -e "s/,/\n/g" | sort | uniq > stats/$spec\_common_targets
        cat data/masterTable/${spec}_masterTable_refined.gtf | grep -w transcript | grep -v SID.....C.... | grep SID.....P.... | grep -v "target \"\"" | awk -F"\"" '{print $6}'| perl -pi -e "s/,/\n/g" | sort | uniq > stats/$spec\_preCap_targets
echo step1

        #----- # of merged targets detected-----#
        cat data/targets/${spec}_CLS3_targetDesign_mergedRegions.gtf | awk -F"\"" '{print $2}'| grep -Ff stats/${spec}_postCap_targets - > stats/${spec}_targetRegions_postCapOnly
        cat data/targets/${spec}_CLS3_targetDesign_mergedRegions.gtf | awk -F"\"" '{print $2}'| grep -Ff stats/${spec}_common_targets - > stats/${spec}_targetRegions_commonOnly
        cat data/targets/${spec}_CLS3_targetDesign_mergedRegions.gtf | awk -F"\"" '{print $2}'| grep -Ff stats/${spec}_preCap_targets - > stats/${spec}_targetRegions_preCapOnly

	#----get the rest of the commons
	cat stats/${spec}_targetRegions_postCapOnly stats/${spec}_targetRegions_preCapOnly stats/${spec}_targetRegions_commonOnly | sort | uniq -c | awk '{if($1>1)print $2;}' > stats/${spec}_restCOMMON
	#19631
	cat stats/${spec}_restCOMMON stats/${spec}_targetRegions_commonOnly | sort | uniq > stats/${spec}_targetRegions_allCOMMON

	#--plus other already extracted common:
	common_targets_count=`cat stats/${spec}_targetRegions_allCOMMON | wc -l`
	#19971 OK
	preCap_targets_count=`cat stats/${spec}_targetRegions_postCapOnly stats/${spec}_targetRegions_allCOMMON | grep -c -v -Ff - stats/${spec}_targetRegions_preCapOnly`
	#706
	postCap_targets_count=`cat stats/${spec}_targetRegions_preCapOnly stats/${spec}_targetRegions_allCOMMON | grep -c -v -Ff - stats/${spec}_targetRegions_postCapOnly`
	#45001

	allRegions=`wc -l data/targets/${spec}_CLS3_targetDesign_mergedRegions.gtf | awk -F" " '{print $1}'`
	detected=`cat stats/${spec}_targetRegions_preCapOnly stats/${spec}_targetRegions_commonOnly stats/${spec}_targetRegions_postCapOnly | sort | uniq | wc -l`

	###-----PRINT into txt files-----#
	echo step2
        allRegions=`wc -l data/targets/${spec}_CLS3_targetDesign_mergedRegions.gtf | awk -F" " '{print $1}'`

        echo -e "$spec\tpre-capture\t$preCap_targets_count\t$allRegions"|awk '{print $0"\t"($3/$4)*100}'>>stats/MergedTargetRegions.proportionAllTogetherACTUAL #stats/$spec\.MergedTargetRegions.proportionAllTogetherACTUAL
	echo -e "$spec\tcommon\t$common_targets_count\t$allRegions"|awk '{print $0"\t"($3/$4)*100}'>>stats/MergedTargetRegions.proportionAllTogetherACTUAL #stats/$spec\.MergedTargetRegions.proportionAllTogetherACTUAL
	echo -e "$spec\tpost-capture\t$postCap_targets_count\t$allRegions"|awk '{print $0"\t"($3/$4)*100}'>>stats/MergedTargetRegions.proportionAllTogetherACTUAL #stats/$spec\.MergedTargetRegions.proportionAllTogetherACTUAL
	echo -e "$spec\tallDetected\t$detected\t$allRegions"|awk '{print $0"\t"($3/$4)*100}'>>stats/MergedTargetRegions.proportionAllTogetherACTUAL #stats/$spec\.MergedTargetRegions.proportionAllTogetherACTUAL
done
