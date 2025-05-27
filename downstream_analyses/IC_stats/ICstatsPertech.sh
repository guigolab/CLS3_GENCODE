for spec in Hv3 Mv2; do
	for tech in pacBioSII ont; do
		echo -e "spec\ttech\ttissue\tcapture\tstats\tcateg2"  > stats/$spec\_refined_$tech\_ALL_anchIC_stats
		echo -e "spec\ttech\ttissue\tcapture\tstats\tcateg2"  > stats/$spec\_refined_$tech\_NI_anchIC_stats
	done

cTs=data/masterTable/${spec}_splicedmasterTable_refined.gtf
cTu=data/masterTable/${spec}_unsplicedmasterTable_refined.gtf

while read meta;
do
        tissue=`echo $meta| awk -F" " '{print $2}'`
        scode=`echo $meta| awk -F" " '{print $1}'` 
        #SIDHBrAOC0301
	stage=`echo $meta| awk -F" " '{print $3}'`
	c1=`echo $scode|cut -b 1-8`
	c2=`echo $scode|cut -b 10-13`
	cap=`echo $scode|cut -b 9`
	tis=`echo $scode|cut -b 5-7`
	tech=`echo $meta| awk -F" " '{print $4}'`
	
	if [[ $tis == "BrA" ]];then     categ2="2"
	elif [[ $tis == "BrE" ]];then categ2="1"
	elif [[ $tis == "HeE" ]];then categ2="4"
	elif [[ $tis == "HeA" ]];then categ2="5"
	elif [[ $tis == "LiE" ]];then categ2="7"
	elif [[ $tis == "LiA" ]];then categ2="8"
	elif [[ $tis == "WbE" ]];then categ2="10"
	elif [[ $tis == "WbA" ]];then categ2="11"
	elif [[ $tis == "TeA" ]];then categ2="13"
	elif [[ $tis == "PlP" && $spec == "Hv3" ]] || [[ $tis == "TpA" && $spec == "Mv2" ]];then categ2="15"
	elif [[ $tis == "TpA" && $spec == "Hv3" ]];then categ2="17"
	elif [[ $tis == "CpA" ]];then categ2="19"
	else categ2="0"
	fi
	
	if [[ $cap == "C" ]]; then
		other="P"
        	stat_post=`cat ${cTs} ${cTu} | awk '{if($3 == "transcript")print $0}' | grep -c $scode`
        	stat_pre=`cat ${cTs} ${cTu} | awk '{if($3 == "transcript")print $0}' | grep -c $c1$other$c2`
        	stat_common=`cat ${cTs} ${cTu} | awk '{if($3 == "transcript")print $0}' | grep $scode| grep -c $c1$other$c2`
        	ACTUALstat_post=$(($stat_post-$stat_common))
		ACTUALstat_pre=$(($stat_pre-$stat_common))

		echo $scode $categ2 $tissue
		echo -e "$spec\t$tech\t$stage$tissue\tpre-capture\t$ACTUALstat_pre\t$categ2\n$spec\t$tech\t$stage$tissue\tcommon\t$stat_common\t$categ2\n$spec\t$tech\t$stage$tissue\tpost-capture\t$ACTUALstat_post\t$categ2" >> stats/$spec\_refined_$tech\_ALL_anchIC_stats
##NI --> Novel Intergenic
		stat_NI_post=`cat ${cTs} ${cTu} | awk '{if($3 == "transcript")print $0}' | awk '{if ($20 ~ /Intergenic/) print $0}' | grep -c $scode`
                stat_NI_pre=`cat ${cTs} ${cTu} | awk '{if($3 == "transcript")print $0}' | awk '{if ($20 ~ /Intergenic/) print $0}'| grep -c $c1$other$c2`
                stat_NI_common=`cat ${cTs} ${cTu} | awk '{if($3 == "transcript")print $0}' | awk '{if ($20 ~ /Intergenic/) print $0}' | grep $scode| grep -c $c1$other$c2`
		ACTUALstat_NI_post=$(($stat_NI_post-$stat_NI_common))
		ACTUALstat_NI_pre=$(($stat_NI_pre-$stat_NI_common))
		echo -e "$spec\t$tech\t$stage$tissue\tpre-capture\t$ACTUALstat_NI_pre\t$categ2\n$spec\t$tech\t$stage$tissue\tcommon\t$stat_NI_common\t$categ2\n$spec\t$tech\t$stage$tissue\tpost-capture\t$ACTUALstat_NI_post\t$categ2" >> stats/$spec\_refined_$tech\_NI_anchIC_stats

	fi
done < data/metadata/${spec}_metadata.tsv
done
