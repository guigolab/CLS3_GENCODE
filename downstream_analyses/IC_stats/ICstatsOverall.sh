#!/bin/bash
set -e
set -o pipefail

#----remove old data if any
rm -rf data/ plots/ stats/

#----download required files
mkdir -p data/ data/masterTable/ data/metadata/ stats/ plots/
echo -e "#### created needed directories ####"

for spec in Hv3 Mv2; do #download CLS transcript (IC) files
      wget https://github.com/guigolab/gencode-cls-master-table/releases/download/GencodeCLS_v3.0/${spec}_splicedmasterTable_refined.gtf.gz -O data/masterTable/${spec}_splicedmasterTable_refined.gtf.gz
      wget https://github.com/guigolab/gencode-cls-master-table/releases/download/GencodeCLS_v3.0/${spec}_unsplicedmasterTable_refined.gtf.gz -O data/masterTable/${spec}_unsplicedmasterTable_refined.gtf.gz
      gunzip data/masterTable/${spec}*

       #download metadata
       wget https://github.com/guigolab/gencode-cls-master-table/releases/download/GencodeCLS_v1.0/${spec}_metadata.tsv.gz -O data/metadata/${spec}_metadata.tsv.gz
       #zcat data/metadata/${spec}_metadata.tsv.gz | cut -f1 | cut -b1-7 | sort | uniq > data/metadata/${spec}_samples
       echo -e "downloaded the needed files $spec"
       gunzip data/metadata/${spec}_metadata.tsv.gz


#----calculations
echo -e "spec\ttech\ttissue\tcapture\tstats\tcateg2" > stats/$spec\_refined_overall_ALL_anchIC_stats
echo -e "spec\ttech\ttissue\tcapture\tstats\tcateg2" > stats/$spec\_refined_overall_N_anchIC_stats
echo -e "spec\tcateg\tcapture\tstats" > stats/$spec\_refined_overall_ALL_anchIC_capstats
echo -e "spec\tcateg\tcapture\tstats" > stats/$spec\_refined_overall_N_anchIC_capstats

cTs=$currnoHiss/IntronChainMT/${spec}_splicedmasterTable_refined.gtf
cTu=$currnoHiss/IntronChainMT/${spec}_unsplicedmasterTable_refined.gtf
      
while read meta;
do
        tissue=`echo $meta| awk -F" " '{print $2}'`
        scode=`echo $meta| awk -F" " '{print $1}'`
        #SIDHBrAOC0301
        stage=`echo $meta| awk -F" " '{print $3}'`
        c0=`echo $scode|cut -b 1-7`
        c1=`echo $scode|cut -b 1-8`
        c2=`echo $scode|cut -b 10-13`
        cap=`echo $scode|cut -b 9`
        echo -e "$c1 $cap $c2"
        tech=`echo $meta| awk -F" " '{print $4}'`
        tis=`echo $scode|cut -b 5-7`
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


	if [[ $cap == "C" && $tech == "ont" ]]; then
		other="P"
		#tech=`echo $meta| awk -F" " '{print $4}'`
		stat_post=`cat ${cTs} ${cTu} | awk '{if($3 == "transcript")print $0}' | grep -c $c0\.$cap$c2`
        	stat_pre=`cat ${cTs} ${cTu} | awk '{if($3 == "transcript")print $0}' | grep -c $c0\.$other$c2`
        	stat_common=`cat ${cTs} ${cTu} | awk '{if($3 == "transcript")print $0}' | grep $c0\.$cap$c2 | grep -c $c0\.$other$c2`
        	ACTUALstat_post=$(($stat_post-$stat_common))
		ACTUALstat_pre=$(($stat_pre-$stat_common))	
		echo -e "$spec\tpacBio+ONT\t$stage$tissue\tpre-capture\t$ACTUALstat_pre\t$categ2\n$spec\tpacBio+ONT\t$stage$tissue\tcommon\t$stat_common\t$categ2\n$spec\tpacBio+ONT\t$stage$tissue\tpost-capture\t$ACTUALstat_post\t$categ2" >> stats/${spec}_refined_overall_ALL_anchIC_stats
echo "Overall $tissue done"
##NI --> Novel Intergenic
#		stat_NI_post=`cat ${cTs} ${cTu} | awk '{if($3 == "transcript")print $0}' | awk '{if ($20 ~ /Intergenic/) print $0}'| grep -c $c0\.$cap$c2`
#                stat_NI_pre=`cat ${cTs} ${cTu} | awk '{if($3 == "transcript")print $0}' | awk '{if ($20 ~ /Intergenic/) print $0}' | grep -c $c0\.$other$c2`
#                stat_NI_common=`cat ${cTs} ${cTu} | awk '{if($3 == "transcript")print $0}' | awk '{if ($20 ~ /Intergenic/) print $0}' | grep $c0\.$cap$c2 | grep -c $c0\.$other$c2`
#		ACTUALstat_NI_post=$(($stat_NI_post-$stat_NI_common))
#		ACTUALstat_NI_pre=$(($stat_NI_pre-$stat_NI_common))
#		echo -e "$spec\tpacBio+ONT\t$stage$tissue\tpre-capture\t$ACTUALstat_NI_pre\t$categ2\n$spec\tpacBio+ONT\t$stage$tissue\tcommon\t$stat_NI_common\t$categ2\n$spec\tpacBio+ONT\t$stage$tissue\tpost-capture\t$ACTUALstat_NI_post\t$categ2" >> stats/$spec\_refined_overall_NI_anchIC_stats
#echo "NI $tissue done"
##N --> Novel ###
               stat_N_post=`cat ${cTs} ${cTu} | awk '{if($3 == "transcript")print $0}' | awk '{if ($20 ~ /Intergenic/ || $20 ~ /Intronic/ || $20 ~ /Extends/ || $20 ~ /Antisense/ || $20 ~ /revIntronic/ || $20 ~ /runOn/) print $0}' | grep -c $c0\.$cap$c2`
                stat_N_pre=`cat ${cTs} ${cTu} | awk '{if($3 == "transcript")print $0}' | awk '{if ($20 ~ /Intergenic/ || $20 ~ /Intronic/ || $20 ~ /Extends/ || $20 ~ /Antisense/ || $20 ~ /revIntronic/ || $20 ~ /runOn/) print $0}' | grep -c $c0\.$other$c2`
                stat_N_common=`cat ${cTs} ${cTu} | awk '{if($3 == "transcript")print $0}' | awk '{if ($20 ~ /Intergenic/ || $20 ~ /Intronic/ || $20 ~ /Extends/ || $20 ~ /Antisense/ || $20 ~ /revIntronic/ || $20 ~ /runOn/) print $0}' | grep $c0\.$cap$c2 | grep -c $c0\.$other$c2`
               ACTUALstat_N_post=$(($stat_N_post-$stat_N_common))
               ACTUALstat_N_pre=$(($stat_N_pre-$stat_N_common))
               echo -e "$spec\tpacBio+ONT\t$stage$tissue\tpre-capture\t$ACTUALstat_N_pre\t$categ2\n$spec\tpacBio+ONT\t$stage$tissue\tcommon\t$stat_N_common\t$categ2\n$spec\tpacBio+ONT\t$stage$tissue\tpost-capture\t$ACTUALstat_N_post\t$categ2" >> stats/$spec\_refined_overall_N_anchIC_stats
echo "N $tissue done"


	fi
done < $currnoHiss/${spec}_sampleCoding


for div in ALL N;
do	if [[ $div == "ALL" ]]; then
        	novelty() {
                grep SID
                }
        elif [[ $div == "N" ]]; then
                novelty() {
                awk '{if ($20 ~ /Intergenic/ || $20 ~ /Intronic/ || $20 ~ /Extends/ || $20 ~ /Antisense/ || $20 ~ /revIntronic/ || $20 ~ /runOn/) print $0}'
                }
        fi
	for cap in preCap common postCap;
        do	if [[ $cap == "preCap" ]]; then
                capC() {
                grep SID.....P.... | grep -v SID.....C....
                }
        elif [[ $cap == "common" ]]; then
                capC() {
                grep SID.....C.... | grep SID.....P....
                }
        elif [[ $cap == "postCap" ]]; then
                capC() {
                grep SID.....C.... | grep -v SID.....P....
                }
        fi
	stat=`cat ${cTs} ${cTu} | awk '{if($3 == "transcript")print $0}'|novelty| capC | awk -F"\"" '{print $2}'| sort |uniq|wc -l`
	echo -e "$spec\t$div\t$cap\t$stat" >> stats/$spec\_refined_overall_${div}_anchIC_capstats
	done
done
done
