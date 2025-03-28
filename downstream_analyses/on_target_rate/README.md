# Plotting on-target rate (OTR) for the control ERCC spike-ins to check the efficacy of the capture/target panel

## merge across tissues for allTissues calculations
```
bash mergeAcrossTissues.sh
```

## Calculate OTR per tissue/sample:
```
bash onTargetRate.sh
```

## prepare input file to plot OTR
```
for tech in pacBio ont; do
        for spec in CapTrap_H CapTrap_M; do
                while read tissue; do
                        #assign categories used to define breaks for clubbing together tissues.
                        if [[ $tissue == "Brain"* ]];then     category="2" 
                        elif [[ $tissue == "EmbBrain"* ]];then category="1"
                        elif [[ $tissue == "EmbHeart"* ]];then category="4"
                        elif [[ $tissue == "Heart"* ]];then category="5"
                        elif [[ $tissue == "EmbLiver"* ]];then category="7"
                        elif [[ $tissue == "Liver"* ]];then category="8"
                        elif [[ $tissue == "iPSC"* ]] || [[ $tissue == "EmbSC"* ]];then category="10"
                        elif [[ $tissue == "WBlood"* ]];then category="11"
                        elif [[ $tissue == "Testis"* ]];then category="13"
                        elif [[ $tissue == "Placenta"* && $spec == "CapTrap_H" ]] || [[ $tissue == "Tpool"* && $spec == "CapTrap_M" ]];then category="15"
                        elif [[ $tissue == "Tpool"* && $spec == "CapTrap_H" ]] || [[ $tissue == "allTissues"* && $spec == "CapTrap_M" ]];then category="17"
                        elif [[ $tissue == "Cpool"* ]];then category="19"
                        elif [[ $tissue == "allTissues" && $spec == "CapTrap_H" ]];then category="21"
                        else category="0"
                        fi

                        for calc in OTRnonERCC OTRwithERCC OTRonlyERCC;do
                                if [[ $calc == "OTRnonERCC" ]];then col="cut -f12"; #echo 11$col; 
                                elif [[ $calc == "OTRwithERCC" ]]; then col="cut -f13"; #echo 12$col; 
                                elif [[ $calc == "OTRonlyERCC" ]]; then col="cut -f14"; #echo 13$col; 
                                fi;
                                echo "$spec $tech $tissue"      
                                pre=`cat ${tech}*${spec}pre*_${tissue}*_OTR| $col`; post=`cat ${tech}*${spec}v*_${tissue}*_OTR| $col`;
                                Enrichment=`echo "scale=1;$post/$pre" | bc| awk '{print $0}'`

                                echo -e "$calc\t$tech\t$spec\t$tissue\tpre\t$pre\t$Enrichment"| awk -v categ=$category '{print $0"\t"$6*100"\t"categ}' >> stats/Enrichment_phase3
                                echo -e "$calc\t$tech\t$spec\t$tissue\tpost\t$post\t$Enrichment"| awk -v categ=$category '{print $0"\t"$6*100"\t"categ}' >> stats/Enrichment_phase3
                        done
                done < <(ls ${tech}*${spec}*_OTR* | awk -F"_" '{print $4}' | awk -F"\." '{print $1}' |sort|uniq)
        done
done
```

## Plot the results
```
Rscript OTR.R
Rscript OTR_ERCConly.R
```
