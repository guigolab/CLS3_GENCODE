# Plotting on-target rate (OTR) for the control ERCC spike-ins to check the efficacy of the capture/target panel

> Calculate on-target rate per tissue/sample:
#-------calculate OTR per tissue/sample
```
qsub onTargetRate.sh
```

#------merge across tissues for allTissues calculations
```
qsub mergeAcrossTissues.sh
```

#-------calculate OTR across all tissues/samples per capture type
```
qsub onTargetRateallTissues.sh
```

#-------prepare input file to plot OTR
```
echo -e "EnrichCateg\ttech\tspec\ttissue\tcap\tOTR\tEnrichment\tOTRperc" > Enrichment_phase3
for tech in pacBio ont; do
        for spec in CapTrap_H CapTrap_M; do
                while read tissue; do

                        for calc in OTRnonERCC OTRwithERCC OTRonlyERCC;do
                                if [[ $calc == "OTRnonERCC" ]];then col="cut -f12"; #echo 11$col; 
                                elif [[ $calc == "OTRwithERCC" ]]; then col="cut -f13"; #echo 12$col; 
                                elif [[ $calc == "OTRonlyERCC" ]]; then col="cut -f14"; #echo 13$col; 
                                fi;

                                pre=`cat ${tech}*${spec}pre*_${tissue}*_OTR| $col`; post=`cat ${tech}*${spec}v*_${tissue}*_OTR| $col`;
                                Enrichment=`echo "scale=1;$post/$pre" | bc| awk '{print $0}'`

                                echo -e "$calc\t$tech\t$spec\t$tissue\tpre\t$pre\t$Enrichment"| awk '{print $0"\t"$6*100}' >> Enrichment_phase3
                                echo -e "$calc\t$tech\t$spec\t$tissue\tpost\t$post\t$Enrichment"| awk '{print $0"\t"$6*100}' >> Enrichment_phase3
                        done
                done < <(ls ${tech}*${spec}*_OTR* | awk -F"_" '{print $4}' | awk -F"\." '{print $1}' |sort | uniq)
        done
done
```

#------------
```
Rscript OTR.R
```
