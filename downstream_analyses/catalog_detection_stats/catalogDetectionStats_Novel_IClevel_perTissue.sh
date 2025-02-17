#!/bin/bash

out=$currnoHiss/plots/TargetAnalaysis/stats/

#delete the output txt if already exists
rm stats/*.MergedTargetRegions.proportionDetected.Novel.Novel_perTissue_*

# Submit the job and wait for its completion
qsub -sync y catalogDetectionStatsNovel_perTissue.sh

# After the job finishes, continue with the following commands
for spec in Hv3 Mv2; do
        if [[ $spec == *Hv3* ]]; then
                code="H"
        else    code="M"
        fi
     cat ${out}/${spec}.MergedTargetRegions.proportionDetected.Novel_perTissue_SID${code}* | egrep -v "sncRNA|allCatalogs|uncaptured" >> ${out}/${spec}.refined_MergedTargetRegions.proportionDetected.Novel_perTissue
     echo -e "spec\tclass\tregionType\tdetectedRegions\tallRegions\tdetectedTMs\tdetectedICs\tdetectedLoci\ttissue\tpropDetectedRegions" > ${out}/${spec}.refined_MergedTargetRegions.proportionDetected.Novel_AllCatalogs
     cat ${out}/${spec}.MergedTargetRegions.proportionDetected.Novel_perTissue_SID${code}* | grep allCatalogs >> ${out}/${spec}.refined_MergedTargetRegions.proportionDetected.Novel_AllCatalogs
     echo -e "spec\tclass\tregionType\tdetectedRegions\tallRegions\tdetectedTMs\tdetectedICs\tdetectedLoci\ttissue\tpropDetectedRegions" > ${out}/${spec}.refined_MergedTargetRegions.proportionDetected.Novel_uncaptured
     cat ${out}/${spec}.MergedTargetRegions.proportionDetected.Novel_perTissue_SID${code}* | grep uncaptured >> ${out}/${spec}.refined_MergedTargetRegions.proportionDetected.Novel_uncaptured


done

source activate forPlotting

Rscript catalogDetectionStatsIntergenic_matrix.R #perTissue matrix plot

Rscript catalogDetectionStats_overallMatrix.R #overall matrices
