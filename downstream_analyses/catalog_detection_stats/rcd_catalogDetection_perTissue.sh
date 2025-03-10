#!/bin/bash
set -e
set -o pipefail

#----data generation
#source activate forPlotting
bash rcd_catalogDetection_perTissue_calc.sh
echo "calculations done"


#----collate into a final stats file
for spec in Hv3 Mv2; do
        if [[ $spec == *Hv3* ]]; then
                code="H"
        else    code="M"
        fi
     cat stats/${spec}.MergedTargetRegions.proportionDetected.perTissue_SID${code}* | egrep -v "sncRNA|allCatalogs|uncaptured" >> stats/${spec}.MergedTargetRegions.proportionDetected_perTissue
     echo -e "spec\tclass\tregionType\tdetectedRegions\tallRegions\tdetectedTMs\tdetectedICs\tdetectedLoci\ttissue\tpropDetectedRegions" > stats/${spec}.MergedTargetRegions.proportionDetected_AllCatalogs
     cat stats/${spec}.MergedTargetRegions.proportionDetected.perTissue_SID${code}* | grep allCatalogs >> stats/${spec}.MergedTargetRegions.proportionDetected_AllCatalogs
     echo -e "spec\tclass\tregionType\tdetectedRegions\tallRegions\tdetectedTMs\tdetectedICs\tdetectedLoci\ttissue\tpropDetectedRegions" > stats/${spec}.MergedTargetRegions.proportionDetected_uncaptured
     cat stats/${spec}.MergedTargetRegions.proportionDetected.perTissue_SID${code}* | grep uncaptured >> stats/${spec}.MergedTargetRegions.proportionDetected_uncaptured
	
     echo "$spec calculations done"
done

Rscript rcd_catalogDetectionStats_matrix.R #perTissue matrix plot
echo main matrix plot created
Rscript rcd_catalogDetectionStats_overallMatrix.R  #overall matrices
echo overall matrices created DONE
