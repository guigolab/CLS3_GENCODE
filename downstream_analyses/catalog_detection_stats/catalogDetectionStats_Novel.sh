#!/bin/bash

#remove any existing result file
rm $currnoHiss/plots/TargetAnalaysis/stats/*.MergedTargetRegions.proportionDetected.Novel*

# Submit the job and wait for its completion
qsub -sync y catalogDetectionStats_Novel_IClevel.sh

# After the job finishes, continue with the following
cat $currnoHiss/plots/TargetAnalaysis/stats/Hv3.MergedTargetRegions.proportionDetected.Novel* | grep -v sncRNA >> $currnoHiss/plots/TargetAnalaysis/stats/Hv3.refined_MergedTargetRegions.proportionDetected.Novel_AllTissues

cat $currnoHiss/plots/TargetAnalaysis/stats/Mv2.MergedTargetRegions.proportionDetected.Novel* | grep -v sncRNA >> $currnoHiss/plots/TargetAnalaysis/stats/Mv2.refined_MergedTargetRegions.proportionDetected.Novel_AllTissues


source activate forPlotting

Rscript catalogDetectionStats_Novel_ICs.R
