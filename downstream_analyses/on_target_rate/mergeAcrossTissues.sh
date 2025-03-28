#!/bin/bash
set -e
set -o pipefail

####---------------merge bams from all tissues; 4 final bams per tech needed; HpreCap, Hv3, MpreCap, Mv2 

##--create directories
mkdir -p data/bam stats/ plots/ out/
##---download bam files
#making files available pending

bamFiles=data/bam/

for tech in pacBioSII-Cshl-CapTrap ont-Crg-CapTrap; do
       for cap in HpreCap Hv3 MpreCap Mv2; do
               samtools merge --threads 4 -o out/${tech}_${cap}_0+_allTissues.bam ${bamFiles}/${tech}_${cap}*bam
               echo -e "${tech}_${cap}_0+_allTissues.bam generated"
       done
done
