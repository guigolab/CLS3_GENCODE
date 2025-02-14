####---------------merge bams from all tissues; 4 final bams per tech needed; HpreCap, Hv3, MpreCap, Mv2 

#!/usr/bin/bash
#$ -N OTRm
#$ -t 1-8
#$ -tc 8
#$ -q rg-el7,long-centos79
#$ -cwd
#$ -l virtual_free=32G,h_rt=6:00:00
#$ -pe smp 1
#$ -P prj006070
#$ -V
#$ -e errOTRm.$TASK_ID
#$ -o outOTRm.$TASK_ID

bamFiles=/users/rg/gkaur/LyRic_gencode_copy/output_SPLIT/mappings/longReadMapping
out=/users/project/gencode_006070/gkaur/CapTrap-CLSfigs/OTRresults

file=$(ls -1 $bamFiles | grep ".bam" | grep -v bai | awk -F"_" '{print $1"_"$2}' | sort | uniq | sed -n ${SGE_TASK_ID}p)

	samtools merge -o ${out}/${file}_0+_allTissues.bam ${bamFiles}/$file*bam
	echo -e "done for $file"

