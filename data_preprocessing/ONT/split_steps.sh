#!/usr/bin/bash 
#$ -N splitJointpp
#$ -t 1-88
#$ -tc 44
#$ -q rg-el7,long-centos79,short-centos79
#$ -cwd
#$ -P prj006070
#$ -l virtual_free=64G,h_rt=7:00:00
#$ -pe smp 8
#$ -V
#$ -e logsJoint_3/$TASK_ID.err
#$ -o logsJoint_3/$TASK_ID.out

input=/users/project/gencode_006070/gkaur/masterTablev2/noHiSS/table/artifacts/duplex_tools/
###mkdir $input/INPUT_FASTQS/ $input/INPUT_FILES/
###ln -s /users/rg/gkaur/LyRic_gencode_copy/fastqs/*fastq.gz /users/project/gencode_006070/gkaur/masterTablev2/noHiSS/table/artifacts/duplex_tools/INPUT_FILES/
###mkdir /users/project/gencode_006070_no_backup/gkaur/splitSamples_JointCadapter/ /users/project/gencode_006070_no_backup/gkaur/splitSamples_ONT_JointCadapter_1step/
file=$(ls -1  /users/project/gencode_006070/gkaur/masterTablev2/noHiSS/table/artifacts/duplex_tools/INPUT_FILES/*.gz | sed -n ${SGE_TASK_ID}p)

#echo $file

sample=`basename $file`
echo $sample
input=/users/project/gencode_006070/gkaur/masterTablev2/noHiSS/table/artifacts/duplex_tools/
###output=/users/project/gencode_006070_no_backup/gkaur/splitSamples_JointCadapter/
output_allDefault=/users/project/gencode_006070_no_backup/gkaur/splitSamples_ONT_JointCadapter_1step/
output_allDefault_2=/users/project/gencode_006070_no_backup/gkaur/splitSamples_ONT_JointCadapter_2step/

###mkdir $input/INPUT_FASTQS/$sample
mkdir $output/$sample/
###ln -s $input/INPUT_FILES/$sample $input/INPUT_FASTQS/$sample/$sample

##step1	#split using ONT seq adapter
conda activate forPlotting
###duplex_tools split_on_adapter --threads 4 --adapter_type ONT_sequencing_adapter $input/INPUT_FASTQS/$sample/ $output/$sample/split1/ Native

##step2 #split using CapTrapSeq adapter
###duplex_tools split_on_adapter --threads 4 --n_bases_to_mask_tail 0 --n_bases_to_mask_head 0 --degenerate_bases 0 --edit_threshold 26 --adapter_type CapTrap_joint $output/$sample/split1/ $output/$sample/split2/ PCR

##stepXX check the combimation of ONTadapter + Joint CapTrap in PCR more: this treats the ONTseqAdapter in Native mode and JointCapTrapAdapter in PCR mode:
##seems to be making the similar splits base on test runs 
duplex_tools split_on_adapter --threads 8 --adapter_type ONT_sequencing_adapter+CapTrapSeqJoint $input/INPUT_FASTQS/$sample/ ${output_allDefault}/$sample/ PCR

#split the resukt of splits from collated adapter using the CapTrapSea adapter to see if there are multisplits getting away (like an example seen)
duplex_tools split_on_adapter --threads 8 --n_bases_to_mask_tail 0 --n_bases_to_mask_head 0 --degenerate_bases 0 --edit_threshold 26 --adapter_type CapTrap_joint ${output_allDefault}/$sample/ ${output_allDefault_2}/$sample/ PCR
echo -e "done #####$sample#####"
