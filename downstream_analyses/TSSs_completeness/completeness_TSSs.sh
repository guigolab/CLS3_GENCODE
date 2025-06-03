#get novel TMs -> antisense, extends, extends, intergenic
set -e
set -o pipefail

#----remove old data if any
rm -rf data/ plots/ stats/

#----download required files
mkdir -p data/masterTable/ data/proCapScores/ data/recount/ stats/ stats/filtered_ICs/ stats/filtered_ICtables/ plots/ data/endSupport/ data/genomes/ data/cage/

echo -e "#### created needed directories ####"

	#--TSS sets
	baselink="https://github.com/guigolab/cls3-final-files-description/releases/download/data.files.v1/"
	wget ${baselink}/novelTSS_anchTMs_loci.TSS.bed -O data/novelTSS_anchTMs_loci.TSS.bed
	wget ${baselink}/annotation_proteincoding_genes.TSS.bed -O data/proteincodingTSS_anchTMs_loci.TSS.bed
        wget ${baselink}/annotation_lncRNA_genes.TSS.bed -O data/lncRNATSS_anchTMs_loci.TSS.bed
	wget ${baselink}/decoy.TSS.bed -O data/decoyTSS_anchTMs_loci.TSS.bed
	#--genome hg38 sizes
	wget ${baselink}/hg38.sorted.genome -O data/genomes/hg38.sorted.genome
	#--human cage support
	wget ${baselink}/hg38_fair+new_CAGE_peaks_phase1and2.bed -O data/cage/hg38_fair+new_CAGE_peaks_phase1and2.bed
	#--proCapNet scores
	wget ${baselink}/proCapScores.zip -O data/proCapScores/proCapScores.zip
	unzip data/proCapScores/proCapScores.zip
echo "downloaded needed data files"

	
echo -e "spec\tcateg\tcount\tdataset\tproCapfilter\ttotal\tperc" > stats/completenessSupport_novelTSSs.txt
echo -e "spec\tcateg\tcount\tdataset\tproCapfilter\ttotal\tperc" > stats/completenessSupport_decoyTSSs.txt

#----Hv3----#
cat data/novelTSS_anchTMs_loci.TSS.bed | awk '{print $1"_"$2"_"$3"_"$6}' | sort | uniq > data/Hv3_novelTSS_anchTMs_loci.TSS.bed
cat data/decoyTSS_anchTMs_loci.TSS.bed | awk '{print $1"_"$2"_"$3"_"$6}' | sort | uniq > data/Hv3_decoyTSS_anchTMs_loci.TSS.bed

for set in proteincoding lncRNA;do
	echo -e "spec\tcateg\tcount\tdataset\tproCapfilter\ttotal\tperc" > stats/completenessSupport_${set}TSSs.txt
	cat data/${set}TSS_anchTMs_loci.TSS.bed | awk '{print $1"_"$2"_"$3"_"$6}' | sort | uniq > data/Hv3_${set}TSS_anchTMs_loci.TSS.bed
done

#-----------#

for type in novel decoy proteincoding lncRNA;
do
     for spec in Hv3;
     do
	total=`cat data/${spec}_${type}TSS_anchTMs_loci.TSS.bed | awk '{print $1"_"$2"_"$3"_"$6}' | sort | uniq | wc -l`
        #cage supported TSSs
	cat data/${type}TSS_anchTMs_loci.TSS.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t0\t"$6}' > data/endSupport/${type}TSSs_Hv3.5pEnds.bed

	cat data/endSupport/${type}TSSs_Hv3.5pEnds.bed | sort -T /tmp/ -k1,1 -k2,2n -k3,3n  | bedtools slop -s -l 50 -r 50 -i stdin -g data/genomes/hg38.sorted.genome  | bedtools intersect -u -s -a stdin -b data/cage/hg38_fair+new_CAGE_peaks_phase1and2.bed | awk '{print $1"_"$2"_"$3"_"$6}' | sort | uniq > data/endSupport/support_cage_${type}TSSs

	#proCapNet supported TSSs       
        if [[ $type == "novel" ]]; then
		dir=data/proCapScores/novel/
	elif [[ $type == "proteincoding" ]]; then
		dir=data/proCapScores/proteincoding_v27/
	elif [[ $type == "lncRNA" ]]; then
		dir=data/proCapScores/lncrna_v27/
	else #[[ $type == "decoy" ]]
		dir=data/proCapScores/decoys/
	fi

	grep -Ff data/${spec}_${type}TSS_anchTMs_loci.TSS.bed ${dir}/support_proCap_scoreMTE5_100bpW > ${dir}/support_proCap_scoreMTE5_100bpW_${type}TSSs

#---calculate completeness
for proCapfilter in MTE5;
        do
                dataset="cage"
                s5_support=`cat data/endSupport/support_cage_${type}TSSs | wc -l`
                sno_support=$((total - s5_support))
		#sno_support=`grep -v -Ff $currnoHiss/endSupport/support_cage_${type}TSSs $currnoHiss/data/${spec}_${type}TSS_anchTMs_loci.TSS.bed  | wc -l`
                echo -e "$spec\t5_support\t$s5_support\t$dataset\t$proCapfilter\t$total"|awk '{print $0"\t"$3/$6*100}' >> stats/completenessSupport_${type}TSSs.txt
                echo -e "$spec\tno_support\t$sno_support\t$dataset\t$proCapfilter\t$total"|awk '{print $0"\t"$3/$6*100}' >> stats/completenessSupport_${type}TSSs.txt

		dataset="proCap"
                s5_support=`cat ${dir}/support_proCap_scoreMTE5_100bpW_${type}TSSs | wc -l`
                sno_support=$((total - s5_support))
                #sno_support=`grep -v -Ff ${proCap}/support_proCap_scoreMTE5_100bpW_${type}TSSs $currnoHiss/data/${spec}_${type}TSS_anchTMs_loci.TSS.bed | wc -l`
                echo -e "$spec\t5_support\t$s5_support\t$dataset\t$proCapfilter\t$total"|awk '{print $0"\t"$3/$6*100}' >> stats/completenessSupport_${type}TSSs.txt
                echo -e "$spec\tno_support\t$sno_support\t$dataset\t$proCapfilter\t$total"|awk '{print $0"\t"$3/$6*100}' >> stats/completenessSupport_${type}TSSs.txt

      		dataset="CageProCap"
                s5_support=`cat ${dir}/support_proCap_scoreMTE5_100bpW_${type}TSSs data/endSupport/support_cage_${type}TSSs | sort | uniq | wc -l`
                sno_support=$((total - s5_support))
                #sno_support=`grep -v -Ff <(cat ${proCap}/support_proCap_scoreMTE5_100bpW_${type}TSSs $currnoHiss/endSupport/support_cage_${type}TSSs) $currnoHiss/data/${spec}_${type}TSS_anchTMs_loci.TSS.bed | sort | uniq | wc -l`
                echo -e "$spec\t5_support\t$s5_support\t$dataset\t$proCapfilter\t$total"|awk '{print $0"\t"$3/$6*100}' >> stats/completenessSupport_${type}TSSs.txt
                echo -e "$spec\tno_support\t$sno_support\t$dataset\t$proCapfilter\t$total"|awk '{print $0"\t"$3/$6*100}' >> stats/completenessSupport_${type}TSSs.txt

	done
done

done
