#!/usr/bin/bash

set -e
set -o pipefail

#----remove old data if any
rm -rf data/ plots/ stats/

#----download required files
mkdir -p data/ data/masterTable/ data/targets/ data/metadata/ data/annotation/ stats/ plots/
echo -e "#### created needed directories ####"

for spec in Hv3 Mv2; do
	#----download target catalog gtfs
        wget https://github.com/guigolab/gencode-cls-master-table/releases/download/Supplementary/${spec}_CLS3_targetDesign_mergedRegions.gtf.gz -O data/targets/${spec}_CLS3_targetDesign_mergedRegions.gtf.gz
        gunzip data/targets/${spec}_CLS3_targetDesign_mergedRegions.gtf.gz

	#----download reference annotations
	if [[ $spec == "Hv3" ]]; then
                v="v27"
        else	v="vM16"
        fi
	wget https://github.com/guigolab/cls3-final-files-description/releases/download/simplified.annotation.v1/${spec}.gencode${v}.simplified_biotypes.gtf.gz -O data/annotation/${spec}.gencode${v}.simplified_biotypes.gtf.gz
	gunzip data/annotation/${spec}.gencode${v}.simplified_biotypes.gtf.gz

	#----download CLS anchored models (anchTM) files
        wget https://github.com/guigolab/gencode-cls-master-table/releases/download/GencodeCLS_v3.0/${spec}_masterTable_refined.gtf.gz -O data/masterTable/${spec}_masterTable_refined.gtf.gz
       gunzip data/masterTable/${spec}*
done

#----data generation
#source activate gffC

#Hv3     detected        176435  23814   17491   nonIntersecting miTranscriptome 13.4973
echo -e "spec\ttype\ttotalALL\ttotal\tcount\tcateg\tclass\tperc\tperct" > stats/targs_biotypes

for spec in Hv3 Mv2;
do
	for type in detected all;	
	do	
		if [[ $spec == "Hv3" ]]; then
		v="v27"
	else	v="vM16"
	fi
	ALLtotal=`cat data/targets/${spec}_CLS3_targetDesign_mergedRegions.gtf | wc -l`

	if [[ $type == "all" ]]; then
		cat data/targets/${spec}_CLS3_targetDesign_mergedRegions.gtf > stats/tmp.targets.gtf
		bedtools intersect -wo -split -s -a stats/tmp.targets.gtf -b <(grep -w exon data/annotation/${spec}.gencode${v}.simplified_biotypes.gtf | grep -v ERCC) | awk '{print $10"\t"$24}' | sort --parallel=4 | uniq > stats/bdi_${v}_${type}
	else
		echo -e "detected!!!"
		cat data/masterTable/$spec\_masterTable_refined.gtf | grep -w transcript | grep -v "target \"\"" | awk -F"\"" '{print $6}' | grep . | perl -pi -e "s/,/\n/g" | sort --parallel=4 -T /tmp/ | uniq | grep -Ff - data/targets/${spec}_CLS3_targetDesign_mergedRegions.gtf > stats/tmp.targets.gtf
		bedtools intersect -wo -split -s -a stats/tmp.targets.gtf -b <(grep -w exon data/annotation/${spec}.gencode${v}.simplified_biotypes.gtf | grep -v ERCC) | awk '{print $10"\t"$24}' | sort --parallel=4 | uniq > stats/bdi_${v}_${type}
	fi
			
		for class in NONCODE miTranscriptome fantomCat gencodeLncRna refSeq bigTranscriptome CMfinderCRSs phyloCSF GWAScatalog UCE fantomEnhancers VISTAenhancers; # allCatalogs;
		do
	        	echo -e "$spec\t$class"
	        	if [[ $class == "allCatalogs" ]]; then
	                	get="cat"
				total=`cat data/targets/${spec}_CLS3_targetDesign_mergedRegions.gtf | wc -l`
	        	else
	                	get="grep $class"
				total=`cat stats/tmp.targets.gtf | eval "$get" | wc -l`
	       		fi

			ALLtotal=`grep ${class} data/targets/${spec}_CLS3_targetDesign_mergedRegions.gtf | wc -l`	

			#--------intersected v27 targets
			intersected=`cat stats/bdi_${v}_${type} | eval "$get" | cut -f1 | sort | uniq | wc -l`
			multiIntersects=`cat stats/bdi_${v}_${type} | eval "$get" | cut -f1,2 | cut -f1 | sort | uniq -c |awk '{if($1!="1")print $2;}' | sort | uniq | wc -l`
			#-------not intersecting v27 targets
			NonIntersected=$((total - intersected))

			echo -e "$spec\t$type\t$ALLtotal\t$total\t$NonIntersected\tnonIntersecting\t$class"| awk '{print $0"\t"$5/$4*100"\t"$5/$3*100}' >> stats/targs_biotypes
			#-- added miRNA, miscRNA, rRNA and multiIntersects to a new category "Others"
			#echo -e "$spec\t$type\t$ALLtotal\t$total\t$multiIntersects\tmultiIntersects\t$class"| awk '{print $0"\t"$5/$4*100"\t"$5/$3*100}' >> stats/targs_biotypes
			#----single intersects
			cat stats/bdi_${v}_${type} | eval "$get" | cut -f1,2 |cut -f1 | uniq -c | awk '{if($1=="1")print $2;}' |  grep -Ff - stats/bdi_${v}_detected | cut -f2 | perl -pi -e "s/misc_RNA|miRNA|rRNA/Others/g" | sort | uniq -c | perl -pi -e "s/\"|;//g" | awk -v type=$type -v At=$ALLtotal -v total=$total -v class=$class -v spec=$spec -v multiIntersects=$multiIntersects '{if($2~"Others")print spec"\t"type"\t"At"\t"total"\t"$1+multiIntersects"\t"$2"\t"class"\t"($1+multiIntersects)/total*100"\t"($1+multiIntersects)/At*100; else print spec"\t"type"\t"At"\t"total"\t"$1"\t"$2"\t"class"\t"$1/total*100"\t"$1/At*100;}' >> stats/targs_biotypes
		done
	done
done
