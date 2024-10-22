# ---------------------------------------------------
# Read length
# ---------------------------------------------------

for capture in Mv2 MpreCap Hv3 HpreCap
do
    while read lab cap tmp sample
    do
        path="${lab}_${cap}_${tmp}_${sample}.gff.gz" #Path to the read mapping GFF
        details=$(echo -e "$lab\t$cap\t$tmp\t$sample")
        if [ -f $path ]; then
            zcat $path | egrep -v "SIRVome_isoforms|ERCC" | awk -F"\t" -v OFS="\t" -v det="$details" '{ split($9, ids, "\""); print $5-$4, ids[2], det }'
        fi
    done < samples.gencode.${capture}.tsv > gencode.readlength.${capture}.tsv
done

for capture in Mv2 MpreCap Hv3 HpreCap
do
    awk -F'\t' '{key = "\"" $2 "\"" "\t" $3 "\t" $4 "\t" $5 "\t" $6; sums[key]+=$1} END {for (key in sums) print key "\t" sums[key]}' gencode.readlength.$capture.tsv > gencode.readlength.$capture.finite.tsv
done

for i in gencode.readlength.*finite.tsv; do Rscript length_distribution.R $i; done
