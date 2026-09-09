working=/home/basiek/Projects/CapTrap-CLS.2026/upsetplot

cd $working

#cp /home/basiek/Pobrane/upSet/.temp3/* /home/basiek/Projects/CapTrap-CLS.2026/upsetplot

cp /home/basiek/Pobrane/upSet/Data/H_catalogs/*.gtf.gz .

cp chess3_1_2_GRCh38_lncRNA.gtf.gz old/
cp lncRNA_LncBookv2_0_GRCh38.gtf.gz old/ 

mv chess3_1_2_GRCh38_lncRNA.gtf.gz chess.gtf.gz
mv lncRNA_LncBookv2_0_GRCh38.gtf.gz lncBookv2_0.gtf.gz

for lab in chess lncBookv2_0
do
    zcat $lab.gtf.gz | /home/basiek/Projects/CapTrap-CLS.2026/upsetplot/gff2bed_full.pl - | awk 'BEGIN{FS=OFS="\t"}$6!="."' | awk '$1 ~ /^chr[0-9XY]{1,2}$/ {print $0}' | grep "\S" | /home/basiek/Projects/CapTrap-CLS.2026/upsetplot/bed12togff - > $lab.noBL.gtf
done

#zcat ../../Data/H_catalogs/chess3_1_2_GRCh38_lncRNA.gtf.gz | $utils/gff2bed_full.pl - | awk 'BEGIN{FS=OFS="\t"}$6!="."' | awk '$1 ~ /^chr[0-9XY]{1,2}$/ {print $0}' | grep "\S" | $utils/bed12togff - > chess3_1_2.hg38.gtf
#zcat ../../Data/H_catalogs/lncRNA_LncBookv2_0_GRCh38.gtf.gz | $utils/gff2bed_full.pl - | awk 'BEGIN{FS=OFS="\t"}$6!="."' | awk '$1 ~ /^chr[0-9XY]{1,2}$/ {print $0}' | grep "\S" | $utils/bed12togff - > lncBookv2_0.hg38.gtf

# -------------------------------------------------
# other
# -------------------------------------------------

echo -e "bigTranscriptome\nfantomCat\ngencode20+\nmiTranscriptome\nNONCODE\nrefSeq" > other.tsv

while read lab
do
    echo $lab
    zcat $lab.gtf.gz | awk 'BEGIN{FS=OFS="\t"}$7!="."' | awk '$1 ~ /^chr[0-9XY]{1,2}$/ {print $0}' | grep "\S" | awk '$3=="exon"' > $lab.noBL.gtf
done < other.tsv

# -------------------------------------------------
# lncRNA merge - exclude chrM, keep unspliced
# -------------------------------------------------
#cat NONCODE.hg38.gtf refSeq.hg38.gtf miTranscriptome.hg38.gtf gencode20+.hg38.gtf bigTranscriptome.hg38.gtf fantomCat.hg38.gtf | $utils/gff2bed_full.pl - | sort | uniq | grep -vFw "ENST00000610119.1" | awk 'BEGIN{FS=OFS="\t"}$6!="."' | awk '$1 ~ /^chr[0-9XY]{1,2}$/ {print $0}' | grep "\S" | $utils/bed12togff - | $utils/sortgff - | $utils/tmerge --exonOverhangTolerance 8 --tmPrefix IC - | grep '\S' > H_concat_catalogs.hg38.gtf


zcat NONCODE.gtf.gz refSeq.gtf.gz miTranscriptome.gtf.gz gencode20+.gtf.gz bigTranscriptome.gtf.gz fantomCat.gtf.gz | /home/basiek/Projects/CapTrap-CLS.2026/upsetplot/gff2bed_full.pl - | sort | uniq | grep -vFw "ENST00000610119.1" | awk 'BEGIN{FS=OFS="\t"}$6!="."' | awk '$1 ~ /^chr[0-9XY]{1,2}$/ {print $0}' | grep "\S" | /home/basiek/Projects/CapTrap-CLS.2026/upsetplot/bed12togff - | /home/basiek/Projects/CapTrap-CLS.2026/upsetplot/sortgff - > all1one.gtf

/home/basiek/Projects/CapTrap-CLS.2026/upsetplot/tmerge/tmerge --exonOverhangTolerance 8 all1one.gtf | grep '\S' > lncRNAmerge.tmp.gtf

./check_fix_gtf.sh lncRNAmerge.tmp.gtf lncRNAmerge.fix.gtf

bedtools intersect -s -wao -a lncRNAmerge.fix.gtf -b lncRNAmerge.fix.gtf | /home/basiek/Projects/CapTrap-CLS.2026/upsetplot/buildLoci.pl - > lncRNAmerge.buildLoci.tmp.gtf

./check_fix_gtf.sh lncRNAmerge.buildLoci.tmp.gtf lncRNAmerge.buildLoci.gtf

gzip lncRNAmerge.buildLoci.gtf


#zcat NONCODE.gtf.gz refSeq.gtf.gz miTranscriptome.gtf.gz gencode20+.gtf.gz bigTranscriptome.gtf.gz fantomCat.gtf.gz | grep -v ENST00000610119.1 | grep -v ENST00000471935.1|  awk '!/^#/' | awk '$7!="."' | sort -k1,1 -k4,4n > all1one.gtf

#./filter_spliced_only.sh all1one.gtf all1one.spliced.gtf

#/home/basiek/Projects/CapTrap-CLS.2026/upsetplot/tmerge/tmerge --exonOverhangTolerance 8 all1one.spliced.gtf | awk '$7!="."' > lncRNAmerge.tmp.gtf

#./check_fix_gtf.sh lncRNAmerge.tmp.gtf lncRNAmerge.fix.gtf

#bedtools intersect -s -wao -a lncRNAmerge.fix.gtf -b lncRNAmerge.fix.gtf | /home/basiek/Projects/CapTrap-CLS.2026/upsetplot/buildLoci.pl - > lncRNAmerge.buildLoci.tmp.gtf

#./check_fix_gtf.sh lncRNAmerge.buildLoci.tmp.gtf lncRNAmerge.buildLoci.gtf

#gzip lncRNAmerge.buildLoci.gtf

# -----------------------------
# buildLoci
# -----------------------------

# get samples
for file in `ls *gtf.gz`
do
    lab=`basename $file| awk -F ".gtf.gz" '{print $1}'`
    echo $lab
done | grep -v gen27 | grep -v gen47 > samples.bl.tsv

while read lab
do
    echo $lab
    bedtools intersect -s -wao -a $lab.noBL.gtf -b $lab.noBL.gtf | /home/basiek/Projects/CapTrap-CLS.2026/upsetplot/buildLoci.pl - > $lab.buildLoci.tmp.gtf
    ./check_fix_gtf.sh $lab.buildLoci.tmp.gtf $lab.buildLoci.gtf
    gzip $lab.buildLoci.gtf 
done < samples.bl.tsv

    
for lab in gen27 gen47
do
    cp $lab.gtf.gz $lab.buildLoci.gtf.gz
done 

# -----------------------------
# get samples
# -----------------------------

for file in `ls *.buildLoci.gtf.gz`
do
    lab=`basename $file| awk -F ".buildLoci.gtf.gz" '{print $1}'`
    echo $lab
done > samples.tsv

# -----------------------------
# generate exon records in bed6
# -----------------------------

while read lab
do
zcat $lab.buildLoci.gtf.gz | gawk -v cat="$lab" 'BEGIN{FS=OFS="\t"} $3=="exon"{
  match($9, /gene_id "([^"]+)"/, m)
  print $1, $4-1, $5, cat"|"m[1], 0, $7
}' | sort -k1,1 -k2,2n > $lab.exons.bed
done < samples.tsv

# check
while read lab
do
    echo $lab
    zcat $lab.buildLoci.gtf.gz| awk '$3=="exon"'| wc -l
    cat $lab.exons.bed| wc -l
done < samples.tsv

# geneNb
while read lab
do
    echo $lab
    zcat $lab.buildLoci.gtf.gz| awk '$3=="exon" {print $10}'|sort| uniq | wc -l
done < samples.tsv


# -----------------------------
# merge exons
# -----------------------------

while read lab
do
    ./merge_exons_per_gene.sh $lab.exons.bed $lab.footprints.bed
done < samples.tsv

# -----------------------------
# consensus_matrix.tsv
# -----------------------------

#pip install pandas networkx --break-system-packages


python3 build_consensus_upset_matrix.refrmtd.py \
  --catalogs bigTranscriptome=bigTranscriptome.footprints.bed \
             CHESS_3.1.2=chess.footprints.bed \
             FantomCat=fantomCat.footprints.bed \
             GENCODE27=gen27.footprints.bed \
             GENCODE47=gen47.footprints.bed \
             "GENCODE20+=gencode20+.footprints.bed" \
             lncBookv2.0=lncBookv2_0.footprints.bed \
             MiTranscriptome=miTranscriptome.footprints.bed \
             NONCODE=NONCODE.footprints.bed \
             refSeq=refSeq.footprints.bed \
             lncRNA-merge=lncRNAmerge.footprints.bed \
  --min-frac 0.95 \
  --outdir results
  
python3 generate_upset_plot.trimm.sort.size.bars.py \
  --matrix results/consensus_matrix.tsv \
  --catalogs CHESS_3.1.2 lncBookv2.0 bigTranscriptome "GENCODE20+" refSeq \
             FantomCat MiTranscriptome NONCODE lncRNA-merge GENCODE27 GENCODE47 \
  --out results/upset_final_fixed.bars.pdf \
  --max-subset-rank 20 \
  --sort-categories-by input \
  --width 14 --height 9 \
  --totals-plot-elements 6 \
  --totals-scale-max 120000 \
  --totals-bar-height 0.5 \
  --totals-axis-style clean \
  --title "Consensus lncRNA loci across annotation catalogs"
  
  
:<<'END'
python3 generate_upset_plot.trimm.sort.size.py \
  --matrix results/consensus_matrix.tsv \
  --catalogs CHESS_3.1.2 lncBookv2.0 bigTranscriptome "GENCODE20+" refSeq \
             FantomCat MiTranscriptome NONCODE lncRNA-merge GENCODE27 GENCODE47 \
  --out results/upset_final_sized.pdf \
  --max-subset-rank 30 \
  --sort-categories-by input \
  --width 14 --height 9 \
  --title "Consensus lncRNA loci across annotation catalogs"
  

python3 generate_upset_plot.trimm.sort.py \
  --matrix results/consensus_matrix.tsv \
  --catalogs CHESS_3.1.2 lncBookv2.0 "GENCODE20+" bigTranscriptome refSeq \
             FantomCat MiTranscriptome NONCODE lncRNA-merge GENCODE27 GENCODE47 \
  --out results/upset_final.pdf \
  --max-subset-rank 30 \
  --sort-categories-by input \
  --title "Consensus lncRNA loci across annotation catalogs"
  
  
python3 generate_upset_plot.trim.py \
  --matrix results/consensus_matrix.tsv \
  --catalogs GENCODE47 GENCODE27 bigTranscriptome FantomCat "GENCODE20+" \
             MiTranscriptome NONCODE refSeq CHESS_3.1.2 lncBookv2.0 lncRNA-merge \
  --out results/upset_final.pdf \
  --max-subset-rank 30 \
  --title "Consensus lncRNA loci across annotation catalogs"
  
    
    
    
    
  
pip install upsetplot matplotlib "pandas<2.2" --break-system-packages

python3 generate_upset_plot.py \
  --matrix results/consensus_matrix.tsv \
  --catalogs gen47 gen27 bigTranscriptome fantomCat "gencode20+" \
             miTranscriptome NONCODE refSeq chess lncBookv2_0 \
  --out results/upset_final.png \
  --title "Consensus lncRNA loci across annotation catalogs"
  
python3 build_consensus_upset_matrix.refrmtd.py \
  --catalogs GENCODE47=gen47.footprints.bed GENCODE27=gen27.footprints.bed \
             lncRNAmerge=merge.footprints.bed NONCODE=noncode.footprints.bed \
             MiTranscriptome=mitranscriptome.footprints.bed RefSeq=refseq.footprints.bed \
  --min-frac 0.6 \
  --outdir results
  
END  
