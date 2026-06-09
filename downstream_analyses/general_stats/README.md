To obtain tagged "v47-CLS3" mapping table. 
The mapping table is further tagged with source sample(s), and source target catalog(s) to aid in the downstream analyses.
Included further mapped sample and target source for the models (anchTMs) from an older version of the release.

```
smT="../../IntronChainMT/Hv3_splicedmasterTable_refined.gtf" 
umT="../../IntronChainMT/Hv3_unsplicedmasterTable_refined.gtf"

python tagv47-CLSmap.anchIC.SID.py <(cat ${smT} ${umT}) <(grep -v geneID_v47 data/tmp.mainChr.novel.transcripts.map) oldIDs.gtf > data/tmp.mainChr.novel.tagged.map
```

To obtain novel TSSs subsets tagged with anchICs and SIDs
```
clsset="novel.transcripts"

python tagTSSs.anchIC.SID.py <(cat ${smT} ${umT}) data/tmp.${clsset}.TSS.bed data/tmp.mainChr.${clsset}.map oldIDs.gtf > data/${clsset}.TSS.anchIC-SIDtagged.bed
```

