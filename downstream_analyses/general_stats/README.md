To obtain tagged "v47-CLS3" mapping table. 
The mapping table is further tagged with source sample(s), and source target catalog(s) to aid in the downstream analyses.

```
smT="../../IntronChainMT/Hv3_splicedmasterTable_refined.gtf" 
umT="../../IntronChainMT/Hv3_unsplicedmasterTable_refined.gtf"

python $currnoHiss/revision/TSS/tagv47-CLSmap.anchIC.SID.py <(cat ${smT} ${umT}) <(grep -v geneID_v47 $currnoHiss/revision/TSS/data/tmp.mainChr.novel.transcripts.map) $currnoHiss/revision/TSS/oldIDs.gtf > data/tmp.mainChr.novel.tagged.map
```

