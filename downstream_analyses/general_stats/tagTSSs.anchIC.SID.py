import sys

#script to add the anchIC, SID tags to the EBI mappings
#--usage
#python tagTSSs.anchIC.SID.py <(cat ${smT} ${umT}) data/${clsset}.TSS.bed data/tmp.mainChr.${clsset}.map oldIDs.gtf > data/${clsset}.TSS.anchIC-SIDtagged.bed
#${smT} ${umT} -> spliced and unspliced masterTables

table=sys.argv[1] 
tsstable=sys.argv[2] 
MAPPINGtable=sys.argv[3]
oldTable=sys.argv[4]

mT = open(table, "r")
tssT = open(tsstable,"r")
mapT= open(MAPPINGtable,"r")
oldT = open(oldTable,"r")

anchSIDmap={}
anchTargetmap={}
#TSSanchSIDmap={}
TSSanchSIDtargmap={}
ENSTanchMap={}

for line1 in oldT:
    line1=line1.strip()
    info = line1.split('"')
    oldID=info[1] + "_OLD" 
    info[5] = info[5].strip()
    anchSIDmap[oldID] = info[17] 
    #anchSIDmap[oldID] = ",".join(oldID + "." + sid for sid in info[17].split(","))
    
    #print(oldID,info[17],anchSIDmap[oldID])
    anchTargetmap[oldID] = info[5]

#print(anchSIDmap)

for line in mT:
    if "\ttranscript\t" in line:
        line=line.strip()
        info = line.split('"') 
        info[5] = info[5].strip()
        anchSIDmap[info[1]] = info[17] 
        #anchSIDmap[info[1]] = ",".join(info[1] + "." + sid for sid in info[17].split(","))
            
        anchTargetmap[info[1]] = info[5]
#print(anchSIDmap)

for mapping in mapT:
    mapInfo=mapping.split('\t') #if ("anchIC" in ID or "anchUC" in ID) and not "OLD" in ID:
#    if "anchIC" in mapInfo[3] or "anchUC" in mapInfo[3]: 
    ENSTanchMap[mapInfo[1]] = mapInfo[3]    #E1 = A1,A2 
    #print(mapInfo[1],mapInfo[3])   OK all OLDs
#    else:
#        ENSTanchMap[mapInfo[1]] = "NA(oldanchTM/UNMAPPED)"


for line in tssT:
        line=line.strip("\n")
        info=line.split("\t")
        ENSTs=info[8].split(",")
        #print(ENSTs)
        anchICs = []
        SIDs = []
        targets = []

        for ENST in ENSTs:
            ENST=ENST.split(".")[0]
            #get all anchICs hence SIDs 
            
            currentAnchICs=ENSTanchMap[ENST]
            
            anchICs.append(currentAnchICs)
            eachanchIC=currentAnchICs.split(",")      

            for anchIC in eachanchIC:
                anchIC = anchIC.strip()
            #    print(anchIC, anchIC in anchSIDmap)
                if "UNMAPPED" not in anchIC and anchIC != "":
                    if "OLDread" in anchIC:
                        anchIC = anchIC.split("_OLD")[0]

                    #print(anchIC, anchIC in anchSIDmap)
                    SIDs.append(anchSIDmap[anchIC])
                    targets.append(anchTargetmap[anchIC])
                    #print(SIDs,"current")        
                else:
                    #print(anchIC)   #488 empty, 1876 UNMAPPED
                    SIDs.append("NA(UNMAPPED)") #or empty
                    targets.append("NA(UNMAPPED)")
                    #print(anchIC) #unclean IDs
           
        #anchUC_000000002822  
        anchICs = list(set(x for x in ",".join(anchICs).split(",") if x)) 
        SIDs = list(set(x for x in ",".join(SIDs).split(",") if x))
        targets = list(set(x for x in ",".join(targets).split(",") if x))
        
        completeTSS = line #".".join([info[0], info[1], info[2], info[3], info[5]])
        TSSanchSIDtargmap[completeTSS] = ",".join(anchICs) + "\t" + ",".join(SIDs) + "\t" + (",".join(targets).strip() or "NA(noOverlappingTargets)")


for key, value in TSSanchSIDtargmap.items():
    print(f"{key}\t{value}")
