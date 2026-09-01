import sys

#script to add the anchIC, SID tags to the EBI mappings
#--usage
#python tagv47-CLSmap.anchIC.SID.py <(cat ${smT} ${umT}) <(grep -v geneID_v47 v47-CLS3mapping_status_anchTMs.txt) > TEST
#${smT} ${umT} -> spliced and unspliced masterTables

table=sys.argv[1] 
MAPPINGtable=sys.argv[2]
oldTable=sys.argv[3]

mT = open(table, "r")
mapT= open(MAPPINGtable,"r")
oldT = open(oldTable,"r")

anchSIDmap={}
ENSTanchSIDmap={}
anchTargetmap={}

for line1 in oldT:
    line1=line1.strip()
    info = line1.split('"')
    oldID=info[1] + "_OLD"
    info[5] = info[5].strip()
    anchSIDmap[oldID] = info[17]
    #anchSIDmap[oldID] = ",".join(oldID + "." + sid for sid in info[17].split(","))
    #print(oldID,info[17],anchSIDmap[oldID])
    anchTargetmap[oldID] = info[5]

for line in mT:
    if "\ttranscript\t" in line:
        line=line.strip()
        info = line.split('"')
        anchSIDmap[info[1]] = info[17]
        #anchSIDmap[info[1]] = ",".join(info[1] + "." + sid for sid in info[17].split(","))
        anchTargetmap[info[1]] = info[5]

for line in mapT:
        line=line.strip("\n")
        info=line.split("\t")
        ENSTs=info[1]
        #print(ENSTs) 

        #currentAnchICs=info[3]
        currentAnchICs = ",".join(sorted(set(info[3].split(","))))
        eachanchIC=currentAnchICs.split(",")    
        SIDs = []
        targets = []
        for anchIC in eachanchIC:
            anchIC = anchIC.strip()
            if "UNMAPPED" not in anchIC and anchIC != "":
            #if "anchTM" not in anchIC and "UNMAPPED" not in anchIC and not anchIC.endswith("_OLD") and anchIC != "":
                if "OLDread" in anchIC:
                    anchIC = anchIC.split("_OLD")[0]
                #print(anchIC)  #get all clean mappable IDs
                SIDs.append(anchSIDmap[anchIC])
                targets.append(anchTargetmap[anchIC])
                #print(SIDs,"current")
            else:
                SIDs.append("NA(oldID/UNMAPPED)") #or empty
                targets.append("NA(oldID/UNMAPPED)")
                #print(anchIC) #unclean IDs
       
        SIDs = list(set(x for x in ",".join(SIDs).split(",") if x))
        targets = list(set(x for x in ",".join(targets).split(",") if x))

        completemapInfo = line #".".join([info[0], info[1], info[2], info[3], info[5]])
        ENSTanchSIDmap[completemapInfo] = currentAnchICs + "\t" + ",".join(SIDs) + "\t" + (",".join(targets).strip() or "NA(noOverlappingTargets)")

for key, value in ENSTanchSIDmap.items():
    print(f"{key}\t{value}")
