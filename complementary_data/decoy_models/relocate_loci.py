import sys
import gzip
import random
import pickle
import os.path

intergenic = open(sys.argv[1])
annotation = gzip.open(sys.argv[2])

if ( not os.path.exists("output/intergenic.pickle")):
    print("Gathering intergenic area")
    intergenic_array = []
    for line in intergenic:
        chr = line.split()[0]
        start = int(line.split()[1])
        end = int(line.split()[2])

        intergenic_array.append((chr,start,end))
    intergenic.close()
    pickle.dump(intergenic_array, open("output/intergenic.pickle","wb"))
else:
    intergenic_array = pickle.load(open("output/intergenic.pickle","rb"))



if ( not os.path.exists("output/cls_loci.pickle")):
    print("Gathering Loci structure")
    annotation_dc = {}
    for line in annotation:
        line = line.decode("utf-8")
        locus = line.split("\"")[1]
        tm = line.split("\"")[3]

        if line.split()[2] == "transcript":
            if locus not in annotation_dc.keys():
                annotation_dc[locus] = {tm:{"length":abs(int(line.split()[4]) - int(line.split()[3])) + 1, "IC" : [], "starts" : [], "start":int(line.split()[3]), "end":int(line.split()[4])}}
            else:
                annotation_dc[locus][tm] = {"length":abs(int(line.split()[4]) - int(line.split()[3])) + 1, "IC" : [], "starts" : [], "start":int(line.split()[3]), "end":int(line.split()[4])}
            prev = int(line.split()[3])
        else:
            annotation_dc[locus][tm]["IC"].append(abs(int(line.split()[4]) - int(line.split()[3])))
            annotation_dc[locus][tm]["starts"].append(int(line.split()[3]) - prev)
            prev = int(line.split()[4])

    annotation.close()
    pickle.dump(annotation_dc, open("output/cls_loci.pickle","wb"))
else:
    annotation_dc = pickle.load(open("output/cls_loci.pickle","rb"))

print("Randomising the genomic coordinates")

intron_chain = {}
loci = list(annotation_dc.keys())
random.shuffle(loci)

for cls in loci:
    #Get only spliced TMs per locus
    tms = [i for i in annotation_dc[cls].keys() if len(annotation_dc[cls][i]["IC"]) > 1 ]

    if len(tms) == 0:
        continue
    
    length = max( [annotation_dc[cls][i]["end"] for i in tms] ) - min( [annotation_dc[cls][i]["start"] for i in tms]) + 1
    putative_regions_index = [ x for x in range(0,len(intergenic_array)) if (intergenic_array[x][2] - intergenic_array[x][1]) > length ]
    if ( len(putative_regions_index) == 0):
        continue

    intergenic_region = random.choice(putative_regions_index)

    new_tss = random.randint(intergenic_array[intergenic_region][1],intergenic_array[intergenic_region][2]-length)

    #Updating the intergenic space to avoid overlaps in relocation
    intergenic_array[intergenic_region] = (intergenic_array[intergenic_region][0], intergenic_array[intergenic_region][1], new_tss - 1)
    intergenic_array.append((intergenic_array[intergenic_region][0], new_tss+length+1, intergenic_array[intergenic_region][2]))
 
    tmp = [ (x, annotation_dc[cls][x]["start"]) for x in tms ] 
    sorted_tms = [ anchtm[0] for anchtm in sorted(tmp, key=lambda tup: tup[1]) ]

    leftmost = annotation_dc[cls][sorted_tms[0]]["start"]

    intron_chain[cls] = {}
    for chain in sorted_tms:
        tss = new_tss + (annotation_dc[cls][chain]["start"]-leftmost)

        ic = (tss,)
        for exonN in range(0,len(annotation_dc[cls][chain]["IC"])):
            ic += (ic[-1]+annotation_dc[cls][chain]["IC"][exonN],)
            if exonN != len(annotation_dc[cls][chain]["IC"])-1:
                ic += (ic[-1]+annotation_dc[cls][chain]["starts"][exonN+1],)

        intron_chain[cls][chain] = (intergenic_array[intergenic_region][0], ic)
    
pickle.dump(intron_chain, open("../output/random.pickle","wb"))


print("Storing to file...")
random_replicates = open("output/random_replicates_locirelocation.gtf", "w")

for locus in intron_chain.keys():
    for tm,coordinates in intron_chain[locus].items():
        random_replicates.write(coordinates[0]+"\trandom\ttranscript\t"+str(coordinates[1][0])+"\t"+str(coordinates[1][-1])+"\t0\t.\t.\tgene_id \""+locus+"\"; transcript_id \""+tm+"\";\n")

        for i in range(0,len(coordinates[1]),2):
            random_replicates.write(coordinates[0]+"\trandom\texon\t"+str(coordinates[1][i])+"\t"+str(coordinates[1][i+1])+"\t0\t.\t.\tgene_id \""+locus+"\"; transcript_id \""+tm+"\";\n")

random_replicates.close()
