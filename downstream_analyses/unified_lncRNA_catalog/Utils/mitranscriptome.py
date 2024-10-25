#!/usr/bin/python

__author__="BU"


"""This script converts mitranscriptome to gencode-like gtf"""

# Usage python mitranscriptome.py -i mitranscriptome.gtf -o output.gtf 

from os import sys
import getopt
import heapq
from optparse import OptionParser

# ----------------
## OPTION PARSING
# ----------------

parser = OptionParser()
parser.add_option("-i",  help="MiTranscriptome annotation in gtf")
parser.add_option("-o",  help='Output gtf file')

options, args = parser.parse_args()


# Store input and output file

infile=''
outfile=''

# Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"i:o:")

# Loop over arguments and assign input to infile and output to outfile 
# o == option
# a == argument passed to the o

#if myopts == []:
#    print("Usage: %s -i input -o output" % sys.argv[0])

for o, a in myopts:
    if o == '-i':
        infile=a
    if o == '-o':
        outfile=a

print ("Input file: %s \nOutput file: %s" % (infile,outfile))

#-------------------------------------------------------------------------------------------------
# Defining Functions
#-------------------------------------------------------------------------------------------------

def mitrans():
   
    gtf={}
    tx2gene={}

    gff=open(infile, "r")
    for line in gff:
        if not line.startswith('#'):
            split=line.split("\t")
            if split[2]=="transcript":
                #print 'TX', split
                chrN=split[0]
                add=split[1]
                typ=split[2]
                start=split[3]
                stop=split[4]
                dot1=split[5]
                strand=split[6]
                dot2=split[7]
                
                spliter=split[8].split(";")
                #print spliter

                for item in spliter:
                    if item.startswith('tcat'): 
                        transcript_biotype=item.split('"')[1]
                        #print transcript_biotype
                    if item.startswith(' gene_id'): 
                        gid=item.split('"')[1]
                        gene_id=gid+"_"+gid
                    if item.startswith(' transcript_id'): 
                        tid=item.split('"')[1]
                        transcript_id=tid+"_"+tid

               
                        record=(chrN+"\t"+add+"\t"+typ+"\t"+start+"\t"+stop+"\t"+dot1+"\t"+strand+"\t"+dot2+"\t"+"gene_id"+" "+'"'+gene_id+'"'"; "+"transcript_id"+" "+'"'+transcript_id+'"'"; "+"transcript_biotype"+" "+'"'+transcript_biotype+'"'";")
           
                        if gene_id not in gtf:
                            gtf[gene_id]=[record]
                        else:
                            gtf[gene_id].append(record)  

            if split[2]=="exon":
                chrN=split[0]
                add=split[1]
                typ=split[2]
                start=split[3]
                stop=split[4]
                dot1=split[5]
                strand=split[6]
                dot2=split[7] 

                spliter=split[8].split(";")
                for item in spliter:
                    if item.startswith('exon_number'): 
                        enb=item.split('"')[1]
                        #print enb
                    if item.startswith(' gene_id'): 
                        gid=item.split('"')[1]
                        gene_id=gid+"_"+gid
                    if item.startswith(' transcript_id'): 
                        tid=item.split('"')[1]
                        transcript_id=tid+"_"+tid
                        exon_id=transcript_id+"."+enb

        
               
                        record=(chrN+"\t"+add+"\t"+typ+"\t"+start+"\t"+stop+"\t"+dot1+"\t"+strand+"\t"+dot2+"\t"+"gene_id"+" "+'"'+gene_id+'"'"; "+"transcript_id"+" "+'"'+transcript_id+'"'"; "+"exon_id"+" "+'"'+exon_id+'"'";")
           
                        if gene_id not in gtf:
                            gtf[gene_id]=[record]
                        else:
                            gtf[gene_id].append(record)  


    f=open(outfile, "w")

    for key in gtf:
        for v in gtf[key]:
            f.write('%s\n' % v)
    f.close() 

mitrans()
