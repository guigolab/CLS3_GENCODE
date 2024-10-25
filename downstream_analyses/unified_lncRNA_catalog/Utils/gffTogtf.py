#!/usr/bin/python

__author__="BU"


"""This script converts gff to gtf"""

# Usage python gff3Togtf.py -i annotation.gff3 -o output.gtf 

from os import sys
import getopt
import heapq
from optparse import OptionParser

# ----------------
## OPTION PARSING
# ----------------

parser = OptionParser()
parser.add_option("-i",  help="Annotation in gff")
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

def gffTogtf():
   
    gtf={}
    tx2gene={}
    reg2chr={}

    gff=open(infile, "r")
    for line in gff:
        if not line.startswith('#'):
            split=line.split("\t")
            #print split
            if split[2]=="region":
                reg=split[0]
                #print reg
                spliter=split[8].split(";")
                #print spliter
                for item in spliter:           
                    if item.startswith('chromosome'): 
                        chrT=item.split('=')[1].split('\n')[0]
                        chrN='chr'+chrT
                        #print chrN
                          
                        reg2chr[reg]=chrN
             
            if split[2]=="gene":
                #print 'Gene', split
                regs=split[0]
                if regs in reg2chr:
                    chrN=reg2chr[regs]   
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
                    if item.startswith('ID'): 
                        ids=item.split('=')[1].split('\n')[0]
                        
                    if item.startswith('Dbxref'): 
                        core=item.split('=')[1].split(',')[0].split(':')[1]
                        #print core
                        gene_id='GeneID:'+core+"_"+ids
                        #print gene_id
  
                    if item.startswith('gene_biotype'): 
                        gene_biotype=item.split('=')[1].split('\n')[0]
                
                        record=(chrN+"\t"+add+"\t"+typ+"\t"+start+"\t"+stop+"\t"+dot1+"\t"+strand+"\t"+dot2+"\t"+"gene_id"+" "+'"'+gene_id+'"'"; "+"transcript_id"+" "+'"'+gene_id+'"'"; "+"gene_biotype"+" "+'"'+gene_biotype+'"'";")
           
                        if gene_id not in gtf:
                            gtf[gene_id]=[record]
                        else:
                            gtf[gene_id].append(record)  


            if split[2]=="transcript":
                #print 'TX', split
                regs=split[0]
                if regs in reg2chr:
                    chrN=reg2chr[regs]
                add=split[1]
                typ=split[2]
                start=split[3]
                stop=split[4]
                dot1=split[5]
                strand=split[6]
                dot2=split[7]
                transcript_biotype='NA'
                
                spliter=split[8].split(";")
                for item in spliter:
                    if item.startswith('ID'): 
                        txids=item.split('=')[1].split('\n')[0]
                        #print transcript_id
                    if item.startswith('Parent'): 
                        ids=item.split('=')[1].split(',')[0]
                        #print gene_id
                    if item.startswith('Dbxref'): 
                        core=item.split('=')[1].split(',')[0].split(':')[1]
                        #print core
                        gene_id='GeneID:'+core+"_"+ids
                        #print gene_id
                    if item.startswith('transcript_id'): 
                        txcore=item.split('=')[1].split('\n')[0]
                        #print core
                        transcript_id=txcore+"_"+txids
                        #print transcript_id

               
                        record=(chrN+"\t"+add+"\t"+typ+"\t"+start+"\t"+stop+"\t"+dot1+"\t"+strand+"\t"+dot2+"\t"+"gene_id"+" "+'"'+gene_id+'"'"; "+"transcript_id"+" "+'"'+transcript_id+'"'"; "+"transcript_biotype"+" "+'"'+transcript_biotype+'"'";")
           
                        if gene_id not in gtf:
                            gtf[gene_id]=[record]
                        else:
                            gtf[gene_id].append(record)  

                        tx2gene[transcript_id]=gene_id

            if split[2]=="lnc_RNA":
                #print 'TX', split
                regs=split[0]
                if regs in reg2chr:
                    chrN=reg2chr[regs]
                add=split[1]
                typ=split[2]
                start=split[3]
                stop=split[4]
                dot1=split[5]
                strand=split[6]
                dot2=split[7]
                transcript_biotype='lncRNA'
                
                spliter=split[8].split(";")
                for item in spliter:
                    if item.startswith('ID'): 
                        txids=item.split('=')[1].split('\n')[0]
                        #print transcript_id
                    if item.startswith('Parent'): 
                        ids=item.split('=')[1].split(',')[0]
                        #print gene_id
                    if item.startswith('Dbxref'): 
                        core=item.split('=')[1].split(',')[0].split(':')[1]
                        #print core
                        gene_id='GeneID:'+core+"_"+ids
                        #print gene_id
                    if item.startswith('transcript_id'): 
                        txcore=item.split('=')[1].split('\n')[0]
                        #print core
                        transcript_id=txcore+"_"+txids
                        #print transcript_id
               
                        record=(chrN+"\t"+add+"\t"+typ+"\t"+start+"\t"+stop+"\t"+dot1+"\t"+strand+"\t"+dot2+"\t"+"gene_id"+" "+'"'+gene_id+'"'"; "+"transcript_id"+" "+'"'+transcript_id+'"'"; "+"transcript_biotype"+" "+'"'+transcript_biotype+'"'";")
           
                        if gene_id not in gtf:
                            gtf[gene_id]=[record]
                        else:
                            gtf[gene_id].append(record)

                        tx2gene[transcript_id]=gene_id     

            if split[2]=="exon":
                regs=split[0]
                if regs in reg2chr:
                    chrN=reg2chr[regs]   
                add=split[1]
                typ=split[2]
                start=split[3]
                stop=split[4]
                dot1=split[5]
                strand=split[6]
                dot2=split[7] 

                spliter=split[8].split(";")
                for item in spliter:
                    if item.startswith('ID'): 
                        exon_id=item.split('=')[1]
                        #print exon_id
                    if item.startswith('Parent'): 
                        txids=item.split('=')[1].split('\n')[0]
                    if item.startswith('transcript_id'): 
                        txcore=item.split('=')[1].split('\n')[0]
                        #print core
                        transcript_id=txcore+"_"+txids
                        if transcript_id in tx2gene:
                            gene_id=tx2gene[transcript_id]
        
               
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

gffTogtf()
