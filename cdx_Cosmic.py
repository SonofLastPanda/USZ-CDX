#Finding the locations of the cdx mutations in Cosmic file.
#Cosmic file is available in:https://cancer.sanger.ac.uk/cosmic/download

import gzip
import csv

Cosmic_f=gzip.open("CosmicGenomeScreenMutantExport.tsv.gz")
cdx_genes=open("cdx_genes.txt")
cdx_cosmic_file=open("cdx_cosmic_common.tsv",'w')
cdx_cosmic_file.write('{a}\t{b}\t{c}\t\n'.format(a="Chromosome",b="Start_pos",c="End_pos"))
cdx=[]
def takeString(word):
    word.split()
    return word[0]

def remENS(word):
    index=word.find("_")
    word=word[0:index]
    return word

def ChromosomeInfo(word):
    DoubleDotInd=word.find(":")
    Chromosome=word[0:DoubleDotInd]
    Dash=word[DoubleDotInd:].find("-")
    Start_pos=word[DoubleDotInd:][1:Dash]
    End_pos=word[DoubleDotInd:][Dash+1:]
    list=[]
    list.append(Chromosome)
    list.append(Start_pos)
    list.append(End_pos)
    return list
for line in cdx_genes:
    line=line.strip('\n')
    if( " " in line):
        line=takeString(line)
    cdx.append(line)

for line in csv.reader(Cosmic_f, dialect="excel-tab"):
    if(line[0].find("_ENS") != -1):
        app=remENS(line[0])
    else:
        app=line[0]
    if(app in cdx):
        C_info=ChromosomeInfo(line[23])
        cdx_cosmic_file.write('{a}\t{b}\t{c}\t\n'.format(a=C_info[0],b=C_info[1],c=C_info[2]))
