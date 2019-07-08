#Aim: is to detect PCAWG mutation samples with cdx genes that are not in Cosmic
import gzip
import csv

Cosmic_f=gzip.open("CosmicGenomeScreenMutantExport.tsv.gz")
cdx_genes=open("cdx_genes.txt")
cdx_cosmic_file=open("cdx_cosmic_common.tsv",'w')
cdx_cosmic_file.write('{a}\t{b}\t{c}\t\n'.format(a="Chromosome",b="Start_pos",c="End_pos"))
cdx=[]
def takeString(word):
    word.split()
    print(word)
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
    print("#1 " + line[0] + " " + line[23])
    if(line[0].find("_ENS") != -1):
        app=remENS(line[0])
        print("#2 " + line[0] + " " + app)
    else:
        app=line[0]
        print("#3 " + line[0] + " " + app)
    if(app in cdx):
        C_info=ChromosomeInfo(line[23])
        cdx_cosmic_file.write('{a}\t{b}\t{c}\t\n'.format(a=C_info[0],b=C_info[1],c=C_info[2]))
        print("#4" + app + " " + line[23] + " " + str(C_info))
    else:
        print("#5" + app + " " + line[23])
