#PCAWG version for finding the patients with only non-cdx mutations.
import gzip
import re
import csv

cdx_genes_file = "cdx_genes.txt" #the file that contains cdx genes
c_file=open(cdx_genes_file)
cdx_genes=[]
NotExistPCAWG=open("nexistgenes_barcode.txt",'w') #output file to store the patients and their information.

def takeString(word):
    fil_word=""
    for a in word:
        if(a!=" "):
            fil_word=fil_word+a
        else:
            break
    return fil_word

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


def remENS(word):
    index=word.find("ENS")
    word=word[0:index]
    return word

for line in c_file:
    line=line.strip('/n')
    if( " " in line):
        line=takeString(line)
    cdx_genes.append(line[0:-1])


c_file.close()
def IsCosmic(word1,word2,word3):
    b=False
    if(cosmic_dict.get(word1,'None')!='None'):
        for list in cosmic_dict[word1]:
            if(list[0]==word2 and list[1]==word3):
                b=True
    return b


ignored={} #dictionary to store ignored barcodes
genes={} #dictionary for storing the genes
var_type={} #dictionary for storing variation types 
f=gzip.open("October_2016_whitelist_2583.snv_mnv_indel.maf.gzip") #PCAWG file
i=1


for line in csv.reader(f, dialect="excel-tab"):
    if(i==1):
        i+=1
        continue
    else:
        app=line[0] #the gene in the line
        if(app != "Unknown"):
            if(ignored.get(line[12],'None')== 'None'): #if the patient barcode is not ignored yet
                if(app not in cdx_genes): #if the gene is not a cdx gene
                    if(genes.get(line[12],'None') == 'None'): #if the info. of the patient is not created yet, create.
                        list_app=[app]
                        list_var=[line[5]]
                        genes[line[12]]=list_app
                        var_type[line[12]]=list_var
                    else: #if it is already created, append the existing one.
                        genes[line[12]].append(app)
                        var_type[line[12]].append(line[5])
                else: #if the gene is a cdx gene, ignore the patient barcode and remove the info. from the other dictionaries.
                    ignored.update({line[12]:line[12]})
                    genes.pop(line[12],'None')
                    var_type.pop(line[12],'None')
                   

i_list=[]
NotExistPCAWG.write('{a}\t{b}\t{c}\t{d}\t\n'.format(a="Tumor_Sample_Barcode",b="Mutated Genes",c="Var_Type",d="Vus_number"))
w_line=""
for key in genes.keys():
    w_line=key+"\t"
    for value in genes[key]:
        w_line+=value+","
    w_line+="\t"
    for val in var_type[key]:
        w_line+=val+","
    w_line+="\t"
    w_line+=str("0")
    w_line=w_line.decode('string_escape')
    NotExistPCAWG.write('{}\n'.format(w_line))
