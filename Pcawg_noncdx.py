#PCAWG version for finding the patients with only non-cdx mutations.
import gzip
import re
import csv

cdx_genes_file = "cdx_genes.txt"
c_file=open(cdx_genes_file)
cdx_genes=[]
NotExistPCAWG=open("nexistgenes_barcode.txt",'w')
#Bar_File=open("Barcodes.txt",'w')
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
'''
cosmic_dict={}
Cosmic_f=gzip.open("CosmicGenomeScreenMutantExport.tsv.gz")

for l in csv.reader(Cosmic_f, dialect="excel-tab"):
    C_info=ChromosomeInfo(l[23])
    if(cosmic_dict.get(C_info[0],'None')=='None'):
        cosmic_dict[C_info[0]]=[[C_info[1],C_info[2]]]
    else:
        cosmic_dict[C_info[0]].append([C_info[1],C_info[2]])
census=open("cancer_gene_census.txt")

def Role(word):
    if("oncogene" in word):
        return "oncogene"
    elif("TSG" in word):
        return "TSG"
    else:
        return "fusion"
#tsg_list=[]

for line in census:
    line=line.strip('\r\n')
    line=line.split('\t')
    gene=line[0]
    role=line[14]
    if(Role(role)=='TSG'):
        tsg_list.append(gene)
'''

ignored={}
out_data={}
genes={}
var_type={}
f=gzip.open("October_2016_whitelist_2583.snv_mnv_indel.maf.gzip")
#f=gzip.open("pcawg_dummy.maf.gz")
i=1


for line in csv.reader(f, dialect="excel-tab"):
    if(i==1):
        i+=1
        continue
    else:
        app=line[0]
        if(app != "Unknown"):
        #Bar_File.write('{}\n'.format(line[12]))
            if(ignored.get(line[12],'None')== 'None'):
                if(app not in cdx_genes):
                    #out_data.update({line[12]:line[12]})
                    #print(line[12])
                    if(genes.get(line[12],'None') == 'None'):
                        list_app=[app]
                        list_var=[line[5]]
                        genes[line[12]]=list_app
                        var_type[line[12]]=list_var
                    else:
                        genes[line[12]].append(app)
                        var_type[line[12]].append(line[5])
                else:
                    #if(IsCosmic(line[1],line[2],line[3])==True or (app in tsg_list and (line[5].find("Frame_Shift") != -1 or line[5].find("Nonsense")!= -1))):
                    ignored.update({line[12]:line[12]})
                    #out_data.pop(line[12],'None')
                    genes.pop(line[12],'None')
                    var_type.pop(line[12],'None')
                    '''
                    else:
                        if(genes.get(line[12],'None') == 'None'):
                            list_app=[app]
                            list_var=[line[5]]
                            genes[line[12]]=list_app
                            var_type[line[12]]=list_var
                        else:
                            genes[line[12]].append(app)
                            var_type[line[12]].append(line[5])
                    '''

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
