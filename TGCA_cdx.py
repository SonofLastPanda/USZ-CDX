#TCGA version
import gzip
import re
import csv
import sys
csv.field_size_limit(sys.maxsize)
cdx_genes_file = "cdx_genes.txt"
c_file=open(cdx_genes_file)
cdx_genes=[]
NotExistTCGA=open("nexistgenes_barcode.txt",'w')
#Bar_File=open("Barcodes.txt",'w')
def takeString(word):
    fil_word=""
    for a in word:
        if(a!=" "):
            fil_word=fil_word+a
        else:
            break
    return fil_word

def IsCosmic(word1,word2,word3):
    b=False
    if(cosmic_dict.get(word1,'None')!='None'):
        for list in cosmic_dict[word1]:
            if(list[0]==word2 and list[1]==word3):
                b=True
    return b


def remENS(word):
    index=word.find("ENS")
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


for line in c_file:
    line=line.strip('\r\n')
    if( " " in line):
        line=takeString(line)
    cdx_genes.append(line)

c_file.close()
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
tsg_list=[]

for line in census:
    line=line.strip('\r\n')
    line=line.split('\t')
    gene=line[0]
    role=line[14]
    if(Role(role)=='TSG'):
        tsg_list.append(gene)


ignored={}
out_data={}
genes={}
var_type={}
f=gzip.open("TCGA-All.maf.gz")
#f=gzip.open("tcga_dummy.maf.gz")
i=1


for line in csv.reader(f, dialect="excel-tab"):
    if(i==1):
        i+=1
        continue
    else:
        app=line[0]
        if(app.find("raw_maf")!=-1):
            ind=app.find(':')
            app=app[ind+1:]
        if(app != "Unknown"):
        #Bar_File.write('{}\n'.format(line[15]))
            if(ignored.get(line[15],'None')== 'None'):
                if(app not in cdx_genes):
                    if(genes.get(line[15],'None') == 'None'):
                        list_app=[app]
                        list_var=[line[8]]
                        genes[line[15]]=list_app
                        var_type[line[15]]=list_var
                    else:
                        genes[line[15]].append(app)
                        var_type[line[15]].append(line[8])
                else:
                    ignored.update({line[15]:line[15]})
                    genes.pop(line[15],'None')
                    var_type.pop(line[15],'None')
                    '''
                    else:
                        if(genes.get(line[15],'None') == 'None'):
                            list_app=[app]
                            list_var=[line[8]]
                            genes[line[15]]=list_app
                            var_type[line[15]]=list_var
                        else:
                            genes[line[15]].append(app)
                            var_type[line[15]].append(line[8])
                    '''
NotExistTCGA.write('{a}\t{b}\t{c}\t{d}\t\n'.format(a="Tumor_Sample_Barcode",b="Mutated Genes",c="Var_Type",d="Vus_number"))
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
    NotExistTCGA.write('{}\n'.format(w_line))
