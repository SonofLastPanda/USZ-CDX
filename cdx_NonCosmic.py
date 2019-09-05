#Script to extend the list of patients with the ones only have VUS mutations.
import gzip
import re
import csv

cdx_genes_file = "cdx_genes.txt"
non_cdx_PCAWG=open("nexistgenes_barcode.txt")
Cosmic_cdx=open("cdx_cosmic_common.tsv")
#Cosmic_cdx=open("cdx_cosmic_common_dummy.tsv")
c_file=open(cdx_genes_file)
cdx_genes=[]
NotExistPCAWG=open("nexistgenes_barcode.txt",'a')
#NotExistPCAWG=open("nexistgenes_barcode_dummy.txt",'w')
#Bar_File=open("Barcodes.txt",'w')
def takeString(word):
    word=word.split()
    return word[0]


def IsCosmic(word1,word2,word3):
    b=False
    for list in cosmic_cdx:
        if(list[0]==word1):
            if(list[1]==word2 and list[2]==word3):
                b=True
    return b


for line in c_file:  #cdx gene list
    if( " " in line):
        line=takeString(line)
    line=line.strip('\n')
    cdx_genes.append(line)



cosmic_cdx=[]
a=0
for line in Cosmic_cdx:
    if(a!=0):
        line=line.strip('\n')
        line=line.split('\t')
        list=[line[0],line[1],line[2]]
        cosmic_cdx.append(list)
    a+=1
c_file.close()

non_cdx=[]
a=0
for row in non_cdx_PCAWG:
    if(a!=0):
        row=row.split()
        non_cdx.append(row[0])
    a+=1


ignored={}
out_data={}
genes={}
var_type={}
f=gzip.open("October_2016_whitelist_2583.snv_mnv_indel.maf.gzip")
#f=gzip.open("bbb.maf.gz")
i=1


for line in csv.reader(f, dialect="excel-tab"):
    if(i==1):
        i+=1
        continue
    else:
        app=line[0]
        if(app != "Unknown"):
            if(line[12] not in non_cdx_PCAWG):
                #Bar_File.write('{}\n'.format(line[12]))
                if(ignored.get(line[12],'None')== 'None'):
                    if(app in cdx_genes):
                        if(IsCosmic(line[1],line[2],line[3])==False):
                            if(genes.get(line[12],'None') == 'None'):
                                list_app=[app]
                                list_var=[line[5]]
                                genes[line[12]]=list_app
                                var_type[line[12]]=list_var
                            else:
                                genes[line[12]].append(app)
                                var_type[line[12]].append(line[5])
                        else:
                            ignored.update({line[12]:line[12]})
                            #out_data.pop(line[12],'None')
                            genes.pop(line[12],'None')
                            var_type.pop(line[12],'None')
                    else:
                        ignored.update({line[12]:line[12]})
                        genes.pop(line[12],'None')
                        var_type.pop(line[12],'None')
#NotExistPCAWG.write('{a}\t{b}\t{c}\t{d}\t\n'.format(a="Tumor_Sample_Barcode",b="Mutated Genes",c="Var_Type",d="Vus_number"))
w_line=""
for key in genes.keys():
    w_line=key+"\t"
    for value in genes[key]:
        w_line+=value+","
    w_line+="\t"
    for val in var_type[key]:
        w_line+=val+","
    w_line+="\t"
    w_line+=str(len(genes[key]))
    w_line=w_line.decode('string_escape')
    NotExistPCAWG.write('{}\n'.format(w_line))
