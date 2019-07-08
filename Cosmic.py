#A script for identifying mutations in PCAWG data that has no mutated gene that are among cdx genes
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
                    ignored.update({line[12]:line[12]})
                    #out_data.pop(line[12],'None')
                    genes.pop(line[12],'None')
                    var_type.pop(line[12],'None')

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
