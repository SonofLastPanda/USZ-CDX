#aim is to find the number of occurences oncogene, tumor suppressor gene and fusion genes in each sample

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

width=0.25
non_exist_bar=open("nexistgenes_barcode.txt")
census=open("cancer_gene_census.txt")
out=open("Sample_Gene_Info.tsv",'w')
bar_gene={}
onco_gene={}
tsg={}
fusion={}
null={}


def Role(word):
    if("oncogene" in word):
        return "oncogene"
    elif("TSG" in word):
        return "TSG"
    else:
        return "fusion"


for line in non_exist_bar:
    line=line.strip('\n')
    line=line.split('\t')
    bar=line[0]
    onco_gene[bar]=[]
    tsg[bar]=[]
    fusion[bar]=[]
    null[bar]=[]
    gene_list=line[1].split(',')
    bar_gene[bar]=gene_list[:-1]

for line in census:
    line=line.strip('\r\n')
    line=line.split('\t')
    gene=line[0]
    for key in bar_gene.keys():
        if(gene in bar_gene[key]):
            role=Role(line[14])
            if(role=="oncogene"):
                onco_gene[key].append(gene)
            elif(role=="TSG"):
                tsg[key].append(gene)
            elif(role=="fusion"):
                fusion[key].append(gene)
            else:
                null[key].append(gene)

out.write('{a}\t{b}\t{c}\t{d}\t{e}\t\n'.format(a="Tumor_Sample_Barcode",b="Oncogene_Number",c="TSG_Number",d="Fusion_number",e="Null"))

for key in bar_gene.keys():
    out.write('{a}\t{b}\t{c}\t{d}\t{e}\t\n'.format(a=key,b=str(len(onco_gene[key])),c=str(len(tsg[key])),d=str(len(fusion[key])),e=str(len(null[key]))))

out.close()

list_onco=[]
list_tsg=[]
list_fus=[]
#list_null=[]
for key in bar_gene.keys():
    list_onco.append(len(onco_gene[key]))
    list_tsg.append(len(tsg[key]))
    list_fus.append(len(fusion[key]))
    #list_null.append(len(null[key]))

r1 = np.arange(len(list_onco))
r2 = [x + width for x in r1]
r3 = [x + width for x in r2]
#r4 = [x + width for x in r3]



plt.bar(r1, list_onco, color='#EC182D', width=width, edgecolor='white', label='Oncogene')
plt.bar(r2, list_tsg, color='#183BEC', width=width, edgecolor='white', label='TSG')
plt.bar(r3, list_fus, color='#216237', width=width, edgecolor='white', label='Fusion')
#plt.bar(r4,list_null,color='#AEC7AF', width=width, edgecolor='white', label='NULL')


plt.xlabel('Tumor_Sample_Barcode', fontweight='bold')
plt.ylabel('Number of Occurence',fontweight='bold')
plt.xticks([r + width for r in range(len(list_onco))], range(1,25),rotation=90)

plt.legend()
plt.show()
