import matplotlib
import matplotlib.pyplot as plt
import numpy as np

width=0.25
non_exist_bar=open("nexistgenes_barcode.txt")
census=open("cancer_gene_census.csv")
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

i=1
for line in non_exist_bar:
    if(i==1):
        i+=1
    else:
        line=line.strip('\n')
        line=line.split('\t')
        bar=line[0]
        onco_gene[bar]=[]
        tsg[bar]=[]
        fusion[bar]=[]
        null[bar]=[]
        gene_list=line[1].split(',')
        bar_gene[bar]=gene_list[:-1]
i=1
for line in census:
    if(i==1):
        i+=1
    else:
        line=line.strip('\n')
        line=line.split(',')
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

out.write('{a}\t{b}\t{c}\t{d}\t{e}\t{f}\t{g}\t{h}\t{i}\t\n'.format(a="Tumor_Sample_Barcode",b="Oncogene",c="Oncogene_Number",d="TSG",e="TSG_Number",f="Fusion",g="Fusion_number",h="Null",i="Null_number"))
#w=""
for key in bar_gene.keys():
    out.write(key+'\t')
    for gene in onco_gene[key]:
        out.write(gene+',')
    out.write('\t')
    out.write(str(len(onco_gene[key])))
    out.write('\t')
    for gene in tsg[key]:
        out.write(gene+',')
    out.write('\t')
    out.write(str(len(tsg[key])))
    out.write('\t')
    for gene in fusion[key]:
        out.write(gene+',')
    out.write('\t')
    out.write(str(len(fusion[key])))
    out.write('\t')
    for gene in null[key]:
        out.write(gene+',')
    out.write('\t')
    out.write(str(len(null[key])))
    out.write('\n')
    '''
    w+=(key+'\t')
    for gene in onco_gene[key]:
        w+=(gene+',')
    w=w[:-1]
    w+='\t'
    w+=str(len(onco_gene[key]))
    for gene in tsg[key]:
        w+=(gene+',')
    w=w[:-1]
    w+='\t'
    w+=(str(len(tsg[key])))
    for gene in fusion[key]:
        w+=(gene+',')
    w=w[:-1]
    w+='\t'
    w+=(str(len(fusion[key])))
    for gene in null[key]:
        w+=(gene+',')
    w=w[:-1]
    w+='\t'
    w+=(str(len(null[key])))
    w+='\n'
    out.write(w)
    '''
#    out.write('{a}\t{b}\t{c}\t{d}\t{e}\t{f}\t{g}\t{h}\t{i}\t\n'.format(a=key,b=onco_gene[key],c=str(len(onco_gene[key])),d=tsg[key],e=str(len(tsg[key])),f=fusion[key],g=str(len(fusion[key])),h=null[key],i=str(len(null[key]))))
'''
out.close()

list_onco=[]
list_tsg=[]
list_fus=[]
list_null=[]
i=1
for key in bar_gene.keys():
    if(i<20):
        i+=1
        list_onco.append(len(onco_gene[key]))
        list_tsg.append(len(tsg[key]))
        list_fus.append(len(fusion[key]))
        list_null.append(len(null[key]))
    else:
        break

r1 = np.arange(len(list_onco))
r2 = [x + width for x in r1]
r3 = [x + width for x in r2]
r4 = [x + width for x in r3]



plt.bar(r1, list_onco, color='#EC182D', width=width, edgecolor='white', label='Oncogene')
plt.bar(r2, list_tsg, color='#183BEC', width=width, edgecolor='white', label='TSG')
plt.bar(r3, list_fus, color='#216237', width=width, edgecolor='white', label='Fusion')
plt.bar(r4,list_null,color='#AEC7AF', width=width, edgecolor='white', label='NULL')


plt.xlabel('Tumor_Sample_Barcode', fontweight='bold')
plt.ylabel('Number of Occurence',fontweight='bold')
plt.xticks([r + width for r in range(len(list_onco))], range(1,25),rotation=90)

plt.legend()
plt.show()
'''
