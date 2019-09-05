import os
import csv
import urllib2
import matplotlib.pyplot as plt
import numpy as np
import operator

class Gene:
    name=[]
    freq=""

    def __init__(self, name,freq):
        self.name = [name]
        self.freq=freq


def Gene_Freq(path):
    gene_freq={}
    for dir in dir_list:
        id=open("/Users/erkin/Projects/TCGA/Barcodes/"+dir+"/Sample_Kegg_Pathway.txt")
        for line in id:
            line=line.strip('\n')
            line=line.split('\t')
            if(path==line[0]):
                gene_l=line[3].split(',')[:-1]
                for g in gene_l:
                    if(gene_freq.get(g,'None')=='None'):
                        gene_freq[g]=1
                    else:
                        gene_freq[g]+=1
        id.seek(0,0)
        id.close()
    return gene_freq

clinical=open("Barcode_Project.txt")
project_dict={}
for line in clinical:
    line=line.strip('\n')
    line=line.split('\t')
    b=line[0]
    p=line[1]
    project_dict[b]=p


id_pathway=open("pathway.txt")
non_exist_bar=open("nexistgenes_barcode.txt")
all_gene={}
for line in non_exist_bar:
    line=line.strip('\n')
    line=line.split('\t')
    gene_list=line[1].split(',')
    for gene in gene_list[:-1]:
        if(all_gene.get(gene,'None')=='None'):
            all_gene[gene]=1
        else:
            all_gene[gene]+=1


url="http://rest.kegg.jp/list/"
dir_list=os.listdir("/Users/erkin/Projects/TCGA/Barcodes/")
dir_list.remove(".DS_Store")

path_name={}

name=open("Pathway Name.txt")
for line in name:
    line=line.strip('\n')
    line=line.split('\t')
    path_name[line[0]]=line[1]

path_map={}
gene_map={}
name_map={}
entrez=open("uniprot-yourlist_M201909028471C63D39733769F8E060B506551E123F05953-filtered-rev--.tab")

i=1
for line in csv.reader(entrez, dialect="excel-tab"):
    if(i==1):
        i+=1
    else:
        entrez_id=line[9]
        mapped_Uni=line[8]
        name=line[0]
        entry=line[1]
        #gene_map[mapped_Uni]=entrez_id
        #name_map[name]=entry
        gene_map[entrez_id]=mapped_Uni
        name_map[entry]=name

dir_list=os.listdir("/Users/erkin/Projects/TCGA/Barcodes/")
dir_list.remove(".DS_Store")

'''
for dir in dir_list:
    os.chdir("/Users/erkin/Projects/TCGA/Barcodes/"+dir)
    kegg=open("Kegg.txt",'w')
    file=open(dir+'.txt')
    for line in file:
        line=line.strip('\n')
        gene=line
        if(name_map.get(gene,'None')!='None'):
            print('a')
            uni=name_map[gene]
            if(gene_map.get(uni,'None')!= 'None'):
                print('b')
                entrez=gene_map[uni]
                kegg.write(entrez+'\n')
    kegg.close()
os.chdir("/Users/erkin/Projects/TCGA/")
'''
count=0

for dir in dir_list:
    count+=1
    id=open("/Users/erkin/Projects/TCGA/Barcodes/"+dir+"/Kegg.txt")
    id_pathway=open("/Users/erkin/Projects/TCGA/pathway.txt")
    out=open("/Users/erkin/Projects/TCGA/Barcodes/"+dir+"/Sample_Kegg_Pathway.txt",'w')
    id_list=[]
    path_dict={}
    for line in id:
        line=line.strip('\r\n')
        id_list.append(line)

    for line in id_pathway:
        line=line.strip('\n')
        line=line.split('\t')
        line_id=line[0][4:]
        if(line_id in id_list):
            if(path_dict.get(line[1],'None')=='None'):
                g=Gene(line_id,1)
                path_dict[line[1]]=g
            else:
                path_dict[line[1]].freq+=1
                path_dict[line[1]].name.append(line_id)
    id.close()
    id_pathway.seek(0,0)
    id_pathway.close()
    out.write("Pathway ID"+'\t'+"Pathway Name"+'\t'+"Occurrence in Sample"+'\t'+"Genes"+'\n')
    for key in path_dict.keys():
        #html=urllib2.urlopen(url+key)
        #h=html.read()
        #h=h.strip('\n')
        #h=h.split('\t')
        #path_name[h[0]]=h[1]
        out.write(key+'\t'+path_name[key]+'\t'+str(path_dict[key].freq)+'\t')
        p_list=[]
        for i in path_dict[key].name:
            if(gene_map[i] not in p_list):
                out.write(gene_map[i]+',')
                p_list.append(gene_map[i])
        out.write('\n')


path_all={}
path_n={}
#path_freq={}
for dir in dir_list:
    path=open("/Users/erkin/Projects/TCGA/Barcodes/"+dir+"/Sample_Kegg_Pathway.txt")
    i=1
    #gene_freq={}
    #project_id=project_dict[dir]
    for line in path:
        if(i==1):
            i+=1
            continue
        else:
            line=line.split('\t')
            path_dir=line[0]
            gene=line[3].split(',')[:-1]
            if(path_all.get(path_dir,'None')=='None'):
                path_all[path_dir]=1
                path_n[path_dir]=gene
            else:
                path_all[path_dir]+=1
                path_n[path_dir]+=gene

out_all=open("All_Samples_Pathways.txt",'w')
out_all.write("Pathway ID"+'\t'+"Pathway Name"+'\t'+"Occurrence"+'\n')

a=0
path_list = sorted(path_all.items(), key=lambda kv: kv[1])
for key in reversed(path_list):
    if(a<10):
        out_2=open("Pathway_Gene_Id_"+key[0][8:]+".txt",'w')
        #list.append(path_all[key])
        d=Gene_Freq(key[0])
        gene_freq_sorted = sorted(d.items(), key=lambda kv: kv[1])
        for i in reversed(gene_freq_sorted):
            for k in gene_map.keys():
                if(i[0]==gene_map[k]):
                    out_2.write(k+'\t'+str(i[1])+'\n')
                    break
        a+=1
        out_2.close()

    out_all.write(key[0]+'\t'+path_name[key[0]]+'\t'+ str(key[1])+'\n')

'''
x=np.arange(len(list))
fig, ax = plt.subplots()
plt.bar(x, list)
plt.xticks(x, path_all.keys(),rotation=90)
plt.show()

out_path=open("Pathway Name.txt",'w')
out_path.write("Pathway ID"+'\t'+"Pathway Name"+'\n')

for k in path_name.keys():
    out_path.write(k+'\t'+path_name[k]+'\n')
'''
