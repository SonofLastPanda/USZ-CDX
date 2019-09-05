#Creates the list of entrez ids of the genes mutated for each patient and creates a list for the 10 most frequent affected pathway.
import os
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
        id=open("/Users/erkin/Projects/PCAWG/Barcodes/"+dir+"/Sample_Kegg_Pathway.txt")
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
dir_list=os.listdir("/Users/erkin/Projects/PCAWG/Barcodes")
dir_list.remove(".DS_Store")
#dir_list.remove("pathway.txt")
url="http://rest.kegg.jp/list/"
path_name={}
gene_map={}

name=open("Pathway Name.txt")
for line in name:
    line=line.strip('\n')
    line=line.split('\t')
    path_name[line[0]]=line[1]


for dir in dir_list:
    id=open("/Users/erkin/Projects/PCAWG/Barcodes/"+dir+"/Kegg.txt")
    id_pathway=open("/Users/erkin/Projects/PCAWG/pathway.txt")
    out=open("/Users/erkin/Projects/Barcodes/"+dir+"/Sample_Kegg_Pathway.txt",'w')
    map=open("/Users/erkin/Projects/PCAWG/Barcodes/"+dir+"/Uniprot.tab")
    for line in map:
        line=line.strip('\r\n')
        line=line.split('\t')
        if(line[8].find(',')==-1):
            gene_map[line[8]]=line[0]
        else:
            line[8]=line[8].strip('"')
            g_l=line[8].split(',')
            for g in g_l:
                gene_map[g]=line[0]
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
    for key in path_dict.keys():
        #out.write(key+'\t'+path_name[key]+'\t'+str(path_dict[key].freq)+'\t')
        p_list=[]
        for i in path_dict[key].name:
            if(gene_map[i] not in p_list):
                #out.write(gene_map[i]+',')
                p_list.append(gene_map[i])
        #out.write('\n')


path_all={}
path_n={}
#path_freq={}
for dir in dir_list:
    path=open("/Users/erkin/Projects/PCAWG/Barcodes"+dir+"/Sample_Kegg_Pathway_id.txt")
    i=1
    #gene_freq={}
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


list=[]
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


x=np.arange(len(list))
fig, ax = plt.subplots()
plt.bar(x, list)
plt.xticks(x, path_all.keys(),rotation=90)
plt.show()


    
