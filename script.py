import os

non_exist_bar=open("nexistgenes_barcode.txt")

#file=open("All_Genes.txt",'w')
all_gene=[]
i=1
for line in non_exist_bar:
    if(i==1):
        i+=1
    else:
        line=line.strip('\n')
        line=line.split('\t')
        gene_list=line[1].split(',')
        os.mkdir("/Users/erkin/Projects/TCGA/Barcodes/"+line[0])
        os.chdir("/Users/erkin/Projects/TCGA/Barcodes/"+line[0])
        file=open(line[0]+".txt",'w')
        for gene in gene_list[:-1]:
            file.write(gene+'\n')
        file.close()
'''
out=open("All_Genes.txt",'w')
i=1
gene_list=[]
for line in non_exist_bar:
    if(i==1):
        i+=1
    else:
        line=line.strip('\n')
        line=line.split('\t')
        genes=line[1][:-1]
        genes=genes.split(',')
        for gene in genes:
            if(gene not in gene_list):
                out.write(gene+'\n')
                gene_list.append(gene)
'''
