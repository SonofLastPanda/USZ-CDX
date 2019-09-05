#Creates a new folder for each patient and opens a txt file with the name of mutated genes, for later use in Kegg pathway mapping.
import os

non_exist_bar=open("nexistgenes_barcode.txt")
#
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
        os.mkdir("/Users/erkin/Projects/PCAWG/Barcodes/"+line[0])
        os.chdir("/Users/erkin/Projects/PCAWG/Barcodes/"+line[0])
        file=open(line[0]+".txt",'w')
        for gene in gene_list[:-1]:
            file.write(gene+'\n')
        file.close()
