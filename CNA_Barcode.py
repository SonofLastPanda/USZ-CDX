#Eliminates the patients that have 6 copies of any of the genes.
import gzip
import csv

f=gzip.open("all_samples.consensus_CN.by_gene.170214.txt.gz")
#f=gzip.open("dummy_1.by_gene.170214.txt.gz")
barcodes=open("nexistgenes_barcode.txt")
out=open("CNA_barcode.txt",'w')
index_list=[]
bar_gene={}
for line in barcodes:
    line=line.strip('\n')
    line=line.split('\t')
    bar=line[0]
    gene_list=line[1].split(',')
    bar_gene[bar]=gene_list[:-1]
i=1
head=[]
for line in csv.reader(f, dialect="excel-tab"):
    if(i==1):
        ind=0
        for bar in line:
            if(bar in bar_gene.keys()):
                index_list.append(ind)
            ind+=1
        for bar in line[1:]:
            if bar not in bar_gene.keys():
                line.remove(bar)

        head=line
        out.write('\t'.join(line)+'\n')
        i+=1
    else:
        gene=line[0]
        list=[gene]
        a=False
        b=False
        for barcode in bar_gene.keys():
            if(gene in bar_gene[barcode]):
                a=True
                break
        if(a==True):
            if(gene == "ERBB2"):
                for ind in range(len(index_list)) :
                    if(line[index_list[ind]][0:1] != 'N'):
                        if(int(line[index_list[ind]][0:1]) >= 5):
                            b=True
                            break
            else:
                for ind in range(len(index_list)) :
                    if(line[index_list[ind]][0:1] != 'N'):
                        #print(head[ind]+ "  "+line[index_list[ind]])
                        if(int(line[index_list[ind]][0:1]) >= 6):
                            b=True
                            break
            if(b):
                for index in range(len(index_list)) :
                    list.append(line[index_list[index]])
                out.write('\t'.join(list)+"\n")

barcodes.seek(0,0)
file=open('CNA_barcode.txt')
i=0
for line in csv.reader(file,dialect="excel-tab"):
    if(i==0):
        i+=1
    else:
        if(line[0]=="ERBB2"):
            for i in range(1,len(line)):
                if(line[i][0:1] != 'N'):
                    if(int(line[i][0:1])>= 5):
                        bar_gene.pop(head[i],'None')
        else:
            for i in range(1,len(line)):
                if(line[i][0:1] != 'N'):
                    if(int(line[i][0:1])>= 6):
                        bar_gene.pop(head[i],'None')
out=open("nexistgenes_barcode_CNA.txt",'w')
out.write('{a}\t{b}\t{c}\t{d}\t\n'.format(a="Tumor_Sample_Barcode",b="Mutated Genes",c="Var_Type",d="Vus_number"))

for l in barcodes:
    l=l.split('\t')
    if(l[0] in bar_gene.keys()):
        out.write('\t'.join(l))
