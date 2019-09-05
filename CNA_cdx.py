#Checks whether the number of the copies of the cdx genes in patients are more than 6.
import gzip
import csv


def takeString(word):
    word.strip('\n')
    fil_word=""
    for a in word:
        if(a!=" "):
            fil_word=fil_word+a
        else:
            break
    ind=word.find('#')
    hash_word=word[ind+1:]
    word_list=[fil_word,hash_word]
    return word_list

barcodes=open("nexistgenes_barcode.txt")
cdx=open("cdx_genes.txt")
f=gzip.open("all_samples.consensus_CN.by_gene.170214.txt.gz")
out=open("CNA_cdx_barcode.txt",'w')
cdx_genes=[]
bar_gene=[]
index_list=[]
for line in barcodes:
    line=line.strip('\n')
    line=line.split('\t')
    bar=line[0]
    bar_gene.append(bar)

for line in cdx:
    line=line.strip('\n')
    if( " " in line):
        line=takeString(line)
        cdx_genes.append(line[0])
        cdx_genes.append(line[1])
    else:
        cdx_genes.append(line)

i=1
for line in csv.reader(f, dialect="excel-tab"):
    if(i==1):
        ind=0
        for bar in line:
            if(bar in bar_gene):
                index_list.append(ind)
            ind+=1
        for bar in line[1:]:
            if bar not in bar_gene:
                line.remove(bar)
        out.write('\t'.join(line)+'\n')
        i+=1
    else:
        gene=line[0]
        list=[gene]
        a=False
        if(gene in cdx_genes):
            a=True
        if(a==True):
            for index in index_list :
                list.append(line[index])
            out.write('\t'.join(list)+"\n")
