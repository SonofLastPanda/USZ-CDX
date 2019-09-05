#eliminates the patients with gene fusions
import csv
barcodes=open("nexistgenes_barcode.txt")
pcawg=open("pcawg-3.gene.fusions.V1.tsv")
out=open("fusion_barcode.txt",'w')
bar_gene=[]
barcode_dict={}
for line in barcodes:
    line=line.strip('\n')
    line=line.split('\t')
    bar=line[0]
    bar_gene.append(bar)
    barcode_dict[bar]=[]

i=1
for line in csv.reader(pcawg, dialect="excel-tab"):
    if(i==1):
        i+=1
    else:
        bar_line=line[24]
        if(bar_line in bar_gene):
            fus_gene=line[2]
            barcode_dict[bar_line].append(fus_gene)

out.write('{a}\t{b}\t{c}\n'.format(a="Tumor_Sample_Barcode",b="Fusion Genes",c="# of Fusion Genes"))
w=""

for key in barcode_dict.keys():
    out.write(key)
    out.write('\t')
    for gene in barcode_dict[key]:
        out.write(gene+',')
    out.write('\t')
    out.write(str(len(barcode_dict[key])))
    out.write('\n')
out.close()

barcodes.seek(0,0)
out_2=open("nexistgenes_barcode_fusion.txt",'w')
#out_2.write('{a}\t{b}\t{c}\t{d}\t\n'.format(a="Tumor_Sample_Barcode",b="Mutated Genes",c="Var_Type",d="Vus_number"))

for line in barcodes:
    line=line.split('\t')
    if(len(barcode_dict[line[0]])==0):
        out_2.write('\t'.join(line))
