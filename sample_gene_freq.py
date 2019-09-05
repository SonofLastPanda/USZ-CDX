
barcodes=open("nexistgenes_barcode.txt")

gene_dict={}
census=open("cancer_gene_census.csv")
all_gene={}
i=1

def Role(word):
    if("oncogene" in word):
        return "oncogene"
    elif("TSG" in word):
        return "TSG"
    else:
        return "fusion"

for line in census:
    if(i==1):
        i+=1
    else:
        line=line.strip('\n')
        line=line.split(',')
        gene=line[0]
        print(gene)
        role=Role(line[14])
        if(role=="oncogene"):
            all_gene[gene]="Oncogene"
        elif(role=="TSG"):
            all_gene[gene]="TSG"
        elif(role=="fusion"):
            all_gene[gene]="Fusion"
        else:
            all_gene[gene]="Null"


for line in barcodes:
    line=line.strip('\n')
    line=line.split('\t')
    genes=line[1][:-1]
    genes=genes.split(',')
    for gene in genes:
        if(gene_dict.get(gene,'None')=='None'):
            gene_dict[gene]=1
        else:
            gene_dict[gene]+=1

out=open("Frequently_Mutated_Genes.txt",'w')
out.write("Gene"+'\t'+"Number of Occurrence of Mutations"+'\n')
list = sorted(gene_dict.items(), key=lambda kv: kv[1])

for key in reversed(list):
    if(all_gene.get(key[0],'None')=='None'):
        out.write(key[0]+'\t'+str(key[1])+'\t'+"Null"+'\n')
    else:
        out.write(key[0]+'\t'+str(key[1])+'\t'+all_gene[key[0]]+'\n')
