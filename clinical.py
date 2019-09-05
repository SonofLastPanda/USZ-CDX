import gzip
import csv
import sys
csv.field_size_limit(sys.maxsize)
'''
barcodes=open("nexistgenes_barcode.txt")
f=gzip.open("TCGA-All.maf.gz")
clinical=open("Barcode_Project.txt",'w')
bar_dict={}

bar_list=[]
for l in csv.reader(barcodes, dialect="excel-tab"):
    barcode=l[0]
    bar_list.append(barcode)
for line in csv.reader(f, dialect="excel-tab"):
    barcode=line[15]
    if(line[-1].find('TCGA')!=-1 and len(line[-1])<10):
        project=line[-1]
    else:
        project="Unknown"
    if(barcode in bar_list):
        if(bar_dict.get(barcode,'None')=='None'):
            bar_dict[barcode]=project
            clinical.write(barcode+'\t'+project+'\n')

'''
clinical=open("Barcode_Project.txt")
project_dict={}
project_all={}
tcga=gzip.open("TCGA-All.maf.gz")

for line in csv.reader(tcga, dialect="excel-tab"):
    pro=line[-1]
    if(pro != "Unknown"):
        if(project_all.get(pro,'None')=='None'):
            project_all[pro]=1
        else:
            project_all[pro]+=1
for line in clinical:
    line=line.strip('\n')
    line=line.split('\t')
    barcode=line[0]
    project=line[1]
    if(project != "Unknown"):
        if(project_dict.get(project,'None')=='None'):
            project_dict[project]=1
        else:
            project_dict[project]+=1

out=open("Project_Frequency.txt",'w')
out.write("Project"+'\t'+"Occurrence"+'\t'+"Occurence in TCGA"+'\t'+"Frequency(%)"+'\n')
list = sorted(project_dict.items(), key=lambda kv: kv[1])

for key in reversed(list):
    out.write(key[0]+'\t'+str(key[1])+'\t'+str(project_all[key[0]])+'\t'+str((float(key[1]/float(project_all[key[0]])))*100)+'\n')
