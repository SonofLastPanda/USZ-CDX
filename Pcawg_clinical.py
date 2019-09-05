#Script to evaluate the clinical data.
import csv
import gzip

class Donor:
    age=0
    sex=""
    cancer_type=""

    def __init__(self, age,sex,type):
        self.age = age
        self.sex =sex
        self.cancer_type=type

def CancerType(line):
    ind=line[0].find('-')
    full=line[0]
    return(full[:ind])

PCAWG=gzip.open("October_2016_whitelist_2583.snv_mnv_indel.maf.gzip")
barcodes=open("nexistgenes_barcode.txt")
PCAWG_donor=open("pcawg_donor_clinical_August2016_v7.tsv")

bar_l=[]
for line in barcodes:
    line=line.strip('\n')
    line=line.split('\t')
    bar=line[0]
    bar_l.append(bar)

bar_dict={}
don_dict={}
for line in csv.reader(PCAWG, dialect="excel-tab"):
    bar_line=line[12]
    donor_id=line[-1]
    if(bar_line in bar_l):
        if(bar_dict.get(donor_id,'None')=='None'):
            bar_dict[donor_id]=bar_line
        else:
            continue
for line in csv.reader(PCAWG_donor, dialect="excel-tab"):
    donor_id=line[2]
    age=line[10]
    sex=line[5]
    type=CancerType(line)
    if(donor_id in bar_dict.keys()):
        if(don_dict.get(donor_id,'None')=='None'):
            d=Donor(age,sex,type)
            don_dict[donor_id]=d

out=open("PCAWG_Clinical.txt",'w')

out.write('{a}\t{b}\t{c}\t{d}\t{e}\n'.format(a="donor_id",b="Tumor_Sample_Barcode",c="Age_at_diagnosis",d="Donor_Sex",e="Cancer Type"))

for key in don_dict.keys():
    out.write(key+'\t'+bar_dict[key]+'\t'+don_dict[key].age+'\t'+don_dict[key].sex+'\t'+don_dict[key].cancer_type+'\n')
