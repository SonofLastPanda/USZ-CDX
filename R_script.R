setwd("~/Projects")
library(fields)
library(pathview)
x <-read.table("CNA_cdx_barcode.txt",sep="\t", row.names=1, header=T)
image.plot(as.matrix(x))
colnames(x)[1]


mydata<-read.table("Pathway_Gene_Id_04360.txt",sep='\t',row.names=1,header=T)
l=as.matrix(mydata)
pv.out <- pathview(gene.data = l[, 1]/max(l[,1]), pathway.id = "05200", species = "hsa", gene.idtype ="entrez", out.suffix = "gse16873", kegg.native = T)
