#script to remove contigs from the data set which are part of the contaminant contigs
#written by Linda Wampach, November 2016 and adjusted by Anna Heintz-Buschart, June 2017
options("encoding"="utf-8")

#this must be set:
evil_dir <- "./" #/path/to/the/directory/with/the/results/of/contaminant/removing/binning/"

args<-commandArgs(TRUE)
LIB <- args[1] #sample name

evil <- readRDS(paste0(LIB,evil_dir,"evilContigs.RDS"))

ann_dir <- paste0(LIB,"/Analysis/annotation/")
ann <- read.delim("annotation.filt.gff",header=F,stringsAsFactors=F)

evil.1 <- gsub(paste0(LIB,"-"),"",evil)

contigs.to.keep <- setdiff(ann[,1],evil.1)
ann.NOevil <- ann[ann[,1] %in% contigs.to.keep & ann[,1] %in% contig.length$V1[contig.length$V2>=1000],]

write.table(ann.NOevil,paste0(LIB,".annotation.filt.NOevilPooled.contigs1000.gff"),sep="\t",quote=F,col.names=F,row.names=F)
