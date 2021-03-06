#script to remove contigs from the data set which are part of the contaminant contigs
#written by Linda Wampach, November 2016

options("encoding"="utf-8")
library(Biostrings)

args<-commandArgs(TRUE)
LIB <- args[1] #sample name

ann_dir <- paste0(LIB,"/Analysis/annotation/")
gff <- read.delim(paste0(ann_dir,paste0(LIB,".annotation.filt.NOevilPooled.gff")),header=F,stringsAsFactors=F)

gff.short <- gsub(";.*","",gff[,9])
gff.short.2 <- gsub("ID=","",gff.short)
gff.new <- cbind(gff,gff.short.2)

prokka <- readAAStringSet(paste0(ann_dir,"prokka.faa"), format="fasta",use.names=TRUE)
ID <- gsub(" .*","",names(prokka))

prokka.NOevil <- prokka[gsub(" .*","",names(prokka)) %in% gff.new[,10],]

writeXStringSet(prokka.NOevil,paste0(ann_dir,LIB,".prokka.NOevil.faa"),append=FALSE,format="fasta")
