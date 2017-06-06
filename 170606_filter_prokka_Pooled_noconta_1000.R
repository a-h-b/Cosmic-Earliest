#script to remove contigs from the data set which are either part of the contaminant contigs
# or below 1000 in length (we don't know where they are from)
#written by Linda Wampach, November 2016 and adjusted by Anna Heintz-Buschart June 2017

options("encoding"="utf-8")
library(Biostrings)

args<-commandArgs(TRUE)
LIB <- args[1] #sample name

ann_dir <- paste0(LIB,"/Analysis/annotation/")
gff <- read.delim(paste0(ann_dir,paste0(LIB,".annotation.filt.NOevilPooled.gff")),header=F,stringsAsFactors=F)

gff.short <- gsub(";.*","",gff[,9])
gff.short.2 <- gsub("ID=","",gff.short)
gff.new <- cbind(gff,gff.short.2)

contig.length <- read.delim(paste0(LIB,"/Analysis/mg.assembly.length.txt"),header=F,stringsAsFactors=F)
gff.new.long <- gff.new[gff.new[,1] %in% contig.length$V1[contig.length$V2>=1000],]

prokka <- readAAStringSet(paste0(ann_dir,"prokka.faa"), format="fasta",use.names=TRUE)
ID <- gsub(" .*","",names(prokka))

prokka.NOevil <- prokka[gsub(" .*","",names(prokka)) %in% gff.new.long[,10],]

writeXStringSet(prokka.NOevil,paste0(ann_dir,LIB,".prokka.NOevil.contigs1000.faa"),append=FALSE,format="fasta")
