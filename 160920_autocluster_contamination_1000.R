#R script to remove contaminating contigs
# criteria in line 57-65 need to be adjusted for other data sets
# written by Anna Heintz-Buschart, September 2016

# arguments from command
args<-commandArgs(TRUE)
LIB <- args[1] #sample name
pk <- as.numeric(args[2])
nn <- as.numeric(args[3])

#libraries
library(caTools)
library(fpc)
library(FNN)
library(RColorBrewer)

#read coordinates
contigInfo <- read.delim(paste0(LIB,"_Contamination.length.tsv"),stringsAsFactors=F) #length and names #change this!

colnames(contigInfo)[1:2] <- c("contig","length")
contigCoord <- read.delim(paste0(LIB,"_concat.contigs.1000.rRNAcut.fa_5mer_clr.coords"),header=F,stringsAsFactors=F,sep=",") #change this!
contigNames <- read.delim(paste0(LIB,"_concat.contigs.1000.rRNAcut.names.txt"),header=F,stringsAsFactors=F) #change this! 
contigNames <- gsub("|","-",contigNames[,1],fixed=T)
contigNames <- gsub(">","",contigNames,fixed=T)
rownames(contigCoord) <- contigNames
colnames(contigCoord) <- c("x","y")
contigInfo <- merge(contigInfo,contigCoord,by.x=1,by.y=0)

contigInfo$conta <- grepl("^Contamination",contigInfo$contig) #check

# auto clustering
coco <- c(which(colnames(contigInfo)=="x"),which(colnames(contigInfo)=="y"))
skn <- sort(knn.dist(contigInfo[,coco],pk)[,pk])
sdkn <- runsd(skn,10)
est <- sort(skn)[min(which(sdkn>quantile(sdkn,0.975)&skn>median(skn)))]
write.table(t(c("scan","reachabilityDist")),paste("reachabilityDistanceEstimates",pk,nn,"tsv",sep="."),row.names=F,col.names=F,quote=F,sep="\t")
write.table(t(c("first",est)),paste("reachabilityDistanceEstimates",pk,nn,"tsv",sep="."),row.names=F,col.names=F,quote=F,sep="\t",append=T)
cdb <- dbscan(contigInfo[,coco],est,pk)

pdf(paste(paste0(LIB,"_scatterPlotConta1"),pk,nn,"pdf",sep="."))
plot(contigInfo[,coco],pch=16,cex=0.25,ann=F,axes=F)
j <- 1
cdbTab <- data.frame("cluster"=names(table(cdb$cluster)),"contigs"=0,"totallength"=0,"contalength"=0,stringsAsFactors=F)
for(i in names(table(cdb$cluster))) {
  points(contigInfo[cdb$cluster==i,coco],pch=16,cex=0.4,col=c("grey",colorRampPalette(brewer.pal(12,"Paired"))(length(unique(cdb$cluster))-1))[j])
  points(contigInfo[cdb$cluster==i&contigInfo$conta,coco],pch=8,cex=0.7,col=c("grey",colorRampPalette(brewer.pal(12,"Paired"))(length(unique(cdb$cluster))-1))[j])
  cdbTab$contigs[cdbTab$cluster==i] <- length(which(cdb$cluster==i))
  cdbTab$totallength[cdbTab$cluster==i] <- sum(contigInfo$length[cdb$cluster==i])
  cdbTab$contalength[cdbTab$cluster==i] <- sum(contigInfo$length[cdb$cluster==i&contigInfo$conta])
  j<-j+1
}
box()
dev.off()
write.table(cdbTab,paste("clusterContaScan",pk,nn,"tsv",sep="."),sep="\t",row.names=F,quote=F)
save.image("clusterTestWS.Rdata")

# save results of evaluation - these values will have to be adjusted for other datasets
totConta <- sum(contigInfo$length[contigInfo$conta])
contigsConta <- contigInfo$contig[!contigInfo$conta&
                                    cdb$cluster %in% cdbTab$cluster[(cdbTab$contalength/totConta)>0.0001&cdbTab$totallength<10000000]&
                                    cdb$cluster != 0]
contigsSafe <- contigInfo$contig[!contigInfo$conta&
                                   cdb$cluster %in% cdbTab$cluster[cdbTab$contalength==0]]
contigsContaBlob <- contigInfo$contig[!contigInfo$conta&
                                    cdb$cluster %in% c(0,cdbTab$cluster[(cdbTab$contalength>0&(cdbTab$contalength/totConta)<=0.0001)|cdbTab$totallength>=10000000])]

# make some plots
pdf(paste(paste0(LIB,"_scatterPlotConta2"),pk,nn,"pdf",sep="."))
#clusters
plot(contigInfo[,coco],pch=16,cex=0.25,ann=F,axes=F)
j <- 1
for(i in cdbTab$cluster) {
  if(cdbTab$totallength[cdbTab$cluster==i]-cdbTab$contalength[cdbTab$cluster==i]>0){
    points(contigInfo[cdb$cluster==i&!contigInfo$conta,coco],pch=16,cex=0.4,col=c("grey",colorRampPalette(brewer.pal(12,"Paired"))(length(which(cdbTab$totallength-cdbTab$contalength>0))-1))[j])
    points(contigInfo[cdb$cluster==i&contigInfo$conta,coco],pch=8,cex=0.7,col=c("grey",colorRampPalette(brewer.pal(12,"Paired"))(length(which(cdbTab$totallength-cdbTab$contalength>0))-1))[j])
    j<-j+1
  }
}
box()
#clusters indicating contaminant contigs in black and "evil contigs" in magenta
plot(contigInfo[,coco],pch=16,cex=0.25,ann=F,axes=F)
j <- 1
for(i in cdbTab$cluster) {
  if(cdbTab$totallength[cdbTab$cluster==i]-cdbTab$contalength[cdbTab$cluster==i]>0){
    points(contigInfo[cdb$cluster==i&!contigInfo$conta,coco],pch=16,cex=0.4,col=c("grey",colorRampPalette(brewer.pal(12,"Paired"))(length(which(cdbTab$totallength-cdbTab$contalength>0))-1))[j])
    points(contigInfo[cdb$cluster==i&contigInfo$conta,coco],pch=8,cex=0.7,col=c("grey",colorRampPalette(brewer.pal(12,"Paired"))(length(which(cdbTab$totallength-cdbTab$contalength>0))-1))[j])
  j<-j+1
  }
}
points(contigInfo[contigInfo$conta,coco],pch=16,cex=0.25)
points(contigInfo[contigInfo$contig %in% contigsConta,coco],pch=4,cex=0.7,col="magenta")
legend("bottomleft",pch=4,cex=0.7,legend="contamination contigs",col="magenta",bty="n")
box()
#clusters indicating contaminant contigs in black and "safe contigs" in green
plot(contigInfo[,coco],pch=16,cex=0.25,ann=F,axes=F)
j <- 1
for(i in cdbTab$cluster) {
  if(cdbTab$totallength[cdbTab$cluster==i]-cdbTab$contalength[cdbTab$cluster==i]>0){
    points(contigInfo[cdb$cluster==i&!contigInfo$conta,coco],pch=16,cex=0.4,col=c("grey",colorRampPalette(brewer.pal(12,"Paired"))(length(which(cdbTab$totallength-cdbTab$contalength>0))-1))[j])
    points(contigInfo[cdb$cluster==i&contigInfo$conta,coco],pch=8,cex=0.7,col=c("grey",colorRampPalette(brewer.pal(12,"Paired"))(length(which(cdbTab$totallength-cdbTab$contalength>0))-1))[j])
    j<-j+1
  }
}
points(contigInfo[contigInfo$conta,coco],pch=16,cex=0.25)
points(contigInfo[contigInfo$contig %in% contigsSafe,coco],pch=4,cex=0.7,col="green")
legend("bottomleft",pch=4,cex=0.7,legend="safe contigs",col="green",bty="n")
box()
#clusters indicating contaminant contigs in black and questionable contigs in cyan
plot(contigInfo[,coco],pch=16,cex=0.25,ann=F,axes=F)
j <- 1
for(i in cdbTab$cluster) {
  if(cdbTab$totallength[cdbTab$cluster==i]-cdbTab$contalength[cdbTab$cluster==i]>0){
    points(contigInfo[cdb$cluster==i&!contigInfo$conta,coco],pch=16,cex=0.4,col=c("grey",colorRampPalette(brewer.pal(12,"Paired"))(length(which(cdbTab$totallength-cdbTab$contalength>0))-1))[j])
    points(contigInfo[cdb$cluster==i&contigInfo$conta,coco],pch=8,cex=0.7,col=c("grey",colorRampPalette(brewer.pal(12,"Paired"))(length(which(cdbTab$totallength-cdbTab$contalength>0))-1))[j])
    j<-j+1
  }
}
points(contigInfo[contigInfo$conta,coco],pch=16,cex=0.25)
points(contigInfo[contigInfo$contig %in% contigsContaBlob,coco],pch=4,cex=0.7,col="cyan")
legend("bottomleft",pch=4,cex=0.7,legend="problematic contigs",col="cyan",bty="n")
box()
dev.off()

saveRDS(contigsConta,"evilContigs.RDS")
saveRDS(contigsContaBlob,"contaBlobContigs.RDS")
saveRDS(contigsSafe,"safeContigs.RDS")


