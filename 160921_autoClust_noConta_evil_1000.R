#R script to bin contigs after removal of contaminant contigs
#call this from the Analysis directory of IMP
# written by Anna Heintz-Buschart and adjusted by Linda Wampach, September 2016

#this need to be set:
contadir <- "./" #/path/to/the/directory/containing/the/results/of/binning/to/remove/contaminants/"
cutdir <- "./" # "/path/to/the/directory/containing/the/results/of/removal/of/rRNAs/and/VizBin/"

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
library(scales)
library(diptest)
library(mixtools)
library(gclus)

#my functions
find.cutoff <- function(data,k=2,proba=0.5) {
  model <- try(normalmixEM(x=data, k=k,maxrestarts=50))
  if(class(model)!="try-error"){
    i <- which.min(model$mu)
    f <- function(x) {
      proba - (model$lambda[i]*dnorm(x, model$mu[i], model$sigma[i]) /(model$lambda[1]*dnorm(x, model$mu[1], model$sigma[1]) + model$lambda[2]*dnorm(x, model$mu[2], model$sigma[2])))
    }
    lower <- min(model$mu)
    upper <- max(model$mu)
    if(f(lower)*f(upper)>0){
      return(median(data))
    }else{
      return(uniroot(f=f, lower=lower, upper=upper)$root)  # Careful with division by zero if changing lower and upper
    }
  }else{
    return(median(data))
  }
}
muClus <- function(clusterName,cRes,resFile){
  covcol <- which(colnames(contigInfo)=="MG_depth")
  dInfo <- contigInfo[cRes$cluster==clusterName,]
  if(dip.test(log10(1+dInfo[,covcol][dInfo$essentialGene!="notEssential"]))$p.value<0.05){
    cuto <- find.cutoff(log10(1+dInfo[,covcol][dInfo$essentialGene!="notEssential"]))
    write.table(t(c(clusterName,cuto)),resFile,append=T,row.names=F,col.names=F,quote=F,sep="\t")
    subs1 <- dInfo[log10(1+dInfo[,covcol])<cuto,]
    num1 <- length(unlist(sapply(subs1$essentialGene[subs1$essentialGene!="notEssential"],function(x) unlist(strsplit(x,split=";")))))
    uni1 <- length(unique(unlist(sapply(subs1$essentialGene[subs1$essentialGene!="notEssential"],function(x) unlist(strsplit(x,split=";"))))))
    if(uni1==0) {
      cRes$cluster[cRes$contig %in% subs1$contig] <- paste(gsub("D","E",cRes$cluster[cRes$contig %in% subs1$contig]),1,sep=".")
    } else if(num1/uni1<=1.2){
      cRes$cluster[cRes$contig %in% subs1$contig] <- paste(gsub("D","C",cRes$cluster[cRes$contig %in% subs1$contig]),1,sep=".")
    } else {
      cRes$cluster[cRes$contig %in% subs1$contig] <- paste(cRes$cluster[cRes$contig %in% subs1$contig],1,sep=".")
      cRes <- muClus(unique(cRes$cluster[cRes$contig %in% subs1$contig]),cRes,resFile)
    }
    subs2 <- dInfo[log10(1+dInfo[,covcol])>=cuto,]
    num2 <- length(unlist(sapply(subs2$essentialGene[subs2$essentialGene!="notEssential"],function(x) unlist(strsplit(x,split=";")))))
    uni2 <- length(unique(unlist(sapply(subs2$essentialGene[subs2$essentialGene!="notEssential"],function(x) unlist(strsplit(x,split=";"))))))
    if(uni2==0) {
      cRes$cluster[cRes$contig %in% subs2$contig] <- paste(gsub("D","E",cRes$cluster[cRes$contig %in% subs2$contig]),2,sep=".")
    } else if(num2/uni2<=1.2){
      cRes$cluster[cRes$contig %in% subs2$contig] <- paste(gsub("D","C",cRes$cluster[cRes$contig %in% subs2$contig]),2,sep=".")
    } else {
      cRes$cluster[cRes$contig %in% subs2$contig] <- paste(cRes$cluster[cRes$contig %in% subs2$contig],2,sep=".")
      cRes <- muClus(unique(cRes$cluster[cRes$contig %in% subs2$contig]),cRes,resFile)
    }
  } else {
    cRes$cluster[cRes$contig %in% dInfo$contig] <- gsub("^.","B",cRes$cluster[cRes$contig %in% dInfo$contig])
  }
  return(cRes)
}

#load WS from IMP
load("results/MG_results.Rdat")

#coordinates
cutrRNA_coords <- read.delim(paste0(cutdir,LIB,"/",LIB,"_concat.contigs.1000.rRNAcut.fa_5mer_clr.coords"),header=F,sep=',')
colnames(cutrRNA_coords) <- c("x","y")
Names <- read.delim(paste0(cutdir,LIB,"/",LIB,"_concat.contigs.1000.rRNAcut.names.txt"),header=F,stringsAsFactors=F) 
colnames(Names) <- "contig"
contigInfo <- data.frame("contig"=Names,cutrRNA_coords,stringsAsFactors=F)

# remove all >
contigInfo$contig <- gsub("[>]","",contigInfo$contig)

# Remove evil contigs from contigInfo
evil_contigs <- readRDS(paste0(contadir,LIB,"/","evilContigs.RDS"))
evil_contigs <- gsub(paste0(LIB,"-"),"",evil_contigs)

contigInfo <- contigInfo[!contigInfo$contig %in% evil_contigs,]

# merge to coverage (metaG) 
contigInfo <- merge(contigInfo,MG.depth,by=1,all.x=T)
contigInfo$MG_depth[is.na(contigInfo$MG_depth)] <- 0

#read essential genes (essential genes were searched using findEssentialGenesInPredictions.sh)
essGenes <- read.table(paste("annotation/prokka.hmm.orfs.essential.hits",sep=""),header=F,stringsAsFactors=F,quote="",strip.white=T,skip=3,colClasses=c("character","NULL","character","NULL","numeric",rep("NULL",14)))
dupliEss <- unique(names(table(essGenes$V1)[table(essGenes$V1)>1]))
dupliID <- dupliEss
for(dupGenes in dupliEss) dupliID[which(dupliID==dupGenes)] <- essGenes$V3[essGenes$V1==dupGenes][which.min(essGenes$V5[essGenes$V1==dupGenes])]
essGenes <- rbind(essGenes[!(essGenes$V1 %in% dupliEss),-3],data.frame("V1"=dupliEss,"V3"=dupliID,stringsAsFactors=F))
colnames(essGenes) <- c("gene","essentialGene")

#read link Genes-Contigs
linkGC <- read.delim(paste("annotation/prokka.bed",sep=""),stringsAsFactors=F,header=F,colClasses=c("character",rep("NULL",2),"character",rep("NULL",2)))
colnames(linkGC) <- c("contig","gene")
linkGC <- linkGC[linkGC$contig %in% contigInfo$contig,]
essGenes <- merge(essGenes,linkGC,by.x=1,by.y=2)
essAnno <- aggregate(essGenes$essentialGene,list(essGenes$contig),function(x) paste(x,sep=";",collapse=";"))
colnames(essAnno) <- c("contig","essentialGene")
contigInfo <- merge(contigInfo,essAnno,by.x=1,by.y=1,all.x=T)
contigInfo$essentialGene[is.na(contigInfo$essentialGene)] <- "notEssential"

#tidy up
rm(list=c("essGenes","dupliEss","dupliID","linkGC"))

# auto clustering
coco <- c(which(colnames(contigInfo)=="x"),which(colnames(contigInfo)=="y"))

#estimate reachability distance (based on number of neighbouting points)
skn <- sort(knn.dist(contigInfo[,coco],pk)[,pk])
sdkn <- runsd(skn,10)
est <- sort(skn)[min(which(sdkn>quantile(sdkn,0.975)&skn>median(skn)))]

#document estimated reachability distances
write.table(t(c("scan","reachabilityDist")),paste("reachabilityDistanceEstimates_NOevil_contigs",pk,nn,"tsv",sep="."),row.names=F,col.names=F,quote=F,sep="\t")
write.table(t(c("first",est)),paste("reachabilityDistanceEstimates_NOevil_contigs",pk,nn,"tsv",sep="."),row.names=F,col.names=F,quote=F,sep="\t",append=T)
#DBscan
cdb <- dbscan(contigInfo[,coco],est,pk)
#count contigs in clusters and document clustering steps
pdf(paste("scatterPlot1_NOevil_contigs",pk,pk,"pdf",sep="."))
plot(contigInfo[,coco],pch=16,cex=0.25,ann=F,axes=F)
j <- 1
cdbTab <- data.frame("cluster"=names(table(cdb$cluster)),"contigs"=0,"numberEss"=0,"uniqueEss"=0,stringsAsFactors=F)
for(i in names(table(cdb$cluster))) {
  points(contigInfo[cdb$cluster==i,coco],pch=16,cex=0.4,col=c("grey",colors(distinct=T))[j])
  cdbTab$contigs[cdbTab$cluster==i] <- table(cdb$cluster)[names(table(cdb$cluster))==i]
  cdbTab$numberEss[cdbTab$cluster==i] <- length(unlist(sapply(contigInfo$essentialGene[contigInfo$essentialGene!="notEssential"&cdb$cluster==i],function(x) unlist(strsplit(x,split=";")))))
  cdbTab$uniqueEss[cdbTab$cluster==i] <- length(unique(unlist(sapply(contigInfo$essentialGene[contigInfo$essentialGene!="notEssential"&cdb$cluster==i],function(x) unlist(strsplit(x,split=";"))))))
  j<-j+1
}
box()
dev.off()
system(paste("mkdir cluster",pk,nn,"Files",sep="."))
write.table(cdbTab,paste("cluster",pk,nn,"Files/clusterFirstScan_NOevil_contigs.tsv",sep="."),sep="\t",row.names=F,quote=F)
write.table(t(c("clusterName","cutoff")),paste("bimodalClusterCutoffs_NOevil_contigs",pk,nn,"tsv",sep="."),sep="\t",row.names=F,col.names=F,quote=F)

#assign clusters to contigs 
### names are "N" for contigs recognized as noise
#### "E" for clusters without essential genes
#### "C" for clusters with less than or equal to 20% duplicated essential genes
#### "D" for clusters with more than 20% duplicated essential genes
clusterRes <- data.frame("contig"=contigInfo$contig,"cluster"="x",stringsAsFactors=F)
clusterRes$cluster[cdb$cluster==0] <- "N"

emptyClus <- cdbTab$cluster[cdbTab$uniqueEss==0&cdbTab$cluster!=0]
clusterRes$cluster[cdb$cluster %in% emptyClus] <- paste("E",cdb$cluster[cdb$cluster %in% emptyClus],sep="")

closeClus <- cdbTab$cluster[cdbTab$numberEss/cdbTab$uniqueEss<=1.2&cdbTab$cluster!=0]
clusterRes$cluster[cdb$cluster %in% closeClus] <- paste("C",cdb$cluster[cdb$cluster %in% closeClus],sep="")

duClus <- cdbTab$cluster[cdbTab$numberEss/cdbTab$uniqueEss>1.2&cdbTab$cluster!=0]
clusterRes$cluster[cdb$cluster %in% duClus] <- paste("D",cdb$cluster[cdb$cluster %in% duClus],sep="")

#cut clusters by metagenomic coverage depth, document cut-offs
write.table(t(c("cluster","cutoff")),paste("bimodalClusterCutoffs_NOevil_contigs",pk,nn,"tsv",sep="."),sep="\t",col.names=F,row.names=F,quote=F)
for(ds in unique(clusterRes$cluster[grep("D",clusterRes$cluster)])){
  clusterRes <- muClus(ds,clusterRes,paste("bimodalClusterCutoffs",pk,nn,"tsv",sep="."))
}

#second iteration of clustering on clusters with more than 20% duplicated essential genes 
### - reachability estimates are based on nth nearest neighbour (independent of number of neighbouring points)
### - neighbouring points are increased by 2
pk2 <- pk + 2
pdf(paste("scatterPlots2_NOevil_contigs",pk,nn,"pdf",sep="."))
for(bb in unique(clusterRes$cluster[grep("B",clusterRes$cluster)])){
  bbInfo <- contigInfo[clusterRes$cluster==bb,]
  if(nrow(bbInfo)>nn){  
    sknBB <- sort(knn.dist(bbInfo[,coco],nn)[,nn])
    if(length(sknBB>10)){
      sdknBB <- runsd(sknBB,10)
      plot(sdknBB,sknBB)
      plot(sdknBB)
      abline(h=quantile(sdknBB,0.975))
      estBB <- sort(sknBB)[min(which(sdknBB>quantile(sdknBB,0.975)&sknBB>median(sknBB)))]
      if(is.na(estBB)){
        estBB <- sort(sknBB)[max(sdknBB)]
        print("estBB could not be found")
      }
      plot(sknBB)
      abline(h=estBB)
      
      write.table(t(c(bb,estBB)),paste("reachabilityDistanceEstimates_NOevil_contigs",pk,nn,"tsv",sep="."),row.names=F,col.names=F,quote=F,sep="\t",append=T)
      BBcdb <- dbscan(bbInfo[,coco],estBB,pk2)
      
      plot(bbInfo[,coco],pch=16,cex=0.25,ann=F)
      j <- 1
      BBcdbTab <- data.frame("cluster"=names(table(BBcdb$cluster)),"contigs"=0,"numberEss"=0,"uniqueEss"=0,stringsAsFactors=F)
      for(i in names(table(BBcdb$cluster))) {
        points(bbInfo[BBcdb$cluster==i,coco],pch=16,cex=0.4,col=c("grey",colors(distinct=T))[j])
        BBcdbTab$contigs[BBcdbTab$cluster==i] <- table(BBcdb$cluster)[names(table(BBcdb$cluster))==i]
        BBcdbTab$numberEss[BBcdbTab$cluster==i] <- length(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";")))))
        BBcdbTab$uniqueEss[BBcdbTab$cluster==i] <- length(unique(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";"))))))
        j<-j+1
      }
      write.table(BBcdbTab,paste("cluster",pk,nn,"Files/blob1Cluster_NOevil_contigs",bb,"tsv",sep="."),sep="\t",row.names=F,quote=F)
      
      BBclusterRes <- data.frame("contig"=bbInfo$contig,"cluster"="x",stringsAsFactors=F)
      BBclusterRes$cluster[BBcdb$cluster==0] <- "N"
      clusterRes$cluster[clusterRes$contig %in% BBclusterRes$contig[BBclusterRes$cluster=="N"]] <- "N"
      
      BBemptyClus <- BBcdbTab$cluster[BBcdbTab$uniqueEss==0&BBcdbTab$cluster!=0]
      BBclusterRes$cluster[BBcdb$cluster %in% BBemptyClus] <- paste("E",BBcdb$cluster[BBcdb$cluster %in% BBemptyClus],sep="")
      BBemptyCont <- BBclusterRes$contig[BBcdb$cluster %in% BBemptyClus]
      clusterRes$cluster[clusterRes$contig %in% BBemptyCont] <- paste(gsub("B","E",bb),
                                                                      gsub("E","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBemptyCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                      sep=".")
      
      BBcloseClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss<=1.2&BBcdbTab$cluster!=0]
      BBclusterRes$cluster[BBcdb$cluster %in% BBcloseClus] <- paste("C",BBcdb$cluster[BBcdb$cluster %in% BBcloseClus],sep="")
      BBcloseCont <- BBclusterRes$contig[BBcdb$cluster %in% BBcloseClus]
      if(length(BBcloseCont>0)) print(paste("new cluster in",bb,"\n"))
      clusterRes$cluster[clusterRes$contig %in% BBcloseCont] <- paste(gsub("B","C",bb),
                                                                      gsub("C","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBcloseCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                      sep=".")
      
      BBduClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss>1.2&BBcdbTab$cluster!=0]
      BBclusterRes$cluster[BBcdb$cluster %in% BBduClus] <- paste("D",BBcdb$cluster[BBcdb$cluster %in% BBduClus],sep="")
      BBduCont <- BBclusterRes$contig[BBcdb$cluster %in% BBduClus]
      uniD <- unique(paste(gsub("B","D",bb),
                           gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                           sep="."))
      clusterRes$cluster[clusterRes$contig %in% BBduCont] <- paste(gsub("B","D",bb),
                                                                   gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                   sep=".")
      for(bds in uniD){
        clusterRes <- muClus(bds,clusterRes,paste("bimodalClusterCutoffs",pk,nn,"tsv",sep="."))
      }
    }else{
      print(paste("Cluster",bb,"too small for knn-approach.\n"))
    }
  }else{
    print(paste("Cluster",bb,"too small for knn-approach.\n"))
  }
}
dev.off()

#third iteration of clustering on clusters with more than 20% duplicated essential genes
### - reachability estimates are based on nth nearest neighbour (independent of number of neighbouring points)
### - neighbouring points are increased by 2 again
pk4 <- pk + 4
pdf(paste("scatterPlots3_NOevil_contigs",pk,nn,"pdf",sep="."))
for(bb in unique(clusterRes$cluster[grep("B",clusterRes$cluster)])){
  bbInfo <- contigInfo[clusterRes$cluster==bb,]
  if(nrow(bbInfo)>nn){  
    sknBB <- sort(knn.dist(bbInfo[,coco],nn)[,nn])
    if(length(sknBB>10)){
      sdknBB <- runsd(sknBB,10)
      plot(sdknBB,sknBB)
      plot(sdknBB)
      abline(h=quantile(sdknBB,0.975))
      estBB <- sort(sknBB)[min(which(sdknBB>quantile(sdknBB,0.975)&sknBB>median(sknBB)))]
      if(is.na(estBB)){
        estBB <- sort(sknBB)[max(sdknBB)]
        print("estBB could not be found")
      }
      plot(sknBB)
      abline(h=estBB)
      
      write.table(t(c(bb,estBB)),paste("reachabilityDistanceEstimates_NOevil_contigs",pk,nn,"tsv",sep="."),row.names=F,col.names=F,quote=F,sep="\t",append=T)
      BBcdb <- dbscan(bbInfo[,coco],estBB,pk4)
      
      plot(bbInfo[,coco],pch=16,cex=0.25,ann=F)
      j <- 1
      BBcdbTab <- data.frame("cluster"=names(table(BBcdb$cluster)),"contigs"=0,"numberEss"=0,"uniqueEss"=0,stringsAsFactors=F)
      for(i in names(table(BBcdb$cluster))) {
        points(bbInfo[BBcdb$cluster==i,coco],pch=16,cex=0.4,col=c("grey",colors(distinct=T))[j])
        BBcdbTab$contigs[BBcdbTab$cluster==i] <- table(BBcdb$cluster)[names(table(BBcdb$cluster))==i]
        BBcdbTab$numberEss[BBcdbTab$cluster==i] <- length(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";")))))
        BBcdbTab$uniqueEss[BBcdbTab$cluster==i] <- length(unique(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";"))))))
        j<-j+1
      }
      write.table(BBcdbTab,paste("cluster",pk,nn,"Files/blob2Cluster_NOevil_contigs",bb,"tsv",sep="."),sep="\t",row.names=F,quote=F)
      
      BBclusterRes <- data.frame("contig"=bbInfo$contig,"cluster"="x",stringsAsFactors=F)
      BBclusterRes$cluster[BBcdb$cluster==0] <- "N"
      clusterRes$cluster[clusterRes$contig %in% BBclusterRes$contig[BBclusterRes$cluster=="N"]] <- "N"
      
      BBemptyClus <- BBcdbTab$cluster[BBcdbTab$uniqueEss==0&BBcdbTab$cluster!=0]
      BBclusterRes$cluster[BBcdb$cluster %in% BBemptyClus] <- paste("E",BBcdb$cluster[BBcdb$cluster %in% BBemptyClus],sep="")
      BBemptyCont <- BBclusterRes$contig[BBcdb$cluster %in% BBemptyClus]
      clusterRes$cluster[clusterRes$contig %in% BBemptyCont] <- paste(gsub("B","E",bb),
                                                                      gsub("E","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBemptyCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                      sep=".")
      
      BBcloseClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss<=1.2&BBcdbTab$cluster!=0]
      BBclusterRes$cluster[BBcdb$cluster %in% BBcloseClus] <- paste("C",BBcdb$cluster[BBcdb$cluster %in% BBcloseClus],sep="")
      BBcloseCont <- BBclusterRes$contig[BBcdb$cluster %in% BBcloseClus]
      if(length(BBcloseCont>0)) print(paste("new cluster in",bb,"\n"))
      clusterRes$cluster[clusterRes$contig %in% BBcloseCont] <- paste(gsub("B","C",bb),
                                                                      gsub("C","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBcloseCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                      sep=".")
      
      BBduClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss>1.2&BBcdbTab$cluster!=0]
      BBclusterRes$cluster[BBcdb$cluster %in% BBduClus] <- paste("D",BBcdb$cluster[BBcdb$cluster %in% BBduClus],sep="")
      BBduCont <- BBclusterRes$contig[BBcdb$cluster %in% BBduClus]
      uniD <- unique(paste(gsub("B","D",bb),
                           gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                           sep="."))
      clusterRes$cluster[clusterRes$contig %in% BBduCont] <- paste(gsub("B","D",bb),
                                                                   gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                   sep=".")
      for(bds in uniD){
        clusterRes <- muClus(bds,clusterRes,paste("bimodalClusterCutoffs",pk,nn,"tsv",sep="."))
      }
    }else{
      print(paste("Cluster",bb,"too small for knn-approach.\n"))
    }
  }else{
    print(paste("Cluster",bb,"too small for knn-approach.\n"))
  }
}
dev.off()

#fourth iteration of clustering on clusters with more than 20% duplicated essential genes
### - reachability estimates are based on nth nearest neighbour (independent of number of neighbouring points)
### - neighbouring points are increased by 2 again
pk6 <- pk + 6
pdf(paste("scatterPlots4_NOevil_contigs",pk,nn,"pdf",sep="."))
for(bb in unique(clusterRes$cluster[grep("B",clusterRes$cluster)])){
  bbInfo <- contigInfo[clusterRes$cluster==bb,]
  if(nrow(bbInfo)>nn){  
    sknBB <- sort(knn.dist(bbInfo[,coco],nn)[,nn])
    if(length(sknBB>10)){
      sdknBB <- runsd(sknBB,10)
      plot(sdknBB,sknBB)
      plot(sdknBB)
      abline(h=quantile(sdknBB,0.975))
      estBB <- sort(sknBB)[min(which(sdknBB>quantile(sdknBB,0.975)&sknBB>median(sknBB)))]
      if(is.na(estBB)){
        estBB <- sort(sknBB)[max(sdknBB)]
        print("estBB could not be found")
      }
      plot(sknBB)
      abline(h=estBB)
      
      write.table(t(c(bb,estBB)),paste("reachabilityDistanceEstimates_NOevil_contigs",pk,nn,"tsv",sep="."),row.names=F,col.names=F,quote=F,sep="\t",append=T)
      BBcdb <- dbscan(bbInfo[,coco],estBB,pk6)
      
      plot(bbInfo[,coco],pch=16,cex=0.25,ann=F)
      j <- 1
      BBcdbTab <- data.frame("cluster"=names(table(BBcdb$cluster)),"contigs"=0,"numberEss"=0,"uniqueEss"=0,stringsAsFactors=F)
      for(i in names(table(BBcdb$cluster))) {
        points(bbInfo[BBcdb$cluster==i,coco],pch=16,cex=0.4,col=c("grey",colors(distinct=T))[j])
        BBcdbTab$contigs[BBcdbTab$cluster==i] <- table(BBcdb$cluster)[names(table(BBcdb$cluster))==i]
        BBcdbTab$numberEss[BBcdbTab$cluster==i] <- length(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";")))))
        BBcdbTab$uniqueEss[BBcdbTab$cluster==i] <- length(unique(unlist(sapply(bbInfo$essentialGene[bbInfo$essentialGene!="notEssential"&BBcdb$cluster==i],function(x) unlist(strsplit(x,split=";"))))))
        j<-j+1
      }
      write.table(BBcdbTab,paste("cluster",pk,nn,"Files/blob3Cluster_NOevil_contigs",bb,"tsv",sep="."),sep="\t",row.names=F,quote=F)
      
      BBclusterRes <- data.frame("contig"=bbInfo$contig,"cluster"="x",stringsAsFactors=F)
      BBclusterRes$cluster[BBcdb$cluster==0] <- "N"
      clusterRes$cluster[clusterRes$contig %in% BBclusterRes$contig[BBclusterRes$cluster=="N"]] <- "N"
      
      BBemptyClus <- BBcdbTab$cluster[BBcdbTab$uniqueEss==0&BBcdbTab$cluster!=0]
      BBclusterRes$cluster[BBcdb$cluster %in% BBemptyClus] <- paste("E",BBcdb$cluster[BBcdb$cluster %in% BBemptyClus],sep="")
      BBemptyCont <- BBclusterRes$contig[BBcdb$cluster %in% BBemptyClus]
      clusterRes$cluster[clusterRes$contig %in% BBemptyCont] <- paste(gsub("B","E",bb),
                                                                      gsub("E","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBemptyCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                      sep=".")
      
      BBcloseClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss<=1.2&BBcdbTab$cluster!=0]
      BBclusterRes$cluster[BBcdb$cluster %in% BBcloseClus] <- paste("C",BBcdb$cluster[BBcdb$cluster %in% BBcloseClus],sep="")
      BBcloseCont <- BBclusterRes$contig[BBcdb$cluster %in% BBcloseClus]
      if(length(BBcloseCont>0)) print(paste("new cluster in",bb,"\n"))
      clusterRes$cluster[clusterRes$contig %in% BBcloseCont] <- paste(gsub("B","C",bb),
                                                                      gsub("C","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBcloseCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                      sep=".")
      
      BBduClus <- BBcdbTab$cluster[BBcdbTab$numberEss/BBcdbTab$uniqueEss>1.2&BBcdbTab$cluster!=0]
      BBclusterRes$cluster[BBcdb$cluster %in% BBduClus] <- paste("D",BBcdb$cluster[BBcdb$cluster %in% BBduClus],sep="")
      BBduCont <- BBclusterRes$contig[BBcdb$cluster %in% BBduClus]
      uniD <- unique(paste(gsub("B","D",bb),
                           gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                           sep="."))
      clusterRes$cluster[clusterRes$contig %in% BBduCont] <- paste(gsub("B","D",bb),
                                                                   gsub("D","",unlist(sapply(clusterRes$contig[clusterRes$contig %in% BBduCont],function(x)BBclusterRes$cluster[BBclusterRes$contig==x]))),
                                                                   sep=".")
      for(bds in uniD){
        clusterRes <- muClus(bds,clusterRes,paste("bimodalClusterCutoffs",pk,nn,"tsv",sep="."))
      }
    }else{
      print(paste("Cluster",bb,"too small for knn-approach.\n"))
    }
  }else{
    print(paste("Cluster",bb,"too small for knn-approach.\n"))
  }
}
dev.off()

#assignment of final names to clusters
### based on completeness
#### "P": more than 100/109 essential genes, less than 115 essential genes in total (<14% duplicated genes, >92% complete)
#### "G": more than 71/109 essential genes (>65% complete), less than 20% in duplicate
#### "O": more than 51/109 essential genes (>47% complete), less than 20% in duplicate
#### "L": more than 31/109 essential genes (>28% complete), less than 20% in duplicate
#### "C": at least 1/109 essential genes (>= 1% complete), less than 20% in duplicate
#### "E": no essential genes
#### "B": at least 1/109 essential genes (>= 1% complete), at least 20% in duplicate
#### "N": noise
for(clus in grep("C",unique(clusterRes$cluster),value=T)){
  uni <- length(unlist(sapply(contigInfo$essentialGene[contigInfo$essentialGene!="notEssential"&clusterRes$cluster==clus],function(x) unlist(strsplit(x,split=";")))))
  ess <- length(unique(unlist(sapply(contigInfo$essentialGene[contigInfo$essentialGene!="notEssential"&clusterRes$cluster==clus],function(x) unlist(strsplit(x,split=";"))))))
  if(uni>100&ess<115){
    clusterRes$cluster[clusterRes$cluster == clus] <- gsub("C","P",clus)
  } else if(uni>71){
    clusterRes$cluster[clusterRes$cluster == clus] <- gsub("C","G",clus)
  } else if(uni>51){
    clusterRes$cluster[clusterRes$cluster == clus] <- gsub("C","O",clus)
  } else if(uni>31){
    clusterRes$cluster[clusterRes$cluster == clus] <- gsub("C","L",clus)
  }
}

#save results
write.table(clusterRes,paste("contigs2clusters_NOevil_contigs",pk,nn,"tsv",sep="."),sep="\t",row.names=F,quote=F)
save(list=c("contigInfo","coco","pk","nn","essAnno","LIB","args","find.cutoff","muClus","clusterRes"),file=paste("clusteringWS_NOevil_contigs",pk,nn,"Rdata",sep="."))
saveRDS(clusterRes,paste("contigs2clusters_NOevil_contigs",pk,nn,"RDS",sep="."))

print("Saved tables. Plotting final map now.")
#plot final map
essPal <- colorRampPalette(brewer.pal(11,"Spectral"))(111)[111:1]
pdf(paste("finalClusterMap_NOevil_contigs",pk,nn,"pdf",sep="."),width=6,height=5,pointsize=8)
layout(matrix(1:2,nrow=1,ncol=2),c(5/6,1/6))
par(mar=rep(1,4),pty="s")
plot(contigInfo[,coco],pch=16,cex=0.25,ann=F,axes=F,bty="o")
cI <- data.frame("essNum"=vector(mode="numeric",length=length(unique(clusterRes$cluster))),
                 "medX"=vector(mode="numeric",length=length(unique(clusterRes$cluster))),
                 "medY"=vector(mode="numeric",length=length(unique(clusterRes$cluster))),
                 "high"=vector(mode="numeric",length=length(unique(clusterRes$cluster))),
                 stringsAsFactors=F)
i <- 1
for(clus in unique(clusterRes$cluster)){
  uni <- length(unique(unlist(sapply(contigInfo$essentialGene[contigInfo$essentialGene!="notEssential"&clusterRes$cluster==clus],function(x) unlist(strsplit(x,split=";"))))))
  cI$essNum[i] <- uni
  rownames(cI)[i] <- clus
  cI$medX[i] <- median(contigInfo[clusterRes$cluster==clus,coco[1]])
  cI$medY[i] <- median(contigInfo[clusterRes$cluster==clus,coco[2]])
  cI$high[i] <- max(contigInfo[clusterRes$cluster==clus,coco[2]])
  rm(uni)
  i <- i +1
}
cI2 <- cI[rownames(cI)!="N",][order(cI$essNum[rownames(cI)!="N"]),]
cI2 <- rbind(cI[rownames(cI)=="N",],cI2)
for(i in 1:nrow(cI2)){
  clus <- rownames(cI2)[i]
  if(clus=="N") pc <- "grey" else if(grepl("E",clus)) pc <- "black" else if(grepl("B",clus)) pc <- brewer.pal(8,"Set3")[8] else{
    pc <- essPal[cI2$essNum[i]]
  }
  points(contigInfo[clusterRes$cluster==clus,coco],cex=0.3,pch=16,col=pc)
  if(cI2$essNum[i]>15&clus!="N") text(cI2$medX[i],cI2$high[i]+0.5,labels=clus,cex=1.0,font=2,col=pc,pos=3)
  rm(pc)
}
mtext(LIB,3,0)
par(mar=c(3,3,3,1),pty="m")
plotcolors(cbind(c(colorRampPalette(brewer.pal(11,"Spectral"))(109),"black"),c(colorRampPalette(brewer.pal(11,"Spectral"))(109),"black")))
hi <- hist(cI$essNum,breaks=0:109,plot=F)
a <- apply(cbind(hi$counts,hi$breaks[-length(hi$breaks)]+1),1,function(x) lines(c(0.5,(x[1]*1.9/max(hi$counts[-1]))+0.5),c(rep(x[2],2))))
abline(h=c(101.5,71.5,51.5,31.5,1.5),col=c(rep("grey20",4),"grey80"),lty=2)
axis(1,seq(from=0.5,to=2.4,length.out=3),labels=c(0,max(hi$counts[-1])/2,max(hi$counts[-1])))
axis(2,seq(from=1,to=11000/100,length.out=6),labels=c(0,20,40,60,80,100),las=1)
mtext("% completeness",2,2)
mtext("clusters",1,2)
dev.off()
