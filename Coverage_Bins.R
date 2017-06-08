# R script to calculate the coverage per bin 
# written by Linda Wampach and Anna Heintz-Buschart, January 2017

# arguments from command
args<-commandArgs(TRUE)
LIB <- args[1] # sample name

load("clusteringWS_NOevil_Pooled_contigs.10.4.Rdata")

PG <- read.delim(paste0(LIB,"/","Binning/clusterFiles/","PG.bins"),header=F,sep=',')
PG.t <- as.vector(t(PG))

Cov.bin <- as.vector(c())
for (i in 1:length(PG.t)){
  Cov.bin[i] <- median(contigInfo$MG_depth[contigInfo$contig %in% clusterRes$contig[clusterRes$cluster==PG.t[i]]])
}

Phylolist <- list()
for (i in 1:nrow(PG)){
setwd(paste0(LIB,"/","Binning/clusterFiles/Phylophlan_output/",LIB,"_cluster.",PG[i,],".faa/"))
pat <- list.files(pattern="conf.txt")
Phylo<- read.delim(pat,header=F,sep='\t')
Phylolist[[i]] <- Phylo # add it to your list
}
Phylo_data <- do.call(rbind, Phylolist)

Phylo_bin <- gsub("cluster.","",Phylo_data$V1)
Phylo <- cbind(Phylo_bin,as.matrix(Phylo_data$V2))

Cov.data <- as.data.frame(cbind(PG,Cov.bin))
Final.data <- merge(Cov.data,Phylo,by.x="V1",by.y="Phylo_bin",sort=F,all=T)
colnames(Final.data) <- c("Bin","Coverage","Phylophlan taxonomy")

setwd(paste0(LIB,"/Analysis/"))
write.table(Final.data,paste0(LIB,"_Coverage.Taxonomy.Bin.PG.txt"),row.names = FALSE,col.names = FALSE)

Cov.20 <- Final.data[as.numeric(Final.data$Coverage)>=20,]
write.table(Cov.20,paste0(LIB,"_Coverage20.Taxonomy.Bin.txt"),row.names = FALSE,col.names = FALSE)

