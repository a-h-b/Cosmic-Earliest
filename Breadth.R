# R script to calculate the coverage per bin 
# written by Linda Wampach, February 2017

# arguments from command
args<-commandArgs(TRUE)
readslib <- args[1] # sample name - reads
bin <- args[2] # sample name - bin name
binID <- args[3] # bin name
lib <- args[4] # sample name

G <- read.delim(paste0('genomeSize.cluster.',lib,'.',binID,'.txt'),sep='',header=F)
G <- G[1,1]
Z <- read.delim(paste0(readslib,'.',bin,'.zeros_sample.txt'),sep='',header=F)
Z <- Z[1]
  
B <- (G-Z)/G*100
write.table(B,paste0(readslib,'.',bin,'.breadth.txt'),row.names = FALSE,col.names = FALSE)

