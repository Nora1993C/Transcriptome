arg <- commandArgs(T)
if(length(arg) != 2){
        cat("Argument: Data_File Out_Dir\n")
        quit('no')
}


dataFile <- arg[1];
outDir <- arg[2];

#dataFile <- "total.Utest.oneSig.data.txt";
#outDir <- "./";

library(ecodist);
FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
SamID <- FileName[1];
Test <- FileName[2];
Level <- FileName[3];

Data <- t(read.table(dataFile,header=T,row.names = 1));

#Methods <- c("difference");
Methods <- c( "euclidean", "bray-curtis", "manhattan","jaccard","difference","sorensen","gower","modgower10", "modgower2");

for (MID in Methods){
Type <- fixdmat(distance(Data, method=MID));
rownames(Type) <- rownames(Data);
colnames(Type) <- rownames(Data);
write.table(Type,file=paste(outDir,"/",SamID,".",Test,".",Level,".",MID,".txt",sep=""),quote = F,sep="\t",col.names=NA);
}

#Corrlations <- c("kendall");
Corrlations <- c("pearson", "kendall", "spearman");
for (CID in Corrlations){
Type <-1- cor(t(Data),method=CID);
write.table(Type,file=paste(outDir,"/",SamID,".",Test,".",Level,".",CID,".txt",sep=""),quote = F,sep="\t",col.names=NA);
}
