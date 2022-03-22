arg <- commandArgs(T)
if(length(arg) != 4){
        cat("Argument: Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}

dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- as.numeric(arg[4]);

#dataFile <- "total.data.txt";
#groupFile <- "total-grouping.info";
#outDir <- "./";
#grp_col <- 2;


SamID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][1];

grpInfo <- read.table(groupFile,header=T);
data <- data.frame(read.table(dataFile,header=T,sep ="\t",row.names=1));
Data <- data[,as.character(grpInfo[,1])];
Data <- Data[rowSums(Data)!=0,]

Pr <- c()
model <- aov(as.numeric(Data[1,]) ~ grpInfo[,2])
adj.p <- matrix(NA,ncol=nrow(TukeyHSD(model)[[1]]),nrow=nrow(Data))
rownames(adj.p) <- rownames(Data)
colnames(adj.p) <- rownames(TukeyHSD(model)[[1]])

for (i in 1:nrow(Data))
{
    model <- aov(as.numeric(Data[i,]) ~ grpInfo[,2])
    Pr[i] <- summary(model)[[1]][1,5]
    adj.p[i,] <- TukeyHSD(model)[[1]][,4]
}
write.csv(adj.p,file=paste(outDir,"/",SamID,".aov.TukeyHSD.csv",sep=""));


Data.sig <- Data[Pr < 0.05,]
Data.sig <- Data.sig[!is.na(Data.sig[,1]),]
write.table(Data.sig,file=paste(outDir,"/",SamID,".aov.Sig.data.txt",sep=""),quote = F,sep="\t",col.names=NA);

Data.sig <- Data[p.adjust(Pr,method='fdr') < 0.05,]
Data.sig <- Data.sig[!is.na(Data.sig[,1]),]
write.table(Data.sig,file=paste(outDir,"/",SamID,".aov.fdr.data.txt",sep=""),quote = F,sep="\t",col.names=NA);

