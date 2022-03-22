library("clusterProfiler")
library("DOSE")
library("org.Rn.eg.db")   ### "org.Mm.eg.db" for Mouse, if for Human should be "org.Hs.eg.db"
library("pathview")

arg <- commandArgs(T)
if(length(arg) != 2){
    cat("Argument: Data_File Out_Dir\n")
    quit('no')
}

dataFile <- arg[1];
outDir <- arg[2];

#dataFile <- '6_sig_gene/CLO2.aov.Sig.data.txt'
#outDir <- '13_Enrichment'

FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
SamID <- FileName[1];
MID <- FileName[2];
MID2 <- FileName[3];
SamID <- paste(SamID,MID,MID2,sep='.')

data <- read.table(dataFile,header=F,row.names=1,sep='\t')

gene <- rownames(data);
gene.df <- bitr(gene, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = org.Rn.eg.db)
gene.id <- gene.df[,2];
write.csv(gene.df,paste(outDir,"/",SamID,".gene.list.csv",sep=""),row.names=F);


pdf(paste(outDir,"/",SamID,".KEGGenriched.barplot.pdf",sep=""),width=8,height=5);
ekk <- enrichKEGG(gene=gene.id,"rno",pvalueCutoff=0.1,pAdjustMethod = "fdr");
write.csv(summary(ekk),paste(outDir,"/",SamID,".KEGGenriched.csv",sep=""),row.names=F);
barplot(ekk,drop=TRUE,showCategory = 15);
dev.off()

pdf(paste(outDir,"/",SamID,".GOenriched.barplot.pdf",sep=""),width=8,height=5);
ekk2 <- enrichGO(gene=gene.id,"org.Rn.eg.db",pvalueCutoff=0.05,pAdjustMethod = "fdr");
write.csv(summary(ekk2),paste(outDir,"/",SamID,".GOenriched.csv",sep=""),row.names=F);
barplot(ekk2,drop=TRUE,showCategory = 15);
dev.off()

