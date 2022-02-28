library(DESeq2)
library(genefilter)
library(ggplot2)

source("/project/Transcriptome/bin/3_my_PCA.R")


arg <- commandArgs(T)
if(length(arg) != 4){
        cat("Argument: Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}


dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- as.numeric(arg[4]);


#dataFile <- "Counts_sum.txt";
#groupFile <- "Grouping-info.txt";
#outDir <- "./";
#grp_col <- 2;

exLim <- function(values, ratio) {
  lim.min <- min(values)
  lim.max <- max(values)
  lim.mid <- (lim.max-lim.min)/2 + lim.min
  lim.min.new <- lim.mid - (lim.max-lim.min)*(1+ratio)/2
  lim.max.new <- lim.mid + (lim.max-lim.min)*(1+ratio)/2
  lims <- list(min=lim.min.new, max=lim.max.new)
  lims
}

FC.cutoff <- 2.0

tag <- "IPAH"
res.path <- outDir
setwd(res.path)
reads.cnt.tbl <- read.table(dataFile,
                            stringsAsFactors=FALSE,
                            header=TRUE)#, sep="\t")

grpInfo <- read.table(groupFile,header=T,check.names = FALSE)
rownames(grpInfo) <- grpInfo[,1]

groupname <- c();
for(i in 1:length(colnames(reads.cnt.tbl)))
{
        groupname <- append(groupname,as.character(grpInfo[grep(paste("^",colnames(reads.cnt.tbl)[i],"$",sep=""),grpInfo[,1]),grp_col]));
}


colData <- as.data.frame(cbind(groupname=groupname))

rownames(reads.cnt.tbl) <- reads.cnt.tbl[ ,1]
reads.cnt.tbl <- reads.cnt.tbl[ , -1]

rownames(colData) <- names(reads.cnt.tbl)
#print(colData)

cds <- DESeqDataSetFromMatrix(countData = reads.cnt.tbl,
                              colData = colData,
                              design = ~ groupname)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
vsd = varianceStabilizingTransformation(cds)
vsd.exp <- assay(vsd)
write.table(vsd.exp, file="vsd_exp.txt", sep="\t", quote=FALSE)

cds <- DESeq(cds)
#print(resultsNames(cds))

#library(Rtsne) # Load package
#set.seed(5)
#tsne_out <- Rtsne(unique(t(vsd.exp)), perplexity=2, theta=0.0, dims=3) # Run TSNE
#tsne.out.y <- as.data.frame(tsne_out$Y)
#names(tsne.out.y) <- c("Y1", "Y2", "Y3")

pdf(file=paste("PCA_", tag, "_DESeq2.pdf", sep=""), width=8, height=8)
myplot <- plotPCA.sny(vsd, intgroup = c("groupname"), returnData=TRUE)

names <- rownames(colData)

percentVar <- round(100 * attr(myplot$d12, "percentVar"))
gg.pca <- ggplot(myplot$d12, aes(PC1, PC2, color=groupname))+#, pch=families)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      geom_text(aes(label=rownames(colData), hjust=0.5, vjust=1), size=3)
print(gg.pca)

percentVar <- round(100 * attr(myplot$d13, "percentVar"))
gg.pca <- ggplot(myplot$d13, aes(PC1, PC3, color=groupname))+#, pch=families)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC3: ",percentVar[2],"% variance")) +
      geom_text(aes(label=rownames(colData), hjust=0.5, vjust=1), size=3)
print(gg.pca)

percentVar <- round(100 * attr(myplot$d23, "percentVar"))
gg.pca <- ggplot(myplot$d23, aes(PC2, PC3, color=groupname))+#, pch=families)) +
      geom_point(size=3) +
      xlab(paste0("PC2: ",percentVar[1],"% variance")) +
      ylab(paste0("PC3: ",percentVar[2],"% variance")) +
      geom_text(aes(label=rownames(colData), hjust=0.5, vjust=1), size=3)
print(gg.pca)

#y1.lims <- exLim(tsne.out.y$Y1, 0.2)
#y2.lims <- exLim(tsne.out.y$Y2, 0.2)
#y3.lims <- exLim(tsne.out.y$Y3, 0.2)

#myplot.tsne <- ggplot(tsne.out.y, aes(x=Y1, y=Y2, color=groupname))+#, pch=families)) +
#  geom_point(size=3) +
#  geom_text(aes(label=rownames(colData),hjust=0.5, vjust=1), size=3) +
#  scale_x_continuous(limit=c(y1.lims$min, y1.lims$max)) +
#  scale_y_continuous(limit=c(y2.lims$min, y2.lims$max))
#print(myplot.tsne)

#myplot.tsne <- ggplot(tsne.out.y, aes(x=Y1, y=Y3, color=groupname))+#, pch=families)) +
#  geom_point(size=3) +
#  geom_text(aes(label=rownames(colData),hjust=0.5, vjust=1), size=3) +
#  scale_x_continuous(limit=c(y1.lims$min, y1.lims$max)) +
#  scale_y_continuous(limit=c(y3.lims$min, y3.lims$max))
#print(myplot.tsne)

#myplot.tsne <- ggplot(tsne.out.y, aes(x=Y2, y=Y3, color=groupname))+#, pch=families)) +
#  geom_point(size=3) +
#  geom_text(aes(label=rownames(colData),hjust=0.5, vjust=1), size=3) +
#  scale_x_continuous(limit=c(y2.lims$min, y2.lims$max)) +
#  scale_y_continuous(limit=c(y3.lims$min, y3.lims$max))
#print(myplot.tsne)

#tsne.out.y <- cbind(tsne.out.y, colData)
#myplot.tsne <- ggplot(tsne.out.y, aes(x=Y1, y=Y2, color=groupname))+#, pch=families)) +
#  geom_point(size=3) +
#  geom_text(aes(label=rownames(tsne.out.y),hjust=0.5, vjust=1), size=3) +
#  facet_grid(families~.) +
#  scale_x_continuous(limit=c(y1.lims$min, y1.lims$max)) +
#  scale_y_continuous(limit=c(y2.lims$min, y2.lims$max))
#print(myplot.tsne)

dev.off()

print(colData(cds))
#cds.1 <- DESeq(cds, test="LRT", reduced=~groupname)
#print(resultsNames(cds.1))
#res.LRT.1 <- results(cds.1)
#write.table(res.LRT.1, file="res_LRT_families.txt", sep="\t", quote=FALSE)

cds.2 <- DESeq(cds, test="Wald")#, reduced=~families)
print(resultsNames(cds.2))
res.LRT.2 <- results(cds.2)
write.table(res.LRT.2, file="res_Wald_groupname.txt", sep="\t", quote=FALSE)

