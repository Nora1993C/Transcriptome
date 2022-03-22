library(dendsort)
library("gplots")
library(RColorBrewer)

arg <- commandArgs(T)
if(length(arg) != 2){
        cat("Argument: Data_File Out_Dir\n")
        quit('no')
}

dataFile <- arg[1];
outDir <- arg[2];

#dataFile <- "sig_types_genes.txt_exp.txt";
#outDir <- "./";

col.cut.k <- 2
row.cut.k <- 2

RdBu = rev(brewer.pal(11, name="RdBu"))
RdYlBu = rev(brewer.pal(11, name="RdYlBu"))


SamID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][1]
MID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][2]
LID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][3]
SamID <- paste0(SamID,'.',MID,'.',LID)

exp.df <- read.table(dataFile, stringsAsFactors=F, header=T, row.names=1)

my.heatmap <- function(x, exp.df, row.cut.k, col.cut.k){
  print(x)
  cor.method <- x[1]
  hclust.method <- x[2]
  dendsort.method <- x[3]
  hc.r <- hclust(as.dist(1-cor(t(exp.df), method=cor.method)), method=hclust.method)
  gr.r <- cutree(hc.r, k=row.cut.k)
  hc.c <- hclust(as.dist(1-cor(exp.df, method=cor.method)), method=hclust.method)
  gr.c <- cutree(hc.c, k=col.cut.k)
  dend.r <- as.dendrogram(hc.r)
  dend.c <- as.dendrogram(hc.c)
  
#  cluster.exp.df <- cbind(clusterID=gr.r, exp.df)
#  hc.r.order <- order.dendrogram(dendsort(dend.r, type = dendsort.method))
#  cluster.exp.df <- cluster.exp.df[hc.r.order, ]
#  write.table(cluster.exp.df, file="clustered_sig_types.txt", sep="\t",quote=FALSE)

  heatmap.2(as.matrix(exp.df),
           col=RdBu,
           dendrogram ="both",
           Rowv=rev(dendsort(dend.r, type = dendsort.method)), Colv=dendsort(dend.c, type = dendsort.method),
           scale="row", labRow=rownames(exp.df), labCol=colnames(exp.df),
           cexRow=0.9, cexCol=0.9, margins=c(10, 10),
           RowSideColors=RdYlBu[gr.r], ColSideColors=RdYlBu[gr.c],
           symm = F, key = T, keysize =1, trace="none", density.info="none",
           main=paste0("1-cor, ", cor.method, " ", hclust.method, " ", dendsort.method))
}

cor.methods <- c("pearson")
hclust.methods <- c("ward.D2", "complete", "average")
dendsort.methods <- c("min", "average")
comb.vars <- expand.grid(cor.methods, hclust.methods, dendsort.methods)

pdf(file=paste0(SamID,"sig_types_genes_parameters.pdf"), width=6, height=10)
res <- apply(comb.vars, 1, my.heatmap, exp.df=exp.df, row.cut.k=row.cut.k, col.cut.k=col.cut.k)
dev.off()
