library(dendsort)
library("gplots")
library(RColorBrewer)
source("/project/Transcriptome/bin/heatmap3.R")

arg <- commandArgs(T)
if(length(arg) != 4){
        cat("Argument:  Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}

dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- as.numeric(arg[4]);

#dataFile <- "sig_types_genes.txt_exp.txt";
#outDir <- "./";
#groupFile <- "Grouping-info.txt";
#grp_col <- 2;

SamID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][1]
MID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][2]
LID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][3]
SamID <- paste0(SamID,'.',MID,'.',LID)

grpInfo <- read.table(groupFile,header=T,check.names = FALSE)
groupname <- grpInfo[,1]


col.cut.k <- 2
row.cut.k <- 2

RdBu = rev(brewer.pal(11, name="RdBu"))
RdYlBu = rev(brewer.pal(11, name="RdYlBu"))
anno1 = brewer.pal(2, name="Dark2")
anno2 = brewer.pal(4, name="Set2")

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

  cluster.exp.df <- cbind(clusterID=gr.r, exp.df)
  hc.r.order <- order.dendrogram(dendsort(dend.r, type = dendsort.method))
  cluster.exp.df <- cluster.exp.df[hc.r.order, ]
  write.table(cluster.exp.df, file=paste0(SamID,"clustered_sig_types.txt"), sep="\t",quote=FALSE)
#  saveRDS(hc.r, file="hc.r.RDS")
#  saveRDS(hc.c, file="hc.c.RDS")

  ## try multiple column annotations
  col.anno <- cbind(anno2[as.factor(groupname)], RdYlBu[gr.c])
  row.anno <- as.matrix(t(RdYlBu[gr.r]))
#  print(col.anno)
#  print(row.anno)

  heatmap.2(as.matrix(exp.df),
           col=RdBu,
           dendrogram ="both",
           Rowv=rev(dendsort(dend.r, type = dendsort.method)), Colv=dendsort(dend.c, type = dendsort.method),
           scale="row", labRow="", labCol=colnames(exp.df),
           cexCol=0.9, margins=c(10, 5),
           RowSideColors=RdYlBu[gr.r], ColSideColors=RdYlBu[gr.c],
           symm = F, key = T, keysize =1, trace="none", density.info="none",
           main=paste0("1-cor, ", cor.method, " ", hclust.method, " ", dendsort.method))

   heatmap.3(as.matrix(exp.df),
             col=RdBu,
             dendrogram ="both",
            Rowv=rev(dendsort(dend.r, type = dendsort.method)), Colv=dendsort(dend.c, type = dendsort.method),
            scale="row", labRow="", labCol=colnames(exp.df),
            cexCol=0.9, margins=c(10, 5),
            RowSideColors=row.anno, ColSideColors=col.anno,
            # RowSideColors=col.anno, ColSideColors=RdYlBu[gr.c],
            symm = F, key = T, keysize =1, trace="none", density.info="none",
            main=paste0("1-cor, ", cor.method, " ", hclust.method, " ", dendsort.method))

    heatmap.3(as.matrix(exp.df),
              col=RdBu,
              dendrogram ="both",
             Rowv=rev(dendsort(dend.r, type = dendsort.method)), Colv=dendsort(dend.c, type = dendsort.method),
             scale="row", labRow=rownames(exp.df), labCol=colnames(exp.df),
             cexCol=0.9, margins=c(10, 5),
             RowSideColors=row.anno, ColSideColors=col.anno,
             # RowSideColors=col.anno, ColSideColors=RdYlBu[gr.c],
             symm = F, key = T, keysize =1, trace="none", density.info="none",
             main=paste0("1-cor, ", cor.method, " ", hclust.method, " ", dendsort.method))

}


 cor.methods <- c("pearson")
 hclust.methods <- c("average")
 dendsort.methods <- c("average")
comb.vars <- expand.grid(cor.methods, hclust.methods, dendsort.methods)

pdf(file=paste0(SamID,"sig_types_genes_nteraction.pdf"), width=6, height=10)
res <- apply(comb.vars, 1, my.heatmap, exp.df=exp.df, row.cut.k=row.cut.k, col.cut.k=col.cut.k)
dev.off()
#saveRDS(res[[1]], file="hc_heatmap2.RDS")
