#args<-commandArgs(TRUE)

library(DESeq2)

countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
colData <- read.table("phenodata.tab", sep="\t", row.names=1,header=T)
##check
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))
##make DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ Group)dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj), ]


