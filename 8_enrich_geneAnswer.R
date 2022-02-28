#! /usr/bin/env Rscript

# Usage: Rscript enrich_pathway_geneAnswer.R tag gene_list.yaml
# Example of yaml:
#     down: down_genes.txt
#     up: up_genes.txt
library(KEGG.db)
library(reactome.db)
library(GO.db)
library(yaml)
library(org.Hs.eg.db)
library(GeneAnswers)
library(ggplot2)
library(reshape2)
library(plyr)

arg <- commandArgs(T)
if(length(arg) != 2){
        cat("Argument: Data_File Out_Dir\n")
        quit('no')
}

dataFile <- arg[1];
outDir <- arg[2];

e2s = toTable(org.Hs.egSYMBOL)
tag <- c("Control")

config <- yaml.load_file(dataFile)
species.anno <- "org.Hs.eg.db"
pvalueT <- 0.1

readGeneList <- function(gene.list.file.name, e2s){
  genes <- read.table(gene.list.file.name, header=FALSE, stringsAsFactors=FALSE)[,1]
  genes.entrez <- e2s$gene_id[e2s$symbol %in% genes]
}

output.genes <- function(gAn.ins, tag, categoryType){
  cat.funs <- c(GO="topGOGenes", DOLITE="topDOLITEGenes", KEGG="topPATHGenes",reactome.path="topREACTOME.PATHGenes")
  # cat.funs <- as.character(data.frame(cat=c("GO", "DOLITE", "KEGG", "reactome.path"),
  #   funs=c("topGOGenes", "topDOLITEGenes", "topPATHGenes", "topREACTOME.PATHGenes")))
#  print(cat.funs)
  fun.name <- cat.funs[categoryType]
  FUN <- match.fun(fun.name)
  pathways.num <- length(gAn.ins@genesInCategory)
  genes.num=max(unlist(llply(gAn.ins@genesInCategory,length)))
  # print(genes.num)
  sub.name <- attr(gAn.ins,"name")
  FUN(gAn.ins, orderby='pvalue', top=pathways.num,topGenes=genes.num,
    fileName=paste0(tag,"_", sub.name,"_", categoryType, "_genes.txt"), file=TRUE)
}

output.go.genes <- function(gAn.ins, tag){
  genes.num <- length(gAn.ins@genesInCategory)
  # print(genes.num)
  sub.name <- attr(gAn.ins,"name")
  topGOGenes(gAn.ins, orderby='pvalue', top=genes.num,
    fileName=paste0(tag,"_", sub.name,"_go_genes.txt"), file=TRUE)
}

output.dolite.genes <- function(gAn.ins, tag){
  genes.num <- length(gAn.ins@genesInCategory)
  # print(genes.num)
  sub.name <- attr(gAn.ins,"name")
  topDOLITEGenes(gAn.ins, orderby='pvalue', top=genes.num,
    fileName=paste0(tag,"_", sub.name,"_dolite_genes.txt"), file=TRUE)
}

output.kegg.genes <- function(gAn.ins, tag){
  genes.num <- length(gAn.ins@genesInCategory)
  # print(genes.num)
  sub.name <- attr(gAn.ins,"name")
  topPATHGenes(gAn.ins, orderby='pvalue', top=genes.num,
    fileName=paste0(tag,"_", sub.name,"_kegg_genes.txt"), file=TRUE)
}

output.reactome.genes <- function(gAn.ins, tag){
  genes.num <- length(gAn.ins@genesInCategory)
  # print(genes.num)
  sub.name <- attr(gAn.ins,"name")
  topREACTOME.PATHGenes(gAn.ins, orderby='pvalue', top=genes.num,
    fileName=paste0(tag,"_", sub.name,"_reactome_genes.txt"), file=TRUE)
}

gAnAnalyze <- function(categoryType, tag, genes.list, species.anno, pvalueT){
    gAn <- lapply(genes.list, geneAnswersBuilder,
                       species.anno, categoryType=categoryType, pvalueT=pvalueT)
    gAn.cluster <- getConceptTable(gAn, items='geneNum')
    write.table(-log10(gAn.cluster$IndexTable),
      file=paste0(tag, "_enrich_", categoryType, ".txt"),
      sep="\t", quote=FALSE, row.names=FALSE)

    g.N <- gAn.cluster$CategoriesTable
    g.N$terms <- row.names(g.N)

    ## keep the order of the terms
    g.N$terms <- factor(g.N$terms, levels=rev(g.N$terms))

    g.N <- g.N[!g.N$terms=="Genes / Group", ]
    g.N.long <- melt(g.N)
    terms.p <- as.data.frame(gAn.cluster$IndexTable)
    terms.p$terms <- row.names(terms.p)
    terms.p.long <- melt(terms.p)
    g.N.long$p.val <- -log10(terms.p.long$value)
    pdf(file=paste0(tag, "_enrich_", categoryType, "_heatmap.pdf"), height=9, width=16)

    nonZero.x <- subset(g.N.long, value>0)
    g.out <- ggplot(nonZero.x, aes(variable, terms)) +
      geom_point(aes(colour=value, size=p.val)) +
      scale_color_gradient(low="black", high="steelblue") +
      xlab("Groups") + ylab(categoryType) +
      labs(colour="Gene numbers", size="-Log10(P-val)") +
      theme(panel.background=element_blank()) +
      theme(axis.text.x = element_text(colour = "black", angle=45, hjust=1),
        axis.text.y = element_text(colour = "black"))
#    print(g.out)
    dev.off()

    write.table(g.N.long,
      file=paste0(tag, "_enrich_", categoryType, "_long_form.txt"),
      sep="\t", quote=FALSE, row.names=FALSE)
}

gAnAnalyze.lev3 <- function(categoryType, tag, genes.list, species.anno, pvalueT){
    gAn <- lapply(genes.list, geneAnswersBuilder,
                       species.anno, categoryType=categoryType, pvalueT=pvalueT,
                       level=3)
#    print(categoryType)
    for(i in 1:length(gAn)){
      attr(gAn[[i]],"name") <- names(gAn)[i]
    }
    if(categoryType %in% c("GO", "DOLITE", "KEGG", "reactome.path")){
      try.out.genes <- llply(gAn, output.genes, tag=tag, categoryType=categoryType)
    }
    # if(categoryType=="GO"){
    #   try.out.genes <- llply(gAn, output.go.genes, tag=tag)
    # }
    # if(categoryType=="DOLITE"){
    #   try.out.genes <- llply(gAn, output.dolite.genes, tag=tag)
    # }
    # if(categoryType=="KEGG"){
    #   try.out.genes <- llply(gAn, output.kegg.genes, tag=tag)
    # }
    # if(categoryType=="reactome.path"){
    #   try.out.genes <- llply(gAn, output.reactome.genes, tag=tag)
    # }

    gAn.cluster <- getConceptTable(gAn, items='geneNum')
    write.table(-log10(gAn.cluster$IndexTable),
      file=paste0(tag, "_enrich_", categoryType, "_lev3.txt"),
      sep="\t", quote=FALSE, row.names=TRUE)

    g.N <- gAn.cluster$CategoriesTable
    g.N$terms <- row.names(g.N)

    ## keep the order of the terms
    g.N$terms <- factor(g.N$terms, levels=rev(g.N$terms))

    g.N <- g.N[!g.N$terms=="Genes / Group", ]
    g.N.long <- melt(g.N)
    terms.p <- as.data.frame(gAn.cluster$IndexTable)
    terms.p$terms <- row.names(terms.p)
    terms.p.long <- melt(terms.p)
    g.N.long$p.val <- -log10(terms.p.long$value)
    pdf(file=paste0(tag, "_enrich_", categoryType, "_lev3_heatmap.pdf"), height=9, width=10)

    nonZero.x <- subset(g.N.long, value>0)
    g.out <- ggplot(nonZero.x, aes(variable, terms)) +
      geom_point(aes(colour=value, size=p.val)) +
      scale_color_gradient(low="black", high="steelblue") +
      xlab("Groups") + ylab(categoryType) +
      labs(colour="Gene numbers", size="-Log10(P-val)") +
      theme(panel.background=element_blank()) +
      theme(axis.text.x = element_text(colour = "black", angle=45, hjust=1),
        axis.text.y = element_text(colour = "black"))
   # print(g.out)
    dev.off()

    write.table(g.N.long,
      file=paste0(tag, "_enrich_", categoryType, "_lev3_long_form.txt"),
      sep="\t", quote=FALSE, row.names=FALSE)
}

genes.list <- list()
genes.list <- lapply(config, function(x) {y <- readGeneList(x, e2s); y})

#categoryTypes <- c("GO", "GO.BP", "GO.CC", "GO.MF", "DOLITE", "KEGG", "reactome.path")

categoryTypes <- c("GO", "KEGG")

lapply(categoryTypes, gAnAnalyze.lev3, tag=tag, genes.list=genes.list, species.anno=species.anno, pvalueT=pvalueT)
