# #! /usr/bin/env Rscript

args <- commandArgs(TRUE)
genes.list <- read.table(args[1], sep="\t", stringsAsFactors=FALSE, quote="")[,1]

file.name <- args[2]
chip.tbl <- read.table(file.name, sep="\t", stringsAsFactors=FALSE, quote="",header=TRUE)


# comb.tbl <- join_all(list(chip.tbl, input.tbl), "id")
# log2.enrich.tbl <- log2((comb.tbl[,2:5] + 1)/ (comb.tbl[,6:9] + 1))
# rownames(log2.enrich.tbl) <- comb.tbl$id

res.tbl <- chip.tbl[match(genes.list, rownames(chip.tbl)), ]
write.table(res.tbl, file=paste0(args[1], "_exp.txt"), quote=FALSE, sep="\t")

