arg <- commandArgs(T)
if(length(arg) != 3){
	cat("Argument: dataDir GroupInfo Out_Dir\n")
	quit('no')
}
GroupInfo <- arg[2]
dataDir <- Sys.glob(arg[1])
Out_Dir <- arg[3]

library(ballgown)
library(genefilter)
a <- read.table(GroupInfo)
bg <- ballgown(samples=dataDir, samplePattern = 'Sample', meas = 'all')
bg_filt <- subset(bg,'rowVars(texpr(bg)) > 0.1', genomesubset=TRUE)
gene_expression <- gexpr(bg_filt) 
write.csv(gene_expression, paste0(Out_Dir,'/gene_expression_fpkm.csv')) 
transcripts_expression <- texpr(bg_filt) 
write.csv(transcripts_expression,paste0(Out_Dir,'/transcripts_expression_fpkm.csv'))
