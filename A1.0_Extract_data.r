arg <- commandArgs(T)
if(length(arg) != 3){
        cat("Argument: Data_File Group_File Out_File\n")
        quit('no')
}

dataFile <- arg[1];
groupFile <- arg[2];
Out_File <- arg[3];

#dataFile <- "5_gene_data/total.data.txt";
#groupFile <- "bin/LPS-grouping.info";
#Out_File <- "5_gene_data/LPS.data.txt";

grpInfo <- read.table(groupFile,header=T);
data <- data.frame(read.table(dataFile,header=T,sep ="\t",row.names=1));
Data <- data[,as.character(grpInfo[,1])];

write.table(Data,file=Out_File,quote = F,sep="\t",col.names=NA)


