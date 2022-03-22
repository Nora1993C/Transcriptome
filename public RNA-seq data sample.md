# RNA-seq

## 数据处理流程：

* 数据资源下载，参考基因组及参考转录组
  * gtf, genome.fa

* 质控， fastqc, multiqc
  * trimmomatic, cutadapt, trim-galore

* 比对
  * STAR, HISAT2, TOPHAT2, BOWTIE2, BWA, SUBREAD
* 计数
  * featureCounts, htseq-counts, bedtools
* normalization 归一化，差异分析等
  * DEseq2, edgeR, limma(voom)

## 软件安装

## conda软件安装

```shell
conda create -y --name rna python=2
source activate rna
conda install -y trimmomatic cutadapt trim_galore
conda install -y bowtie2 hisat2 star
conda install -y subread tophat htseq bedtools deeptools
conda install -y salmon multiqc
#conda deactivate
```

## 数据下载

### 下载SRA的数据

```shell
cat id｜while read id; do (nohup prefetch $id &);done
ps -ef|grep prefetch|awk '{print $2}'|while read id; do kill $id;done

ls *.sra|while read id ; do (nohup fastq-dump $id --gzip --split-3 -O ./ &);done
```

### 质控

```shell
ls *.gz|xargs fastqc -t 10
multiqc ./ #路径仅含有fastqc产生的gz和html
```

### 过滤低质量的reads和去掉接头

```shell
#双端trim_galore -q 25 --phred33 --length 50 -e 0.1 --stringency 3 --paired -o ./clean $fq1 $fq2
#单端trim_galore -q 25 --phred33 --length 50 -e 0.1 --stringency 3 -o ./1_clean/ PRJNA407106/SRR6032289.1.fastq.gz
ls *.gz|while read id ; do (nohup trim_galore $id -q 25 --phred33 --length 36 -e 0.1 --stringency 3 -o ../1_clean &);done
```

### 比对

```shell
### 1. build HISAT2 index
### Note: Only need to run 1 time, and the index files located in 
#Requirement:
#1) HISAT2; 
#2) reference genome fasta file, which can be downloaded from UCSC
#3) reference GTF file, I suggest use Gencode annotation file

genome="mm10"
echo $genome

HISAT2genomeDir="/project/Ref_database/reference/hisat2/mm10"
[ ! -d $STARgenomeDir ] && mkdir $STARgenomeDir

fasta="/project/Ref_database/reference/genome/mm10/mm10.fa"
gtf="/project/Ref_database/reference/rnaseqc/gencode.vM16.annotation.gtf"

##Just for once, have already done
#extract_exons.py $gtf > $HISAT2genomeDir/genome.exon
#extract_splice_sites.py $gtf > $HISAT2genomeDir/genome.ss
hisat2-build -p 120 --ss $HISAT2genomeDir/genome.ss --exon $HISAT2genomeDir/genome.exon $fasta $HISAT2genomeDir/mm10

## HISAT2 alignment
### number of threads for HISAT
nThreadsSTAR=30
#Change to working directory
[ ! -d $Data/2_hisat2_align ] && mkdir $Data/2_hisat2_align
workDIR=$Data/2_hisat2_align

cd $workDIR
for i in $(ls $Data/1_filted_fastq/*_1_val_1.fq);
do
	read1=$(basename $i)  # fastq file for read1
	read2=${read1%_1_val_1.fq}_2_val_2.fq 
	hisat2 -p $nThreadsSTAR -x $HISAT2genomeDir/mm10 -1 $Data/1_filted_fastq/$read1 -2 $Data/1_filted_fastq/$read2 -S $workDIR/${read1%_1_val_1.fq}.sam
	#hisat2 -p $nThreadsSTAR -x $HISAT2genomeDir/mm10 -U $Data/1_filted_fastq/$read1 -S $workDIR/${read1%_1_val_1.fq}.sam  单端
done

###
subjunc -T 5 -i index -r R1.fastq -R R2.fastq -o temp.subjunc.bam
bowtie2 -p 10 -x index -1 R1.fastq -2 R2.fastq -S temp.bowtie.sam
bwa men -t 5 -M index R1.fastq R2.fastq > temp.bwa.sam

###
ls *.sam|while read id; do (samtools sort -O sam -@ 5 -o $(basename $id '.sam').bam $id);done   ##-@ 线程 -O filetype
#ls *.bam|xargs -i samtools index {}
#ls *.bam|xargs -i samtools flagstat -@ 10 {} >
ls *.bam|while read id; do ( samtools flagstat $id > $(basename $id '.bam').flagstat);done 

cat *flagstat|awk '{print $1}'|paste - - - - - - - - - - - - -   #
```

### 计数 count

```shell
featureCounts -T 5 -p -t exon -g gene_id -a /project/Ref_database/reference/rnaseqc/gencode.vM16.annotation.gtf -o all.id.txt *.bam
#解释：-T 1：线程数为1；-p：表示数据为paired-end，双末端测序数据；-t exon表示 feature名称为exon；-g gene_id表示meta-feature名称为gene_id(ensembl名称)；-a $gtf 表示输入的GTF基因组注释文件；-o all.id.txt 设置输出文件类型；*.sort.bam是被分析的bam文件。featureCounts支持通配符*

# 使用multiqc对featureCounts统计结果进行可视化。
multiqc all.id.txt.summary

# 只保留all.id.txt文件的【基因名】和【样本counts】
cat all.id.txt | cut -f1,7- > counts.txt
```

### StringTie

```shell
### Do mapping for the sequence reads to the reference genomes
Data=/project/CYY/RNA-seq;
genome="mm10"

[ ! -d $Data/3_stringtie ] && mkdir $Data/3_stringtie
gtf="/project/Ref_database/reference/rnaseqc/gencode.vM16.annotation.gtf"
nThreads=80
for i in $(ls $Data/2_hisat2_align/*.bam);
do
	ID=$(basename $i)
	stringtie $i -o $Data/3_stringtie/${ID%.bam}.gtf -p $nThreads -G $gtf -A $Data/3_stringtie/${ID%.bam}._abund.tab -B -e

done

Data=/project/CYY/RNA-seq;

genome="mm10"

[ ! -d $Data/3_stringtie ] && mkdir $Data/3_stringtie

gtf="/project/Ref_database/reference/rnaseqc/gencode.vM16.annotation.gtf"
cd $Data/3_stringtie;
cp $Data/bin/gtflist.txt $Data/3_stringtie
python $Data/bin/prepDE.py -i $Data/3_stringtie/gtflist.txt -g $Data/3_stringtie/gene_count.csv -t $Data/3_stringtie/transcript.csv
python $Data/bin/getFPKM.py -i $Data/3_stringtie/gtflist.txt -g $Data/3_stringtie/gene_fpkm.csv -t $Data/3_stringtie/transcript.fpkm.csv


cd $Data/3_stringtie;
ls $Data/3_stringtie/*gtf > mergelist.txt
stringtie --merge -p 100 -G $gtf -o stringtie_merged.gtf mergelist.txt
```



## 差异分析

* DEGs

```R
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(rlist)
library(corrplot)

options(stringsAsFactors = F)
ref <- readxl::read_excel('样本信息 1.xlsx',col_names = F)
rawdata <- read.table('all.id.txt',header = T)

meta <- as_tibble(rawdata[,1:6])
countData <- rawdata[,c(7:ncol(rawdata))]
rownames(countData) <- rawdata$Geneid
names <- data.frame(cbind(id=paste0(ref$...1,'.bam'),name=ref$...2))
rownames(names) <- names$id
colnames(countData) <- names[colnames(countData),2]

pdf('heatmap.pdf')
pheatmap(scale(cor(log2(countData+1))))
dev.off()

nodePar <- list(lab.cex=0.6,pch=c(NA,19),
                cex=0.6,col='#ef7a82')
hc = hclust(dist(t(log2(countData+1))))
pdf('hclust.pdf',w=4,h=5)
plot(as.dendrogram(hc),nodePar=nodePar,horiz=T)
dev.off()

###DEGs
colData <- data.frame(id=colnames(countData),
                      grpname=as.factor(list.rbind(strsplit(colnames(countData),'_'))[,1]))

dds <- DESeqDataSetFromMatrix(countData, colData, design= ~ grpname)
dds <- DESeq(dds)
```

### 选择目的基因的FPKM分析差异

```R
library(tidyverse)
library(pheatmap)
library(rlist)
library(ggplot2)
options(stringsAsFactors = F)
ref <- readxl::read_excel('样本信息 1.xlsx',col_names = F)
rawdata <- read.csv('3_stringtie/gene_fpkm.csv',header = T,row.names = 1)

names <- data.frame(name=ref$...2,
                    grpname=as.factor(list.rbind(strsplit(ref$...2,'_'))[,1]))
rownames(names) <- ref$...1
colnames(countData) <- names[colnames(countData),2]

count <- as_tibble(rawdata) %>% mutate(ID=list.rbind(strsplit(rownames(rawdata),'|',fixed = T))[,1]) %>% mutate(gene=list.rbind(strsplit(rownames(rawdata),'|',fixed = T))[,2])

data <- count %>% gather(sample,FPKM,-ID,-gene) %>% left_join(ref,by=c('sample'='...1')) %>% select(-...3)

Select <- data[data$gene %in% c('Neu1','Neu2','Neu3'),]
#Select <- data[grep('Fgl',data$gene),] %>% add_row(Select)

plot <- Select %>% select(-ID,-sample) %>% spread(...2,FPKM) %>% as.data.frame()
plot.data <- plot[,-1]
rownames(plot.data) <- plot$gene
plot.data <- plot.data[rowSums(plot.data) != 0,]
names <- within(names,{grpname <- factor(grpname,levels=c('Sham','Isc1h','Rep3h','Rep6h','Rep12h','Rep24h'))})
names <- names[order(names$name),]
names <- names[c(26:30,1:5,16:25,6:15),]
plot.data <- plot.data[,names$name]
gg <- as.data.frame(names$grpname)
rownames(gg) <- names$name
colnames(gg) <- 'Groups'
ann_cols <- list(Groups = c('#4980c9','#feefc4','#e57335','#c93835','#7C1435','#06052d'))
names(ann_cols$Groups) <- levels(gg$Groups)
pheatmap(as.matrix(plot.data),scale = 'row',cluster_rows = F,cluster_cols = F,annotation_col = gg,annotation_colors = ann_cols,fontsize_col = 6,gaps_col=cumsum(table(gg$Groups)),filename = 'Neu.heatmap.pdf',width = 6,h=1.8)
```

![Neu.heatmap](/Users/caiyuanyuan/Documents/RNA-seq/Neu.heatmap.png)

```
###################
neu1 <- data %>% filter(gene == 'Neu1') %>% select(-ID,-sample) 
neu1 <- neu1 %>% mutate(grp=as.factor(list.rbind(strsplit(neu1$...2,'_',fixed = T))[,1]))

neu1 <- within(neu1,{grp <- factor(grp,levels=c('Sham','Isc1h','Rep3h','Rep6h','Rep12h','Rep24h'))})

reflect <- tibble(grp=levels(neu1$grp)[2:6],
                     h=c(1,3,6,12,24)) 

stat <- neu1 %>% filter(grp %in% reflect$grp) %>% left_join(reflect)
R = cor(x=stat$FPKM,y=stat$h,method = 'spearman')
P = cor.test(x=stat$FPKM,y=stat$h,method = 'spearman')$p.value

pdf('neu1.pdf',w=4,h=3)
ggplot(data=neu1,aes(grp, FPKM)) + 
  labs(x="", y="Expression (fpkm)",title = 'Neu1') + 
  geom_boxplot(fill=c('#4980c9','#feefc4','#e57335','#c93835','#7C1435','#06052d'),outlier.shape = NA, width=0.5) +
  #scale_fill_manual(values = c('#4980c9','#feefc4','#e57335','#c93835','#7C1435','#06052d')) +
  theme_bw() +
  annotate("text", x = 2.3, y = 23, label=paste0('Spearman correlation with stage:\nR=',signif(R,2),', p=',signif(P,2)), size=3)+
  theme(axis.line = element_line(color="grey20"),
        axis.text.x = element_text(color = "black",size = 12,angle = 45,hjust = 1),
        plot.title = element_text(hjust = .5, face = "bold.italic",size = 15),
        panel.grid = element_blank())
dev.off()
```

![neu1](/Users/caiyuanyuan/Documents/RNA-seq/neu1.png)

