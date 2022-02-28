### Do mapping for the sequence reads to the reference genomes

Data=/project/CYY/LQY;

### 1. build HISAT2 index
### Note: Only need to run 1 time, and the index files located in 
#Requirement:
#1) HISAT2; 
#2) reference genome fasta file, which can be downloaded from UCSC
#3) reference GTF file, I suggest use Gencode annotation file

genome="mm10"
echo $genome

HISAT2genomeDir="/project/Ref_database/reference/hisat2/mm10"
[ ! -d $HISAT2genomeDir ] && mkdir $HISAT2genomeDir

fasta="/project/Ref_database/reference/genome/mm10/mm10.fa"
gtf="/project/Ref_database/reference/rnaseqc/gencode.vM16.annotation.gtf"

##Just for once, have already done
#extract_exons.py $gtf > $HISAT2genomeDir/genome.exon
#extract_splice_sites.py $gtf > $HISAT2genomeDir/genome.ss
#hisat2-build -p 120 --ss $HISAT2genomeDir/genome.ss --exon $HISAT2genomeDir/genome.exon $fasta $HISAT2genomeDir/mm10

## HISAT2 alignment
### number of threads for HISAT
nThreadsSTAR=100

#Change to working directory
[ ! -d $Data/2_hisat2_align ] && mkdir $Data/2_hisat2_align
workDIR=$Data/2_hisat2_align

cd $workDIR
for i in $(ls $Data/1_filted_fastq/*_1_val_1.fq);
do
	read1=$(basename $i)  # fastq file for read1
	read2=${read1%_1_val_1.fq}_2_val_2.fq 
	hisat2 -p $nThreadsSTAR -x $HISAT2genomeDir/mm10 -1 $Data/1_filted_fastq/$read1 -2 $Data/1_filted_fastq/$read2 -S $workDIR/${read1%_1_val_1.fq}.sam
done

