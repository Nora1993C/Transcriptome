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
[ ! -d $STARgenomeDir ] && mkdir $STARgenomeDir

fasta="/project/Ref_database/reference/genome/mm10/mm10.fa"
gtf="/project/Ref_database/reference/rnaseqc/gencode.vM16.annotation.gtf"

##Just for once, have already done
#extract_exons.py $gtf > $HISAT2genomeDir/genome.exon
#extract_splice_sites.py $gtf > $HISAT2genomeDir/genome.ss
hisat2-build -p 120 --ss $HISAT2genomeDir/genome.ss --exon $HISAT2genomeDir/genome.exon $fasta $HISAT2genomeDir/mm10
