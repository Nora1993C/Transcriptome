### Do mapping for the sequence reads to the reference genomes

Data=/project/CYY/LQY;

genome="mm10"

[ ! -d $Data/3_stringtie ] && mkdir $Data/3_stringtie

gtf="/project/Ref_database/reference/rnaseqc/gencode.vM16.annotation.gtf"

nThreadsSTAR=30

for i in $(ls $Data/2_hisat2_align/*_sorted.bam);
do
	ID=$(basename $i)
	echo stringtie $i \-o $Data/3_stringtie/${ID%_sorted.bam}.gtf \-p $nThreadsSTAR \-G $gtf \-A $Data/3_stringtie/${ID%_sorted.bam}._abund.tab \-B \-e > $Data/3_stringtie/${ID%_sorted.bam}_stringtie.sh
	qsub $Data/3_stringtie/${ID%_sorted.bam}_stringtie.sh
done

