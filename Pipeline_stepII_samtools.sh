### Do mapping for the sequence reads to the reference genomes

Data=/project/CYY/LQY;

genome="mm10"

fasta="/project/Ref_database/reference/genome/mm10/mm10.fa"
gtf="/project/Ref_database/reference/rnaseqc/gencode.vM16.annotation.gtf"

#Change to working directory
workDIR=$Data/2_hisat2_align

cd $workDIR

#samtools
for i in $(ls $workDIR/*.sam);
do
	ID=$(basename $i);
	echo samtools view -bS $i \> $workDIR/${ID%.sam}.bam >> ${ID%.sam}_samtools.sh
	echo samtools sort $workDIR/${ID%.sam}.bam \> $workDIR/${ID%.sam}_sorted.bam >> ${ID%.sam}_samtools.sh
	qsub ${ID%.sam}_samtools.sh
done


