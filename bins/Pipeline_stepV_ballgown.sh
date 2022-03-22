### Do mapping for the sequence reads to the reference genomes

Data=/project/CYY/LQY;

[ ! -d $Data/4_ballgown ] && mkdir $Data/4_ballgown

gtf="/project/Ref_database/reference/rnaseqc/gencode.vM16.annotation.gtf"

for i in $(ls $Data/2_hisat2_align/*_sorted.bam);
do
	ID=$(basename $i)
	stringtie -p 120 -G $Data/3_stringtie/stringtie_merged.gtf -e -B -o $Data/4_ballgown/${ID%_sorted.bam}/${ID%_sorted.bam}.gtf $i
done
