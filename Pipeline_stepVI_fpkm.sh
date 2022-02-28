### Do mapping for the sequence reads to the reference genomes

Data=/project/Transcriptome/ZYP;

[ ! -d $Data/5_gene_count ] && mkdir $Data/5_gene_count

gtf="/project/Ref_database/reference/rnaseqc/gencode.vM16.annotation.gtf"

for i in $(ls $Data/2_hisat2_align/*_sorted.bam);
do
	ID=$(basename $i)
#	mv $Data/4_ballgown/${ID%_sorted.bam} $Data/4_ballgown/Sample_${ID%_sorted.bam}
done

ls $Data/4_ballgown/ > $Data/bin/Group.info
Rscript ballgown.r $Data/4_ballgown/ $Data/bin/Group.info $Data/5_gene_count

featureCounts -T 64 -t exon -g gene_name -s 1 -a $gtf -o $Data/5_gene_count/Join_mRNA_counts.txt $Data/2_hisat2_align/*_sorted.bam
Rscript $Data/bin/0_calFPKM.R $Data/5_gene_count/Join_mRNA_counts.txt

