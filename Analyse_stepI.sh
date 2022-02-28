Data=/project/Transcriptome/LY;

### 1. get gene data

[ ! -d $Data/5_gene_data ] && mkdir $Data/5_gene_data
#grep "Geneid" $Data/4_reads_count/Counts_sum.txt > $Data/5_gene_data/total.data.txt
#less $Data/4_reads_count/Counts_sum_rpkm.txt >> $Data/5_gene_data/total.data.txt
#perl -pi -e "s///" $Data/5_gene_data/total.data.txt;

for i in $(ls $Data/bin/*.info)
do
    ID=$(basename $i);
    Rscript $Data/bin/A1.0_Extract_data.r $Data/5_gene_data/total.data.txt $i $Data/5_gene_data/${ID%-*.*}.data.txt
done

for i in $(ls $Data/5_gene_data/*.data.txt);
do
    perl -pi -e "s/^\t/Geneid\t/g" $i;
done

### 1. get sig gene [anova]
[ ! -d $Data/6_sig_gene ] && mkdir $Data/6_sig_gene

for i in $(ls $Data/5_gene_data/*.data.txt);
do
    ID=$(basename $i);
    Rscript $Data/bin/A1.1_Extract_Sig_anova.r $i $Data/bin/${ID%.*.*}-grouping.info $Data/6_sig_gene 2;

done

for i in $(ls $Data/6_sig_gene/*Sig.data.txt);
do
    perl -pi -e "s/^\t/ID\t/g" $i;
done

###Plot Volcano figure
[ ! -d $Data/7_Volcano ] && mkdir $Data/7_Volcano;
Rscript $Data/bin/A1.2_Plot_Volcano.r $Data/5_gene_data/CLO2.data.txt $Data/bin/CLO2-grouping.info $Data/7_Volcano 2 LO2 Control;
Rscript $Data/bin/A1.2_Plot_Volcano.r $Data/5_gene_data/CLPS.data.txt $Data/bin/CLPS-grouping.info $Data/7_Volcano 2 LPS Control;
Rscript $Data/bin/A1.2_Plot_Volcano.r $Data/5_gene_data/LO2M.data.txt $Data/bin/LO2M-grouping.info $Data/7_Volcano 2 LO2DMM LO2;
Rscript $Data/bin/A1.2_Plot_Volcano.r $Data/5_gene_data/LPSM.data.txt $Data/bin/LPSM-grouping.info $Data/7_Volcano 2 LPSDMM LPS;
