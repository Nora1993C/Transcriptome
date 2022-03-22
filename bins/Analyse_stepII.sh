Data=/project/Transcriptome/LY;

### Calculate Distance for VIP_Sig Data
[ ! -d $Data/8_Sig_Dis ] && mkdir $Data/8_Sig_Dis;
for i in $(ls $Data/6_sig_gene/*.data.txt)
do
       ID=$(basename $i);
       Rscript $Data/bin/A2.1_Calculate_Distance.r $i $Data/8_Sig_Dis;
done

###Get PCA Plot
[ ! -d $Data/9_Sig_PCA ] && mkdir $Data/9_Sig_PCA;
for i in $(ls $Data/8_Sig_Dis/*.txt)
do
        ID=$(basename $i);
        Rscript $Data/bin/A2.2_Plot_pca.r $i $Data/bin/${ID%.*.*.*.*}-grouping.info $Data/9_Sig_PCA 2
done

###Get MDS Plot
[ ! -d $Data/10_Sig_MDS ] && mkdir $Data/10_Sig_MDS;
for i in $(ls $Data/8_Sig_Dis/*.txt)
do
       ID=$(basename $i);
       Rscript $Data/bin/A2.3_Plot_mds.r $i $Data/bin/${ID%.*.*.*.*}-grouping.info $Data/10_Sig_MDS 2
done

###Get Heatmap figure

[ ! -d $Data/11_Sig_Heatmap ] && mkdir $Data/11_Sig_Heatmap;
for i in $(ls $Data/6_sig_gene/*data.txt)
do
        ID=$(basename $i);
        Rscript $Data/bin/A2.4_Plot_heatmap.r $i $Data/bin/${ID%.*.*.*.*}-grouping.info $Data/11_Sig_Heatmap 2
done
