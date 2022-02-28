Data=/project/Transcriptome/LY;
#### 1. Heatmap 
#[ ! -d $Data/12_Heatmap ] && mkdir $Data/12_Heatmap
#cd $Data/12_Heatmap;
#for i in $(ls $Data/6_sig_gene/*data.txt)
#do
#	ID=$(basename $i);
#	Rscript $Data/bin/6_E_hclust_heatmap_01_parameters.R $i $Data/12_Heatmap;
#	Rscript $Data/bin/7_F_hclust_heatmap_02_table.R $i $Data/bin/${ID%.*.*.*.*}-grouping.info $Data/12_Heatmap 2;
#done

#### 2. Enrichment analysis
[ ! -d $Data/13_Enrichment ] && mkdir $Data/13_Enrichment
cd $Data/13_Enrichment;
for i in $(ls $Data/6_sig_gene/*data.txt)
do
	ID=$(basename $i);
	Rscript $Data/bin/A3.1_Enrich_analysis.r $i $Data/13_Enrichment;
done
#for i in $(ls $Data/12_Heatmap/*clustered_sig_types.txt)
#do
#	cat $i | awk '{if($2==1) print $1}' > $Data/13_Enrichment/${ID%clustered_sig_types.txt}clustered_sig_genes_cluster1.txt
#	cat $i | awk '{if($2==2) print $1}' > $Data/13_Enrichment/${ID%clustered_sig_types.txt}clustered_sig_genes_cluster2.txt
#done

#cp $Data/bin/enrich_geneAnswer.yaml $Data/13_Enrichment
#Rscript $Data/bin/8_enrich_geneAnswer.R $Data/13_Enrichment/enrich_geneAnswer.yaml $Data/13_Enrichment/
