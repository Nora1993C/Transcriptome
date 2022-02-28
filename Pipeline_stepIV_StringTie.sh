### Do mapping for the sequence reads to the reference genomes

Data=/project/CYY/LQY;

genome="mm10"

[ ! -d $Data/3_stringtie ] && mkdir $Data/3_stringtie

gtf="/project/Ref_database/reference/rnaseqc/gencode.vM16.annotation.gtf"

cp $Data/bin/gtflist.txt $Data/3_stringtie
python $Data/bin/prepDE.py -i $Data/3_stringtie/gtflist.txt -g $Data/3_stringtie/gene_count.csv -t $Data/3_stringtie/transcript.csv

ls $Data/3_stringtie/*gtf > mergelist.txt
stringtie --merge -p 100 -o stringtie_merged.gtf mergelist.txt
