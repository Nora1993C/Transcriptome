for i in $(ls *.sra);
do
        fastq-dump --split-3 $i;
done
