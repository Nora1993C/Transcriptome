#!/bin/bash

Data=/project/CYY/LQY;
TG=/project/bin/trim_galore;

# Trim FASTQs for quality & adaptors. FASTQC output generated also

######################################### Trim low quality bases and adaptors from reads
###Example of trim_galore: cutadapt, fastqc
### ./trim_galore -q 25 --stringency 5 --dont_gzip --fastqc --retain_unpaired -r1 31 -r2 31 --length 30 -o ./ --paired
###  -phred33 -a adaptor1 -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT *_R1.fastq.gz *_R2.fastq.gz 


cd $Data;
[ ! -d 1_filted_fastq ] && mkdir 1_filted_fastq;
echo "STARTING TRIM_GALORE";

for i in $(ls $Data/0_raw_fastq/*_1.fastq);
do
	$TG/trim_galore -q 20 --trim1 --paired --fastqc $i ${i%_1.fastq}_2.fastq -o $Data/1_filted_fastq;
done
