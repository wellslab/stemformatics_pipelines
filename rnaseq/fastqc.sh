#!/bin/bash
RUN_PATH=$1
cd $RUN_PATH
for file in $(ls $RUN_PATH)
do
    SAMPLE=`basename $file`
    fastqc -t 5 ${SAMPLE} -o ../qc/fastqc_out
done

# usage ./fastqc.sh /scratch/yf0/sb0894/rnaseq/datasets/7387_2.0/FASTQ
