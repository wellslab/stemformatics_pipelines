#!/bin/bash
RUN_PATH=$1
cd $RUN_PATH
for file in $(ls $RUN_PATH)
do
    SAMPLE=`basename $file`
    fastqc -t 5 ${SAMPLE} -o ./fastqc_out
done


# usage ./fastqc.sh /scratch/yf0/sb0894/7387_test/FASTQ_Files
