#!/bin/bash

readarray files < $1
cores=$2
for ((i=0;i<${#files[@]};i++));
do
        name=${files[i]//[$'\t\r\n']}
        infile1="${name}_R1.fastq.gz"
        infile2="${name}_R2.fastq.gz"
        #forout=(${infile//'_'/ })
        outfile1="${name}_R1_q20trimmed.fastq.gz"
        outfile2="${name}_R2_q20trimmed.fastq.gz"
        echo $infile1
        cutadapt --cores=$cores -q 20 -m 10 -o $outfile1 -p $outfile2 $infile1 $infile2
done
