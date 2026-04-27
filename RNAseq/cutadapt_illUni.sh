#!/bin/bash

readarray files < $1
cores=$2
for ((i=0;i<${#files[@]};i++));
do
        name=${files[i]//[$'\t\r\n']}
        infile1="${name}_1.fq.gz"
        infile2="${name}_2.fq.gz"
        #forout=(${infile//'_'/ })
        outfile1="${name}_1_trimmed.fq.gz"
        outfile2="${name}_2_trimmed.fq.gz"
        echo $infile1
        cutadapt --cores=$cores -a "AGATCGGAAGAG" -A "AGATCGGAAGAG" -m 130 -o $outfile1 -p $outfile2 $infile1 $infile2
done
