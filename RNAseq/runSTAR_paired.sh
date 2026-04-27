#!/usr/bin/bash

N=$(wc -l $1 | cut -d' ' -f1)
FILES=$(cat $1)
INPUTS=(${FILES//' '/ })
#echo $FILES
for ((i=0;i<$N;i++));
do
	ulimit -n 10000
        f1="${INPUTS[i]}_1_trimmed.fq.gz"
	f2="${INPUTS[i]}_2_trimmed.fq.gz"
	#echo $f1 $f2
	STAR \
	--outReadsUnmapped Fastx \
	--genomeDir /netscratch/dep_psl/grp_rgo/taklee/genomes/Hvulgare/GP/concat_Plecuc1pilon_HvGP_BPGv2/HvBPGv2Plecuc1pol_gindex_gtf_149 \
	--readFilesCommand zcat \
	--readFilesIn $f1 $f2 \
	--outSAMtype BAM SortedByCoordinate \
	--limitBAMsortRAM 1200000000 \
	--runThreadN 32 \
	--outFileNamePrefix bam/${INPUTS[i]}_ 

done

