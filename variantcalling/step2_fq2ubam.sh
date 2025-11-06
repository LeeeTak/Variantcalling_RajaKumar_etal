#!/usr/bin/bash

indir=/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/bbsplit
outdir=/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/preprocess_files_Leifsonia_split

parallel -j 32 \
    'in1={};in2=${in1/R1/R2};
     out='"$outdir"'/$(basename "$in1" .fastq.gz | sed "s/_R1//").bam;
     rg=$(zcat $in1 | head -1 | cut -f 2,4 -d":");
     IFS="_" read -a sn <<< $(basename "$in1");
     sampname=${sn[0]}_${sn[1]}_${sn[2]};
     picard FastqToSam \
	 FASTQ=$in1\
	 FASTQ2=$in2\
	 OUTPUT=$out\
	 READ_GROUP_NAME=$rg\
	 SAMPLE_NAME=$sampname\
	 LIBRARY_NAME=pcucumerina\
	 PLATFORM=illumina' ::: $indir/*_hybrid_assembly_R1.fastq.gz
