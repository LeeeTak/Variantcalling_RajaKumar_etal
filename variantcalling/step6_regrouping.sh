#!/usr/bin/bash

dir=preprocess_files_Leifsonia_split
parallel -j 32 \
    'f={}; name=$(basename "$f" .bam | sed "s/_Plecuc1_hybrid_assembly_markilluadpaters_mergedbam_bwaalignsortedrmdup//");
    out='"$dir"'/$name"_hybrid_bwa_sortedrmdupregrouped.bam";
    picard AddOrReplaceReadGroups \
	I=$f\
	O=$out\
	RGID=$name\
	RGLB=$name\
	RGPL=ILLUMINA\
	RGPU=$name\
	RGSM=$name' ::: $dir/*_bwaalignsortedrmdup.bam
