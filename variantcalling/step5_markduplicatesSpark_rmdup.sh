#!/usr/bin/bash
#from here, the data are grouped in samples. So the input is a list of bam files that come from the same sample. The .bam files from the same sample but from different sequence lanes are merged to a single bam.

inputdir=/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/preprocess_files_Leifsonia_split
ls ${inputdir}/*mergedbam_bwaalign.bam > step4_input_filelist
FILES=$(cat step4_input_filelist)
for F in $FILES;
do
    inlist=${F/.bam/.list}
    ls ${F} > $inlist
    outf=${F/.bam/sortedrmdup.bam}
    met=${F/.bam/sortedrmdup.met}
    gatk MarkDuplicatesSpark \
	-I ${inlist} \
	-O $outf\
	-M $met\
	--allow-multiple-sort-orders-in-input \
	--remove-sequencing-duplicates \
	--conf 'spark.executor.cores=32' 

done
