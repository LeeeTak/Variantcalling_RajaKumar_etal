#!/usr/bin/bash

indir=/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/preprocess_files_Leifsonia_split

parallel -j 32 \
   picard MarkIlluminaAdapters \
   I={}\
   O={.}_markilluadpaters.bam\
   M={.}_mark_metrics.txt\
   TMP_DIR={.}_temp/ ::: $indir/*.bam

