#!/usr/bin/bash

genome=/biodata/dep_psl/common/culture_collections/fungal_genomes/Plecuc1_AssemblyScaffolds_Repeatmasked.fasta
ancestral=/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/preprocessed_inputfiles/5808_A_sorteddedup.bam

pilon \
    --genome $genome \
    --frags $ancestral \
    --output Plecuc1_hybrid_assembly \
    --vcf \
    --changes \
    --threads 12
