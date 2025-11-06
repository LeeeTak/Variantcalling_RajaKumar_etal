#!/usr/bin/bash
#the input reference genome fasta file requires an index file
#simply create the index file with: samtools faidx ref.fasta
REF=/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/Plecuc1_hybrid_assembly.fasta
#samtools faidx $REF
#input bam files need to be indexed... again
dir=/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/from_hybrid_preprocessed_inputfiles
outdir=hybridgenome_first_variantcall_rm_softclipped_chimeric_multimapped
mkdir -p $outdir

map_len=130
parallel -j 32 \
    'f={}; name=$(basename "$f" .bam | sed "s/_sortedrmdupregrouped_rm_scmap'"$map_len"'_chimeric_multi//");
    output='"$outdir"'/$name"_rm_scmap'"$map_len"'_chim_multi_freebayes.vcf";
    samtools index -M $f;
    freebayes -f '"$REF"' -C 10 -p 1 $f > $output' ::: $dir/*_hybrid_sortedrmdupregrouped_rm_scmap${map_len}_chimeric_multi.bam

