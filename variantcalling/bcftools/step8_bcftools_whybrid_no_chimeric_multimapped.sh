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
parallel -j 40 \
    'f={}; name=$(basename "$f" .bam | sed "s/_sortedrmdupregrouped_rm_scmap'"$map_len"'_chimeric_multi//"); output='"$outdir"'/$name"_rm_scmap'"$map_len"'_chim_multi_bcftools.vcf.gz";
    bcftools mpileup --max-depth 10000 -Ou -f '"$REF"' $f | bcftools call -mv -Oz -o $output --ploidy 1 --threads 2' ::: $dir/*_hybrid_sortedrmdupregrouped_rm_scmap${map_len}_chimeric_multi.bam
