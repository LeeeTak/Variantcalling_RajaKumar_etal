#!/usr/bin/bash
#the input reference genome fasta file requires an index file
#simply create the index file with: samtools faidx ref.fasta
REF=/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/Plecuc1_hybrid_assembly.fasta
#samtools faidx $REF
#input bam files need to be indexed... again
dir=/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/preprocess_files_Leifsonia_split
outdir=hybridgenome_first_INDELcall_rm_chimeric_multimapped
outpileup=mpileups
mkdir -p $outdir
mkdir -p $outpileup
#map_len=130
#-q 30: minimum mapping quality for alignment
#-Q 20: minimum base quality
parallel --tmpdir /netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/varscan/tmp -j 32 \
    'f={}; name=$(basename "$f" .bam | sed "s/_bwa_sortedrmdupregrouped_rm_scmap10_chimeric_multi//");
    outputindel='"$outdir"'/$name"_rm_chim_multim_varscan_mapq60baseq30dp10varfreq0.1p0.05minread2_indel.vcf";
    pileup='"$outpileup"'/$name".pileup";
    samtools index $f;
    samtools mpileup -f '"$REF"' -q 60 -Q 30 $f > $pileup;
    varscan mpileup2indel $pileup --p-value 0.05 --min-var-freq 0.1 --min-coverage 10 --min-avg-qual 30 --min-reads2 2 --output-vcf 1 > $outputindel' ::: $dir/Ath_G08_T3_hybrid_bwa_sortedrmdupregrouped_rm_scmap10_chimeric_multi.bam


#samtools mpileup -f '"$REF"' -q 30 -Q 20 $f > $pileup; mpileup already generated
#removed --min-reads 1, default is 2 anyway
