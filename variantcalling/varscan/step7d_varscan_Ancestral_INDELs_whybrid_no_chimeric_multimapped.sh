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
parallel --tmpdir tmp -j 1 \
    'f={}; name=$(basename "$f" .bam | sed "s/_bwa_sortedrmdupregrouped_rm_scmap10_chimeric_multi//");
    outputindel='"$outdir"'/$name"_rm_chim_multim_varscan_dp5_varfreq0.001p0.5minread1_indel.vcf";
    samtools index $f;
    pileup='"$outpileup"'/$name".pileup";
    samtools mpileup -f '"$REF"' -q 1 -Q 1 $f > $pileup;
    varscan mpileup2indel $pileup --p-value 0.5 --strand-filter 0 --min-var-freq 0.001 --min-coverage 5 --min-reads2 1 --min-avg-qual 1 --output-vcf 1 > $outputindel' ::: $dir/Ancestral_hybrid_bwa_sortedrmdupregrouped_rm_scmap10_chimeric_multi.bam


#samtools mpileup -f '"$REF"' -q 30 -Q 20 $f > $pileup; mpileup already generated
#removed --min-reads 1, default is 2 anyway
