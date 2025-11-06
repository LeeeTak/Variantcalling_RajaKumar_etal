#!/usr/bin/bash
#the input reference genome fasta file requires an index file
#simply create the index file with: samtools faidx ref.fasta
REF=/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/Plecuc1_hybrid_assembly.fasta
#samtools faidx $REF
#input bam files need to be indexed... again
dir=/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/preprocess_files_Leifsonia_split
#-q 40: minimum mapping quality for alignment
#-Q 30: minimum base quality
softclips=(150)
for map_len in ${softclips[@]}; do
    outdir=hybridgenome_first_SNPcall_rm_chimeric_multimapped_sc${map_len}
    outpileup=mpileups_forsnps_sc${map_len}
    mkdir -p $outdir
    mkdir -p $outpileup
    parallel --tmpdir /netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/varscan/tmp -j 32 \
	'f={}; name=$(basename "$f" .bam | sed "s/_bwa_sortedrmdupregrouped_rm_scmap'"${map_len}"'_chimeric_multi//");
	outputsnp='"$outdir"'/$name"_rm_chim_multim_sc'"${map_len}"'_varscan_mapq60baseq30dp10varfreq0.1p0.05minread2_snp.vcf";
	pileup='"$outpileup"'/$name".pileup";
	samtools index $f;
	samtools mpileup -f '"$REF"' -q 60 -Q 30 $f > $pileup;
	varscan mpileup2snp $pileup --p-value 0.05 --min-var-freq 0.1 --min-coverage 10 --min-reads2 2 --min-avg-qual 30 --output-vcf 1 > $outputsnp;
	' ::: $dir/Ath_G08_T3_hybrid_bwa_sortedrmdupregrouped_rm_scmap${map_len}_chimeric_multi.bam
done
wait

#removed --min-reads 1, default is 2 anyway
#samtools mpileup -f '"$REF"' -q 40 -Q 30 $f > $pileup;
