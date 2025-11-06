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
    parallel --tmpdir tmp -j 1 \
	'f={}; name=$(basename "$f" .bam | sed "s/_bwa_sortedrmdupregrouped_rm_scmap${map_len}_chimeric_multi//");
	outputsnp='"$outdir"'/$name"_rm_chim_multim_sc${map_len}_varscan_dp5_varfreq0.001p0.5minread1_snp.vcf";
	pileup='"$outpileup"'/$name".pileup";
	samtools index $f;
	samtools mpileup -f '"$REF"' -q 1 -Q 1 $f > $pileup;
	varscan mpileup2snp $pileup --p-value 0.5 --strand-filter 0 --min-var-freq 0.001 --min-coverage 2 --min-reads2 1 --min-avg-qual 1 --output-vcf 1 > $outputsnp;
	' ::: $dir/Ancestral_hybrid_bwa_sortedrmdupregrouped_rm_scmap${map_len}_chimeric_multi.bam
done
wait

#removed --min-reads 1, default is 2 anyway
#samtools mpileup -f '"$REF"' -q 40 -Q 30 $f > $pileup;
