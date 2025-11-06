#!/usr/bin/bash
#the input reference genome fasta file requires an index file
#simply create the index file with: samtools faidx ref.fasta
REF=/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/Plecuc1_hybrid_assembly.fasta
#samtools faidx $REF
#input bam files need to be indexed... again
dir=/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/contaminant_removed_preprocess_files
outdir=hybridgenome_rm_contaminants_softclipped_chimeric_multimapped
mkdir -p $outdir
#5809_Z_hybrid_bbduk_bwa_sortedrmdupregrouped_rm_scmap130_chimeric_multi.bam
for map_len in 130; do
    parallel -j 32 \
	'f={}; name=$(basename "$f" .bam | sed "s/_sortedrmdupregrouped_rm_scmap'"$map_len"'_chimeric_multi//"); output='"$outdir"'/$name"_rm_contaminant_scmap'"$map_len"'_chim_multi_haplotypecaller.vcf.gz";
	samtools index $f;
	gatk HaplotypeCaller \
	    -R '"$REF"'\
	    -I $f\
	    -O $output\
	    -ploidy 1' ::: $dir/*_hybrid_bbduk_bwa_sortedrmdupregrouped_rm_scmap${map_len}_chimeric_multi.bam &
	done
	wait
