#!/usr/bin/bash
#input files are adapter marked bams
#the input genome needs to be indexed with the following command
#bowtie2-build Plecuc1_AssemblyScaffolds_Repeatmasked.fasta Plecuc1_AssemblyScaffolds_Repeatmasked
#and also a dictionary is required
#run creategenomedict.sh

genome=/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/Plecuc1_hybrid_assembly.fasta
tmp=${F/.bam/_temp}
indir=/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/preprocess_files_Leifsonia_split

parallel -j 1 \
    'f={}; name=$(basename "$f" .bam | sed "s/_markilluadapters//"); unmapped='"$indir"'/$name.bam;
    output='"$indir"'/$name"_mergedbam_bwaalign.bam";
    picard SamToFastq \
	I=$f \
	FASTQ=/dev/stdout \
	CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
	TMP_DIR='"$indir"'/temp | \
    bwa mem -M -t 32 -p '"$genome"' /dev/stdin | \
    picard MergeBamAlignment \
	ALIGNED_BAM=/dev/stdin \
	UNMAPPED_BAM=$unmapped \
	OUTPUT=$output \
	R='"$genome"' CREATE_INDEX=true ADD_MATE_CIGAR=true \
	CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
	INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
	PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
	TMP_DIR={.}_temp' ::: $indir/*_markilluadpaters.bam
