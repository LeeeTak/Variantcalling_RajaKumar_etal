#!/usr/bin/bash
FILES=$(cat /netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/sampleids.all)
for F in $FILES;
do
    infile=${F}_hybrid_rm_scmap130_chim_multi_haplotypecaller_Ancestralfilt_repeatfilt.vcf
    input=filt2_repeatfiltered/$infile
    outdir=filt3_hardfiltered
    mkdir -p $outdir
    snps=$outdir/${infile/.vcf/.snps.vcf}
    snpsfilt=${snps/snps/snpsfiltered}
    indels=$outdir/${infile/.vcf/.indels.vcf}
    indelsfilt=${indels/indels/indelsfiltered}
    bothfilt=${indelsfilt/indelsfiltered/allfiltered}
    excluded=${bothfilt/allfiltered/excluded}

    gatk SelectVariants \
	-V $input \
	-select-type SNP \
	-O $snps \

    gatk SelectVariants \
	-V $input \
	-select-type INDEL \
	-O $indels \

    gatk VariantFiltration \
	-V $snps \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "SOR > 3.0" --filter-name "SOR3" \
	-filter "FS > 60.0" --filter-name "FS60" \
	-filter "MQ < 40.0" --filter-name "MQ40" \
	-filter "DP < 10.0" --filter-name "DP10" \
	-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
	-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
	-O $snpsfilt \

    gatk VariantFiltration \
	-V $indels \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "FS > 200.0" --filter-name "FS200" \
	-filter "DP < 10.0" --filter-name "DP10" \
	-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
	-O $indelsfilt \

    gatk MergeVcfs \
    	-I $snpsfilt \
	-I $indelsfilt \
	-O $bothfilt \

    gatk SelectVariants \
	-V $bothfilt \
	-O $excluded \
	--exclude-filtered



done
