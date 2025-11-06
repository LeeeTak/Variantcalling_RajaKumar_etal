#!/usr/bin/bash

FILES=$(cat samples.list)
mkdir -p filt3_hardfiltered
for F in $FILES;
do
    infile=${F}_hybrid_rm_scmap130_chim_multi_freebayes_Ancestralfilt_repeatfilt.vcf
    input=filt2_repeatfiltered/$infile
    filtered=${infile/.vcf/.filtered.vcf}
    excluded=${filtered/filtered/excluded}

    gatk VariantFiltration \
	-V $input \
	-filter "QUAL < 1.0" --filter-name "QUAL1" \
	-filter "DP < 10.0" --filter-name "DP10" \
	-O filt3_hardfiltered/$filtered \

    gatk SelectVariants \
	-V filt3_hardfiltered/$filtered \
	-O filt3_hardfiltered/$excluded \
	--exclude-filtered



done
#-filter "QUAL < 1.0" --filter-name "QUAL1" \
#-filter "SRP < 1.0" --filter-name "SRP10" \
#-filter "SAP < 1.0" --filter-name "SAP10" \
#-filter "EPP < 1.0" --filter-name "EPP10" \
#-filter "DP < 10.0" --filter-name "DP10" \
#-filter "SAF < 1" --filter-name "SAF1" \
#-filter "SAR < 1" --filter-name "SAR1" \
#-filter "RPR < 1" --filter-name "RPR1" \
#-filter "RPL < 1" --filter-name "RPL1" \
