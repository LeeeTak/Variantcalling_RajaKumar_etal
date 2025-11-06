#!/usr/bin/bash

dir=preprocess_files_Leifsonia_split

for len in 10 150; do #150; do
    parallel -j 40 \
	'samtools view -h {} | awk '"'"'function cigMapLen(c,a,i,l,s) {
				      split(c,a,"[MIDNSHP=X]",s);
				      for(i in a) if(s[i] ~ /[MDN=X]/) l+=a[i];
				      return l;
				  }
				  ($1 ~ /^@/) || ($12 !~ /SA:Z:/ && $7 == "=" && cigMapLen($6) > '"$len"')'"'"' | samtools view -b -o {.}_rm_scmap'"${len}"'_chimeric_multi.bam -' ::: ${dir}/*hybrid_bwa_sortedrmdupregrouped.bam &
			      done
			      wait
