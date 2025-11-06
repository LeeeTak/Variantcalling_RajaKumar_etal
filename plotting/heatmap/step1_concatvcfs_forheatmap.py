#!/usr/bin/python
import os,sys
import subprocess
import glob
indir = "filt2_repeatfiltered_sc150"
files = glob.glob(f"{indir}/*repeatfilt.vcf")

for f in files:
    zipcmd = f"bgzip -k {f}"
    subprocess.run(zipcmd,shell=True)
    gzipped = f"{f}.gz"
    indexcmd = f"tabix -p vcf {gzipped}"
    subprocess.run(indexcmd,shell=True)

mergecmd = f"bcftools merge {indir}/*repeatfilt.vcf.gz -o {indir}/repeatfilt_sc150_splitgenome_VERYrelaxedAnc_strictercall.merged.vcf"
subprocess.run(mergecmd,shell=True)

