#!/usr/bin/python
import os,sys
import glob
import subprocess
#Left-align and normalize indels, check if REF alleles match the reference, split multiallelic sites into multiple rows; recover multiallelics from multiple rows.

ref = "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/Plecuc1_hybrid_assembly.fasta"
indir = "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/varscan/filt2_repeatfiltered_sc150"
files = glob.glob(f"{indir}/*_hybrid_rm_chim_multim_sc150_varscan_mapq60baseq30dp10varfreq0.1p0.05minread2_Ancestralfilt_repeatfilt.vcf.gz")
outdir = "normalised"
subprocess.run(f"mkdir -p {outdir}",shell=True)
subprocess.run(f"mkdir -p {outdir}/mSly",shell=True)
subprocess.run(f"mkdir -p {outdir}/Ath",shell=True)
subprocess.run(f"mkdir -p {outdir}/Hvu",shell=True)
subprocess.run(f"mkdir -p {outdir}/system",shell=True)
for f in files:
    fa = f.split("/")[-1]
    outfile = f"{outdir}/{fa.split('_')[0]}/{fa.replace('.vcf.gz','.norm.vcf.gz')}"
    cmd= f"bcftools norm -f {ref} -m -both -Oz -o {outfile} {f}"
    subprocess.run(cmd,shell=True)
    cmd2 = f"tabix -p vcf {outfile}"
    subprocess.run(cmd2,shell=True)
