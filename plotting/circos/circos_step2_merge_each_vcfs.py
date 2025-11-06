#!/usr/bin/python
import os,sys
import subprocess

workdir = "normalised"
outd = f"{workdir}/merged"
subprocess.run(f"mkdir -p {outd}",shell=True)
for c in ["Ath","mSly","Hvu","system"]:
    out = f"{outd}/{c}_merged.vcf.gz"
    cmd=f"bcftools merge -m none -Oz -o {out} {workdir}/{c}/*.norm.vcf.gz"
    subprocess.run(cmd,shell=True)
    cmd2 = f"tabix {out}"
    subprocess.run(cmd2,shell=True)

