#!/usr/bin/python
import os,sys
import glob
import subprocess

f = sorted(glob.glob("G2T1_r3*R1*.fastq.gz"))
genomedir = "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/bismark_genome"
for r1 in f:
    r2 = r1.replace("R1","R2")
    prfx = r1.split("_")[0]+"_"+r1.split("_")[1]
    cmd = f"bismark --parallel 24 --non_directional --genome {genomedir} -1 {r1} -2 {r2}"
    if not os.path.exists(r1.replace(".fastq.gz","_bismark_bt2_PE_report.txt")):
        subprocess.run(cmd,shell=True)
    #print(cmd)


