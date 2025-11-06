#!/usr/bin/python
import os,sys
import subprocess
import glob

ref1 = "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid_with_Leifsonia/Plecuc1_hybrid_assembly.fasta"
ref2 = "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid_with_Leifsonia/Leifsonia_sp_T36S-04.fasta"

for infile in glob.glob("/netscratch/dep_psl/common/1323_Ram_to_Tak/concatenated/*_concat.fastq.gz"):
    sampname = "Ancestral"
    if sampname not in infile:
        sampa = infile.split("/")[-1].split("_")
        sampname = f"{sampa[0]}_{sampa[1]}_{sampa[2]}"
    if "toothpick" not in sampname:
        outhybrid = sampname+"_"+ref1.split("/")[-1].replace(".fasta",".fastq")
        outother = sampname+"_"+ref2.split("/")[-1].replace(".fasta",".fastq")
        o1 = outhybrid.replace(".fastq","_R1.fastq")
        o2 = outhybrid.replace(".fastq","_R2.fastq")
        if not os.path.exists(f"{o1}.gz"):
            r1 = infile
            r2 = infile.replace("R1","R2")
            cmd = f"bbsplit.sh ref={ref1},{ref2} in={r1} in2={r2} interleaved=False basename={sampname}_%.fastq"
            subprocess.run(cmd,shell=True)
            cmd_splitfile = f"reformat.sh in={outhybrid} out1={o1} out2={o2}"
            subprocess.run(cmd_splitfile,shell=True)
            compress = f"pigz -p 32 {sampname}_*R*.fastq"
            subprocess.run(compress,shell=True)
            rmcmd = f"rm {outhybrid} {outother}"
            subprocess.run(rmcmd,shell=True)
        else:
            print("file already processed")
