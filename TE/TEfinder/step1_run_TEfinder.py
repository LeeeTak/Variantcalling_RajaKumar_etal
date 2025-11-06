#!/usr/bin/python
import os,sys
import subprocess
import glob

bamdir = f"/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/preprocess_files_Leifsonia_split"
inputbam = sorted(glob.glob(f"{bamdir}/*_rm_scmap10_chimeric_multi.bam"))
reference = "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/Plecuc1_hybrid_assembly.fasta"
TE = "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/repeatmodeler/repeatmasker_out/Plecuc1_hybrid_TEs.gtf"
TElist = "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/repeatmodeler/repeatmasker_out/Plecuc1_hybrid_TElist.txt"
workdir = ""
for bam in inputbam:
    if "toothpick" not in bam:
        workdir = "Ancestral"
        bama = bam.split("/")[-1].split("_")
        if bama[0] != "Ancestral":
            workdir = f"{bama[0]}_{bama[1]}_{bama[2]}"
        cmd = f"./TEfinder -alignment {bam} -fa {reference} -gtf {TE} -te {TElist} -out GTF -intermed yes -picard picard.jar -workingdir {workdir}"
        subprocess.run(cmd,shell=True)
