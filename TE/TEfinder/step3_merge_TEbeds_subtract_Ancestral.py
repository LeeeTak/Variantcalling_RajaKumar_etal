#!/usr/bin/python
import os,sys
import subprocess

fai = "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/Plecuc1_hybrid_assembly.fasta.fai"
ancestral = "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/TEfinder-1.0.1/Ancestral/TEinsertions.bed"
sorta = f"bedtools sort -faidx {fai} -i {ancestral} > Ancestral.sorted.bed"
subprocess.run(sorta,shell=True)
for c in ["Ath","mSly","Hvu","system"]:
    outfile = f"{c}.TE.union.bed"
    cmd = f"cat {c}*TE.bed | bedtools sort -faidx {fai} -i - | bedtools merge -i - > {outfile}"
    subprocess.run(cmd,shell=True)
    sortedbed = outfile.replace(".bed",".sorted.bed")
    sorts = f"bedtools sort -faidx {fai} -i {outfile} > {sortedbed}"
    subprocess.run(sorts,shell=True)
    subtracted = sortedbed.replace(".sorted.bed",".sorted.afilt.bed")
    subtract = f"bedtools subtract -a {sortedbed} -b Ancestral.sorted.bed > {subtracted}"
    subprocess.run(subtract,shell=True)
