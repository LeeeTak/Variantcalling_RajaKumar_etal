#!/usr/bin/python

import os,sys
import subprocess
import glob

genomei = "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/Plecuc1_hybrid_assembly.fasta.fai"

def reheader_sort(f,flag):
    if flag == "all":
        reheader = f.replace(".vcf.gz","_reheadered.vcf.gz")
        cmd1 = f"bcftools reheader -f {genomei} -o {reheader} {f}"
        subprocess.run(cmd1,shell=True)
        fsorted = reheader.replace(".vcf.gz","_sorted.vcf.gz")
        cmd2 = f"bcftools sort {reheader} -o {fsorted}"
        subprocess.run(cmd2,shell=True)
        cmd3 = f"tabix -p vcf {fsorted}"
        subprocess.run(cmd3,shell=True)
        return(fsorted)
    elif flag == "sorting":
        fsorted = f.replace(".vcf.gz","_sorted.vcf.gz")
        cmd2 = f"bcftools sort {f} -o {fsorted}"
        subprocess.run(cmd2,shell=True)
        cmd3 = f"mv {fsorted} {f}"
        subprocess.run(cmd3,shell=True)
        cmd4 = f"tabix -p vcf {f}"
        subprocess.run(cmd4,shell=True)

softclipped = ["150"]#["140","150"]
indeldir = "hybridgenome_first_INDELcall_rm_chimeric_multimapped"

for sc in softclipped:
    concatdir = f"hybridgenome_first_CONCATcall_rm_chimeric_multimapped_sc{sc}"
    subprocess.run(f"mkdir -p {concatdir}",shell=True)
    snpdir = f"hybridgenome_first_SNPcall_rm_chimeric_multimapped_sc{sc}"
    snpcalls = glob.glob(f"{snpdir}/Ath_G08_T3*.vcf")#*snp.vcf")
    for snp in snpcalls:
        indel = snp.replace(snpdir,indeldir).replace("_snp.vcf","_indel.vcf").replace(f"_sc{sc}","")
        concat = snp.replace(snpdir,concatdir).replace("_snp","")
        cmd1 = f"bgzip -c {snp} > {snp}.gz"
        cmd2 = f"bgzip -c {indel} > {indel}.gz"
        cmd3 = f"tabix -p vcf {snp}.gz"
        cmd4 = f"tabix -p vcf {indel}.gz"
        subprocess.run(cmd1,shell=True)
        subprocess.run(cmd2,shell=True)
        subprocess.run(cmd3,shell=True)
        subprocess.run(cmd4,shell=True)
        #reheader all the files
        snpfsorted = reheader_sort(f"{snp}.gz","all")
        indelfsorted = reheader_sort(f"{indel}.gz","all")

        cmd5 = f"bcftools concat -n {snpfsorted} {indelfsorted} -o {concat}.gz"
        if not os.path.isfile(f"{concat}.gz"):
            subprocess.run(cmd5,shell=True)
        #have to sort again
        reheader_sort(f"{concat}.gz","sorting")
