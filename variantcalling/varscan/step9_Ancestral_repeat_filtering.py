#!/usr/bin/python

import subprocess
import os,sys
import re
import glob
tool = "varscan"
suffix = f"{tool}_mapq60baseq30dp10varfreq0.1p0.05minread2.vcf.gz"
sys.path.insert(0, '/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/'+tool+'/scripts')
from further_variant_filtering import statsgen,get_rank,merge_chkrep_vcfs,getflank,filter_Ancestral
Qth = 30
DPth = 10
workingpath = "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/"+tool
ancestralcall = "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/varscan/hybridgenome_first_CONCATcall_rm_chimeric_multimapped_sc150/Ancestral_hybrid_rm_chim_multim_sc150_varscan_dp5_varfreq0.001p0.5minread1.vcf.gz"
softclipped = ["150"]#["140","150"]

for sc in softclipped:
    calledpath = f"{workingpath}/hybridgenome_first_CONCATcall_rm_chimeric_multimapped_sc{sc}"
    repeatfiltpath = f"{workingpath}/filt2_repeatfiltered_sc{sc}"
    subprocess.run("mkdir -p "+repeatfiltpath,shell=True)
    ancfiltpath = f"{workingpath}/filt1_Ancestralfiltered_sc{sc}"
    subprocess.run("mkdir -p "+ancfiltpath,shell=True)
    firstinputfiles = glob.glob(f"{calledpath}/Ath_G08_T3*{suffix}")
    firstinputs = ""
    if len(firstinputfiles) > 0:
        firstinputs = f"{calledpath}/*{suffix}"
        unzipping = "parallel bgzip -k -d {} ::: "+firstinputs
        print(1)
        subprocess.run(unzipping,shell=True)
    else:
        suffix2 = suffix.replace(".vcf.gz",".vcf")
        firstinputs = f"{calledpath}/*{suffix2}"
        zipping = "parallel bgzip -k {} ::: "+firstinputs
        print(2)
        subprocess.run(zipping,shell=True)

    #ancestralcall = f"{calledpath}/Ancestral_hybrid_rm_chim_multim_sc150_varscan_dp2_varfreq0.001p0.5minread1.vcf.gz"
    if not os.path.exists(f"{ancestralcall}.tbi"):
        indexing = f"parallel tabix -p vcf {{}} ::: {ancestralcall}"
        subprocess.run(indexing,shell=True)

    def plothist(infile):
        plothistcmd = "Rscript /netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/"+tool+"/scripts/hist.R "+infile
        subprocess.run(plothistcmd,shell=True)


    #statsfile = statsgen(calledpath,Qth,DPth,tool+"_Firstcall")
    #plothist(statsfile)
    filter_Ancestral(firstinputfiles,ancfiltpath,ancestralcall)
    #unzipping2 = "parallel bgzip -k -d {} ::: "+ancfiltpath+"/*.gz"
    #ancfiltstatsfile = statsgen(ancfiltpath,Qth,DPth,tool+"_Ancestralfilt")
    #plothist(ancfiltstatsfile)
    zipping = "parallel bgzip -k {} ::: "+ancfiltpath+"/*.vcf"
    subprocess.run(zipping,shell=True)
    indexing = "parallel tabix -p vcf {} ::: "+ancfiltpath+"/*.vcf.gz"
    subprocess.run(indexing,shell=True)
    merge_chkrep_vcfs(ancfiltpath,repeatfiltpath)
    #repfiltstats = statsgen(repeatfiltpath,Qth,DPth,tool+"_Repfilt")
    #plothist(repfiltstats)
