#!/usr/bin/python

import subprocess
import os,sys
import re
import glob
tool = "haplotypecaller"
sys.path.insert(0, '/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/'+tool+'/accessory_scripts')
from further_variant_filtering import statsgen,get_rank,merge_chkrep_vcfs,getflank,filter_Ancestral

Qth = 30
DPth = 10
workingpath = "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/"+tool
calledpath = workingpath+"/hybridgenome_first_variantcall_rm_softclipped130_chimeric_multimapped"
repeatfiltpath = workingpath+"/filt2_repeatfiltered"
subprocess.run("mkdir -p "+repeatfiltpath,shell=True)
ancfiltpath = workingpath+"/filt1_Ancestralfiltered"
subprocess.run("mkdir -p "+ancfiltpath,shell=True)
firstinputfiles = glob.glob(calledpath+"/*_"+tool+".vcf.gz")
firstinputs = ""
if len(firstinputfiles) > 0:
    firstinputs = calledpath+"/*_"+tool+".vcf.gz"
    unzipping = "parallel bgzip -k -d {} ::: "+firstinputs
    subprocess.run(unzipping,shell=True)
else:
    firstinputs = calledpath+"/*_"+tool+".vcf"
    zipping = "parallel bgzip -k {} ::: "+firstinputs
    subprocess.run(zipping,shell=True)

ancestralcall = calledpath+"/5808_A_hybrid_rm_scmap130_chim_multi_haplotypecaller.vcf.gz"


def plothist(infile):
    plothistcmd = "Rscript /netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/"+tool+"/scripts/hist.R "+infile
    subprocess.run(plothistcmd,shell=True)

#statsfile = statsgen(calledpath,Qth,DPth,tool+"_Firstcall")
#plothist(statsfile)
filter_Ancestral(calledpath,ancfiltpath,ancestralcall)
unzipping2 = "parallel bgzip -k -d {} ::: "+ancfiltpath+"/*.gz"
#ancfiltstatsfile = statsgen(ancfiltpath,Qth,DPth,tool+"_Ancestralfilt")
#plothist(ancfiltstatsfile)
zipping = "parallel bgzip -k {} ::: "+ancfiltpath+"/*.vcf"
subprocess.run(zipping,shell=True)
indexing = "parallel tabix -p vcf {} ::: "+ancfiltpath+"/*.vcf.gz"
subprocess.run(indexing,shell=True)
merge_chkrep_vcfs(ancfiltpath,repeatfiltpath)
repfiltstats = statsgen(repeatfiltpath,Qth,DPth,tool+"_Repfilt")
plothist(repfiltstats)
