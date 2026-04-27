#!/usr/bin/python
import os,sys
import subprocess
import glob

inputf = sorted(glob.glob("*context*txt.gz"))
for f in inputf:
    context = f.split("_")[0]
    samp = f.split("_")[2]+"_"+f.split("_")[3]
    bedgraph = f"{samp}_{context}_bedgraph"
    cmd = f"bismark2bedGraph --CX --dir . -o {bedgraph} {f}"
    #subprocess.run(cmd,shell=True)
    outf = f"{samp}_{context}_bedgraph.gz.bismark.cov.gz"
    chng = outf.replace("_bedgraph.gz","")
    mvcmd = f"mv {outf} {chng}"
    if not os.path.exists(f"{bedgraph}.gz"):
        subprocess.run(cmd,shell=True)
    if not os.path.exists(chng):
        subprocess.run(mvcmd,shell=True)


