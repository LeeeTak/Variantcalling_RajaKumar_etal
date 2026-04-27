#!/usr/bin/python
import os,sys
import glob
import subprocess
bams = sorted(glob.glob("*pe.bam"))

for b in bams:
    if not os.path.exists(b.replace(".bam",".deduplicated.bam")):
        cmd = f"deduplicate_bismark {b}"
        subprocess.run(cmd,shell=True)
        #print(cmd)
