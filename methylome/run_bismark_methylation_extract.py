#!/usr/bin/python

import subprocess
from pathlib import Path
import glob

def run(cmd, cwd=None):
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)


GENOME_FOLDER = "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/bismark_genome"
BAMS = sorted(glob.glob("*deduplicated.bam"))
for bam in BAMS:
    run([
        "bismark_methylation_extractor",
        "--paired-end",
        "--gzip",
        "--comprehensive",
        "--bedGraph",
        "--counts",
        "--parallel",str(24),
        "--cytosine_report",
        "--genome_folder", GENOME_FOLDER,
        bam
    ])
