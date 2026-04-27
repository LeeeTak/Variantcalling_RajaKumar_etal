#!/usr/bin/python

import subprocess
from pathlib import Path
import glob

def run(cmd, cwd=None):
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)

outdir = "methylation_extract_ignore10bp"
subprocess.run(f"mkdir -p {outdir}",shell=True)

GENOME_FOLDER = "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/bismark_genome"
BAMS = sorted(glob.glob("G2T1_r3*deduplicated.bam"))
for bam in BAMS:
    run([
        "bismark_methylation_extractor",
        "--paired-end",
        "--ignore",str(10),
        "--ignore_r2",str(10),
        "--gzip",
        "--comprehensive",
        "--bedGraph",
        "--counts",
        "--parallel",str(24),
        "--cytosine_report",
        "--genome_folder", GENOME_FOLDER,
        "--output_dir", outdir,
        bam
    ])
