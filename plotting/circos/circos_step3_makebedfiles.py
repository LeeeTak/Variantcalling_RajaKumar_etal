#!/usr/bin/python
import os,sys
import glob
import subprocess
import shlex

indir = "normalised/merged"
infiles = glob.glob(f"{indir}/*vcf.gz")
outdir = "bedfiles"
subprocess.run(f"mkdir -p {outdir}",shell=True)

def outname(infile,typ):
    outputname = outdir+"/"+infile.split("/")[-1].replace(".vcf.gz",f".{typ}.bed")
    return(outputname)

def points(inputfile):
    outputfile= outname(inputfile,"points")
    cmd = f"bcftools query -f '%CHROM\t%POS0\t%POS\n' {inputfile} | sort -k1,1 -k2,2n -u > {outputfile}"
    subprocess.run(cmd,shell=True)

def intervals(inputfile):
    outputfile= outname(inputfile,"intervals")
    cmd = f"bcftools query -f '%CHROM\t%POS\t%REF\n' {inputfile} | awk 'BEGIN{{OFS=\"\t\"}}{{start=$2-1; end=$2+length($3)-1; print $1, start, end}}' | sort -k1,1 -k2,2n -u > {outputfile}"
    subprocess.run(cmd,shell=True)

def pointsNcounts(inputfile):
    outputfile= outname(inputfile,"pointsNcounts")
    cmd = r"""bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' {infile} \
            | awk 'BEGIN{{OFS="\t"}} {{ n_alt=0; for(i=5;i<=NF;i++){{ gt=$i; if(gt!="." && gt!="./." && gt!=".|." && gt!="0/0" && gt!="0|0") n_alt++; }} print $1, $2-1, $2, n_alt }}' \
            | sort -k1,1 -k2,2n > {outfile}""".format(
    infile=shlex.quote(inputfile),
    outfile=shlex.quote(outputfile)
    )
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


for f in infiles:
    #points(f)
    intervals(f)
    #pointsNcounts(f)
