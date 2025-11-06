#!/usr/bin/python
import os,sys

conditions = ["Ath","mSly","Hvu","system"]
tra = ["T1","T2","T3","T4","T5"]
gen = ["G02","G04","G06","G08","G10","G12"]

for c in conditions:
    for g in gen:
        for t in tra:
            samp = f"{c}_{g}_{t}"
            bed = f"{samp}/TEinsertions.bed"
            bedout = f"bedfiles/{samp}_TE.bed"
            with open(bed,"r") as bedfile, open(bedout,"w+") as bedoutfile:
                for l in bedfile:
                    la=l.strip().split("\t")
                    if "track name" in la[0]:
                        bedoutfile.write(l.strip()+"\n")
                    else:
                        laa=la[6].split(";")
                        f=int(laa[0].replace("FR=",""))
                        r=int(laa[1].replace("RR=",""))
                        meanreads = (f+r)/2
                        if meanreads > 10:
                            if "weak" not in l.strip():
                                bedoutfile.write(l.strip()+"\n")
