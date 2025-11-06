#!/usr/bin/python
import os,sys

i=0
for l in open(sys.argv[1],"r"):
    l=l.strip()
    if not l.startswith("##"):
        la=l.strip().split("\t")
        if la[0] == "#CHROM":
            print("variant\t"+"\t".join(la[9:]))
        else:
            i+=1
            variant = la[0]+":"+la[1]
            #if i < 10:
            #    variant="variant000"+str(i)
            #elif i < 100:
            #    variant="variant00"+str(i)
            #elif i < 1000:
            #    variant="variant0"+str(i)
            #else:
            #    variant="variant"+str(i)
            per = [variant]
            for v in la[9:]:
                if v.startswith('.'):
                    per.append(str(0))
                else:
                    per.append(v.split(":")[6].replace("%",""))
            print("\t".join(per))




