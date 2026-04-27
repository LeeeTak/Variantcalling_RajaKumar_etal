#!/usr/bin/python
import os,sys

annots = {la[0]:la[1:4]+la[24:] for la in (l.strip().split("\t") for l in open("/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/snpEff/from_previous_rnaseq/withJGIname_hybridplecuc1_RNAseq_lfc_annotated_table_generenamed_notintranscript.txt","r"))}
dummy = ['NA']*len(annots["polished_locus_id"])
for f in open(sys.argv[1],"r"):
    f=f.strip()
    outf = open(f.replace("diff.txt","DEGs_onlyJGI.tsv"),"w+")
    i=0
    degs = {}
    outf.write("locusID"+"\t"+"log2FC"+"\t"+"FDR adjusted p-value"+"\t"+"\t".join(x for x in annots["polished_locus_id"])+"\n")
    for l in open(f,"r"):
        la = l.strip().split("\t")
        if i > 0:
            g = la[0]
            lfc = 0
            adjp = 1
            if la[2] != "NA":
                lfc = float(la[2])
            if la[6] != "NA":
                adjp = float(la[6])
            if abs(lfc) >= 1 and adjp <= 0.05:
                degs[la[0]] = (lfc,adjp)
        i=i+1
    s_degs = sorted(degs.items(), key=lambda x: x[1],reverse=True)
    for x in s_degs:
        if x[0] in annots:
            outf.write(x[0]+"\t"+str(x[1][0])+"\t"+str(x[1][1])+"\t"+"\t".join(annots[x[0]])+"\n")


