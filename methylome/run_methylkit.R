#!/usr/bin/Rscript
library(methylKit)
library(genomation)
inputpath <- "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/methylome/bismark/methylation_extract_ignore10bp"
outdir <- "methylkit_differentialmeth"
genebed <- readTranscriptFeatures("/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/funannotate/predict_results/Plectosphaerella_cucumerina_Plecuc1.renamed.bed",remove.unusual=FALSE)
sampsNreps <- function(fl){
    bn <- basename(fl)
    parts <- strsplit(bn, "_")
    samples    <- vapply(parts, `[`, character(1), 1)
    replicates <- vapply(parts, `[`, character(1), 2)
    return(paste(samples,replicates,sep="_"))
}

allsamps <- sampsNreps(filelist)
#meth <- methRead(as.list(filelist),
#           sample.id=as.list(allsamps),
#           assembly="plecuc1",
#           treatment=c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4),
#           context="CpG",
#           mincov = 10,
#	   pipeline='bismarkCoverage'
#           )
#meth_merged <- unite(meth, destrand=FALSE)
#print(meth)

#pdf("methylome_clustering_correlation.pdf",width=15,height=15)
#clusterSamples(meth_merged,dist="correlation",method="ward.D",plot=TRUE)
#dev.off()

#pdf("methylome_clustering_euclidean.pdf",width=15,height=15)
#clusterSamples(meth_merged,dist="euclidean",method="ward.D",plot=TRUE)
#dev.off()

#pdf("methylome_PCA.pdf",width=9,height=9)
#PCASamples(meth_merged)
#dev.off()

context <- c("CpG","CHH","CHG")
samples <-c("G2T1","G2T5","G12T1","G10T5")
coverages <- c(15,30,40)
for (cov in coverages){
    for (s in samples){
	collated_table <- data.frame()
	tablename <- paste(outdir,paste(s,"_differential_methylation_cov",cov,".tsv",sep=""),sep="/")
	for (c in context){
	    apatt <- paste("^Anc_r\\d+_",c,"\\.bismark\\.cov\\.gz$",sep="")
	    Anc <- list.files(path=inputpath,
			      pattern = apatt,
			      full.names = TRUE)
	    cpatt <- paste("^",s,"_r\\d+_",c,"\\.bismark\\.cov\\.gz$",sep="")
	    comp <- list.files(path=inputpath,
			      pattern = cpatt,
			      full.names = TRUE)
	    infiles <- c(Anc,comp)
	    compsamps <- sampsNreps(infiles)
	    meth_comp <- methRead(as.list(infiles),
		   sample.id=as.list(compsamps),
		   assembly="plecuc1",
		   treatment=c(0,0,0,1,1,1),
		   context=c,
		   mincov = cov,
		   pipeline='bismarkCoverage'
		   )
	    comp_merged <- unite(meth_comp,destrand=FALSE)
	    diff <- calculateDiffMeth(comp_merged,adjust="BH")
	    #diff_filtered <- getMethylDiff(diff,difference=2,pvalue=1,type="all")
	    diffAnn <- annotateWithGeneParts(as(diff,"GRanges"),genebed)
	    members <- getMembers(diffAnn)          # columns: prom exon intron
	    m <- members
	    m[is.na(m)] <- 0
	    part_multi <- apply(m, 1, function(x) {
				    hits <- names(x)[x > 0]
				    hits <- sub("^prom$", "promoter", hits)      # nicer label
				    if (length(hits) == 0) "intergenic" else paste(hits, collapse = ",")
		   })
	    tss <- getAssociationWithTSS(diffAnn)
	    tss$gene_parts <- part_multi[tss$target.row]
	    tss$context <- rep(c,dim(tss)[1])
	    tss$qvalue <- diff$qvalue
	    tss$methdiff <- diff$meth.diff
	    tss$position <- diff$start
	    tss$pvalue <- diff$pvalue
	    tss$chr <- diff$chr
	    collated_table <- rbind(collated_table,tss[c("chr","position","context","feature.name","feature.strand","gene_parts","methdiff","pvalue","qvalue")])
	    }
	outtable <- collated_table[order(collated_table$pvalue),]
	write.table(outtable,file=tablename,row.names=F,sep="\t",quote=F)
    }
}

