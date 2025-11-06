#!/usr/bin/Rscript

library(ggpubr)
library(ggplot2)
Args <- commandArgs(TRUE)

df <- read.csv(Args[1],header=T,sep="\t")
df$value <- as.numeric(df$value)
for (i in unique(df$id)){
    preprefix <- gsub(".txt","",Args[1])
    prefix <- gsub("QUAL_DP_","",preprefix)
    outf <- paste(prefix,i,"hist.pdf",sep="_")
    pdf(outf,width=12,height=9)
    dat <- df[df$id == i,]
    b <- dim(dat)[1]
    histplot <- ggplot(dat, aes(x=value,fill=id)) +
    geom_histogram(bins=150)
    if (i == "DP") {
	histwline <- histplot + geom_vline(aes(xintercept=10))
    } else if (i == "QUAL") {
	histwline <- histplot + geom_vline(aes(xintercept=30)) }
    print(histwline)
    dev.off()
}
