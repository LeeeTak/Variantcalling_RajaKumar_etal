#!/usr/bin/Rscript
#script written by Tak Lee
#usage: Rscript compheatmap_forvariants.R

library(RColorBrewer)
library(dplyr)
library(stringr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)
perthreshold <- 10

#read data
infile <- "repeatfilt_sc150_splitgenome_VERYrelaxedAnc_strictercall_readper_table_orderedbypos.txt"
rawdat <- read.table(infile,sep='\t',header=T,row.names=1,check.names=F)
dat <- rawdat[c("Ath_G02_T1","Ath_G04_T1","Ath_G06_T1","Ath_G08_T1","Ath_G10_T1","Ath_G12_T1","Ath_G02_T2","Ath_G04_T2","Ath_G06_T2","Ath_G08_T2","Ath_G10_T2","Ath_G12_T2","Ath_G02_T3","Ath_G04_T3","Ath_G06_T3","Ath_G08_T3","Ath_G10_T3","Ath_G12_T3","Ath_G02_T4","Ath_G04_T4","Ath_G06_T4","Ath_G08_T4","Ath_G10_T4","Ath_G12_T4","Ath_G02_T5","Ath_G04_T5","Ath_G06_T5","Ath_G08_T5","Ath_G10_T5","Ath_G12_T5","mSly_G02_T1","mSly_G04_T1","mSly_G06_T1","mSly_G08_T1","mSly_G10_T1","mSly_G12_T1","mSly_G02_T2","mSly_G04_T2","mSly_G06_T2","mSly_G08_T2","mSly_G10_T2","mSly_G12_T2","mSly_G02_T3","mSly_G04_T3","mSly_G06_T3","mSly_G08_T3","mSly_G10_T3","mSly_G12_T3","mSly_G02_T4","mSly_G04_T4","mSly_G06_T4","mSly_G08_T4","mSly_G10_T4","mSly_G12_T4","mSly_G02_T5","mSly_G04_T5","mSly_G06_T5","mSly_G08_T5","mSly_G10_T5","mSly_G12_T5","Hvu_G02_T1","Hvu_G04_T1","Hvu_G06_T1","Hvu_G08_T1","Hvu_G10_T1","Hvu_G12_T1","Hvu_G02_T2","Hvu_G04_T2","Hvu_G06_T2","Hvu_G08_T2","Hvu_G10_T2","Hvu_G12_T2","Hvu_G02_T3","Hvu_G04_T3","Hvu_G06_T3","Hvu_G08_T3","Hvu_G10_T3","Hvu_G12_T3","Hvu_G02_T4","Hvu_G04_T4","Hvu_G06_T4","Hvu_G08_T4","Hvu_G10_T4","Hvu_G12_T4","Hvu_G02_T5","Hvu_G04_T5","Hvu_G06_T5","Hvu_G08_T5","Hvu_G10_T5","Hvu_G12_T5","system_G02_T1","system_G04_T1","system_G06_T1","system_G08_T1","system_G10_T1","system_G12_T1","system_G02_T2","system_G04_T2","system_G06_T2","system_G08_T2","system_G10_T2","system_G12_T2","system_G02_T3","system_G04_T3","system_G06_T3","system_G08_T3","system_G10_T3","system_G12_T3","system_G02_T4","system_G04_T4","system_G06_T4","system_G08_T4","system_G10_T4","system_G12_T4","system_G02_T5","system_G04_T5","system_G06_T5","system_G08_T5","system_G10_T5","system_G12_T5")]
dat[dat < perthreshold] <- 0
dat <- dat[rowSums(dat) > 0,]
blacklists <- c("scaffold_3_pilon:5389696","scaffold_6_pilon:84121","scaffold_3_pilon:5788746","scaffold_3_pilon:4067561","scaffold_5_pilon:3302727","scaffold_6_pilon:1734032","scaffold_3_pilon:5608240")
dat_blacklisted <- dat[!(rownames(dat) %in% blacklists),]

for (chk in c("all","blacklist"))
{
    if (chk == "all"){
	dat <- dat
    }
    else if (chk == "blacklist"){
	     dat <- dat_blacklisted
    }
    impacts <- read.csv("/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/snpEff/varscan/varscan_variant_impacts.txt",sep="\t")
    high <- impacts[impacts$impact == "HIGH",]$variant
    moderate <- impacts[impacts$impact == "MODERATE",]$variant
    modifier <- impacts[impacts$impact == "MODIFIER",]$variant
    low <- impacts[impacts$impact == "LOW",]$variant
    intergenic <- rownames(dat[!rownames(dat) %in% c(high,moderate,modifier,low),])
    cn <- colnames(dat)
    parts <- do.call(rbind, strsplit(cn, "_"))
    cond_raw <- parts[,1]
    gen <- parts[,2]
    traj <- parts[,3]
    cond <- sub("^system$", "Sys", cond_raw)  # map 'system' -> 'Sys' for legend

    varcount <- matrix(nrow=length(colnames(dat)),ncol=6)
    colnames(varcount) <- c("INTERGENIC","MODIFIER","LOW","MODERATE","HIGH","ALL")
    rownames(varcount) <- colnames(dat)

    make_mat <- function(mat,glist,impct){
	subdat <- drop_na(dat[glist,])
	subdat[subdat > 0] <- 1
	mat[,impct] <- colSums(subdat)
	return(mat)
    }

    varcount <- make_mat(varcount,high,"HIGH")
    varcount <- make_mat(varcount,moderate,"MODERATE")
    varcount <- make_mat(varcount,low,"LOW")
    varcount <- make_mat(varcount,modifier,"MODIFIER")
    varcount <- make_mat(varcount,intergenic,"INTERGENIC")
    varcount <- make_mat(varcount,c(high,moderate,low,modifier,intergenic),"ALL")


    colsplit <- c()
    coltitle <- c()
    colgap <- c()
    samps <- 20
    reps <- 6
    j<-1
    for (i in seq(1,samps,1)){
	colsplit <- c(colsplit,rep(i,reps))
	if (j == 5){
	    colgap <- c(colgap,3.0)
	}
	else{
	    colgap <- c(colgap,1.5)
	}
	if (j > 5){
	    j<-1
	}
	coltitle <- c(coltitle,j)
	j<-j+1
    }
    colgap <- head(colgap,-1)

    # base colors per condition
    base_cols <- c(Ath = "#008b45", mSly = "#8b3626", Hvu = "#8b864e", Sys = "#999999")
    gen_levels <- c("G02","G04","G06","G08","G10","G12")

    # build per-condition generation shades (white -> base; G12 == base color), don't need shades anymore. Just shades of grey
    col_grad <- c()
    for (k in names(base_cols)) {
	pal <- c("grey45","grey60","grey75","grey80","grey90","#FAFAFA") 
	#pal <- colorRampPalette(c("white", base_cols[k]))(length(gen_levels))
	#pal <- colorRampPalette(c("white", base_cols[k]))(length(gen_levels))
	names(pal) <- paste(k, gen_levels, sep = "_")  # e.g., "Ath_G02", ...
	col_grad <- c(col_grad, pal)
    }
    traj_num <- sub("^T", "", traj)

    vartype <- read.csv("snp_ins_dels.txt",sep="\t",row.names=1)
    vartypefilt <- vartype[rownames(dat),]

    #RNAseq data
    rnaseqdat <- read.csv("/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/snpEff/varscan/DEs_prioritizing_dist.txt",sep="\t",header=T,check.names=F,row.names=1)
    rnaseqdatfilt <- rnaseqdat[rownames(dat),]
    rownames(rnaseqdatfilt) <- rownames(dat)
    rnaseqdatfilt[is.na(rnaseqdatfilt)] <- 0
    rnaseqborder <- rownames(rnaseqdatfilt)[rowSums(rnaseqdatfilt) > 0]
    rnaseqdatfilt[rnaseqdatfilt > 0] <- "up"
    rnaseqdatfilt[rnaseqdatfilt < 0] <- "down"
    rnaseqdatfilt[rnaseqdatfilt == 0] <- "notDE"

    #annotations
    ha <- HeatmapAnnotation(
			    variant_count = anno_barplot(
							 #varcount,
							 varcount[,"ALL"],
							 gp=gpar(fill="black",col=NA),#c("#c3b091","#1e90ff","#00ff00","#ffa500","#ff1493"),col=NA),
							 height=unit(18,"mm")),
			    generation = anno_simple(
						     paste(cond, gen, sep = "_"),
						     col= col_grad,
						     border = TRUE,
						     gp = gpar(col = "black", lwd = 0.1),
						     height = unit(3, "mm")
						     ),
			    annotation_name_side = "left",
			    annotation_name_gp   = gpar(fontsize = 8)
    )
    ba <- HeatmapAnnotation(
			    Trajectory=anno_block(gp=gpar(fill="black"),
						  labels = coltitle,
						  labels_gp = gpar(col="white",fontsize=8),
						  ),
			    annotation_name_side = "right",
			    height = unit(3,"mm")
			    )
    sc_annotation <- rowAnnotation(scaffold=anno_block(gp=gpar(fill="black"),
							  labels = c("1","2","3","4","5","6","7"),
							  labels_gp = gpar(col="white",fontsize=8),
							  width=unit(3,"mm"),
							  ),
				   annotation_label="scaffold",
				   annotation_name_side="bottom"
    )

    type_annotation <- rowAnnotation(Ath= anno_simple(
							rnaseqdatfilt$AthDE,
							col=c("up"="black","notDE"="white","down"="white"),
							width = unit(2,"mm"),
							gp = gpar(col = "black", lwd = 0.1)
							),
				     ` `  = anno_empty(width = unit(0.5, "mm"),border=F),
				     mSly= anno_simple(
							rnaseqdatfilt$mSlyDE,
							col=c("up"="black","notDE"="white","down"="white"),
							width = unit(2,"mm"),
							gp = gpar(col = "black", lwd = 0.1)
							),
				     `  `  = anno_empty(width = unit(0.5, "mm"),border=F),
				     Hvu= anno_simple(
							rnaseqdatfilt$HvuDE,
							col=c("up"="black","notDE"="white","down"="white"),
							width = unit(2,"mm"),
							gp = gpar(col = "black", lwd = 0.1)
							),
				     `   `  = anno_empty(width = unit(0.5, "mm"),border=F),
				     type= anno_simple(
						       vartypefilt,
						       col=c("SNP"="#b788c2","Ins"="#7abf95","Del"="#1A5b77"),
						       width = unit(2,"mm"),
						       gp = gpar(col = "black", lwd = 0.1)
						       ),
				     annotation_name_gp   = gpar(fontsize = 7)
    )


    rn <- rownames(dat)
    scaf_id <- sub("^scaffold_(\\d+)_pilon:.*", "\\1", rn)
    rowsplit <- factor(scaf_id, levels = as.character(sort(unique(as.numeric(scaf_id)))))
    rowtitle <- c()
    for (i in seq(1,length(unique(scaf_id)),1)){
	rowtitle <- c(rowtitle,"")
    }

    #colors <- colorRamp2(c(0,50,100),c("grey45", "yellow","red"))
    colors <- colorRamp2(c(0,100),c("white", "#c3b091"))
    pinks <- colorRamp2(c(0, 100), c("white", "#ff1493"))
    pinks2 <- colorRamp2(c(0, 100), c("#ededed", "#ff1493"))
    limes <- colorRamp2(c(0,100),c("white","#00ff00"))
    limes2 <- colorRamp2(c(0, 100), c("#ededed", "#00ff00"))
    skyblue <- colorRamp2(c(0,100),c("white","#1e90ff"))
    skyblue2 <- colorRamp2(c(0,100),c("#ededed","#1e90ff"))
    oranges <- colorRamp2(c(0,100),c("white","#ffa500"))
    oranges2 <- colorRamp2(c(0,100),c("#ededed","#ffa500"))
    outfile <- paste(gsub(".txt","_heatmap_",infile),chk,".pdf",sep="")
    pdf(outfile,height=9,width=16)
    heatmap <-Heatmap(
	    as.matrix(dat),
	    col=colors,
	    cluster_columns = F,
	    cluster_rows = F,
	    clustering_distance_rows = "euclidean",
	    clustering_method_rows= "ward.D",
	    #this is used to separate the dendrogram based on the best number of clusters
	    show_row_names = F,
	    show_column_names = F,
	    column_split = colsplit,
	    column_title = paste("variant heatmap sorted by position, filtered by ",perthreshold,"% allele frequency",sep=""),
	    row_split = rowsplit,
	    row_title = "Variants",#rowtitle,
	    rect_gp = gpar(col = "#E5E4E2", lwd = 0.5),
	    top_annotation = ha,
	    bottom_annotation = ba,
	    left_annotation = sc_annotation,
	    right_annotation = type_annotation,
	    gap = unit(1.5, "mm"),
	    column_gap = unit(colgap,"mm"),
	    #row_km = 7,
	    use_raster = F,
	    border=T,
	    show_heatmap_legend = F,
	    #heatmap_legend_param = list(title="percentage read support"),
	    cell_fun = function(j, i, x, y, w, h, col) {
		rn <- rownames(dat)[i]
		bcol <- "#E5E4E2"
		linew <- 0.5
		if (rn %in% high) {
		    v <- dat[i, j]
		    if (rn %in% rnaseqborder)
		    {
			grid.rect(x,y,w,h, gp = gpar(fill = pinks(v), col = bcol,lwd=linew))
		    }
		    else
		    {
			grid.rect(x,y,w,h, gp = gpar(fill = pinks(v), col = bcol,lwd=linew))
		    }
		}
		if (rn %in% moderate){
		    v <- dat[i,j]
		    if (rn %in% rnaseqborder)
		    {
			grid.rect(x,y,w,h,gp = gpar(fill = oranges(v), col = bcol,lwd=linew))
		    }
		    else
		    {
			grid.rect(x,y,w,h,gp = gpar(fill = oranges(v), col = bcol,lwd=linew))
		    }
		}
		if (rn %in% low){
		    v <- dat[i,j]
		    if (rn %in% rnaseqborder)
		    {
			grid.rect(x,y,w,h,gp = gpar(fill = limes(v), col = bcol,lwd=linew))
		    }
		    else
		    {
			grid.rect(x,y,w,h,gp = gpar(fill = limes(v), col = bcol,lwd=linew))
		    }
		}
		if (rn %in% modifier){
		    v <- dat[i,j]
		    if (rn %in% rnaseqborder)
		    {
			grid.rect(x,y,w,h,gp = gpar(fill = skyblue(v), col = bcol,lwd=linew))
		    }
		    else
		    {
			grid.rect(x,y,w,h,gp = gpar(fill = skyblue(v), col = bcol,lwd=linew))
		    }
		}
	    },
	    )
    blacks = colorRamp2(c(0,100),c("white", "black"))
    lg_per <- Legend(col_fun=blacks,title="allele frequency",direction="horizontal")
    lg_impact <- Legend(labels=c("High","Moderate","Low","Modifier","intergenic"),title = "strongest variant impact on genes",legend_gp = gpar(fill= c("#ff1493","#ffa500","#00ff00","#1e90ff","#c3b091")))
    #lg_condition <- Legend(labels=c("Ath","mSly","Hvu","Sys"),title="conditions",legend_gp=gpar(fill=c("#008b45","#8b3626","#8b864e","#999999")))
    lg_gen <- Legend(labels=c("2","4","6","8","10","12"),title="generation",legend_gp=gpar(fill=c("grey45","grey60","grey75","grey80","grey90","#FAFAFA")),border="black",direction="horizontal")
    lg_type <- Legend(labels=c("SNP","Ins","Del"),title="type",legend_gp=gpar(fill=c("#b788c2","#7abf95","#1A5b77")))
    lg_DE <- Legend(labels=c("induced","not induced"),title="differential expression",legend_gp=gpar(fill=c("black","white")),border="black")
    draw(heatmap,annotation_legend_list = list(lg_per,lg_impact,lg_gen,lg_DE,lg_type),padding = unit(c(10, 10, 10, 10), "mm"),ht_gap = unit(0.25, "cm"))
    dev.off()
}
