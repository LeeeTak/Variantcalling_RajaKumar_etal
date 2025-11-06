#i!/usr/bin/Rscript
library(plot3D)
pca3d_vector <- function(table, outstem = "pca3d",
                         label_quantile = 0.75,  #quartile for labeling
                         theta = 50, phi = 25) {
    stopifnot(is.matrix(table) || is.data.frame(table))
    pca <- prcomp(table)#, center = TRUE, scale. = TRUE)
    pc1 <- pca$x[,1]
    pc2 <- pca$x[,2]
    pc3 <- pca$x[,3]
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    percentVar <- round(100 * percentVar)
    pcadat <- data.frame(PC1 = pc1, PC2 = pc2, PC3=pc3)
    split_names <- strsplit(rownames(matdat), "_")
    pcadat$condition <- sapply(split_names, function(x) x[1])
    pcadat$generation <- sapply(split_names, function(x) sub("G0", "", x[2]))
    splitnames2 <- pcadat$generation
    pcadat$generation <- as.numeric(sub("G","",splitnames2))
    pcadat$generation <- ordered(pcadat$generation)
    pcadat$trajectory <- sapply(split_names, function(x) sub("trajectory", "T", x[3]))
    pcadat["Ancestral",]$generation <- 6
    pcadat["Ancestral",]$trajectory <- "Ancestral"
    pcadat$label <- rownames(pcadat)
    anc_idx <- which(pcadat$label == "Ancestral")
    color_map <- c("system"="#999999","Ath"="#008b45","Hvu"="#8b864e","mSly"="#8b3626","Ancestral"="black")
    point_col <- ifelse(pcadat$condition %in% names(color_map),
			color_map[pcadat$condition], "#555555")
    names(point_col) <- rownames(pcadat)
    text_col <- point_col[!names(point_col) %in% c("Ancestral")]
    label_show <- rep("", nrow(pcadat))
    anc_idx <- which(pcadat$label == "Ancestral")
    ax <- pcadat$PC1[anc_idx]
    ay <- pcadat$PC2[anc_idx]
    az <- pcadat$PC3[anc_idx]
    pcadat$dist_from_anc <- sqrt((pcadat$PC1 - ax)^2 + (pcadat$PC2 - ay)^2 + (pcadat$PC3 - az)^2)
    thr <- stats::quantile(pcadat$dist_from_anc, 0.75, na.rm = TRUE)
    label_show <- ifelse(
			 pcadat$dist_from_anc >= thr | rownames(pcadat) == "Ancestral", pcadat$label, ""
    )

    xl <- paste0("PC1: ", percentVar[1], "% variance")
    yl <- paste0("PC2: ", percentVar[2], "% variance")
    zl <- paste0("PC3: ", percentVar[3], "% variance")
    draw_plot <- function() {
    scatter3D(x = pcadat$PC1, y = pcadat$PC2, z = pcadat$PC3,
	      xlab = xl, ylab = yl, zlab = zl,
	      col = point_col, pch = 16, cex = 1,colvar=NULL,
	      colkey = FALSE, theta = theta, phi = phi,
	      bty = "b2", # filled panels; change to "b2" for simple box
	      ticktype = "detailed")

    #just overplotting Ancestral to ensure visibility
    if (length(anc_idx) == 1) {
	points3D(pcadat$PC1[anc_idx], pcadat$PC2[anc_idx], pcadat$PC3[anc_idx],
		 add = TRUE, pch = 16, cex = 3, col = "black",colvar=NULL)
    }

    groups <- split(pcadat, interaction(pcadat$condition, pcadat$trajectory, drop = TRUE))
    for (g in groups) {
	if (nrow(g) < 1) next
	g <- g[order(g$generation), , drop = FALSE]
	if (length(anc_idx) == 1) {
	    x <- c(pcadat$PC1[anc_idx], g$PC1)
	    y <- c(pcadat$PC2[anc_idx], g$PC2)
	    z <- c(pcadat$PC3[anc_idx], g$PC3)
	} else {
	    x <- g$PC1; y <- g$PC2; z <- g$PC3
	}
	col_line <- if (g$condition[1] %in% names(color_map)) color_map[g$condition[1]] else "#555555"
	lines3D(x, y, z, add = TRUE, col = col_line, lwd = 0.1,colvar=NULL)
    }

    dx <- 0.02 * diff(range(pcadat$PC1, na.rm = TRUE))
    dy <- 0.02 * diff(range(pcadat$PC2, na.rm = TRUE))
    idx_lab <- which(label_show != "")
    idx_anc <- intersect(idx_lab, which(pcadat$label == "Ancestral"))
    idx_oth <- setdiff(idx_lab, idx_anc)
    text_col <- point_col[names(point_col) %in% label_show[idx_oth]]

    if (length(idx_oth) > 0) {
	text3D(pcadat$PC1[idx_oth]+dx, pcadat$PC2[idx_oth]+dy, pcadat$PC3[idx_oth],
	       labels = label_show[idx_oth],colvar=NULL,colkey=FALSE,
         add = TRUE, cex = 0.8, col = text_col)
    }
    if (length(idx_anc) == 1) {
	text3D(pcadat$PC1[idx_anc]+dx, pcadat$PC2[idx_anc]+dy, pcadat$PC3[idx_anc],
	       labels = label_show[idx_anc],colvar=NULL,
	       add = TRUE, cex = 1.2, col = "black")
    }
    leg_levels <- unique(pcadat$condition)
    leg_cols <- ifelse(leg_levels %in% names(color_map), color_map[leg_levels], "#555555")
    legend("topright", legend = leg_levels, col = leg_cols, pch = 16, pt.cex = 1.2, bty = "n")
  }

    pdf(paste0(outstem, ".pdf"), width = 16, height = 16)
    draw_plot()

    svg(paste0(outstem, ".svg"), width = 16, height = 16)
    draw_plot()
    dev.off()

}
dat<-read.table("/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/varscan/variants_heatmap/repeatfilt_sc150_splitgenome_VERYrelaxedAnc_strictercall_readper_table_orderedbypos.txt",header=T,row.names=1)
dat$Ancestral <- rep(0,length(dat[1]))
blacklists <- c("scaffold_3_pilon:5389696","scaffold_6_pilon:84121","scaffold_3_pilon:5788746","scaffold_3_pilon:4067561","scaffold_5_pilon:3302727","scaffold_6_pilon:1734032","scaffold_3_pilon:5608240")
dat <- dat[!(rownames(dat) %in% blacklists),]
tdat <- t(dat)
matdat <- matrix(as.numeric(tdat),ncol=length(colnames(tdat)))
rownames(matdat) <- rownames(tdat)
colnames(matdat) <- colnames(tdat)

matdat_binary <- matdat
matdat_binary[is.na(matdat_binary)] <- 0
matdat_binary[matdat_binary != 0] <- 1

for (tht in seq(0,360,30)){
    outfile = paste0("PCA_3D_by_presence_noblacklist_colorbycond_namefurther_th",tht)
    pca3d_vector(matdat_binary, outstem = outfile,theta=tht)
}
