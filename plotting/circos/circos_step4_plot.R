#!/usr/bin/Rscript
library("circlize")

fai_path <- "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/Plecuc1_hybrid_assembly.fasta.fai"
indir <- "bedfiles"
bed_ath    <- paste0(indir,"/Ath_merged.intervals.bed")
bed_msly   <- paste0(indir,"/mSly_merged.intervals.bed")
bed_hvu    <- paste0(indir,"/Hvu_merged.intervals.bed")
bed_system <- paste0(indir,"/system_merged.intervals.bed")

bed_genes <- "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/funannotate/predict_results/Plectosphaerella_cucumerina_Plecuc1.renamed.scaffoldrenamed_geneonly.bed"
bed_repeats <- "/netscratch/dep_psl/grp_rgo/taklee/Ram_variantcalling/hybrid/repeatmodeler/repeatmasker_out/Plecuc1_hybrid_assembly.fasta.out.bed"
out_pdf <- "circos_variants.pdf"

# Colors
cols <- c(ath = "#008b45", msly = "#8b3626", hvu = "#8b864e", system = "grey60")

read_bed3 <- function(path) {
    x <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
		  quote = "", comment.char = "")
    if (ncol(x) < 3) stop("BED must have at least 3 columns: chr, start, end: ", path)
    colnames(x)[1:3] <- c("chr","start","end")
    x$start <- as.integer(x$start)
    x$end   <- as.integer(x$end)
    # drop malformed rows
    x <- x[!is.na(x$start) & !is.na(x$end) & x$end > x$start, ]
}

add_tick_track <- function(df, col, track_height = 0.03,label=NULL,label.col=NULL) {
    circos.genomicTrack(df[, c("chr","start","end")],
			ylim = c(0, 1),
			track.height = track_height,
			bg.border = "black",
			bg.lwd = 0.3,
			panel.fun = function(region, value, ...) {
			    circos.genomicRect(region,
					       ybottom = 0, ytop = 1,col = col, border =col,lwd=0.0)
			    if (!is.null(label)) {
				circos.text(CELL_META$xcenter, 0.5, labels = label,facing = "bending", niceFacing = TRUE,
					    col = if (is.null(label.col)) col else label.col,cex = 0.6)
			    }
			})
}
#add tracks just for genic regions 
add_gene_track <- function(df, col, track_height = 0.02) {
    circos.genomicTrack(df[, c("chr","start","end")],
			ylim = c(0, 1),
			track.height = track_height,
			bg.border = "black",
			bg.lwd = 0.3,
			panel.fun = function(region, value, ...) {
			    circos.genomicRect(region,
					       ybottom = 0, ytop = 1,
					       col = col, border =NA,lwd=0.0)
			})
}


cs <- read.table(fai_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
cs <- cs[, 1:2]
colnames(cs) <- c("chr","length")

athbed    <- read_bed3(bed_ath)
mslybed   <- read_bed3(bed_msly)
hvubed    <- read_bed3(bed_hvu)
systembed <- read_bed3(bed_system)
genebed <- read_bed3(bed_genes)
repeatbed <- read_bed3(bed_repeats)

keep_chr <- cs$chr
athbed    <- athbed[athbed$chr %in% keep_chr, ]
mslybed   <- mslybed[mslybed$chr %in% keep_chr, ]
hvubed    <- hvubed[hvubed$chr %in% keep_chr, ]
systembed <- systembed[systembed$chr %in% keep_chr, ]
genebed <- genebed[genebed$chr %in% keep_chr, ]
repeatbed <- repeatbed[repeatbed$chr %in% keep_chr, ]

athbed$chr    <- factor(athbed$chr,    levels = cs$chr)
mslybed$chr   <- factor(mslybed$chr,   levels = cs$chr)
hvubed$chr    <- factor(hvubed$chr,    levels = cs$chr)
systembed$chr <- factor(systembed$chr, levels = cs$chr)
genebed$chr <- factor(genebed$chr, levels = cs$chr)
repeatbed$chr <- factor(repeatbed$chr, levels = cs$chr)
pdf(out_pdf, width = 9, height = 9)
circos.clear()
circos.par(
	   cell.padding = c(0, 0, 0, 0),
	   track.margin = c(0.003, 0.003),
	   start.degree = 90,
	   gap.after = rep(0.75, nrow(cs))
)

# Initialize
circos.initialize(factors = cs$chr, xlim = cbind(0, cs$length))

# Outer label/axis track: show chromosome names + ticks
circos.track(ylim = c(0, 1), bg.col=NA,bg.border = NA, track.height = 0.05,
	     panel.fun = function(x, y) {
		 sector <- CELL_META$sector.index
		 circos.axis(h = "bottom", labels.cex = 0.35, major.tick.length = 0.06,col = "black", labels.col = "black")
		 circos.text(CELL_META$xcenter, CELL_META$ycenter,
			     labels = sector,
			     facing = "bending.outside", niceFacing = TRUE,
			     cex = 0.5, adj = c(0.5, 2.0),col="black")
	     })

add_gene_track(genebed,"blue")
add_gene_track(repeatbed,"#fc9803")
add_tick_track(athbed,cols["ath"])
add_tick_track(mslybed,cols["msly"])
add_tick_track(hvubed,cols["hvu"])
add_tick_track(systembed, cols["system"])

# Legend
legend("topleft",
       legend = c("genes","repeatregion","Ath", "mSly", "Hvu", "System"),
       col = c("blue","#fc9803",unname(cols[c("ath","msly","hvu","system")])),
       lty = 1, lwd = 3, bty = "n")

circos.clear()
dev.off()

#zooming in to short scaffolds
scaffolds_zoom <- cs$chr[grepl("^scaffold_(8|9|10|11).*pilon$", cs$chr, ignore.case = TRUE)]
if (length(scaffolds_zoom) == 0) {
    scaffolds_zoom <- intersect(cs$chr, c("scaffold_8_pilon","scaffold_9_pilon","scaffold_10_pilon","scaffold_11_pilon"))
}
if (length(scaffolds_zoom) == 0) stop("No zoom scaffolds matched.")

cs_zoom <- cs[cs$chr %in% scaffolds_zoom, , drop = FALSE]

genebed_zoom <- genebed[genebed$chr %in% scaffolds_zoom, ]
repeatbed_zoom <- repeatbed[repeatbed$chr %in% scaffolds_zoom, ]
athbed_zoom <- athbed[athbed$chr %in% scaffolds_zoom, ]
mslybed_zoom <- mslybed[mslybed$chr %in% scaffolds_zoom, ]
hvubed_zoom <- hvubed[hvubed$chr %in% scaffolds_zoom, ]
systembed_zoom <- systembed[systembed$chr %in% scaffolds_zoom, ]

#matching orders
genebed_zoom$chr <- factor(genebed_zoom$chr,    levels = cs_zoom$chr)
repeatbed_zoom$chr <- factor(repeatbed_zoom$chr,    levels = cs_zoom$chr)
athbed_zoom$chr <- factor(athbed_zoom$chr,    levels = cs_zoom$chr)
mslybed_zoom$chr <- factor(mslybed_zoom$chr,   levels = cs_zoom$chr)
hvubed_zoom$chr <- factor(hvubed_zoom$chr,    levels = cs_zoom$chr)
systembed_zoom$chr <- factor(systembed_zoom$chr, levels = cs_zoom$chr)
out_pdf_zoom <- "circos_variants_scaffold_8-11.pdf"
pdf(out_pdf_zoom, width = 9, height = 9)
circos.clear()
circos.par(
	   cell.padding = c(0, 0, 0, 0),
	   track.margin = c(0.003, 0.003),
	   start.degree = 120,
	   gap.after = c(rep(1, nrow(cs_zoom) - 1), 360 - 30)
)

circos.initialize(factors = cs_zoom$chr, xlim = cbind(0, cs_zoom$length))

circos.track(ylim = c(0, 1), bg.col = NA, bg.border = NA, bg.lwd = 0.5, track.height = 0.05,
	     panel.fun = function(x, y) {
		 sector <- CELL_META$sector.index
		 xlim <- CELL_META$xlim
		 brks_mb <- pretty(c(xlim[1], xlim[2]) / 1e6)
		 major_at <- brks_mb * 1e6
		 major_at <- major_at[major_at >= xlim[1] & major_at <= xlim[2]]
		 major_at <- unique(c(xlim[1], major_at, xlim[2]))
		 keep <- major_at >= xlim[1] & major_at <= xlim[2]
		 labs <- paste0(round(major_at/1e6, 2), " Mb")
		 circos.axis(h = "bottom",
			     major.at = major_at,
			     labels = labs,
			     labels.cex = 0.35,
			     labels.col = "black",
			     col = "black",
			     labels.niceFacing = TRUE,
			     direction = "outside",
			     major.tick.length = 0.06)
		 circos.text(CELL_META$xcenter, CELL_META$ycenter,
			     labels = sector,
			     facing = "bending.inside", niceFacing = TRUE,
			     cex = 0.5, adj = c(0.5, -2.0), col = "black")
	     })

add_gene_track(genebed_zoom,"blue")
add_gene_track(repeatbed_zoom,"#fc9803")
add_tick_track(athbed_zoom,cols["ath"],label = "Ath",label.col = "white")
add_tick_track(mslybed_zoom,cols["msly"])
add_tick_track(hvubed_zoom,cols["hvu"])
add_tick_track(systembed_zoom, cols["system"])

circos.clear()
dev.off()
