# Plot smoothed curves for published genes of interest
# Edited 3/14/22

library(SingleCellExperiment)
library(tidyverse)
library(Seurat)
library(princurve)
library(TrajectoryUtils)
library(slingshot)
library(tradeSeq)
library(RColorBrewer)

setwd("/proj/jmsimon/Song/SPLITseq/Analysis_052021/")

load("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_slingshot_pseudotime_lineages_tradeSeq_DEconditions.Rdata")

set.seed(123)
sce = as.SingleCellExperiment(subset,assay = "RNA")
sce.cts = as.SingleCellExperiment(subset,assay = "RNA")

goi = c("Gfap","Apoe","Sox11","Eomes","Aldoc","Stmn1","Hopx","Fgfr3","Ccnd2","Clu","Sox9","Hmgn2","Neurod1","Nnat","Dcx")
sce <- fitGAM(counts = counts(sce), sds = exp.sling, conditions = factor(colData(sce)$trt), genes = goi)


# Create color vector for color-coding output plots
exp.col.pre = rep("black",length(exp.ids.pre))
names(exp.col.pre) = names(exp.ids.pre)

#Set cluster color vectors
col.15.yfp <- "#007800"			# Dark green
col.15.chr2 <- "#00DC00"		# Light green
col.16.yfp <- "#000078"		# Dark blue
col.16.chr2 <- "#0000DC"		# Light blue
col.25.yfp <- "#780000"		    # Dark red
col.25.chr2 <- "#dc0000"		# Light red
col.2.yfp <- "dark gray"		    # Dark red
col.2.chr2 <- "light gray"		# Light red

exp.col.pre[intersect(grep("RGL/Astrocyte",exp.ids.pre,perl=TRUE),grep("yfp",exp.gt.pre,perl=TRUE))] <- col.15.yfp
exp.col.pre[intersect(grep("RGL/Astrocyte",exp.ids.pre,perl=TRUE),grep("chr2",exp.gt.pre,perl=TRUE))] <- col.15.chr2
exp.col.pre[intersect(grep("Neuroblast",exp.ids.pre,perl=TRUE),grep("yfp",exp.gt.pre,perl=TRUE))] <- col.16.yfp
exp.col.pre[intersect(grep("Neuroblast",exp.ids.pre,perl=TRUE),grep("chr2",exp.gt.pre,perl=TRUE))] <- col.16.chr2
exp.col.pre[intersect(grep("RGL Like",exp.ids.pre,perl=TRUE),grep("yfp",exp.gt.pre,perl=TRUE))] <- col.25.yfp
exp.col.pre[intersect(grep("RGL Like",exp.ids.pre,perl=TRUE),grep("chr2",exp.gt.pre,perl=TRUE))] <- col.25.chr2
exp.col.pre[intersect(grep("Early GC",exp.ids.pre,perl=TRUE),grep("yfp",exp.gt.pre,perl=TRUE))] <- col.2.yfp
exp.col.pre[intersect(grep("Early GC",exp.ids.pre,perl=TRUE),grep("chr2",exp.gt.pre,perl=TRUE))] <- col.2.chr2


pdf("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_slingshot_pseudotime_lineages_tradeSeq_conditions_smoothers_goi_v2.pdf")
for(i in 1:length(goi)) {
	show(plotSmoothers(sce, assays(sce)$counts, gene = goi[i], border = FALSE, pointCol = paste(sce.cts$ident,sce.cts$trt,sep="_")) +
	ggtitle(goi[i]) +
	scale_color_manual(breaks=c("15_yfp","15_chr2","16_yfp","16_chr2","25_yfp","25_chr2","2_yfp","2_chr2"),values=c("#007800","#00DC00","#000078","#0000DC","#780000","#dc0000","dark gray","light gray")))
}
dev.off()
