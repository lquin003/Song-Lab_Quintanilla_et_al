library(rlang,lib="/proj/jmsimon/Rlibs36")
library(vctrs,lib="/proj/jmsimon/Rlibs36")
library(backports,lib="/proj/jmsimon/Rlibs36")
library(cli,lib="/proj/jmsimon/Rlibs36")
library(data.table,lib="/proj/jmsimon/Rlibs36")
library(presto,lib="/proj/jmsimon/Rlibs36")
library(Seurat,lib="/proj/jmsimon/Rlibs36")
library(rstudioapi,lib="/proj/jmsimon/Rlibs36")
library(withr,lib="/proj/jmsimon/Rlibs36")
library(tidyverse,lib="/proj/jmsimon/Rlibs36")
library(ggplot2,lib="/proj/jmsimon/Rlibs36")
library(sctransform,lib="/proj/jmsimon/Rlibs36")
library(farver,lib="/proj/jmsimon/Rlibs36")
library(usethis,lib="/proj/jmsimon/Rlibs36")
library(labeling,lib="/proj/jmsimon/Rlibs36")
library(Matrix,lib="/proj/jmsimon/Rlibs36")
library(stringdist,lib="/proj/jmsimon/Rlibs36")
library(circlize,lib="/proj/jmsimon/Rlibs36")
library(ComplexHeatmap,lib="/proj/jmsimon/Rlibs36")
library(princurve,lib="/proj/jmsimon/Rlibs36")
library(slingshot,lib="/proj/jmsimon/Rlibs36")
library(umap,lib="/proj/jmsimon/Rlibs36")

setwd("/proj/jmsimon/Song/SPLITseq/Analysis_052021/")

dg.integrated = readRDS("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered.rds")

subset = subset(dg.integrated,idents = c(15,16,25,2))

# Get most informative genes of these clusters
vargenes<- presto::wilcoxauc(subset, 'seurat_clusters', seurat_assay = 'SCT')
top_vargenes = top_markers(vargenes, n = 50, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)
all_markers<- top_vargenes %>%
	select(-rank) %>% 
	unclass() %>% 
	stack() %>%
	pull(values) %>%
	unique() %>%
	.[!is.na(.)]


# Retrieve integrated data for subsetted object
exp.pre = GetAssayData(subset,assay="integrated")

# Filter for just the marker genes
exp.pre.filter = exp.pre[rownames(exp.pre) %in% all_markers,]

# Run UMAP
exp.pc.pre = umap::umap(t(as.matrix(exp.pre.filter)))
exp.pc = exp.pc.pre$layout

idents = c(15,16,25,2)
names(idents) = c("RGL/Astrocyte","Neuroblast","RGL Like","Early GC")

exp.ids.pre = as.factor(names(idents[match(as.numeric(as.character(Idents(subset)[rownames(exp.pc)])),idents)]))
names(exp.ids.pre) = names(Idents(subset)[rownames(exp.pc)])

exp.gt.pre = subset$trt[rownames(exp.pc)]

# Run slingshot
exp.sling = slingshot(exp.pc, exp.ids.pre,start.clus='RGL/Astrocyte',end.clus='Early GC')
print(exp.sling)

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


# Plot linear and curve trajectories superimposed on PCA plot of values
pdf("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_slingshot_pseudotime_lineage_all_plot.pdf")
plot(exp.pc, col = exp.col.pre, pch=19, cex=0.5, asp=1)
plot(exp.sling, type="curves", add=T, show.constraints = T)
dev.off()
