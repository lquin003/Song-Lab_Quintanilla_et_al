library(data.table,lib="/proj/jmsimon/Rlibs36")
library(presto,lib="/proj/jmsimon/Rlibs36")
library(Rcpp,lib="/proj/jmsimon/Rlibs36")
library(rlang,lib="/proj/jmsimon/Rlibs36")
library(vctrs,lib="/proj/jmsimon/Rlibs36")
library(backports,lib="/proj/jmsimon/Rlibs36")
library(cli,lib="/proj/jmsimon/Rlibs36")
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

#Set working directory
setwd("/proj/jmsimon/Song/SPLITseq/Analysis_052021/")

# Read in zUMIs data
s1 = readRDS(file="/proj/jmsimon/Song/SPLITseq/Song_DG_SPLITseq.dgecounts.rds")
s2 = readRDS(file="/proj/jmsimon/Song/SPLITseq/Song_DG_SPLITseq_SL2.dgecounts.rds")

# Rename cells to include sample number, and divide barcodes
s1.cnames = gsub(colnames(s1$umicount$inex$all),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="Song_DG_SL1_zUMIs_\\1-\\2-\\3",perl=T)
s2.cnames = gsub(colnames(s2$umicount$inex$all),pattern="([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})",replacement="Song_DG_SL2_zUMIs_\\1-\\2-\\3",perl=T)

colnames(s1$umicount$inex$all) = s1.cnames
colnames(s2$umicount$inex$all) = s2.cnames

# Collapse transcripts to genes
tx2gene = read.table("/proj/jmsimon/genomeAnnotation/gencode.vM26.annotation.lookup.txt",header=F,sep="\t",col.names=c("tx","gene"))

s1.txId = rownames(s1$umicount$inex$all)
s1.geneId <- as.vector(tx2gene$gene[match(s1.txId, tx2gene$tx)])
s1.tx.grp = t(sparse.model.matrix(~ 0 + s1.geneId))
s1.summarized <- s1.tx.grp %*% s1$umicount$inex$all  
rownames(s1.summarized) <- rownames(s1.summarized) %>% str_replace_all(".+.geneId","")

s2.txId = rownames(s2$umicount$inex$all)
s2.geneId <- as.vector(tx2gene$gene[match(s2.txId, tx2gene$tx)])
s2.tx.grp = t(sparse.model.matrix(~ 0 + s2.geneId))
s2.summarized <- s2.tx.grp %*% s2$umicount$inex$all  
rownames(s2.summarized) <- rownames(s2.summarized) %>% str_replace_all(".+.geneId","")


# Import data into Seurat objects
s1.seurat = CreateSeuratObject(s1.summarized,project="Song_DG_SL1")
s2.seurat = CreateSeuratObject(s2.summarized,project="Song_DG_SL2")

# Check raw dimensions
dim(s1.seurat)
[1] 37008 28707
dim(s2.seurat)
[1] 37046 26818


# Combine into one seurat object for now
combined = merge(x=s1.seurat, y=s2.seurat)


# Look up barcodes (last 8, this orientation) to get sample name
# Searches for a hamming distance within 1, but for this dataset, it seems to be all exact matches
# If no match, cell name will be prepended with "TRASH" for subsequent filtering
# Matches will have cell name prepended with GTrep-ampBC

bclookup=read.table("bcLookup.txt",header=F,sep="\t")

ids = colnames(combined@assays$RNA)
ids_l8 = substr(ids,37,44)


newNames = ids
for(i in 1:length(ids_l8)) {
	ham = stringdist(ids_l8[i],bclookup$V1,method="hamming")
	ct = length(ham[ham<2])
	if(ct==0){
		newNames[i] = paste0("TRASH_",newNames[i])
	}
	else if(ct > 1){
		print(paste0("Barcode ",ids_l8[i]," had ",ct," matches within a Hamming distance of 1"))
		newNames[i] = paste0("TRASH_",newNames[i])
	}
	else{
		bestMatch = which.min(ham)
		matchingBC = as.character(bclookup[bestMatch,1])
		newNames[i] = paste0(as.character(bclookup[bestMatch,2]),"_",newNames[i])
	}
}
combined = RenameCells(combined,new.names = newNames)


# Parse out sample metadata, drop TRASH cells
cnames = colnames(combined)
names = as.character(gsub("(.+)_(Song_DG_SL[[:digit:]])_(zUMIs)_(.+)","\\1",cnames,perl=T))
names = as.character(gsub("TRASH","TRASH-TRASH",names,perl=T))
names = as.character(gsub("yfp_","yfp-",names,perl=T))
names = as.character(gsub("chr2_","chr2-",names,perl=T))

trt = as.character(gsub("(.+)-(.+)","\\1",names,perl=T))
reps = as.character(gsub("(.+)-(.+)","\\2",names,perl=T))

combined = AddMetaData(combined,metadata=reps,col.name="rep")
combined = AddMetaData(combined,metadata=trt,col.name="trt")
combined = AddMetaData(combined,metadata=names,col.name="trtrep")

combined = subset(combined, subset = trt!="TRASH") 

> dim(combined)
[1] 39469 51326


as.data.frame(table(names))

         names Freq
1    chr2-Rep1 4653
2    chr2-Rep2 4225
3    chr2-Rep3 4308
4    chr2-Rep4 4905
5    chr2-Rep5 9536
6  TRASH-TRASH 4199
7     yfp-Rep1 6950
8     yfp-Rep2 4130
9     yfp-Rep3 3683
10    yfp-Rep4 5033
11    yfp-Rep5 3903

combined <- PercentageFeatureSet(object = combined, pattern = "^mt-", col.name = "percent.mt")

combined.filtered <- subset(combined, subset = nCount_RNA > 1250 & nFeature_RNA > 750 & percent.mt < 10)

> dim(combined.filtered)
[1] 39469 21647

as.data.frame(table(combined.filtered$trtrep))
        Var1 Freq
1  chr2-Rep1 3019
2  chr2-Rep2 1670
3  chr2-Rep3 2485
4  chr2-Rep4 2320
5  chr2-Rep5 2081
6   yfp-Rep1 2555
7   yfp-Rep2 1325
8   yfp-Rep3 1469
9   yfp-Rep4 2902
10  yfp-Rep5 1821

all.list <- SplitObject(combined.filtered, split.by = "trtrep")

for (i in 1:length(x = all.list)) {
	all.list[[i]] <- PercentageFeatureSet(object = all.list[[i]], pattern = "^mt-", col.name = "percent.mt")
	all.list[[i]] <- SCTransform(all.list[[i]], vars.to.regress = c("percent.mt"), verbose = FALSE,return.only.var.genes=FALSE)
}


#Save progress
save.image("Song_DG_SPLITseq_zUMIs_combined_scTransform_052021.Rdata")










# Ran separate job with 400GB RAM, it used X
library(rlang,lib="/proj/jmsimon/Rlibs36")
library(vctrs,lib="/proj/jmsimon/Rlibs36")
library(backports,lib="/proj/jmsimon/Rlibs36")
library(cli,lib="/proj/jmsimon/Rlibs36")
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

#Set working directory
setwd("/proj/jmsimon/Song/SPLITseq/Analysis_052021/")


load("Song_DG_SPLITseq_zUMIs_combined_scTransform_052021.Rdata")
options(future.globals.maxSize=5242880000)		# This was needed otherwise I got an error about something being too large/above limits. This particular value is 5000*1024^2, based off of this page: https://stackoverflow.com/questions/40536067/how-to-adjust-future-global-maxsize-in-r
dg.features <- SelectIntegrationFeatures(object.list = all.list, nfeatures = 5000)
all.list <- PrepSCTIntegration(object.list = all.list, anchor.features = dg.features, verbose = TRUE)
dg.anchors <- FindIntegrationAnchors(object.list = all.list, normalization.method = "SCT", anchor.features = dg.features, verbose = TRUE)
dg.integrated <- IntegrateData(anchorset = dg.anchors, normalization.method = "SCT", verbose = TRUE)

#Save image at point of object integration
saveRDS(dg.integrated,"Song_DG_SPLITseq_zUMIs_combined_scTransform_052021_integrated.rds")



# Load back into interactive session
dg.integrated = readRDS("Song_DG_SPLITseq_zUMIs_combined_scTransform_052021_integrated.rds")





#Run PCA
dg.integrated <- RunPCA(dg.integrated, verbose = FALSE, npcs = 100)
dg.integrated <- RunUMAP(dg.integrated, dims = 1:100)
dg.integrated <- FindNeighbors(dg.integrated, dims = 1:100, verbose = FALSE)
dg.integrated <- FindClusters(dg.integrated, verbose = FALSE, resolution = 1.5, algorithm=2)


saveRDS(dg.integrated,"Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_res1.5.rds")
#dg.integrated = readRDS("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered.rds")



#Clustering plots
pdf("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_res1.5.pdf",width=7)
DimPlot(dg.integrated,reduction = "umap", label = TRUE)
DimPlot(dg.integrated,reduction = "umap", group.by = "trt")
DimPlot(dg.integrated,reduction = "umap", group.by = "trtrep", cols = c("#285F2F","#3B9244","#71C163","#A1D4A1","#D3DFCC","#9F384A","#AE5C61","#BE7D7D","#CFA19F","#E3CAC8"))
dev.off()

#Comparison plots
pdf("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_comparison_res1.5.pdf",width=40)
DimPlot(dg.integrated,reduction = "umap", label = TRUE, split.by = "trtrep")
dev.off()


# QC plots
#Make violin plots
features = c("nCount_RNA","nFeature_RNA","percent.mt")
plotwidth = 100
plotheight = 15
pdf("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_vlnplot_QC_res1.5.pdf",width=plotwidth,height=plotheight)
for (i in features) {
    print(VlnPlot(dg.integrated, features = i, split.by = "trtrep"))
}
dev.off()


#Make violin plots split just by treatment alone
features = c("nCount_RNA","nFeature_RNA","percent.mt")
plotwidth = 100
plotheight = 15
pdf("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_vlnplot_QC_splitByTreatment_res1.5.pdf",width=plotwidth,height=plotheight)
for (i in features) {
    print(VlnPlot(dg.integrated, features = i, split.by = "trt"))
}
dev.off()


# Compute cell proportions by cluster by replicate

trt <- unique(dg.integrated@meta.data$trt)
rep <- unique(dg.integrated@meta.data$rep)

comp.prop <- vector()
names.prop <- vector()
for (i in 1:length(trt)) {
        #Subset dg.integrated object
        sub <- dg.integrated[,dg.integrated$trt == trt[i]]
        
        #Generate proportion table by replicate
        sub.table <- prop.table(table(Idents(sub), sub$rep), margin = 2)
        
        #Create table to add to comp (code necessary when not all clusters are in subset)
        sub.final <- array(0, dim=c(length(levels(dg.integrated)),length(unique(sub$rep))))
        rownames(sub.final) <- levels(dg.integrated)
        colnames(sub.final) <- c(1:length(unique(sub$rep)))
        sub.final[rownames(sub.final) %in% rownames(sub.table),] <- sub.table
        
        #Add to compilation and colnames vector
        comp.prop <- cbind(comp.prop,sub.final)
        names.prop <- c(names.prop, paste0(trt[i],"_",colnames(sub.table)))
}
colnames(comp.prop) <- names.prop

write.table(comp.prop,"Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_proportions_replicates_res1.5.txt",col.names=NA,row.names=TRUE,quote=FALSE,sep="\t")


# import using tidy's readr function
tbl = as_tibble(rownames_to_column(as.data.frame(comp.prop),var = "Cluster"))

# turn data matrix into tidy tibble
tibble = tbl %>% gather(key="Sample",value="Proportion",-Cluster) %>% separate(Sample,sep="_",into=c("Sample","Replicate"),convert=T)

# boxplot
pdf("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_proportions_replicates_boxplots_res1.5.pdf",width=15,height=10)
tibble %>%
    mutate(Cluster = fct_relevel(Cluster, unique(tibble$Cluster))) %>%
    mutate(Sample = fct_relevel(Sample, "yfp", "chr2")) %>%
    ggplot(aes(fill = Sample, x= Sample, y=Proportion)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(pch = 21, position = position_jitterdodge(jitter.width=0.9)) +
    facet_wrap(~Cluster, scales="free") +
    xlab("")
dev.off()







# Broke this into separate script, took X hrs to run with 100GB RAM
library(rlang,lib="/proj/jmsimon/Rlibs36")
library(vctrs,lib="/proj/jmsimon/Rlibs36")
library(backports,lib="/proj/jmsimon/Rlibs36")
library(cli,lib="/proj/jmsimon/Rlibs36")
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

#Set working directory
setwd("/proj/jmsimon/Song/SPLITseq/Analysis_052021/")

dg.integrated = readRDS("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_res1.5.rds")

markers <- FindAllMarkers(dg.integrated, assay="SCT", slot="scale.data", only.pos=T, logfc.threshold = 0.25)
write.table(markers,"Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_markers_SCT_res1.5.txt",quote=F,sep="\t",col.names=NA)





goi = read_tsv("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_markers_SCT_res1.5.txt") %>% 
    group_by(cluster) %>%
    arrange(p_val_adj, .by_group = T) %>%
    slice_head(n = 10) %>%
    pull(gene)

goi.uniq = unique(goi)

goi.uniq2 = unique(c("Hopx","Cd38","Fam227b","Kcng4","Pamr1","Dnaic1","Grin2c","Ppp1r3g","AI464131","Mertk","Gm266","Cox4i2","Rgs5","Rhcg","Vnn1","Lpar1","Hopxos","1500015O10Rik","Wnt8b","Cacng5","E330013P04Rik","Cdk1","Cdkn3","Aurkb","Cdca3","Pbk","Plk1","Top2a","Ccnb1","Ccnb2","Spc25","Hs3st1","Tac2","Cldn3","Mfng","Mpped1","Neurod4","Slc17a6","Elavl4","Tbr1","Eomes","Rhov","Kcnf1","Gabra5","Npas4","Ctxn2","Mtus2","Cntnap5a","Shisa7","Plk5","Caly","Insc","Tac2","Syt2","Gal","Pdlim4","Calb2","Rbp1","Ddah2","Plin2","Draxin","Stc1","Sdk2","Rasd1","Fxyd7","Il16","Dsp","Rspo3","Prickle2","A330050F15Rik","Smim3","Rgs4","Trpc6","BC030499","Ntng1","Rasgrp1","Nos1","Snhg11","Smad3","Adarb2","Ipcef1","Cst3","Fos","Jun","Trem2","Cx3cr1",goi.uniq))


pdf("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_dotplot_res1.5.pdf",width=50,height=10)
DotPlot(dg.integrated,features=rev(goi.uniq2),assay="SCT",cols = c("blue","red")) + RotatedAxis()
dev.off()



# Now extract just the GC clusters and call markers individually to get a better sense of what subtype each are

dg.cells = c(0,1,2,3,4,5,7,8,10,11,14,22,26,30,32)
dg.subset = subset(dg.integrated, subset = seurat_clusters %in% dg.cells)
dg.markers <- FindAllMarkers(dg.subset, assay="SCT", slot="scale.data", only.pos=T, logfc.threshold = 0.25)

write.table(dg.markers,"Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_markers_DGcellsOnly_SCT.txt",quote=F,sep="\t",col.names=NA)


# Write out average expression tables for all clusters
dg.avg = AverageExpression(dg.integrated)
write.table(dg.avg$SCT,"Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_averages_SCT_res1.5.txt",quote=F,sep="\t",col.names=NA)

// 
// # Save progress for today
// 
// saveRDS(dg.integrated,"Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_withMarkers.rds")
// save.image("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_withMarkers.Rdata")
// 
// 
// # Make feature plot of marker genes
// DefaultAssay(dg.integrated) = "SCT"
// #pdf("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_featureplot_MKI67_DCX_SOX2.pdf",width=18)
// FeaturePlot(object = dg.integrated, features = c("Mki67","Dcx","Rgs4","Stc1","Hopx","Trpc6"),ncol=3)
// #dev.off()



luismarkers = unique(c("Aldoc","Apoe","Aqp4","Atp1a2","Clu","Csf1","Gfap","Gja1","S100b","Slc1a2","Cst3","Hopx","Hes5","Nestin","Gfap","S100b","Sox2","Vim","Sox6","Gli1","Ascl1","Vim","Neurog2","Sox2","Cdk1","Cdkn3","Aurkb","Cdca3","Pbk","Plk1","Top2a","Ccnb1","Ccnb2","Spc25","Eomes","Sox5","Neurod1","Calb2","Eomes","Neurod6","Tbr1","Mex3b","Cdk5r1","Sema3c","Sema5a","Cntnap5a","Kcnf1","Neurod2","Dcx","Dsp","Slc4a4","Camk4","Prox1","hectd2","stxbp6","camk2","Tenm1","Pdzd2","Stc1","Npy1r","Neurod2","Snhg11","Smad3","Sv2c","Zfp365","Kcnip3","Nbl1","Krt2","Map2k1","Rin1","Trpc6","Trpm3","camk2a","Prox1","Neurod2","Mbp","Mog","Plp1","Mal","Car2","Cryab","Mobp","Olig1","SV2b","Cntn6","CCK","Calb2","Grm8","Chga","Csf2rb2","Drd2","Gal","Glp1r","Prrx1","Amigo2","Nrxn1","Cplx2","Ncdn","Olfm1","Synpr","Sema5a","Cpne4","Ly6e","Npy2r","Atp2b1","Hpca","Spink8","Arpc2","Ptn","Cck","Tspan13","Plk2","Itm2b","Pex5l","Neurod6","Ctss","Cst3","Ctsd","P2ry12","C1qa","Csf1r","C1qc","Cx3cr1","C1qb","Hexb","Zfp36","Ptprc","Dcn","Nnat","Timp2","Cpne7","Ndufc2","Nptxr","Cpe","Rspo3","Pde1a","Grp","Bace2","Fibcd1","Camk2d","Cnr1","Grik1","Gad1","Gad2","Slc32a1","Dlx6","Dlx5","Nrxn3","Grik1","Erbb4","Gm14204","Qpct","Dlx6os1","Six3","Flt1","Rgs5","Sptbn1","Ptprb","Rn18s","Abcb1a","Sgms1","Utrn","Slco1a4","Adgrl4","Fn1","Calml4","Ak7","Ccdc153","1110017D15Rik","Dynlrb2","1500015O10Rik","Tmem212","Fam216b","Ttr","Hdc","Lbp","Spag16","Acta2","Tpm1","Tpm2","Crip1","Tagln","Apod","Vtn","Atp1a2","Malat1","Itih5","Dcn"))
luis.markers = unique(luismarkers[luismarkers %in% rownames(dg.integrated)])

pdf("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_dotplot_res1.5_luisMarkers.pdf",width=30,height=10)
DotPlot(dg.integrated,features=rev(luis.markers),assay="SCT",cols = c("blue","red")) + RotatedAxis()
dev.off()



dup.dg.integrated = dg.integrated

clustOrder = rev(c(14, 12, 28, 27, 23, 32, 26, 4, 3, 2, 5, 0, 6, 1, 8, 10, 20, 30, 7, 34, 35, 16, 17, 18, 21, 22, 19, 15, 13, 36, 24, 31, 9, 11, 29, 33, 25, 37))
Idents(dup.dg.integrated) <- factor(Idents(dup.dg.integrated), levels = clustOrder)

pdf("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_dotplot_res1.5_luisMarkers_ordered.pdf",width=30,height=10)
DotPlot(dup.dg.integrated,features=rev(luis.markers),assay="SCT",cols = c("blue","red")) + RotatedAxis()
dev.off()





# Trying to resolve whether clusters 14 and 28 in the res1.5 data are heterogeneous mixes containing astrocytes and RGL cells

astroRGLgenes = c("Abhd3","Id3","Ncan1","Acsl6","Sepp1","Pla2g7","Slc4a4","Bcan","Htra4","Aqp4","Scg3","Gpr37l1","Cspg5","Atp1b2","Gja1","Glul","Clu","Sparcl1","Atp1a2","Aldoc","Hopx","Pbx1","Xist","Btg1","Thrsp","Notch2","Gm10561","Itih3","Rpl23a-ps3","Jun","Tfap2c","Emx1","Hopxos","Btg2","Frzb","2700094K13Rik","E330013P04Rik","Ascl1","Sox4","Riiad1")

dg.integrated = readRDS("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_res1.5.rds")
dg.sub = subset(dg.integrated, subset = seurat_clusters %in% c(14,28))
DoHeatmap(dg.sub,features = astroRGLgenes, assay = "SCT", slot = "scale.data")


# Use *extremely* fast method for calling markers, in the package "presto"

vargenes<- presto::wilcoxauc(dg.integrated, 'seurat_clusters', seurat_assay = 'SCT')
top_vargenes = top_markers(vargenes, n = 50, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)

all_markers<- top_vargenes %>%
	select(-rank) %>% 
	unclass() %>% 
	stack() %>%
	pull(values) %>%
	unique() %>%
	.[!is.na(.)]


# Create clustered heatmap using ComplexHeatmap

mat<- dg.sub@assays$SCT[rownames(dg.sub@assays$SCT) %in% c(astroRGLgenes,all_markers), ] %>% as.matrix()
## scale the rows
mat<- t(scale(t(mat)))

cluster_anno<- dg.sub@meta.data$seurat_clusters
quantile(mat, c(0.1, 0.95))
col_fun = circlize::colorRamp2(c(-1, 0, 3), c("white", "white", "red"))

Heatmap(mat, name = "Expression",  
        column_split = factor(cluster_anno),
        cluster_columns = TRUE,
        show_column_dend = TRUE,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 8),
        column_gap = unit(0.5, "mm"),
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 7),
        column_title_rot = 90,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
        show_column_names = FALSE)





# Iterate clustering resolution to evaluate astro/RGL clusters

dg.integrated = readRDS("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered.rds")

pdf("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_iterateRes.pdf")
dg.integrated <- FindClusters(dg.integrated, verbose = FALSE, resolution = 1.5, algorithm=2)
DimPlot(dg.integrated,reduction = "umap", label = TRUE)
dg.integrated <- FindClusters(dg.integrated, verbose = FALSE, resolution = 1.75, algorithm=2)
DimPlot(dg.integrated,reduction = "umap", label = TRUE)
dg.integrated <- FindClusters(dg.integrated, verbose = FALSE, resolution = 2, algorithm=2)
DimPlot(dg.integrated,reduction = "umap", label = TRUE)
dg.integrated <- FindClusters(dg.integrated, verbose = FALSE, resolution = 2.25, algorithm=2)
DimPlot(dg.integrated,reduction = "umap", label = TRUE)
dg.integrated <- FindClusters(dg.integrated, verbose = FALSE, resolution = 2.5, algorithm=2)
DimPlot(dg.integrated,reduction = "umap", label = TRUE)
dg.integrated <- FindClusters(dg.integrated, verbose = FALSE, resolution = 2.75, algorithm=2)
DimPlot(dg.integrated,reduction = "umap", label = TRUE)
dg.integrated <- FindClusters(dg.integrated, verbose = FALSE, resolution = 3, algorithm=2)
DimPlot(dg.integrated,reduction = "umap", label = TRUE)
dev.off()



# Focus on res=2.0
dg.integrated <- FindClusters(dg.integrated, verbose = FALSE, resolution = 2, algorithm=2)






# Remade on 12/21/21 to capture PDF
dg.integrated = readRDS("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered.rds")

vargenes<- presto::wilcoxauc(dg.integrated, 'seurat_clusters', seurat_assay = 'SCT')
top_vargenes = top_markers(vargenes, n = 50, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)
all_markers<- top_vargenes %>%
	select(-rank) %>% 
	unclass() %>% 
	stack() %>%
	pull(values) %>%
	unique() %>%
	.[!is.na(.)]

#dg.integrated <- BuildClusterTree(dg.integrated, reorder = F, verbose = T, assay="SCT", features = all_markers)
#PlotClusterTree(dg.integrated)

dg.avg = AverageExpression(dg.integrated)
dg.sct.cor = cor(dg.avg$SCT[all_markers,])

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
col_fun(seq(-1, 1))

pdf("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_prestoMarkers_SCT_pearson_cor_heatmap.pdf",height=8,width=10)
ComplexHeatmap::Heatmap(dg.sct.cor,name="Pearson correlation",column_names_side = "top",row_names_side = "left",col = col_fun, clustering_distance_columns=function(x) as.dist(1-cor(t(x))), clustering_distance_rows=function(x) as.dist(1-cor(t(x))), clustering_method_columns = "average", clustering_method_rows = "average")
dev.off()








# DE testing
library(rlang,lib="/proj/jmsimon/Rlibs36")
library(vctrs,lib="/proj/jmsimon/Rlibs36")
library(backports,lib="/proj/jmsimon/Rlibs36")
library(cli,lib="/proj/jmsimon/Rlibs36")
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

#Set working directory
setwd("/proj/jmsimon/Song/SPLITseq/Analysis_052021/")

dg.integrated = readRDS("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered.rds")

for (i in as.numeric(levels(dg.integrated$seurat_clusters))){	
	# Instantiate all variables to make sure no carry-overs from previous clusters
	names = 0
	nullnames = 0
	ctl1names = 0
	diff = 0
	sig.up = 0
	sig.dn = 0
	sig.up.gp = 0
	sig.dn.gp = 0
	sig.up.gp.flat = 0
	sig.dn.gp.flat = 0

	print(paste0("Working on cluster ",i))

	names = colnames(dg.integrated)[dg.integrated$seurat_clusters==i]
	nullnames = names[grepl("chr2",names)]
	ctl1names = names[grepl("yfp",names)]

	diff = FindMarkers(dg.integrated, ident.1 = nullnames, ident.2=ctl1names, assay="SCT", slot="scale.data")

	write.table(diff,paste0("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_chr2_vs_yfp_cluster",i,"_DEresults_061821.txt"),quote=F,sep="\t",col.names=NA)

}

for (i in as.numeric(levels(dg.integrated$seurat_clusters))){		
	# Instantiate all variables to make sure no carry-overs from previous clusters
	names = 0
	nullnames = 0
	ctl1names = 0
	diff = 0
	sig.up = 0
	sig.dn = 0
	sig.up.gp = 0
	sig.dn.gp = 0
	sig.up.gp.flat = 0
	sig.dn.gp.flat = 0

	print(paste0("Working on cluster ",i))

	diff = read.table(paste0("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_chr2_vs_yfp_cluster",i,"_DEresults_061821.txt"),header=T,sep="\t",row.names=1)

	sig.up = rownames(diff[diff$p_val_adj<0.05 & diff$avg_diff > 0,])
	sig.dn = rownames(diff[diff$p_val_adj<0.05 & diff$avg_diff < 0,])

	if(length(sig.up)>=5 & length(sig.dn)>=5) {
		q = list(sig.up, sig.dn)
		names(q) = c(paste0("Cluster",i,"_UP"),paste0("Cluster",i,"_DN"))
		gp = gost(q,significant=TRUE,organism = "mmusculus",evcodes=TRUE,correction_method="fdr",sources=c("GO:BP","GO:MF","GO:CC","KEGG","REAC"),multi_query=FALSE)
		gp.flat = as_tibble(gp$result) %>%
			select(-parents,-evidence_codes)	
		write.table(gp.flat,file=paste0("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_chr2_vs_yfp_cluster",i,"_gProfiler_061821.txt"),quote=F,sep="\t",col.names=NA)
	} else if(length(sig.up)>=5) {
		q = list(sig.up)
		names(q) = paste0("Cluster",i,"_UP")
		gp = gost(q,significant=TRUE,organism = "mmusculus",evcodes=TRUE,correction_method="fdr",sources=c("GO:BP","GO:MF","GO:CC","KEGG","REAC"),multi_query=FALSE)
		gp.flat = as_tibble(gp$result) %>%
			select(-parents,-evidence_codes)	
		write.table(gp.flat,file=paste0("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_chr2_vs_yfp_cluster",i,"_gProfiler_061821.txt"),quote=F,sep="\t",col.names=NA)
	} else if(length(sig.dn)>=5) {
		q = list(sig.dn)
		names(q) = paste0("Cluster",i,"_DN")
		gp = gost(q,significant=TRUE,organism = "mmusculus",evcodes=TRUE,correction_method="fdr",sources=c("GO:BP","GO:MF","GO:CC","KEGG","REAC"),multi_query=FALSE)
		gp.flat = as_tibble(gp$result) %>%
			select(-parents,-evidence_codes)	
		write.table(gp.flat,file=paste0("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_chr2_vs_yfp_cluster",i,"_gProfiler_061821.txt"),quote=F,sep="\t",col.names=NA)
	} else {
		# skip
	}
}




# Make umap split by treatment

pdf("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_splitByTreatment.pdf",width=12)
DimPlot(dg.integrated,reduction = "umap", split.by = "trt",label=T)
dev.off()



# Read in DE output files and make barplot of number of DE genes per cluster

counts = as.data.frame(matrix(data=NA,nrow=length(unique(levels(dg.integrated@meta.data$seurat_clusters))),ncol=2,byrow=T))
colnames(counts) = c("Upregulated","Downregulated")
rownames(counts) = unique(levels(dg.integrated@meta.data$seurat_clusters))
counts = cbind(counts,"Group" = c("GC","GC","GC","GC","GC","GC","Excitatory","GC","GC","GABA","GC","GC","Excitatory","Excitatory","GC","NSC","NSC","Mossy","GABA","GABA","GABA","Mossy","GC","GABA","Mossy","NSC","GC","Oligo","Excitatory","Endothelial","GC","VLMC","GC","DG Meis2","Excitatory","Microglia","Astrocyte","Cajal-Retzius","GABA","Neuron","GABA","Excitatory","Ependymal"))

for (i in unique(levels(dg.integrated@meta.data$seurat_clusters))) {
	file = paste0("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_chr2_vs_yfp_cluster",i,"_RNA_scaledata_DEresults_061821.txt")
	if(file.exists(file)) {
		upct = length(read_tsv(file,col_names=c("Gene","p_val","avg_diff","pct.1","pct.2","p_val_adj"),skip=1) %>%
			filter(p_val < 0.01) %>%
			filter(avg_diff > 0) %>%
			pull(Gene))
		dnct = length(read_tsv(file,col_names=c("Gene","p_val","avg_diff","pct.1","pct.2","p_val_adj"),skip=1) %>%
			filter(p_val < 0.01) %>%
			filter(avg_diff < 0) %>%
			pull(Gene))
		counts[i,1] = upct
		counts[i,2] = dnct		
	}
}

rownames_to_column(counts,var="Cluster") %>% 
	as_tibble() %>%
	pivot_longer(cols=contains("regulated"),names_to="Direction",values_to="Count") %>%
	mutate(
		Count = case_when(
			Direction=="Downregulated" & !is.na(Count) ~ -1*Count,
			!is.na(Count) ~ 1*Count,
			T ~ 0
		)
	) %>%
	mutate(Cluster = as.numeric(Cluster)) %>%
	ggplot(aes(x=Cluster,y=Count,fill=Group)) +
	geom_col() +
	ylab("Downregulated count    Upregulated count") +
	geom_hline(yintercept=0)




# Now plot sum of total by principal cell type	
	
rownames_to_column(counts,var="Cluster") %>% 
	as_tibble() %>%
	pivot_longer(cols=contains("regulated"),names_to="Direction",values_to="Count") %>%
	mutate(
		Count = case_when(
			Direction=="Downregulated" & !is.na(Count) ~ 1*Count,
			!is.na(Count) ~ 1*Count,
			T ~ 0
		)
	) %>%
	mutate(Cluster = as.numeric(Cluster)) %>%
	group_by(Group) %>%
	summarize(Sum = sum(Count)) %>%
	ggplot(aes(x=reorder(Group,-Sum),y=Sum,fill=Group)) +
	geom_col() +
	ylab("Total DEGs")
