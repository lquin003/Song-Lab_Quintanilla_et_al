assign(".lib.loc", "/nas/longleaf/home/lquin003/R/alternate_libs/", envir = environment(.libPaths))


library(ggplot2,lib.loc = "/nas/longleaf/rhel8/apps/r/4.1.0/lib64/R/library/")
library(RColorBrewer,lib.loc = "/nas/longleaf/rhel8/apps/r/4.1.0/lib64/R/library/")
library(stringr ,lib.loc = "/nas/longleaf/rhel8/apps/r/4.1.0/lib64/R/library/")
library(data.table,lib.loc = "/nas/longleaf/rhel8/apps/r/4.1.0/lib64/R/library/")
library(tidyr,lib.loc = "/nas/longleaf/rhel8/apps/r/4.1.0/lib64/R/library/")
library(pdftools,lib.loc = "/nas/longleaf/rhel8/apps/r/4.1.0/lib64/R/library/")
library(readxl,lib.loc = "/nas/longleaf/rhel8/apps/r/4.1.0/lib64/R/library/")
library(UpSetR,lib.loc = "/nas/longleaf/home/lquin003/R/x86_64-pc-linux-gnu-library/4.1")
library(patchwork,lib.loc = "/nas/longleaf/home/lquin003/R/x86_64-pc-linux-gnu-library/4.1")
library(scales,lib.loc = "/nas/longleaf/rhel8/apps/r/4.1.0/lib64/R/library/")
library(tibble,lib.loc = "/nas/longleaf/rhel8/apps/r/4.1.0/lib64/R/library/")
library(forcats,lib.loc = "/nas/longleaf/rhel8/apps/r/4.1.0/lib64/R/library/")
library(dplyr)
library(janitor)
library(ggdendro)
library(circlize,lib.loc = "/nas/longleaf/rhel8/apps/r/4.1.0/lib64/R/library/")
library(ComplexHeatmap,lib.loc = "/nas/longleaf/rhel8/apps/r/4.1.0/lib64/R/library/")
library(ggrepel)
library(gprofiler2)
library(Seurat,lib.loc = "/nas/longleaf/rhel8/apps/r/4.1.0/lib64/R/library/")
library(Matrix,lib.loc = "/nas/longleaf/rhel8/apps/r/4.1.0/lib64/R/library/")
library(presto,lib.loc = "/nas/longleaf/home/lquin003/R/x86_64-pc-linux-gnu-library/4.1")
library(cowplot,lib.loc = "/nas/longleaf/home/lquin003/R/x86_64-pc-linux-gnu-library/4.1")


#### Begin Analysis after preprocessing steps this script Assumes Seurat Object has been created, cells have been renamed,
#### clusters have been scaled, split, integrated between control and treatment, regressed by top vars, ran PCA, ran UMAP,
#### Found neighbors found clusters and have iterated through for optimum cluster resolution object should have SCT, RNA, and integrated as assays ###

#set output for graphs and other files use this paste0 output before each file name to use this
output <- "/nas/longleaf/home/lquin003/Documents/Output/"

#import cluster names from word file copy and paste
cluster_names <- read.csv("/nas/longleaf/home/lquin003/Seurat/Cluster_Annotations_R2.txt")
print(cluster_names)
cluster_names <- cluster_names[1:43,]
cluster_names <- as.character(cluster_names)
#### Start analysis for DEG's from here. #object is 11.5 GB, make sure to have enough memory if using interactive R studio ####
dg.cpm = readRDS(file ="/proj/jsonglab/projects/Luis/single_cell_analysis/Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_wrs2_deg_cpm_new.rds")


# Assign order for future Vln plots contorl and treatment
dg.cpm$trt <- factor(dg.cpm$trt, levels = c("yfp","chr2"))
dg.cpm$trt

# Add cluster names to identties 
names(cluster_names) <- levels(dg.cpm)
dg.cpm <- RenameIdents(dg.cpm, cluster_names)


########## Correlation Plot from Jeremy 
#dg.integrated = readRDS("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered.rds")

library(presto)

vargenes<- presto::wilcoxauc(dg.cpm, 'seurat_clusters', seurat_assay = 'SCT')
top_vargenes = top_markers(vargenes, n = 50, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)
all_markers<- top_vargenes %>%
  dplyr::select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

#dg.integrated <- BuildClusterTree(dg.integrated, reorder = F, verbose = T, assay="SCT", features = all_markers)
#PlotClusterTree(dg.integrated)

dg.avg = AverageExpression(dg.cpm)
dg.sct.cor = cor(dg.avg$SCT[all_markers,])

col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
col_fun(seq(-1, 1))

fill=list(1:43)
#df_anno = data.frame(cluster_names,fill) %>% rename("colors"="X1.43")

#df_anno

#ha=HeatmapAnnotation(cell_class =df_anno$colors)



cor_heatmap = ComplexHeatmap::Heatmap(dg.sct.cor,
                                      name="Pearson\nCorrelation",heatmap_height = unit(20,"cm"),heatmap_width = unit(25,"cm"),
                                      column_names_side = "top",
                                      row_names_side = "left",
                                      col = col_fun, 
                                      show_column_dend = F,
                                      show_column_names = F,
                                      border = T,
                                      rect_gp = gpar(col= "black"),
                                      #top_annotation = ha,
                                      clustering_distance_columns=function(x) as.dist(1-cor(t(x))), 
                                      clustering_distance_rows=function(x) as.dist(1-cor(t(x))), 
                                      clustering_method_columns = "average",
                                      heatmap_legend_param = list(direction="horizontal",title_position = "topcenter"), 
                                      clustering_method_rows = "average")


heatmap_row_order = ComplexHeatmap::row_order(ComplexHeatmap::draw(cor_heatmap))-1


pdf(paste0(output,"Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_prestoMarkers_SCT_pearson_cor_heatmap__dim_50_4.0.pdf"),height=8,width=12)
draw(cor_heatmap, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()


p1= grid.grabExpr(draw(cor_heatmap, heatmap_legend_side = "bottom", annotation_legend_side = "bottom"))



########### Dot plot visualization of All genes Dot plot has been changed to Violin Plot ######### Figure 5B
feature_list <- c("Sox2","Gfap","Aqp4","Clu","Hopx","Cd74","Csf1R","Plp1","Gad2","Gad1","Dpp10","Dcx","Prox1","Sv2b",
                  "Flt1","Grik1","Ctss","Cst3","C1qa","Pex5l","Arhgap12","Fibcd1","Cpne4","Ly6e","Npy2R",
                  "MBP", "Plp1","Mobp","Pdgfra","Vcan","Cspg4","Amigo2","Cntn5","Calb2","Calb1",
                  "Foxg1","Nes","Olig","Nrxn1","Reln","Dkk3","Prkcb")

total_feature_list <- c("Gfap","Apoe","Aqp4","Clu","Hopx","Cdh5","Foxj1","Dsp","Gad2","Gad1","Slc17a7","Calb1","Prox1","Sv2b","Cntn6","Grm8","Cst3","Pex5l","Arhgap12","Mobp","Pdgfra","Snap25","Foxp2","Cacng5","Htr2c","Abi3bp","Ccnd2","Bhlhe22","Gabra5","Prrx1")

total_feature_list2 = c("Mapk8","Itgb1","Itgb5","Itgb8","Fads2","Sorl1")

total_dp <- dg.cpm %>% 
  DotPlot(feature=total_feature_list2,cluster.idents = T, assay = "SCT",cols = c("grey","red","blue"))+
  labs(x="Genes",y="Clusters")+
  theme(legend.position = "none",
        axis.line = element_blank())
total_dp

pdf(paste0(output,"total_dotplot_SCT_intracellular_receptors.pdf") ,height =8 ,width =18)
total_dp
dev.off()

#astrocyte_genes = c("Sox9","Etv4","Sall3","Grin2c","Kcng4","Aqp4","Gfap","Hopx","S100b","mki67","Pcna","Sox2","Tbr2","Mcm2","Dcx","Ncam1","Mxd3","Ascl1","Neurog2","Rfp4","Lockd","Eomes","Calb2","Neurod4","Slc17a6")
#Total Features Vln plot stacked

total_feature_vlnplot <- dg.cpm %>%
  VlnPlot(object =, assay = "SCT",slot="data",sort = T,features = total_feature_list,stack = T,flip = T)+ #facet_wrap(~total_feature_list,strip.position = "left") +
  NoLegend() +
  plot_annotation() & theme(strip.text = element_text(face = "italic"),
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.line = element_blank(), axis.text.y = element_blank(),
                            axis.text.x = element_text(angle = 90,vjust = .5,hjust = 1),
                            axis.title.y.right = element_blank(),
                            axis.text.y.right = element_blank())

pdf(paste0(output,"total_vlnplot_SCT.pdf") ,height =8 ,width =18)
total_feature_vlnplot
dev.off()



#### Attempt at heatmap row order ### it works
my_levels= heatmap_row_order
dg.cpm@active.ident = dg.cpm$seurat_clusters
Idents(dg.cpm) <- factor(Idents(dg.cpm), levels= my_levels)

side_total_feature_vlnplot <- dg.cpm %>%
  VlnPlot(object =, assay = "SCT",slot = "data",features = total_feature_list ,stack = T)+ 
  NoLegend() + scale_y_discrete(limits=rev) + plot_annotation() & theme(strip.text = element_text(face = "italic"),
                                                                        axis.title = element_blank(),
                                                                        axis.text.x = element_blank(),
                                                                        axis.ticks.x = element_blank())

side_total_feature_vlnplot

## Combine with complex heatmap
p1 = grid.grabExpr(draw(cor_heatmap, heatmap_legend_side = "top", annotation_legend_side = "top"))
p2 = side_total_feature_vlnplot

#output for figure 5B
pdf(paste0(output,"Cor_plot_vln_plotmerge_fig_5b.pdf"),height =9,width = 20)
cowplot::plot_grid(p1,p2,nrow=1, axis="r",rel_heights = c(.9,1),rel_widths = c(1,1))
dev.off()




#### Figure 5 C UMAP Plots #### 
################# UMap Split by treatement colors
ach$trtrep <- factor(ach$trtrep,levels=c("yfp-Rep1",
                                         "yfp-Rep2",
                                         "yfp-Rep3",
                                         "yfp-Rep4",
                                         "yfp-Rep5",
                                         "chr2-Rep1",
                                         "chr2-Rep2",
                                         "chr2-Rep3",
                                         "chr2-Rep4",
                                         "chr2-Rep5"))


ctrl_vs_trtp <- UMAPPlot(ach,group.by="trtrep")+
  scale_color_manual(labels=c("C1","C2","C3","C4",'C5',"T1","T2","T3","T4","T5"),
                     values = c( "yfp-Rep1"="green",
                                 "yfp-Rep2"="green1",
                                 "yfp-Rep3"="green2",
                                 "yfp-Rep4"="green3",
                                 "yfp-Rep5"="green4",
                                 "chr2-Rep1"="deepskyblue",
                                 "chr2-Rep2"="deepskyblue1",
                                 "chr2-Rep3"="deepskyblue2",
                                 "chr2-Rep4"="deepskyblue3",
                                 "chr2-Rep5"="deepskyblue4"))+
  labs(title="Biological Replicates") +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank())

pdf(file= paste0(output,"controlvstreatUMAP.pdf"), width = 6, height=6)
ctrl_vs_trtp
dev.off()


cell_location_UMAP =  DimPlot(ach, group.by ="cell_location") +labs(title = "Cell Location")+scale_color_manual(values = c("blue3","darkorange","red3","purple3")) +theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())

cell_location_UMAP


cell_type_UMAP = DimPlot(ach, group.by ="cell_type") +labs(title = "Cell Type")+scale_color_manual(values = c("blue3","darkorange","red3","purple3"))  +theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
cell_type_UMAP


pdf(paste0(output,"Representative_UMAPs.pdf"), width=18,height=6)
ctrl_vs_trtp+cell_location_UMAP+cell_type_UMAP
dev.off()

###Supplemental Figure 6B UMAP with names####

#Make sure to rename clusters accordingly
names(cluster_names) <- levels(dg.cpm)
dg.cpm <- RenameIdents(dg.cpm, cluster_names)

#determine length of unique clusters
colorcount <- length(unique(dg.cpm$seurat_clusters))
# function for color ramp
getPalette = colorRampPalette(brewer.pal(12, "Set3"))

## Make UMAP plot that has cluster labels split by treatment
umap_ach <- DimPlot(dg.cpm,label=F) + NoLegend()+
  scale_color_manual(values = getPalette(colorcount))+
  plot_annotation() & theme(axis.line = element_blank(),
                            axis.title = element_blank(),
                            axis.ticks = element_blank(),
                            axis.text = element_blank())

#add labels to the UMAP plot and use repel to remove label overlap

pdf(paste0(output,"Supplemental Fig 6B UMAP with Names.pdf"),width = 10,height = 10)
LabelClusters(umap_ach, id="ident",repel = T,size=4)
dev.off()

#### Start of analysis for Figure 6 ####

#### Find DEGs using RNA assay and counts per million on Seurat object ####
setwd("/proj/jsonglab/projects/Luis/single_cell_analysis/Output/")
dg.integrated = readRDS("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_wrs2_deg.rds")
dg.sct = SCTransform(dg.integrated, vars.to.regress = "Percent.mt",verbose =T)

dg.cpm = NormalizeData(dg.integrated,normalization.method = "RC", scale.factor = 1e6, verbose=T, assay="RNA")
dg.cpm <- ScaleData(dg.cpm, assay="RNA")

for (i in as.numeric(levels(dg.cpm$seurat_clusters))){
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
  
  names = colnames(dg.cpm)[dg.cpm$seurat_clusters==i]
  nullnames = names[grepl("chr2",names)]
  ctl1names = names[grepl("yfp",names)]
  if((length(nullnames)>0) & (length(ctl1names)>0)) {
    diff = FindMarkers(dg.cpm, ident.1 = nullnames, ident.2=ctl1names, assay="RNA", slot="scale.data",logfc.threshold = 0.1)
    write.table(diff,paste0("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_chr2_vs_yfp_cluster",i,"_RNA_scaledata_DEresults_061821.txt"),quote=F,sep="\t",col.names=NA)
  }
  
}

# Save normalized RDS file for analysis
saveRDS(dg.cpm, file ="/proj/jsonglab/projects/Luis/single_cell_analysis/Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_wrs2_deg_cpm_new.rds")

#### Start analysis for DEG's from here. #object is 11.5 GB, make sure to have enough memory if using interactive R studio ####
dg.cpm = readRDS(file ="/proj/jsonglab/projects/Luis/single_cell_analysis/Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_wrs2_deg_cpm_new.rds")

# Assign order for future Vln plots contorl and treatment
dg.cpm$trt <- factor(dg.cpm$trt, levels = c("yfp","chr2"))
dg.cpm$trt

# Add cluster names to identties 
names(cluster_names) <- levels(dg.cpm)
dg.cpm <- RenameIdents(dg.cpm, cluster_names)


######## Import DEG Data CSV Import and Cluster Name ###### ########
## Read in differential expression data from directory for WRS method 2
setwd("/nas/longleaf/home/lquin003/Seurat/DE_TestingRNA/")
dedirectory2 <- setwd("/nas/longleaf/home/lquin003/Seurat/DE_TestingRNA/")

files_list <- list.files(dedirectory2)

#take only files that have DEresults
deonlyfiles2 <- str_sort(files_list[grep("scaledata_DEresults",x = files_list,ignore.case = T)],numeric=T)
number_of_files <- length(deonlyfiles2)
print(number_of_files)

#Read cluster numbers from file list and determine cluster numbers
cluster_numbers <- str_match(deonlyfiles2,regex("cluster*.\\d"))
cluster_number_only <- matrix(as.numeric((gsub("\\D","",cluster_numbers))))

#import cluster names from word file copy and paste
cluster_names <- read.csv("/nas/longleaf/home/lquin003/Seurat/Cluster_Annotations_R2.txt")
print(cluster_names)
cluster_names <- cluster_names[1:43,]
cluster_names <- as.character(cluster_names)

#Empty Cluster names
cluster_names_mat <- matrix()
for(i in seq_along(cluster_names)){
  if(as.numeric(gsub("\\D","",unlist(str_extract(deonlyfiles2[i],"cluster*.\\d")))) %in% as.numeric(str_extract(cluster_names[i],".\\d?"))){
    cluster_names_mat[i] <- cluster_names[i]
  }  else{
    deonlyfiles2<- append(deonlyfiles2,NA,after=i-1)
  }
}

#remove NA used in list extension to prevent issues later on
cluster_names_mat <-  na.exclude(cluster_names_mat)
deonlyfiles2 <-  na.exclude(deonlyfiles2)
print(cluster_names_mat)
print(deonlyfiles2)

#read in each file by searching names of file, then add cluster number and cluster name.
#For all other clusters other than the initial continue to add to original data frame
for (i in 1:number_of_files){
  if(i == 1){
    datatemp <- read.csv(deonlyfiles2[i],sep = "",header = T,row.names = NULL)
    datatemp$cluster_number <- cluster_number_only[i]
    datatemp$cluster_name <- str_replace(cluster_names_mat[i],".\\d?\\s","")
    data <- datatemp
  } 
  #If it's not the first file then merge the list with past list to create large data frame
  else
    if(i>1 ){
      datatemp <- read.csv(deonlyfiles2[i],sep = "",header = T,row.names = NULL)
      datatemp$cluster_number <- cluster_number_only[i]
      datatemp$cluster_name <- str_replace(cluster_names_mat[i],".\\d?\\s","")
      data <- rbind(data,datatemp)
    }
}

data <- as.data.frame(data)

#rename dataframe first row as gene
colnames(data)[1] <- "gene"

#subset only significant data less than .05 adjusted p_value
significant_data2 <- na.omit(data[data$p_val < .01,])

## Remove RNA genes and Xist a sex linked genes from analysis
rna_genes<- c("Rn18s","Meg3","Snhg11","Xist")
significant_data2 <- significant_data2[!grepl(paste(rna_genes,collapse ="|"),significant_data2$gene),]

#rename data
wrs_method2_data <- significant_data2
wrs_method2_data$analysis_method <- "wrs_2"

str(wrs_method2_data)

#DEG set 2
wrs_method2_data
write.csv(wrs_method2_data, paste0(output,"wrs2_DEG.csv"))

###### DE Heatmap for top 5 genes and bottom 5 genes average Diff ###### #####################
Granule_cells <- c(0,1,2,3,4,5,7,8,10,11,14,26,30,32)
GABA <- c(18,19,20,23,40)
Mossy <- c(17,21,24)
Subiculum <- c(9,33,38,39,41)
CA_Excitatory <- c(6,12,13,22,28,34)
Neurogenic <- c(15,25,16)
Glia <- c(27,31,35,36)
Other_notdg <- c(37,42)
Other<- c(29)
RGL <-(15)

#Select each gene category based on clusters and assign cell type. This will be then used when making heatplot
select_cell_categories_in_df_plus_column <- function(data_frame, cell_category1){
  out <- as_tibble(data_frame) %>%
    filter(cluster_number %in% cell_category1)
  add_col <- deparse(substitute(cell_category1))
  out <- out %>% tibble::add_column("cell_type" = add_col)
  return(out)
}

# must provide 2 variables indicating data frame and cell category above
mossy_genes <- select_cell_categories_in_df_plus_column(wrs_method2_data,Mossy)
gc_genes <- select_cell_categories_in_df_plus_column(wrs_method2_data,Granule_cells)
neurogenic_genes <- select_cell_categories_in_df_plus_column(wrs_method2_data,Neurogenic)
gaba_genes <- select_cell_categories_in_df_plus_column(wrs_method2_data,GABA)
glia_genes <- select_cell_categories_in_df_plus_column(wrs_method2_data,Glia)
ca_genes <- select_cell_categories_in_df_plus_column(wrs_method2_data,CA_Excitatory)
sub_genes <- select_cell_categories_in_df_plus_column(wrs_method2_data,Subiculum)
other_notdg <- select_cell_categories_in_df_plus_column(wrs_method2_data,Other_notdg)
other_genes <- select_cell_categories_in_df_plus_column(wrs_method2_data,Other)

# Create tibble by binding all cell type data frames then removing CA excitatory and subiculum genes and then sorting by highest avg_diff
heatmapdf <- rbind(mossy_genes,gc_genes,neurogenic_genes,gaba_genes,glia_genes,sub_genes,ca_genes,other_notdg,other_genes) %>% 
  filter(!cell_type %in% c("CA_Excitatory","Subiculum","Other_notdg")) %>%
  arrange(-avg_diff)


##### Create Bar plot Supplemental Figure included in DEGS to quantify overalap of DEG's across difrerent clusters #####

## List of overlapped genes in data set 
str(heatmapdf)
#tibble [2,676 × 10]
total_deg = 2676

#Find duplicates using janitor package get dupes creates new df Subset any gene that is greater than 1
percent_plot = heatmapdf %>% get_dupes(gene) %>% group_by(gene,cell_type) %>% filter(length(cell_type) >1) 
#Table all DEG's and make data frame for all DEGs
counts_df = as.data.frame(table(heatmapdf$cell_type)) %>% dplyr::rename("Cell_class"="Var1", "Values1"="Freq")
#Table filtered duplicate genes across cell type and make df
counts_df2 =  as.data.frame(table(percent_plot$cell_type))  %>% dplyr::rename("Cell_class"="Var1", "Values2"="Freq")
#Joing both DF's and create column that indicates percentage
deg_percent_by_class = full_join(counts_df,counts_df2,"Cell_class") %>% mutate("percent"=round((Values2/Values)*100,digits = 2))
#arrange df from highest to lowest percent of overlap
deg_percent_by_class =arrange(deg_percent_by_class,-percent)

#Create bar plot that contains all overlapped percentages
bar = ggplot(deg_percent_by_class, aes(x=fct_reorder(Cell_class,-percent) ,y=percent, fill=Cell_class)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label = paste0(percent, "%")),position = position_stack(vjust=0.5),check_overlap = T) +
  labs(y ="Percent of DEG's that overlap Across Cell Class") +
  theme_classic() +
  #labs(title = "Percent of DEG's Overlaping in Cell Class")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        legend.position = "none",) +
  scale_fill_brewer(palette = "Set3",breaks=c("Granule_cells","GABA","Mossy","Glia","Other"))

#Make pdf with barplot
pdf(paste0(output,"Overlapping DEG's Bar Chart.pdf"),height = 7,width = 7)
bar
dev.off()

#write all overlapped genes in CSV table
write.csv(x=overlapped_genes,paste0(output,"all overlapped genes.csv"))
#write perecent overlapped genes in csv table
write.csv(x=percent_plot,paste0(output,"Overlapped_Genes_by_Cell_Class.csv"))



#Create bar plot that contains all overlapped percentages
bar = ggplot(deg_percent_by_class, aes(x=fct_reorder(Cell_class,-percent) ,y=percent, fill=Cell_class)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label = paste0(percent, "%")),position = position_stack(vjust=0.5),check_overlap = T) +
  labs(y ="Percent of DEG's that overlap Across Cell Class") +
  theme_classic() +
  #labs(title = "Percent of DEG's Overlaping in Cell Class")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        legend.position = "none",) +
  scale_fill_brewer(palette = "Set3",breaks=c("Granule_cells","GABA","Mossy","Glia","Other"))


### Create Pie chart ## Make sure to load data frame m before running through pie chart creation
piechart=  m %>% 
  mutate("total"= Up+Dwn) %>%   
  mutate("total_sum" = sum(.[[6]])) %>% 
  mutate("total_percent" = total/total_sum * 100) %>% 
  group_by(cell_class) %>% 
  mutate("cell_type_sum"=sum(total_percent)) %>% 
  ungroup() %>% 
  distinct(cell_class,cell_type_sum) %>%
  mutate(cell_class = fct_reorder(cell_class,-cell_type_sum))%>%
  ggplot(aes(x = "", y = round(cell_type_sum,2), fill = cell_class)) +
  geom_col(color = "black") +
  geom_text(aes(label = round(cell_type_sum,2)),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  labs(title="Percent of DEG's by Cell Class")+
  scale_fill_brewer(palette = "Set3") +
  scale_fill_discrete(breaks = c("Granule Cell","GABA","Mossy","Neurogenic","Glia","Other"))+
  theme_void()

pdf(paste0(output,"Piechart_fig6_degpercent.pdf"),height = 4,width = 4)
piechart
dev.off()

# looking for percent of DEG overlap between 2 GC cells
#percent_plot = heatmapdf %>% filter(cluster_number %in% c(0,7)) %>%  get_dupes(gene) %>% group_by(gene,cell_type) %>% filter(length(cell_type) >1) 

#### Figure 6A and Figure 6B #####
######### Count total DEG's by cluster ############

# count upregualted and downregulated DEGs
upcounts_new = heatmapdf %>% filter(avg_diff >0) %>% group_by(cluster_name) %>% dplyr::count(cluster_name) %>%mutate("up"="up")
# count upregualted and downregulated DEGs
dwncounts_new <- heatmapdf %>% filter(avg_diff < 0) %>% group_by(cluster_name)%>% dplyr::count(cluster_name) %>% mutate("dwn" = "dwn")

#Make DF that contains the ratio of upregualted and wonregualted DEGS
m <- full_join(by = "cluster_name",upcounts_new,dwncounts_new) %>% dplyr::rename("Up"="n.x","Dwn"="n.y") %>% dplyr::select(!c(up,dwn)) 
m = m %>% replace(is.na(m),0) %>% mutate(ratio = ((Up-Dwn)/(Up+Dwn)))

#Select cells for each class to make list using cluster names
gc_list = m[grepl("GC", m$cluster_name),] %>% pull(cluster_name)  
gc_list = gc_list[!grepl("Immature",gc_list)]
mossy_list = m[grepl("Mossy", m$cluster_name),]%>% pull(cluster_name)
gaba_list = m[grepl("^DG.*?GABA", m$cluster_name),]%>% pull(cluster_name)
glia_list = m[grepl("Astrocyte|Microglia|Oligo|OPC",m$cluster_name),] %>% pull(cluster_name)
neurogenic_list = m[grepl("NSC|IPC|Immature",m$cluster_name),] %>% pull(cluster_name)

#add column with cell class to counts df
m =  add_column(m, "cell_class" = case_when(m$cluster_name %in% gc_list ~ "Granule Cell",
                                            m$cluster_name %in% mossy_list ~ "Mossy",
                                            m$cluster_name %in% gaba_list ~ "GABA",
                                            m$cluster_name %in% glia_list ~ "Glia",
                                            m$cluster_name %in% neurogenic_list ~ "Neurogenic",
                                            TRUE~"Other"))




#subset dataframe by cell class and then create bar plots used for ratios
ratio_plot_gc = m %>% filter(cell_class=="Granule Cell") %>% ungroup() %>%  ggplot(aes(x = fct_reorder(cluster_name,ratio),y=ratio))+ geom_bar(stat = "identity",col="black",alpha=.8, aes(fill=ratio<0)) + scale_fill_manual(guide=F,breaks=c(T,F),values = c("steelblue","firebrick"))+theme(axis.title.x = element_blank(),legend.position = "none",panel.background  = element_blank(), axis.text.y= element_text(size = 12), axis.text.x = element_text(angle =0,size = 12)) + labs(x="",title="Ratio of Up and Down DEGs") + coord_flip() +lims(y=c(-1,1))+ geom_hline(yintercept=0)
ratio_plot_gaba = m %>% filter(cell_class=="GABA") %>% ungroup()  %>%  ggplot(aes(x = fct_reorder(cluster_name,ratio),y=ratio))+ geom_bar(stat = "identity",col="black",alpha=.8,aes(fill=ratio<0)) + scale_fill_manual(guide=F,breaks=c(T,F),values = c("steelblue","firebrick"))+theme(axis.title.x = element_blank(),legend.position = "none",panel.background  = element_blank(), axis.text.y= element_text(size = 12), axis.text.x = element_text(angle =0,size = 12)) + labs(x="") + coord_flip()+lims(y=c(-1,1))+ geom_hline(yintercept=0)
ratio_plot_mossy = m %>% filter(cell_class=="Mossy") %>% ungroup() %>%  ggplot(aes(x = fct_reorder(cluster_name,ratio),y=ratio))+ geom_bar(stat = "identity",col="black",alpha=.8,aes(fill=ratio<0)) + scale_fill_manual(guide=F,breaks=c(T,F),values = c("steelblue","firebrick"))+theme(axis.title.x = element_blank(),legend.position = "none",panel.background  = element_blank(), axis.text.y= element_text(size = 12), axis.text.x = element_text(angle =0,size = 12)) + labs(x="") + coord_flip()+lims(y=c(-1,1))+ geom_hline(yintercept=0)
ratio_plot_ng = m %>% filter(cell_class=="Neurogenic") %>%  ggplot(aes(x = fct_reorder(cluster_name,ratio),y=ratio))+ geom_bar(stat = "identity",col="black",alpha=.8,aes(fill=ratio<0)) + scale_fill_manual(guide=F,breaks=c(T,F),values = c("steelblue","firebrick"))+theme(axis.title.x = element_blank(),legend.position = "none",panel.background  = element_blank(), axis.text.y= element_text(size = 12), axis.text.x = element_text(angle =0,size = 12) ) + labs(x="") + coord_flip()+lims(y=c(-1,1))+ geom_hline(yintercept=0)
ratio_plot_glia = m %>% filter(cell_class=="Glia") %>% ungroup() %>%  ggplot(aes(x = fct_reorder(cluster_name,ratio),y=ratio))+ geom_bar(stat = "identity",col="black",alpha=.8,aes(fill=ratio<0)) + scale_fill_manual(guide=F,breaks=c(T,F),values = c("steelblue","firebrick"))+theme(axis.title.x = element_blank(),legend.position = "none",panel.background  = element_blank(), axis.text.y= element_text(size = 12), axis.text.x = element_text(angle =0,size = 12)) + labs(x="") + coord_flip()+lims(y=c(-1,1))+ geom_hline(yintercept=0)
ratio_plot_other = m_other_r = m %>% filter(cell_class=="Other") %>% ungroup()  %>%  ggplot(aes(x = fct_reorder(cluster_name,ratio),y=ratio))+ geom_bar(stat = "identity",col="black",alpha=.8,aes(fill=ratio<0)) + scale_fill_manual(guide=F,breaks=c(T,F),values = c("steelblue","firebrick"))+theme(axis.title.x = element_blank(),legend.position = "none",panel.background  = element_blank(), axis.text.y= element_text(size = 12), axis.text.x = element_text(angle =0,size = 12)) + labs(x="") + coord_flip()+lims(y=c(-1,1))+ geom_hline(yintercept=0)



#Figure 6B bar pltos stacked ###(Removed from Manuscript) #### 
pdf(paste0(output,"DEG ratio plots.pdf"),height = 10,width = 5)
(ratio_plot_gc+ratio_plot_gaba+ratio_plot_mossy+ratio_plot_ng+ratio_plot_glia+ratio_plot_other) +plot_layout(heights = c(3.2,1.5,1,1,1.5,.5),ncol=1)
dev.off()


# CSV with coutns 
write.csv(m,paste0(output,"cluster_counts_deg_new_merged.csv"))


#multiply DF by negative number to show downregulated DEGs in negative direction
m$Dwn = m$Dwn*-1

#Pivot long for graphing
m =  pivot_longer(m,c(Up,Dwn))

#set order for up and down
m$name= factor(m$name,levels = c("Up","Dwn"))

#Filter by cell class for DEG Counts Fig 6A
m_gc = m %>% filter(cell_class=="Granule Cell") %>% ungroup() %>% mutate("cluster_name" = fct_reorder(cluster_name, value,.fun =max ,.desc = F))
m_gaba = m %>% filter(cell_class=="GABA") %>% ungroup() %>% mutate("cluster_name" = fct_reorder(cluster_name, value,.fun =max ,.desc = F))
m_mossy = m %>% filter(cell_class=="Mossy") %>% ungroup() %>% mutate("cluster_name" = fct_reorder(cluster_name, value,.fun =max ,.desc = F))
m_neurogenic = m %>% filter(cell_class=="Neurogenic") %>% ungroup() %>% mutate("cluster_name" = fct_reorder(cluster_name, value,.fun =max ,.desc = F))
m_glia = m %>% filter(cell_class=="Glia") %>% ungroup() %>% mutate("cluster_name" = fct_reorder(cluster_name, value,.fun =max ,.desc = F))
m_other = m %>% filter(cell_class=="Other") %>% ungroup() %>% mutate("cluster_name" = fct_reorder(cluster_name, value,.fun =max ,.desc = F))

# plot each cell class as a seperate bar plot
gc_p = ggplot(m_gc,aes(x=cluster_name,y=value,fill=name))+ geom_bar(stat = "identity",alpha=.8,col="black")+labs(x="",y="",title = "Number of Up and Down DEGs",fill="Direction \nof DEG")+coord_flip()+theme(panel.background = element_blank() ,axis.text.x = element_text(size = 12),axis.text.y = element_text(size=12),axis.title.y = element_blank(),legend.position = "none")+ scale_fill_manual(values = c("firebrick","steelblue")) +lims(y=c(-250,250))+ geom_hline(yintercept=0)
gaba_p = ggplot(m_gaba,aes(x=cluster_name,y=value,fill=name))+ geom_bar(stat = "identity",alpha=.8,col="black")+labs(x="",y="",fill="Direction \nof DEG")+coord_flip()+theme(panel.background = element_blank() ,axis.text.x = element_text(size = 12),axis.text.y = element_text(size=12),axis.title.y = element_blank(),legend.position = "none")+ scale_fill_manual(values = c("firebrick","steelblue")) +lims(y=c(-250,250)) + geom_hline(yintercept=0)
mossy_p = ggplot(m_mossy,aes(x=cluster_name,y=value,fill=name))+ geom_bar(stat = "identity",alpha=.8,col="black")+labs(x="",y="",fill="Direction \nof DEG")+coord_flip()+theme(panel.background = element_blank() ,axis.text.x = element_text(size = 12),axis.text.y = element_text(size=12),legend.position = "none", axis.title.y = element_blank())+ scale_fill_manual(values = c("firebrick","steelblue")) +lims(y=c(-250,250)) + geom_hline(yintercept=0)
neurogenic_p = ggplot(m_neurogenic,aes(x=cluster_name,y=value,fill=name))+ geom_bar(stat = "identity",alpha=.8,col="black")+labs(x="",y="",fill="Direction \nof DEG")+coord_flip()+theme(panel.background = element_blank() ,axis.text.x = element_text(size = 12),axis.text.y = element_text(size=12),axis.title.y = element_blank(),legend.position = "none")+ scale_fill_manual(values = c("firebrick","steelblue")) +lims(y=c(-250,250)) + geom_hline(yintercept=0)
glia_p = ggplot(m_glia,aes(x=cluster_name,y=value,fill=name))+ geom_bar(stat = "identity",alpha=.8,col="black")+labs(x="",y="",fill="Direction \nof DEG")+coord_flip()+theme(panel.background = element_blank() ,axis.text.x = element_text(size = 12),axis.text.y = element_text(size=12),axis.title.y = element_blank(),legend.position = "none")+ scale_fill_manual(values = c("firebrick","steelblue")) +lims(y=c(-250,250)) + geom_hline(yintercept=0)
other_p = ggplot(m_other,aes(x=cluster_name,y=value,fill=name))+ geom_bar(stat = "identity",alpha=.8,col="black")+labs(x="",y="",fill="Direction \nof DEG")+coord_flip()+theme(panel.background = element_blank() ,axis.text.x = element_text(size = 12),axis.text.y = element_text(size=12),axis.title.y = element_blank(),legend.position = "none")+ scale_fill_manual(values = c("firebrick","steelblue")) +lims(y=c(-250,250)) + geom_hline(yintercept=0)

#make pdf output of all graphs
pdf(paste0(output,"DEG bar plots.pdf"),height = 10,width = 5)
(gc_p+ gaba_p+ mossy_p+ neurogenic_p+ glia_p+ other_p)  +plot_layout(heights = c(3.2,1.5,1,1,1.5,.5),ncol=1)
dev.off()




################# DE Heatmap for top 5 genes and bottom 5 genes average Diff ###########################
#The following section was used to create heatmpas for Figure 6C and Sup fig 10 it includes all DEG heatmaps

########## Function to Create Heatmap ###########
#heatdf = genes selected only for certain cell types or the entire data set
#class = subsetting of data based on either "cluster_name","cluster_number", or "cell_type"
#number_of_genes_ = number of genes to select for heatmap usually set to 5
#title = string text for title of Heatplot
#Our function has dendograms off
#heatmap height and width is in CM 
make_heatmap = function(heatdf,class,number_of_genes, title, values, heatmap_height,heatmap_width){
  subsets <- c("gene",eval(class),"avg_diff")
  col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  col_fun(seq(-3, 3))
  
  df<- bind_rows(
    heatdf %>% 
      group_by(.data[[class]]) %>%
      filter(avg_diff >0) %>%
      top_n(number_of_genes,-p_val),
    
    heatdf %>% 
      group_by(.data[[class]]) %>%
      filter(avg_diff < 0) %>%
      top_n(number_of_genes,-p_val)
  )
  
  if (values == "up"){
    df = df %>% filter(avg_diff >0)
    col_fun = circlize::colorRamp2(c(0, 2, 4), c("white", "red1", "red4"))
    col_fun(seq(-3, 3))
  } else if(values == "down"){
    df = df %>% filter(avg_diff <0)
    col_fun = circlize::colorRamp2(c(-4, -2, 0), c("blue4", "blue1", "white"))
    col_fun(seq(-3, 3))
  }else if (values == "both"){
    df = df
    col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    col_fun(seq(-3, 3))
  }
  
  # Convert selected genes to matrix for heatmap
  mat_ud = df %>% ungroup() %>%
    dplyr::select(all_of(subsets)) %>%
    pivot_wider(names_from = all_of(class), values_from = avg_diff, values_fill = 0, values_fn = mean) %>% 
    tibble::column_to_rownames(var="gene") %>% as.matrix() 
  #column_order_heat = df %>% arrange(cell_type) %>% pull(cluster_name) %>% unique()
  # row_order_heat = df %>% arrange(p_val) %>% pull(gene) %>% unique()
  column_order_heat = c("DG-Trpc6 GC", "DG-Kcnip4 GC", "DG-Kcnc2-GABA", "DG-SST-GABA", "DG-Prrx1-Mossy","DG-Sv2b-Mossy","Astrocyte","Oligo","OPC/VLMC","Microglia/Macrophage","Endothelial")
  
  #h2 = HeatmapAnnotation()

  heatmap_plot <- Heatmap(mat_ud, name="Average\nDifference", col = col_fun, column_title = title,height = unit(heatmap_height,"cm"),width = unit(heatmap_width,"cm"),
                          #row_split = 9,
                          border = T,
                          rect_gp = gpar(col= "black"),
                          column_order = column_order_heat,
                          #row_order = row_order_heat,
                          row_title_side = "left",
                          show_row_dend = F,
                          #show_row_names = F,
                          show_column_dend = F,
                          #column_dend_reorder = F,
                          column_names_side = 'bottom',
                          #show_row_names = F,
                          row_names_side = "left",
                          row_names_gp = gpar(fontsize=8),
                          heatmap_legend_param = list(direction="horizontal",title_position = "topcenter"))
  
  return(draw(heatmap_plot, heatmap_legend_side = "top", annotation_legend_side = "top"))
  
}

# Selected highest responders from our bar plots
selected_cell_names = c("DG-Trpc6 GC", "DG-Kcnip4 GC", "DG-Kcnc2-GABA", "DG-SST-GABA", "DG-Prrx1-Mossy","DG-Sv2b-Mossy","Astrocyte","Oligo","OPC/VLMC","Microglia/Macrophage","Endothelial")

selected_cells = heatmapdf %>% filter(cluster_name%in%selected_cell_names) 

make_heatmap(heatdf = selected_cells,class = "cluster_name",number_of_genes = 5,title = "",values = "up",heatmap_height = ,heatmap_width = )

heat_down = make_heatmap(heatdf = selected_cells, class = "cluster_name",number_of_genes = 5,title = "",values = "down")

selected_cell_names= c("DG-NSC")


## Figure 6 C heatmaps didn't use combined
pdf(paste0(output,"DEGs in clusters of interest_both.pdf"),height = 6,width = 15 )
heat_both =  make_heatmap(heatdf = selected_cells,class = "cluster_name",number_of_genes = 15,title = "DEGs in NSCs",values = "both",heatmap_height =.5 ,heatmap_width =15 )
dev.off()

## Figure 6 C heatmaps 
pdf(paste0(output,"DEGs in clusters of interest_up.pdf"),height = 10,width = 4)
make_heatmap(heatdf = selected_cells,class = "cluster_name",number_of_genes = 5,title = "",values = "up",heatmap_height =15 ,heatmap_width =4.2)
dev.off()

## Figure 6 C heatmaps 
pdf(paste0(output,"DEGs in clusters of interest_down.pdf"),height = 10,width = 4)
make_heatmap(heatdf = selected_cells,class = "cluster_name",number_of_genes = 5,title = "",values = "down",heatmap_height =15 ,heatmap_width =4.2)
dev.off()


### Supplemental heat maps for all genes Sup figure 10
pdf(paste0(output,"GABA Heatmap.pdf"),height =10 ,width =4)
make_heatmap(gaba_genes,"cluster_name",3,"",values = "both",heatmap_height = 15,heatmap_width = 2.75)
dev.off()

pdf(paste0(output,"Mossy Heatmap.pdf"),height =10 ,width =4)
make_heatmap(mossy_genes,"cluster_name",5,"",values = 'both',heatmap_height = 15, heatmap_width = 2.5)
dev.off()

pdf(paste0(output,"GC Heatmap.pdf"),height =12 ,width =4)
make_heatmap(gc_genes,"cluster_name",3,"",values = "both" ,heatmap_height = 20,heatmap_width = 5.4)
dev.off()

pdf(paste0(output,"Neurogenic Heatmap.pdf"),height =10 ,width =4)
make_heatmap(neurogenic_genes,"cluster_name",3,"",values="both",heatmap_height = 15,heatmap_width = 2.5)
dev.off()

pdf(paste0(output,"Glia Heatmap.pdf"),height =10 ,width =4)
make_heatmap(glia_genes,"cluster_name",3,"",values='both',heatmap_height = 15,heatmap_width = 2.5)
dev.off()

pdf(paste0(output,"Other Heatmap.pdf"),height =10 ,width =4)
make_heatmap(other_genes,"cluster_name",3,"",values = 'both',heatmap_height =5 ,heatmap_width =.8)
dev.off()

#potential figure for both
pdf(file = paste0(output,"Guidance_cues_heatmap.pdf"),width = 6, height=7)
guidance_cues_heatmap = heatmapdf %>% filter(grepl("Tenm|Cntn|Sema|Nrx|Flr|Eph|Cdh|Slit|Unc",x = gene)) %>% make_heatmap(class = "cluster_name",number_of_genes = 5,title = "Development-Like Cues",values = "both")
dev.off()



## Goterm analysis for seperatedd upregulated and downregulated DEGs 
got_df_up = selected_cell_df %>% filter(avg_diff >0) %>%  pull(gene) %>% gost( 
  organism ="mmusculus", ordered_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
  measure_underrepresentation = FALSE, evcodes = TRUE, 
  user_threshold = 0.05, correction_method = "g_SCS", 
  custom_bg = NULL, numeric_ns = "", sources = "GO:BP", as_short_link = FALSE)

got_df_down = selected_cell_df %>% filter(avg_diff < 0) %>% pull(gene) %>% gost( 
  organism ="mmusculus", ordered_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
  measure_underrepresentation = FALSE, evcodes = TRUE, 
  user_threshold = 0.05, correction_method = "g_SCS", 
  custom_bg = NULL, numeric_ns = "", sources = "GO:BP", as_short_link = FALSE)

got_df_up$result$query = str_replace(got_df_up$result$query,pattern = "query_1",replacement = "up")

got_df_down$result$query = str_replace(got_df_down$result$query,pattern = "query_1",replacement = "down")

got_df_merged = rbind(got_df_up$result,got_df_down$result) %>% tibble()

got_df_merged = add_column(.data= got_df_merged, "log10"= -log10(got_df_merged$p_value),.after = "p_value")


### Used for goeterm analysis of entire data set New go term analysis perfomred 3-14-22 using Gprofiler 2. This go term analysis was used for Figure 6D
setwd("/nas/longleaf/home/lquin003/Seurat/goterms/")
dedirectory2 <- setwd("/nas/longleaf/home/lquin003/Seurat/goterms/")
files_list <- list.files(dedirectory2)

Granule_cells <- c(0,1,2,3,4,5,7,8,10,11,14,26,30,32)
GABA <- c(18,19,20,23,40)
Mossy <- c(17,21,24)
Subiculum <- c(9,33,38,39,41)
CA_Excitatory <- c(6,12,13,22,28,34)
Neurogenic <- c(15,25,16)
Glia <- c(27,31,35,36)
Other_notdg <- c(37,42)
Other<- c(29)
RGL <-(15)

#take only files that have DEresults
goonlyfiles <- str_sort(files_list[grep("03_14_22",x = files_list,ignore.case = T)],numeric=T)

# forloop to read in all data from all Gpfofiler files 
goonlydf <- data.frame()
for (i in seq_along(goonlyfiles)) {
  goonlytemp <- read.delim(goonlyfiles[i],header = T,row.names = NULL)
  goonlydf <- rbind(goonlydf,goonlytemp)
}


#Filtering of GO Terms
#Add column where p value is converted to -log10 for graphing
#Convert to tibble
#add column including cluster names after query
#remove GO:CC or cell components from the go terms list
#remove CA Excitatory Cluster Go terms since they are considered contamination


goonlydf <- add_column(goonlydf, "log10"= -log10(goonlydf$p_value),.after = "p_value") %>% 
  as_tibble() %>%
  add_column(.after ="query", "cluster_number" = as.numeric(str_extract(goonlydf$query, "\\d+"))) %>%
  filter(source %in% c("GO:BP") & significant == TRUE) %>%  
  filter(!cluster_number %in% c(CA_Excitatory,Subiculum))

# name each cluster based on their cluster number
goonlydf <- goonlydf %>% mutate("cluster_name" = case_when(
  cluster_number %in% 0 ~ "0 DG-Kcnip4 GC", cluster_number %in% 1 ~ "1 DG-Stxbp6 GC", cluster_number %in% 2 ~ "2 DG-Mast4 GC", cluster_number %in% 3 ~ "3 DG-Hectd2 GC", cluster_number %in% 4 ~ "4 DG-Ntng1 GC",
  cluster_number %in% 5 ~ "5 DG-Csgalnact1 GC", cluster_number %in% 6 ~ "6 CA3-Excitatory", cluster_number %in% 7 ~ "7 DG-Trpc6 GC", cluster_number %in% 8 ~ "8 DG-Cdh12 GC",
  cluster_number %in% 9 ~ "9 SB-Dpp10-Neuron", cluster_number %in% 10 ~ "10 DG-Map4k GC", cluster_number %in% 11 ~ "11 DG-Slc4a4 GC", cluster_number %in% 12 ~ "12 CA3-Subtype",
  cluster_number %in% 13 ~ "13 CA1-Pex5l",cluster_number %in% 14 ~ "14 DG-Cdh13-GC", cluster_number %in% 15 ~ "15 DG-NSC", cluster_number %in% 16 ~ "16 DG-Immature GC",
  cluster_number %in% 17 ~ "17 DG-Mossy-Prrx1", cluster_number %in% 18 ~ "18 DG-Kcnc2-GABA", cluster_number %in% 19 ~ "19 DG-Cnr1-GABA", cluster_number %in% 20 ~ "20 DG-Nos-GABA",
  cluster_number %in% 21 ~ "21 DG-Mossy-Sv2b", cluster_number %in% 22 ~ "22 CA1-Unc5d", cluster_number %in% 23 ~ "23 DG-SST-GABA" , cluster_number %in% 24 ~ "24 DG-Mossy-Slit2",
  cluster_number %in% 25 ~ "25 DG-Ipc/Neuroblast", cluster_number %in% 26 ~ "26 DG-Negr1-GC", cluster_number %in% 27 ~ "27 Oligo", cluster_number %in% 28 ~ "28 CA2/CA3-Excitatory-1",
  cluster_number %in% 29 ~ "29 Endothelial", cluster_number %in% 30 ~ "30 DG-Sema6d-GC", cluster_number %in% 31 ~ "31 OPC/VLMC", cluster_number %in% 32 ~ "32 DG-Arc-GC",
  cluster_number %in% 33 ~ "33 SB-Foxp2-GABA", cluster_number %in% 34 ~ "34 CA2/CA3 Excitatory-2", cluster_number %in% 35 ~ "35 Microglia/Macrophage", cluster_number %in% 36 ~ "36 Astrocyte",
  cluster_number %in% 37 ~ "37 Cajal-Retzius-GABA", cluster_number %in% 38 ~ "38 SB-Tox-GABA", cluster_number %in% 39 ~ "39 SB-Excitatory-Neuron", cluster_number %in% 40 ~ "40 DG-PV-GABA",
  cluster_number %in% 41 ~ "41 SB-Meis2-GABA", cluster_number %in% 42 ~ "42 Ependymal"  
))


#print go only df to check data table
goonlydf

# name each cluster based on their cluster number
write.csv(goonlydf,paste0(output,"DE GO Terms Method2_Filtered for gobp only new DEGs.csv"))


#subset smaller DF for each cell type cluster class
select_go <- function(cell_type,data_frame){
  selection <- filter(data_frame, cluster_number %in% cell_type)
  return(selection)
}

#subset each data set and appened cell_type which will be used for graphing later
gc_goterms <- select_go(Granule_cells,goonlydf) %>% mutate("cell_type"="Granule Cell")
mossy_goterms <- select_go(Mossy,goonlydf) %>% mutate("cell_type"="Mossy")
neurogenic_goterms <- select_go(Neurogenic,goonlydf) %>% mutate("cell_type"="Neurogenic")
gaba_goterms <- select_go(GABA,goonlydf) %>% mutate("cell_type"="GABA")
glia_goterms <- select_go(Glia,goonlydf) %>% mutate("cell_type"="Glia")
other_goterms <- select_go(Other,goonlydf) %>% mutate("cell_type"="Other")

wrs2go_data <- rbind(gc_goterms,gaba_goterms,mossy_goterms,neurogenic_goterms,glia_goterms,other_goterms) %>% as_tibble()

### Go term plots ###
wrs2go_data = wrs2go_data %>% filter(!grepl("urinary",x = term_name))

selected_cell_numbers = c(7,0,18,23,17,21)
selected_cell_names = c("7 DG-Trpc6 GC", "0 DG-Kcnip4 GC", "18 DG-Kcnc2-GABA", "23 DG-SST-GABA", "17 DG-Mossy-Prrx1","21 DG-Mossy-Sv2b")

selected_glia_clusters = c(Glia,Other)
wrs2go_data %>% filter(selected_glia_clusters %in% cluster_number)


#wrs2go_data$cell_type <- factor(wrs2go_data$cell_type, levels = c("Granule Cell","GABA","Mossy","Neurogenic","Glia","Other"))

wrs2go_data = wrs2go_data %>% filter(cluster_number %in% selected_glia_clusters)

wrs2go_data = wrs2go_data %>% filter(cluster_number %in% c(Mossy,Granule_cells,GABA), !cluster_number %in% selected_cell_numbers)

#wrs2go_data$cluster_name = factor(wrs2go_data$cluster_name, levels =selected_cell_names)

write.csv(wrs2go_data,paste0(output,"Glia_cell_GO_terms.csv"))

### GO Term Function Plot Description of Terms ###
#go_df = go term data set wrs2go_data
#Up or down = string of up or down go terms
#column_to_group = column used to select data from df either cell type cluster number or cell name
#number_of_genes = number of genes to select
#axis_text_size = numeric value for axis text size
#nrow = number of rows for graph usually 1 or 2
#to_filter = cell types that wish to be filtered, enter empty string or none as a string if you don't want to filter cell types out



make_GO_Plot = function(go_df, up_or_down, column_to_group, number_of_genes ,bar_color,sig_line_col, axis_text_size, numrow, to_filter){
  
  #Theme for Go term plot_ edit if you need to change parameters
  my_theme <- theme(strip.background = element_rect(fill=NA),
                    strip.text = element_text(size = 8),
                    strip.placement = "inside",
                    axis.text.y = element_text(size=axis_text_size), 
                    panel.background = element_blank())
  
  #Plot by Cluster
  plot_df <- go_df %>%
    filter(grepl(up_or_down, query)) %>%
    #filter(!cluster_number %in% to_filter) %>%
    group_by(cluster_name) %>%
    top_n(number_of_genes, log10) %>%
    arrange(log10,.by_group = T) %>%
    ungroup()%>%
    distinct(term_name,.keep_all = T)
  
  # Create Bar graph with plot data
  go_plot <- plot_df %>%
    mutate(names = paste(term_id,term_name)) %>%
    mutate(cluster_name = stringr::str_replace(cluster_name, "[:digit:][:digit:]*[:space:]*", "")) %>% 
    ggplot(aes(x=fct_reorder(names, log10,.fun = mean), y=log10))+
    #ggplot(aes(x=names, y=log10))+
    geom_bar(stat="identity",color="black",fill=bar_color,width = .8)+
    coord_flip()+
    geom_hline(yintercept = -log10(.05),linetype= "dashed",col=sig_line_col)+
    my_theme+
    ylim(c(0,14))+
    labs(title ="",x="",y="-log10(FDR)") +
    geom_text(aes(y=10,label=cluster_name),size=6,hjust=-.3)
    #facet_wrap(~.data[[column_to_group]], nrow = numrow)+
    #annotate("segment", x=-Inf,xend = Inf, y=-Inf, yend = -Inf)+
    #annotate("segment",x=-Inf,xend = -Inf,y=-Inf,yend=Inf)#+ 
  #geom_text(data = go_df, mapping= aes(x=fct_reorder(term_name, log10), y=log10, label=cluster_name))
  go_plot
  #write.csv(plot_df,paste0(output,"mainfiggoplot1.csv"))
}

#make_GO_Plot = function(go_df, up_or_down, column_to_group, number_of_genes ,bar_color, y_axis_text_size)
p1 <- make_GO_Plot(go_df = wrs2go_data, up_or_down =  "up", column_to_group = "cluster_name", number_of_genes = 10, bar_color = "pink", sig_line_col = "black", axis_text_size = 16, numrow =  1, to_filter = NA)
p2 <- make_GO_Plot(go_df = wrs2go_data, up_or_down =  "down", column_to_group = "cluster_name", number_of_genes = 10, bar_color = "lightblue",sig_line_col = "black", axis_text_size = 16, numrow =  1, to_filter = NA)
p1

pdf(file =paste0(output,"Figure 6 F.pdf"),width = 12,height = 12)
p1/p2
dev.off()





####### Old figure 6 D
up_or_down = "up"
number_of_genes = 1
go_df = wrs2go_data


#main figure 6D
my_theme <- theme(strip.background = element_rect(fill=NA),
                  strip.text = element_text(size = 16),
                  strip.placement = "inside",
                  axis.text.x = element_text(size = 16),
                  axis.text.y = element_text(size=20), 
                  panel.background = element_blank())

plot_df <- go_df %>%
  filter(grepl(up_or_down, query)) %>%
  #filter(!cluster_number %in% to_filter) %>%
  group_by(cluster_name) %>%
  top_n(number_of_genes, log10) %>%
  arrange(log10,.by_group = T) %>%
  ungroup()%>%
  distinct(term_name,.keep_all = T)

# Create Bar graph with plot data
go_plot1 <- plot_df %>%
  mutate(names = paste(term_id,term_name)) %>% 
  mutate(cluster_name = stringr::str_replace(cluster_name, "[:digit:][:digit:]*[:space:]*", "")) %>% 
  ggplot(aes(x=fct_reorder(names, log10,.fun = mean), y=log10,fill=cluster_name))+
  #ggplot(aes(x=names, y=log10))+
  geom_bar(stat="identity",color="black", fill="pink",width = .8)+
  coord_flip()+
  #geom_hline(yintercept = -log10(.05),linetype= "dashed",col=sig_line_col)+
  my_theme+
  labs(title ="",x="",y="-log10(FDR)")+
  geom_text(aes(y=-Inf,label=cluster_name),size=6,hjust=-.3)

go_plot1


#Plot by Cluster
plot_df2 <- go_df %>%
  filter(grepl("down", query)) %>%
  #filter(!cluster_number %in% to_filter) %>%
  group_by(cluster_name) %>%
  top_n(number_of_genes, log10) %>%
  arrange(log10,.by_group = T) %>%
  ungroup()%>%
  distinct(term_name,.keep_all = T)

# Create Bar graph with plot data
go_plot2 <- plot_df2 %>%
  mutate(names = paste(term_id,term_name)) %>%
  mutate(cluster_name = stringr::str_replace(plot_df2$cluster_name, "[:digit:][:digit:]*[:space:]*", "")) %>% 
  ggplot(aes(x=fct_reorder(names, log10,.fun = mean), y=log10,fill=cluster_name))+
  #ggplot(aes(x=names, y=log10))+
  geom_bar(stat="identity",color="black", fill="lightblue",width = .8)+
  coord_flip()+
  #geom_hline(yintercept = -log10(.05),linetype= "dashed",col=sig_line_col)+
  my_theme+
  labs(title ="",x="",y="-log10(FDR)")+
  geom_text(aes(y=-Inf,label=cluster_name),size=6,hjust=-.3)
go_plot2


pdf(paste0(output,"Figure 6D Go Plot.pdf"),width =16,height = 7)
go_plot1+go_plot2 + plot_layout(ncol=1)
dev.off()


### Figure 7 D Go terms for NSC 

# Create Bar graph with plot data
go_plot1 <- wrs2go_data %>% 
  filter(cluster_number == c(15))%>% 
  mutate(names = paste(term_id,term_name)) %>% 
  mutate(cluster_name = stringr::str_replace(cluster_name, "[:digit:][:digit:]*[:space:]*", "")) %>% 
  ggplot(aes(x=fct_reorder(names, log10,.fun = mean), y=log10,fill=cluster_name))+
  #ggplot(aes(x=names, y=log10))+
  coord_flip()+
  geom_bar(stat="identity",color="black", fill="pink",width = .8)+
  #geom_hline(yintercept = -log10(.05),linetype= "dashed",col=sig_line_col)+
  my_theme+
  labs(title ="",x="",y="-log10(FDR)")+
  geom_text(aes(y=-Inf,label=cluster_name),size=6,hjust=-.3)


pdf(paste0(output,"Figure 6D Go Plot.pdf"),width =8,height = 2)
go_plot1
dev.off()


##### Correlation between Control and Trt 


##### Correlation Scatter Plots ##### Testing to see if correlation with waterfall data exists

dg.cpm_avg_expression_RNA_ctrl = dg.cpm %>% subset(trt == "yfp") %>%  AverageExpression() #*1e6 #might not need to multiply by 1million
dg.cpm_avg_expression_RNA_trt = dg.cpm %>% subset(trt == "chr2") %>%  AverageExpression() #*1e6 #might not need to multiply by 1million

df_ctrl = dg.cpm_avg_expression_RNA_ctrl$SCT %>% as.data.frame() %>% rownames_to_column(var="gene") %>% pivot_longer(cols = !gene)
df_trt =  dg.cpm_avg_expression_RNA_trt$SCT %>% as.data.frame() %>% rownames_to_column(var="gene") %>% pivot_longer(cols = !gene)

df_merged = inner_join(df_ctrl,df_trt,"gene") %>% rename("control" = "value.x", "treatment"= "value.y") 

wf_nsc_merge = inner_join(wf_rowmeans_df,nsc_df) %>% as.tibble() %>% filter(!nsc1==0,!nsc2==0,!early==0,!late==0)

ctrl_vs_trt = round(cor(log10(df_merged$control+1),log10(df_merged$treatment+1) ,method = "pearson"),4)

pairwise_corplot = ggplot(df_merged, aes(x=log10(control+1),y=log10(treatment+1))) +
  geom_point(size=1,alpha=.5,color="black")+
  labs(x="Control", y="Treatment",title = "Control_vs_Treatment")+
  theme(panel.background = element_blank(),axis.title = element_text(size = 15),panel.border = element_rect(fill = "transparent",color = "black",))+
  geom_smooth(method = "lm",se = F,color="red")#+
  #annotate("text", x=1, y=2,color="red" ,label =paste("Correlation =",ctrl_vs_trt))

pairwise_corplot



###Figure 6

all_markers <- FindMarkers(object = dg.cpm, assay = "RNA", group.by = "trt",
                           slot = "scale.data", 
                           test.use = "wilcox", 
                           ident.1 = "chr2", 
                           ident.2 = "yfp",
                           subset.ident = c(Granule_cells,Glia,Mossy,Other,GABA),
                           logfc.threshold = .1)

all_markers = all_markers %>%  rownames_to_column(var="gene") %>%  filter(!gene=="Xist")


#### Function to create Volcano Plot, Default number of top Genes is 10,####
make_volc_p = function(cluster_markers,number_of_genes,ptitle){
  
  volc_df <- cluster_markers %>% 
    #rownames_to_column("gene") %>% 
    as_tibble()%>%
    mutate(Expression = case_when(avg_diff >= log(1.00) & p_val <= .01 ~ "Up_regulated",
                                  avg_diff <= -log2(1.00) & p_val <= .01 ~ "Down_regulated",
                                  TRUE~ "Unchanged"))
  
  top_genes <-  bind_rows(
    volc_df %>% 
      filter(Expression == 'Up_regulated') %>% 
      arrange(p_val_adj, desc(abs(avg_diff))) %>% 
      head(number_of_genes),
    volc_df %>% 
      filter(Expression == 'Down_regulated') %>% 
      arrange(p_val_adj, desc(abs(avg_diff))) %>% 
      head(number_of_genes)
  )      
  
  mytheme <- theme(panel.background = element_blank(),
                   legend.position = "none",
                  text = element_text(size = 30),
                   axis.line = element_line(size=.5),
                   axis.text = element_text(size=12),
                   axis.ticks = element_line(size=.5),
                   axis.title = element_text(size=18),
                  title = element_text(size=14))
  
  volc_p <- ggplot(volc_df, aes(x=avg_diff, y=-log(p_val,10)))+
    geom_point(aes(color=Expression),inherit.aes = T,stat = "identity", size = 1.5)+
    xlab("Average Log Fold Change")+
    ylab(expression("-log"[10]*"(P Value)"))+
    labs(title = ptitle )+
    mytheme+
    geom_vline(xintercept = 0, linetype="longdash")+
    geom_vline(xintercept=log(1.1), linetype="longdash",color="red")+
    geom_vline(xintercept = -log(1.1),linetype="longdash",color="blue")+
    scale_color_manual(values = c("blue2","black","red2"))+
    geom_text_repel(data=top_genes, mapping= aes(x=avg_diff,y=-log(p_val,10), label=gene,size=5))
  volc_p
  return(list(volc_p,volc_df))
}


all_trt_vs_ctrl_volcp_data = make_volc_p(cluster_markers = all_markers,number_of_genes = 15, ptitle = "Ctrl Vs Trt")
all_trt_vs_ctrl_volcp_plot = all_trt_vs_ctrl_volcp_data[[1]]
all_trt_vs_ctrl_volcp_df = all_trt_vs_ctrl_volcp_data[[2]]

pdf(paste0(output,"Main_Figure_Ctrl_vs_trt_volcp.pdf"),height = 9,width = 9)
all_trt_vs_ctrl_volcp_plot
dev.off()



export_sct =  dg.cpm@assays$SCT@data


write.table(export_sct,file = "dg.cpm_SCT_gene_cell_matrix.txt", quote = F,sep = "\t",col.names = NA)

export_counts = dg.cpm@assays$RNA@counts

write.table(export_counts,file = "DG_raw_counts_gene_cell_matrix.txt", quote = F,sep = "\t",col.names = NA)


write.table(export_idents,file="cluster_idents.txt",quote = F, sep = "\t",col.names = NA)
