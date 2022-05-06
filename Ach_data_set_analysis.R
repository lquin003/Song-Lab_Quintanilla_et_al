#uses Seurat Version 3.2.2
assign(".lib.loc", "/nas/longleaf/home/lquin003/R/alternate_libs/", envir = environment(.libPaths))

#r_library_location = "/nas/longleaf/home/lquin003/R/alternate_libs/"
output <- "/nas/longleaf/home/lquin003/Documents/Output/"

#Libray location 
#library(dplyr,lib.loc = r_library_location)
#library(Seurat,lib.loc =r_library_location)
#library(patchwork,lib.loc = r_library_location)
#library(SeuratObject,lib.loc =r_library_location)
#library(readr,lib.loc = r_library_location)
#library(devtools,lib.loc = r_library_location)
#library(presto,lib.loc = r_library_location)

library(dplyr)
library(Seurat)
library(patchwork)
library(readr)
library(devtools)
library(presto)
library(slingshot)

#set output for graphs and other files use this paste0 output before each file name to use this
output <- "/nas/longleaf/home/lquin003/Documents/Output/"
cluster_names <- read.csv("/nas/longleaf/home/lquin003/Seurat/Cluster_Annotations_R2.txt") %>% .[1:43,] %>% as.character()

###### Load the Ach dataset #####
ach <- readRDS(file = "/proj/jsonglab/projects/Luis/single_cell_analysis/Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_.rds")

#rename clusters with cluster names
names(cluster_names) <- levels(ach)
ach <- RenameIdents(ach, cluster_names)
  
###Add metadata to Ach including all cell types
###Add metadata to Ach including all cell types
##Create DF that contains seurat clusters in one column and cell type in the other
assign_cell_location <- function(seurat_object){
  #Cluster assignment using cell type for all clusters
  Granule_cells <- c(0,1,2,3,4,5,7,8,10,11,14,26,30,32)
  GABA <- c(18,19,20,23,40)
  Mossy <- c(17,21,24)
  Subiculum <- c(9,33,38,39,41)
  CA_Excitatory <- c(6,12,13,22,28,34)
  Neurogenic <- c(15,25,16)
  Glia <- c(27,31,35,36)
  Other<- c(29,37,42)
  
  # Create Data frame with cell names and cluster numbers in first col
  df <- data.frame("seurat_clusters" = seurat_object$seurat_clusters, "cell_location"= NA)
  # assign cell types to each cluster according to specified identities above.
  df<-  mutate(df, "cell_location" = case_when(df$seurat_clusters %in% c(Granule_cells,GABA,Neurogenic,Mossy) ~ "Dentate Gyrus",
                                           df$seurat_clusters %in% Subiculum ~ "Subiculum",
                                           df$seurat_clusters %in% c(Glia,Other) ~ "Hippocampus",
                                           df$seurat_clusters %in% CA_Excitatory ~"CA Regions"))
  
  seurat_object <- AddMetaData(seurat_object,df)
  return(seurat_object)
}

ach = assign_cell_location(ach)

assign_cell_type <- function(seurat_object){
  #Cluster assignment using cell type for all clusters
  Granule_cells <- c(0,1,2,3,4,5,7,8,10,11,14,26,30,32)
  GABA <- c(18,19,20,23,40)
  Mossy <- c(17,21,24)
  Subiculum = c(9,33,38,39,41)
  CA_Excitatory <- c(6,12,13,22,28,34)
  Neurogenic <- c(15,16,25)
  Glia <- c(27,31,35,36)
  Other<- c(29,37,42)
# Create Data frame with cell names and cluster numbers in first col
df <- data.frame("seurat_clusters" = seurat_object$seurat_clusters, "cell_type"= NA)
# assign cell types to each cluster according to specified identities above.
df<-  mutate(df, "cell_type" = case_when(df$seurat_clusters %in% Granule_cells ~ "Granule Cell",
                                           df$seurat_clusters %in% Mossy ~ "Mossy",
                                           df$seurat_clusters %in% GABA ~ "GABA",
                                           df$seurat_clusters %in% CA_Excitatory~"CA Excitatory",
                                            df$seurat_clusters %in% Subiculum~"Subiculum",
                                           df$seurat_clusters %in% Neurogenic ~ "Neurogenic",
                                           df$seurat_clusters %in% Glia ~ "Glia",
                                           df$seurat_clusters %in% Other~ "Other"
                                           ))
seurat_object <- AddMetaData(seurat_object,df)
return(seurat_object)
  }

ach = assign_cell_type(ach)

assign_cluster_names = function(seurat_object){
  df  = data.frame("seurat_clusters" = seurat_object$seurat_clusters, "cluster_name"= NA)
  df = mutate(df, "cluster_name" = case_when(
    df$seurat_clusters %in% 0 ~ "0 DG-Kcnip4 GC", df$seurat_clusters %in% 1 ~ "1 DG-Stxbp6 GC", df$seurat_clusters %in% 2 ~ "2 DG-Mast4 GC", df$seurat_clusters %in% 3 ~ "3 DG-Hectd2 GC", df$seurat_clusters %in% 4 ~ "4 DG-Ntng1 GC",
    df$seurat_clusters %in% 5 ~ "5 DG-Csgalnact1 GC", df$seurat_clusters %in% 6 ~ "6 CA3-Excitatory", df$seurat_clusters %in% 7 ~ "7 DG-Trpc6 GC", df$seurat_clusters %in% 8 ~ "8 DG-Cdh12 GC",
    df$seurat_clusters %in% 9 ~ "9 SB-Dpp10-Neuron", df$seurat_clusters %in% 10 ~ "10 DG-Map4k GC", df$seurat_clusters %in% 11 ~ "11 DG-Slc4a4 GC", df$seurat_clusters %in% 12 ~ "12 CA3-Subtype",
    df$seurat_clusters %in% 13 ~ "13 CA1-Pex5l",df$seurat_clusters %in% 14 ~ "14 DG-Cdh13 GC", df$seurat_clusters %in% 15 ~ "15 DG-NSC", df$seurat_clusters %in% 16 ~ "16 DG-Immature GC",
    df$seurat_clusters %in% 17 ~ "17 DG-Mossy-Prrx1", df$seurat_clusters %in% 18 ~ "18 DG-Kcnc2-GABA", df$seurat_clusters %in% 19 ~ "19 DG-Cnr1-GABA", df$seurat_clusters %in% 20 ~ "20 DG-Nos-GABA",
    df$seurat_clusters %in% 21 ~ "21 DG-Mossy-Sv2b", df$seurat_clusters %in% 22 ~ "22 CA1-Unc5d", df$seurat_clusters %in% 23 ~ "23 DG-SST-GABA" , df$seurat_clusters %in% 24 ~ "24 DG-Mossy-Slit2",
    df$seurat_clusters %in% 25 ~ "25 DG-IPC/Neuroblast", df$seurat_clusters %in% 26 ~ "26 DG-Negr1 GC", df$seurat_clusters %in% 27 ~ "27 Oligo", df$seurat_clusters %in% 28 ~ "28 CA2/CA3-Excitatory-1",
    df$seurat_clusters %in% 29 ~ "29 Endothelial", df$seurat_clusters %in% 30 ~ "30 DG-Sema6d GC", df$seurat_clusters %in% 31 ~ "31 OPC/VLMC", df$seurat_clusters %in% 32 ~ "32 DG-Arc GC",
    df$seurat_clusters %in% 33 ~ "33 SB-Foxp2-GABA", df$seurat_clusters %in% 34 ~ "34 CA2/CA3 Excitatory-2", df$seurat_clusters %in% 35 ~ "35 Microglia/Macrophage", df$seurat_clusters %in% 36 ~ "36 Astrocyte",
    df$seurat_clusters %in% 37 ~ "37 Cajal-Retzius-GABA", df$seurat_clusters %in% 38 ~ "38 SB-Tox-GABA", df$seurat_clusters %in% 39 ~ "39 SB-Excitatory-Neuron", df$seurat_clusters %in% 40 ~ "40 DG-PV-GABA",
    df$seurat_clusters %in% 41 ~ "41 SB-Meis2-GABA", df$seurat_clusters %in% 42 ~ "42 Ependymal"
  ))
  seurat_object <- AddMetaData(seurat_object,df)
  return(seurat_object)
}


ach <- assign_cluster_names(ach)

#verify that meta data was inserted assigning a column called cell type
head(ach@meta.data)

saveRDS(ach,file ="/nas/longleaf/home/lquin003/Seurat/Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_LQ_names.rds")

################# Start code from here after saving initial object with names ###############

# read in new ach file with names attached to seurat object slot
ach <- readRDS("/nas/longleaf/home/lquin003/Seurat/Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_LQ_names.rds")

#rename clusters with cluster names
names(cluster_names) <- levels(ach)
ach <- RenameIdents(ach, cluster_names)

#assign factor levels so YFP comes first in graphs to seurat object useful when making plots
ach$trt <- factor(ach$trt, levels = c("yfp","chr2"))

####
# Qualty control Figures ## using ggplot to show bar graphs of different reads #####
c_ncount_p <- ach %>% 
subset(trt =="yfp") %>%
VlnPlot(assay = "RNA", features = c("nCount_RNA"))+NoLegend()+
  ylim(1000,8000) + theme(axis.title = element_blank(),axis.text.x = element_blank())
c_ncount_p


t_ncount_p <- ach %>% 
  subset(trt =="chr2") %>%
  VlnPlot(assay = "RNA", features = c("nCount_RNA"))+NoLegend()+
  ylim(1000,8000)+ theme(axis.title = element_blank(),plot.title = element_blank())
t_ncount_p

c_ncount_p/t_ncount_p

#subset only YFP data
ach_c <- ach %>% subset(trt=="yfp")
ach_c <- bind_cols(
  as_tibble(ach_c[["seurat_clusters"]]),
  as_tibble(ach_c[["nCount_RNA"]]),
  as_tibble(ach_c[["nFeature_RNA"]]),
  as_tibble(ach_c[["percent.mt"]]),
    "trt"="yfp")
#subset only chr2 data
ach_t <- ach %>% subset(trt=="chr2")
ach_t <- bind_cols(
  as_tibble(ach_t[["seurat_clusters"]]),
  as_tibble(ach_t[["nCount_RNA"]]),
  as_tibble(ach_t[["nFeature_RNA"]]),
  as_tibble(ach_t[["percent.mt"]]),
  "trt"="chr2")


# assign color
colorcount <- length(unique(ach_c$seurat_clusters))
# function for color ramp
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))

bptheme1 <- theme(legend.position = "none",panel.background = element_blank(),
                  axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank())

bptheme <- theme(legend.position = "none",panel.background = element_blank())

bp1_c <- ggplot(ach_c,aes(x=seurat_clusters,y=nCount_RNA),fill=seurat_clusters)+
      geom_boxplot(aes(color=seurat_clusters))+ labs (y ="Number of Counts",title ="Control")+
      ylim(2000,8000) + bptheme1 + scale_colour_manual(values = getPalette(colorcount))

bp2_c <- ggplot(ach_c,aes(x=seurat_clusters,y=nFeature_RNA),fill=seurat_clusters)+
  geom_boxplot(aes(color=seurat_clusters))+ labs(y="Number of Features")+
  ylim(1000,4000)+bptheme1 + scale_colour_manual(values = getPalette(colorcount))

bp3_c <- ggplot(ach_c,aes(x=seurat_clusters,y=percent.mt),fill=seurat_clusters)+
  geom_boxplot(aes(color=seurat_clusters))+ labs(x="Cluster Numbers",y="Percent MT") +
  ylim(0,10) + bptheme + scale_colour_manual(values = getPalette(colorcount))

# cell counts of seurat clusters for export 
ach_c_counts <- as.data.frame(table(ach_c$seurat_clusters))


bp1_c
bp2_c
bp3_c

#Control Plot quality contorl
pdf(file = paste0(output,"QC_c.pdf"),width = 8,height = 8)
grid.newpage()
grid.draw(rbind(ggplotGrob(bp1_c),ggplotGrob(bp2_c),ggplotGrob(bp3_c)))
dev.off()

colorcount <- length(unique(ach_t$seurat_clusters))
# function for color ramp
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))

#treatment plot quality contorl
bp1_t <- ggplot(ach_t,aes(x=seurat_clusters,y=nCount_RNA),fill=seurat_clusters)+
  geom_boxplot(aes(color=seurat_clusters))+ labs(y="Number of Counts",title="Treatment")+
  ylim(2000,8000) + bptheme1 + scale_colour_manual(values = getPalette(colorcount))

bp2_t <- ggplot(ach_t,aes(x=seurat_clusters,y=nFeature_RNA),fill=seurat_clusters)+
  geom_boxplot(aes(color=seurat_clusters))+ labs(y="Number of Features")+
  ylim(1000,4000)+bptheme1 + scale_colour_manual(values = getPalette(colorcount))

bp3_t <- ggplot(ach_t,aes(x=seurat_clusters,y=percent.mt),fill=seurat_clusters)+
  geom_boxplot(aes(color=seurat_clusters))+ labs(x="Cluster Numbers",y="Percent MT")+
  ylim(0,10) + bptheme + scale_colour_manual(values = getPalette(colorcount)) 
bp3_t
#### 
ach_t_counts <- as.data.frame(table(ach_t$seurat_clusters))


cell_number = inner_join(by ="Var1",ach_c_counts,ach_t_counts) %>% rename("Cluster_Number"=Var1, "Control"=Freq.x, "ChR2"=Freq.y) 
cell_number = cell_number[heatmap_row_order,]
p1 = tableGrob(cell_number)


pdf(file = paste0(output,"Table_of_counts.pdf") ,height=12, width =4 )
grid.draw(p1)
dev.off()

pdf(file = paste0(output,"QC_t.pdf"),width =8 ,height = 8)
grid.newpage()
grid.draw(rbind(ggplotGrob(bp1_t),ggplotGrob(bp2_t),ggplotGrob(bp3_t)))
dev.off()


pdf(file = paste0(output,"Quality Control plots.pdf"),width =16 ,height = 8)
(bp1_c/bp2_c/bp3_c)-(bp1_t/bp2_t/bp3_t)
dev.off()


#### Count cells in clusters Bar Graphs for Supplemental Figure 9 ####
cluster_cell_counts_yfp <- ach %>% subset(trt == "yfp") %>% Idents() %>% table()
cluster_cell_counts_chr2 <- ach %>% subset(trt == "chr2") %>% Idents() %>% table()


cluster_cell_counts <- table(Idents(ach)) %>% 
  as_tibble(.name_repair = "unique")%>%
  rename(cluster_name = ...1)%>%
  mutate(cluster_numbers = 0:42)

#Assign data frame cell classes and then calculate percentages of cells in each cluster and by cell type
cluster_cell_counts = cluster_cell_counts %>%
  mutate("cell_type" = case_when(cluster_cell_counts$cluster_numbers %in% Granule_cells ~ "Granule Cell",
                                 cluster_cell_counts$cluster_numbers %in% Mossy ~ "Mossy",
                                 cluster_cell_counts$cluster_numbers %in% GABA ~ "GABA",
                                 cluster_cell_counts$cluster_numbers %in% CA_Excitatory~"CA Cells",
                                 cluster_cell_counts$cluster_numbers %in% Subiculum~"Subiculum",
                                 cluster_cell_counts$cluster_numbers %in% Neurogenic ~ "Neurogenic",
                                 cluster_cell_counts$cluster_numbers %in% Glia ~ "Glia",
                                 cluster_cell_counts$cluster_numbers %in% Other~ "Other")) %>%
  filter(!cell_type %in% "CA Excitatory") %>%
  mutate("percent_cluster"= (n/sum(n)*100))%>%
  group_by(cell_type)%>%
  mutate("cell_type_sum" = sum(n)) %>%
  ungroup()%>%
  mutate("percent_cell_type" = (cell_type_sum/(sum(n)) * 100))

## assign order of cell classes for plotting
cluster_cell_counts$cell_type <- factor(cluster_cell_counts,levels = c("Granule Cell","GABA","Mossy","Neurogenic","Glia","Other"))

## Plot of percent of cells in each cell class
percentage_plot <-
  cluster_cell_counts %>%
  distinct(cell_type,percent_cell_type) %>%
  mutate(cell_type = fct_reorder(cell_type,-percent_cell_type))%>%
  ggplot(aes(x="", y= percent_cell_type, fill=cell_type))+
  geom_bar(stat="identity",width = 1)+
  ylim(0,100)+
  theme(panel.background = element_blank(),
        axis.line = element_line(size = .5))+
  scale_fill_brewer(palette = "Set3")+
  #geom_text(aes(x="",y=percent_cell_type,label=paste0(percent_cell_type)))+
  #theme_bw()+
  labs(y="Percentage of Cells from DG Dissection",fill="Cell Type",x="")

print(cluster_cell_counts$percent_cell_type)

# Export Pdf of plot 
pdf(file = paste0(output,"cell_type_barplot.pdf"),width = 3,height = 5)  
percentage_plot
dev.off()

# assign color
colorcount <- length(unique(cluster_cell_counts$cluster_name))
# function for color ramp
getPalette = colorRampPalette(brewer.pal(12, "Set3"))
# plot of percentage of cluster 
percentage_plot_clust <- cluster_cell_counts %>%
  ungroup() %>%
  #mutate(fct_reorder(cluster_name,percent_cluster))%>%
  ggplot(aes(x=fct_reorder(cluster_name,-percent_cluster),y=percent_cluster,fill=cluster_name))+
  geom_bar(stat="identity")+
  #  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,size = 12),
        legend.position = "none", panel.background = element_blank(),
        axis.line = element_line(size = .5))+
  scale_fill_manual(values = getPalette(colorcount))+
  #scale_color_brewer(getPalette(colocount))+
  labs(y="Percent of DG",x="") 
#geom_text(aes(label=percent_cluster),vjust=1.6)
# export as pdf
pdf(file = paste0(output,"barplot_percentage_clust.pdf"),width = 9,height = 6)
percentage_plot_clust
dev.off()

#determine length of unique clusters
colorcount <- length(unique(ach$seurat_clusters))
# function for color ramp
getPalette = colorRampPalette(brewer.pal(12, "Set3"))

## Make UMAP plot that has cluster labels split by treatment
umap_ach <- DimPlot(ach,label=F) + NoLegend()+
  scale_color_manual(values = getPalette(colorcount))+
  plot_annotation() & theme(axis.line = element_blank(),
                          axis.title = element_blank(),
                          axis.ticks = element_blank(),
                          axis.text = element_blank())

#add labels to the UMAP plot and use repel to remove label overlap
LabelClusters(umap_ach, id="ident",repel = T,size=4)

#visualize Umap Plot
umap_ach

pdf(file = paste0(output,"New_UMAP.pdf"),width = 8,height = 8)
LabelClusters(umap_ach, id="ident",repel = T,size=4)
dev.off()


#visualize GABA feature plot
gaba_feature_plot <- ach%>%
subset(seurat_clusters %in% GABA) %>%
FeaturePlot(features = c("Gad1","Gad2","Dpp10","Slc6a1"),order = T,cols = c("grey","goldenrod1","goldenrod4")) + 
  plot_annotation() & theme(axis.text = element_blank(), axis.ticks = element_blank(),axis.title = element_blank())


pdf(paste0(output,"GABA_Feature_Plots.pdf") ,height = 8, width = 8)
gaba_feature_plot
dev.off()

#### GABA Dot Plot For Supplemental ####
gaba_dot_plot <- ach%>%
  subset(seurat_clusters %in% GABA) %>%
  DotPlot(assay = "SCT",features = c("Gad1","Gad2","Dpp10","Slc6a1","Slc6a13","Slc6a11","Vip","Pvalb","Sst","Nos1","Npy","Cck","Calb2","Reln") ,cols = c("grey","goldenrod1","goldenrod4"),cluster.idents = T)+
  NoLegend()+ theme(axis.text.x = element_text(angle = 90,vjust = .5,hjust = 1),
                    axis.title.y = element_blank(),
                    axis.title.x = element_blank())
  
gaba_dot_plot

pdf(paste0(output,"GABA_dot_Plots.pdf") ,height = 5, width = 5)
gaba_dot_plot
dev.off()

#### GABA vln plot for supplemental figures
make_stacked_violin_plot = function(s_object,assay_to_use,slot_to_use,gene_list,title,remove_x,sort){
  violin_plot = s_object %>%
    VlnPlot(object =, assay = assay_to_use,slot=,sort = sort,features =gene_list, stack = T,flip = T) + NoLegend() + labs(title=title) +
    plot_annotation() & theme(strip.text = element_text(face = "italic"),
                              axis.title.x = element_blank(),
                              axis.title.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.line = element_blank(), axis.text.y = element_blank(),
                              axis.text.x = element_text(angle = 90,vjust = .5,hjust = 1),
                              axis.title.y.right = element_blank(),
                              axis.text.y.right = element_blank())
  
  if(remove_x == T){
    violin_plot + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
  } else{
    violin_plot
  }
}

gaba_genes = c("Gad1","Gad2","Slc6a1","Slc6a11","Slc32a1","Vip","Cnr1","Pvalb","Sst","Nos1","Npy","Cck","Cnr1","Reln","Lamp5")

gaba_vln_plot = ach %>% subset(seurat_clusters %in% GABA) %>% make_stacked_violin_plot(assay_to_use = "SCT",slot_to_use = "data",gene_list = gaba_genes,title = "GABA Genes",remove_x = F,sort = F)

gaba_vln_plot


pdf(paste0(output,"GABA_Vln_Plot.pdf") ,height = 6, width = 4)
gaba_vln_plot
dev.off()

### GABA Feature List
gaba_feature_list <- c("Gad2","Gad1","Slc6a11","Slc6a1")
#### GABA Feature plot 
gaba_fp <- ach %>% subset(seurat_clusters %in% GABA) %>%
  FeaturePlot(feature=gaba_feature_list,order = T, cols=c("grey","darkgoldenrod1","darkgoldenrod4"))+
  labs(x="Genes",y="Clusters") + plot_annotation() & theme(axis.text = element_blank(),axis.ticks = element_blank(),axis.title = element_blank())

pdf(paste0(output,"GABA_Feature_plot.pdf"),height = 6,width = 6)
gaba_fp
dev.off()


### GC Feature List
gc_feature_list <- c("Prox1","Calb1","","")
#### GC Feature plot 
gc_fp <- ach %>% subset(seurat_clusters %in% GABA) %>%
  FeaturePlot(feature=gaba_feature_list,order = T, cols=c("grey","darkgoldenrod1","darkgoldenrod4"))+
  labs(x="Genes",y="Clusters") + plot_annotation() & theme(axis.text = element_blank(),axis.ticks = element_blank(),axis.title = element_blank())

pdf(paste0(output,"GABA_Feature_plot.pdf"),height = 6,width = 6)
gab_fp
dev.off()







#GC Vln plot Subplot for GC's
gc_vln_plot <- ach%>%
  subset(seurat_clusters %in% Granule_cells) %>%
  VlnPlot(object =, assay = "SCT", sort = T,features = c("Prox1","Negr1","Unc5d","Stxbp6","Csgalnact1","Mast4","Hectd2","Slc4a4","Cdh12","Cdh13","Mapk4"),
          #cols = c("red1","red2","red3","red4","red1","red2","blue","blue","blue","blue","blue","blue","blue","blue"),
          stack = T,flip = T) + NoLegend() + plot_annotation() & theme(strip.text = element_text(face = "italic"),
                                         axis.title.x = element_blank(),
                                         axis.title.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.line = element_blank(), axis.text.y = element_blank(),
                                         axis.text.x = element_text(angle = 90,vjust = .5,hjust = 1),
                                         axis.title.y.right = element_blank(),
                                         axis.text.y.right = element_blank()) 

pdf(paste0(output,"gc_Vln_Plot.pdf") ,height = 6, width = 5)
gc_vln_plot
dev.off()

#Mossy Cell Clustrers Vln plot stacked
mossy_feature_vlnplot <- ach %>%
  subset(seurat_clusters %in% Mossy)%>%
  VlnPlot(object =, assay = "SCT",slot="data",sort = T,features = c("Cntn6","Grm8","Prrx1","Slit2","Pex5l","Cpne4","Nrxn1","Sv2b","Pdzd2"), stack = T,flip = T)+
  NoLegend() + 
  plot_annotation() & 
  theme(strip.text = element_text(face = "italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90,vjust = .5,hjust = 1),
        axis.title.y.right = element_blank(),
        axis.text.y.right = element_blank()) 

pdf(paste0(output,"mossy_vlnplot_SCT.pdf") ,height =8 ,width =6)
mossy_feature_vlnplot
dev.off()












## Neuroblast Cluster Identification  ##
ach%>% subset(seurat_clusters == 16) %>%
Seurat::DoHeatmap()



c_16_markers_SCT <- FindMarkers(object = ach, assay = "SCT", slot = "scale.data", test.use = "wilcox", ident.1 = "16 Neuroblast", min.pct = .25)

c_15_markers_SCT <- FindMarkers(object = ach, assay = "SCT", slot = "scale.data", test.use = "wilcox", ident.1 = "15 RGL 1", min.pct = .25)


### CA1/Mossy/CA2/3 feature plots
ca_feature_plot <- DotPlot(ach,features = c("Cntn6","Dkkl1","Fgf1","Ephb6","Drd2", "Gal", "Glp1r", "Grm8", "Nmb","Prrx1","Rac3","Vat1","Stxbp6","Igfbp5","Pex5l","Arhgap12","Fibcd1","Cpne4","Nrxn1","Ly6e","Np2yr","Sv2b","Calb2","Pdzd2","Prox1","Gad1","Gad2","Dpp10"),
        assay = "SCT",
        cluster.idents = T,
        cols = c("grey","goldenrod1","goldenrod4"))

ca_feature_plot


########### Dot plot visualization of All genes Dot plot has been changed to Violin Plot ######### Figure 5B
feature_list <- c("Sox2","Gfap","Aqp4","Clu","Hopx","Cd74","Csf1R","Plp1","Gad2","Gad1","Dpp10","Dcx","Prox1","Sv2b",
                  "Flt1","Grik1","Ctss","Cst3","C1qa","Pex5l","Arhgap12","Fibcd1","Cpne4","Ly6e","Npy2R",
                  "MBP", "Plp1","Mobp","Pdgfra","Vcan","Cspg4","Amigo2","Cntn5","Calb2","Calb1",
                  "Foxg1","Nes","Olig","Nrxn1","Reln","Dkk3","Prkcb")

total_feature_list <- c("Gfap","Apoe","Aqp4","Clu","Hopx","Gad2","Gad1","Slc17a7","Calb1","Prox1","Sv2b","Cntn6","Grm8","Cst3","Pex5l","Arhgap12","Mobp","Pdgfra","Snap25","Foxp2","Cacng5","Htr2c","Abi3bp","Ccnd2","Bhlhe22","Gabra5","Prrx1")

total_dp <- ach %>% 
  DotPlot(feature=total_feature_list,cluster.idents = T, assay = "SCT",cols = c("grey","red","blue"))+
  labs(x="Genes",y="Clusters")+
  theme(legend.position = "none",
        axis.line = element_blank())
total_dp

pdf(paste0(output,"total_dotplot_SCT.pdf") ,height =8 ,width =18)
total_dp
dev.off()

#astrocyte_genes = c("Sox9","Etv4","Sall3","Grin2c","Kcng4","Aqp4","Gfap","Hopx","S100b","mki67","Pcna","Sox2","Tbr2","Mcm2","Dcx","Ncam1","Mxd3","Ascl1","Neurog2","Rfp4","Lockd","Eomes","Calb2","Neurod4","Slc17a6")
#Total Features Vln plot stacked

total_feature_vlnplot <- ach %>%
  VlnPlot(object =, assay = "SCT",slot="data",sort = T,features = total_feature_list ,stack = T,flip = T)+ #facet_wrap(~total_feature_list,strip.position = "left") +
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
ach@active.ident = ach$seurat_clusters
Idents(ach) <- factor(Idents(ach), levels= my_levels)

side_total_feature_vlnplot <- ach %>%
  VlnPlot(object =, assay = "SCT",slot="data",features = total_feature_list ,stack = T)+ 
  NoLegend() + scale_y_discrete(limits=rev) + plot_annotation() & theme(strip.text = element_text(face = "italic"),
                                           axis.title = element_blank(),
                                           axis.text.x = element_blank(),
                                           axis.ticks.x = element_blank())

side_total_feature_vlnplot

## Combine with complex heatmap
p1 = grid.grabExpr(draw(cor_heatmap, heatmap_legend_side = "top", annotation_legend_side = "top"))
p2 = side_total_feature_vlnplot

pdf(paste0(output,"Cor_plot_vln_plotmerge.pdf"),height =8,width = 20)
plot_grid(p1,p2,nrow=1, axis="r",rel_heights = c(.9,1),rel_widths = c(1,1))
dev.off()


##
total_dp_ctrl <- ach %>% subset(trt == "yfp") %>%
  DotPlot(feature=feature_list,cluster.idents = T)+
  labs(x="Genes",y="Clusters",title = "YFP")+
  theme(#legend.position = "none",
    axis.line = element_blank())
total_dp_ctrl

total_dp_trt <- ach %>% subset(trt == "chr2") %>%
  DotPlot(feature=feature_list,cluster.idents = T)+
  labs(x="Genes",y="Clusters",title = "chr2")+
  theme(#legend.position = "none",
    axis.line = element_blank())
total_dp_trt

### GABA Feature List
gaba_feature_list <- c("Gad2","Gad1","Slc6a11","Slc")
#### GABA plot 
gab_dp <- ach %>% subset(seurat_clusters == GABA) %>%
DotPlot(feature=gaba_feature_list,cluster.idents = T, cols=c("grey","darkgoldenrod1","darkgoldenrod4"))+
  labs(x="Genes",y="Clusters")+
  theme(legend.position = "none",
        axis.line = element_blank())
gab_dp

### GC Plot
gc_feature_list <- c("Prox1","hectd2","Camk2a","Camk4","Neurod2","Camk2d","Syt1","Grin1","Snap25")
  gc_dp <- ach %>% subset(seurat_clusters == Granule_cells) %>%
  DotPlot(feature=gc_feature_list,cluster.idents = T, cols=c("grey","seagreen1","seagreen4"))+
  labs(x="Genes",y="Clusters")+
  theme(legend.position = "none",
        axis.line = element_blank())
gc_dp

# Mossy Plot
mc_feature_list <- c("Sv2b","Glur2","Glur3","Calb2","Slit2","Galntl6","Htr2a")
mc_dp <- ach %>% subset(seurat_clusters == Mossy) %>%
  DotPlot(feature=mc_feature_list,cluster.idents = T, cols=c("grey","deepskyblue1","deepskyblue4"))+
  labs(x="Genes",y="Clusters")+
  theme(legend.position = "none",
        axis.line = element_blank())
mc_dp

#### Neurogenic Plot
ng_feature_list <- c("Sox2","Gfap","Clu","Vim","Sox11","Cdk1","Tbr2","Hopx","Dcx","Notch1","Calb2","Calb1","Prox1")
ng_dp <- ach %>% subset(seurat_clusters == c(15,25,16,36)) %>%
  DotPlot(feature=ng_feature_list, split.by = "trt" , cluster.idents = T, cols=c("grey","green1","green4"))+
  labs(x="Genes",y="Clusters")+
  theme(legend.position = "none",
        axis.line = element_blank())
ng_dp

#### Glia Plot
gl_feature_list <- c("Ctss","Cst3","C1qa","Gfap")
gl_dp <- ach %>% subset(seurat_clusters == Glia) %>%
  DotPlot(feature=gl_feature_list,cluster.idents = T,cols=c("grey","purple1","purple4"))+
  labs(x="Genes",y="Clusters")+
  theme(legend.position = "none",
        axis.line = element_blank())
gl_dp

#### Other Plot
o_feature_list <- c("Ctss","Cst3","C1qa","Gfap")
o_dp <- ach %>% subset(seurat_clusters == Glia) %>%
  DotPlot(feature=gl_feature_list,cluster.idents = T, cols=c("grey","purple1","purple4"))+
  labs(x="Genes",y="Clusters")+
  theme(legend.position = "none",
        axis.line = element_blank())
o_dp


#### Location Feature Plot ###########
#### Location Feature Plots for proliferation markers ####
FeaturePlot(ach,features = "Dcx",cols=c("grey","green1","green4"), order = T, label = T,split.by = "trt",repel = T)

FeaturePlot(ach,features = "Sox2",cols=c("grey","red1","red4"),order = T, label = T,split.by = "trt", repel = T)

FeaturePlot(ach,features = "Cpeb4",cols=c("grey","blue1","blue4"),order = T,label = T,split.by = "trt", repel = T)

#Location feature plot for manuscript 
total_fplot <- ach %>% 
#  subset(seurat_clusters == c(15,25)) %>%
  FeaturePlot(features =c("Hopx","Gad2","Prox1","Gfap","Pex5l","Sv2b","Mobp","Ctss"),order = T, ncol = 4 ,min.cutoff = 0,cols = c("grey","coral1","coral2","magenta1","magenta2","magenta3"))+
  plot_annotation() &
  theme(axis.text = element_blank(),axis.line = element_blank(),axis.ticks = element_blank(), axis.title = element_blank())
total_fplot

pdf(file=paste0(output,"Featuremap_location Plot.pdf"),width = 12,height = 6)
total_fplot
dev.off()

#### Subset just Neurogenic Clusters and make violin plots
ach$trt <- factor(ach$trt, levels = c("yfp","chr2"))

ng_vp1 <- ach %>% 
  subset(seurat_clusters == c(15,25)) %>%
  VlnPlot(features ="Mt3",split.by = "trt")+
  scale_fill_discrete(type = c("springgreen2","deepskyblue2"))
ng_vp1

ng_vp2 <- ach %>% 
  subset(seurat_clusters == c(15,25)) %>%
  VlnPlot(features ="Hopx",split.by = "trt")+
  scale_fill_discrete(type = c("springgreen2","deepskyblue2"))
ng_vp2

ng_vp3 <- ach %>% 
  subset(seurat_clusters == c(15,25)) %>%
  VlnPlot(features ="Dbi",split.by = "trt")+
  scale_fill_discrete(type = c("springgreen2","deepskyblue2"))
ng_vp3

ng_vp4 <- ach %>% 
  subset(seurat_clusters == c(15,25)) %>%
  VlnPlot(features ="Grin1",split.by = "trt")+
  scale_fill_discrete(type = c("springgreen2","deepskyblue2"))
ng_vp4

ng_vp1+ng_vp2+ng_vp3

######### Identity of RGL cluster DIving or not Jessberger paper
ng_vp1_i <- dg.cpm %>% 
  #subset(seurat_clusters %in% c(15,25,36)) %>%
  VlnPlot(features = c("Mt3","Hopx","Dbi","Gfap") ,assay = "SCT", slot="data",stack = T)#,cols = c("springgreen2","deepskyblue2","darkorchid"))
ng_vp1_i

pdf(file = paste0(output,"neurogenic jessberger 3 genes.pdf"),width = 8,height = 6)
ng_vp1_i
dev.off()

ng_vp_div <- ach %>% 
  subset(seurat_clusters == c(15,25)) %>% 
  VlnPlot(features = c("Cox8a","Thrsp","Apoe","Cst3"),cols =c("springgreen2","deepskyblue2"))+ NoLegend()
ng_vp_div

#### Genes from Jessberger paper main text as violin plots in both of our RGL clusters ####
ng_vp1_all <- ach %>% 
  VlnPlot(features = c("Mt3","Hopx","Dbi"),ncol=1)
ng_vp1_all

pdf(file = paste0(output,"neurogenic jessberger 3 genes all clust.pdf"), width = 12,height = 18)
ng_vp1_all 
dev.off()

##### Main Figures Analysis #####

### load in genes list

# Load in Neurotransmitter Binding genes list
nt_genes <- read_delim("/nas/longleaf/home/lquin003/Seurat/neurotransmitter_binding_genes.txt",col_names = F,delim = "\t") %>%
  arrange(X1) %>% dplyr::select(X1,X3) 
# Load in Ach signaling genes list
ach_genes <- read_delim("/nas/longleaf/home/lquin003/Seurat/acetylcholine_binding_genes.txt",col_names = F,delim = "\t")%>%
  arrange(X1) %>% dplyr::select(X1,X3)
# Load in Calcium Signaling Genes list
ca_sig_genes <- read_delim("/nas/longleaf/home/lquin003/Seurat/calcium_signaling_genes.txt",col_names = F,delim = "\t") %>%  as_tibble %>% 
  dplyr::select(X1) %>% filter(!X1 %in% c("Chrm3","Chrna7","Grin2a","Grin2b","Grm1","Grm5"))

# load learning and mem genes
learn_n_mem_genes = read_delim("/nas/longleaf/home/lquin003/Seurat/learning_and_memory.txt",col_names = F,delim = "\t")%>%
  arrange(X1) %>% dplyr::select(X1,X2)

#load cadherin genes 
cadherin_genes = read_delim("/nas/longleaf/home/lquin003/Seurat/cadherins_list.txt",col_names = F,delim = "\t")%>%
  arrange(X1) %>% dplyr::select(X1,X2)

#Load in cancer genes
cancer_genes = read_delim("/nas/longleaf/home/lquin003/Seurat/cancer_gene_list.txt",col_names = F,delim = "\t")%>%
  arrange(X1) %>% dplyr::select(X1,X2)


#Select NT genes against our DEG list excluding CA clusters
nt_match <- heatmapdf %>% 
  filter(gene %in% nt_genes$X1,!cluster_number %in%CA_Excitatory)
#Select ach genes against our DEG list excluding CA clusters
ach_match <- heatmapdf %>% 
  filter(gene %in% ach_genes$X1, !cluster_number %in% CA_Excitatory)
#Select Calcium Signaling genes against our DEG list excluding CA clusters
ca_match <- heatmapdf %>%
  filter(gene %in% ca_sig_genes$X1, !cluster_number%in%CA_Excitatory)
#Select learn and mem genes against our DEG list and exclude CA clusters
learn_match =  heatmapdf %>%
  filter(gene %in% learn_n_mem_genes$X1, !cluster_number%in%CA_Excitatory)
## Cadherin match
cadherin_match = heatmapdf %>%
  filter(gene %in% cadherin_genes$X1, !cluster_number%in%CA_Excitatory)
## Cadherin match
cancer_match = heatmapdf %>%
  filter(gene %in% cancer_genes$X1, !cluster_number%in%CA_Excitatory)


#### Heat map cluster ####
#Parameters for Coloring Heatmap
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-3, 3))

#select DF and subset based off select columns 
subsets <- c("gene","cluster_name","avg_diff")

#Create matrix for heatmap that contains subsets
mat_ud <-ach_match %>% select(all_of(subsets)) %>%
  as_tibble() %>% #distinct(gene) %>%
  pivot_wider(names_from = cluster_name, values_from = avg_diff, values_fill = 0)%>% 
  tibble::column_to_rownames(var="gene") %>% as.matrix()

### Compex heatmap attempt at Clustering
h_ach_match_byclust <- Heatmap(mat_ud, name="Average\nDifference", col = col_fun, column_title = "Acetylcholine Receptors",
                        border = T,
                        rect_gp = gpar(col= "black"),
                        row_title_side = "left",
                        show_row_dend = F,
                        show_column_dend = F,
                        column_dend_reorder = T,
                        column_names_gp = gpar(fontsize=11),
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize=11))
#View heatplot
h_ach_match_byclust

pdf(file = paste0(output,"ach_downregulation_heatmap_byclust.pdf"),width = 4, height=3)
h_ach_match_byclust
dev.off()


#select DF and subset based off select columns 
subsets <- c("gene","cell_type","avg_diff")

#Create matrix for heatmap that contains subsets
mat_ud <-ach_match %>% select(all_of(subsets)) %>%
  as_tibble() %>% 
  # distinct(gene) %>%
  pivot_wider(names_from = cell_type, values_from = avg_diff, values_fill = 0)%>% 
  tibble::column_to_rownames(var="gene") %>% as.matrix()

### Compex heatmap attempt at Clustering
h_ach_match_bycelltype <- Heatmap(mat_ud, name="Average\nDifference", col = col_fun, column_title = "Acetylcholine Receptors",
                               border = T,
                               rect_gp = gpar(col= "black"),
                               column_order = c("Granule_cells","Mossy","Other"),
                               # Not changing GC name to remove dash because the column order changes
                               column_labels = c("Granule Cells","Mossy","Other"),
                               row_title_side = "left",
                               show_row_dend = F,
                               show_column_dend = F,
                               column_dend_reorder = T,
                               column_names_gp = gpar(fontsize=11),
                               row_names_side = "left",
                               row_names_gp = gpar(fontsize=11))
#View heatplot
h_ach_match_bycelltype

pdf(file = paste0(output,"ach_downregulation_heatmap_bycell_type.pdf"),width = 4, height=3)
h_ach_match_bycelltype
dev.off()



## Calcium Signaling Heatmap
###labeling colors for heatmap
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-3, 3))
#select DF and subset based off select columns 
subsets <- c("gene","cell_type","avg_diff")

#Create matrix for heatmap that contains subsets
mat_ud <-ca_match %>% group_by(cell_type) %>% distinct(gene,.keep_all = T) %>% ungroup() %>%
  select(all_of(subsets)) %>%
  as_tibble() %>% 
  pivot_wider(names_from = cell_type, values_from = avg_diff, values_fill = 0)%>% 
  tibble::column_to_rownames(var="gene") %>% as.matrix()

### Compex heatmap attempt at Clustering
heat_ca_sig_by_cell_type <- Heatmap(mat_ud, name="Average\nDifference", col = col_fun, column_title = "Calcium Signaling Genes",
                                  border = T,
                                  rect_gp = gpar(col= "black"),
                                  column_order = c("Granule_cells","GABA","Mossy","Neurogenic","Glia","Other"),
                                  # Not changing GC name to remove dash because the column order changes
                                  #column_labels = c("Granule cells","GABA","Mossy","Neurogenic","Glia","Other"),
                                  row_title_side = "left",
                                  show_row_dend = F,
                                  show_column_dend = F,
                                  column_dend_reorder = T,
                                  column_names_gp = gpar(fontsize=11),
                                  row_names_side = "left",
                                  row_names_gp = gpar(fontsize=10))
#View heatplot
heat_ca_sig_by_cell_type

pdf(file = paste0(output,"ca_sig_heatmap_bycell_type.pdf"),width = 3, height=7)
heat_ca_sig_by_cell_type
dev.off()


## Calcium Signaling Heatmap
###labeling colors for heatmap
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-3, 3))
#select DF and subset based off select columns 
subsets <- c("gene","cluster_name","avg_diff")

#Create matrix for heatmap that contains subsets

ca_row_names = ca_match %>% arrange(p_val) %>% distinct(gene) %>% pull(gene)

column_order_heat = ca_match %>% arrange(cell_type) %>% pull(cluster_name) %>% unique()

#ca_row_names <- ca_match %>% arrange(gene) %>% dplyr::select(gene) %>% distinct()
ca_col_names <- ca_match %>% group_by(cell_type) %>% arrange(cell_type) %>% ungroup() %>% dplyr::select(cluster_name) %>% distinct()

mat_ud <-ca_match %>% group_by(cluster_name) %>% distinct(gene,.keep_all = T) %>% ungroup() %>%
  dplyr::select(all_of(subsets)) %>%
  as_tibble() %>% 
  pivot_wider(names_from = cluster_name, values_from = avg_diff, values_fill = 0)%>% 
  tibble::column_to_rownames(var="gene") %>% as.matrix()

#Create right annotation of p_value
pvalue_col_fun = colorRamp2(c(0, 2, 3), c("green", "white", "red")) 
row_ha = rowAnnotation(p_val = )

### Compex heatmap attempt at Clustering
heat_ca_sig_by_cluster <- Heatmap(mat_ud, name="Average\nDifference", col = col_fun, column_title = "Calcium Signaling Genes by Cluster",
                                    border = T,
                                    rect_gp = gpar(col= "black"),
                                    column_order = ca_col_names$cluster_name,
                                    row_order = ca_row_names,
                                    row_title_side = "left",
                                    show_row_dend = F,
                                    show_column_dend = F,
                                    column_dend_reorder = T,
                                    column_names_gp = gpar(fontsize=11),
                                    row_names_side = "left",
                                    heatmap_legend_param = list(at=c(-2,-1,0,1,2),direction="horizontal",title_position = "topcenter"),
                                    row_names_gp = gpar(fontsize=10))
                                    
#View heatplot
draw(heat_ca_sig_by_cluster, heatmap_legend_side = "top", annotation_legend_side = "top")

#

#heat_ca_sig_by_cluster

pdf(file = paste0(output,"ca_sig_heatmap_bycluster.pdf"),width = 6, height=7)
draw(heat_ca_sig_by_cluster, heatmap_legend_side = "top", annotation_legend_side = "top")
dev.off()

## Neurotransmitter Heatmap
#Create matrix for heatmap that contains subsets
mat_ud <-nt_match %>% group_by(cell_type) %>% distinct(gene,.keep_all = T) %>% ungroup() %>%
  select(all_of(subsets)) %>%
  as_tibble() %>% 
  pivot_wider(names_from = cell_type, values_from = avg_diff, values_fill = 0)%>% 
  tibble::column_to_rownames(var="gene") %>% as.matrix()

### Compex heatmap attempt at Clustering
heat_nt_binding_celltype <- Heatmap(mat_ud, name="Average\nDifference", col = col_fun, column_title = "Neurotransmitter Binding Genes",
                                    border = T,
                                    rect_gp = gpar(col= "black"),
                                    column_order = c("Granule_cells","GABA","Mossy","Neurogenic","Glia", "Other"),
                                    # Not changing GC name to remove dash because the column order changes
                                    #column_labels = c("Granule cells","GABA","Mossy","Neurogenic","Glia","Other"),
                                    row_title_side = "left",
                                    show_row_dend = F,
                                    show_column_dend = F,
                                    column_dend_reorder = T,
                                    column_names_gp = gpar(fontsize=11),
                                    row_names_side = "left",
                                    row_names_gp = gpar(fontsize=10))
#View heatplot
heat_nt_binding_celltype

pdf(file = paste0(output,"nt_binding_heatmap_bycell_type.pdf"),width = 3, height=3)
heat_nt_binding_celltype
dev.off()

# Neurotransmitter Heatmap
subsets <- c("gene","cluster_name","avg_diff")
#Create matrix for heatmap that contains subsets by cluster number


nt_row_names <- nt_match %>% arrange(desc(gene)) %>% select(gene) %>% distinct()
nt_col_names <- nt_match %>% group_by(cell_type) %>% arrange(cell_type) %>% ungroup() %>% select(cluster_name) %>% distinct()

mat_ud <-nt_match %>% #group_by(cluster_name) %>% distinct(gene,.keep_all = T) %>% ungroup() %>%
  select(all_of(subsets)) %>%
  as_tibble() %>% 
  pivot_wider(names_from = cluster_name, values_from = avg_diff, values_fill = 0) %>% 
  tibble::column_to_rownames(var="gene") %>% as.matrix()

### Compex heatmap attempt at Clustering
heat_nt_binding_cluster <- Heatmap(mat_ud, name="Average\nDifference", col = col_fun, column_title = "Neurotransmitter Binding Genes",
                                    border = T,
                                    rect_gp = gpar(col= "black"),
                                    #column_order = c("Granule_cells","GABA","Mossy","Neurogenic","Other"),
                                    # Not changing GC name to remove dash because the column order changes
                                    column_order = nt_col_names$cluster_name,
                                    row_order = nt_row_names$gene,
                                    row_title_side = "left",
                                    show_row_dend = F,
                                    show_column_dend = F,
                                    column_dend_reorder = T,
                                    column_names_gp = gpar(fontsize=11),
                                    row_names_side = "left",
                                    row_names_gp = gpar(fontsize=10))
#View heatplot
heat_nt_binding_cluster

pdf(file = paste0(output,"nt_binding_heatmap_cluster.pdf"),width = 5, height=5)
heat_nt_binding_cluster
dev.off()

# Figure 5 C 
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



#Neurogenic Feature Plots
factor(ach_ng, levels = c(25,15,16))

VlnPlot(ach_ng, features = "Sox2", split.by = "trt",split.plot = T)

RidgePlot(ach_ng, features = c("Aldoc","Sox2" ), split.by = "trt")
DoHeatmap(ach_ng, features = c("Aldoc","Sox2" ))

### Find different markers in Cluster 15 and 25 and 16 ###
#When finding DEG's we should use CPM or counts per million as an assay instead of SCT or RNA
#These normalization methods aren't consistent with those used in the past

cluster_15_markers_RNA <- FindMarkers(object = ach, assay = "RNA", slot = "data", test.use = "wilcox", ident.1 = "15 RGL 1", ident.2 = "25 RGL Like", min.pct = .25)
cluster_15_markers_SCT <- FindMarkers(object = ach, assay = "SCT", slot = "data", test.use = "wilcox", ident.1 = "15 RGL 1", ident.2 = "25 RGL Like", min.pct = .25)

cluster_15_markers_trt <- FindMarkers(object = ach, 
                                      assay = "SCT",
                                      group.by = "trt",
                                      slot = "data", 
                                      test.use = "wilcox", 
                                      ident.1 = "chr2", 
                                      ident.2 = "yfp", 
                                      subset.ident = "25 RGL Like",
                                      min.pct = .25)

cluster_25_markers_trt <- FindMarkers(object = ach, 
                                      assay = "SCT",
                                      group.by = "trt",
                                      slot = "data", 
                                      test.use = "wilcox", 
                                      ident.1 = "chr2", 
                                      ident.2 = "yfp", 
                                      subset.ident = "25 RGL Like",
                                      min.pct = .25)




cluster_15_markers2 <- cluster_15_markers_SCT %>% 
  rownames_to_column("gene") %>% 
  as_tibble()%>%
  mutate(Expression = case_when(avg_log2FC >= log(2) & p_val_adj <= .05 ~ "Up_regulated",
                          avg_log2FC <= -log2(2) & p_val_adj <= .05 ~ "Down_regulated",
                          TRUE~ "Unchanged"))

top <- 15
top_genes <-  bind_rows(
  cluster_15_markers2 %>% 
    filter(Expression == 'Up_regulated') %>% 
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>% 
    head(top),
  cluster_15_markers2 %>% 
    filter(Expression == 'Down_regulated') %>% 
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>% 
    head(top)
)      

mytheme <- theme(panel.background = element_blank(),
      legend.position = "none",
      axis.line = element_line(size=.5),
      axis.text = element_text(size=12),
      axis.ticks = element_line(size=.5),
      axis.title = element_text(size=18))

c15_volc_p <- ggplot(cluster_15_markers2, aes(x=avg_log2FC, y=-log(p_val_adj,10)))+
  geom_point(aes(color=Expression),inherit.aes = T,stat = "identity", size = 1.5)+
  xlab("Average Log Fold Change")+
  ylab(expression("-log"[10]*"(P Value)"))+
  mytheme+
  scale_color_manual(values = c("blue2","black","red2"))+
  geom_text_repel(data=top_genes, mapping= aes(avg_log2FC,-log(p_val_adj), label=gene,size=3))

c15_volc_p


c15_volc_p + geom_label_repel(data = top_genes,
                              mapping = aes(avg_log2FC, -log(p_val_adj), label=gene), 
                              label.size = 0,
                              point.size = 5,
                              box.padding = .35,
                              point.padding = .5,)

ach %>%
  subset(seurat_clusters %in% c(15,25))%>%
  VlnPlot(features = c("Tenm2","Plp1","Sngh11","Meg3"), assay = "SCT",split.by = "trt")

####

cluster_25_markers_trt <- FindMarkers(object = ach, 
                                      assay = "SCT",
                                      group.by = "trt",
                                      slot = "data", 
                                      test.use = "wilcox", 
                                      ident.1 = "chr2", 
                                      ident.2 = "yfp", 
                                      subset.ident = "25 RGL Like",
                                      min.pct = .25)


ach %>%
  subset(seurat_clusters %in% c(15,25))%>%
  VlnPlot(features = "Ctnnd2", assay = "RNA",split.by = "trt")

DotPlot(dg.cpm, assay = "SCT", features = c("Cdh13"),cluster.idents = T ,cols = c("grey","goldenrod1","goldenrod4"))

c_15_expression_counts <- cluster_15_markers2 %>%
  count(Expression) 


cdh_filt = wrs_method2_data %>% filter(gene %in% c("Cdh6","Cdh9","Cdh10","Cdh2","Cdh13","Cdh3","Cdh8"))
write.csv(paste0(output,cluster_15_markers),file ="Cluster15_Markers_vs_Cluster25_Markers.csv")
#####
cluster_25_markers <- FindMarkers(object = ach, ident.1 = 25, ident.2 = 15, min.pct = .25)

### Different markers between 16 and 25
cluster_16_markers <- FindMarkers(object = ach, ident.1 = 16, ident.2 = 25, min.pct = .25)

cluster_15_markers_trt <- FindMarkers(object = ach_trt, ident.1 = 15, ident.2 = 25, min.pct = .25)

write.csv(paste0(output,cluster_16_markers),file ="Cluster16_Markers_vs_Cluster25_Markers.csv")



########## Correlation Plot from Jeremy 
# Remade on 12/21/21 to capture PDF
dg.integrated = readRDS("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered.rds")

vargenes<- presto::wilcoxauc(ach, 'seurat_clusters', seurat_assay = 'SCT')
top_vargenes = top_markers(vargenes, n = 500, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)
all_markers<- top_vargenes %>%
  dplyr::select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

#dg.integrated <- BuildClusterTree(dg.integrated, reorder = F, verbose = T, assay="SCT", features = all_markers)
#PlotClusterTree(dg.integrated)

dg.avg = AverageExpression(ach)
dg.sct.cor = cor(dg.avg$SCT[all_markers,])

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
col_fun(seq(-1, 1))

pdf("Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_prestoMarkers_SCT_pearson_cor_heatmap.pdf",height=8,width=10)
ComplexHeatmap::Heatmap(dg.sct.cor,name="Pearson correlation",column_names_side = "top",row_names_side = "left",col = col_fun, clustering_distance_columns=function(x) as.dist(1-cor(t(x))), clustering_distance_rows=function(x) as.dist(1-cor(t(x))), clustering_method_columns = "average", clustering_method_rows = "average")
dev.off()



