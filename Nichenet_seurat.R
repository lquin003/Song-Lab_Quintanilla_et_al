#assign(".lib.loc", "/nas/longleaf/home/lquin003/R/niche_net_libs/", envir = environment(.libPaths))

library(Seurat,lib.loc = "/nas/longleaf/rhel8/apps/r/4.1.0/lib64/R/library") # please update to Seurat V4
library(nichenetr)
library(tidyverse)
library(SeuratObject)
library(circlize)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(cowplot)
library(tidyr)
library(ggpubr)
library(purrr)

#set output for graphs and other files use this paste0 output before each file name to use this
output <- "/nas/longleaf/home/lquin003/Nicenet/Output_N/"

#Read in seurat object
dg.cpm = readRDS(file ="/proj/jsonglab/projects/Luis/single_cell_analysis/Song_DG_SPLITseq_zUMIs_scTransform_052021_integrated_clustered_wrs2_deg_cpm_new.rds")

# Assign order for future Vln plots contorl and treatment
dg.cpm$trt <- factor(dg.cpm$trt, levels = c("yfp","chr2"))
dg.cpm$trt

seuratObj = dg.cpm
#Remove duplicate Seurat object to save memory
rm(dg.cpm)
gc()

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

seuratObj = assign_cluster_names(seuratObj)

##### Read in Niche Net Ligand-target piror model, lignd receptor network and weighthed integrated networks ####

ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5]

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

#### Convert from Human to Mouse Symbols ####
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

head(weighted_networks$lr_sig)



##### Define Sender/Niche and Receiver cells cell population present in your expression data and determine which genes are expressed in both populations ####
# Load in Cluster names



cluster_names <- read.csv("/nas/longleaf/home/lquin003/Seurat/Cluster_Annotations_R2.txt")
print(cluster_names)
cluster_names <- cluster_names[1:43,]
cluster_names <- as.vector(cluster_names)


#Cluster assignment using cell type for all clusters
Granule_cells <- c(0,1,2,3,4,5,7,8,10,11,14,26,30,32)
GABA <- c(18,19,20,23,40)
Mossy <- c(17,21,24)
Subiculum <- c(9,33,38,39,41)
CA_Excitatory <- c(6,12,13,22,28,34)
Neurogenic <- c(25,16)
Glia <- c(27,31,35,36)
Other<- c(29,37,42)
Endothelial = (29)
RGL= 15


##Create list of sender and receivers
sender = c(Granule_cells,GABA,Mossy,Glia,Endothelial)
receiver = c(RGL)
title = "GC Niche Signaling"
  
seuratObj@active.ident = seuratObj$seurat_clusters
seuratObj@active.assay = "RNA"
    
## receiver
receiver = receiver
list_expressed_genes_receiver = receiver %>% unique() %>% lapply(get_expressed_genes,seuratObj, pct = 0.10, assay_oi="SCT")
expressed_genes_receiver = list_expressed_genes_receiver %>% unlist() %>% unique()

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## Sender
sender_celltypes = sender

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.10, assay_oi="SCT") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

##### Define gene set of interest: These are genes in the receiver/target cell ppulation that are potentially affected by ligands expressed by interacting cells #####
seurat_obj_receiver= subset(seuratObj, idents = receiver)

#select cells based on column you wish to split data by
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["trt"]])

#treatment condition
condition_oi = "chr2"
#control condition
condition_reference = "yfp" 

DE_table_receiver = FindMarkers(object = seurat_obj_receiver,assay= "RNA" ,slot = "scale.data", ident.1 = condition_oi, ident.2 = condition_reference, logfc.threshold =0.1 ) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val <= 0.01 & pct.1 >= .1) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]


##### Define a set of potential ligands: Thse are ligands that are expressed by the "sender/Niche" cell population and bind a (putative) receptor expressed by the receiver/target population #####
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = dplyr::intersect(ligands,expressed_genes_sender)
expressed_receptors = dplyr::intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()



##### Perform nichenet ligand activity analysis: rank potential ligands based on the presence of their target genes in gene set of itnerest compared to background genes #####
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson))) %>% mutate(cell_class = title)

##DF of best ligand activites sorted as person
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

#Change idents before making graphs
names(cluster_names) <- levels(seuratObj)
seuratObj <- RenameIdents(seuratObj, cluster_names)

#p1 = DotPlot(seuratObj,assay = "SCT" ,features = best_upstream_ligands %>% rev(),cols = "RdYlBu") + RotatedAxis() +labs(x="Upstream Ligands") +theme(axis.title.y = element_blank(),axis.text.x = element_text(angle=90)) + NoLegend() + labs(title = title)

#pdf(paste0(output, title,"upstream ligand.pdf"),height = 10,width = 10)
#p1
#dev.off()

##### Infer receptors and top predicted target genes of ligands that are top-ranked in ligand activity analysis #####
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
##DF of active ligand targers

#Edit cutoff to include all top 20 N for receptor ligands and targeted genes
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.3)
order_ligands = dplyr::intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()


order_targets = active_ligand_target_links_df$target %>% unique() %>% dplyr::intersect(rownames(active_ligand_target_links)) %>% make.names()

rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

#p2= vis_ligand_target %>% make_heatmap_ggplot("", "", color = "red",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "red", breaks = c(0,0.0045,0.0090)) +geom_tile(color="black",lwd=.75,linetype=1)

p2 = vis_ligand_target %>% make_threecolor_heatmap_ggplot("","", low_color = "blue",mid_color = "white", mid = median(vis_ligand_target), high_color = "orange",legend_position = "top", x_axis_position = "top", legend_title = "Predicted Target \nPotential") +
  geom_tile(color="black",lwd=.75,linetype=1) + 
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.text.x = element_text(face="italic")) +  
  scale_fill_gradient2(low = "whitesmoke",  high = "orange",breaks=c(0,.004,.007))

p2+coord_fixed()
##Added code for Ligands from Best Ligands list for strict receptor pairs

lr_network_top = lr_network %>% filter(from %in% order_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% order_ligands & to %in% best_upstream_receptors)


lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% dplyr::intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% dplyr::intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

vis_ligand_receptor_network = vis_ligand_receptor_network %>% t() 
vis_ligand_receptor_network = vis_ligand_receptor_network[order_ligands,]

#ligand receptor network plot
p3 = vis_ligand_receptor_network %>% make_threecolor_heatmap_ggplot("","", low_color = "white",mid_color = "white", mid = mean(vis_ligand_receptor_network), high_color = "darkgreen",legend_position = "top", x_axis_position = "top", legend_title = "Receptor Interaction Potnetial") +geom_tile(color="black",lwd=.75,linetype=1) + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())

p3

##### Receptors of Top-ranked ligands ####
#lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
#best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()


#lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)



#lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
#lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

#dist_receptors = dist(lr_network_top_matrix, method = "binary")
#hclust_receptors = hclust(dist_receptors, method = "ward.D2")
#order_receptors = hclust_receptors$labels[hclust_receptors$order]

#dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
#hclust_ligands = hclust(dist_ligands, method = "ward.D2")
#order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

#order_receptors = order_receptors %>% dplyr::intersect(rownames(lr_network_top_matrix))
#order_ligands_receptor = order_ligands_receptor %>% dplyr::intersect(colnames(lr_network_top_matrix))

#vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
#rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
#colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

#ligand receptor network plot
#p3= vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands of Cholinergic Signaling from Niche",paste("Receptors of Ligands from Cholinergic Signaling in", title), color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")+geom_tile(color="black",lwd=.75,linetype=1)
#p3


## Receptors of top ranked Ligands but using only bona fide ligand-receptor interactions from documented literature
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

#
lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% dplyr::intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% dplyr::intersect(colnames(lr_network_top_matrix_strict))


vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()

#ligand receptor network strict
p4 = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("","" ,color = "green", x_axis_position = "top",legend_title = "Interaction Potential") +geom_tile(color="black",lwd=.75,linetype=1)
p4 +coord_fixed()

#### Add log fold change Information of Ligands from sender cells ####

#### Follow Up analysis 2 Visualize expression of Top-predicted ligands and their target genes in a combined heatmap

ligand_pearson_matrix = ligand_activities %>%  mutate(zscore = ((pearson-(mean(pearson)))/sd(pearson))) %>% dplyr::select(zscore) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

vis_ligand_pearson = ligand_pearson_matrix[order_ligands,] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")

#ligand pearson plot
p5 = vis_ligand_pearson %>% make_heatmap_ggplot("","", color = "royalblue",legend_position = "top", x_axis_position = "bottom", legend_title = "Zscore of \nLigand Activity")+theme(axis.title.x = element_blank(),axis.text.x = element_blank()) +geom_tile(color="black",lwd=.75,linetype=1)+ theme(legend.key.size = unit(.5,"cm"),axis.ticks.x = element_blank(),legend.spacing.x = unit(.1,"cm"),legend.text = element_text(size=8),legend.title = element_text(size=8))

p5

### Aad log fold change to Sender heat map plot
names(cluster_names) <- levels(seuratObj)
seuratObj <- RenameIdents(seuratObj, cluster_names)

sender_celltypes = sender
sender_celltypes = cluster_names[sender_celltypes+1]

# this uses a new nichenetr function - reinstall nichenetr if necessary!
DE_table_all = sender_celltypes %>% lapply(get_lfc_celltype, seurat_obj = seuratObj, condition_colname = "trt", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10,celltype_col ="cluster_name") %>% reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% dplyr::select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% dplyr::select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% str_replace("X\\d*_","") %>% make.names()

p6 = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("","", low_color = "blue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.x = element_text(size=8))+geom_tile(color="black",lwd=.75,linetype=1)
p6

#write.csv(ligand_activities,paste(output,"Ligand_activities_",title,".csv"))

# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% dplyr::select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

# ligand expression Seurat dotplot
order_ligands_adapted = order_ligands

names(cluster_names) <- levels(seuratObj)
seuratObj <- RenameIdents(seuratObj, cluster_names)

p7 = DotPlot(seuratObj %>% subset(seurat_clusters %in% sender),assay = "SCT" ,features = order_ligands_adapted, cols = "RdBu") +coord_flip()  + theme(legend.text = element_text(size = 10),axis.text.y = element_blank(), panel.border = element_rect(color="black",size=1), axis.ticks.y = element_blank() ,axis.title.y = element_blank(),axis.title.x = element_blank() ,axis.text.x = element_text(size = 8, angle=90, hjust=.5, vjust=.1), legend.position = "top",legend.title = (element_text(size = 8)),legend.key.size = unit(.5,"cm")) +scale_y_discrete(position = "right") + guides(color = guide_colorbar(title.position = "top",title = "Average Expression"))
p7
p7 + coord_fixed()




p5 = p5 + coord_fixed()
p7 = p7 + coord_fixed()
p3 = p3  + coord_fixed()
p2 = p2 + coord_fixed()
p6 = p6 + coord_fixed()

pdf(paste0(output,"Niche to RGL plots_nichenet 1.pdf"),width = 22, height = 6)
(p5 + p7 + p3 + p2) +plot_layout(widths = c(.2,3,3,1),nrow = 1) 
dev.off()

pdf(paste0(output,"Niche_plots_nichenet 2.pdf"),width = 10, height = 6)
p6
dev.off()

pdf("GC_plots_nichenet.pdf",width = 18, height = 6)
(p5 + p7 + p3+ p2+ p6)+plot_layout(widths = c(.3,3,3,2,2), nrow = 1) 
dev.off()

pdf("GABA_plots_nichenet.pdf",width = 18, height = 6)
(p5 + p7 + p3+ p2+ p6)+plot_layout(widths = c(.3,3,3,2,2), nrow = 1) 
dev.off()

pdf("Mossy_plots_nichenet.pdf",width = 18, height = 6)
(p5 + p7 + p3+ p2+ p6)+plot_layout(widths = c(.3,3,3,2,2), nrow = 1) 
dev.off()

pdf("Glia_plots_nichenet.pdf",width = 18, height = 6)
(p5 + p7 + p3+ p2+ p6)+plot_layout(widths = c(.3,3,3,2,2), nrow = 1) 
dev.off()

pdf("Other_plots_nichenet.pdf",width = 18, height = 6)
(p5 + p7 + p3+ p2+ p6)+plot_layout(widths = c(.3,3,3,2,2), nrow = 1) 
dev.off()



