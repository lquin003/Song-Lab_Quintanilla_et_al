#assign(".lib.loc", "/nas/longleaf/home/lquin003/R/niche_net_libs/", envir = environment(.libPaths))

library(Seurat) # please update to Seurat V4
library(nichenetr,lib.loc = "/nas/longleaf/home/lquin003/R/niche_net_libs/")
library(tidyverse)
library(SeuratObject)
library(dplyr)
library(circlize)
library(patchwork)
library(RColorBrewer)
library(cowplot)
library(ggpubr)


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

head(weighted_networks_lr)

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
title = "nsc Niche Signaling"

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

DE_table_receiver = FindMarkers(object = seurat_obj_receiver,assay= "RNA" ,slot = "scale.data", ident.1 = condition_oi, ident.2 = condition_reference, logfc.threshold =0.1) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val <= 0.01 & pct.1 >=.1) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

## Sender list of genes identified using getxpressed genes function
list_expressed_genes_gc = Granule_cells %>% unique() %>% lapply(get_expressed_genes,seuratObj, pct = 0.10, assay_oi="SCT")
expressed_genes_gc = list_expressed_genes_gc %>% unlist() %>% unique()
str(expressed_genes_gc)

list_expressed_genes_mossy = Mossy %>% unique() %>% lapply(get_expressed_genes,seuratObj, pct = 0.10, assay_oi="SCT")
expressed_genes_mossy = list_expressed_genes_mossy %>% unlist() %>% unique()
str(expressed_genes_mossy)

list_expressed_genes_GABA = GABA %>% unique() %>% lapply(get_expressed_genes,seuratObj, pct = 0.10, assay_oi="SCT")
expressed_genes_GABA = list_expressed_genes_GABA %>% unlist() %>% unique()
str(expressed_genes_GABA)

list_expressed_genes_glia = Glia %>% unique() %>% lapply(get_expressed_genes,seuratObj, pct = 0.10, assay_oi="SCT")
expressed_genes_glia = list_expressed_genes_glia %>% unlist() %>% unique()
str(expressed_genes_glia)

#list_expressed_genes_ng = Neurogenic %>% unique() %>% lapply(get_expressed_genes,seuratObj, pct = 0.10, assay_oi="SCT")
#expressed_genes_ng = list_expressed_genes_ng %>% unlist() %>% unique()
#str(expressed_genes_ng)

list_expressed_genes_other = Other %>% unique() %>% lapply(get_expressed_genes,seuratObj, pct = 0.10, assay_oi="SCT")
expressed_genes_other = list_expressed_genes_GABA %>% unlist() %>% unique()
str(expressed_genes_other)

# Identify sender cell types of interest in order to determine how each is communicating with target or receiver cell
ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands_gc = dplyr::intersect(ligands,expressed_genes_gc)
expressed_ligands_mossy = dplyr::intersect(ligands,expressed_genes_mossy)
expressed_ligands_GABA = dplyr::intersect(ligands,expressed_genes_GABA)
expressed_ligands_glia = dplyr::intersect(ligands,expressed_genes_glia)
#expressed_ligands_ng = dplyr::intersect(ligands,expressed_genes_ng)
expressed_ligands_other = dplyr::intersect(ligands,expressed_genes_other)
expressed_ligands = union(expressed_ligands_gc,expressed_ligands_mossy) %>% union(expressed_ligands_GABA) %>% union(expressed_ligands_glia) %>% union(expressed_ligands_other)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = dplyr::intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
head(potential_ligands)

#Predict ligands using nichenet function
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities %>% arrange(-pearson) 

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

best_upstream_ligands %>% dplyr::intersect(expressed_ligands_gc)
best_upstream_ligands %>% dplyr::intersect(expressed_ligands_mossy)
best_upstream_ligands %>% dplyr::intersect(expressed_ligands_GABA)
#best_upstream_ligands %>% dplyr::intersect(expressed_ligands_ng)
best_upstream_ligands %>% dplyr::intersect(expressed_ligands_glia)
best_upstream_ligands %>% dplyr::intersect(expressed_ligands_other)

# lot of overlap between both cell types in terms of expressed ligands
# therefore, determine which ligands are more strongly expressed in which of the two
#Get expression from seurat object 
ligand_expression_tbl = tibble(
  ligand = best_upstream_ligands,
  gc_expression = GetAssayData(object = subset(seuratObj,seurat_clusters %in% Granule_cells),assay="SCT",slot="data")[best_upstream_ligands,] %>% rowMeans(),
  mossy_expression = GetAssayData(object = subset(seuratObj,seurat_clusters %in% Mossy),assay="SCT",slot="data")[best_upstream_ligands,] %>% rowMeans(),
  GABA_expression = GetAssayData(object = subset(seuratObj,seurat_clusters%in% GABA),assay="SCT",slot="data")[best_upstream_ligands,] %>% rowMeans(),
  glia_expression = GetAssayData(object = subset(seuratObj,seurat_clusters%in% Glia),assay="SCT",slot="data")[best_upstream_ligands,] %>% rowMeans(),
  #ng_expression = GetAssayData(object = subset(seuratObj,seurat_clusters%in% Neurogenic),assay="SCT",slot="data")[best_upstream_ligands,] %>% rowMeans(),
  other_expression = GetAssayData(object = subset(seuratObj,seurat_clusters%in% Other),assay="SCT",slot="data")[best_upstream_ligands,] %>% rowMeans()
)

# Identify cell type specific ligands by determining which cell type has the highest expression of ligand
gc_specific_ligands = ligand_expression_tbl %>% filter(gc_expression > mossy_expression & gc_expression > GABA_expression & gc_expression > glia_expression & gc_expression > other_expression) %>% pull(ligand)
mossy_specific_ligands = ligand_expression_tbl %>% filter(mossy_expression > gc_expression & mossy_expression > GABA_expression & mossy_expression > glia_expression & mossy_expression > other_expression) %>% pull(ligand)
gaba_specific_ligands = ligand_expression_tbl %>% filter(GABA_expression > mossy_expression & GABA_expression > gc_expression & GABA_expression > glia_expression & GABA_expression > other_expression) %>% pull(ligand)
glia_specific_ligands = ligand_expression_tbl %>% filter(glia_expression > mossy_expression & glia_expression > GABA_expression & glia_expression > gc_expression & glia_expression > other_expression) %>% pull(ligand)
#ng_specific_ligands = ligand_expression_tbl %>% filter(ng_expression > mossy_expression & ng_expression > GABA_expression & ng_expression > gc_expression & ng_expression > other_expression & ng_expression > gc_expression) %>% pull(ligand)
other_specific_ligands = ligand_expression_tbl %>% filter(other_expression > mossy_expression & other_expression > GABA_expression & other_expression > glia_expression & other_expression > gc_expression) %>% pull(ligand)
general_ligands = setdiff(best_upstream_ligands,c(gc_specific_ligands,mossy_specific_ligands,gaba_specific_ligands,glia_specific_ligands,
                                                  #,ng_specific_ligands
                                                    other_specific_ligands))

# Create tibble that has information order of ligands is important to match them correctly
ligand_type_indication_df = tibble(
  ligand_type = c(rep("GC-specific", times = gc_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length()),
                  rep("Mossy-specific", times = mossy_specific_ligands %>% length()),
                  rep("GABA-specific", times = gaba_specific_ligands %>% length()),
                  rep("glia-specific", times = glia_specific_ligands %>% length()),
                  #rep("ng-specific", times = ng_specific_ligands %>% length()),
                  rep("other-specific", times = other_specific_ligands %>% length())),
  ligand = c(gc_specific_ligands,general_ligands,mossy_specific_ligands,gaba_specific_ligands,glia_specific_ligands,
             #ng_specific_ligands,
             other_specific_ligands))

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "RGL") %>% inner_join(ligand_type_indication_df) %>% na.omit()# if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type

cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.67)

active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

##Identify ligands and targets to remove
ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

#remove 1
circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

#Add 3 colors to match ligand type data frame
grid_col_ligand =c("GC-specific" = "turquoise3",
                   "General" = "orange",
                   "Mossy-specific" = "darkorchid3",
                   #"ng-specific" = "darkolivegreen",
                   "GABA-specific" = "goldenrod3",
                   "glia-specific" = "turquoise3",
                   "other-specific" = "steelblue3")

grid_col_target = c("RGL" = "red1")


grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)

grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

#remove 2 suffix
circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
#remove 3 suffix
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand,target, weight)


ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

transperancy2= circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% select(all_of(c("target","transparency")))

target_order = circos_links$target %>% unique()
ligand_order = c(gc_specific_ligands,general_ligands,mossy_specific_ligands,gaba_specific_ligands,glia_specific_ligands,
                 #ng_specific_ligands,
                 other_specific_ligands) %>% c(paste(.," ")) %>% dplyr::intersect(circos_links$ligand)
order = c(ligand_order,target_order)


#repare the circos visuaizaiton define the gaps between the different segments
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

gaps = c(
  # width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "GC-specific") %>% distinct(ligand) %>% nrow() -1)),
  #width_different_cell,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  #width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "glia-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "ng-specific") %>% distinct(ligand) %>% nrow() -1)),
  #width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "other-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Mossy-specific") %>% distinct(ligand) %>% nrow() -1)), 
  #width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "GABA-specific") %>% distinct(ligand) %>% nrow() -1)), 
  width_different_cell,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "RGL") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
)

circos.par(gap.degree = gaps,track.height=3)
chordDiagram(links_circle, directional = 1,order = order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col, transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1.5,col = "black")
}, bg.border = NA) #

circos.clear()

svg("ligand_target_circos_RGL_genes.svg", width = 10, height = 10)
circos.par(gap.degree = gaps,track.height =3)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transperancy, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1.5)
}, bg.border = NA) #
circos.clear()
dev.off()


## png 
##   2

#second Circos plot showing ligand receptor interactions instead of ligand gene interactions
##Visualize Ligand-receptor interactions of prioritized ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

#get the weights of the ligand-receptor interactions as used in the NicheNet model
## Possibly edit this next line since they used human and not mouse genes
#Import list of best upstream receptors and 
lr_network_top_df = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors) %>% rename(ligand = from, receptor = to)

lr_network_top_df = lr_network_top_df %>% mutate(receptor_type = "RGL") %>% inner_join(ligand_type_indication_df)

#Add 3 colors to match ligand type data frame
grid_col_ligand =c("GC-specific" = "turquoise1",
                   "General" = "orange",
                   "Mossy-specific" = "darkorchid3",
                   "GABA-specific" = "goldenrod3",
                   "glia-specific" = "turquoise3",
                   #"ng-specific" = "darkolivegreen",
                   "other-specific" = "steelblue3")

grid_col_receptor = c("RGL" = "firebrick4")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_receptor = tibble(receptor_type = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)


circos_links = lr_network_top_df %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as receptor!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_receptor)
 

links_circle = circos_links %>% dplyr::select(ligand,receptor, weight)

#takes color from ligand and color ligand type column
ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
#Takes ligand color and sets names as ligands
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
# Takes receptor from receptor type column
receptor_color = circos_links %>% distinct(receptor,color_receptor_type)
#Changes names to color form recpetor grid column
grid_receptor_color = receptor_color$color_receptor_type %>% set_names(receptor_color$receptor)

#creats a color df combined of both ligand and receptor colors
grid_col =c(grid_ligand_color,grid_receptor_color)


# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

receptor_order = circos_links$receptor %>% unique()
ligand_order = c(gc_specific_ligands,general_ligands,mossy_specific_ligands,gaba_specific_ligands,glia_specific_ligands,other_specific_ligands) %>% c(paste(.," ")) %>% dplyr::intersect(circos_links$ligand)
order = c(ligand_order,receptor_order)

width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_receptor = 15
width_same_cell_same_receptor_type = 0.5


gaps = c(
  # width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "GC-specific") %>% distinct(ligand) %>% nrow() -1)),
  #width_different_cell,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  #width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Mossy-specific") %>% distinct(ligand) %>% nrow() -1)), 
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "GABA-specific") %>% distinct(ligand) %>% nrow() -1)), 
  width_ligand_receptor,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "glia-specific") %>% distinct(ligand) %>% nrow() -1)), 
  width_ligand_receptor,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "ng-specific") %>% distinct(ligand) %>% nrow() -1)), 
  #width_ligand_receptor,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "other-specific") %>% distinct(ligand) %>% nrow() -1)), 
  width_ligand_receptor,
  rep(width_same_cell_same_receptor_type, times = (circos_links %>% filter(receptor_type == "RGL") %>% distinct(receptor) %>% nrow() -1)),
  width_ligand_receptor
)

plot_title = paste0(output,"ligand_receptor_circos_in_RGL.svg")

svg(filename = "ligand_receptor_RGL.svg",width = 10, height = 10)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col, transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = .8)
}, bg.border = NA) #

circos.clear()
dev.off()

