# Song-Lab Quintanilla et al. 
This Repo contains the code used to process and analyze snRNA-seq data using Split-seq from Quintanilla et al.

The following scripts include preprocessing, normalization, differential gene expression analysis, pseudotime analysis and Nichenet analysis.

For questions reach out to juansong@email.unc.edu or 01luisquintanilla@gmail.com

## List of Codes and brief explanation

Ach_data_set_analysis.R - Used to create and annotate initial clusters includes explorative data not included in manuscript.

Cluster_Annotationsn_R2.txt - Name of assigned cluster identities using zUMIs.

DG_SPLITseq_zUMIS_combined_seurat.R - Code used after generating gene x cell matrix from zUMIS to rename cells according to barcodes.

Nichenet_seurat.R - Mirrors nichenet analysis from Nichenet Repo using a Seurat Object from our study.

circos_plot_my_data.R - Creates circle plots using Nichenet generated ligand-receptor associations

combined_DEG_scripts.3.31.22.R - Contains the most comperhensive script used to generate most figures within the manuscript. It includes DEG analysis.

slinghsot.R - Used to generate pseudotime figures in the manuscript using Rv3.6.

slingshot_smoothers.R - Used to generate pseudotime plots of different genes using Rv4.1
