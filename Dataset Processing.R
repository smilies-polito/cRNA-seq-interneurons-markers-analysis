library(monocle3)
library(Seurat)
library(SeuratObject)
library(SeuratWrappers)
library(dplyr)
library(Matrix)






allen_human_acc <- read.csv("C:/Users/loren/single cell/Data/Allen/Human/ACC/human_ACC_2018-10-04_exon-matrix.csv")
gene_allen_acc <- read.csv("C:/Users/loren/single cell/Data/Allen/Human/ACC/human_ACC_2018-10-04_genes-rows.csv")
rownames(allen_human_acc) <- gene_allen_acc$gene
allen_human_acc[1] <- NULL

allen_human_AccSEU <- CreateSeuratObject(allen_human_acc)
allen_human_AccCDS <- as.cell_data_set(allen_human_AccSEU)

allen_human_AccCDS <- detect_genes(allen_human_AccCDS)
allen_human_AccCDS <- estimate_size_factors(allen_human_AccCDS)
allen_human_AccCDS <- preprocess_cds(allen_human_AccCDS, method = "LSI")
allen_human_AccCDS <- reduce_dimension(allen_human_AccCDS, reduction_method = 'UMAP', preprocess_method = "LSI")
allen_human_AccCDS <- cluster_cells(allen_human_AccCDS, resolution = 2.5e-4)

plot_cells(allen_human_AccCDS)

rowData(allen_human_AccCDS)$gene_short_name <- rownames(rowData(allen_human_AccCDS))
plot_cells(allen_human_AccCDS, genes = "GAD2")
saveRDS(allen_human_AccCDS, "C:/Users/loren/single cell/Processed CDS/allen_human_ACC")


######################################


allen_human <- read.csv("C:/Users/loren/single cell/Data/Allen/Human/human_MTG_2018-06-14_exon-matrix.csv")
gene_allen <-  read.csv("C:/Users/loren/single cell/Data/Allen/Human/human_MTG_2018-06-14_genes-rows.csv")
cell_allen <-  read.csv("C:/Users/loren/single cell/Data/Allen/Human/human_MTG_2018-06-14_samples-columns.csv")


rownames(allen_human) <- gene_allen$gene
allen_human[1] <- NULL

allen_humanSEU <- CreateSeuratObject(allen_human)

allen_humanCDS <- as.cell_data_set(allen_humanSEU)

allen_humanCDS <- detect_genes(allen_humanCDS)
allen_humanCDS <- estimate_size_factors(allen_humanCDS)
allen_humanCDS <- preprocess_cds(allen_humanCDS, method = "LSI")
allen_humanCDS <- reduce_dimension(allen_humanCDS, reduction_method = 'UMAP', preprocess_method = "LSI")
allen_humanCDS <- cluster_cells(allen_humanCDS, resolution = 1e-5)

cell_allen <- cell_allen[-1,]
colData(allen_humanCDS)$clus <- cell_allen$cluster


# this commands are for plotting the results
plot_cells(allen_humanCDS)

# this command puts the genes names in a partciular column of the rowData to allow the visualization of gene expression
rowData(allen_humanCDS)$gene_short_name <- rownames(rowData(allen_humanCDS))
plot_cells(allen_humanCDS, genes = "GAD2")

saveRDS(allen_human_AccCDS, "C:/Users/loren/single cell/Processed CDS/allen_human_MTG")

############


allen_human_visp <- read.csv("C:/Users/loren/single cell/Data/Allen/Human/viusual cortex Visp/human_VISp_2018-10-04_exon-matrix.csv")
gene_allen_V <-  read.csv("C:/Users/loren/single cell/Data/Allen/Human/viusual cortex Visp/human_VISp_2018-10-04_genes-rows.csv")
cell_allen_v <-  read.csv("C:/Users/loren/single cell/Data/Allen/Human/viusual cortex Visp/human_VISp_2018-10-04_samples-columns.csv")


rownames(allen_human_visp) <- gene_allen_V$gene
allen_human_visp[1] <- NULL


allen_human_VispSEU <- CreateSeuratObject(allen_human_visp)
allen_human_VispCDS <- as.cell_data_set(allen_human_VispSEU)

allen_human_VispCDS <- detect_genes(allen_human_VispCDS)
allen_human_VispCDS <- estimate_size_factors(allen_human_VispCDS)
allen_human_VispCDS <- preprocess_cds(allen_human_VispCDS, method = "LSI")
allen_human_VispCDS <- reduce_dimension(allen_human_VispCDS, reduction_method = 'UMAP', preprocess_method = "LSI")
allen_human_VispCDS <- cluster_cells(allen_human_VispCDS, resolution = 1e-3)

plot_cells(allen_human_VispCDS)

rowData(allen_human_VispCDS)$gene_short_name <- rownames(rowData(allen_human_VispCDS))
plot_cells(allen_human_VispCDS, genes = "GAD2")

saveRDS(allen_human_VispCDS, "C:/Users/loren/single cell/Processed CDS/allen_human_VISP")

#########################################################


allen_brainCDS <- as.cell_data_set(allen_brain)


allen_brainCDS <- detect_genes(allen_brainCDS)
allen_brainCDS <- estimate_size_factors(allen_brainCDS)
allen_brainCDS <- preprocess_cds(allen_brainCDS, method = 'LSI', norm_method = 'none')
allen_brainCDS <- reduce_dimension(allen_brainCDS, reduction_method = 'UMAP', preprocess_method = "LSI")
allen_brainCDS <- cluster_cells(allen_brainCDS, resolution = 1e-4)

rowData(allen_brainCDS)$gene_short_name <- rownames(rowData(allen_brainCDS))
plot_cells(allen_brainCDS)
plot_cells(allen_brainCDS, genes = "Adarb2")
plot_cells(allen_brainCDS, color_cells_by = "partition")
top_allen_brain <- top_markers(allen_brainCDS)

saveRDS(allen_brainCDS, "C:/Users/loren/single cell/Processed CDS/allen_mouse")


#######################################


SCP2 <- readMM("C:/Users/loren/single cell/Data/Single cell portal/Mouse/2/expression_matrix.mtx")
SCP2_features <- read.table("C:/Users/loren/single cell/Data/Single cell portal/Mouse/2/genes.tsv")
SCP2_barcodes <- read.table("C:/Users/loren/single cell/Data/Single cell portal/Mouse/2/barcodes1.tsv")

#SCP2_metadata <- read.csv("C:/Users/loren/single cell/Data/Single cell portal/Mouse/2/Book1.csv")

#SCP2_metadata <- SCP2_metadata[-1,]
#SCP2_metadata <- SCP2_metadata[-1,]
colnames(SCP2_features)[1] <- "gene_short_name" 
colnames(SCP2) <- SCP2_barcodes$V1
row.names(SCP2) <- SCP2_features$gene_short_name

SCP2CDS <- new_cell_data_set(SCP2)

SCP2CDS <- detect_genes(SCP2CDS)
SCP2CDS <- estimate_size_factors(SCP2CDS)
SCP2CDS <- preprocess_cds(SCP2CDS, method = "LSI")

SCP2CDS <- reduce_dimension(SCP2CDS, preprocess_method = "LSI")
SCP2CDS <- cluster_cells(SCP2CDS, resolution = 1e-3)

rowData(SCP2CDS)$gene_short_name <- rownames(rowData(SCP2CDS))
#colData(SCP2CDS)$class <- SCP2_metadata$Column5
plot_cells(SCP2CDS)
plot_cells(SCP2CDS, genes = "Gad1")
#plot_cells(SCP2CDS, color_cells_by = "class")


saveRDS(SCP2CDS, "C:/Users/loren/single cell/Processed CDS/SCP2")

############################################


SCP3 <- readMM("C:/Users/loren/single cell/Data/Single cell portal/Mouse/ORB/gene_sorted-matrix.mtx")
SCP3_features <- read.table("C:/Users/loren/single cell/Data/Single cell portal/Mouse/ORB/genes.tsv")
SCP3_barcodes <- read.table("C:/Users/loren/single cell/Data/Single cell portal/Mouse/ORB/barcodes.tsv")

#SCP3_metadata <- read.csv("C:/Users/loren/single cell/Data/Single cell portal/Mouse/2/Book1.csv")

#SCP3_metadata <- SCP3_metadata[-1,]
#SCP3_metadata <- SCP3_metadata[-1,]
colnames(SCP3_features)[2] <- "gene_short_name" 
colnames(SCP3) <- SCP3_barcodes$V1
row.names(SCP3) <- SCP3_features$gene_short_name

SCP3CDS <- new_cell_data_set(SCP3)

SCP3CDS <- detect_genes(SCP3CDS)
SCP3CDS <- estimate_size_factors(SCP3CDS)
SCP3CDS <- preprocess_cds(SCP3CDS, method = "LSI")

SCP3CDS <- reduce_dimension(SCP3CDS, preprocess_method = "LSI")
SCP3CDS <- cluster_cells(SCP3CDS, resolution = 1e-3)

rowData(SCP3CDS)$gene_short_name <- rownames(rowData(SCP3CDS))
#colData(SCP3CDS)$class <- SCP3_metadata$Column5
plot_cells(SCP3CDS)
plot_cells(SCP3CDS, genes = "Gad1")
#plot_cells(SCP3CDS, color_cells_by = "class")


saveRDS(SCP3CDS, "C:/Users/loren/single cell/Processed CDS/SCP3")


#################################

SCP5 <- readMM("C:/Users/loren/single cell/Data/Single cell portal/Mouse/Visp/gene_sorted-matrix.mtx")
SCP5_features <- read.table("C:/Users/loren/single cell/Data/Single cell portal/Mouse/Visp/genes.tsv")
SCP5_barcodes <- read.table("C:/Users/loren/single cell/Data/Single cell portal/Mouse/Visp/barcodes.tsv")

#SCP5_metadata <- read.csv("C:/Users/loren/single cell/Data/Single cell portal/Mouse/2/Book1.csv")

#SCP5_metadata <- SCP5_metadata[-1,]
#SCP5_metadata <- SCP5_metadata[-1,]
colnames(SCP5_features)[2] <- "gene_short_name" 
colnames(SCP5) <- SCP5_barcodes$V1
row.names(SCP5) <- SCP5_features$gene_short_name

SCP5CDS <- new_cell_data_set(SCP5)

SCP5CDS <- detect_genes(SCP5CDS)
SCP5CDS <- estimate_size_factors(SCP5CDS)
SCP5CDS <- preprocess_cds(SCP5CDS, method = "LSI")

SCP5CDS <- reduce_dimension(SCP5CDS, preprocess_method = "LSI")
SCP5CDS <- cluster_cells(SCP5CDS, resolution = 2.5e-3)

rowData(SCP5CDS)$gene_short_name <- rownames(rowData(SCP5CDS))
#colData(SCP5CDS)$class <- SCP5_metadata$Column5
plot_cells(SCP5CDS)
plot_cells(SCP5CDS, genes = "Gad1")
#plot_cells(SCP5CDS, color_cells_by = "class")


saveRDS(SCP5CDS, "C:/Users/loren/single cell/Processed CDS/SCP5")


#####################


SCP6 <- readMM("C:/Users/loren/single cell/Data/Single cell portal/Mouse/MOP/gene_sorted-matrix.mtx")
SCP6_features <- read.table("C:/Users/loren/single cell/Data/Single cell portal/Mouse/MOP/genes.tsv")
SCP6_barcodes <- read.table("C:/Users/loren/single cell/Data/Single cell portal/Mouse/MOP/barcodes.tsv")

#SCP6_metadata <- read.csv("C:/Users/loren/single cell/Data/Single cell portal/Mouse/2/Book1.csv")

#SCP6_metadata <- SCP5_metadata[-1,]
#SCP6_metadata <- SCP5_metadata[-1,]
colnames(SCP6_features)[2] <- "gene_short_name" 
colnames(SCP6) <- SCP6_barcodes$V1
row.names(SCP6) <- SCP6_features$gene_short_name

SCP6CDS <- new_cell_data_set(SCP6)

SCP6CDS <- detect_genes(SCP6CDS)
SCP6CDS <- estimate_size_factors(SCP6CDS)
SCP6CDS <- preprocess_cds(SCP6CDS, method = "LSI")

SCP6CDS <- reduce_dimension(SCP6CDS, preprocess_method = "LSI")
SCP6CDS <- cluster_cells(SCP6CDS, resolution = 2.5e-3)

rowData(SCP6CDS)$gene_short_name <- rownames(rowData(SCP6CDS))
#colData(SCP6CDS)$class <- SCP6_metadata$Column5
plot_cells(SCP6CDS)
plot_cells(SCP6CDS, genes = "Gad2")
#plot_cells(SCP6CDS, color_cells_by = "class")

saveRDS(SCP6CDS, "C:/Users/loren/single cell/Processed CDS/SCP6")


##########################################


indata <- readMM("C:/Users/loren/single cell/Data/Snare/GSE126074_AdBrainCortex_SNAREseq_cDNA.counts.mtx") 
cellinfo <- read.table("C:/Users/loren/single cell/Data/Snare/GSE126074_AdBrainCortex_SNAREseq_cDNA.barcodes.tsv")
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"
peakinfo <- read.table("C:/Users/loren/single cell/Data/Snare/GSE126074_AdBrainCortex_SNAREseq_cDNA.genes.tsv")
names(peakinfo) <- "gene_short_name"
row.names(peakinfo) <- peakinfo$gene_short_name
rownames(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)


cds_rna <- suppressWarnings(new_cell_data_set(indata, cell_metadata = cellinfo, gene_metadata = peakinfo))
cds_rna <- detect_genes(cds_rna)
cds_rna <- estimate_size_factors(cds_rna)
cds_rna <- preprocess_cds(cds_rna, method = "LSI")

cds_rna <- reduce_dimension(cds_rna, reduction_method = 'UMAP', 
                            preprocess_method = "LSI")
cds_rna = cluster_cells(cds_rna, resolution=1e-3)

plot_cells(cds_rna,reduction_method = 'UMAP')
plot_cells(cds_rna,genes = "Gad1")


saveRDS(cds_rna, "C:/Users/loren/single cell/Processed CDS/snare_rna")

###################################