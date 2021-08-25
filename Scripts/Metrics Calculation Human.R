
library(monocle3)
library(Seurat)
library(SeuratObject)
library(SeuratWrappers)
library(dplyr)
library(Matrix)
library(aricode)

############################# Functions Definitions ###############

specificity_matrix <- function(agg_expr_matrix, cores=1){
  specificity_mat <-
    pbmcapply::pbmclapply(row.names(agg_expr_matrix),
                          FUN = function(x) {
                            agg_exprs = as.numeric(agg_expr_matrix[x,])
                            agg_exprs = makeprobsvec(agg_exprs)
                            perfect_spec_matrix = diag(ncol(agg_expr_matrix))
                            sapply(1:ncol(agg_expr_matrix), function(col_idx) {
                              1 - JSdistVec(agg_exprs,
                                            perfect_spec_matrix[,col_idx])
                            })
                          }, mc.cores=cores,
                          ignore.interactive = TRUE)
  specificity_mat = do.call(rbind, specificity_mat)
  colnames(specificity_mat) = colnames(agg_expr_matrix)
  row.names(specificity_mat) = row.names(agg_expr_matrix)
  return(specificity_mat)
  #
}

JSdistVec <- function (p, q)
{
  JSdiv <- shannon.entropy((p + q)/2) - (shannon.entropy(p) +
                                           shannon.entropy(q)) * 0.5
  JSdiv[is.infinite(JSdiv)] <- 1
  JSdiv[JSdiv < 0] <- 0
  JSdist <- sqrt(JSdiv)
  JSdist
}

makeprobsvec <- function(p) {
  phat <- p/sum(p)
  phat[is.na(phat)] = 0
  phat
}

makeprobs <- function(a) {
  colSums<-apply(a,2,sum)
  b <- Matrix::t(Matrix::t(a)/colSums)
  b[is.na(b)] = 0
  b
}

shannon.entropy <- function(p) {
  if (min(p) < 0 || sum(p) <=0)
    return(Inf)
  p.norm <- p[p>0]/sum(p)
  -sum(log2(p.norm)*p.norm)
}

gini.index <- function(p){
  l <- length(p)
  s <- summary(p)
  s <- (s/l)^2
  G <- 1-sum(s)
  return(G)
}

score_table <- function(CDS, cell_group_df)
{
  cluster_binary_exprs = as.matrix(aggregate_gene_expression(CDS, cell_group_df = cell_group_df, norm_method = "binary"))
  
  cluster_mean_exprs = as.matrix(aggregate_gene_expression(CDS, 
                                                           cell_group_df = cell_group_df, norm_method = "size_only"))
  
  cluster_mean_exprs_value <- cluster_mean_exprs/cluster_binary_exprs
  cluster_mean_exprs_value[cluster_mean_exprs_value == 'NaN'] <- 0
  
  cluster_spec_mat = specificity_matrix(cluster_mean_exprs)
  cluster_marker_score_mat = as.matrix(cluster_binary_exprs * 
                                         cluster_spec_mat)
  
  
  cluster_marker_score_table = tibble::rownames_to_column(as.data.frame(cluster_marker_score_mat))
  cluster_marker_score_table = tidyr::gather(cluster_marker_score_table, "cell_group", "marker_score", -rowname)
  
  cluster_spec_table = tibble::rownames_to_column(as.data.frame(cluster_spec_mat))
  cluster_spec_table = tidyr::gather(cluster_spec_table, "cell_group", "specificity", -rowname)
  
  cluster_expr_table = tibble::rownames_to_column(as.data.frame(cluster_mean_exprs))
  cluster_expr_table = tidyr::gather(cluster_expr_table, "cell_group", "mean_expression", -rowname)
  
  cluster_expr_val_table = tibble::rownames_to_column(as.data.frame(cluster_mean_exprs_value))
  cluster_expr_val_table = tidyr::gather(cluster_expr_val_table, "cell_group", "mean_expression_value", -rowname)
  
  cluster_fraction_expressing_table = tibble::rownames_to_column(as.data.frame(cluster_binary_exprs))
  cluster_fraction_expressing_table = tidyr::gather(cluster_fraction_expressing_table, "cell_group", "fraction_expressing", -rowname)
  
  
  cluster_marker_score_table$specificity = cluster_spec_table$specificity
  cluster_marker_score_table$mean_expression = cluster_expr_table$mean_expression
  cluster_marker_score_table$fraction_expressing = cluster_fraction_expressing_table$fraction_expressing
  cluster_marker_score_table$mean_expression_value = cluster_expr_val_table$mean_expression_value
  
  return(cluster_marker_score_table)
  
}


################



############# One needs to load a preprocessed CDS of a HUMAN dataset
human_CDS <- readRDS("C:/Users/loren/single cell/Processed CDS/allen_human_ACC")
human_CDS@colData@listData[["clusters"]] <- human_CDS@clusters@listData[["UMAP"]][["clusters"]]
human_CDS@colData@listData[["partition"]] <- human_CDS@clusters@listData[["UMAP"]][["partitions"]]


############# one needs to select the interneurons cells
human_CDS_sub <- choose_cells(human_CDS)




#GG <- human_CDS[row.names(subset(rowData(human_CDS),gene_short_name %in% c("GAD1", "GAD2"))),]
#GG <- cluster_cells(GG, resolution = 1e-4)
#plot_cells(GG, color_cells_by = "partition")
#plot_cells(human_CDS, genes = "SLC17A7")
#GG@colData@listData[["clusters"]] <- GG@clusters@listData[["UMAP"]][["clusters"]]
#GG@colData@listData[["partition"]] <- GG@clusters@listData[["UMAP"]][["partitions"]]

###############
###############   FIRST LEVEL - GAD1, GAD2  ################
###############

GG <- human_CDS[row.names(subset(rowData(human_CDS),gene_short_name %in% c("GAD1", "GAD2"))),]

plot_percent_cells_positive(GG, group_cells_by = "clusters")

Num_expr <- GG@rowRanges@elementMetadata@listData[["num_cells_expressed"]]


############# here one needs to manually set all the four cluster derived from the clustering processing
############# as one cluster, beacuse we want to calculate metrics on all the interneurons as a whole
############# if one emploies the processed dataset provided, one can look at the file "Processed datasets cluster numbers" to know the clusters numbers.
GG@clusters@listData[["UMAP"]][["clusters"]] <- recode_factor(GG@clusters@listData[["UMAP"]][["clusters"]], 
                                                              "20" = "4",
                                                              "7" = "4",
                                                              "10" = "4",
                                                              "1" = "4"
)


cell_group_df <- data.frame(row.names = row.names(colData(GG)), 
                            cell_id = row.names(colData(GG)))
cell_group_df$cell_group <- tryCatch({  clusters(GG, reduction_method = "UMAP")}, error = function(e) {
  NULL
})

cluster_marker_score_table <- score_table(GG, cell_group_df)

Geni <- agg_mat@Dimnames[[1]]
association <- as.data.frame(Geni)
############### here one needs to manually insert the cluster name set before, and associates the markers to the respective cluster 
association$cluster <- c(4,4)
Gain<-0
Impurity_Reduction<-0


################ here it calculates the IG and IR for all the markers

for (x in 1:2) {
  
  colname <- agg_mat@Dimnames[[1]][x]
  cell_group_df[colname] <- agg_mat[x,]
  
  clus_numb <- association$cluster[association$Geni == colname]
  
  ############binarize the class in two categories the ones of the cluster we want to study and all the others together
  binary_class <- cell_group_df
  binary_class$cell_group[binary_class$cell_group %in% clus_numb ] <- clus_numb
  binary_class$cell_group[binary_class$cell_group != clus_numb ] <- 2
  
  ############### divide the cells if they express or not the marker  
  Yes <- subset(binary_class, binary_class[colname]==1)
  No <- subset(binary_class, binary_class[colname]==0)
  
  NY <- dim(Yes)[1]
  NN <- dim(No)[1]
  tot <- dim(cell_group_df)[1]
 
  ######### calculate all the gini indexes and consequentially the impurity reduction
  Gini <- gini.index(binary_class$cell_group)  
  GiniYes <- gini.index(Yes$cell_group)
  GiniNo <- gini.index(No$cell_group)
  Impurity_Reduction[x] <- Gini - NY*GiniYes/tot -NN*GiniNo/tot
  
  ######### calculate all the entropy and consequentially the information gain
  ENTROPY <- entropy(binary_class$cell_group,binary_class[[colname]])$U
  EntrYes <- entropy(Yes$cell_group, as.vector(Yes[[colname]]))$U
  EntrNo <- entropy(No$cell_group, as.vector(No[[colname]]))$U
  Gain[x] <- ENTROPY - EntrYes*NY/tot - EntrNo*NN/tot
}

association$gain <- Gain
association$Imp_red <- Impurity_Reduction

Gad_HUMAN <- association
Gad_scores_HUMAN <- cluster_marker_score_table
Gad_scores_HUMAN_spec <- Gad_scores_HUMAN %>%
  top_n(2, marker_score)

##################
################### SECOND LEVEL - ADARB2 LHX6 HTR3A
###################


ALH <- human_CDS_sub[row.names(subset(rowData(human_CDS_sub),gene_short_name %in% c("ADARB2", "HTR3A", "LHX6"))),]

############# again here one needs to manually set all the four clusters derived from the clustering processing
############# as two clusters, beacuse we want to calculate metrics on MGE and CGE derived cells respectively

ALH@clusters@listData[["UMAP"]][["clusters"]] <-recode_factor(ALH@clusters@listData[["UMAP"]][["clusters"]], 
                                                              "7" = "4",
                                                              "20" = "4",
                                                              "10" = "1"
)




cell_group_df <- data.frame(row.names = row.names(colData(ALH)), 
                            cell_id = row.names(colData(ALH)))
cell_group_df$cell_group <- tryCatch({  clusters(ALH, reduction_method = "UMAP")}, error = function(e) {
  NULL
})

cluster_marker_score_table <- score_table(ALH, cell_group_df)

agg_mat <- normalized_counts(ALH, norm_method="binary",
                             pseudocount=1)

Geni <- agg_mat@Dimnames[[1]]
association <- as.data.frame(Geni)
association$cluster <- c(1,1,4)


for (x in 1:3) {
  
  colname <- agg_mat@Dimnames[[1]][x]
  cell_group_df[colname] <- agg_mat[x,]
  
  clus_numb <- association$cluster[association$Geni == colname]
  
  binary_class <- cell_group_df
  binary_class$cell_group[binary_class$cell_group %in% clus_numb ] <- clus_numb
  binary_class$cell_group[binary_class$cell_group != clus_numb ] <- 2
  
  Yes <- subset(binary_class, binary_class[colname]==1)
  No <- subset(binary_class, binary_class[colname]==0)
  
  NY <- dim(Yes)[1]
  NN <- dim(No)[1]
  tot <- dim(cell_group_df)[1]
  
  Gini <- gini.index(binary_class$cell_group)  
  GiniYes <- gini.index(Yes$cell_group)
  GiniNo <- gini.index(No$cell_group)
  Impurity_Reduction[x] <- Gini - NY*GiniYes/tot -NN*GiniNo/tot
  
  
  ENTROPY <- entropy(binary_class$cell_group,binary_class[[colname]])$U
  EntrYes <- entropy(Yes$cell_group, as.vector(Yes[[colname]]))$U
  EntrNo <- entropy(No$cell_group, as.vector(No[[colname]]))$U
  Gain[x] <- ENTROPY - EntrYes*NY/tot - EntrNo*NN/tot
}

association$gain <- Gain
association$Imp_red <- Impurity_Reduction

MGE_HUMAN <- association
MGE_score_HUMAN <- cluster_marker_score_table
MGE_score_spec_HUMAN <- cluster_marker_score_table
MGE_score_spec_HUMAN <- MGE_score_spec_HUMAN %>%
  top_n(3, marker_score)


##################
################### THIRD LEVEL - PVALB SST VIP LAMP5
###################


PSVL <- human_CDS_sub[row.names(subset(rowData(human_CDS_sub),gene_short_name %in% c("PVALB", "SST", "VIP", "LAMP5"))),]

PSVL@clusters@listData[["UMAP"]][["clusters"]] <-recode_factor(PSVL@clusters@listData[["UMAP"]][["clusters"]], 
                                                              "20"="7"
)
plot_cells(PSVL, genes = "LAMP5")

cell_group_df <- data.frame(row.names = row.names(colData(PSVL)), 
                            cell_id = row.names(colData(PSVL)))
cell_group_df$cell_group <- tryCatch({  clusters(PSVL, reduction_method = "UMAP")}, error = function(e) {
  NULL
})

plot_percent_cells_positive(PSVL , group_cells_by = "clusters")

cluster_marker_score_table <- score_table(PSVL, cell_group_df)


agg_mat <- normalized_counts(PSVL, norm_method="binary",
                             pseudocount=1)

Geni <- agg_mat@Dimnames[[1]]
association <- as.data.frame(Geni)
association$cluster <- c(10,7,4,1)

for (x in 1:4) {
  
  colname <- agg_mat@Dimnames[[1]][x]
  cell_group_df[colname] <- agg_mat[x,]
  
  clus_numb <- association$cluster[association$Geni == colname]
  
  binary_class <- cell_group_df
  binary_class$cell_group[binary_class$cell_group %in% clus_numb ] <- clus_numb
  binary_class$cell_group[binary_class$cell_group != clus_numb ] <- 2
  
  Yes <- subset(binary_class, binary_class[colname]==1)
  No <- subset(binary_class, binary_class[colname]==0)
  
  NY <- dim(Yes)[1]
  NN <- dim(No)[1]
  tot <- dim(cell_group_df)[1]
  
  
  Gini <- gini.index(binary_class$cell_group)  
  GiniYes <- gini.index(Yes$cell_group)
  GiniNo <- gini.index(No$cell_group)
  Impurity_Reduction[x] <- Gini - NY*GiniYes/tot -NN*GiniNo/tot
  
  ENTROPY <- entropy(binary_class$cell_group,binary_class[[colname]])$U
  EntrYes <- entropy(Yes$cell_group, as.vector(Yes[[colname]]))$U
  EntrNo <- entropy(No$cell_group, as.vector(No[[colname]]))$U
  Gain[x] <- ENTROPY - EntrYes*NY/tot - EntrNo*NN/tot
}

association$gain <- Gain
association$Imp_red <- Impurity_Reduction


PSVL_HUMAN <- association
PSVL_score_HUMAN <- cluster_marker_score_table 
PSVL_score_HUMAN <- PSVL_score_HUMAN %>%
  top_n(4, marker_score)
############################


############# put all together in one table
INF_HUMAN = rbind(GAD_HUMAN, MGE_HUMAN, PSVL_HUMAN)
score_HUMAN <- rbind(Gad_scores_HUMAN_spec, MGE_score_spec_HUMAN, PSVL_score_HUMAN)


final_scores_HUMAN <- merge.data.frame(INF_HUMAN, score_HUMAN, by.x = c("Geni", "dataset", "cluster"), by.y = c("rowname","dataset", "cell_group"))



############### plotting

ggplot(final_scores_HUMAN, aes(x = specificity , y = fraction_expressing)) + geom_point(aes(color = dataset, size = mean_expression)) + geom_text(aes(label = Geni),hjust=-0.2, vjust=-0.2, size = 5) +
  labs(x = "Specificity", y = "Fraction Expressing", size = "Mean Expression", color = "Dataset", title = "Specificity-Fraction Plot", legend.text=element_text(size=13)) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5), legend.title=element_text(size=18), legend.text=element_text(size=16), axis.title = element_text(size=18), axis.text = element_text(size=13))


ggplot(final_scores_HUMAN, aes(x = gain , y = Imp_red)) + geom_point(aes(color = dataset)) + geom_text(aes(label = Geni),hjust=-0.2, vjust=-0.2, size = 5) +
  labs(x = "Information Gain", y = "Impurity Reduction", color = "Dataset", title = "IG-IR Plot", legend.text=element_text(size=13)) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5), legend.title=element_text(size=18), legend.text=element_text(size=16), axis.title = element_text(size=18), axis.text = element_text(size=13))


ggplot(final_scores_HUMAN, aes(x =  marker_score , y = gain)) + geom_point(aes(color = dataset)) + geom_text(aes(label = Geni),hjust=-0.2, vjust=-0.2, size = 5) +
  labs(x = "Marker Score", y = "Information Gain", color = "Dataset", title = "Marker Score-IG Plot", legend.text=element_text(size=13)) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5), legend.title=element_text(size=18), legend.text=element_text(size=16), axis.title = element_text(size=18), axis.text = element_text(size=13))


ggplot(final_scores_HUMAN, aes(x =  marker_score , y = log(mean_expression))) + geom_point(aes(color = dataset)) + geom_text(aes(label = Geni),hjust=-0.2, vjust=-0.2, size = 5) +
  labs(x = "Marker Score", y = "Mean Expression", color = "Dataset", title = "Marker Score-Mean Expression Plot", legend.text=element_text(size=13)) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5), legend.title=element_text(size=18), legend.text=element_text(size=16), axis.title = element_text(size=18), axis.text = element_text(size=13))

