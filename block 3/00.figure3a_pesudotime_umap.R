args <- commandArgs(T)

library(Seurat)
library(ggplot2)
library(dplyr)
library(monocle3)
library(patchwork)
library(sctransform)

rds <- readRDS(args[1]) #myonucleus snRNA-seq rds
type1 <- rds[,rds$anno_0713 %in% c("Type I","Type I-Aging","Type I-NMJ 1","Type I-NMJ 2")]
type1 <- type1[,type1$integrated_snn_res.2 %in% c(2,4,7,8,10,12,13,14,15,16,17,18,19,21,26,29)]
type2 <- rds[,rds$anno_0713 %in% c("Type II","Type II-Aging","Type II-ENOX1+")]
type2 <- type2[,type2$integrated_snn_res.2 %in% c(0,1,3,5,6,11,22,23,24)]

data_t1 <- GetAssayData(type1, assay = 'RNA', slot = 'counts')
cell_metadata_t1 <- type1@meta.data
gene_annotation_t1 <- data.frame(gene_short_name = rownames(data_t1))
rownames(gene_annotation_t1) <- rownames(data_t1)
cds_t1 <- new_cell_data_set(data_t1,
                         cell_metadata = cell_metadata_t1,
                         gene_metadata = gene_annotation_t1)
cds_t1 <- preprocess_cds(cds_t1, num_dim = 50)
cds_t1 <- reduce_dimension(cds_t1,preprocess_method = "PCA")
cds.embed_t1 <- cds_t1@int_colData$reducedDims$UMAP
int.embed_t1 <- Embeddings(type1, reduction = "umap")
int.embed_t1 <- int.embed_t1[rownames(cds.embed_t1),]
cds_t1@int_colData$reducedDims$UMAP <- int.embed_t1
cds_t1 <- cluster_cells(cds = cds_t1, reduction_method = "UMAP")
cds_t1 <- learn_graph(cds_t1, use_partition = TRUE)

get_earliest_principal_node <- function(cds, time_bin=c(53)){
  cell_ids <- which(colData(cds)[, "integrated_snn_res.7"] == time_bin)
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds_t1 <- order_cells(cds_t1, root_pr_nodes=get_earliest_principal_node(cds_t1))
p<-plot_cells(cds_t1,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=0.3,
           trajectory_graph_segment_size=0.5)+NoLegend()

pseudotime <- as.data.frame(cds_t1@principal_graph_aux$UMAP$pseudotime)
colnames(pseudotime) = 'pseu_type1'
pseudotime$rank_type1=cut_number(pseudotime$pseu_type1,n=100,labels=F)

data_t2 <- GetAssayData(type2, assay = 'RNA', slot = 'counts')
cell_metadata_t2 <- type2@meta.data
gene_annotation_t2 <- data.frame(gene_short_name = rownames(data_t2))
rownames(gene_annotation_t2) <- rownames(data_t2)
cds_t2 <- new_cell_data_set(data_t2,
                         cell_metadata = cell_metadata_t2,
                         gene_metadata = gene_annotation_t2)
cds_t2 <- preprocess_cds(cds_t2, num_dim = 50)
cds_t2 <- reduce_dimension(cds_t2,preprocess_method = "PCA")
cds.embed_t2 <- cds_t2@int_colData$reducedDims$UMAP
int.embed_t2 <- Embeddings(type2, reduction = "umap")
int.embed_t2 <- int.embed_t2[rownames(cds.embed_t2),]
cds_t2@int_colData$reducedDims$UMAP <- int.embed_t2
cds_t2 <- cluster_cells(cds_t2)
cds_t2 <- learn_graph(cds_t2)

get_earliest_principal_node <- function(cds, time_bin=13){
  cell_ids <- which(colData(cds)[, "integrated_snn_res.7"] == time_bin)
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds_t2_two <- order_cells(cds_t2, root_pr_nodes=get_earliest_principal_node(cds_t2))
p<-plot_cells(cds_t2_two,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=0.3,
           trajectory_graph_segment_size=0.5)+NoLegend()

pseudotime <- as.data.frame(cds_t2@principal_graph_aux$UMAP$pseudotime)
colnames(pseudotime) = 'pseu_type2'
pseudotime$rank_type2=cut_number(pseudotime$pseu_type2,n=100,labels=F)

rds$pseu_type1 <- type1$pseu_type1
rds$rank_type1 <- type1$rank_type1
rds$pseu_type2 <- type2$pseu_type2
rds$rank_type2 <- type2$rank_type2
rds$rank <- rds$rank_type1
rds$pseu <- rds$pseu_type1
rds$pseu[match(rownames(type2@meta.data), rownames(rds@meta.data))] <- type2$pseu_type2
rds$rank[match(rownames(type2@meta.data), rownames(rds@meta.data))] <- type2$rank_type2

p <- FeaturePlot(rds,features = "pseu", order = T) + scale_colour_gradientn(colours = rev(viridis(7)))
