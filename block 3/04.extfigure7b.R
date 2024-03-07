args <- commandArgs(T)

library(doParallel)
library(BuenColors)
library(FigR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(ComplexHeatmap)
library(networkD3)
require(htmlwidgets)
library(cowplot)
library(ggrepel)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(viridis)
get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  
get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
######################### type1
atac_sub_cor<-read.csv(args[1],row.names=1) #type1 dorc atac correlation
metadata<-meta[colnames(atac_sub_cor),]
metadata = subset(metadata,select = "rank_new_label")
#metadata$rank_type1= factor(metadata$rank_type1, levels = 1:100)
metadata$rank_new_label=factor(metadata$rank_new_label,levels = c("1_type1","2_type1","3_type1","4_type1","5_type1","6_type1","7_type1","8_type1","9_type1","10_type1"))
ann_colors = list(
  rank_new_label = c("1_type1"="#440154FF","2_type1"="#482878FF","3_type1"="#3E4A89FF","4_type1"="#31688EFF","5_type1"="#26828EFF",
"6_type1"="#1F9E89FF","7_type1"="#35B779FF","8_type1"="#6DCD59FF","9_type1"="#B4DE2CFF","10_type1"="#FDE725FF"))


bk = seq(0.1, 1,by = 0.01)
col_length = length(bk)
if_1 = inferno(col_length/2)
if_2 = rev(mako(col_length/2))
col= append(if_2,if_1)

q<-pheatmap(atac_sub_cor,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = F,
         show_colnames = F,
         annotation_colors=ann_colors,
         clustering_method="ward.D2",
         annotation_col=metadata,
         breaks=  bk,
         na_col = NA,
         legend = F,
         annotation_legend = F,
         color = col)
atac_sub_cor<-atac_sub_cor[q$tree_row$order,q$tree_col$order]
atac_sub_cor_utri = get_upper_tri(atac_sub_cor)
q1<-pheatmap(atac_sub_cor_utri,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         annotation_colors=ann_colors,
         annotation_col=metadata,
         breaks=  bk,
         na_col = NA,
         legend = F,
         annotation_legend = F,
         color = col)
atac_sub_cor_ltri = get_lower_tri(atac_sub_cor)
q2<-pheatmap(atac_sub_cor_ltri,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         annotation_colors=ann_colors,
         annotation_col=metadata,
         breaks=  bk,
         na_col = NA,
         legend = F,
         annotation_legend = F,
         color = col)

rna_sub_cor<-read.csv("./output/type1_dorc_rna_1000_cor.csv",row.names=1)
rna_sub_cor<-rna_sub_cor[q$tree_row$order,q$tree_col$order]

bk = seq(0.5, 1,by = 0.01)
col_length = length(bk)
if_1 = inferno(col_length/2)
if_2 = rev(mako(col_length/2))
col= append(if_2,if_1)

p<-pheatmap(rna_sub_cor,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         annotation_colors=ann_colors,
         #clustering_method="ward.D2",
         annotation_col=metadata,
         breaks=  bk,
        na_col = NA,
        legend = F,
        annotation_legend = F,color = col)
rna_sub_cor_ltri = get_lower_tri(rna_sub_cor)
p1<-pheatmap(rna_sub_cor_ltri,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         annotation_colors=ann_colors,
         annotation_col=metadata,
         breaks=  bk,
         na_col = NA,
         legend = F,
         annotation_legend = F,
         color = col)

rna_sub_cor_utri = get_upper_tri(rna_sub_cor)
p2<-pheatmap(rna_sub_cor_utri,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         annotation_colors=ann_colors,
         annotation_col=metadata,
         breaks=  bk,
         na_col = NA,
         legend = F,
         annotation_legend = F,
         color = col)

######################## type2
atac_sub_cor<-read.csv("type2_dorc_atac_1000_cor.csv",row.names=1)
metadata<-meta[colnames(atac_sub_cor),]
metadata = subset(metadata,select = "rank_new_label")
#metadata$rank_type1= factor(metadata$rank_type1, levels = 1:100)
metadata$rank_new_label=factor(metadata$rank_new_label,levels = c("1_type2","2_type2","3_type2","4_type2","5_type2","6_type2","7_type2","8_type2","9_type2","10_type2"))
ann_colors = list(
  rank_new_label = c("1_type2"='#000004FF',"2_type2"='#1B0C42FF',"3_type2"='#4B0C6BFF',"4_type2"='#781C6DFF',"5_type2"='#A52C60FF',
"6_type2"='#CF4446FF',"7_type2"='#ED6925FF',"8_type2"='#FB9A06FF',"9_type2"='#F7D03CFF',"10_type2"='#FCFFA4FF'))

bk = seq(0.1, 1,by = 0.01)
col_length = length(bk)
if_1 = inferno(col_length/2)
if_2 = rev(mako(col_length/2))
col= append(if_2,if_1)
q<-pheatmap(atac_sub_cor,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = F,
         show_colnames = F,
         annotation_colors=ann_colors,
         clustering_method="average",
         annotation_col=metadata,
         breaks=  bk,
         na_col = NA,
         #legend = F,
         annotation_legend = F,
         color = col)

atac_sub_cor<-atac_sub_cor[q$tree_row$order,q$tree_col$order]

atac_sub_cor_utri = get_upper_tri(atac_sub_cor)

q1<-pheatmap(atac_sub_cor_utri,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         annotation_colors=ann_colors,
         annotation_col=metadata,
         breaks=  bk,
         na_col = NA,
         #legend = F,
         annotation_legend = F,
         color = col)

 atac_sub_cor_ltri = get_lower_tri(atac_sub_cor)

 q2<-pheatmap(atac_sub_cor_ltri,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         annotation_colors=ann_colors,
         annotation_col=metadata,
         breaks=  bk,
         na_col = NA,
         legend = F,
         annotation_legend = F,
         color = col)
ggsave('type2_lowertricangle_atac.png',q2, dpi  = 600, width = 4, height = 4)

rna_sub_cor<-read.csv("type2_dorc_rna_1000_cor.csv",row.names=1)
rna_sub_cor<-rna_sub_cor[q$tree_row$order,q$tree_col$order]


bk = seq(0.5, 1,by = 0.01)
col_length = length(bk)
if_1 = inferno(col_length/2)
if_2 = rev(mako(col_length/2))
col= append(if_2,if_1)

p<-pheatmap(rna_sub_cor,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         annotation_colors=ann_colors,
         annotation_col=metadata,
         breaks=  bk,
         na_col = NA,
         #legend = F,
         annotation_legend = F,
         color = col)

rna_sub_cor_utri = get_upper_tri(rna_sub_cor)
p1<-pheatmap(rna_sub_cor_utri,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         annotation_colors=ann_colors,
         annotation_col=metadata,
         breaks=  bk,
         na_col = NA,
         legend = F,
         annotation_legend = F,
         color = col)

rna_sub_cor_ltri = get_lower_tri(rna_sub_cor)
p2<-pheatmap(rna_sub_cor_ltri,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,
         annotation_colors=ann_colors,
         annotation_col=metadata,
         breaks=  bk,
         na_col = NA,
         legend = F,
         annotation_legend = F,
         color = col)
