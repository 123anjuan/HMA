args <- commandArgs(T)

library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggplot2)
library(cowplot)

type1_atac_smooth<-read.csv(args[1],row.names=1) #type1 dorc ATAC smooth
type1_rna_smooth<-read.csv("args[2],row.names=1) #type1 dorc RNA smooth
type1_atac_rna_smooth<-type1_atac_smooth-type1_rna_smooth

type2_atac_smooth<-read.csv(args[3],row.names=1) #type2 dorc ATAC smooth
type2_rna_smooth<-read.csv(args[4],row.names=1) #type2 dorc RNA smooth
type2_atac_rna_smooth<-type2_atac_smooth-type2_rna_smooth

gene<-c("JUND", "JUN", "FOS", "JUNB", "RUNX1", "EGR1", "FOSL2", "STAT3")
gene<-intersect(gene,row.names(type1_atac_smooth))
type1_atac_smooth_sub<-type1_atac_smooth[gene,]
type1_rna_smooth_sub<-type1_rna_smooth[gene,]
type1_atac_rna_smooth_sub<-type1_atac_smooth_sub-type1_rna_smooth_sub

type2_atac_smooth_sub<-type2_atac_smooth[gene,]
type2_rna_smooth_sub<-type2_rna_smooth[gene,]
type2_atac_rna_smooth_sub<-type2_atac_smooth_sub-type2_rna_smooth_sub

hmcols=colorRampPalette(brewer.pal(n = 7, name ="YlGnBu"))(100)
bks <- seq(-0.1, 1.1, length.out = length(hmcols))
p1<-pheatmap::pheatmap(type1_atac_smooth_sub,
                      cluster_rows = F,
                     #cutree_rows=k_means, 
                     breaks = bks,
                     color = hmcols,
                     cluster_cols = F,
                      border_color=NA,
                         show_rownames =T,
                   # annotation_row=split_matrix
                     show_colnames = F)

hmcols=colorRampPalette(brewer.pal(n = 9, name ="RdPu"))(100)
bks <- seq(-0.1, 1.1, length.out = length(hmcols))
p2<-pheatmap::pheatmap(type1_rna_smooth_sub,
                      cluster_rows = F,
                     #cutree_rows=k_means, 
                     breaks = bks,
                     color = hmcols,
                     cluster_cols = F,
                      border_color=NA,
                         show_rownames =T,
                   # annotation_row=split_matrix
                     show_colnames = F)

hmcols=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlGn")))(100)
bks <- seq(-1.1, 1.1, length.out = length(hmcols))
p3<-pheatmap::pheatmap(type1_atac_rna_smooth_sub,
                      cluster_rows = F,
                     #cutree_rows=k_means, 
                     breaks = bks,
                     color = hmcols,
                     cluster_cols = F,
                      border_color=NA,
                         show_rownames =T,
                   # annotation_row=split_matrix
                     show_colnames = F)

hmcols=colorRampPalette(brewer.pal(n = 7, name ="YlGnBu"))(100)
bks <- seq(-0.1, 1.1, length.out = length(hmcols))
p4<-pheatmap::pheatmap(type2_atac_smooth_sub,
                      cluster_rows = F,
                     #cutree_rows=k_means, 
                     breaks = bks,
                     color = hmcols,
                     cluster_cols = F,
                      border_color=NA,
                         show_rownames =T,
                   # annotation_row=split_matrix
                     show_colnames = F)

hmcols=colorRampPalette(brewer.pal(n = 9, name ="RdPu"))(100)
bks <- seq(-0.1, 1.1, length.out = length(hmcols))
p5<-pheatmap::pheatmap(type2_rna_smooth_sub,
                      cluster_rows = F,
                     #cutree_rows=k_means, 
                     breaks = bks,
                     color = hmcols,
                     cluster_cols = F,
                      border_color=NA,
                         show_rownames =T,
                   # annotation_row=split_matrix
                     show_colnames = F)

hmcols=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlGn")))(100)
bks <- seq(-1.1, 1.1, length.out = length(hmcols))
p6<-pheatmap::pheatmap(type2_atac_rna_smooth_sub,
                      cluster_rows = F,
                     #cutree_rows=k_means, 
                     breaks = bks,
                     color = hmcols,
                     cluster_cols = F,
                      border_color=NA,
                         show_rownames =T,
                   # annotation_row=split_matrix
                     show_colnames = F)
