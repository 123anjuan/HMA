args <- commandArgs(T)

library(CellChat)
library(Seurat)
library(ggplot2)
library(ggridges)
library(dplyr)
library(igraph)
library(data.table)
library(patchwork)
library(viridisLite)
library(viridis)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(showtext)
library(future)
library(NMF)
library(ComplexHeatmap)
library(forcats)
options(stringsAsFactors = FALSE)

rds = readRDS(args[1]) 
ide = as.character(unique(rds$anno_0629))

# remove unwanted cell types
ide = setdiff(ide, c('Specialized MF','Erythrocyte','Mast cell'))
rds = subset(rds, idents = ide)

# CellChat
CellChatDB <- CellChatDB.human 
rds_young = subset(rds, cells = WhichCells(rds, expression = `age_pop` == 'young_pop'))
rds_old = subset(rds, cells = WhichCells(rds, expression = `age_pop` == 'old_pop'))
cellchat_young <- createCellChat(object = rds_young, assay = 'RNA')
cellchat_young@DB <- CellChatDB
cellchat_young <- subsetData(cellchat_young)
cellchat_young <- identifyOverExpressedGenes(cellchat_young)
cellchat_young <- identifyOverExpressedInteractions(cellchat_young)
cellchat_young <- projectData(cellchat_young, PPI.mouse)
cellchat_young <- computeCommunProb(cellchat_young)
cellchat_young <- computeCommunProbPathway(cellchat_young)
cellchat_young <- aggregateNet(cellchat_young)
cellchat_young <- netAnalysis_computeCentrality(cellchat_young)
cellchat_old <- createCellChat(object = rds_old, assay = 'RNA',  group.by = "anno_0629")
cellchat_old@DB <- CellChatDB
cellchat_old <- subsetData(cellchat_old)
cellchat_old <- identifyOverExpressedGenes(cellchat_old)
cellchat_old <- identifyOverExpressedInteractions(cellchat_old)
cellchat_old <- projectData(cellchat_old, PPI.mouse)
cellchat_old <- computeCommunProb(cellchat_old)
cellchat_old <- computeCommunProbPathway(cellchat_old)
cellchat_old <- aggregateNet(cellchat_old)
cellchat_old <- netAnalysis_computeCentrality(cellchat_old)

object.list <- list(Young = cellchat_young, Old = cellchat_old)
color_key = c("#4A308E","#FFEC8B","#66CDAA","#DB3931","#AD5551","#D309A8",
	      	"#117FA5","#437C06","#974A99","#FFA54F","#F08080","#291D99")
names(color_key) = c('Adipocyte', 'EC', 'FAP', 'Lymphocyte', 'MuSC',
		      'Myeloid cell', 'Pericyte', 'Schwann cell', 'SMC', 'Tenocyte', 'Type I', 'Type II')

age_group_color <- c(Young = "#50BA47", Old = "#A15DBA")
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# count number of interactions
reform_list <- subsetCommunication(cellchat,slot.name = "net",signaling = NULL)
complex_ref <- readRDS('CellChatDB.human.complex.gene_list.RDS')
inter_df = list()
for(n in names(reform_list)){
	    reform_list[[n]]$age_group = n
    tmp = merge(reform_list[[n]], complex_ref, by.x = 3, by.y = 2, all.x = T)
        inter_df[[n]] = merge(tmp, complex_ref, by.x = 4, by.y = 2, all.x = T)[,c(3,4,2,13,1,14,5:12)]
        inter_df[[n]]$gene_name.x[is.na(inter_df[[n]]$gene_name.x)] = inter_df[[n]]$ligand[is.na(inter_df[[n]]$gene_name.x)]
	    inter_df[[n]]$gene_name.y[is.na(inter_df[[n]]$gene_name.y)] = inter_df[[n]]$receptor[is.na(inter_df[[n]]$gene_name.y)]
	    inter_df[[n]]$ID = paste0(inter_df[[n]]$source,"->",inter_df[[n]]$target,"&",inter_df[[n]]$interaction_name,"&",inter_df[[n]]$gene_name.x,"_",inter_df[[n]]$gene_name.y)     
	        colnames(inter_df[[n]])[c(4,6)] = c("gn1", "gn2")
	        inter_df[[n]] =inter_df[[n]][,c(15,1:14)]
}

Young =     rowSums(data.matrix(table(inter_df[["Young"]]$source, inter_df[["Young"]]$target)))
Old =       rowSums(data.matrix(table(inter_df[["Old"]]$source, inter_df[["Old"]]$target)))
df_y = data.frame(Young)
df_y$cell_type = rownames(df_y)
rownames(df_y) = NULL
df_o = data.frame(Old)
df_o$cell_type = rownames(df_o)
rownames(df_o) = NULL

CCI_count_all = merge(df_y,df_o,by="cell_type",all=T)
rownames(CCI_count_all) = CCI_count_all$cell_type
CCI_count_all = CCI_count_all[,c(2:3,1)]
CCI_count_all[CCI_count_all$Young=="0", "Young"] = NA
CCI_count_all[CCI_count_all$Old=="0", "Old"] = NA
CCI_count_all = CCI_count_all[complete.cases(CCI_count_all),]
data_all_Y <- select(CCI_count_all,Young,cell_type)
data_all_Y$age_group <- "Young"
colnames(data_all_Y) <- c("Number","cell_type","age_group")
data_all_O <- select(CCI_count_all,Old,cell_type)
data_all_O$age_group <- "Old"
colnames(data_all_O) <- c("Number","cell_type","age_group")
df_mer = rbind(data_all_Y,data_all_O)
celltype_all = unique(CCI_count_all$cell_type)
