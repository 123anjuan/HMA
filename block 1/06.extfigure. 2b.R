#correlation
args <- commandArgs(T)
library(psych)
library(ggplot2)
library(Seurat)

seRNA = readRDS(args[1])
DefaultAssay(seRNA) = 'RNA'
Idents(seRNA) = seRNA$celltype
seRNA
seRNA <- FindVariableFeatures(seRNA)
genesUse <- VariableFeatures(object = seRNA)

proj =  loadArchRProject(showLogo = F,path = 'ATAC')#ATAC ArchRobject
gene_score = getMatrixFromProject(proj,useMatrix = 'GeneScoreMatrix')
GS = gene_score$GeneScoreMatrix
GS = assays(gene_score)$GeneScoreMatrix
rownames(GS) = rowData(gene_score)$name
gene_score <- GS[rownames(GS) %in% genesUse,]
mat <- log(gene_score + 1)

snATAC <- CreateSeuratObject(
                counts = mat,
                assay = 'GeneScore',
                project = 'ATAC',
                min.cells = 1,
                meta.data = as.data.frame(proj@cellColData)
        )
snATAC <- ScaleData(snATAC, verbose = FALSE)

Idents(snATAC) = snATAC$celltype
averages_ATAC<-AverageExpression(snATAC)$GeneScore
averages<-AverageExpression(seRNA)$RNA
averages_1 = averages[rownames(averages_ATAC),]
averages_1 = averages_1[,colnames(averages_ATAC)]
averages_ATAC = averages_ATAC[,order(colnames(averages_ATAC))]
averages_1 = averages_1[,colnames(averages_ATAC)]

V1_Scaled <- t(scale(t(averages_ATAC), center = T, scale=T))
V2_Scaled <- t(scale(t(averages_1), center = T, scale=T))
cortest_psy <- corr.test(V1_Scaled, V2_Scaled, method = "pearson",adjust = "fdr")
p = pheatmap(cortest_psy$r, cluster_rows = F, border_color = NA,cluster_cols = F, show_rownames = T, show_colnames = T,color = colorRampPalette(c("#999999","white","#ef8a62"))(100),angle_col = 90)

