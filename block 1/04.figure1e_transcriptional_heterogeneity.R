args <- commandArgs(T)
library(Seurat)
library(Matrix)
library(pheatmap)
library(sfsmisc)
library(MASS)
library(hopach)
library(showtext)
showtext_auto()
library(egg)
set.seed(1)

library(ggplot2)
library(dplyr)
library(cowplot)
library(ggsci)

rds = readRDS(args[1])#snRNA rds

rds$celltype = as.character(rds$celltype)
DefaultAssay(rds) <- 'RNA'
celltypes <- unique(rds@meta.data$celltype)

#---------------
#celltype
DefaultAssay(rds) <- 'RNA'
rds@meta.data$celltype = rds@meta.data$Final_annotation
celltypes <- unique(rds@meta.data$celltype)


getEuclideanDistance <- function(type){
  print(paste("Working on", type))
  tmp <- subset(rds, celltype == as.character(type))
  expr <- tmp@assays$RNA@counts
expr = round(expr)  
print(range(expr))
  gc()
  zeros <- which(Matrix::rowSums(expr) == 0)
   expr <- expr[rownames(expr) %ni% rownames(expr)[zeros],]
Down_Sample_Matrix <-function (expr_mat) {
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
      prob <- min_lib_size/sum(x)
      return(unlist(lapply(x, function(y) {
        rbinom(1, y, prob)
      })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
  }
  ds_expr <- Down_Sample_Matrix(expr)
  nsample <- min(table(tmp@meta.data$age_group)[c("Adult", "Old")])
      if(nsample >300){
  nsample = 300
  }

      if(nsample < 10){
    print("Not enough cells")
    return(NULL)
  } 
  print(nsample)  
Adult_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$age_group == "Adult")], nsample) 
Old_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$age_group == "Old")], nsample)
ds_expr_r <- ds_expr[, c(Adult_r, Old_r)]
    getLowCVgenes <- function(matr){
      means <- Matrix::rowMeans(matr)
      bins <- quantile(means, c(seq(from = 0, to = 1, length = 11)),na.rm=TRUE)
      mean_bin <- unlist(lapply(means, function(x) min(which(bins >= x))))
      asplit <- split(names(means), mean_bin)
      genes <- unique(unlist(lapply(asplit[setdiff(names(asplit), c("1", "11"))], function(x){
        coef_var <- apply(matr, 1, function(x) sd(x)/mean(x))
        bottom10percent <- names(head(sort(coef_var), round(10*length(coef_var))))
      })))
      genes
    }
    genes <- getLowCVgenes(ds_expr_r)
    calcEuclDist <- function(matr, Adult, Old){
    tmp <- data.matrix(sqrt(matr[genes, Adult]))
    mean <- rowMeans(sqrt(matr[genes, Adult]))
    d_Adult <- distancevector(t(tmp), mean , d="euclid")
    names(d_Adult) <- Adult
    
    mean <- rowMeans(sqrt(matr[genes, Old]))
    d_Old <- distancevector(t(tmp), mean , d="euclid")
    names(d_Old) <- Old
    
    
  
    list(Adult = d_Adult,Old = d_Old)
  }
  ds <- calcEuclDist(matr = ds_expr_r, Adult = Adult_r,Old = Old_r)
  ds
}


#----------------------------------------
res <- lapply(celltypes, function(x) getEuclideanDistance(x))
names(res) <- celltypes

res_original <- res

# Calculate mean differences and p-values ####
diffs <- unlist(lapply(res_original, function(x) log2(mean(x[[2]]) / mean(x[[1]]))))
pvals <- unlist(lapply(res_original, function(x) wilcox.test(x[[1]], x[[2]])$p.value))
adj_pvals <- p.adjust(pvals, method = "BH")
sizes <- (-log10(adj_pvals))
sizes[which(sizes < 1)] <- 1
sizes[which(sizes > 4)] <- 4
sizes <- sizes * 0.75
farben <- rep("grey", length(adj_pvals))
farben[which(adj_pvals < 0.05)] <- "purple"

# Calculate mean differences and p-values ####
diffs <- unlist(lapply(res_original, function(x) log2(mean(x[[2]]) / mean(x[[1]]))))
pvals <- unlist(lapply(res_original, function(x) wilcox.test(x[[1]], x[[2]])$p.value))
adj_pvals <- p.adjust(pvals, method = "BH")
sizes <- (-log10(adj_pvals))
sizes[which(sizes < 1)] <- 1
sizes[which(sizes > 4)] <- 4
sizes <- sizes * 0.75
farben <- rep("grey", length(adj_pvals))
farben[which(adj_pvals < 0.05)] <- "purple"

ord <- order(names(res))
par(mar = c(15,5,2,5))
boxplot(do.call(c, res[ord]), las = 2, outline = F, col = c("#50BA47", "#A15DBA"), ylab = "Transcriptional noise", xaxt = 'n') +
axis(1, at = seq(from = 1.5, to = 22.5, by = 2), names(res)[ord], las = 2)
