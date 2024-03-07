#sn/scRNA-seq clustering
library(Seurat)
library(dplyr)
library(future)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 50000*1024^2)

rds.list <- dir()[grep("*.rds",dir())]

for (i in seq(1,length(rds.list),1)) {
	    assign(paste("rds", i ,sep="."), readRDS(rds.list[i]))
        assign(paste("name", i ,sep="."), levels(get(paste("rds",i,sep="."))@meta.data$orig.ident))
}

rds_all <- as.list(mget(ls(pattern = "rds.")[1:length(rds.list)]))

ifnb.list <- lapply(X = rds_all, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#SCTransform
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)

ifnb.list <- lapply(X = ifnb.list, FUN = RunPCA, features = features)

obj.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
    anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
obj.combined.sct <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT", dims = 1:30) 

obj.combined.sct <- RunPCA(obj.combined.sct, verbose = FALSE)
obj.combined.sct <- RunUMAP(obj.combined.sct, reduction = "pca", dims = 1:30)
obj.combined.sct <- FindNeighbors(obj.combined.sct, reduction = "pca", dims = 1:30)
obj.combined.sct <- FindClusters(obj.combined.sct, resolution = 0.8)

DefaultAssay(obj.combined) <- "RNA"
