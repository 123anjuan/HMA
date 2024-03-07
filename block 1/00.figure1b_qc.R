#sn/scRNA-seq QC
args <- commandArgs(T)

library(Seurat)
library(SoupX)
library(DoubletFinder)
library(dplyr)

soupx<-function(FilterMatrix,RawMatrix){
	toc <- Read10X(FilterMatrix,gene.column =1)
    tod <- Read10X(RawMatrix,gene.column =1)
    tod <- tod[rownames(toc),]
    all <- toc
    all <- CreateSeuratObject(all)
    all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
    all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
    all.genes <- rownames(all)
    all <- ScaleData(all, features = all.genes)
    all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = F)
    all <- FindNeighbors(all, dims = 1:30)
    all <- FindClusters(all, resolution = 0.5)
    all <- RunUMAP(all, dims = 1:30)
    matx <- all@meta.data
    sc = SoupChannel(tod, toc)
    sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
    sc = autoEstCont(sc)
    out = adjustCounts(sc)
    return(out)
}

doubletfinder<-function(Count,SampleID,Sample){
	data <- Read10X(Count)
	obj <- CreateSeuratObject(counts=data,project=SampleID,min.cells = 3, min.features = 200)
	obj <- NormalizeData(obj)
	obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
	obj <- ScaleData(obj)
	obj <- RunPCA(obj)
	obj <- RunUMAP(obj, dims = 1:10)
	sweep.res.list_obj <- paramSweep_v3(obj, PCs = 1:10, sct = FALSE)
	sweep.stats_obj <- summarizeSweep(sweep.res.list_obj, GT = FALSE)
	bcmvn_obj <- find.pK(sweep.stats_obj)
	mpK<-as.numeric(as.vector(bcmvn_obj$pK[which.max(bcmvn_obj$BCmetric)]))
	annotations <- obj@meta.data$orig.ident
	homotypic.prop <- modelHomotypic(annotations)
	nExp_poi <- round(0.075*ncol(obj@assays$RNA))
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	obj_filterDouble <- doubletFinder_v3(obj, PCs = 1:40, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
	colnames(obj_filterDouble@meta.data)[length(obj_filterDouble@meta.data)] <- "DoubleScore"
	obj_singlet <- subset(obj_filterDouble,cells=WhichCells(obj_filterDouble,expression=`DoubleScore`!="Doublet"))
	obj_singlet@meta.data <- select(obj_singlet@meta.data,orig.ident,nCount_RNA,nFeature_RNA)
	obj_singlet$sample <- Sample
	return(obj_singlet) 
}

filter_count<-args[1]
raw_count<-args[2]
sampleid<-args[3]
sample<-args[4]

clean_matrix<-soupx(filter_count,raw_count)

singlet<-doubletfinder(clean_matrix,sampleid,sample)
saveRDS(singlet,args[5])
