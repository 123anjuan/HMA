#snATAC-seq QC and clustering
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)

addArchRThreads(threads = 20) 
addArchRGenome("mm10")

args <- dir()[grep("gz$", dir())]
name <- gsub(".fragments.tsv.gz", "",args )
inputFiles <- args
names(inputFiles) <- name
names(inputFiles)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = name,
  filterTSS = 8, 
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

library(Seurat)
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, 
    knnMethod = "UMAP", 
    LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "./",
  copyArrows = F 
)

proj <- filterDoublets(proj, filterRatio = 2)

col <- rownames(proj@cellColData)
coln <- unlist(lapply(X = col, FUN = function(x) {return(strsplit(x, split = "#")[[1]][[1]])}))
head(coln)

coln1 <- unlist(lapply(X = coln, FUN = function(x) {return(strsplit(x, split = "-")[[1]][[2]])}))
proj$sample <- coln1

proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI1", 
    iterations = 5, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI1",
    name = "Harmony",
    groupBy = "Sample",
)

proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "Harmony", 
    name = "harmony_UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

proj <- addClusters(
    input = proj,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "harmony_Clusters",
    resolution = 2,
    maxClusters = NULL,force=T
)

seRNA <- readRDS("hM-sn_scRNA-seq_celltype.rds")

proj <- addImputeWeights(
  ArchRProj = proj,
  reducedDims = "Harmony",
  dimsToUse = NULL,
  scaleDims = NULL,
  corCutOff = 0.75,
  td = 3,
  ka = 4,
  sampleCells = 5000,
  nRep = 2,
  k = 15,
  epsilon = 1,
  useHdf5 = TRUE,
  randomSuffix = FALSE,
  threads = getArchRThreads(),
  seed = 1,
  verbose = TRUE,
  logFile = createLogFile("addImputeWeights")
)

proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = seRNA,
    reduction = "cca",
    addToArrow = FALSE,
    sampleCellsATAC = 10000,
    groupRNA = "celltype",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

pathToMacs2 <- "/home/an_juan/miniconda3/bin/macs2"
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "predictedGroup_Un")
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "predictedGroup_Un", 
    pathToMacs2 = pathToMacs2
)

proj <- addPeakMatrix(proj,force = T)
proj <- addMotifAnnotations(proj, motifSet = "cisbp", name = "Motif",force = T)

saveArchRProject(ArchRProj = proj, outputDirectory = "./", load = T, overwrite = F)
