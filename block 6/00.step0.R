##Step 0: find markerpeaks
library(ArchR)
library(Seurat)
addArchRThreads(threads = 6)
addArchRGenome("hg38")
proj<-loadArchRProject(path = ./, force = T, showLogo = F)

proj<- addGroupCoverages(ArchRProj = proj, groupBy = "celltype",
)

proj<- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "celltype",
  force=T
)
 markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "celltype",
  testMethod = "wilcoxon"
)

#
peak<-getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")

for (i in 1:length(markerList)){
  differtest<-markerList[[i]]
  differtest<-as.data.frame(differtest)
  differtest<-differtest[,c('seqnames','start','end')]
  print(head(differtest))
  name = names(markerList)[i]
    print(name)
  write.table(differtest,file=paste0("differtest_peak_",i,".txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
}
