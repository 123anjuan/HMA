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

###step4:Enrichment.pl
#!/usr/bin/perl -w
use strict;
open IN, "phenotype.tsv" || die "$!";
while(<IN>){
  next if($. == 1);
  chomp;
  my($a, $b) = split(/\s+/, $_);
  `python ldsc.py --h2-cts $b --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. --out $a --ref-ld-chr-cts ATAC_Disease.ldcts --w-ld-chr weights_hm3_no_hla/weights.`;
}

library(tidyverse)
library(ggplot2)
res_files=list.files(pattern=".cell_type_results.txt")

data_list=lapply(res_files, function(x){read.csv(file=x,header=T, sep="\t")})

data_list_gwas=Map(cbind, data_list, gwas = gwas_list)

data.df=do.call(rbind,data_list_gwas)

bigplot=ggplot(data.df, aes(x=gwas,y=-log10(Coefficient_P_value), color=Name)) +
  labs(y="-log10(Coefficient P-value)")+
  geom_point(size=3) + coord_flip()+ scale_color_brewer(palette="Set3")+ theme_bw()

