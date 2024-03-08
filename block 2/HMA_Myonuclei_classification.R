#Load libraries:
library(SingleCellExperiment)
library(Seurat)
library(dplyr)
library(tibble)
library(stringr)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(future)
library(ggplot2)
library(sleepwalk)
library(SCINA)
library(prettydoc)
library(wesanderson)
library(parallel)
library(MASS)
library(sctransform)
library(RColorBrewer)
library(harmony)

# Load opt options or introduce here the input files
output <- opt$output_folder
rds_name <- opt$data_file
dataset_name <- opt$name_dataset

# Load object
new_OM <- readRDS(file = rds_name)
setwd(output)
metadata <- new_OM@meta.data

#Calculate Type I and Type II score to classify each of them. Next, calculate Type IIa and Type IIx scores to further discriminate Type II fiber types
DefaultAssay(new_OM) <- "RNA"
new_OM$scoreI  <- PercentageFeatureSet(object = new_OM, pattern = "TNNT1|MYH7|MYH7B|TNNC1|TNNI1|ATP2A2")
new_OM$scoreII  <- PercentageFeatureSet(object = new_OM, pattern = "TNNT3|MYH1|MYH2|TNNC2|TNNI2|ATP2A1")
new_OM$scoreIIa  <- PercentageFeatureSet(object = new_OM, pattern  = "MYH2|ANKRD2|NDUFA8|MYOM3|CASQ2|HSPB6|RDH11|AIMP1")
new_OM$scoreIIx  <- PercentageFeatureSet(object = new_OM, pattern = "MYH1|MYLK2|ACTN3|MYBPC2|PCYOX1|CAPZA1|CD38|PDLIM7|COBL|TMEM159|HNRNPA1|TFRC")

I <- new_OM$scoreI
II <- new_OM$scoreII
IIa <- new_OM$scoreIIa
IIx <- new_OM$scoreIIx

new_OM$IorII <- I/II
new_OM$IIclass <- IIx/IIa
new_OM$fiber_class <- new_OM@meta.data$age_pop
i <- 1
for(i in 1:length(rownames(new_OM@meta.data))){
  x <- new_OM@meta.data$IorII[i]
  z <- new_OM@meta.data$scoreIIa[i]
  y <- new_OM@meta.data$IIclass[i]
  if(x == "NaN"){
    new_OM@meta.data$fiber_class[i] <- "Unknown MF"
  }else if(x == "Inf"){
    new_OM@meta.data$fiber_class[i] <- "Pure Type I"
  }else if(x > 1){
    if(z < 0.25){
      new_OM@meta.data$fiber_class[i] <- "Type I"
    }else{
      new_OM@meta.data$fiber_class[i] <- "Hybrid I/IIa"
    }
  }else{
    if(y == "NaN"){
      new_OM@meta.data$fiber_class[i] <- "Unknown II"
    }else if(y == "Inf"){
      new_OM@meta.data$fiber_class[i] <- "Pure Type IIx"
    }else if(y == 0){
      new_OM@meta.data$fiber_class[i] <- "Pure Type IIa"
    }else if(y > 1.25){
      new_OM@meta.data$fiber_class[i] <- "Type IIx"
    }else if(y < 0.75){
      new_OM@meta.data$fiber_class[i] <- "Type IIa"
    }else if(y < 1.25 & y > 0.75){
      new_OM@meta.data$fiber_class[i]  <- "Hybrid IIx/IIa"
    }else{
      new_OM@meta.data$fiber_class[i] <- "Unknown II-subtype"
    }
    
  }
}
x <- Idents(new_OM)
Idents(new_OM) <- new_OM$fiber_class
cell_number <- table(Idents(new_OM), new_OM$sample)
write.table(cell_number, file = "Myonuclei_class_V5_number.txt", sep = "\t")
Idents(new_OM) <- x

#8. Save RDS object

saveRDS(new_OM, file="Myonuclei_class_V5.rds")


