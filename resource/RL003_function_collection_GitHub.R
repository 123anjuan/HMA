combineSmallWnnClusters <- function(object=NA,resolutions=NA,graphNames=NA,minClusterSize=100) {
  for (graphName in graphNames) {
    comparison <- gsub("wsnn_","",graphName)
    for (resolution in resolutions) {
      object@meta.data[,paste0(graphName,"_res.",resolution,"_combined")] <- as.numeric(as.character(object@meta.data[,paste0(graphName,"_res.",resolution)]))
      for (myDistance in 1:5) {
        for (counter in 1:10) { 
          for (i in names(table(object@meta.data[,paste0(graphName,"_res.",resolution,"_combined")])[table(object@meta.data[,paste0(graphName,"_res.",resolution,"_combined")])<minClusterSize])) {
            for (cell in rownames(object@meta.data[object@meta.data[,paste0(graphName,"_res.",resolution,"_combined")]==i,])) {
              nearestNeighbor <- object@neighbors[[paste0("weighted.nn_",comparison)]]@cell.names[object@neighbors[[paste0("weighted.nn_",comparison)]]@nn.idx[object@neighbors[[paste0("weighted.nn_",comparison)]]@cell.names==cell,myDistance]]
              object@meta.data[cell,paste0(graphName,"_res.",resolution,"_combined")] <- object@meta.data[nearestNeighbor,paste0(graphName,"_res.",resolution,"_combined")]
            }
          }
        }
      }
    }
  }
  return(object)
}

compareGenotype <- function(ref_gt=NA,souporcell_output_dir=NA,sample_id=NA) {
  if (!file.exists(paste0(souporcell_output_dir,"/",sample_id,"/cluster_genotypes.vcf.gz"))) { print("souporcell vcf.gz not found") } else {
    print(sample_id)
    dgts <- as.data.frame(ref_gt$GT,stringsAsFactors = F)
    cell_data <- load_GT_vcf(paste0(souporcell_output_dir,"/",sample_id,"/cluster_genotypes.vcf.gz"),na.rm = F)
    cgts <- as.data.frame(cell_data$GT[!duplicated(rownames(cell_data$GT)),],stringsAsFactors = F)
    cgts <- cgts[,apply(cgts,2,function(x) sum(!is.na(x)))!=0]
    cgts <- cgts[rownames(cgts)%in%rownames(dgts),]
    cgts <- cgts[order(rownames(cgts)),]
    dgts <- dgts[rownames(dgts)%in%rownames(cgts),]
    dgts <- dgts[order(rownames(dgts)),]
    if (any(rownames(cgts)!=rownames(dgts))) { print("ref and soc vcf do not overlap") }
    #cor(cgts,as.data.frame(dgts[,patientList$NB_sample[patientList$CITE_fileName==i & !is.na(patientList$NB_sample)]]),use = "pairwise.complete.obs")
    myCors <- cor(cgts,as.data.frame(dgts),use="pairwise.complete.obs")
    return(myCors)
  }
}

runDoubletFinderOnSouporcellOutput <- function(object=NA) {
  object@meta.data[,c("pANN","DF.classifications","doubletFinder_params","df_classification_onSinglets")] <- NULL
  DefaultAssay(object) <- "RNA"
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
  object <- ScaleData(object)
  object <- RunPCA(object)
  object <- RunUMAP(object, dims = 1:15)
  object <- FindNeighbors(object, dims = 1:15)
  object <- FindClusters(object)
  
  ## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
  sweep.res.list_object <- paramSweep_v3(object, PCs = 1:15, sct = FALSE)
  gt.calls <- ifelse(object@meta.data$status=="doublet","Doublet","Singlet")   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
  expectedNewDoublets <- (sum(gt.calls=="Doublet")/length(gt.calls))*(1/length(unlist(str_split(unique(object@meta.data$pool_patients),";")))) #assuming that doublets from the same donor has the same frequency as multidonor doublets
  if (length(gt.calls)>10000) { gt.calls <- gt.calls[1:10000] } # Somehow this doesn't work with more than 10000 known doublets
  sweep.stats_object <- summarizeSweep(sweep.res.list_object, GT = TRUE, GT.calls = gt.calls)
  bcmvn_object <- try(find.pK(sweep.stats_object))
  
  bestPK <- as.numeric(as.character(bcmvn_object$pK[bcmvn_object$BCmetric==max(bcmvn_object$BCmetric)]))
  #bestPK <- 0.175 #Median for the ones that don't work
  annotations <- object@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(expectedNewDoublets*nrow(object@meta.data))  ## Assuming 5% extra ( with 4 genotypes we should already be quite OK )
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  object <- doubletFinder_v3(object, PCs = 1:15, pN = 0.25, pK = bestPK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  object@meta.data$DF.classifications <- object@meta.data[,grep("DF.classifications_",colnames(object@meta.data))]
  object@meta.data$pANN <- object@meta.data[,grep("pANN_",colnames(object@meta.data))]
  object@meta.data$df_classification_onSinglets <- ifelse(object@meta.data$pANN>=object@meta.data$pANN[object@meta.data$status!="doublet"][order(object@meta.data$pANN[object@meta.data$status!="doublet"],decreasing = T)][nExp_poi],"doublet","singlet")
  object@meta.data$doubletFinder_params <- paste0("pK",bestPK,";pN0.25;nExp",nExp_poi.adj,";",expectedNewDoublets*100,"%")
  
  return(object)
}

multiModal_processing <- function(object=NA,gex=T,adt=T,sct=T,gexAdtWnn=T,sctAdtWnn=T,doHarmony=T,npca=30,regress_cellcycle_gex=T,makeFinalWnnUmap=T,doFreshSct=F) {
  gc()
  if (gex) {
    # RNA assay
    DefaultAssay(object) <- 'RNA'
    if (regress_cellcycle_gex) {
      object <- NormalizeData(object)
      object <- CellCycleScoring(object, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
      object <- FindVariableFeatures(object) %>% ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) %>% RunPCA(reduction.name = "pca_RNA")
      gc()
    } else {
      object <- NormalizeData(object) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(reduction.name = "pca_RNA")
      gc()
    }
    object <- RunUMAP(object,reduction = "pca_RNA", dims = 1:npca,reduction.name = "umapBeforeHarmony_RNA",reduction.key='umapBeforeHarmony_RNA_')
    gc()
    if (doHarmony) { 
      object <- RunHarmony(object, group.by.vars = "orig.ident",assay.use = "RNA",reduction = "pca_RNA",reduction.save = "harmony_RNA")
      object <- RunUMAP(object, reduction = "harmony_RNA", dims = 1:npca,reduction.name = "umapAfterHarmony_RNA",reduction.key='umapAfterHarmony_RNA_')
      gc()
    }
  }
  
  # ADT assay
  if (adt) {
    DefaultAssay(object) <- 'ADT'
    VariableFeatures(object) <- rownames(object[["ADT"]])
    object <- NormalizeData(object, normalization.method = 'CLR', margin = 2) %>% ScaleData() %>% RunPCA(reduction.name = "pca_ADT")
    object <- RunUMAP(object,reduction = "pca_ADT", dims = 1:npca,reduction.name = "umapBeforeHarmony_ADT",reduction.key='umapBeforeHarmony_ADT_')
    gc()
    if (doHarmony) { 
      object <- RunHarmony(object, group.by.vars = "orig.ident",assay.use = "ADT",reduction = "pca_ADT",reduction.save = "harmony_ADT")
      object <- RunUMAP(object, reduction = "harmony_ADT", dims = 1:npca,reduction.name = "umapAfterHarmony_ADT",reduction.key='umapAfterHarmony_ADT_')
      gc()
    }
  }
  
  # Run WNN with RNA + ADT
  if (gexAdtWnn) {
    if (doHarmony) {
      object <- FindMultiModalNeighbors(
        object, reduction.list = list("harmony_RNA", "harmony_ADT"), 
        dims.list = list(1:npca, 1:npca), modality.weight.name = "RNA.weight_rnaAdt",
        knn.graph.name = "wknn_rnaAdt",
        snn.graph.name = "wsnn_rnaAdt",
        weighted.nn.name = "weighted.nn_rnaAdt"
      )
    } else {
      object <- FindMultiModalNeighbors(
        object, reduction.list = list("pca_RNA", "pca_ADT"), 
        dims.list = list(1:npca, 1:npca), modality.weight.name = "RNA.weight_rnaAdt",
        knn.graph.name = "wknn_rnaAdt",
        snn.graph.name = "wsnn_rnaAdt",
        weighted.nn.name = "weighted.nn_rnaAdt"
      )
    }
    gc()
    if (makeFinalWnnUmap) {
      object <- RunUMAP(object, nn.name = "weighted.nn_rnaAdt", reduction.name = "wnn.umap_rnaAdt", reduction.key = "wnnUMAP_rnaAdt_")
      gc()
    }
  }
  
  # SCT assay
  if (sct) {
    if (doFreshSct) {
      suppressWarnings(object <- SCTransform(object, assay = "RNA" , verbose = F,conserve.memory=T))
      gc()
    }
    DefaultAssay(object) <- 'SCT'
    VariableFeatures(object[["SCT"]]) <- rownames(object[["SCT"]]@scale.data)
    object <- RunPCA(object,reduction.name = "pca_SCT")
    object <- RunUMAP(object,reduction = "pca_SCT", dims = 1:npca,reduction.name = "umapBeforeHarmony_SCT",reduction.key='umapBeforeHarmony_SCT_')
    gc()
    if (doHarmony) { 
      object <- RunHarmony(object, group.by.vars = "orig.ident",assay.use = "SCT",reduction = "pca_SCT",reduction.save = "harmony_SCT")
      gc()
      object <- RunUMAP(object, reduction = "harmony_SCT", dims = 1:npca,reduction.name = "umapAfterHarmony_SCT",reduction.key='umapAfterHarmony_SCT_')
      gc()
    }
  }
  # Run WNN with SCT + ADT
  if (sctAdtWnn) {
    object <- FindMultiModalNeighbors(
      object, reduction.list = list("harmony_SCT", "harmony_ADT"), 
      dims.list = list(1:npca, 1:npca), modality.weight.name = "RNA.weight_sctAdt",
      knn.graph.name = "wknn_sctAdt",
      snn.graph.name = "wsnn_sctAdt",
      weighted.nn.name = "weighted.nn_sctAdt"
    )
    gc()
    if (makeFinalWnnUmap) {
      object <- RunUMAP(object, nn.name = "weighted.nn_sctAdt", reduction.name = "wnn.umap_sctAdt", reduction.key = "wnnUMAP_sctAdt_")
      gc()
    }
  }
  return(object)
  gc()
}

downloadScData <- function(cite=NA, gex=NA, bcr=NA, tcr=NA, overwrite=F, alignment="cellranger", out_dir=NA) {
  #download GEX CITE BCR and TCR data from iRODS
  out_dir <- gsub("(.*)/$","\\1",out_dir)
  if (!is.na(cite) & !is.na(gex)) { cat("supply either a gex or cite sample, not both\n") } else {
    if (is.na(cite) & !is.na(gex)) { 
      sample_id <- gex 
      if (!dir.exists(paste0(out_dir,"/gex"))) { dir.create(paste0(out_dir,"/gex")) }
      new_dir <- paste0(out_dir,"/gex/",sample_id)
      dir.create(new_dir)
    }
    if (!is.na(cite) & is.na(gex)) { 
      sample_id <- cite 
      if (!dir.exists(paste0(out_dir,"/cite"))) { dir.create(paste0(out_dir,"/cite")) }
      new_dir <- paste0(out_dir,"/cite/",sample_id)
      dir.create(new_dir)
    }
    if ((!is.na(cite) | !is.na(gex)) & alignment=="cellranger") {
      doesDataExist <- tryCatch(system(paste0("ils /archive/HCA/10X/",sample_id,"/outs/filtered_feature_bc_matrix")))
      if (doesDataExist==0) {
        if (!overwrite) { 
          system(paste0("iget -r /archive/HCA/10X/",sample_id,"/outs/filtered_feature_bc_matrix /archive/HCA/10X/",sample_id,"/outs/analysis /archive/HCA/10X/",sample_id,"/outs/raw_feature_bc_matrix ",new_dir)) 
        } else {
          system(paste0("iget -fr /archive/HCA/10X/",sample_id,"/outs/filtered_feature_bc_matrix /archive/HCA/10X/",sample_id,"/outs/analysis /archive/HCA/10X/",sample_id,"/outs/raw_feature_bc_matrix ",new_dir)) 
        }
      } else { cat(paste(sample_id,"does not seem to exist\n")) }
    }
  }
  if (!is.na(bcr)) {
    if (!dir.exists(paste0(out_dir,"/bcr"))) { dir.create(paste0(out_dir,"/bcr")) }
    new_dir <- paste0(out_dir,"/bcr/",bcr)
    dir.create(new_dir)
    doesDataExist <- tryCatch(system(paste0("ils /archive/HCA/10X-VDJ/",bcr,"/outs/filtered_contig.fasta")))
    if (doesDataExist==0) {
      if (!overwrite) { 
        system(paste0("iget -r /archive/HCA/10X-VDJ/",bcr,"/outs/filtered_contig.fasta /archive/HCA/10X-VDJ/",bcr,"/outs/filtered_contig_annotations.csv ",new_dir))
      } else {
        system(paste0("iget -fr /archive/HCA/10X-VDJ/",bcr,"/outs/filtered_contig.fasta /archive/HCA/10X-VDJ/",bcr,"/outs/filtered_contig_annotations.csv ",new_dir))
        
      }
    } else { cat(paste(bcr,"BCR file does not seem to exist\n")) }
  }
  if (!is.na(tcr)) {
    if (!dir.exists(paste0(out_dir,"/tcr"))) { dir.create(paste0(out_dir,"/tcr")) }
    new_dir <- paste0(out_dir,"/tcr/",tcr)
    dir.create(new_dir)
    doesDataExist <- tryCatch(system(paste0("ils /archive/HCA/10X-VDJ/",tcr,"/outs/filtered_contig.fasta")))
    if (doesDataExist==0) {
      if (!overwrite) { 
        system(paste0("iget -r /archive/HCA/10X-VDJ/",tcr,"/outs/filtered_contig.fasta /archive/HCA/10X-VDJ/",tcr,"/outs/filtered_contig_annotations.csv ",new_dir))
      } else {
        system(paste0("iget -fr /archive/HCA/10X-VDJ/",tcr,"/outs/filtered_contig.fasta /archive/HCA/10X-VDJ/",tcr,"/outs/filtered_contig_annotations.csv ",new_dir))
      }
    } else { cat(paste(bcr,"TCR file does not seem to exist\n")) }
  }
}

processCiteSamples <- function(sample=NA, SoupX_rna=T, SoupX_adt=T, save_raw=F, doSct=T, data_dir=NA, min_cells=0, min_features=200, add_sample_name_to_cellbarcode=T) {
  if (!dir.exists(paste0(data_dir,"/",sample,"/filtered_feature_bc_matrix/"))) { print(paste0(sample,": data not found")) } else {
    rawData = Read10X(data.dir = paste0(data_dir,"/",sample,"/raw_feature_bc_matrix/"))
    filData = Read10X(data.dir = paste0(data_dir,"/",sample,"/filtered_feature_bc_matrix/"))
    clusters <- read.csv(paste0(data_dir,"/",sample,"/analysis/clustering/graphclust/clusters.csv"))
    
    if (SoupX_rna) {
      tod = rawData[['Gene Expression']]
      toc = filData[['Gene Expression']]
      sc = SoupChannel(tod = tod,toc = toc)
      sc = setClusters(sc, setNames(clusters$Cluster, clusters$Barcode))
      sc = autoEstCont(sc)
      out = Matrix::drop0(adjustCounts(sc))
      filData[["RNA"]] <- out
    } else { filData[["RNA"]] <- filData[['Gene Expression']] }
    if (SoupX_adt) {
      # Do SoupX on ADT with parameters used by Krzysztof and LFH
      tod = rawData[['Antibody Capture']]
      toc = filData[['Antibody Capture']]
      sc = SoupChannel(tod = tod,toc = toc)
      sc = setClusters(sc, setNames(clusters$Cluster, clusters$Barcode))
      soupQuantile = 0.25
      tfidfMin = 0.2
      try(sc <- autoEstCont(sc, soupQuantile=soupQuantile, tfidfMin=tfidfMin, forceAccept=TRUE,doPlot=F,verbose = F))
      message = try(adjustCounts(sc))
      #we really don't want this to fail though. so keep dropping the parameters and trying again
      #(once tfidMin hits 0.05 the next iteration puts it at 0, so just give up)
      while ((class(message) == 'try-error') && (tfidfMin > 0.05))
      {
        soupQuantile = soupQuantile - 0.05
        tfidfMin = tfidfMin - 0.05
        try(sc <- autoEstCont(sc, soupQuantile=soupQuantile, tfidfMin=tfidfMin, forceAccept=TRUE,doPlot=F,verbose = F))
        message = try(adjustCounts(sc))
      }
      # message = try(sc <- autoEstCont(sc, soupQuantile=0.25, tfidfMin=0.2, forceAccept=TRUE))
      # if (class(message) == 'try-error') {
      #   sc = autoEstCont(sc, soupQuantile=0.05, tfidfMin=0.05, forceAccept=TRUE)
      # }
      if (class(message) == 'try-error') { 
        cat(paste0("SoupX on ADT failed for ",sample,", raw ADT data returned instead\n")) 
        filData[["ADT"]] <- filData[['Antibody Capture']]
        soupxAdtFailed <- T
      } else {
        out = Matrix::drop0(adjustCounts(sc))
        filData[["ADT"]] <- out
        soupxAdtFailed <- F
      }
    } else { filData[["ADT"]] <- filData[['Antibody Capture']] }
    
    fil <- CreateSeuratObject(counts = filData[["RNA"]],min.cells = min_cells, min.features = min_features,project = sample,assay = "RNA")
    fil[["ADT"]] <- CreateAssayObject(counts = filData[["ADT"]][,colnames(fil)])
    if (save_raw) {
      fil[["ADT_raw"]] <- CreateAssayObject(counts = filData[['Antibody Capture']][,colnames(fil)])
      fil[["RNA_raw"]] <- CreateAssayObject(counts = filData[['Gene Expression']][,colnames(fil)])
    }
    if (doSct) {
      suppressWarnings(fil <- SCTransform(fil, assay = "RNA" , verbose = F))
    }
    if (add_sample_name_to_cellbarcode) { 
      fil <- RenameCells(fil,add.cell.id = sample)
    }
    if (SoupX_adt) {
      fil@meta.data$soupxOnAdt <- ifelse(soupxAdtFailed,"Fail","Pass")
    }
    return(fil)
  }
}

# Below is written by or adapted from Natsuhiko Kumasaka
col.rb <-
  c("#053061", "#073568", "#0A3A70", "#0D4077", "#10457F", "#134B86", 
    "#15508E", "#185696", "#1B5B9D", "#1E61A5", "#2166AC", "#246AAE", 
    "#286FB0", "#2B74B3", "#2F78B5", "#327DB7", "#3581BA", "#3986BC", 
    "#3C8ABE", "#408FC1", "#4494C3", "#4C99C6", "#549EC9", "#5CA3CB", 
    "#64A8CE", "#6CADD1", "#74B2D3", "#7CB7D6", "#84BCD9", "#8CC1DC", 
    "#93C5DE", "#9AC9E0", "#A0CCE2", "#A7CFE4", "#ADD2E5", "#B3D6E7", 
    "#BAD9E9", "#C0DCEB", "#C6DFED", "#CDE3EE", "#D2E5F0", "#D6E7F0", 
    "#DAE9F1", "#DEEBF2", "#E1EDF3", "#E5EEF3", "#E9F0F4", "#EDF2F5", 
    "#F1F4F5", "#F5F6F6", "#F7F5F4", "#F7F2EF", "#F8EFEA", "#F9EDE6", 
    "#F9EAE1", "#FAE7DC", "#FAE4D7", "#FBE1D2", "#FCDECD", "#FCDCC8", 
    "#FCD7C2", "#FBD2BB", "#FACCB4", "#F9C7AD", "#F8C1A6", "#F7BC9F", 
    "#F7B799", "#F6B192", "#F5AC8B", "#F4A684", "#F1A07E", "#EE9978", 
    "#EB9273", "#E88B6E", "#E58468", "#E27D63", "#DF765E", "#DC6F58", 
    "#D96853", "#D6614E", "#D35A4A", "#CF5246", "#CB4B43", "#C8443F", 
    "#C43D3C", "#C03539", "#BD2E35", "#B92732", "#B6202E", "#B2182B", 
    "#AB1529", "#A31328", "#9C1027", "#940E26", "#8C0C25", "#850923", 
    "#7D0722", "#760421", "#6E0220", "#67001F")

drawDendrogram <-
  function(tree, a, b, transpose=F){
    m = tree$merge
    o = tree$order
    h = tree$height/max(tree$height)
    
    # number of leavs
    K = length(o)
    # number of nodes
    N = nrow(m)
    
    x = (matrix(match(-m, o), N)-1)/(K-1)
    y = matrix(rep(0, N*2), N)
    
    #plot(1,1,xlim=c(0,1),ylim=c(0,1),type="n")
    if(transpose==F){
      for(i in 1:N){
        if(is.na(x[i,1])){
          x[i,1] = mean(x[m[i,1],])
          y[i,1] = h[m[i,1]]
        }
        if(is.na(x[i,2])){
          x[i,2] = mean(x[m[i,2],])
          y[i,2] = h[m[i,2]]
        }
        lines(x[i,c(1,1,2,2)]*b[1]+a[1], c(y[i,1],h[i],h[i],y[i,2])*b[2]+a[2])
      }
    }else{
      for(i in 1:N){
        if(is.na(x[i,1])){
          x[i,1] = mean(x[m[i,1],])
          y[i,1] = h[m[i,1]]
        }
        if(is.na(x[i,2])){
          x[i,2] = mean(x[m[i,2],])
          y[i,2] = h[m[i,2]]
        }
        lines(c(y[i,1],h[i],h[i],y[i,2])*b[2]+a[2], x[i,c(1,1,2,2)]*b[1]+a[1])
      }
    }
  }

Dotplot <-
  function(X, ltsr, col=col.rb, cex=1, zlim=c(-1,1), xlab="", ylab="", lab=rep("",ncol(X)), SORT=c(F,F), beta, border=F, add=F, xlim=c(0,ncol(X)), tick=F, side=2, w=0.03, ww=50, labx=NULL, laby=NULL,srt=0,cex.axis=1){
    dotplot=function(x,y,z,col,add,zlim,axes,xlab,ylab,xlim=range(x),ylim=range(y),ltsr,cex=1){
      z=ceiling((z-zlim[1])/diff(range(zlim))*99)+1
      ltsr=-log10((1-ltsr)*2)*4
      par(xpd=NA)
      plot(rep(x,length(y)),rep(y,rep(length(x),length(y))),cex=sqrt(ltsr)*cex, pch=20, col=col[z],
           xlim=range(x)+c(0,0.5),axes=axes,xlab=xlab,ylab=ylab)
      par(xpd=F)
    }
    
    X[X>zlim[2]]=zlim[2]
    X[X<zlim[1]]=zlim[1]
    ltsr[ltsr>0.9999]=0.9999
    N0=N=ncol(X)
    M0=M=nrow(X)
    M1=N1=0
    
    if(SORT[1]){
      dendx = hclust(dist(t(X)))
      X = X[, dendx$ord]
      ltsr = ltsr[, dendx$ord]
    }
    if(SORT[2]){
      dendy = hclust(dist(X))
      X = X[dendy$ord,]
      ltsr = ltsr[dendy$ord,]
    }
    if(!is.null(labx)){
      X=rbind(X,array(NA,c(ncol(labx)+1, N)))
      rownames(X)[(M+2):(M+ncol(labx)+1)]=colnames(labx)
      M1=ncol(labx)+1
      M=M0+M1
    }
    if(!is.null(laby)){
      X=cbind(array(NA,c(M, ncol(laby)+1)),X)
      colnames(X)[(1):(ncol(laby))]=colnames(laby)
      N1=ncol(laby)+1
      N=N0+N1
    }
    
    if(ncol(X)>1){
      dotplot(1:N-0.5, 1:M-0.5, t(X[nrow(X):1,,drop=F]),col=col,zlim=zlim,axes=F,xlab=xlab,ylab=ylab,add=add, xlim=xlim,ltsr=t(ltsr[nrow(ltsr):1,,drop=F]), cex=cex)
    }else{
      dotplot(1:2-0.5, 1:M-0.5, t(cbind(X,NA)[nrow(X):1,,drop=F]),col=col,zlim=zlim,axes=F,xlab=xlab,ylab=ylab,add=add, xlim=xlim+c(0,1),ltsr=t(cbind(ltsr,NA)[nrow(X):1,,drop=F]), cex=cex)
    }
    if(!is.null(labx)){
      labx=cbind(NA,matrix(as.numeric(as.factor(unlist(lapply(labx,as.character)))),N0))
      if(SORT[1]){image(N1:N, 1:M1-0.5, (labx[dendx$ord,ncol(labx):1]),col=c(col24,col24)[1:max(labx,na.rm=T)], add=T)
      }else{image(N1:N, 1:M1-0.5, (labx[,ncol(labx):1]),col=c(col24,col24)[1:max(labx,na.rm=T)], add=T)}
    }
    if(!is.null(laby)){
      laby=cbind(matrix(as.numeric(as.factor(unlist(lapply(laby,as.character)))),M0),NA)
      if(SORT[2]){ image(1:N1-0.5, M1:M, t(laby[rev(dendy$ord),]),col=c(col24,col24)[1:max(laby,na.rm=T)], add=T)
      }else{image(1:N1-0.5, M1:M, t(laby[nrow(laby):1,]),col=c(col24,col24)[1:max(laby,na.rm=T)], add=T)}
    }
    
    if(border){
      for(i in (-1):(N)){segments(i+0.5,0.5,i+0.5,M+0.5)}
      for(i in (-1):(M)){segments(0.5,i+0.5,N+0.5,i+0.5)}
    }
    mtext(lab,1,at=1:N)
    if(tick){
      flag=!is.na(rownames(X))
      par(xpd=NA)
      segments(N,c(M:1-0.5)[flag],N+0.1,c(M:1-0.5)[flag])
      segments(N+0.1,c(M:1-0.5)[flag],N+1, c(M:1-0.5)[flag]*w + (rank((M:1-0.5)[flag])/sum(flag)*M-ww)*(1-w))
      text(N+1, c(M:1-0.5)[flag]*w + (rank((M:1-0.5)[flag])/sum(flag)*M-ww)*(1-w), rownames(X)[flag], pos=4, offs=0.1, cex=0.75)
      par(xpd=F)
    }else{
      mtext(rownames(X),at=M:1-0.5,side,las=2,line=0,cex = cex.axis)
    }
    par(xpd=NA)
    text(1:ncol(X)-0.2,rep(-.5,ncol(X)),colnames(X),pos=2,srt=srt,offs=0,cex = cex.axis)
    if(SORT[1]) drawDendrogram(dendx, c(N1+0.5,M), c(N0-1,3))
    if(SORT[2]) drawDendrogram(dendy, c(M-0.5,N), c(-(M0-1),3), T)
    #Gauge
    ybase=-1
    rect(rep(N+1+3-0.25,100),(seq(100)-1)/100*6+10+ybase,rep(N+1+3+0.25,100),seq(100)/100*6+10+ybase,border=NA,col=col.rb)
    text(N+5.5,16.7+ybase,"Fold change",offs=0.5)
    text(rep(N+2+3,3), c(0.02,0.5,0.97)*6+10+ybase, c(paste("<1/",exp(zlim[2]),sep=""),"1",paste(">",exp(zlim[2]),sep="")),pos=4,offs=-0.5)
    #LTSR
    points(rep(N+1+3,5),1:5*1.2+ybase,cex=sqrt(-log10((1-c(0.501,0.9,0.99,0.999,0.9999))*2)*4)*cex,pch=20)
    text(rep(N+2+3,5),  1:5*1.2+ybase,c("0.5","0.9","0.99","0.999",">0.9999"),pos=4,offs=-0.3)
    text(N+3+3,7.2+ybase,"LTSR",offs=0.5)
    par(xpd=F)
  }

Dotplot_forCorHeatmap <-
  function(X, ltsr, col=col.rb, cex=1, zlim=c(-1,1), xlab="", ylab="", lab=rep("",ncol(X)), SORT=c(F,F), beta, border=F, add=F, xlim=c(0,ncol(X)), tick=F, side=2, w=0.03, ww=50, labx=NULL, laby=NULL,measure="r(s)",cex.axis=1,srt=0){
    dotplot=function(x,y,z,col,add,zlim,axes,xlab,ylab,xlim=range(x),ylim=range(y),ltsr,cex=1){
      z=ceiling((z-zlim[1])/diff(range(zlim))*99)+1
      ltsr=-log10((1-ltsr)*2)*4
      par(xpd=NA)
      plot(rep(x,length(y)),rep(y,rep(length(x),length(y))),cex=sqrt(ltsr)*cex , pch=20, col=col[z],
           # plot(rep(x,length(y)),rep(y,rep(length(x),length(y))),cex=sqrt(ltsr)*cex , pch=20, col=col[z],
           xlim=range(x)+c(0,0.5),axes=axes,xlab=xlab,ylab=ylab)
      par(xpd=F)
    }
    
    X[X>zlim[2]]=zlim[2]
    X[X<zlim[1]]=zlim[1]
    ltsr[ltsr>0.999]=0.999
    N0=N=ncol(X)
    M0=M=nrow(X)
    M1=N1=0
    
    if(SORT[1]){
      dendx = hclust(dist(t(X)))
      X = X[, dendx$ord]
      ltsr = ltsr[, dendx$ord]
    }
    if(SORT[2]){
      dendy = hclust(dist(X))
      X = X[dendy$ord,]
      ltsr = ltsr[dendy$ord,]
    }
    if(!is.null(labx)){
      X=rbind(X,array(NA,c(ncol(labx)+1, N)))
      rownames(X)[(M+2):(M+ncol(labx)+1)]=colnames(labx)
      M1=ncol(labx)+1
      M=M0+M1
    }
    if(!is.null(laby)){
      X=cbind(array(NA,c(M, ncol(laby)+1)),X)
      colnames(X)[(1):(ncol(laby))]=colnames(laby)
      N1=ncol(laby)+1
      N=N0+N1
    }
    
    if(ncol(X)>1){
      dotplot(1:N-0.5, 1:M-0.5, t(X[nrow(X):1,,drop=F]),col=col,zlim=zlim,axes=F,xlab=xlab,ylab=ylab,add=add, xlim=xlim,ltsr=t(ltsr[nrow(ltsr):1,,drop=F]), cex=cex)
    }else{
      dotplot(1:2-0.5, 1:M-0.5, t(cbind(X,NA)[nrow(X):1,,drop=F]),col=col,zlim=zlim,axes=F,xlab=xlab,ylab=ylab,add=add, xlim=xlim+c(0,1),ltsr=t(cbind(ltsr,NA)[nrow(X):1,,drop=F]), cex=cex)
    }
    if(!is.null(labx)){
      labx=cbind(NA,matrix(as.numeric(as.factor(unlist(lapply(labx,as.character)))),N0))
      if(SORT[1]){image(N1:N, 1:M1-0.5, (labx[dendx$ord,ncol(labx):1]),col=c(col24,col24)[1:max(labx,na.rm=T)], add=T)
      }else{image(N1:N, 1:M1-0.5, (labx[,ncol(labx):1]),col=c(col24,col24)[1:max(labx,na.rm=T)], add=T)}
    }
    if(!is.null(laby)){
      laby=cbind(matrix(as.numeric(as.factor(unlist(lapply(laby,as.character)))),M0),NA)
      if(SORT[2]){ image(1:N1-0.5, M1:M, t(laby[rev(dendy$ord),]),col=c(col24,col24)[1:max(laby,na.rm=T)], add=T)
      }else{image(1:N1-0.5, M1:M, t(laby[nrow(laby):1,]),col=c(col24,col24)[1:max(laby,na.rm=T)], add=T)}
    }
    
    if(border){
      for(i in (-1):(N)){segments(i+0.5,0.5,i+0.5,M+0.5)}
      for(i in (-1):(M)){segments(0.5,i+0.5,N+0.5,i+0.5)}
    }
    mtext(lab,1,at=1:N)
    if(tick){
      flag=!is.na(rownames(X))
      par(xpd=NA)
      segments(N,c(M:1-0.5)[flag],N+0.1,c(M:1-0.5)[flag])
      segments(N+0.1,c(M:1-0.5)[flag],N+1, c(M:1-0.5)[flag]*w + (rank((M:1-0.5)[flag])/sum(flag)*M-ww)*(1-w))
      text(N+1, c(M:1-0.5)[flag]*w + (rank((M:1-0.5)[flag])/sum(flag)*M-ww)*(1-w), rownames(X)[flag], pos=4, offs=0.1, cex=0.75)
      par(xpd=F)
    }else{
      mtext(rownames(X),at=M:1-0.5,side,las=2,line=0,cex = cex.axis)
    }
    par(xpd=NA)
    text(1:ncol(X)-0.2,rep(-.5,ncol(X)),colnames(X),pos=2,srt=srt,offs=0,cex = cex.axis)
    if(SORT[1]) drawDendrogram(dendx, c(N1+0.5,M), c(N0-1,3))
    if(SORT[2]) drawDendrogram(dendy, c(M-0.5,N), c(-(M0-1),3), T)
    #Gauge
    ybase=-1
    rect(rep(N+1+3-0.25,100),(seq(100)-1)/100*6+10+ybase,rep(N+1+3+0.25,100),seq(100)/100*6+10+ybase,border=NA,col=col.rb)
    text(N+5.5,16.7+ybase,measure,offs=0.5)
    text(rep(N+2+3,3), c(0.02,0.5,0.97)*6+10+ybase, c(paste("-",zlim[2],sep=""),"0",zlim[2]),pos=4,offs=-0.5)
    #LTSR
    points(rep(N+1+3,4),1:4*1.2+ybase,cex=sqrt(-log10((1-c(0.501,0.9,0.99,0.999))*2)*4)*cex,pch=20)
    text(rep(N+2+3,4),  1:4*1.2+ybase,c("0.5","0.1","0.01","<0.001"),pos=4,offs=-0.3)
    text(N+3+3,7.2+ybase,"P value",offs=0.5)
    par(xpd=F)
  }

dplot <-
  function(x,y=NULL,z,dims=c(100,100),col=colorRampPalette(c("#053061", "white", "#67001F"))(100),cex=1){
    if(!is.null(y)){x=cbind(x,y)}
    x1=ceiling((x[,1]-min(x[,1],na.rm=T)+1e-20)/(diff(range(x[,1],na.rm=T))+1e-20)*dims[1])
    x2=ceiling((x[,2]-min(x[,2],na.rm=T)+1e-20)/(diff(range(x[,2],na.rm=T))+1e-20)*dims[2])
    z=unlist(lapply(split(z,paste(x1,x2)),mean,na.rm=T))
    x=apply(t(matrix(unlist(strsplit(names(z)," ")),2)),2,as.numeric)
    gplot(x,,z,col=col,sym=F,cex=cex)
  }
gplot <-
  function(x, y=NULL, z, col=colorRampPalette(c("#053061", "white", "#67001F"))(100), symm=T, zlim=NULL, gauge=T, ...){
    x=x[order(z),]
    z=z[order(z)]
    discretise<-function(x,k=100,symm=T,zlim=NULL){
      if(!is.null(zlim)){
        x[x<zlim[1]]=zlim[1]
        x[x>zlim[2]]=zlim[2]
      }else{
        zlim = range(x,na.rm=T)
      }
      if(symm){
        x=c(x,-x)
        x=floor((x-zlim[1])/diff(zlim)*(k));
        x[x%in%k]=(k-1);
        x[1:(length(x)/2)]+1
      }else{
        x=floor((x-zlim[1])/diff(zlim)*(k)); 
        x[x%in%k]=(k-1); 
        x+1
      }
    }
    plot(x, y, col=col[discretise(z,symm=symm,zlim=zlim)], pch=20, axes=T, xlab="", ylab="",...)
    if(gauge){
      yy = seq(max(x[,2])-diff(range(x[,2]))*0.3, max(x[,2]), diff(range(x[,2]))*0.3/100)
      xx1= min(x[,1])
      xx2= min(x[,1])+diff(range(x[,1]))*0.05
      rect(rep(xx1,100),yy[-101],rep(xx2,100),yy[-1],col=col,border=NA)
      rect(xx1,min(yy),xx2,max(yy),col=NA,lwd=3)
      #par(xpd=NA);text(rep(mean(c(xx1,xx2)),2),yy[c(1,101)],c("low","high"),pos=c(1,3))
      par(xpd=NA);text(rep(mean(c(xx1,xx2)),2),yy[c(1,101)],c(round(min(z),2),round(max(z),2)),pos=c(1,3))
    }
  }

# obtaining posterior mean and
# local false sign rate (lfsr)

getCondVal <- function(res.prop.ranef, id, ncells, nfactors=2){
  tmp = data.frame(res.prop.ranef)[data.frame(res.prop.ranef)[[1]]==id,]
  if(length(grep(":",tmp$grp))==0){
    cnam = matrix(as.character(tmp$term),ncells)[1,]
    rnam = matrix(as.character(tmp$grp),ncells)[,1] 
  }else if(nfactors==2){
    cnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[1,],ncells)[1,]
    rnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[2,],ncells)[,1]
  }else if(nfactors==3){
    cnam1 = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[1,],ncells)[1,]
    cnam2 = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[2,],ncells)[1,]
    cnam = paste(cnam1,cnam2,sep=":")
    rnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[3,],ncells)[,1]
  }
  condval = matrix(tmp$condval,ncells)
  condsd  = matrix(tmp$condsd, ncells)
  rownames(condval)=rownames(condsd)=rnam
  colnames(condval)=colnames(condsd)=cnam
  lfsr = pnorm(condval,0,condsd)
  lfsr[lfsr>0.5]=1-lfsr[lfsr>0.5]
  list(condval=condval, lfsr=lfsr)
}

getCondVal=function(res.prop.ranef, id, ncells, nfactors=2, celltype){
  tmp = data.frame(res.prop.ranef)[data.frame(res.prop.ranef)[[1]]==id,]
  if(length(grep(":",tmp$grp))==0){
    cnam = matrix(as.character(tmp$term),ncells)[1,]
    rnam = matrix(as.character(tmp$grp),ncells)[,1]
  }else if(nfactors==2){
    a=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[1,])
    b=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[2,])
    tmp=tmp[match(paste(a[rep(1:length(a),rep(length(b),length(a)))],b[rep(1:length(b),length(a))],sep=":"),tmp$grp),]
    tmp$grp=paste(a[rep(1:length(a),rep(length(b),length(a)))],b[rep(1:length(b),length(a))],sep=":")
    cnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[1,],ncells)[1,]
    rnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[2,],ncells)[,1]
  }else if(nfactors==3){
    a=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[1,])
    b=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[2,])
    C=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[3,])
    ab=paste(a[rep(1:length(a),rep(length(b),length(a)))],b[rep(1:length(b),length(a))],sep=":")
    abc=paste(ab[rep(1:length(ab),rep(length(C),length(ab)))],C[rep(1:length(C),length(ab))],sep=":")
    tmp=tmp[match(abc,tmp$grp),]
    tmp$grp=abc
    cnam1 = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[1,],ncells)[1,]
    cnam2 = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[2,],ncells)[1,]
    cnam = paste(cnam1,cnam2,sep=":")
    rnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[3,],ncells)[,1]
  }
  condval = matrix(tmp$condval,ncells)
  condsd  = matrix(tmp$condsd, ncells)
  condval[is.na(condval)]=0
  condsd[is.na(condsd)]=1
  rownames(condval)=rownames(condsd)=rnam
  condval=condval[match(celltype,rnam),]
  condsd =condsd[match(celltype,rnam),]
  colnames(condval)=colnames(condsd)=cnam
  lfsr = pnorm(condval,0,condsd)
  lfsr[is.na(lfsr)]=1
  lfsr[lfsr>0.5]=1-lfsr[lfsr>0.5]
  list(condval=condval,  lfsr=lfsr)
}

Forest <-
  function(x,ord=T,labs=rownames(x),xlim=c(-1,1.5)){
    x=x[rownames(x)!="theta.Celltype.(Intercept)",]
    if(ord){labs=labs[rev(order(x[,1]))]; x=x[rev(order(x[,1])),]}
    labs=gsub("Celltype.","",gsub("X10x","10x_kit",gsub("Sample","residual",gsub(":Celltype.\\(Intercept\\)","",gsub("theta.","",rownames(x))))))
    x=rbind(x[labs!="residual",],x[labs=="residual",])
    labs=c(labs[labs!="residual"],labs[labs=="residual"])
    n=nrow(x)
    x=cbind(x,x[,1]-x[,2]*1.96,x[,1]+x[,2]*1.96);
    plot(x[,1],n:1,type="n",axes=F,xlab="",ylab="",xlim=xlim)
    points(x[,1],n:1,pch=15)
    axis(1)
    #segments(xlim[1],(n:1)[x[,3]<xlim[1]],xlim[1]+0.05,(n:1)[x[,3]<xlim[1]]+0.1)
    #segments(xlim[1],(n:1)[x[,3]<xlim[1]],xlim[1]+0.05,(n:1)[x[,3]<xlim[1]]-0.1)
    #x[x[,3]<xlim[1],3]=xlim[1]
    segments(x[,3],n:1,x[,4],n:1)
    abline(v=0,lty=2)
    par(xpd=NA)
    text(rep(xlim[1],nrow(x)),n:1,labs,pos=2)
    x
  }
