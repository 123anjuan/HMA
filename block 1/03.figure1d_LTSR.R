#LTSR
args <- commandArgs(T)
source('RL003_function_collection_GitHub.R')
library(Seurat)
library(ggplot2)
library(lme4)
library(showtext)
showtext_auto()
library(egg)

meta = read.csv(args[1],row.names=1) #the matedata of merged snRNA,scRNA and snATAC GeneScore RDS
Y = table(meta$patient_id_sample,meta$Final_annotation)

mymetadata <- meta[!duplicated(meta$patient_id_sample),]##library
mymetadata$Gender <- mymetadata$gender
mymetadata$Ethnicity <- mymetadata$Country
mymetadata$Batch <- mymetadata$tech
mymetadata$AgeGroup <- mymetadata$age_pop  #age_pop


metadata <- mymetadata[,c("Age","Batch","Gender","Ethnicity","patient_id_sample")] 

metadata$Age =as.numeric(metadata$Age)
range(metadata$Age)
mean(metadata$Age)

metadata$Age = scale(as.numeric(metadata$Age))##Centralization, non-standardization
range(metadata$Age)
mean(metadata$Age)


Y <- Y[rownames(Y)%in%metadata$patient_id_sample,]
dim(Y)
# number of samples / number of cell types
nsamples = nrow(Y)
ncells = ncol(Y)
# repeating the meta data table by the number of cell types
metadataExp=cbind(metadata[rep(match(rownames(Y),as.character(metadata$patient_id_sample)),ncells),],Celltype=rep(colnames(Y),rep(nsamples,ncells)))
head(metadataExp)
head(Y)

res.prop=glmer(I(c(Y))~Age
+(1|Celltype)              
+(1|patient_id_sample)
+(1|Gender)
+(1|Ethnicity)
+(1|Batch)
               
+(Age-1|Celltype)     
+(1|patient_id_sample:Celltype)
+(1|Gender:Celltype)
+(1|Ethnicity:Celltype)
+(1|Batch:Celltype),
#+(1|AgeGroup:Celltype),

family=poisson,data=metadataExp,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

# standard errors of standard deviations (squre root of the variance parameters)]
summary(res.prop)
library('numDeriv')

devfun = update(res.prop, devFunOnly=T)
pars = getME(res.prop, c("theta","fixef"))
hess = hessian(devfun, unlist(pars))
sdse.prop = data.frame(sd=unlist(pars), se=sqrt(diag(solve(hess))))

# posterior means and their standard deviations
res.prop.ranef = ranef(res.prop)

# Forest plot
rownames(sdse.prop)[rownames(sdse.prop)=="theta.Celltype.(Intercept)"] <- "Residual"
par(mar=c(3,6,1,1),mgp=c(1.2,0.5,0))
Forest(sdse.prop[grep("(Celltype|Residual)",rownames(sdse.prop)),],xlim=c(0,2.5))

# Forest plot
rownames(sdse.prop)[rownames(sdse.prop)=="theta.Celltype.(Intercept)"] <- "Residual"
par(mar=c(3,6,1,1),mgp=c(1.2,0.5,0))
Forest(sdse.prop[grep("(Celltype|Residual)",rownames(sdse.prop)),],xlim=c(0,2.5))

data = getCondVal(res.prop.ranef,"Celltype",ncells,celltype=colnames(Y))
data1 = getCondVal(res.prop.ranef,"Batch:Celltype",ncells,celltype=colnames(Y))
data2 = getCondVal(res.prop.ranef,"Gender:Celltype",ncells,celltype=colnames(Y))
data3 = getCondVal(res.prop.ranef,"Ethnicity:Celltype",ncells,celltype=colnames(Y))


postmean = cbind(
    NA,
    data[[1]],
   NA,
    data1[[1]],
     NA,
    data2[[1]],
     NA
 #   data3[[1]]
    
)
lfsr = cbind(
    NA,
    data[[2]], 
    NA,
    data1[[2]],
     NA,
    data2[[2]],
     NA
 #   data3[[2]]
)
# Dotplot
postmean_oldAgeGroupsPlusSeverity <- postmean
lfsr_oldAgeGroupsPlusSeverity <- lfsr
myClust <- hclust(dist(postmean_oldAgeGroupsPlusSeverity*(1-lfsr_oldAgeGroupsPlusSeverity)),method = "complete")$order
postmean_oldAgeGroupsPlusSeverity <- postmean_oldAgeGroupsPlusSeverity[myClust,]
lfsr_oldAgeGroupsPlusSeverity <- lfsr_oldAgeGroupsPlusSeverity[myClust,]
par(mar=c(13,12,8,10))
Dotplot(postmean_oldAgeGroupsPlusSeverity, SORT=c(F,F),zlim=c(log(1/3),log(3)),ltsr=1-lfsr_oldAgeGroupsPlusSeverity, cex=0.8,srt=90,cex.axis = .8)

postmean_oldAgeGroupsPlusSeverity = postmean_oldAgeGroupsPlusSeverity[,-3]
lfsr_oldAgeGroupsPlusSeverity= lfsr_oldAgeGroupsPlusSeverity[,-3]
myClust = order(rownames(postmean_oldAgeGroupsPlusSeverity))
myClust
postmean_oldAgeGroupsPlusSeverity <- postmean_oldAgeGroupsPlusSeverity[myClust,]
lfsr_oldAgeGroupsPlusSeverity <- lfsr_oldAgeGroupsPlusSeverity[myClust,]
par(mar=c(13,12,8,10))
Dotplot(postmean_oldAgeGroupsPlusSeverity, SORT=c(F,F),zlim=c(log(1/3),log(3)),ltsr=1-lfsr_oldAgeGroupsPlusSeverity, cex=0.8,srt=90,cex.axis = .8)
