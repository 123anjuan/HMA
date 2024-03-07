library(ggplot2)
library(ArchR)
library(Seurat)
addArchRGenome("hg38") #hg19/mm9/mm10(new)
addArchRThreads(threads = 6)
metadata<-read.csv("snRNA_fiber.meta.csv",row.names=1)
type1_pair<-readRDS("rna_atac_pairing.rds")
type1_pair<-as.data.frame(type1_pair)
row.names(type1_pair)<-type1_pair$RNA

type2_pair<-readRDS("rna_atac_pairing.rds")
type2_pair<-as.data.frame(type2_pair)
row.names(type2_pair)<-type2_pair$RNA

type1_meta<-metadata[type1_pair$RNA,]
type1_meta<-type1_meta[order(type1_meta$rank_type1),]
type1_pair<-type1_pair[row.names(type1_meta),]

type2_meta<-metadata[type2_pair$RNA,]
type2_meta<-type2_meta[order(type2_meta$rank_type2),]
type2_pair<-type2_pair[row.names(type2_meta),]

type1_meta$rank_rna<-paste("type1",type1_meta$rank_type1_new,sep="_")
type2_meta$rank_rna<-paste("type2",type2_meta$rank_type2_new,sep="_")
type1_pair$rank_rna<-type1_meta$rank_rna
type2_pair$rank_rna<-type2_meta$rank_rna

type1_meta<-type1_meta[type1_pair$RNA,]
type2_meta<-type2_meta[type2_pair$RNA,]

meta<-rbind(type1_meta,type2_meta)
pair<-rbind(type1_pair,type2_pair)

cell<-pair$ATAC

proj<-loadArchRProject(path = "fiber_0612", force = T, showLogo = F)
proj_sub<-proj[cell, ]
proj_sub@cellColData$rank_rna<-pair$rank_rna
p1 <- plotEmbedding(ArchRProj = proj_sub, colorBy = "cellColData", name = "rank_rna", embedding = "UMAPHarmony",baseSize=10,size = 1,labelAsFactors = F)

col = c(viridis(10),inferno(10))
names(col)<-c(paste("type1_",1:10,sep=""),paste("type2_",1:10,sep=""))


p2 <- plotEmbedding(ArchRProj = proj_sub, colorBy = "cellColData", name = "rank_rna", embedding = "UMAPHarmony",pal=col,
                    size = 0.2,
                    bgWidth = 0,fgColor = "#00000000",bgColor = "#00000000")+ theme_void() +  theme(legend.position="none")+   theme(title=element_blank())
plotPDF(p2, name = "ATAC_RNA_rank.pdf", ArchRProj = proj, addDOC = FALSE,width = 2.5, height = 2.5)

saveArchRProject(ArchRProj=proj_sub,outputDirectory="fiber_proj",overwrite = T,load=F,dropCells = F)

type1 <- meta[meta$rank_type1_new != "NA",]
cell.prop<-as.data.frame(prop.table(table(type1$group, type1$rank_type1_new)))
colnames(cell.prop)<-c("group","rank_type1_new","proportion")

p1<-ggplot(cell.prop,aes(rank_type1_new,proportion,fill=group))+
geom_bar(stat="identity",position="fill")+
ggtitle("")+
theme_bw()+
theme(axis.text.x = element_text(angle = 45),axis.ticks.length=unit(0.1,'cm'),panel.background=element_rect(fill='white'))+
guides(fill=guide_legend(title=NULL))+
scale_fill_manual(values=c('#A15DBA','#50BA47'))+ theme_classic()+
theme(axis.title = element_blank()) +
theme(axis.text.x = element_text(colour = "black", size = rel(2),angle = 45,hjust = 0.8,vjust = 0.8)) +
theme(axis.text.y = element_text(colour = "black", size = rel(2))) 

type2 <- meta[meta$rank_type2_new != "NA",]
cell.prop<-as.data.frame(prop.table(table(type2$group, type2$rank_type2_new)))
colnames(cell.prop)<-c("group","rank_type2_new","proportion")
p2<-ggplot(cell.prop,aes(rank_type2_new,proportion,fill=group))+
geom_bar(stat="identity",position="fill")+
ggtitle("")+
theme_bw()+
theme(axis.text.x = element_text(angle = 45),axis.ticks.length=unit(0.1,'cm'),panel.background=element_rect(fill='white'))+
guides(fill=guide_legend(title=NULL))+
scale_fill_manual(values=c('#A15DBA','#50BA47'))+ theme_classic()+
theme(axis.title = element_blank()) +
theme(axis.text.x = element_text(colour = "black", size = rel(2),angle = 45,hjust = 0.8,vjust = 0.8)) +
theme(axis.text.y = element_text(colour = "black", size = rel(2))) 

ggsave("./plot/barplot_type1_type2_age.pdf",p1 + p2,height = 4,width = 15)
