suppressPackageStartupMessages({
library(RColorBrewer)
library(ArchR)
library(dplyr)
library(ggplot2)
})

proj <-  loadArchRProject("./")

markersTest_com <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "celltype_group",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = "Quie.MuSC_Old",
  bgdGroups = "Quie.MuSC_Adult"
)

motifsUp <- peakAnnoEnrichment(markersTest_com, proj, peakAnnotation = "Motif", cutOff = "FDR <= 0.01 & Log2FC >= 1")
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
df <- df[!grepl("LINE.", df$TF),]

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(20)), ], aes(x = rank, y = mlog10Padj, label = TF), label.r = 0.15, box.padding = 0.25, label.size = 0.5,
        segment.size = 1,
        size = 4.2,
        nudge_x = 2,
        color = "black",max.overlaps=40
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

motifsDown <- peakAnnoEnrichment(markersTest_com, proj, peakAnnotation = "Motif", cutOff = "FDR <= 0.01 & Log2FC <= -1")
df_down <- data.frame(TF = rownames(motifsDown), mlog10Padj = assay(motifsDown)[,1])
df_down <- df_down[order(df$mlog10Padj, decreasing = TRUE),]
df_down$rank <- seq_len(nrow(df_down))
df_down <- df_down[!grepl("LINE.", df$TF),]

ggDo <- ggplot(df_down, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(8)), ], aes(x = rank, y = mlog10Padj, label = TF), label.r = 0.15, box.padding = 0.25, label.size = 0.5,
        segment.size = 1,
        size = 4.2,
        nudge_x = 2,
        color = "black",max.overlaps=60
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
