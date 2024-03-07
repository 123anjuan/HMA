args <- commandArgs(T)

library(Seurat)
library(ggplot2)
library(dplyr)
library(monocle3)
library(patchwork)
library(sctransform)

type1 <- readRDS(args[1]) 
type2 <- readRDS(args[2])

df <- type1@meta.data %>%
  group_by(rank_type1, group) %>%
  summarise(
    per_anno = n()) %>% data.frame()
df$percent <- 0
for (i in seq_along(1:nrow(df))){
    sgnum <- colSums(as.matrix(df[df$rank_type1 == df[i,"rank_type1"],"per_anno"]))
    df$percent[i] <- df$per_anno[i] / sgnum
}
color_age <- c('Adult'='#50BA47','Old'='#A15DBA')
p <- ggplot(data = df, mapping = aes(
  x = `rank_type1`, fill = `group`, y = `percent`))
p1 <- p + geom_col() + theme_classic()+
scale_fill_manual(values=color_age) + theme(axis.title = element_blank()) +
theme(axis.text.x = element_text(colour = "black", size = rel(2),angle = 45,hjust = 0.8,vjust = 0.8)) +
theme(axis.text.y = element_text(colour = "black", size = rel(2))) 
p1

df2 <- type2@meta.data %>%
  group_by(rank_type2, group) %>%
  summarise(
    per_anno = n()) %>% data.frame()
df2$percent <- 0
for (i in seq_along(1:nrow(df2))){
    sgnum <- colSums(as.matrix(df[df2$rank_type2 == df2[i,"rank_type2"],"per_anno"]))
    df2$percent[i] <- df2$per_anno[i] / sgnum
}
p <- ggplot(data = df2, mapping = aes(
  x = `rank_type2`, fill = `group`, y = `percent`))
p2 <- p + geom_col() + theme_classic()+
scale_fill_manual(values=color_age) + theme(axis.title = element_blank()) +
theme(axis.text.x = element_text(colour = "black", size = rel(2),angle = 45,hjust = 0.8,vjust = 0.8)) +
theme(axis.text.y = element_text(colour = "black", size = rel(2))) 
p2
