coding_gene <- read.csv("hsapiens_ensembl_gene.txt",sep = "\t")

exp <- t(as.matrix(GetAssayData(type1, assay = 'RNA', slot = 'data')))
type1_exp = aggregate(exp, FUN = mean, by = list(type1@meta.data$rank))
type1_exp %>% 
  column_to_rownames("Group.1") -> df
df_t <- t(df)
df_t_coding <- intersect(rownames(df_t),coding_gene$SYMBOL)
df_t_gene <- df_t[rownames(df_t) %in% df_t_coding,]
gene_sd<-apply(df_gene,1,sd)
df_sd <- cbind(df_gene,gene_sd)
write.csv(df_sd,"Type1_sd.csv")

t2 <- t(as.matrix(GetAssayData(type2, assay = 'RNA', slot = 'data')))
type2_exp = aggregate(t2, FUN = mean, by = list(type2@meta.data$rank))
type2_exp %>%  column_to_rownames("Group.1") -> df2
df2_t <- t(df2)
df2_t_gene <- df2_t[rownames(df2_t) %in% coding_gene$SYMBOL,]
gene_sd2<-apply(df2_t_gene,1,sd)
df2_sd <- cbind(df2_t_gene,gene_sd2)
write.csv(df2_sd,"Type2_sd.csv")

type1 <- read.csv("Type1_sd.csv",row.names = 1)
type2 <- read.csv("Type2_sd.csv",row.names = 1)

type1_t4000_gene <- rownames(type1)[1:4000]
type2_t4000_gene <- rownames(type2)[1:4000]

h_gene_name <- union(type1_t4000_gene,type2_t4000_gene)

type1_select <- type1[rownames(type1) %in% h_gene_name,]
type2_select <- type2[rownames(type2) %in% h_gene_name,]

type1_select <- type1_select[sort(rownames(type1_select)),]
type2_select <- type2_select[sort(rownames(type2_select)),]

type1_select_t <- t(type1_select)
type2_select_t <- t(type2_select)

reversed_type2 <- type2_select_t[rev(row.names(type2_select_t)), ]

matrix <- rbind(type1_select_t,reversed_type2)
mat_cha <- t(scale(matrix))

mat_cha[mat_cha >= 5] = 5
mat_cha[mat_cha <= -5] = -5
mat_cha_type1 <- mat_cha[,1:100]
mat_cha_type1 <- t(apply(mat_cha_type1,1,function(x){smooth.spline(x,df=3)$y}))
mat_cha_type2 <- mat_cha[,101:200]
mat_cha_type2 <- t(apply(mat_cha_type2,1,function(x){smooth.spline(x,df=3)$y}))
mat_cha_smooth <- cbind(mat_cha_type1,mat_cha_type2)
mat_cha_smooth[mat_cha_smooth >= 5] = 5
mat_cha_smooth[mat_cha_smooth <= -5] = -5

k_means <- 10
p_clust=pheatmap::pheatmap(mat_cha,
                     cutree_rows=k_means, breaks = seq(-5,5, by = 0.1),
                     cluster_cols = F,show_rownames = FALSE,
                     show_colnames = FALSE      )
split_matrix=data.frame(cutree(p_clust$tree_row,k=k_means))


clust_hm2 <- Heatmap(mat_cha_smooth,colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100),
                    left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 2:11),
                                      labels_gp = gpar(col = "white", fontsize = 10))),
                    row_split = split_matrix, 
                    cluster_columns = FALSE,
                    #right_annotation = ha,
                    show_row_names = FALSE,
                    show_column_names = FALSE
                    )
