df_all = data.frame(unique_number=colSums(CCI_count_all[,c(1:2)]))
df_all$age_group = rownames(df_all)
df_all$age_group <- factor(df_all$age_group,levels = c("Young","Old"))
gg_all <- ggplot(df_all, aes(x=age_group, y=unique_number,fill = age_group)) +
	        geom_bar(stat="identity", position=position_dodge()) + 
		        geom_text(aes(label=unique_number), vjust=-0.3, size=5, position = position_dodge(0.9)) + 
			        xlab("") +
				        ylab("Number of total interactions")  + 
					        theme_classic() +
						        scale_fill_manual(values = age_group_color) +
							        theme(axis.text.x = element_text(size = 10, family = "sans" , color = "Black"), 
								                    axis.text.y = element_text(size = 10, family = "sans" , color = "Black")) +
        theme(axis.title.x=element_text(face = "bold", size = 10, family = "sans" , color = "Black"),
	                    axis.title.y=element_text(face = "bold", size = 10, family = "sans" , color = "Black")) +
        theme(legend.position = "none") +
	        theme(plot.title = element_text(family = "sans", face = "bold", color = "Black", size = 20, hjust = 0.5,vjust = 0.5, angle = 0))
	gg_all
	ggsave("CCI_number_all_barplot.pdf",gg_all,width = 3, height = 4, dpi = 1200)

	df_mer$age_group <- factor(df_mer$age_group,levels=c("Young","Old"))
	df_mer %>%
		         mutate(cell_type = fct_reorder(cell_type, Number, .fun='median')) %>%
			          ggplot(aes(x=Number, y=cell_type, group=age_group, color=age_group)) + geom_point(size=7) +
				           scale_color_manual(values = age_group_color) +
					            ylab("") + 
						             xlab("Number of interactions") +
							              theme_bw() + 
								               theme(axis.text.x = element_text(size = 12,  color = "Black"), 
										                    axis.text.y = element_text(size = 12,  color = "Black")) +
         theme(axis.title.x=element_text(face = "bold", size = 20, color = "Black"),
	                      axis.title.y=element_text(face = "bold", size = 20, color = "Black")) + 
         guides(color=guide_legend(title = "Age group")) + 
	          theme(legend.title = element_text(size=20,  color = "Black"), 
			               legend.text = element_text(size=15,  color = "Black"))+
         theme(plot.title = element_text(hjust = 0.5,face = "bold", size = 20,  color = "Black"))

 ggsave("CCI_number_subtype_dotplot.pdf",width = 6, height = 4, dpi = 300)

#---------------------------------------------------------------------------------------------------------------------------------
 fc_interaction = log(cellchat@net$Old$count/  cellchat@net$Young$count, 2)
 bk = c(seq(-2,2,by=0.1))
 pdf('interaction.phm.pdf')
 pheatmap(fc_interaction,
	                 color = colorRampPalette(rev(brewer.pal(n = 7, name ="Spectral")))(length(bk)),cellheight = 15,cellwidth=15,
			                fontsize = 12, border_color=NA, breaks=bk, 
			                legend = TRUE, 
					               show_rownames=TRUE,cluster_cols = F,cluster_rows = F)
 dev.off()
