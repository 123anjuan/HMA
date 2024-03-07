 immune = c('TNF','CXCL','MIF','CCL','TWEAK','IL1','TIGIT','IL6','IFN-II')
 rkn = rankNet(
	         cellchat,
		   slot.name = "netP", , stacked = T,
		     signaling = immune,
		   measure = "count",
		     mode = c("comparison", "single"),
		     comparison = c(1, 2),
		       color.use = age_group_color)
 ggsave('immune_net.pdf',rkn, width = 4, height =4)


 factor = c('WNT','IGF','BMP','NOTCH','EGF','FGF','HGF')
 rkn = rankNet(
	         cellchat,
		   slot.name = "netP", , stacked = T,
		     signaling = factor,
		   measure = "count",
		     mode = c("comparison", "single"),
		     comparison = c(1, 2),
		       color.use = age_group_color)
 ggsave('factor.pdf',rkn, width = 4, height =4)

 ecm_pw_f= c('HSPG','LAMININ','COLLAGEN','FN1','TGFb','RELN','SPP1','MK')
 rkn = rankNet(
	         cellchat,
		   slot.name = "netP", , stacked = T,
		     signaling = ecm_pw_f,
		   measure = "count",
		     mode = c("comparison", "single"),
		     comparison = c(1, 2),
		       color.use = age_group_color)
 ggsave('ecm_pw_f.pdf',rkn, width = 4, height =4)
