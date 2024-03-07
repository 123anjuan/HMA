pathways.show <- c("TGFb", 'ACTIVIN', 'CCL','CXCL', 'IL6','TNF','NOTCH','IGF','BMP','WNT')

 for (f in pathways.show) {
	 pdf(paste0(f,'.pdf'))
	 df = netVisual_aggregate(object.list[[1]], signaling = pathways.show, layout = "circle", color.use = color_key, edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[1]))
	 df = netVisual_aggregate(object.list[[2]], signaling = pathways.show, layout = "circle", color.use = color_key, edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[2]))
	 dev.off()
 }
