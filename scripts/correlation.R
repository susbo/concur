# R script to plot correlations between different read sets.

# Change log:
# 2025-08-20 Edited by SB to remove calculation of read length correlation. This
#            is now read from .csv file instead to (avoid R dependency for pipeline)

library(pheatmap)
library(RColorBrewer)

breaksList = seq(0.6, 1, by = 0.01)
breaksList = c(seq(0, 0.5, by = 0.05),breaksList)
palette = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))

args = commandArgs(trailingOnly=TRUE)
filename = args[1]
out = args[2]

# Read saved correlation from file
cor_matrix = read.csv(paste(out,"/",filename,".correlation.csv",sep=""),row.names=1,header=TRUE,check.names=FALSE)

# Create annotation
ann = matrix(ncol=1,nrow=ncol(cor_matrix))
colnames(ann) = c("Position")
ann_colors = list()
ann[grep("BG",colnames(cor_matrix)),1] = "control"
ann[grep("\\.4",colnames(cor_matrix)),1] = "E"
ann[grep("\\.5",colnames(cor_matrix)),1] = "P"
ann[grep("\\.6",colnames(cor_matrix)),1] = "A"
ann[grep("\\.[123]",colnames(cor_matrix)),1] = "up"
ann[grep("\\.[789]",colnames(cor_matrix)),1] = "down"
ann_colors[["Position"]] = c("black","darkgrey","lightgrey","#339999","purple","orange")
names(ann_colors[["Position"]]) = c("control","up","down","E","P","A")
ann = as.data.frame(ann)
rownames(ann) = colnames(cor_matrix)

pdf(paste(out,"/",filename,".correlation.pdf",sep=""),onefile=FALSE,width=10,height=9)
pheatmap(cor_matrix,border=NA,border_color=NA,cex=1,annotation_row=ann,annotation_col=ann,annotation_colors=ann_colors,clustering_method="average",clustering_distance_rows="correlation",clustering_distance_cols="correlation",color = palette, breaks = breaksList, fontsize=7)
dev.off()

pdf(paste(out,"/",filename,".correlation.order.pdf",sep=""),onefile=FALSE,width=10,height=9)
pheatmap(cor_matrix,border=NA,border_color=NA,cex=1,annotation_row=ann,annotation_col=ann,annotation_colors=ann_colors,cluster_rows = FALSE,cluster_cols=FALSE,color = palette, breaks = breaksList, fontsize=7)
dev.off()
