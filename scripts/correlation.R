library(pheatmap)
library(RColorBrewer)

breaksList = seq(0.6, 1, by = 0.01)
breaksList = c(seq(0, 0.5, by = 0.05),breaksList)
palette = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))

args = commandArgs(trailingOnly=TRUE)
folder = args[1]
filename = args[2]
genome = args[3]
out = args[4]

background = read.csv(paste("data/",genome,".bg.txt",sep=""),sep="\t",header=FALSE)

# List all codons from AAA to TTT
codons = c()
for (i in c("A","C","G","T")) {
  for (j in c("A","C","G","T")) {
    for (k in c("A","C","G","T")) {
      codons = c(codons,paste(i,j,k,sep=""))
    }
  }
}

alldat = matrix(nrow=64,ncol=0)

for (i in 1:9) {
	files=system(paste("ls ",folder,"/",filename,".codon.",i,"*.txt",sep=""),intern=TRUE)
	
	dat = matrix(nrow=64,ncol=0)
	rownames(dat) = codons
	
	for (file in files) {
		info = file.info(file)
		if (info$size == 0) next
		
		name = strsplit(gsub(filename,"",gsub(folder,"",file,fixed=TRUE),fixed=TRUE),"\\.")[[1]][4]
		
		tmp = read.csv(paste(file),header=FALSE,stringsAsFactors=FALSE,strip.white=TRUE)
		tmp2 = matrix(unlist(strsplit(tmp$V1," ")),ncol=2,byrow=TRUE)
		
		# Handle missing codons and set them to zero
		new = matrix(rep(0,64),nrow=64,ncol=1)
		rownames(new) = codons
		new[tmp2[,2],1] = as.numeric(tmp2[,1])
		
		dat = cbind(dat,as.numeric(new))
		colnames(dat)[dim(dat)[2]] = name
	}
	colnames(dat) = paste(colnames(dat),".",i,sep="")
	alldat = cbind(alldat,dat)
}
alldat = cbind(alldat,background[,2])
colnames(alldat)[dim(alldat)[2]] = "BG"

# Create annotation
ann = matrix(ncol=1,nrow=ncol(alldat))
colnames(ann) = c("Position")
ann_colors = list()
ann[grep("BG",colnames(alldat)),1] = "control"
ann[grep("\\.4",colnames(alldat)),1] = "E"
ann[grep("\\.5",colnames(alldat)),1] = "P"
ann[grep("\\.6",colnames(alldat)),1] = "A"
ann[grep("\\.[123]",colnames(alldat)),1] = "up"
ann[grep("\\.[789]",colnames(alldat)),1] = "down"
ann_colors[["Position"]] = c("black","darkgrey","lightgrey","#339999","purple","orange")
names(ann_colors[["Position"]]) = c("control","up","down","E","P","A")
ann = as.data.frame(ann)
rownames(ann) = colnames(alldat)

pdf(paste(out,"/",filename,".correlation.pdf",sep=""),onefile=FALSE,width=10,height=9)
pheatmap(cor(alldat,alldat),border=NA,border_color=NA,cex=1,annotation_row=ann,annotation_col=ann,annotation_colors=ann_colors,clustering_method="average",clustering_distance_rows="correlation",clustering_distance_cols="correlation",color = palette, breaks = breaksList, fontsize=7)
dev.off()
pdf(paste(out,"/",filename,".correlation.order.pdf",sep=""),onefile=FALSE,width=10,height=9)
pheatmap(cor(alldat,alldat),border=NA,border_color=NA,cex=1,annotation_row=ann,annotation_col=ann,annotation_colors=ann_colors,cluster_rows = FALSE,cluster_cols=FALSE,color = palette, breaks = breaksList, fontsize=7)
dev.off()

write.csv(cor(alldat,alldat),paste(out,"/",filename,".correlation.csv",sep=""),quote=FALSE)
