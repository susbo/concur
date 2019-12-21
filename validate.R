args = commandArgs(trailingOnly=TRUE)
out = args[1]
name = args[2]
genome = args[3]

dat = read.csv(paste(out,"/",name,".codons.txt",sep=""),sep="\t",header=TRUE,stringsAsFactors = FALSE,row.names=1)

colors = rep("#33333366",64)

start = c("ATG")
for (s in start) colors[grep(s,rownames(dat))] = "#FF00FFCC"

stop = c("TAG","TAA","TGA")
for (s in stop) colors[grep(s,rownames(dat))] = "#0000FFCC"

common = system(paste("cat data/",genome,".bg.txt | sort -k2,2nr | head -n10 | cut -f1",sep=""),intern=TRUE)
common = common[!common %in% c(start,stop)][1:3]
for (s in common) colors[grep(s,rownames(dat))] = "#00FF00CC"

rare = system(paste("cat data/",genome,".bg.txt | sort -k2,2n | head -n10 | cut -f1",sep=""),intern=TRUE)
rare = rare[!rare %in% c(start,stop)][1:3]
for (s in rare) colors[grep(s,rownames(dat))] = "#FF0000CC"

pdf(paste(out,"/",name,".pdf",sep=""))
plot(1:9,unlist(dat[1,]),type="l",ylim=c(0,max(dat)),lwd=3,col=colors[1],main="Codon frequency per position",ylab="Count",xlab="Position",xaxt="n")
for (i in 2:64) {
	lines(1:9,unlist(dat[i,]),type="l",lwd=3,col=colors[i])
}
legend("topleft",cex=0.8,c("Start codon","Stop codons",paste("3 most common [",paste(common,collapse=','),"]",sep=""),paste("3 most rare [",paste(rare,collapse=','),"]",sep=""),"All other codons"),col=c("#FF00FFCC","#0000FFCC","#00FF00CC","#FF0000CC","#33333366"),lwd=3)
axis(1,at=1:9,labels = c("-5","-4","-3","-2 [E]","-1 [P]","0 [A]","+1","+2","+3"))
dev.off()
