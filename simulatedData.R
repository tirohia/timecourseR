library(gplots)
FileChoice<-file.choose() # choose input file
b <- read.delim(FileChoice, skip=0, sep="\t", as.is=TRUE)
rownames(b)<-b[,1];b<-b[,-1]
b<-log2(b) #optionally log2-transformn the data
k=10
c<-kmeans(b,k)
s<-cbind(c$cluster,b)
colnames(s)[1]<-"CL"
lb1<-length(b[1,])
x<-1:lb1
s.ordered<- s[order(s$CL),]; s.ordered<-data.frame(s.ordered)
pdf(file=paste(date(),"k_means_graphs.pdf") ,width=6,height=6)
for (cycle1 in 1: k) {
  r<- s.ordered[s.ordered$CL==cycle1,-1]; l=length(r[,1])
  print(length(rownames(r)))
  textplot(paste("k =",cycle1))
  textplot(rownames(r))
  for (cycle2 in 1: l-1) {
    plot(x,r[cycle2,],col=cycle1+1,type="l",ylab="",ylim=range(min(r),max(r)) ,cex.axis=0.5)
    par(new=TRUE)
  }
  plot(x,r[l,],col=cycle1+1,type="l",ylab="",ylim=range(min(r),max(r)),main=paste("k =",cycle1) ,cex.axis=0.5)
  abline(v=8.5); abline(v=16.5) #use if you have blocks of chips you want to separate
}
dev.off()