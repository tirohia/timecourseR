setwd("/home/ben/workspace/timeCourse/data")
source("../code/R/tools.R")


counts = read.table( datafile,sep=",", header=TRUE,row.names=1 , stringsAsFactors = FALSE )
counts<-cleanDataNames(counts)
delist<-c("IYO_14535","IYO_13640","IYO_18340","IYO_23120","IYO_03285","IYO_05360","IYO_13220","IYO_17550","IYO_06520","IYO_27715","IYO_27720")
counts<-counts[delist,]

design <- getDataRange(counts,24,120)
design<-design[design$condition=="psa",]
cols<-design[design$time==24,3:5]
collist<-unlist(cols)
countcols<-counts[,design$collist]
#horrid kludge
countcols<-counts[,c(3,41,42,10 ,43, 44,6,8,12,55,56,57,49,50,51)]
colnames(countcols)<-rep(c(24,48,72,96,120),each=3)


par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(0,0, xlim=c(20, 120), ylim=c(0, max(log(counts))),xlab="HPI",ylab="log Reads" ,type="n",main="Differentially expressed psa genes",cex=2,, bty='L')
xl <- seq(min(design$time),max(design$time), (max(design$time) - min(design$time))/1000)
y.leg <- c(11,10.5,10,9.5,9,8.5,8,7.5,7,6.5,6)

for (i in 1:nrow(countcols)) {
  psadf<-data.frame(time=rep(c(24,48,72,96,120),each=3), reads=t(countcols[i,]))
  name<-colnames(psadf)[2]
  colnames(psadf)<-c("time","geneName")
  psadf[,2]<-log(psadf[,2]+0.1)
  loess_fit <- loess(geneName ~ time ,psadf, na.action=na.exclude)
  lines(xl, predict(loess_fit,xl), col = i*3 ,lwd=2) 
  legend(130,y.leg[i],inset=c(-.4,0),legend=rownames(countcols)[i],col = i*3,ncol = 2, cex = 0.75, lwd = 3, text.font = i, text.col = i*3,y.intersp=0,lty = c(1,2))
  #legend(2.8,-1,c("group A", "group B"), pch = c(1,2), lty = c(1,2))
}















y.leg<-y.leg-5
length(y.leg)
for (i in 1:nrow(countcols)) {
  #legend(2,legend=rownames(colcounts), col = 2:3, pch = 1:2, lty = 1, merge = TRUE)   #, trace = TRUE)
  legend(3, y.leg[i],legend=rownames(countcols)[i] , pch = "sSvV", col = c(1, 3))

}
warnings()
design
psadf<-data.frame(time=unique(design$time))
psadf
plot(0,0, xlim=c(20, 120), ylim=c(0, max(log(counts[delist,]))),xlab="HPI",ylab="RPM" ,type="n",main="genes",cex=2)
loess_fit <- loess("IYO_03285" ~ time , counts[,1:5], na.action=na.exclude)
counts


###
notableList<-read.csv("deResults/transcriptionsFactors-foundby3methods.csv",header=FALSE)
notableList<-tail(notableList, n = -12)
delist<-as.vector(head(notableList[,1],12))
delist<-c("Achn103051","Achn103061","Achn103071","Achn103081")

plotGenes(delist,24,120)

#colnames(plotdf)<-c("gene","time")  