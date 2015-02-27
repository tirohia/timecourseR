source("http://bioconductor.org/biocLite.R")
biocLite("betr")
biocLite("GEOquery")
library("GEOquery")
library("betr")
library("affy")
source("../code/tools.R")

setwd ("/home/ben/workspace/timeCourse/data")

datafile<-file.choose()
rawCounts = read.table( datafile,sep=",", header=TRUE,row.names=1 , stringsAsFactors = FALSE )
rawCounts<-cleanDataNames(rawCounts)

design <- getDataRange(rawCounts,24,120)
datacolumns<-getDataCols(design)
counts<-rawCounts[,datacolumns]
scale<-scales[datacolumns]
datacolumns


M3 <- as.matrix(counts[1:39039,])
M3<-log(M3+0.000001)
M3<-scale(M3, center= FALSE, scale = scale)
M3<-scale(M3, center=FALSE, scale=colSums(M3))
length(scale)
head(M3)

counts<-counts+1
reps <- rep(c("rep1","rep2","rep3"),length(design[,1]))
conditions <- rep(design[,1],each=3)
time<- rep(design[,2],each=3)

prob <- betr(eset=M3, cond=conditions,timepoint=time, replicate=reps, alpha=0.01)

prob
delist<-names(head(sort(prob,decreasing=TRUE),12))
delist<-tail(delist, n = -12)
delist
problist<-prob[prob ==1]
length(problist)
problist<-prob[prob >= 0.99]
length(problist)
write.csv(problist,"betrDE.csv")
length(problist)

par(mfrow=c(1,1))
hist(prob, breaks=1000, col="skyblue", border="slateblue", main="")




##needed



geneName
normalize.loess(as.matrix(psadf[,2:4]))
thing<-rlm(log((rep1 + rep2 + rep3)/3) ~ time , psadf)
plot(thing)
with(psadf,{
  ok <- is.finite(time) & is.finite(rep1)
  x <- as.numeric(time)[ok]   ## you reproduce  the error without coercion 
  y <- as.numeric(rep1)[ok]
  data <- list(x = x, y = y)
  loess( y~x,data = data)
})
warnings()

dev.off()

notableList<-read.csv("timecourseNotableGenes.csv",header=FALSE)
notableList<-tail(notableList, n = -24)
delist<-as.vector(head(notableList[,1],12))
delist
warnings()



geneList<-c("Achn241691","Achn014741","Achn025281","Achn107721","Achn110911","Achn120901","Achn124291","Achn126231","Achn241691","Achn241691","Achn267561","Achn290511")
geneList<-c("Achn322901","Achn380061","Achn383541")

plotGenes(geneList,24,120)
