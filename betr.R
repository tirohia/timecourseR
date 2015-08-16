source("http://bioconductor.org/biocLite.R")
library("GEOquery")
library("betr")
library("affy")

source("../code/R/tools.R")
setwd("/home/ben/workspace/timeCourse/data")


subject<-"bacteria"  #should be plant or bacteria
track<-"psa"  #should be psa or control, i.e. treatment

## can't go to zero when looking at bacteria reads - 0 time is just a copy of the plant zero timepoint - bacterial reads 
## are all zero, which screws with DESeq. So go from 1.5 instead. 
counts<-load.data(subject,track,1.5,120) 
design <- getDataRange(1.5,120,track)
colData <- data.frame(row.names=colnames(counts), time=as.factor(rep(design$time,each=3)),condition=as.factor(rep(design$condition,each=3)))
dds <- DESeqDataSetFromMatrix(countData = counts,colData = colData, design =~ time )
dds <- DESeq(dds, test="Wald")
rld<-rlog(dds,blind=FALSE)

counts<-assay(rld) #rld has to come from a deseq analysis of the same experimental design. 

reps <- rep(c("rep1","rep2","rep3"),length(design[,1]))
conditions <- rep(design[,1],each=3)
time<- rep(design[,2],each=3)

prob <- betr(eset=as.matrix(counts), cond=conditions,timepoint=time, replicate=reps, alpha=0.01)

problist<-prob[prob ==1]
length(problist)
problist<-prob[prob >= 0.99]
length(problist)

write.csv(problist,"deResults/betr/psaGenes.csv")

hist(prob, breaks=100, col="skyblue", border="slateblue", main="")

##needed ##why? 
#normalize.loess(as.matrix(psadf[,2:4]))
#thing<-rlm(log((rep1 + rep2 + rep3)/3) ~ time , psadf)
#plot(thing)
#with(psadf,{
#  ok <- is.finite(time) & is.finite(rep1)
#  x <- as.numeric(time)[ok]   ## you reproduce  the error without coercion 
#  y <- as.numeric(rep1)[ok]
#  data <- list(x = x, y = y)
#  loess( y~x,data = data)
#})
#warnings()

#dev.off()

#notableList<-tail(notableList, n = -24)
#delist<-as.vector(head(notableList[,1],12))
#delist

#plotGenes(geneList,24,120)



