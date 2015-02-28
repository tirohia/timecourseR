setwd("/home/ben/workspace/timeCourse/data")

library("plyr")

source("../code/R/tools.R")

datafile="/home/ben/workspace/timeCourse/data/actinidia/cornellGenomeModelTH-TCFrequencyMatrix.csv"
counts = read.table( datafile,sep=",", header=TRUE, row.names=1,stringsAsFactors = FALSE )

counts<-counts1
design <- getDataRange(counts,24,120)
design
datacolumns<-getDataCols(design)
counts<-counts[datacolumns,]
head(counts)
logcounts<-log(counts)

ddply(dat, .(id), summarise, var1 = var(x), var2 = var(y))
?ddply
var(counts)
?var
dim(counts)
for (i in 1:nrow(counts)){
  variance<-var(logcounts[i,])
  mean<-mean(logcounts[i,])
  if (mean >=10 && variance <=10){
    print (c(rownames(counts)[i],mean,variance))
  }
}
rownames(counts)[1]
var(logcounts[1,])
sum(counts[1,])

var(c))
?var