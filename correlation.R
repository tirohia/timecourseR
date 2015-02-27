setwd ("/home/ben/workspace/timeCourse/data")
source("../code/R/tools.R")
library("preprocessCore")


datafile="cornellGenomeModelTH-TCFrequencyMatrix.csv"
psadatafile="psaGenome-TCFrequencyMatrix.csv"

## load data into rawCounts
rawCounts = read.table( datafile,sep=",", header=TRUE,row.names=1 , stringsAsFactors = FALSE )
psaCounts1 = read.table( psadatafile,sep=",", header=TRUE,row.names=1 , stringsAsFactors = FALSE )

rawCounts<-cleanDataNames(rawCounts)
psaCounts1<-cleanDataNames(psaCounts1)

psaNames<-rownames(psaCounts1)

design <- getDataRange(rawCounts,24,120)
datacolumns<-getDataCols(design)
length(datacolumns)
counts<-rawCounts[,datacolumns]
counts<-normalize.quantiles(as.matrix(rawCounts),copy=TRUE)
psaCounts<-normalize.quantiles(as.matrix(psaCounts),copy=TRUE)

row1<-design[design$condition == "psa",]$X1
row2<-design[design$condition == "psa",]$X2
design[design$condition == "psa",]$X3

row<-rawCounts[,row1]
row
row2<-rawCounts[,row1]
row2
cor(row,row2)

rows<-c("X1","X2","X3")
cormat<-data.frame(row.names=rows)
rawCounts<- rawCounts[rowSums(rawCounts[, -1])>0, ]
rawCounts<- log(rawCounts+1)
psaCounts<- psaCounts[rowSums(psaCounts[, -1])>0, ]
psaCounts<- log(psaCounts+1)

i<-6
j<-1
k<-"X1"
l<-"X3"
row1<-design[design$condition == "psa",k]
row2<-design[design$condition == "psa",l]
row1
row2
datarow1<-rawCounts[i,row1]
datarow2<-psaCounts[j,row2]
datarow1
datarow2
cor(datarow1,datarow2)
cormat[k,l]<-cor(datarow1,datarow2)
cormat[k,l]
cormat16<-cormat

is.na(sum(rowSums(cormat16)))

mean(unlist(cormat),na.rm=TRUE)
min(cormat)
max(cormat)

rowSums(cormat)
is.na(cormat)
warnings()
is.na(sum(rowSums(cormat))) 
mean(unlist(cormat),na.rm=TRUE)
?rank
?rownames
counts[1,]
rawCounts[1,]
length(rownames(counts))
length(rownames(rawCounts))
rownames(psaCounts)<-psaNames
rownames(counts)<-rownames(rawCounts)

for (i in rownames(counts)){
  for (j in rownames(psaCounts)){
    cormat<-data.frame(row.names=rows)
    for (k in rows){
      for (l in rows){
        row1<-design[design$condition == "psa",k]
        rowdata<-rawCounts[1,row1]
        if(sum(rowdata) == 0 ) {next}
        row2<-design[design$condition == "psa",l]
        row2<-psaCounts[j,row2]
        if(sum(row2) == 0 ) {next}
        cormat[k,l]<- cor(row1,row2)
      }
    }
    if(sum(rowSums(cormat,na.rm=TRUE)) != 0){
      if( mean(unlist(cormat),na.rm=TRUE) >= 0.9){
        i<-3
        rownames(rawCounts)
        counts[i,3]
        print(rownames(counts[i,]))
        print(j)
        print(mean(unlist(cormat),na.rm=TRUE))
        print(cormat)
      }
      else if(mean(unlist(cormat),na.rm=TRUE) <= -0.9){
        print(i)
        print(j)
        print(mean(unlist(cormat),na.rm=TRUE))
        print(cormat)
      }
    }
  }
}
cormat
?data.frame
