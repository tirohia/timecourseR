library("cycle")
library("Biobase")
library("caret")
setwd ("/home/ben/workspace/timeCourse/data")
datafile<-file.choose()
networkCounts = read.table( datafile,sep=",", header=TRUE,row.names=1 , stringsAsFactors = FALSE )
networkCounts<-cleanDataNames(networkCounts)
design <- getDataRange(rawCounts,24,120)
datacols<-getDataCols(design)
networkCounts<-networkCounts[,datacols]  

#phenotype data for annotated data frame. 
pDataFile <- file.path("pdata.csv")
pData <- read.csv(pDataFile,header=TRUE, sep=",")
phenoData <- new("AnnotatedDataFrame",data=pData, varMetadata=metadata)

#annotated data frame. 
metadata <- data.frame(labelDescription=c("datacol","time","trt","assay"))
colnames(networkCounts)<-rownames(pData)
all(rownames(pData) == colnames(networkCounts))
phenoData <- new("AnnotatedDataFrame", data = pData, varMetadata = metadata)

pData<-data.frame(c(24,48,72,96,124,24,48,72,96,124))
colnames(pData)<-c("datacol","time","trt","assay")
pData
metadata

## Experiment data - complete bollocks atm
experimentData <- new("MIAME",name="Ben Curran", lab="No lab in particular", contact="b.curran@auckland.ac.nz")
?experimentData

#Matrix
M3<-as.matrix(networkCounts)
M3<- M3+1
colnames(M3)<-rownames(pData)
## Expression Set
exampleSet <- ExpressionSet(assayData=M4, phenoData=phenoData,experimentData=experimentData,annotation="hgu95av2")
exampleSet
cycleresults<-ar1analysis(exampleSet)

head(as.vector(cycleresults))
## seriously. What does this do? sigh. Ah. takes every 3 and aggregates? 
M4<-NULL
for (i in 1:39039){
    #print(M3[,i])
    insert<-aggregate(x=(M3[i,]),by = list(index),mean)
    #insert<-c(2,2,2)
    insert.rownames<-NULL
    insert<-insert[,2]
    M4 <- rbind(M4,t(insert))  
}
rownames(M4)<-rownames(networkCounts)
head(M4)
M4psa<-M4[,1:5]
exampleSet<-standardise(exampleSet)

auto.corr <- 0
for (i in 2:dim(exprs(exampleSet))[[2]]){
  auto.corr[i]<-cor(exprs(exampleSet)[,i-1],exprs(exampleSet)[,i])
}
auto.corr

ffdrfourier()
data(yeast)
yeast <- yeast[1:200,]

index <- rep(1:10,3)
index<-sort(index)
zeroVar<-nzv(t(M3), freqCut = 80/20, uniqueCut = 10, saveMetrics = FALSE)
length(zeroVar)
?cor
head(M5)
rownames(M4)
zeroVar
!rownames(M4) %in% zeroVar
M3<-M3[-zeroVar,]
dim(M4Psa)
