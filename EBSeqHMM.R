source("http://bioconductor.org/biocLite.R")
library("EBSeqHMM")

browseVignettes("EBSeqHMM")

design <- getDataRange(rawCounts,24,120,condition="psa")
datacolumns<-getDataCols(design)
length(datacolumns)
counts<-rawCounts[,datacolumns]

M3 <- as.matrix(counts[1:39039,]) #for actinidia

design

assay <- rep(c("rep1","rep2","rep3"),length(design[,1]))
trt <- rep(design[,1],each=3)
times<- rep(design[,2],each=3)
times
conditions <- factor(times, unique(time))
sizes <- MedianNorm(counts)
geneNormData <- GetNormalizedMat(counts, sizes)
geneNormData
PlotExp(geneNormData, conditions, Name="Achn138421")

EBSeqHMMGeneOut <- EBSeqHMMTest(Data=M3, sizeFactors=sizes, Conditions=conditions, UpdateRd=5)
GeneDECalls <- GetDECalls(EBSeqHMMGeneOut, FDR=.05)
toplist<-rownames(head(GeneDECalls,666))
head(GeneDECalls)

GeneDECalls[rownames(GeneDECalls)=="Achn053721"]
            [,2] >= 0.98]

write.csv(toplist,"EBSeqHMM.csv")
par(mfrow=c(4,5))
for (gene in toplist){
  PlotExp(geneNormData, conditions, Name=gene)
}
