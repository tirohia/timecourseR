source("http://bioconductor.org/biocLite.R")
library("preprocessCore")
counts
countsMatrix<-as.matrix(counts)
qnormalized<-normalize.quantiles(countsMatrix,copy=TRUE)
head(qnormalized)
