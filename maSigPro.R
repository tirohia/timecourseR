library("maSigPro")
library("limma")
v <- voom(counts1,design,plot=TRUE)
v

plotMDS(v,xlim=c(-2.5,2.5))
class(v)
as.matrix(v)



fit <- eBayes(lmFit(v,design))
topTable(fit,coef=2)


counts1<-counts
colnames(counts1)<-conditions
par(mfrow=c(1,1))
?voom
head(counts)
conditions<-factor(c(rep("C",15),rep("P",15)))
f <- factor(counts$genes, levels=conditions)
design <- model.matrix(~0+f) 


library("Rsubread")
library(limma)
library(edgeR)

options(digits=2)
targets <- readTargets()
celltype <- factor(targets$CellType)
design <- model.matrix(~conditions)
design
targets