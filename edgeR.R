library("edgeR")

## Selcet data source. 
setwd("/home/ben/workspace/timeCourse/data")

source("../code/R/tools.R")
datafile="/home/ben/workspace/timeCourse/data/actinidia/cornellGenomeModelTH-TCFrequencyMatrix.csv"
datafile="/home/ben/workspace/timeCourse/data/psa/psaFrequencyMatrix.csv"
#datafile=file.choose()


## Load data into DESeq
## Clean up the column names, select the data to be used

counts = read.table( datafile,sep=",", header=TRUE, row.names=1,stringsAsFactors = FALSE )
counts <-cleanDataNames(counts)

## remove ribosomal entries from psa data. 
ribosomalEntriesfile="/home/ben/workspace/models/psa/ribosomalIYO"
ribosomalEntries = read.table( ribosomalEntriesfile,sep=",", header=FALSE,stringsAsFactors = FALSE )
counts<-counts[ ! rownames(counts) %in% ribosomalEntries, ] 


head(counts)
design <- getDataRange(counts,24,120)
design <- design[design$condition =="psa",]
design <- design[c(1,3),]
design

datacolumns<-getDataCols(design)
length(datacolumns)
counts<-counts[,datacolumns]

group <- factor(c(1,1,1,2,2,2))
group <- factor(c(rep(design$time,each=3)))

y <- DGEList(counts=counts,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)

geneList<-data.frame(topTags(et,200))
topTags(et)
write.csv(res,"edgeRDEGenes.csv")
class(geneList)
res=geneList[geneList$logFC < -2,]
res=geneList[geneList$logFC > 1.5,]
dim(res)

res

##


dim(counts)
head(counts)
cplist<-cpm(counts)
head(cplist)
countlist<-cplist[list,]
dim(countlist)
head(countlist)


####






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