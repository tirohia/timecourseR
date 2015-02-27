source("http://bioconductor.org/biocLite.R")
library("Biobase")
library("DESeq2")
library("gplots")
library("vsn")



## Selcet data source. 
setwd("/home/ben/workspace/timeCourse/data")

source("../code/R/tools.R")
%datafile="/home/ben/workspace/timeCourse/data/actinidia/cornellGenomeModelTH-TCFrequencyMatrix.csv"
datafile="/home/ben/workspace/timeCourse/data/psa/psaFrequencyMatrix.csv"
#datafile=file.choose()


## Load data into DESeq
## Clean up the column names, select the data to be used

counts = read.table( datafile,sep=",", header=TRUE, row.names=1,stringsAsFactors = FALSE )
counts <-cleanDataNames(counts)

## remove ribosomal entries from psa data. 
ribosomalEntriesfile="/home/ben/workspace/models/psa/ribosomalIYO"
ribosomalEntries = read.table( ribosomalEntriesfile,sep=",", header=FALSE,stringsAsFactors = FALSE )
counts<- counts[ ! rownames(counts) %in% ribosomalEntries, ] 


## run a D.E. analysis.
design <- getDataRange(rawCounts,24,120,"psa")
design
datacolumns<-getDataCols(design)
counts<-counts[,datacolumns]
colData <- data.frame(row.names=colnames( counts ), time=as.factor(rep(design$time,each=3)))
dds <- DESeqDataSetFromMatrix(countData = counts,colData = colData, design =~ time)
dds <- DESeq(dds)
plotDispEsts(dds, main="Modelling dispersion for gene shrinkage")
sizeFactors<-sizeFactors(dds)


### Heteroscedacity control - for exporting data with the mean-variance relationship smoothed. 

notAllZero<-(rowSums(counts(dds))>0)
##prior to controlling

##relative log control
rld<-rlog(dds,blind=FALSE)
#varaince stabilising control
vsd<-varianceStabilizingTransformation(dds)

## plots of SD-mean relationship.
par(mfrow=c(1,3))
meanSdPlot(log2(counts(dds)[notAllZero,]+1))
meanSdPlot(assay(rld[notAllZero]))
meanSdPlot(assay(vsd[notAllZero]))


## Controlling for library size for export of data
## check this on sample data. 
sizeFactors<-sizeFactors(dds)
counts1<-as.matrix(counts) %*% diag(sizeFactors)
counts1


## get data and design for desired comparson. This is for a mass comparison.
wrapper(24,120,"down")
row.names=colnames( counts )
time=as.factor(design$time)
colData <- data.frame(row.names=colnames( counts ), time=as.factor(design$time))
deExp(counts,design,"up")
dds<-deExp(counts,design,"up",banded="TRUE")

wrapper <- function(time1,time2,direction){
  counts = read.table( datafile,sep=",", header=TRUE, row.names=1,stringsAsFactors = FALSE )
  counts <-cleanDataNames(counts)
  counts<-getRawData(counts,time1,time2,"psa")
  design <- getDesign(counts,time1,time2,"psa")
  deExp(counts,design,direction)
}

## prepare results matrix - this is just the number of DE genes, not their names/read counts etc. 
geneMapUp<-matrix(data=0,nrow=length(time1),ncol=length(time2))
geneMapDown<-matrix(data=0,nrow=length(time1),ncol=length(time2))
rownames(geneMapUp)<-time2
colnames(geneMapUp)<-time1
rownames(geneMapDown)<-time2
colnames(geneMapDown)<-time1

time1<-c(24,48,72,96,120)
time2<-time1

for(j in 1:length(time1)){
  for(i in 1:length(time2)){
    print (c(i,j))
    print(c(time1[i],time2[j]))
    resultCount<-nrow(wrapper(time1[i],time2[j],"up"))
    if (is.integer(resultCount)){
      geneMapUp[i,j]<-resultCount
    }
    else{
      geneMapUp[i,j]<-0
    }
    resultCount<-nrow(wrapper(time1[i],time2[j],"down"))
    if (is.integer(resultCount)){
      geneMapDown[i,j]<-resultCount
    }
    else{
      geneMapDown[i,j]<-0
    }
    #geneMapUp[i,j] <- nrow(wrapper(time1[i],time2[j],"up"))
    #geneMapDown[i,j] <- nrow(wrapper(time1[i],time2[j],"down"))
   }
}
geneMapUp
geneMapDown
write.csv(geneMapUp,"psageneMapLogUp.csv")
write.csv(geneMapDown,"psageneMapLogDown.csv")

palette<-brewer_pal(pal="Greys")(n=9)

heatmap.2( geneMapUp,col = palette,margins = c(5, 10),trace = "none",lhei = c(2, 8),xlab = "Comparison", na.color="blue", main = "psa DE genes, Up",Colv = FALSE,Rowv =FALSE)
heatmap.2( geneMapDown,col = palette,margins = c(5, 10),trace = "none",lhei = c(2, 8),xlab = "Comparison", na.color="blue", main = "psa DE genes, Down",Colv = FALSE,Rowv =FALSE)

library(RColorBrewer)
library("scales")



######################### Random tootling ##############################
datafile<-file.choose()
psaGeneMapUp <-  read.csv( datafile, header=TRUE, row.names=1,stringsAsFactors = FALSE )
psaGeneMapUp

datafile<-file.choose()
actInocGeneMapUp <-  read.csv( datafile, header=TRUE, row.names=1,stringsAsFactors = FALSE )
actInocGeneMapUp


dds4872
res <- results(dds4872)
plotMA(res)
hist(res$padj, breaks=100, col="skyblue", border="slateblue", main="")
res = res[ res$log2FoldChange < -2, ]
res <- na.omit(res)
res=res[res$padj < 0.1,]
res

file="deResults/deIntersections.csv"
intersections=read.csv(file,row.names=1)
intersections
rownames(intersections)
inter<-log(intersections+1)

heatmap.2( as.matrix(inter),col = palette,margins = c(5, 10),trace = "none",lhei = c(2, 8),xlab = "Comparison", na.color="blue", main = "psa DE genes, Up",Colv = FALSE,Rowv =FALSE)
?heatmap.2
)
########################################################################



colData <- data.frame(row.names=colnames( counts ), condition=as.factor(design$condition),time=as.factor(design$time))
colData <- data.frame(row.names=colnames( rawCounts ), time=as.factor(design$time))
colData


use <- res$baseMean > attr(res,"filterThreshold")
table(use)
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE, col = colori, space = 0, main = "", ylab="frequency")
