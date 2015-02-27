source("http://bioconductor.org/biocLite.R")
library("Biobase")
library("DESeq2")
library("gplots")
library("vsn")



## Selcet data source. 
setwd("/home/ben/workspace/timeCourse/data")

source("../code/R/tools.R")
#datafile="/home/ben/workspace/timeCourse/data/actinidia/cornellGenomeModelTH-TCFrequencyMatrix.csv"
datafile="/home/ben/workspace/timeCourse/data/psa/psaFrequencyMatrix.csv"
#datafile=file.choose()

## Load data into DESeq
## Clean up the column names, select the data to be used

counts = read.table( datafile,sep=",", header=TRUE, row.names=1,stringsAsFactors = FALSE )
counts <-cleanDataNames(counts)
design <- getDataRange(rawCounts,24,120,"psa")
datacolumns<-getDataCols(design)
counts<-counts[,datacolumns]
colData <- data.frame(row.names=colnames( counts ), time=as.factor(rep(design$time,each=3)))
dds <- DESeqDataSetFromMatrix(countData = counts,colData = colData, design =~ time)
dds <- DESeq(dds)


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


dds<-deExp(counts,design,dds=TRUE)
time1 <- c(24,48,72,96,120)
time2 <- time1

wrapper <- function(time1,time2,direction){
  counts = read.table( datafile,sep=",", header=TRUE, row.names=1,stringsAsFactors = FALSE )
  counts <-cleanDataNames(counts)
  counts<-getRawData(counts,time1,time2)
  design <- getDesign(counts,time1,time2)
  deExp(counts,design,direction)
}

## prepare results matrix - this is just the number of DE genes, not their names/read counts etc. 
geneMapUp<-matrix(data=0,nrow=length(time1),ncol=length(time2))
geneMapDown<-matrix(data=0,nrow=length(time1),ncol=length(time2))
rownames(geneMapUp)<-time2
colnames(geneMapUp)<-time1
rownames(geneMapDown)<-time2
colnames(geneMapDown)<-time1



for(j in 1:length(time1)){
  for(i in 1:length(time2)){
    print (c(i,j))
    print(c(time1[i],time2[j]))
    geneMapUp[i,j] <- wrapper(time1[i],time2[j],"up")
    geneMapDown[i,j] <- wrapper(time1[i],time2[j],"down")
   }
}

write.csv(geneMapUp,"psageneMapLogUp.csv")
write.csv(geneMapDown,"psageneMapLogDown.csv")

palette<-brewer_pal(pal="Greys")(n=9)

heatmap.2( geneMapUp,col = palette,margins = c(5, 10),trace = "none",lhei = c(2, 8),xlab = "Comparison", na.color="blue", main = "psa DE genes, Up",Colv = FALSE,Rowv =FALSE)
heatmap.2( geneMapDown,col = palette,margins = c(5, 10),trace = "none",lhei = c(2, 8),xlab = "Comparison", na.color="blue", main = "psa DE genes, Down",Colv = FALSE,Rowv =FALSE)

library(RColorBrewer)
library("scales")


dds4872<-wrapper(48,72,"down")
dds7296<-wrapper(72,96,"up")

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


rawCounts <- counts #everything
colnames( rawCounts )
tail(rawCounts)

row.names = colnames( rawCounts )
row.names



#design for 24hr cycle

timeTable = read.table( "lookupsheet.csv",sep=",", header=TRUE, stringsAsFactors = FALSE )
timesLookup <- timeTable[timeTable$time >= 24, ] 
timesLookup 
rawCounts <- counts[,timesLookup] #24hr cycle timepoints
time <- timeTable[,2]
condition <- timeTable[,1]
condition
as.vector()
rep(timeTable["condition"],3)
design = data.frame(
  row.names = colnames( counts ),
  condition=c("Control","Psa","Control","Control","Psa","Control","Psa","Control","Psa","Control","Psa","Control","Control","Control","Control","Control","Psa","Psa","Psa","Psa","Control","Control","Control","Control","Psa","Psa","Psa","Control","Control","Control","Psa","Psa","Psa"),
  time=c("72","24","48","24","72","72","72","48","48","48","72","0","0","0","24","24","24","24","48","48","72","120","120","120","120","120","120","96","96","96","96","96","96"),
  libtype = c("paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end")
)
condition <- rep(design[,1],each=3),
libtype=c(rep("paired-end",30)),
time<- rep(design[,2],each=3)

#design for just control time points
design = data.frame(
  row.names = colnames(rawCounts),
  time= c("12","72","48","24","72","48","48","0","0","0","1.5","1.5","1.5","3","3","3","6","6","6","12","12","24","24","72","120","120","120","96","96","96"),
  condition= c("Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control"),
  libType = c("paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end")
)

colnames(rawCounts)
# design for everything combined. 
design = data.frame (
  row.names = colnames (rawCounts),
  time= c(12,72,24,48,24,72,72,72,48,48,48,72,0,0,0,1.5,1.5,1.5,1.5,1.5,1.5,3,3,3,3,3,3,6,6,6,6,6,6,12,12,12,12,12,24,24,24,24,48,48,72,120,120,120,120,120,120,96,96,96,96,96,96),
  condition= c("Control","Psa","Control","Psa","Control","Control","Control","Psa","Psa","Psa","Control","Control","Control","Control","Control","Control","Control","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Control","Control","Control","Psa","Control","Control","Psa","Psa","Control","Control","Control","Psa","Psa","Psa","Psa","Control","Control","Control","Control","Psa","Control","Psa","Psa","Control","Control","Control","Psa","Psa","Psa","Psa","Control","Psa","Control"),
  libType = c("paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end")
)

?DESeq
colData <- data.frame(row.names=colnames( counts ), condition=as.factor(design$condition),time=as.factor(design$time))
colData <- data.frame(row.names=colnames( rawCounts ), time=as.factor(design$time))
colData

head(counts)
## time reduced LRT
dds <- DESeqDataSetFromMatrix(countData = counts,colData = colData, design =~  time)
dds <- DESeq(dds)
sizeFactors<-sizeFactors(dds)

res <- results(dds7296)
par(mfrow=c(1,1))

plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
res <- na.omit(res)
res = res[ res$log2FoldChange > 2, ]
res = res[ res$log2FoldChange < -2, ]
res=res[res$padj < 0.1,]
res

use <- res$baseMean > attr(res,"filterThreshold")
table(use)
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE, col = colori, space = 0, main = "", ylab="frequency")



write.csv(as.data.frame(res[order(res$padj),]),file="../results/cornell72-96_timecourse_up.csv")
colData

## condition and time, reduced to condition. 
dds <- DESeqDataSetFromMatrix(countData = rawCounts,colData = colData, design =~  condtion)
dds <- DESeq(dds, test="LRT", full=~time+condition, reduced=~time)
res <- results(dds)
plotMA(res)
hist(res$padj, breaks=100, col="skyblue", border="slateblue", main="")
res = res[ res$log2FoldChange < -2, ]
res <- na.omit(res)
res=res[res$padj < 0.1,]
res





ddsLRT<- dds
## run DESeq
dds <- DESeq(dds)
ddsCtrst <- dds[, dds$condition == "Psa"]
colData(ddsCtrst)
ddsCtrst
as.data.frame(colData(ddsCtrst)[,c("0","42")])



ddsLRT <- nbinomLRT(dds, reduced = ~ condition)
resLRT <- results(ddsLRT)


write.csv(as.data.frame(res[order(res$padj),]),file="pfr_manual_models_timecourse_down.csv")
resLRT <- na.omit(resLRT)
res = resLRT[ resLRT$log2FoldChange < -2, ]
res <- na.omit(res)
res=res[resLRT$padj < 0.1,]
res


abline(h=c(2,-2))


dds

## all of the data types generated.
substr(names(mcols(dds)),1,10)

dev.off()
#qucik sanity check
plotDispEsts(dds)

# Dissection of results
compare(dds,"48","72" )
resultSet <- results(dds, contrast=c("condition","Control","Psa"))
res=na.omit(resultSet)
abline(h=c(2,2))
plotMA(res)



compare <- function(resultDataSet,condition1,condition2) {
  resultSet <- results(resultDataSet, contrast=c("condition",condition1,condition2))
  res=na.omit(resultSet)
  par(mfrow=c(1,1))
  
  abline(h=c(2,2))
  
  plotMA(res)
  
  res = res[ res$log2FoldChange > 2, ]
  res=res[res$padj < 0.1,]
  res
  plotMA(res)
#hist(resSig$pval, breaks=100, col="skyblue", border="slateblue", main="")
#write.csv(as.data.frame(resSig[order(resSig$padj),]),file="pfr_manual_models_desq2_sdw_psa_results_up.csv")
  res=na.omit(resultSet)
  res = res[ res$log2FoldChange < -2, ]
  res=res[res$padj < 0.1,]
  res
  plotMA(res)
#hist(resSig$pval, breaks=100, col="skyblue", border="slateblue", main="")
#write.csv(as.data.frame(resSig[order(resSig$padj),]),file="pfr_manual_models_desq2_sdw_psa_results_down.csv")



}








?plotDispEsts

res <- results(dds)
resOrdered <- res[order(res$padj),]
head(resOrdered)
res=na.omit(res)
res



res <- results(dds)
plotMA(res)
res

genesetLists <- resMFCondition["stat"]
genesetLists
write.table(genesetLists, file = "geneSets",append = FALSE, sep = "\t",eol="\n",row.names=TRUE,col.names=TRUE)

            
plotMA(res, main="DESeq2", ylim=c(-5,5))
plotMA(dds, ylim=c(-5,5),main="DESeq2")
?plotMA
countData <- counts(pasillaGenes)
colData <- pData(pasillaGenes)[,c("condition","type")]
dds <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~ condition)



dds <- makeExampleDESeqDataSet()

