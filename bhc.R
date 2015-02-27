install.packages("som")
library("som")

M4<-filtering(M3, lt=20, ut=16000, mmr=3, mmd=200)
M5<-normalize(M4)
?som
foo <- som(M5, xdim=5, ydim=6)
foo <- som(M5, xdim=5, ydim=6, topol="hexa", neigh="gaussian")
par(mfrow=c(1,1))
plot(foo, sdbar=1, ylim=c(-4, 5), color=TRUE, ntik=3, yadj=0.1, xlab="a", ylab="b")
summary(foo)



################# Bayesian heirachical clustering
library("edgeR")
edgeRUsersGuide()
datafile="actinidiaReadCounts.tsv"
rawCounts = read.table( datafile,sep="\t", header=TRUE,row.names=1 , stringsAsFactors = FALSE )
rawCounts<-cleanDataNames(rawCounts)

head(rawCounts)
design <- getDataRange(rawCounts,48,72)
design <- design[design$condition =="psa",]


datacolumns<-getDataCols(design)
length(datacolumns)
counts<-rawCounts[,datacolumns]



group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=counts,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)


?DGEList
geneList<-data.frame(topTags(et,2000))
topTags(et)
write.csv(res,"edgeRDEGenes.csv")
class(geneList)
res=geneList[geneList$logFC < -2,]
dim(res)
