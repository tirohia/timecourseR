inputfile="/home/ben/workspace/timeCourse/code/nail/input/sample_data/sample_expression_data.txt"
nailIn = read.table( inputfile,sep="\t", header=TRUE, row.names=1,stringsAsFactors = FALSE )



nailIn
chip2<-nailIn[,2]/nailIn[,1]
log2(chip2)

sums
sums<-colSums(nailIn)

chip1 <-nailIn[,1]/sums[1]
log2(chip1)
log2(nailIn[,1]/(sum(sums)/6))

log2(nailIn[,1])

log2(nailIn[,2]/sums[2])

log2(nailIn[1,]/(sum(nailIn[1,])/6))
rld<-rlog(dds)
rld
library("vsn")
par(mfrow=c(1,3))
notAllZero <- (rowSums(counts(dds))>0)
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1))
meanSdPlot(assay(rld[notAllZero,]))
mcols(rld,use.names=TRUE)[1:4,1:4]

transformedData<-assays(dds)[["mu"]]
head(assays(dds)[["mu"]])

class(transformedData)
dim(transformedData)
head(counts)
write.csv(transformedData, "nailInput.csv" )

transformedData<-transformedData[candidates,]

candidates<-read.csv("../results/tripleDE.csv",row.names=1,header=TRUE)
candidates<-rownames(candidates)
candidates

