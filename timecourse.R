source("http://bioconductor.org/biocLite.R")
library("timecourse")
library("preprocessCore")

## Select data source. 
setwd("/home/ben/workspace/timeCourse/data")

source("../code/R/tools.R")
#datafile="/home/ben/workspace/timeCourse/data/actinidia/cornellGenomeModelTH-TCFrequencyMatrix.csv"
datafile="/home/ben/workspace/timeCourse/data/psa/psaFrequencyMatrix.csv"
#datafile=file.choose()


counts = read.table( datafile,sep=",", header=TRUE,row.names=1 , stringsAsFactors = FALSE )
counts<-cleanDataNames(rawCounts)

## remove ribosomal entries from psa data. 
ribosomalEntriesfile="/home/ben/workspace/models/psa/ribosomalIYO"
ribosomalEntries = read.table( ribosomalEntriesfile,sep=",", header=FALSE,stringsAsFactors = FALSE )
counts<-counts[ ! rownames(counts) %in% ribosomalEntries, ] 
numberofGenes<- length(rownames(counts))

design <- getDataRange(rawCounts,24,120)
datacolumns<-getDataCols(design)
length(datacolumns)
counts<-counts[,datacolumns]

qnormalized<-normalize.quantiles(as.matrix(counts),copy=TRUE)

M3 <- as.matrix(qnormalized[1:numberofGenes,]) 
M3 <- log(M3+1)
sizeFactors
M3<-scale(M3, center=FALSE, scale=sizeFactors)
M3<-M3[complete.cases(M3), ]

assay <- rep(c("rep1","rep2","rep3"),length(design[,1]))
trt <- rep(design[,1],each=3)
time<- rep(design[,2],each=3)


size <- matrix(3, nrow=numberofGenes, ncol=2) # actinidia

## sanity check
length(trt)
length(assay)
length(time)
dim(size)
dim(counts)
dim(M3)
counts
trt
assay
time
datacolumns
size
design

## and go ...
MB.paired <- mb.long(M3, method="paired", times=5, reps=size, condition.grp=trt, rep.grp=assay, time.grp=time)
MB.paired1 <- mb.long(M3, method="paired", times=5, reps=size)
MB.2D <- mb.long(M3, method="2", times=5, reps=size, condition.grp=trt, rep.grp=assay,time.grp=time)

summary(MB.paired)
297mm x 210mm

genenames=rownames(rawCounts)


##################### BY ranking
## Remember to change times argument denpending on how many time points there are. 
png(file="../results/timecourse1st24HActGenes1-12.png" ,width=1050,height=1485)#,width=1920,height=1080)#
par(mfrow=c(4,3))
for (cycle1 in 13: 24) {
  plotProfile(MB.paired, type="b", gnames=genenames,  pch=c("1","2","3"), xlab="Hour",ranking=cycle1,legloc=c(1000,1000),col=c("green","blue"),lwd=1,cex.main=1.5)
}
dev.off()
source("../code/R/plotProfile.R")

################ By gene id.
delist<-c("IYO_14535","IYO_13640","IYO_18340","IYO_23120","IYO_03285","IYO_05360","IYO_13220","IYO_13640","IYO_17550","IYO_06520","IYO_27715","IYO_27720")

delist<-c("Achn160161","Achn061231","Achn253571","Achn033321","Achn229111","Achn309921","Achn334481","Achn136171","Achn126231","Achn143171","Achn334501","Achn334541","Achn229071","Achn276611","Achn056131")
notableList<-read.table("../results/tripleDE.csv",sep=",")
delist<-notableList[,1]
notableList<-read.csv("timecourseNotableGenes.csv",header=FALSE)
notableList<-tail(notableList, n = -24)
delist<-as.vector(head(notableList[,1],12))

png(file="../results/timecourse24HrcycleActGenesOfInterest25-31.png" ,width=1050,height=1485)#,width=1920,height=1080)#

par(mfrow=c(4,4))
for (geneId in delist) {
  plotProfile(MB.paired, type="b", gnames=genenames,  pch=c("1","2","3"), xlab="Hour",gid=geneId,legloc=c(1000,1000),col=c("green","blue"),lwd=1,cex.main=1.5)
}
dev.off()

?plotProfile
#



### Do a for loop and print out when == Test
?plotProfile
for(cycle1 in 1:1000){
  if (rownames([cycle1,]=="test"){ ## MB.Paired thing. 
    print(cycle1)
  }
}

dev.off()

source("../code/plotProfile.R")


png(file="timecourseGenes2D1-20.png" ,width=1920,height=1080)

par(mfrow=c(1,1))
for (cycle1 in 1: 20) {
  plotProfile(MB.2D, type="b", gnames=genenames,  pch=c("1","2","3"), xlab="Hour",ranking=1)
}
dev.off()

plotProfile(MB.2D,type="b", gnames=genenames,ranking=1)
plotProfile(MB.multi, stats="MB", type="b")





par(mfrow=c(1,1))
for (cycle1 in 1: 50) {
  qqnorm(rawCounts[cycle1,])
  qqline(rawCounts[cycle1,])
}
cor(rawCounts)
cov(rawCounts)
pairs(rawCounts)
e.cor<-eigen(cor(rawCounts))
e.cor
e.cov<-eigen(cov(rawCounts))
e.cov












################### viginette stuff ###############################
head(M3)
head(M2)
M2 <- matrix(rep(14,39044*30), ncol=30)
M2[1:20,] <- t(apply(M2[1:20,],1,sim.data2))
M2[21:39044,] <- t(apply(M2[21:39044,],1,sim.data2, 0)) 
M2

trt <- rep(c("wt","mt"),each=15)
assay <- rep(rep(c("rep1","rep2","rep3"),each=5),2)
size <- matrix(3, nrow=39044, ncol=2)
MB.paired <- mb.long(M2, method="paired", times=5, reps=size, condition.grp=trt, rep.grp=assay)




################## For the simulation of data #########################
SS <- matrix(c( 1e-02, -8e-04, -0.003,  7e-03,  2e-03,
                -8e-04,  2e-02,  0.002, -4e-04, -1e-03,
                -3e-03,  2e-03,  0.030, -5e-03, -9e-03,
                7e-03, -4e-04, -0.005,  2e-02,  8e-04,
                2e-03, -1e-03, -0.009,  8e-04,  7e-02), ncol=5)

library("timecourse")
sim.Sigma <- function(){
  S <- matrix(rep(0,25),ncol=5)
  x <- mvrnorm(n=10, mu=rep(0,5), Sigma=10*SS)
  for(i in 1:10)
    S <- S+crossprod(t(x[i,]))
  
  solve(S)
}

sim.data2 <- function(x, indx=1){
  mu <- rep(runif(1,8,x[1]),5)
  if(indx==1)
    res <- c(as.numeric(t(mvrnorm(n=3, mu=mu+rnorm(5,sd=5), Sigma=sim.Sigma()))),
             as.numeric(t(mvrnorm(n=3, mu=mu+rnorm(5,sd=3.2), Sigma=sim.Sigma()))))
  
  if(indx==0) res <- as.numeric(t(mvrnorm(n=6, mu=mu+rnorm(5,sd=3), Sigma=sim.Sigma())))
  res 
}

M3<-rbind(M3,M2)
dim(M3)
rownames(M2)<-c(rep("test",100))
head(M2)
M2 <- matrix(rep(14,1000*30), ncol=30)
M2 <- matrix(rep(14,100*30), ncol=30)
M2[1:1000,] <- t(apply(M2[1:1000,],1,sim.data2))
M2[21:1000,] <- t(apply(M2[21:1000,],1,sim.data2, 0))
dim(M2)

head(M2,30)
