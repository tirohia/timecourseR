source("http://bioconductor.org/biocLite.R")
setwd("/home/ben/workspace/timeCourse/data")
library("TDARACNE")

subject<-"plant"  #should be plant or bacteria
track<-"control"  #should be psa or control, i.e. treatment

### raw data. 
counts<-load.data(subject,track,0,120)
##normalized data
counts<-assay(rld)
design<-getDataRange(0,120,track)

## the genes of interest. 
datafile<-"/home/ben/workspace/timeCourse/data/deResults/networkGenesMedians.csv"
networkGenes <- read.csv(datafile,row.names=1)
counts<-counts[rownames(networkGenes),]


##Expression set
eset<-createESet(counts,design)
aracneGraph<-TDARACNE(eset,11,"netIRMAon",delta=3,likehood=1.2,norm=2,logarithm=1,thresh=0,ksd=0,0.15);
