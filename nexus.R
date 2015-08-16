setwd("/home/ben/workspace/timeCourse/data")
source("../code/R/tools.R")
source(("../code/R/voomlimma.R"))
source(("../code/R/edgeR.R"))
source(("../code/R/timecourse.R"))

outputImages<-"/home/ben/Google Drive/PhD/writing/Timecourse/images/"
outputDir<-"/home/ben/workspace/timeCourse/data/deResults/"

subject<-"plant"  #should be plant or bacteria
track<-"control"  #should be psa or control, i.e. treatment
counts<-load.data(subject,track,0,120)
edgeRanalysis(counts,subject,track,0,120,outputImages,paste(outputDir,"edgeR/",sep=""))
limmaAnalysis(counts,subject,track,0,120,outputImages,paste(outputDir,"limma/",sep=""))

subject<-"plant"  #should be plant or bacteria
track<-"psa"  #should be psa or control, i.e. treatment
counts<-load.data(subject,track,0,120)
edgeRanalysis(counts,subject,track,0,120,outputImages,paste(outputDir,"edgeR/",sep=""))
limmaAnalysis(counts,subject,track,0,120,outputImages,paste(outputDir,"limma/",sep=""))

subject<-"bacteria"  #should be plant or bacteria
track<-"psa"  #should be psa or control, i.e. treatment
counts<-load.data(subject,track,0,120)
edgeRanalysis(counts,subject,track,0,120,outputImages,paste(outputDir,"edgeR/",sep=""))
limmaAnalysis(counts,subject,track,0,120,outputImages,paste(outputDir,"limma/",sep=""))

subject<-"plant"
track<-"all"
counts<-load.data(subject,track,24,120)
design<-getDataRange(24,120)
timecourseAnalysis(counts,design)

subject<-"bacteria"
track<-"all"
counts<-load.data(subject,track,1.5,120)
head(counts)
design<-getDataRange(1.5,120)
timecourseAnalysis(counts,design)


## counts are loaded