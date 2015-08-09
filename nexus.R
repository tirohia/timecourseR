subject<-"plant"  #should be plant or bacteria
track<-"control"  #should be psa or control, i.e. treatment



setwd("/home/ben/workspace/timeCourse/data")
source("../code/R/tools.R")

counts<-load.data("plant","psa",0,120)
outputImages<-"/home/ben/Google Drive/PhD/writing/Timecourse/images/"
outputDir<-"/home/ben/workspace/timeCourse/data/deResults/"

subject<-"plant"  #should be plant or bacteria
track<-"control"  #should be psa or control, i.e. treatment
edgeRanalysis(counts,subject,track,0,120,outputImages,paste(outputDir,"edgeR/",sep=""))
limmaAnalysis(counts,subject,track,0,120,outputImages,paste(outputDir,"limma/",sep=""))

subject<-"plant"  #should be plant or bacteria
track<-"psa"  #should be psa or control, i.e. treatment
edgeRanalysis(counts,subject,track,0,120,outputImages,paste(outputDir,"edgeR/",sep=""))
limmaAnalysis(counts,subject,track,0,120,outputImages,paste(outputDir,"limma/",sep=""))

subject<-"bacteria"  #should be plant or bacteria
track<-"psa"  #should be psa or control, i.e. treatment
edgeRanalysis(counts,subject,track,0,120,outputImages,paste(outputDir,"edgeR/",sep=""))
limmaAnalysis(counts,subject,track,0,120,outputImages,paste(outputDir,"limma/",sep=""))

## counts are loaded