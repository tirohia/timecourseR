source("http://bioconductor.org/biocLite.R")
library("BHC")
library("pheatmap")
library("dendextend")
library('gplots')

#datafile<-"/home/ben/workspace/tiemcourse/data/networkCandidates/networkCandidates.csv"

## get Medians from file. The file is generated in intersection.R It should consist of 1 column for each timepoint, 
## with log transformed values. 

datafile<-"/home/ben/workspace/timeCourse/data/deResults/networkGenesMedians.csv"
counts <- read.csv(datafile,row.names=1)
#counts<-medians

itemLabels <- rownames(counts)
counts<-as.matrix(counts)
timePoints<-c( 24.0,48.0 ,72.0, 96.0,120.0)

hc2 <- bhc(counts, itemLabels, 0, timePoints, "time-course", numReps=1, noiseMode=0,  verbose=TRUE)

#make hc via bhc. Cut into six using cuttree and then six of those lists of names. 

## output the clusters that bhc finds. 
WriteOutClusterLabels(hc2, "labels.txt", verbose=TRUE)

## output the clusters when the number of clusters is pre-defined (k)
nodelist<-cutree(hc2,k = 6)
write.csv(nodelist[nodelist==1],"nodelist1.csv")
write.csv(nodelist[nodelist==2],"nodelist2.csv")
write.csv(nodelist[nodelist==3],"nodelist3.csv")
write.csv(nodelist[nodelist==4],"nodelist4.csv")
write.csv(nodelist[nodelist==5],"nodelist5.csv")
write.csv(nodelist[nodelist==6],"nodelist6.csv")

## output a dendrogram
png("hc2.png", width=800, height=800)
  plot(hc2)
dev.off()

## output a heatmap with the dendrogram to the side. 
png("heatmap.png", width=800, height=800)
  heatmap.2(counts,trace="none",col = bluered(10),  Rowv = hc2)
dev.off()

## output the details of the lists from earlier. Yes, I know, this could be done better. 
group1<-read.csv("clustering/bhc/nodelist1.csv", row.names=1)
group2<-read.csv("clustering/bhc/nodelist2.csv", row.names=1)
group3<-read.csv("clustering/bhc/nodelist3.csv", row.names=1)
group4<-read.csv("clustering/bhc/nodelist4.csv", row.names=1)
group5<-read.csv("clustering/bhc/nodelist5.csv", row.names=1)
group6<-read.csv("clustering/bhc/nodelist6.csv", row.names=1)



## should be done with a loop/apply. Can't be arsed at the moment. Bloody R and the lack of dynamic variable assignment.   
group1<-rownames(group1)
bac.details<-geneDetails("bacteria", list = group1[grep("IYO",group1)]) 
plant.details<-geneDetails("plant", list = group1[grep("Achn",group1)])
details1<-rbind(plant.details,bac.details)
write.csv(details1,"clustering/bhc/bhc1.csv")

group2<-rownames(group2)
bac.details<-geneDetails("bacteria", list = group2[grep("IYO",group2)]) 
plant.details<-geneDetails("plant", list = group2[grep("Achn",group2)])
details<-rbind(plant.details,bac.details)
write.csv(details,"clustering/bhc/bhc2.csv")

group3<-rownames(group3)
bac.details<-geneDetails("bacteria", list = group3[grep("IYO",group3)]) 
plant.details<-geneDetails("plant", list = group3[grep("Achn",group3)])
details<-rbind(plant.details,bac.details)
write.csv(details,"clustering/bhc/bhc3.csv")

group4<-rownames(group4)
bac.details<-geneDetails("bacteria", list = group4[grep("IYO",group4)]) 
plant.details<-geneDetails("plant", list = group4[grep("Achn",group4)])
details<-rbind(plant.details,bac.details)
write.csv(details,"clustering/bhc/bhc4.csv")

group5<-rownames(group5)
bac.details<-geneDetails("bacteria", list = group5[grep("IYO",group5)]) 
plant.details<-geneDetails("plant", list = group5[grep("Achn",group5)])
details<-rbind(plant.details,bac.details)
write.csv(details,"clustering/bhc/bhc5.csv")

group6<-rownames(group6)
bac.details<-geneDetails("bacteria", list = group6[grep("IYO",group6)]) 
plant.details<-geneDetails("plant", list = group6[grep("Achn",group6)])
details<-rbind(plant.details,bac.details)
write.csv(details,"clustering/bhc/bhc6.csv")

intersect(group6,aGroup)
intersect(group6,bGroup)
intersect(group6,cGroup)
intersect(group6,dGroup)
intersect(group6,eGroup)
intersect(group6,fGroup)
group5
