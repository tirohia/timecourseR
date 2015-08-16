baseDir<-"/home/ben/workspace/timeCourse/data/deResults/"

##DESeq Genes
directory<-paste(baseDir,"deseq",sep="")
filenames <- list.files(directory, pattern="*.csv", full.names=TRUE)
deseqGenes <- vector()
for (file in filenames){
  genes <-read.csv(file)
  deseqGenes <-c(deseqGenes, as.vector(genes[,1]))
}
deseqGenes<-unique(deseqGenes)
length(deseqGenes)

##EdgeR genes
directory<-paste(baseDir,"edgeR",sep="")
filenames <- list.files(directory, pattern="*.csv", full.names=TRUE)
edgeRGenes <- vector()
for (file in filenames){
  genes <-read.csv(file)
  edgeRGenes <-c(edgeRGenes, as.vector(genes[,1]))
}

edgeRGenes<-unique(edgeRGenes)
length(edgeRGenes)

##Limma Genes
directory<-paste(baseDir,"limma",sep="")
filenames <- list.files(directory, pattern="*.csv", full.names=TRUE)
limmaGenes <- vector()
for (file in filenames){
  genes <-read.csv(file)
  limmaGenes <-c(limmaGenes, as.vector(genes[,1]))
}

limmaGenes<-unique(limmaGenes)
length(limmaGenes)

##Betr Genes
#betr uses rlog transformed data. 

betrpsaGenes<-read.csv("/home/ben/workspace/timeCourse/data/deResults/betr/psaGenes.csv")
betrplantGenes<-read.csv("/home/ben/workspace/timeCourse/data/deResults/betr/plantGenes.csv")
betrGenes<-rbind(betrpsaGenes,betrplantGenes)
colnames(betrGenes)=c("gene","probability")
betrGenes<-betrGenes$gene


##timecourse Genes
## timecourse I've used voom to prep the data becuase timecourse is designed for microarray data. 

timecourseGenes<-read.csv("/home/ben/workspace/timeCourse/data/deResults/timecourse/rankedTimecourseGenesVoomComplatecases-farBest.csv",header=FALSE)
timecourseBacGenes<-read.csv("/home/ben/workspace/timeCourse/data/deResults/timecourse/bacteria-all-rankedTimecourseGenes.csv",header=FALSE)
timecourseGenes<-rbind(timecourseGenes,timecourseBacGenes)
timecourseGenes<-timecourseGenes[,1]


## previous list.
previous<-read.csv("/home/ben/workspace/timeCourse/results/tripleDE.csv")
previous<-previous[,1]
length(previous)
##intersection


#combos<-combinations(5,3,c("edgeRGenes","timecourseGenes","deseqGenes","betrGenes","limmaGenes"))

allIntersections<-intersect(intersect(intersect(intersect(edgeRGenes,deseqGenes),betrGenes),timecourseGenes),limmaGenes)
length(allIntersections)

i1<-intersect(intersect(betrGenes,  deseqGenes), edgeRGenes)     
i2<-intersect(intersect(betrGenes,  deseqGenes), limmaGenes)
i3<-intersect(intersect(betrGenes,  deseqGenes), timecourseGenes)
i4<-intersect(intersect(betrGenes,  edgeRGenes), limmaGenes)
i5<-intersect(intersect(betrGenes,  edgeRGenes), timecourseGenes)
i6<-intersect(intersect(betrGenes,  limmaGenes), timecourseGenes)
i7<-intersect(intersect(deseqGenes, edgeRGenes), limmaGenes)
i8<-intersect(intersect(deseqGenes, edgeRGenes), timecourseGenes)
i9<-intersect(intersect(deseqGenes, limmaGenes), timecourseGenes)
i10<-intersect(intersect(edgeRGenes, limmaGenes), timecourseGenes)


all<-c(i1,i2,i3,i4,i5,i6,i8,i9,i10)
all<-unique(all)
length(all)

expandedAll<-c(i1,i2,i3,i4,i5,i6,i8,i9,i10,i7)
expandedAll<-unique(expandedAll)
length(expandedAll)

table(previous %in% unique(all))
table(previous %in% unique(expandedAll))

## this bit is looking at some of the notable ones I've pulled out in previous ones - just to make sure I'm doing the right thing.

c("Achn229111","Achn309921","Achn324811","Achn199471","Achn334501","Achn334481","Achn193721","Achn033321","Achn261701") %in% betrGenes

transcriptionFactors<-c("Achn241691","Achn014741","Achn025281","Achn107721","Achn110911","Achn120901","Achn241691","Achn126231","Achn241691","Achn267561","Achn290511","Achn322901","Achn380061","Achn383541")
transcriptionFactors %in% all
transcriptionFactors %in% expandedAll

responseProteins<-c("Achn241691","Achn014741","Achn025281","Achn107721","Achn110911","Achn120901","Achn241691","Achn126231","Achn241691","Achn267561","Achn290511","Achn322901","Achn380061","Achn383541")
responseProteins %in% all
responseProteins %in% expandedAll

thing %in% all
thing<-c("Achn033321","Achn033321","Achn334541","Achn229111","Achn353281","Achn353831","Achn196241","Achn309921","Achn324811","Achn334501","Achn126231","Achn160161","Achn334481","Achn253571","Achn199471","Achn277471","Achn059341","Achn287861","Achn029231","Achn367821","Achn169421","Achn193721","Achn374841")
thing2<-c("Achn240741","Achn103081","Achn332771","Achn120901","Achn332071","Achn245451","Achn273491","Achn322901","Achn327311")
