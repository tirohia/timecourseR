library("edgeR")
library(RColorBrewer)

subject<-"bacteria"  #should be plant or psa
track<-"psa"  #should be psa or control, i.e. treatment

outputDir<-"/home/ben/Google Drive/PhD/writing/Timecourse/images/"

setwd("/home/ben/workspace/timeCourse/data")
source("../code/R/tools.R")
if (subject=="plant"){
  datafile="/home/ben/workspace/timeCourse/data/actinidia/cornellGenomeModelTH-TCFrequencyMatrix.csv"
}else if (subject=="bacteria"){
  datafile="/home/ben/workspace/timeCourse/data/psa/psaFrequencyMatrix.csv"
}else {
  print ("please specify subject <- plant or bacteria")
}


## load data into rawCounts
counts = read.table( datafile,sep=",", header=TRUE, row.names=1,stringsAsFactors = FALSE )
counts <-cleanDataNames(counts)
design <- getDataRange(counts,1.5,120,track)
datacolumns<-getDataCols(design)
counts<-counts[,datacolumns]

if (subject=="bacteria"){
  ## remove ribosomal entries from psa data. 
  ribosomalEntriesfile="/home/ben/workspace/models/psa/ribosomalIYO"
  ribosomalEntries = read.table( ribosomalEntriesfile,sep=",", header=FALSE,stringsAsFactors = FALSE )
  counts<-counts[ ! rownames(counts) %in% ribosomalEntries, ] 
  counts <- head(counts, -5)
  #tail(counts)
}


####################
### glm. 


#condition=as.factor(c(rep("control",each=15),rep("treated",each=15)))
condition=as.factor(rep("treated",each=27))

time=c(1,2,3,4,5,6,7,8,9)
time=as.factor(c(rep(time,each=3)))

targets <- data.frame(condition=condition,time=time)


#Group <- as.factor(paste(targets$condition,targets$time,sep="."))
#targets<-cbind(targets,Group=Group)

y <- DGEList(counts=counts,group=time)
design <- model.matrix(~time)
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
fit <- glmFit(y,design)

#targets$condition <- relevel(targets$condition, ref="treated")
#design <- model.matrix(~condition + time, data=targets)

total<-vector('character')  #for unique genes
geneMap<-matrix(data=0,nrow=9,ncol=9)
rownames(geneMap)<-paste("time",c(1:9),sep="")
colnames(geneMap)<-paste("time",c(1:9),sep="")


for (i in c(1:9) ){
  for (j in c(1:9)){
    if(i==j){
      geneMap[paste("time",i,sep=""),paste("time",j,sep="")]<-0
      next
    }
    a = c(0,0,0,0,0,0,0,0,0)
    if (i==1){
      a<-replace(a,j,1)
    }
    else if (j==1){
      a<-replace(a,i,-1)
    }
    else{
      a<-replace(a,c(i,j),c(1,-1))
    }
    lrt <- glmLRT(fit, contrast=a)
    tags<-as.data.frame(topTags(lrt,39039,sort.by="logFC"))
    tags<- tags[tags$FDR < 0.05,]
    tagsUp<- tags[tags$logFC > 2,]
    tagsDown<- tags[tags$logFC < -2,]
    total<-c(total,rownames(tagsUp))
    write.csv(rownames(tagsUp),paste("/home/ben/workspace/timeCourse/data/deResults/edgeR/edgeRpsa",i,"-",j,"up.csv",sep=""),row.names=FALSE)
    write.csv(rownames(tagsDown),paste("/home/ben/workspace/timeCourse/data/deResults/edgeR/edgeRpsa",i,"-",j,"down.csv",sep=""),row.names=FALSE)
    
    if (is.integer(nrow(tagsUp))){
      geneMap[paste("time",i,sep=""),paste("time",j,sep="")]<-nrow(tagsUp)
    }
    
  }  
}

###

#get number of unique genes. 

###


if(subject=="plant" && track=="control"){
  filenameUp = paste(outputDir,"edgeRGeneMapControlUp.png",sep="")
  filenameDown = paste(outputDir,"edgeRGeneMapControlDown.png",sep="")
  mainUp<-paste("Actinidia, control - up",sep="")
  mainDown<-paste("Actinidia, control - down",sep="")
  palette<-brewer.pal(9,name="Greens")
}else if (subject=="plant" && track=="psa"){
  filenameUp = paste(outputDir,"edgeRGeneMapInnoculatedUp.png",sep="")
  filenameDown = paste(outputDir,"edgeRGeneMapInnoculatedDown.png",sep="")
  mainUp<-paste("Actinidia, innoculated - up",sep="")
  mainDown<-paste("Actinidia, innoculated - down",sep="")
  palette<-brewer.pal(9,name ="Blues")
}else if (subject=="bacteria" && track=="psa"){
  filenameUp = paste(outputDir,"edgeRGeneMapPsaUp.png",sep="")
  filenameDown = paste(outputDir,"edgeRGeneMapPsaDown.png",sep="")
  mainUp<-paste("Psa, bacterial genes - up",sep="")
  mainDown<-paste("Psa, bacterial genes - down",sep="")
  palette<-brewer.pal(9,name="Reds")
}else{
  print("there should have been an error before here but either the plant/bacteria or psa/control thing hasn't been specified")
}

colnames(geneMap)<-c(1.5,3,6,12,24,48,72,96,120)
rownames(geneMap)<-c(1.5,3,6,12,24,48,72,96,120)

matUp<-geneMap
matUp[upper.tri(geneMap)]<-0
png(filename=filenameUp,width=800)
heatmap.2( log(matUp+1),scale="none",col = palette,margins = c(5, 15),trace = "none",lhei = c(0.75,3),xlab = "Comparison", na.color="white", main = mainUp ,Colv = FALSE,Rowv =FALSE, dendrogram = "none")
dev.off()



#heatmap.2( log(matUp+1),scale="none",col = palette,margins = c(5, 15),trace = "none",lhei = c(0.75,3),xlab = "Comparison", na.color="white", main = mainUp ,Colv = FALSE,Rowv =FALSE, dendrogram = "none")

matDown<-geneMap
matDown[lower.tri(geneMap)]<-0
png(filename=filenameDown,width=800)
heatmap.2( log(matDown+1),scale="none",col = palette,margins = c(5, 15),trace = "none",lhei = c(0.75,3),xlab = "Comparison", na.color="white", main = mainDown,Colv = FALSE,Rowv =FALSE, dendrogram = "none")
dev.off()


print("total number of gene changes:")
print(length(total))
print("total number of unique genes changing:")
print(length(unique(total)))



###################################

#design <- model.matrix(~condition * time, data=targets)
#fit <- glmFit(y, design)
#colnames(fit)

#fit <- glmFit(y, design)
#lrt <- glmLRT(fit, contrast=c(0,-1,0,1))
#tags<-as.data.frame(topTags(lrt,39039,sort.by="logFC"))
#tags<- tags[tags$FDR < 0.05,]
#tagsUp<- tags[tags$logFC > 2,]
#tagsDown<- tags[tags$logFC < -2,]
#dim(tagsUp)
#dim(tagsDown)




