#source("http://bioconductor.org/biocLite.R")
library("limma")
library("Biobase")
library("gplots")
library(RColorBrewer)

subject<-"bacteria"  #should be plant or bacteria
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


## load data into counts
counts = read.table( datafile,sep=",", header=TRUE, row.names=1,stringsAsFactors = FALSE )
counts <-cleanDataNames(counts)
design <- getDataRange(counts,1.5,120,track)
datacolumns<-getDataCols(design)
counts<-counts[,datacolumns]

if (subject=="psa"){
## remove ribosomal entries from psa data. 
  ribosomalEntriesfile="/home/ben/workspace/models/psa/ribosomalIYO"
  ribosomalEntries = read.table( ribosomalEntriesfile,sep=",", header=FALSE,stringsAsFactors = FALSE )
  counts<-counts[ ! rownames(counts) %in% ribosomalEntries, ] 
  counts <- head(counts, -5)
  #tail(counts)
}

####################
### glm. 


time=c("s1","s2","s3","s4","s5","s6","s7","s8","s9")
time=as.factor(c(rep(time,each=3)))
design <- model.matrix(~0+time)
colnames(design) <- c("s1","s2","s3","s4","s5","s6","s7","s8","s9")


dge <- DGEList(counts=counts)
dge <- calcNormFactors(dge)
v <- voom(dge,design,plot=TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)


total<-vector('character')   #for unique genes
geneMap<-matrix(data=0,nrow=9,ncol=9)
rownames(geneMap)<-paste("time",c(1:9),sep="")
colnames(geneMap)<-paste("time",c(1:9),sep="")


for (i in c(1:9) ){
  for (j in c(1:9)){
    if(i==j){
      geneMap[paste("time",i,sep=""),paste("time",j,sep="")]<-0
      next
    }
    cont.matrix <- makeContrasts(
      assign("contrast",paste("s",i,"-s",j,sep="")),
      levels=design
    )
    a = c(0,0,0,0,0,0,0,0,0)
    a<-replace(a,c(i,j),c(1,-1))
    cont.matrix[,1]<-a
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2)
    
    top<-as.data.frame(topTable(fit2,  number=Inf, sort.by="P",adjust="BH"))
    top<- top[top$adj.P.Val < 0.05,]
    topUp<- top[top$logFC > 2,]
    topDown<- top[top$logFC < -2,]
    
    write.csv(rownames(topUp),paste("/home/ben/workspace/timeCourse/data/deResults/limma/limma",subject,i,"-",j,"up.csv",sep=""),row.names=FALSE)
    write.csv(rownames(topDown),paste("/home/ben/workspace/timeCourse/data/deResults/limma/limma",subject,i,"-",j,"down.csv",sep=""),row.names=FALSE)
    # only counting the Up genes as given the reflection, as it runs throug the second half of the (i,j) matrix, it'll switch and the up's will be down
    # this is preventing changes being counted twice.
    total<-c(total,rownames(topUp))
    if (is.integer(nrow(topUp))){
        geneMap[paste("time",i,sep=""),paste("time",j,sep="")]<-nrow(topUp)
    }
    
  }
}



####


if(subject=="plant" && track=="control"){
  filenameUp = paste(outputDir,"limmaGeneMapControlUp.png",sep="")
  filenameDown = paste(outputDir,"limmaGeneMapControlDown.png",sep="")
  mainUp<-paste("Actinidia, control - up",sep="")
  mainDown<-paste("Actinidia, control - down",sep="")
  palette<-brewer.pal(9,name="Greens")
}else if (subject=="plant" && track=="psa"){
  filenameUp = paste(outputDir,"limmaGeneMapInnoculatedUp.png",sep="")
  filenameDown = paste(outputDir,"limmaGeneMapInnoculatedDown.png",sep="")
  mainUp<-paste("Actinidia, innoculated - up",sep="")
  mainDown<-paste("Actinidia, innoculated - down",sep="")
  palette<-brewer.pal(9,name ="Blues")
}else if (subject=="bacteria" && track=="psa"){
  filenameUp = paste(outputDir,"limmaGeneMapPsaUp.png",sep="")
  filenameDown = paste(outputDir,"limmaGeneMapPsaDown.png",sep="")
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

par(mfrow=c(1,1))
heatmap.2( log(matDown+1),scale="none",col = palette,margins = c(5, 15),trace = "none",lhei = c(0.75,3),xlab = "Comparison", na.color="white", main = mainDown,Colv = FALSE,Rowv =FALSE, dendrogram = "none")













