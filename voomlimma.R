#source("http://bioconductor.org/biocLite.R")
library("limma")
library("Biobase")
library("gplots")
library(gtools)
library(RColorBrewer)


limmaAnalysis<-function(counts,subject,track,time1,time2,outputImages,outputDir){
  #condition=as.factor(c(rep("control",each=15),rep("treated",each=15)))
  if (track=="all"){
    design <- getDataRange(time1,time2)
  }else{
    design <- getDataRange(time1,time2,track)
  }
  design <- getDataRange(time1,time2,"psa")
  design
  timeList<-design$time
  numberOfPoints<-length(design$time)
  condition=as.factor(rep("treated",each=numberOfPoints*3))
  time=c(1:length(design$time))
  time=as.factor(c(rep(time,each=3)))
  targets <- data.frame(condition=condition,time=time)

####################
### glm. 

  #time=as.factor(c(rep(time,each=3)))
  design <- model.matrix(~0+time)
  dge <- DGEList(counts=counts)
  dge <- calcNormFactors(dge)
  
  v <- voom(dge,design,plot=TRUE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)


  total<-vector('character')   #for unique genes
  geneMap<-matrix(data=0,nrow=numberOfPoints,ncol=numberOfPoints)
  rownames(geneMap)<-paste("time",c(1:numberOfPoints),sep="")
  colnames(geneMap)<-paste("time",c(1:numberOfPoints),sep="")


  for (i in c(1:numberOfPoints) ){
    for (j in c(1:numberOfPoints)){
      if(i==j){
        geneMap[paste("time",i,sep=""),paste("time",j,sep="")]<-0
        next
      }
      
      cont.matrix <- makeContrasts(time1-time2, levels=design)
      a = c(rep(0,each=numberOfPoints))
      a<-replace(a,c(i,j),c(1,-1))
      cont.matrix[,1]<-a
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2)
      
      top<-as.data.frame(topTable(fit2,  number=Inf, sort.by="P",adjust="BH"))
      top<- top[top$adj.P.Val < 0.05,]
      topUp<- top[top$logFC > 2,]
      topDown<- top[top$logFC < -2,]
      if(!is.null(outputDir)){
          if(!dir.exists(file.path(outputDir))){
            dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
          }
          write.csv(rownames(topUp),paste(outputDir,"/limma",subject,"-",track,i,"-",j,"up.csv",sep=""),row.names=FALSE)
          write.csv(rownames(topDown),paste(outputDir,"/limma",subject,"-",track,i,"-",j,"down.csv",sep=""),row.names=FALSE)
      }
      # only counting the Up genes as given the reflection, as it runs throug the second half of the (i,j) matrix, it'll switch and the up's will be down
      # this is preventing changes being counted twice.
      total<-c(total,rownames(topUp))
      if (is.integer(nrow(topUp))){
          geneMap[paste("time",i,sep=""),paste("time",j,sep="")]<-nrow(topUp)
      }
    }
  }



####
  if(!is.null(outputImages)){
    if(subject=="plant" && track=="control"){
      filenameUp = paste(outputImages,"limmaGeneMapControlUp",time1,"-",time2,".png",sep="")
      filenameDown = paste(outputImages,"limmaGeneMapControlDown",time1,"-",time2,".png",sep="")
      mainUp<-paste("Actinidia, control - up",sep="")
      mainDown<-paste("Actinidia, control - down",sep="")
      palette<-brewer.pal(9,name="Greens")
    }else if (subject=="plant" && track=="psa"){
      filenameUp = paste(outputImages,"limmaGeneMapInnoculatedUp",time1,"-",time2,".png",sep="")
      filenameDown = paste(outputImages,"limmaGeneMapInnoculatedDown",time1,"-",time2,".png",sep="")
      mainUp<-paste("Actinidia, innoculated - up",sep="")
      mainDown<-paste("Actinidia, innoculated - down",sep="")
      palette<-brewer.pal(9,name ="Blues")
    }else if (subject=="bacteria" && track=="psa"){
      filenameUp = paste(outputImages,"limmaGeneMapPsaUp",time1,"-",time2,".png",sep="")
      filenameDown = paste(outputImages,"limmaGeneMapPsaDown",time1,"-",time2,".png",sep="")
      mainUp<-paste("Psa, bacterial genes - up",sep="")
      mainDown<-paste("Psa, bacterial genes - down",sep="")
      palette<-brewer.pal(9,name="Reds")
    }else{
      print("there should have been an error before here but either the plant/bacteria or psa/control thing hasn't been specified")
    }
  
    colnames(geneMap)<-timeList
    rownames(geneMap)<-timeList

    matUp<-geneMap
    matUp[upper.tri(geneMap)]<-0
    png(filename=filenameUp,width=800)
    heatmap.2( log(matUp+1),scale="none",col = palette,margins = c(5, 15),trace = "none",lhei = c(0.75,3),xlab = "Hours post innoculation", na.color="white", main = mainUp ,Colv = FALSE,Rowv =FALSE, dendrogram = "none")
    dev.off()

    matDown<-geneMap
    matDown[lower.tri(geneMap)]<-0
    png(filename=filenameDown,width=800)
    heatmap.2( log(matDown+1),scale="none",col = palette,margins = c(5, 15),trace = "none",lhei = c(0.75,3),xlab = "Hours post innoculation", na.color="white", main = mainDown,Colv = FALSE,Rowv =FALSE, dendrogram = "none")
    dev.off()
  }
    print("total number of gene changes:")
    print(length(total))
    print("total number of unique genes changing:")
    print(length(unique(total)))

    par(mfrow=c(1,1))
  heatmap.2( log(matDown+1),scale="none",col = palette,margins = c(5, 15),trace = "none",lhei = c(0.75,3),xlab = "Hours post innoculation", na.color="white", main = mainDown,Colv = FALSE,Rowv =FALSE, dendrogram = "none")

}











