require("edgeR")
require(RColorBrewer)
require(gplots)

## data should be loaded elsewhere using load.data() from tools.R

####################
### glm. 

edgeRanalysis<-function(counts,subject,track,time1,time2,outputImages,outputDir){
  #condition=as.factor(c(rep("control",each=15),rep("treated",each=15)))
  if (track=="all"){
    design <- getDataRange(time1,time2)
  }else{
    design <- getDataRange(time1,time2,track)
  }
  design <- getDataRange(0,120,"psa")
  timeList<-design$time
  numberOfPoints<-length(design$time)
  condition=as.factor(rep("treated",each=numberOfPoints*3))
  time=c(1:length(design$time))
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
  geneMap<-matrix(data=0,nrow=numberOfPoints,ncol=numberOfPoints)
  rownames(geneMap)<-paste("time",c(1:numberOfPoints),sep="")
  colnames(geneMap)<-paste("time",c(1:numberOfPoints),sep="")


  for (i in c(1:numberOfPoints) ){
    for (j in c(1:numberOfPoints)){
      if(i==j){
        geneMap[paste("time",i,sep=""),paste("time",j,sep="")]<-0
        next
      }
      a = c(rep(0,each=numberOfPoints))
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
      if(!is.null(outputDir)){
        if(!dir.exists(file.path(outputDir))){
          dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
        }
        write.csv(rownames(topUp),paste(outputDir,"/edgeR",subject,"-",track,i,"-",j,"up.csv",sep=""),row.names=FALSE)
        write.csv(rownames(topDown),paste(outputDir,"/edgeR",subject,"-",track,i,"-",j,"down.csv",sep=""),row.names=FALSE)
      }
      if (is.integer(nrow(tagsUp))){
        geneMap[paste("time",i,sep=""),paste("time",j,sep="")]<-nrow(tagsUp)
      }
    }  
  }

###

#get number of unique genes. 
  if(!is.null(outputImages)){
    if(subject=="plant" && track=="control"){
      filenameUp = paste(outputImages,"edgeRGeneMapControlUp",time1,"-",time2,".png",sep="")
      filenameDown = paste(outputImages,"edgeRGeneMapControlDown",time1,"-",time2,".png",sep="")
      mainUp<-paste("Actinidia, control - up",sep="")
      mainDown<-paste("Actinidia, control - down",sep="")
      palette<-brewer.pal(9,name="Greens")
    }else if (subject=="plant" && track=="psa"){
      filenameUp = paste(outputImages,"edgeRGeneMapInnoculatedUp",time1,"-",time2,".png",sep="")
      filenameDown = paste(outputImages,"edgeRGeneMapInnoculatedDown",time1,"-",time2,".png",sep="")
      mainUp<-paste("Actinidia, innoculated - up",sep="")
      mainDown<-paste("Actinidia, innoculated - down",sep="")
      palette<-brewer.pal(9,name ="Blues")
    }else if (subject=="bacteria" && track=="psa"){
      filenameUp = paste(outputImages,"edgeRGeneMapPsaUp",time1,"-",time2,".png",sep="")
      filenameDown = paste(outputImages,"edgeRGeneMapPsaDown",time1,"-",time2,".png",sep="")
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

}
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




