library(stringr)
library("Biobase")
setwd("/home/ben/workspace/timeCourse/data")


timeTable = read.table( "lookupsheet.csv",sep=",", header=TRUE, stringsAsFactors = FALSE )

load.data<-function(subject,track,time1,time2){
  if (subject=="plant"){
    datafile="/home/ben/workspace/timeCourse/data/actinidia/cornellGenomeModelTH-TCFrequencyMatrix.csv"
    
  }else if (subject=="bacteria"){
    datafile="/home/ben/workspace/timeCourse/data/psa/psaFrequencyMatrix.csv"
  }else {
    stop ("please specify subject <- plant or bacteria" )
  }
  
  ## load data into counts
  counts = read.table( datafile,sep=",", header=TRUE, row.names=1,stringsAsFactors = FALSE )
  counts <-cleanDataNames(counts)
  if (track=="all"){
    design <- getDataRange(time1,time2)
  }else{
    design <- getDataRange(time1,time2,track)
  }
  datacolumns<-getDataCols(design)
  counts<-counts[,datacolumns]
  
  if (subject=="bacteria"){
    ## remove ribosomal entries from psa data. 
    ribosomalEntriesfile="/home/ben/workspace/models/psa/ribosomalIYO"
    ribosomalEntries = read.table( ribosomalEntriesfile,sep=",", header=FALSE,stringsAsFactors = FALSE )
    counts<-counts[ ! rownames(counts) %in% ribosomalEntries, ] 
    counts <- head(counts, -5)
    
    ribo<-read.csv("/home/ben/workspace/timeCourse/data/psa/NZv-13invitroRNA-seq.csv")
    riboNames<-counts[grep("ribosom", ribo[,8]), ]
    riboNames<- rownames(riboNames)
    counts<-counts[ ! rownames(counts) %in% riboNames, ] 
    
    #tail(counts)
  }
  return(counts)
}

cleanDataNames <- function(counts){
  #sometimes the csv has an empty column on for some reason. Get rid of it. 
  drops <- c("X")
  counts<-counts[,!(names(counts) %in% drops)]
  #get rid of the stupid X's in the column names
  columnnames <- colnames(counts)
  columnnames <- gsub("X", "", columnnames)
  colnames(counts) <- columnnames
  colnames(counts)
  counts <- counts[order(as.numeric(colnames(counts))) ]
}



getDataRange<-function(time1,time2,condition="all"){
  intable <- timeTable[timeTable$time >=time1,]
  fintable<- intable[intable$time<= time2,]
  #print(fintable[fintable$condtion==condition,])
  
  if (condition=="all" || is.null(condition)){
    return(fintable)
  }
  else if (condition=="psa"){
    return(fintable[fintable$condition=="psa",])
  }
  else if (condition=="control"){
    
    return(fintable[fintable$condition=="control",])
  }
}

getDataCols <- function(design){
  timecolumns<-NULL
  for(j in unique(design$time)){
    timecolumns<-append(timecolumns,design[design$time==j & design$condition=="control",3:5])
    
  }
  for(j in unique(design$time)){
    timecolumns<-append(timecolumns,design[design$time==j & design$condition=="psa",3:5])
  }
  names(timecolumns) <- NULL
  unlist(timecolumns)
}


getRawData <- function (rawdata, time1, time2,condition=NULL){
  #should be error checking in here to check that counts is the full 1-57 data frame. 
  #set1<-c(timeTable[timeTable$time == time1 && timeTable$condition == condition, ]$X1,timeTable[timeTable$time == time1 && timeTable$condition == condition, ]$X2,timeTable[timeTable$time == time1 && timeTable$condition == condition, ]$X3)
  #set2<-c(timeTable[timeTable$time == time2, ]$X1,timeTable[timeTable$time == time2, ]$X2,timeTable[timeTable$time == time2, ]$X3)
  print(head(rawdata))
  timeTable1<-timeTable[timeTable$condition == condition,]
  tableCols<-timeTable1[timeTable1$time == time1, ]
  set1<-unlist(tableCols[,3:5])
  tableCols<-timeTable1[timeTable1$time == time2, ]
  set2<-unlist(tableCols[,3:5])
  rawdata <- rawdata[,c(set1,set2)] 
}

gm_mean = function(x, thing, n,na.rm=TRUE){
  list<-lapply(split(x, ceiling(seq_along(x)/n)), function(x) exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))) 
  unlist(list)
  
  #exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  #  c("thing","thingngn","salnfhd")
}

getMedians <- function (data,time1,time2,condition=NULL){
  ## returns geometric means
  n<-3
  t(apply(data, 1, function(x) gm_mean(x,seq(1, ncol(m), n), n)))
  
  #assumes log transformed data is being passed in. Means I have to take the geometirc rather than the arithmetic mean. 
  
  #timeTable = read.table( "lookupsheet.csv",sep=",", header=TRUE, stringsAsFactors = FALSE )
  
  #empty data frame, setting only rownames
  #mediansDF<-data.frame(gene=rownames(rawData))
  #
  #unique(design$time)
  #thing<-apply(rawData, 1, function(x) tapply(x, rep(seq(1, length(x), n), each=n), mean))
  
  #gm_mean(rawData[1,1:3])
  
  #if (!is.null(condition)){
  #  design<-timeTable[timeTable$time >= time1 & timeTable$time <= time2 & timeTable$condition ==condition,]
  #}
  #else{
  #  design<-timeTable[timeTable$time >= time1 & timeTable$time <= time2,]
  #}
  #for(j in unique(design$time)){
  #  columns <-  design[design$time==j,3:5]
  #  print(columns)
  #  data<-rawData[,unlist(columns)]
  #  M3<-as.matrix(data)
  #  medians<-rowMedians(M3)
  #  #mediansDF[,paste(condition,j,sep="")]<-as.vector(medians)
  #  mediansDF[,paste(j,sep="")]<-as.vector(medians)
  #}
  #mediansDF$gene<-NULL
  #rownames(mediansDF)<-rownames(rawData)
  #return(mediansDF)
}


## Design for DE Seq. Only copes with two time points atm. 
# not sure if it's worth the time to make it cope with different lengths - i.e. all condition 1 vs all condition2
getDesign <- function (counts,time1,time2,factor=NULL){
  time1 <- as.numeric(str_extract_all(time1,"\\(?[0-9,.]+\\)?")[[1]])
  time2 <- as.numeric(str_extract_all(time2,"\\(?[0-9,.]+\\)?")[[1]])
  design = data.frame(
    row.names = colnames(counts),
    time= c(time1,time1,time1,time2,time2,time2),
    libType = c("paired-end","paired-end","paired-end","paired-end","paired-end","paired-end")
  )
  design
}




## Do the D.E expression and return the number of rows. 
# ideally this would also have a parameter write to file, 
# that won't work with the wrapper though, which I'm using to call this. Not ideal.
## Banded is the prefered method according to Anders, Huber et al. 

deExp <- function(counts,design,direction=NULL, banded=FALSE,outputName=FALSE){
  if(length(unique(design[,"time"]))<2){
    return(0)
  }
  else{
    colData <- data.frame(row.names=colnames( counts ), time=as.factor(design$time))
    
    dds <- DESeqDataSetFromMatrix(countData = counts,colData = colData, design =~  time)
    dds <- DESeq(dds)
    res <- results(dds)
    res <- na.omit(res)
    if (banded==FALSE){
      if (direction=="up"){res = res[ res$log2FoldChange > 2, ]}
      else if(direction=="down"){res = res[ res$log2FoldChange < -2, ]}
      res=res[res$padj < 0.05,]
    #write.csv(res,paste("deResults/psa",direction,design[1,"time"],"-",design[4,"time"],".csv",sep="")) #actInoculatedUp
      if(nrow(res)==0){
        return(0)
      }
    }
    else {
      res <- results(dds, lfcThreshold=2, altHypothesis="greaterAbs")
      #contrast=c("condition","Control","Psa")
      res <- na.omit(res)
      
      if (direction=="up"){
        res=res[res$log2FoldChange > 2 ,]
      }
      else if (direction=="down"){
        res=res[res$log2FoldChange < -2,]  
      }
      res<-res[res$padj< 0.05,]
      if(nrow(res)==0){
        return(0)
      }
    }
    res
    #log(nrow(res))
  }
}


createESet<-function(counts,design){
  ## creates an expression set object. 
  ## errors regarding rownames etc are likely due to counts not being derived from the correct design object. 
  ## inputs? Counts is obvious. Log transformed, hetero wotsit stabilsed, geometric medians.  via assay(rld). 
  ## design? ah. counts shouldn't be medians. Should be replicates. 
  exprs <- as.matrix(counts)    
  metadata<-data.frame(lableDescription=c("Innoculation vs Control","timepoint"),row.names=(c("condition","time")))
  print(design)
  design<-data.frame(condition=rep(design$condition,each=3), time=rep(design$time,each=3))
  print(design)
  rownames(design)<-as.character(colnames(counts))
  phenoData <- new("AnnotatedDataFrame", data=design, varMetadata=metadata)
  experimentData <- new("MIAME", name="Ben Curran", lab="PFR",contact="ben.curran@auckland.ac.nz", title="Kiwifruit-psa host-pathogen timecourse", abstract="ExpressionSet for dimensionality reduction", url="www.sbs.acukland.ac.nz",other=list(notes="Created from text files"))
  eSet <- ExpressionSet(assayData=exprs, phenoData=phenoData, experimentData=experimentData) 
}


geneDetails<-function(subject,list,keyword=NULL,printToFile=FALSE){
  ## reads the annotation file and extracts genes in the list. 
  ## if keyword is specified, finds all those genes with that keyword (i.e. "transcription factor") in the annotation. 
  if(subject=="plant"){
    deets<-read.csv("/home/ben/workspace/models/cornell/Kiwifruit_gene_annotation_AHRD.csv")
  }
  else if (subject=="bacteria"){
    deets<-read.csv("/home/ben/workspace/models/psa/Nz13v-duplicatesRemoved.csv")
    deets<-deets[,2:4]
   
  }
  listDetails<-deets[deets[,1] %in% list,]
  if(!is.null(keyword)){
    keywordLocations<-grep(keyword, listDetails[,2], ignore.case = TRUE, perl = TRUE)
    listDetails<-listDetails[keywordLocations,1:2]
    if(printToFile==TRUE){
      write.csv(listDetails,paste("deResults/geneDetails-",length(all),"-",keyword,".csv"))
    }
  }
  else{
    print(listDetails)
    if(printToFile==TRUE){
      write.csv(listDetails,paste("deResults/geneDetails-",length(all),".csv"),row.names=TRUE)
    }
  }
  colnames(listDetails)<-c("gene","annotation")
  listDetails[,1:2]
  
}

removeIsolates<-function(graph){
  isolates <- which(igraph::degree(graph) == 0)
  graph<-delete.vertices(graph, isolates)
  graph
}
