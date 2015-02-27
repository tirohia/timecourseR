library("TDARACNE")

datafile<-file.choose()
networkCounts = read.table( datafile,sep=",", header=TRUE,row.names=1 , stringsAsFactors = FALSE )
networkCounts<-cleanDataNames(networkCounts)
design <- getDataRange(rawCounts,24,120)
datacols<-getDataCols(design)
networkCounts<-networkCounts[,datacols]  

#phenotype data for annotated data frame. 
pDataFile <- file.path("pdata.csv")
pData <- read.csv(pDataFile,header=TRUE, sep=",")
phenoData <- new("AnnotatedDataFrame",data=pData, varMetadata=metadata)

#annotated data frame. 
colnames(networkCounts)<-rownames(pData)
all(rownames(pData) == colnames(networkCounts))
metadata <- data.frame(labelDescription=c("Condition", "Assay","Time"))

phenoData <- new("AnnotatedDataFrame", data = pData, varMetadata = metadata)

## Experiment data - complete bollocks atm
experimentData <- new("MIAME",name="Pierre Fermat", lab="Francis Galton Lab", contact="pfermat@lab.not.exist")

#Matrix
NM3M<-as.matrix(networkCounts)
NM3M<- NM3M+1
colnames(M3)<-rownames(pData)

##Expression set
exampleSet <- ExpressionSet(assayData=NM3M, phenoData=phenoData,experimentData=experimentData,annotation="hgu95av2")
aracneGraph<-TDARACNE(exampleSet,11,"netIRMAon",delta=3,likehood=1.2,norm=2,logarithm=1,thresh=0,ksd=0,0.15);

