source("http://bioconductor.org/biocLite.R")
library("DESeq2")


setwd("/home/ben/workspace/timeCourse/data")

## load data into rawCounts
datafile<-file.choose()
counts <- read.delim(datafile,row.names=1,quote="")
colnames(counts)

# get's rid of that annoying last column with nothing in it that keeps appearing from somewhere. 
drops <- c("X")
counts<-counts[,!(names(counts) %in% drops)]

## get all the columns in the right order to make selecting groups easier.
columnnames <- colnames(counts)
columnnames <- gsub("X", "", columnnames)
colnames(counts) <- columnnames
colnames(counts)
counts <- counts[order(as.numeric(colnames(counts))) ]
head(counts)

#select the data wanted
rawCounts <- counts[,c(1,3,5,13:42)] # extract the 0-24 hour time points
rawCounts <- counts[,c(2:7,12,23,33:57)] #24hr cycle timepoints
rawCounts <- counts[,c(1,3,5:7,11,13:17,25:27,29,30,33:35,40:43,45,48:50,55,57)] #control time points
rawCounts <- counts[,c(3,6,8,10,12,16,17,18,25,26,27,28,29,30,34,37,38,41,42,43,44,49,50,51,55,56,57)] #psa time points
rawCounts <- counts #everything

head(rawCounts)
colnames(rawCounts) <- time
write.csv(rawCounts,file="psaCounts.csv")
row.names = colnames( rawCounts )
row.names
design = data.frame(
  row.names = colnames( rawCounts ),
  condition = c("Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Psa"),
  time = c("24","72","72","48","72","1.5","1.5","1.5","3","3","3","6","6","6","12","12","12","24","24","48","48","120","120","120","96","96","96"),
  libtype = c("paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end","paired-end")
)

colData <- data.frame(row.names=colnames( rawCounts ), time=as.factor(design$time))
dds <- DESeqDataSetFromMatrix(countData = rawCounts,colData = colData, design =~  time)
dds <- DESeq(dds)
res <- results(dds)
plotMA(res)