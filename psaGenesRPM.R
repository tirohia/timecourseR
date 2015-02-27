setwd("/home/ben/workspace/timeCourse/data")

setwd("/home/ben/workspace/timeCourse/data")
source("../code/R/tools.R")
datafile="cornellGenomeModelTH-TCFrequencyMatrix.csv"

readsMapped <- read.csv(datafile, stringsAsFactors=FALSE, row.names=1)
columnnames <- colnames(readsMapped)
columnnames <- gsub("X", "P", columnnames)
colnames(readsMapped) <- columnnames
colnames(readsMapped)

design <- getDataRange(rawCounts,24,120)
datacolumns<-getDataCols(design)
readsMapped<-readsMapped[,datacolumns]

readsMapped<-scale(readsMapped, center=FALSE, scale=sizeFactors)
head(readsMapped)  
scales <- c(37.347786,38.243575,42.958841,37.693273,38.833448,38.767085,39.531887,66.581748,58.79224,36.833425,36.766411,38.43041,41.391057,38.195016,42.786511,43.100056,39.045232,41.279979,38.454703,41.514493,37.64413,37.562333,37.627361,42.407026,42.816286,36.399431,42.436111,41.946002,39.902334,42.970703,39.638666,40.412914,39.782447,41.603102,39.267669,38.407583,44.687265,38.630291,38.819944,43.256377,42.243905,39.799117,38.818551,42.459334,44.892463,29.408106,30.442288,38.035275,33.161491,37.71993,31.988759,35.262211,32.950623,41.12849,36.21799,34.777404,31.711137)

x_name <- "time"
y_name <- "gene"


plot(0, 0, xlim=c(0, 120), ylim=c(0, 500), type="n")
for (i in 1:length(readsMapped[,1])) {
  gene <- as.numeric(as.vector(readsMapped[i,]))
  gene <- t(t(gene)/scales)
  d<- data.frame( time = time,gene =gene )
  d <- d[order(time),]
  loess_fit <- loess( gene~ time , d)
  xl <- seq(min(time),max(time), (max(time) - min(time))/1000)
  if (grepl("Achn",rownames(readsMapped[i,])) ) {
    lines(xl, predict(loess_fit,xl), col = "green",lwd=2)
    points(gene ~ time, d,col="green")
  }
  else{
    lines(xl, predict(loess_fit, xl), col = "blue",lwd=2)
    points(gene ~ time, d,col="blue")
  }
}

warnings()

gene <- as.numeric(as.vector(readsMapped[2,]))
gene <- t(t(gene)/scales)
d<- data.frame( time = time,gene =gene )
d <- d[order(time),]
loess_fit <- loess( gene~ time , d)
xl <- seq(min(time),max(time), (max(time) - min(time))/5)
xl
predict(loess_fit,xl)

?loess
max(readsMapped[,23])

for (i in 1:length(series[,1]) {
  
    series <- readsMapped[readsMapped$time == i, ]
    print(series)
    loess_fit <- loess( point ~ time, series)
    #lines(series$x, predict(loess_fit), col = "blue")
    #lines(startPoint:endPoint, series[i,], col=col, lwd=lwd)
  }
  
  #randomIndex <- sample(1:length(series[,1]), 1)
  #lines(startPoint:endPoint, series[randomIndex,], col=hcol, lwd=hlwd)
  
  if (showPoints) {
    points(startPoint:endPoint, series[randomIndex,], col=hcol, pch=4)
  }
}





################################
#time <- c(24,72,72,48,72,1.5,1.5,1.5,3,3,3,6,6,6,12,12,12,24,24,48,48,120,120,120,96,96,96)
#readsMapped <- readsMapped[,c(3,6,8,10,12,16,17,18,25,26,27,28,29,30,34,37,38,41,42,43,44,49,50,51,55,56,57)]
#scales <- c(37.347786,38.243575,42.958841,37.693273,38.833448,38.767085,39.531887,66.581748,58.79224,36.833425,36.766411,38.43041,41.391057,38.195016,42.786511,43.100056,39.045232,41.279979,38.454703,41.514493,37.64413,37.562333,37.627361,42.407026,42.816286,36.399431,42.436111,41.946002,39.902334,42.970703,39.638666,40.412914,39.782447,41.603102,39.267669,38.407583,44.687265,38.630291,38.819944,43.256377,42.243905,39.799117,38.818551,42.459334,44.892463,29.408106,30.442288,38.035275,33.161491,37.71993,31.988759,35.262211,32.950623,41.12849,36.21799,34.777404,31.711137)
#scales <- scales[c(3,6,8,10,12,16,17,18,25,26,27,28,29,30,34,37,38,41,42,43,44,49,50,51,55,56,57)]
#scales