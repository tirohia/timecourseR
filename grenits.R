source("http://bioconductor.org/biocLite.R")
library("GRENITS")
setwd ("/home/ben/workspace/timeCourse/data")

data(Athaliana_ODE)
head(Athaliana_ODE)

file<-file.choose()
data = read.table( file,sep=",", header=TRUE,row.names=1 , stringsAsFactors = FALSE )

data<-as.matrix(data)
data<-log(data+1)

output.folder <- paste(tempdir(), "testGranitsNet", sep="")
LinearNet(output.folder, data)
analyse.output(output.folder)
dir(output.folder)
