library(timecourse)
## load data into rawCounts
rawCounts = read.table( "subset.csv",sep=",", header=TRUE,row.names=1 , stringsAsFactors = FALSE )

## select timepoints
counts <- rawCounts[,c(2:12,39:57)] #24hr cycle timepoints

#set up parmeters
M3 <- as.matrix(counts[1:366,]) 
size <- matrix(3, nrow=366, ncol=2) 
assay<-c("rep1","rep1","rep1","rep1","rep1","rep2","rep2","rep2","rep1","rep3","rep3","rep2","rep3","rep2","rep3","rep2","rep3","rep3","rep1","rep2","rep3","rep1","rep2","rep3","rep1","rep2","rep3","rep1","rep2","rep3")
trt<-c("Control","Psa","Control","Control","Psa","Control","Psa","Control","Psa","Control","Psa","Control","Control","Psa","Psa","Psa","Psa","Control","Control","Control","Control","Psa","Psa","Psa","Control","Control","Control","Psa","Psa","Psa")
time<-c(72,24,48,24,72,72,72,48,48,48,72,24,24,24,24,48,48,72,120,120,120,120,120,120,96,96,96,96,96,96)
genenames<-rownames(counts)




#run and output top 20.
MB.paired <- mb.long(M3, method="paired", times=5, reps=size, condition.grp=trt, rep.grp=assay, time.grp=time)
par(mfrow=c(5,4))
for (cycle1 in 1: 20) {
  plotProfile(MB.paired, type="b", gnames=genenames,  pch=c("1","2","3"), xlab="Hour",ranking=cycle1)
}


#reset. The first 11 samples are now being listed last in the data frame. 
counts <- rawCounts[,c(39:57,2:12)] #24hr cycle timepoints

#set up parameters
M3 <- as.matrix(counts[1:336,]) 
size <- matrix(3, nrow=336, ncol=2) 
# the first 11 parameters have been moved to the corresponding position - at the end of the parameters lists. 
assay<-c("rep2","rep3","rep2","rep3","rep2","rep3","rep3","rep1","rep2","rep3","rep1","rep2","rep3","rep1","rep2","rep3","rep1","rep2","rep3","rep1","rep1","rep1","rep1","rep1","rep2","rep2","rep2","rep1","rep3","rep3")
trt<-c("Control","Control","Psa","Psa","Psa","Psa","Control","Control","Control","Control","Psa","Psa","Psa","Control","Control","Control","Psa","Psa","Psa","Control","Psa","Control","Control","Psa","Control","Psa","Control","Psa","Control","Psa")
time<-c(24,24,24,24,48,48,72,120,120,120,120,120,120,96,96,96,96,96,96,72,24,48,24,72,72,72,48,48,48,72)

#run and output top 20 ranked genes. 
MB.paired <- mb.long(M3, method="paired", times=5, reps=size, condition.grp=trt, rep.grp=assay, time.grp=time)
par(mfrow=c(4,3))
for (cycle1 in 1: 12) {
  plotProfile(MB.paired, type="b", gnames=genenames,  pch=c("1","2","3"), xlab="Hour",ranking=cycle1)
}

