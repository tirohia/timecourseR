setwd("/home/ben/workspace/timeCourse/data")
source("../code/R/tools.R")

## This should be turned into a function and put in tools.R

subject<-"plant"  #should be plant or bacteria
track<-"psa"  #should be psa or control, i.e. treatment
counts<-load.data(subject,track,1.5,120)
counts<-assay(rld)

#load the data - needs to be the design as the deseq analysis used to generate the rld assay if normalized reads are being used. 
design <- getDataRange(1.5,120,track)
genesOfInterest<-c("Achn212061","Achn095821","Achn357801","Achn168181") #list of genes that you want on the same plot. 
counts<-counts[genesOfInterest,]


# basic plot parameters/outlines. 
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
mainTitle<-"Levels of nominal housekeeping genes in innoculated plants"
xlabel<-"HPI"
ylabel<-"normalized log Reads"

#number of divisions that the loess curve prediction is going to be over. 
xl <- seq(min(design$time),max(design$time), (max(design$time) - min(design$time))/1000)

time<-rep(design$time,each=3)
plot(0,0, xlim=c(0, 120), ylim=c(0, max(reads)),xlab=xlabel,ylab=ylabel ,type="n",main=mainTitle,cex=2,, bty='L')

colour<-1
for(gene in genesOfInterest){
  reads<-counts[rownames(counts)==gene,] 
  df<-data.frame(time=time,reads=reads) #dataframe cosnisting of time point and corresponding reads.
  colnames(df)<-c("time","geneCounts")
  loess_fit <- loess(geneCounts ~ time ,df)
  lines(xl, predict(loess_fit,xl), col = colour ,lwd=2) 
  ##  first two variables are x/y positon of the legend -
  legend(10,1+0.5*colour,colour,inset=c(-.4,0),legend=gene,col = i*3,ncol = 2, cex = 0.75, lwd = 3, text.font = i, text.col = i*3,y.intersp=0,lty = c(1,2))
  colour<-colour+1
}

