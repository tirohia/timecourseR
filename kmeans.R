source("http://bioconductor.org/biocLite.R")
library("timecourse")
library(gplots)

setwd("/home/ben/workspace/timeCourse/data")
datafile="cornellGenomeModelTH-TCFrequencyMatrix.csv"

## load data into rawCounts
datafile<-file.choose()
rawCounts = read.table( datafile,sep=",", header=TRUE, row.names=1,stringsAsFactors = FALSE )
dim(rawCounts)
colnames(rawCounts)
counts<-log(rawCounts+1)



#### log 2 transformed
countsLog2<-log2(rawCounts +1 ) #optionally log2-transformn the data
countsLog2 <- countsLog2[complete.cases(countsLog2),]
nrow(countsLog2)
?kmeans

c<-kmeans(countsLog2,k)
s<-cbind(c$cluster,rawCounts)


for (k in 2:10){
  c<-kmeans(countsLog2,k)
  s<-cbind(c$cluster,rawCounts)
  
  colnames(s)[1]<-"CL"
  lb1<-length(rawCounts[1,])
  x<-1:lb1
  s.ordered<- s[order(s$CL),]; s.ordered<-data.frame(s.ordered)
    
  png(file=paste("k_means_graphs",k,".png", sep="") ,width=1900,height=1080)
  par(mfrow=c(2,2))
  
  if(k > 4) {par(mfrow=c(2,3))}
  if(k > 6) {par(mfrow=c(3,3))}
  if(k > 9) {par(mfrow=c(3,4))}
  
  
  for (cycle1 in 1: k) {
    r<- s.ordered[s.ordered$CL==cycle1,-1]
    l=length(r[,1])
    
    fileConn<-file(paste("../results/kmeans",k,"-",cycle1,".txt",sep=""))
    writeLines(rownames(r), fileConn)
    close(fileConn)
    for (cycle2 in 1: l-1) {
      plot(x,r[cycle2,],col=cycle1+1,type="l",ylab="",ylim=range(min(r),max(r)) ,cex.axis=0.5)
      par(new=TRUE)
    }
    plot(x,r[l,],col=cycle1+1,type="l",ylab="",ylim=range(min(r),max(r)),main=paste("k =",cycle1) ,cex.axis=0.5)
  }

  dev.off()
}



##############################################################################################

#Finding a useful k


k=4
data<-countsLog2[1:1000,]
data<-rawCounts[1:1000,]
data<-counts1
c<-kmeans(data,k)
s<-cbind(c$cluster,data)

colnames(s)[1]<-"CL"
head(s)
lb1<-length(data[1,])

x<-1:lb1
s.ordered<- s[order(s$CL),]; s.ordered<-data.frame(s.ordered)

par(mfrow=c(2,2))

for (cycle1 in 1: k) {
  r<- s.ordered[s.ordered$CL==cycle1,-1]
  l=length(r[,1])
  for (cycle2 in 1: l-1) {
    plot(x,r[cycle2,],col=cycle1+1,type="l",ylab="",ylim=range(min(r),max(r)) ,cex.axis=0.5)
    par(new=TRUE)
  }
  plot(x,r[l,],col=cycle1+1,type="l",ylab="",ylim=range(min(r),max(r)),main=paste("k =",cycle1) ,cex.axis=0.5)
}


### using the elbow method.
par(mfrow=c(1,1))

wss <- (nrow(data)-1)*sum(apply(data,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(data,centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",   ylab="Within groups sum of squares")

### using BIC
library("mclust")
# Run the function to see how many clusters
# it finds to be optimal, set it to search for
# at least 1 model and up 20.
d_clust <- Mclust(as.matrix(t(log(counts+1))), G=1:10)
m.best <- dim(d_clust$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")
# 4 clusters
plot(d_clust)

#### BIC method 2

mydata_for_clustering<- counts
dim(mydata_for_clustering)
dim(counts)
n <- nrow(mydata_for_clustering)
nclus = 7
myclus = kmeans(mydata_for_clustering,centers=nclus)
print(names(myclus))
myclus$cluster[myclus$cluster==3]


library("sfsmisc")
mult.fig(1,main="Simulated data with two clusters")
scol = rainbow(nclus,end=0.8) # select colors from the rainbow

library(cluster)
library("fpc")

plotcluster(counts, myclus$cluster,col=scol[myclus$cluster],pch=2)
clusplot(counts, myclus$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)

kmax = 15 # the maximum number of clusters we will examine; you can change this
totwss = rep(0,kmax) # will be filled with total sum of within group sum squares
kmfit = list() # create and empty list
for (i in 1:10){
  kclus = kmeans(mydata_for_clustering,centers=i,iter.max=20)
  print(kclus$tot.withinss)
  totwss[i] = kclus$tot.withinss
  print(ncol(kclus$centers))
  print(length(kclus$cluster))
  print(nrow(kclus$centers))
  kmfit[[i]] = kclus
}

dim(counts)
kmeansBIC = function(fit){
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return( D + 0.5*log(n)*m*k)
}
kmeansAIC = function(fit){
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + 2*m*k)
}
kmfit[2]$tot.withinss
bic=sapply(kmfit,kmeansBIC)
aic=sapply(kmfit,kmeansAIC)
bic
aic

mult.fig(1,main="Simulated data with two clusters")

v
v = -diff(bic)
nv = length(v)
fom = v[1:(nv-1)]/v[2:nv]
nclus = which.max(fom)+1
cat("The apparent number of clusters is: ",nclus,"\n")
plot(1:kmax,bic[1:kmax],main="Minimizing the BIC",xlab="k", ylab="BIC",type="b")
?plot


#### G clustering
head(counts)
plot(density(counts))
shapiro.test(counts);
qqline(counts, col = 2)
