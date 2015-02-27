source("http://bioconductor.org/biocLite.R")
library("RTN")


data(dt4rtn)
head(dt4rtn$gexp)
head(dt4rtn$gexpIDs)
head(dt4rtn$pheno)
head(dt4rtn$gexp)
head(dt4rtn$gexp)
head(dt4rtn$gexp)
?gexp
rtni <- new("TNI", gexp=dt4rtn$gexp, transcriptionFactors=dt4rtn$tfs[c("PTTG1","E2F2","FOXM1")] )
rtni
rtni<-tni.preprocess(rtni,gexpIDs=dt4rtn$gexpIDs)

rtni<-tni.permutation(rtni)
rtni
rtni<-tni.bootstrap(rtni)
rtni
rtni<-tni.dpi.filter(rtni)
rtni
tni.get(rtni,what="summary")
refnet<-tni.get(rtni,what="refnet")
refnet
tnet<-tni.get(rtni,what="tnet")
tnet
g<-tni.graph(rtni)
g
plot(g)
