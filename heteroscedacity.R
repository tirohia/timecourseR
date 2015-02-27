library("lmtest")
setwd ("/home/ben/workspace/timeCourse/data")
source("../code/R/tools.R")


?bptest()
lm(y~x1+x2)->p
?bptest
?voom
?limma
coeftest(p,vcov=hccm(p))

plot(resid(model)~fitted(model),B)

datafile="../results/tripleDE.csv"

## load data into rawCounts
counts = read.table( datafile,sep=",", header=TRUE,row.names=1 , stringsAsFactors = FALSE )

?bptest()

df1 <-data.frame(condition = factor(c(rep(c("control","psa"), each=3)) ),time=factor(c(rep(design$time, each=3))), counts = c(t(data[i,])))
plot(df1$time,df1$counts)
model<-lm(counts~time , df1)
abline(model,col="blue")
test<-bptest(model)
if(test$p.value <= 0.05){
  print(i)
  print(test$p.value)  
}
model<-update(model, correlation=corAR1())
AICc(model1,model,model2,model3,model4)
?AICc
summary(model3)
summary(model4)
plot(ACF(model),alpha=0.05)
plot(ACF.gls(model4),alpha=0.05)
summary(model)
install.packages("MuMIn")
library("MuMIn")
library("nlme")
?gls
?data.frame
colnames=c("gene","pvalue")
nonConstant<-data.frame()
dim(nonConstant)
model
for (i in triples){
  df1 <-data.frame(condition = factor(c(rep(c("control","psa"), each=3)) ),time=factor(c(rep(design$time, each=3))), counts = c(t(data[i,])))
  model1<-gls(counts~time ,df1)
  test<-bptest(model)
  if(test$p.value <= 0.05){
    nonConstant[i,"pvalue"]<-test$p.value
    print(i)
    print(test$p.value)  
  }
}

model2<-glm(counts~time ,df1,family="poisson")
model3<-glm(counts~time ,df1,family="gaussian")
model4<-lm(counts~time ,df1)
triples
summary(model4)
length(triples)
df1
plot(df1$time,df1$counts)
resid<-model$residuals
df1 <- data.frame(condition = factor(c(rep(c("control","psa"), each=3)) ),time=factor(c(rep(design$time, each=3))), counts = c(t(data[1,])))
resid
triples <- rownames(counts)
triples

meanSDPlot(as.matrix(counts[triples,]))
datafile="cornellGenomeModelTH-TCFrequencyMatrix.csv"
counts = read.table( datafile,sep=",", header=TRUE,row.names=1 , stringsAsFactors = FALSE )
dim(counts[triples,])


design <- getDataRange(rawCounts,24,120)
design
datacols<-getDataCols(design)
datacols
data<-counts[triples,datacols]
data<-log(data+1)
data<-cleanDataNames(data)

cols<-design[design$time==24 & design$condition=="psa",c("X1","X2","X3")]
cols
data

plot(hx, x, type="l", lty=2, xlab="x value", ylab="Density", main="Comparison of t Distributions")
boxplot(mpg~duration,data=data, main="Car Milage Data",      xlab="Number of Cylinders", ylab="Miles Per Gallon")
boxplot(x)

mod<-lm()

df1
require(ggplot2)
require(reshape2)
library("RColorBrewer")

head(data)
design
df1$time<-design$time

datacols
data[2,]
t(data[3,])
?vsn()
meanSdPlot(df1)
meanSdPlot(as.matrix(counts))
df1
?vst
library("vsn")
library("DESeq2")
head(df1)
x<-df$counts
df1
df1 <- data.frame(condition = factor(c(rep(c("control","psa"), each=3)) ),time=factor(c(rep(design$time, each=3))), counts = c(t(data[4,])))
ggplot(df1,aes_string(x="time")) + geom_density(aplha=.3)
ggplot(df1, aes(x=counts,  fill=condition)) + geom_density(alpha=.3) #+scale_fill_manual(values=c("green","blue"))
head(df1)
ggplot(df1, aes(x=time,y=counts, fill=condition)) + geom_histogram(binwidth=.5, position="dodge")
?geom_density
class(df1$time)
