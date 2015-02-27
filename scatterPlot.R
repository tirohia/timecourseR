setwd("/home/ben/workspace/timeCourse/data")

datafile<-file.choose()
readsMapped = read.table( datafile,sep=",", header=TRUE,stringsAsFactors = FALSE)
head(readsMapped)

colnames(readsMapped)
columnnames <- colnames(readsMapped)
columnnames <- gsub("X", "P", columnnames)
colnames(readsMapped) <- columnnames
colnames(readsMapped)


total <- readsMapped$total
total
scales <- total/total[1]
scales

readsMapped$mapped <- readsMapped$mapped*scales
readsMapped$psa <- readsMapped$psa/scales
readsMapped$psa <- (readsMapped$psa/total)*100
readsMapped$mapped <- (readsMapped$mapped/total)*100

partialReads <- readsMapped[readsMapped$condition == "Psa", ]
head(partialReads)

library(car)
scatterplot(psaM ~ time| condition, data=readsMapped,   xlab="Time (HPI)", ylab="Reads/Million", main="Psa levels",reg.line=TRUE, smooth=TRUE)

scatterplot(gene ~ time |condition , data=partialReads,   xlab="Time (HPI)", ylab="Reads", main="Psa overlap gene levels",reg.line=TRUE, smooth=TRUE)

plot(P5 ~time , data=partialReads,   xlab="Time (HPI)", ylab="Reads/Million", main="Psa levels",reg.line=TRUE, smooth=TRUE)


abline(lm(psaM ~ time,data=control))

?glm
psaCount.lm <- mlm(count ~ time ,family = data=psa)

apply()
?apply
max(actM)

counts <- readsMapped$counts
time <- readsMapped$time
actM <- readsMapped$actM
condition <- readsMapped$condition

abline(glm.fit(time,actM + counts))

model <- glm(time ~ actM + condition,family=gaussian, data=readsMapped)
anova(model)
summary(model)
coefficients(psaCount.lm)
summary(psaCount.lm)
psaCount.res <- resid(psaCount.lm)
psaCount.stres <- rstandard(psaCount.lm)

plot(psa$time, psaCount.stres)
abline(0, 0)
abline(lm(count ~ time,data=psa))


readsMapped[readsMapped$condtion == "Psa",]
psa <- subset(readsMapped, condition =="Psa")
control <- subset(readsMapped, condition =="Control")
