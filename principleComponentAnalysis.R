library("car")
library(gclus)


par(mfrow=c(1,1))
counts<-na.omit(rawCounts)
conditions = c("Control","Psa","Control","Psa","Control","Control","Control","Psa","Psa","Psa","Control","Control","Control","Control","Control","Control","Control","Psa","Psa","Psa","Psa","Psa","Psa","Psa","Control","Control","Control","Psa","Control","Control","Psa","Psa","Control","Control","Control","Psa","Psa","Psa","Psa","Control","Control","Control","Control","Psa","Control","Psa","Psa","Control","Psa","Control")
colnames(counts) <- conditions


rld <- rlog(dds)
print(plotPCA(rld, intgroup=c("condition", "time"),ntop=500,col=c("red","blue","brown","springgreen","antiquewhite","cadetblue","coral","hotpink","azure","black","blueviolet","cyan","darkblue","darkgrey","green","chartreuse4","darkseagreen1","gold")))
print(plotPCA(rld, intgroup=c("condition", "time"),ntop=500,col=c("red","red","red","red","red","red","red","red","red","red","blue","blue","blue","blue","blue","blue","blue","blue")))
print(plotPCA(rld, intgroup=c("condition", "time"),ntop=500,col=c("red","black","pink","blue","brown","green","grey","yellow","yellow","cyan","black","pink","blue","brown","green","grey","yellow","cyan")))
plotPCA(vsd, intgroup=c("condition", "time"),ntop=500)
?plotPCA

vsd <- varianceStabilizingTransformation(dds)
vsd
print(plotPCA(vsd, intgroup=c("condition", "time"),ntop=500,col=c("red","blue","brown","springgreen","antiquewhite","cadetblue","coral","hotpink","azure","black","blueviolet","cyan","darkblue","darkgrey","green","chartreuse4","darkseagreen1","gold")))
print(plotPCA(vsd, intgroup=c("condition", "time"),ntop=500,col=c("red","red","red","red","red","red","blue","blue","blue","blue","blue")))
print(plotPCA(vsd, intgroup=c("condition", "time"),ntop=500,col=c("red","black","pink","blue","brown","green","grey","yellow","cyan")))
print(plotPCA(vsd, intgroup=c("condition", "time"),ntop=500,col=c("red","black","pink","blue","brown","green","black","pink","blue","brown","green")))

library("RColorBrewer")
library("gplots")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, time, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))


head(counts)
tail(counts)

CF <- data.frame(counts)
cf <- data.frame(t(CF))
cf
counts <- cf

fit <- princomp(counts, cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)
?biplot

scatterplotMatrix(counts[2:6])
plot(counts$Psa120, counts$Psa24)

sapply(counts[2:57],sd)

standardisedconcentrations <- as.data.frame(scale(counts[1:39000]))
standardisedconcentrations1 <- data.frame(t(na.omit(t(standardisedconcentrations))))

counts.pca <- prcomp(standardisedconcentrations1,retx=TRUE)

summary(counts.pca)
counts.pca$sdev

screeplot(counts.pca, type="lines")
plot(counts.pca$x[,1],counts.pca$x[,2])
text(counts.pca$x[,1],counts.pca$x[,2], counts$Psa120, cex=0.7, pos=4, col="red")
text(counts.pca$x[,1],counts.pca$x[,2], counts$Psa24, cex=0.7, pos=4, col="blue")

my.abs     <- abs(cor(counts[,-1]))
my.colors  <- dmat.color(my.abs)
my.ordered <- order.single(cor(my.wines[,-1]))
cpairs(my.wines, my.ordered, panel.colors=my.colors, gap=0.5)

?prcomp
counts.pca <- prcomp(counts, center=TRUE, scale=TRUE)

screeplot(counts.pca, main="Scree Plot", xlab="Components")
screeplot(counts.pca, main="Scree Plot", type="line" )

load    <- counts.pca$rotation
sorted.loadings <- load[order(load[, 1]), 1]
myTitle <- "Loadings Plot for PC1" 
myXlab  <- "Variable Loadings"
dotplot(sorted.loadings, main=myTitle, xlab=myXlab, cex=1.5, col="red")

sorted.loadings <- load[order(load[, 2]), 2]
myTitle <- "Loadings Plot for PC2"
myXlab  <- "Variable Loadings"
dotplot(sorted.loadings, main=myTitle, xlab=myXlab, cex=1.5, col="red")

biplot(counts.pca, cex=c(1, 0.7))



library(RColorBrewer)
names <- conditions
mylist <- list(counts$Psa120,counts$Psa24,counts$Psa48,counts$Psa96,counts$Psa72)
makeProfilePlot(mylist,names)








makeProfilePlot <- function(mylist,names)
{
  require(RColorBrewer)
  # find out how many variables we want to include
  numvariables <- length(mylist)
  # choose 'numvariables' random colours
  colours <- brewer.pal(numvariables,"Set1")
  # find out the minimum and maximum values of the variables:
  mymin <- 1e+20
  mymax <- 1e-20
  for (i in 1:numvariables)
  {
    vectori <- mylist[[i]]
    mini <- min(vectori)
    maxi <- max(vectori)
    if (mini < mymin) { mymin <- mini }
    if (maxi > mymax) { mymax <- maxi }
  }
  # plot the variables
  for (i in 1:numvariables)
  {
    vectori <- mylist[[i]]
    namei <- names[i]
    colouri <- colours[i]
    if (i == 1) { plot(vectori,col=colouri,type="l",ylim=c(mymin,mymax)) }
    else         { points(vectori, col=colouri,type="l")                                     }
    lastxval <- length(vectori)
    lastyval <- vectori[length(vectori)]
    text((lastxval-10),(lastyval),namei,col="black",cex=0.6)
  }
}

