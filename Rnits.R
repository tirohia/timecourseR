biocLite("knitr")
library("knitr")
library("Rnits")
setwd ("/home/ben/workspace/timeCourse/data")
## load data into rawCounts
datafile<-file.choose()
networkCounts = read.table( datafile,sep=",", header=TRUE,row.names=1 , stringsAsFactors = FALSE )
networkCounts<-cleanDataNames(networkCounts)
design <- getDataRange(rawCounts,24,120)
design
datacols<-getDataCols(design)

networkCounts<-networkCounts[,datacols]  

library(GEOquery)
library(stringr)

gds<-getGEO("GSE4158",AnnotGPL=FALSE)[[1]]
pdata<-pData(gds)
pData

filt<-pdata$characteristics_ch2 %in% names(which((table(pdata$characteristics_ch2 == 2))))
pdata$characteristics_ch2
rownames(dataDF)<-fData(dat)$ID
rownames(dataDF)
rnit <-build.Rnits(dataDF,probedata=probedata,phenodata=phenodata,logscale=TRUE,normmethod="Between")



#######################Venn diagram shit

library("venneuler")

require(plotrix)
draw.circle(4,14,2,border="black",col="white")



library("VennDiagram")
draw.pairwise.venn(8, 10, 4,category=c("p(B)","p(A)"))
?draw.pairwise.venn
grid.newpage();
venn.plot <- draw.pairwise.venn(
  area1 = 7,
  area2 = 7,
  cross.area = 3,
  category = c("p(A)", "p(B)"),
  fill = c("white", "white"),
  lty = "solid",
  cat.default.pos="text",
  #cex = 2,
  cat.cex = 2,
  #cat.pos = c(285, 105),
  #cat.dist = 0.09,
  #cat.just = list(c(-1, -1), c(1, 1)),
  ext.pos = 30,
  ext.dist = -0.05,
  ext.length = 0.85,
  ext.line.lwd = 2,
  ext.text=FALSE
);
grid.draw(venn.plot);


venn.plot <- draw.pairwise.venn(100, 70, 30, c("First", "Second"));
grid.draw(venn.plot);

venn.plot <- draw.pairwise.venn(100, 70, 30, c("First", "Second"), scaled = FALSE);
grid.draw(venn.plot);
grid.newpage();
venn.plot <- draw.pairwise.venn(
  area1 = 100,
  area2 = 70,
  cross.area = 0,
  category = c("First", "Second"),
  cat.pos = c(0, 180),
  euler.d = TRUE,
  sep.dist = 0.03,
  rotation.degree = 45
);

