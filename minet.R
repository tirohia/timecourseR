library("minet")
library("Rgraphviz")
library("RCytoscape")
setwd ("/home/ben/workspace/timeCourse/data")
source("../code/R/tools.R")
library("igraph")
library("graph")


counts<-read.csv("networkCandidateGeneReads.csv",row.names=1,header=TRUE)
candidates<-read.csv("../results/tripleDE.csv",row.names=1,header=TRUE)
candidates<-rownames(candidates)
candidates
medians<-getMedians(counts[candidates,],24,120,condition="psa")
dim(medians)
psaMedians<-as.data.frame(t(medians))
psaMedians<-log(psaMedians+1)
psaMim<-build.mim(psaMedians, estimator = "spearman", disc = "none", nbins = sqrt(NROW(dataset)))
psaNet<-aracne(psaMim,eps=0)
psaNet[is.na(psaNet)] <- 0
psaNet
psaNet<-mrnet( psaMim)


am.graph<-new("graphAM", adjMat=psaNet, edgemode="undirected")

plot(am.graph, attrs = list(node = list(fillcolor = "lightblue"),edge = list(arrowsize=0.5)))


medians<-getMedians(candidates,24,120,condition="control")
controlMedians<-as.data.frame(t(medians))
controlMedians<-log(controlMedians+1)
mim<-build.mim(controlMedians, estimator = "spearman", disc = "none", nbins = sqrt(NROW(dataset)))
controlNet<-aracne(mim,eps=0)
controlNet[is.na(controlNet)] <- 0

am.graph<-new("graphAM", adjMat=controlNet[1800:2000,1800:2000], edgemode="undirected")
plot(am.graph, attrs = list(node = list(fillcolor = "lightblue"),edge = list(arrowsize=0.5)))

##subtract control from psa, see what changes. 
net<-psaNet-controlNet
?pdf()

net[net < 0] <- 0
dim(psaNet)
am.graph<-new("graphAM", adjMat=net[1800:2000,1800:2000], edgemode="undirected")
am.matrix <- as(am.graph,"matrix")
am.nel <-as(am.graph,"graphNEL")
graph <- graph.adjacency(am.matrix, mode="undirected")
class(graph)

plot(am.graph, attrs = list(node = list(fillcolor = "lightblue"),edge = list(arrowsize=0.5)))





### export for cytoscape

am.matrix<-as(am.graph,"graphNEL")
am.matrix
To <- as(am.matrix,"matrix")
To
?write.graph
write.graph(am.matrix,"graph.gml","gml")
g2 <- igraph.from.graphNEL(am.matrix)

myAdjacencyMatrix <- matrix(runif(400),nc=20,nr=20)
myAdjacencyMatrix
g  <- graph.adjacency(To,weighted=TRUE)
df <- get.data.frame(g)
head(df)
write.table(df,"graph.gml")










plot( as(net[1:100,1:100] ,"graphNEL") )
colnames(net)
class(net)
net <- minet(rawCounts, method = "mrnet", estimator = "mi.empirical", disc = "equalwidth", nbins = sqrt(nrow(rawCounts)))
net[1:10,1:10]

net[1:100,1:100]
am.graph



head(candidates)
head(thing)
rownames(medians)
cor(medians[,2:5])
plot(as(net ,"graphNEL"))
poi <- read.csv( "deResults/psa/psaDown1.5-12.csv", header=TRUE, row.names=1,stringsAsFactors = FALSE )
rownames(poi)
