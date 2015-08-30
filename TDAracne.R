source("http://bioconductor.org/biocLite.R")
setwd("/home/ben/workspace/timeCourse/data")
library("TDARACNE")
library("igraph")
source("toCytoscape.R")

######## on sample data
subject<-"plant"  #should be plant or bacteria
track<-"psa"  #should be psa or control, i.e. treatment

### design
design<-getDataRange(24,120,track)
design<-data.frame(condition<-design$condition,time<-design$time)

## the genes of interest. 
datafile<-"/home/ben/workspace/timeCourse/data/networks/genesOfInterest.csv"
networkGenes <- read.csv(datafile,row.names=1)
geneSample<-networkGenes[sample(nrow(networkGenes),20), ]

##Expression set and run.
eset<-createESet(geneSample,design)
graph<-TDARACNE(eset,11,"netIRMAon",delta=3,likehood=1.2,norm=2,logarithm=1,thresh=0,ksd=0,0.15,adj=TRUE,adj =TRUE); #adj to return mim
write.csv(graph,"aracrneAdj.csv") #writes a mutual information matrix to file. 

#read mutual information mimatrix in from file. 
aracneMim<-read.csv("/home/ben/workspace/timeCourse/data/networks/aracneAdj7.csv",row.names=1)
aracneMim<-as.matrix(aracneMim)

graph<-graph_from_adjacency_matrix(aracneMim, mode = "undirected", weighted = NULL, diag = TRUE, add.colnames = NULL, add.rownames = NA)
plot(graph)


## remove isolated nodes. 

isolates <- which(igraph::degree(graph) == 0)
graph<-delete.vertices(graph, isolates)
plot(graph)

##Send graph to cytoscape. 
## graph properties
graph$name = "288 TDAracne"
graph$density = graph.density(graph)
V(graph)$degree <- igraph::degree(graph)
V(graph)$closeness <- closeness(graph)
V(graph)$betweenness <- betweenness(graph)
V(graph)$page_rank <- page.rank(graph)$vector
V(graph)$community <- label.propagation.community(graph)$membership
E(graph)$betweenness <- edge.betweenness(graph)


# Convert to Cytoscape style JSON object
cygraph <- toCytoscape(graph)

# Send it to Cytoscape!
network.url = paste(base.url, "networks", sep="/")
res <- POST(url=network.url, body=cygraph, encode="json")
network.suid = unname(fromJSON(rawToChar(res$content)))







