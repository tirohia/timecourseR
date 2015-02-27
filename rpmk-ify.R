data<-file.choose()

counts
candidates = read.table( datafile,sep=",", header=TRUE, row.names=1,stringsAsFactors = FALSE, fill=TRUE,quote="")
list<-rownames(candidates)


datafile="/home/ben/workspace/timeCourse/data/actinidia/cornellGenomeModelTH-TCFrequencyMatrix.csv"
counts = read.table( datafile,sep=",", header=TRUE, row.names=1,stringsAsFactors = FALSE )

lengthfile<-"/home/ben/workspace/models/cornell/cornellGeneLengths.csv"
lengthfile<-"/home/ben/workspace/timeCourse/data/psa/NZv13invitro.csv"
lengths= read.table( lengthfile,sep=",", header=FALSE, row.names=1,stringsAsFactors = FALSE ,fill=TRUE, quote="")
lengths[,1]
design <- getDataRange(counts,72,72,"psa")
design

listdesign <- as.list(t(design[,3:5]))
listdesign<-unlist(listdesign)
listdesign

counts<-counts[,listdesign]
counts
sizes<-colSums(counts)
sizes
d<-DGEList(counts=as.matrix(counts),lib.size=sizes)
d
d$genes$Length<-lengths[,1]
wibble<-rpkm(d)


dim(wibble)
head(wibble)
max(wibble)

listlist<-unlist(as.list(wibble))
breaklist<-seq(0,1000,5)
breaklist
fac<-factor(cut(listlist,breaks=breaklist))
#Tabulate and turn into data.frame
xout <- as.data.frame(table(fac))
xout
scatterplot(Freq~fac, data=xout,smoother=loessLine)

chisq.test(table(fac)) 

?chisq.test


#Add cumFreq and proportions
#xout <- transform(xout, cumFreq = cumsum(Freq), relative = prop.table(Freq))
xout

attach(mtcars)
plot(wt, mpg, main="Scatterplot Example",  xlab="Car Weight ", ylab="Miles Per Gallon ", pch=19)
wt
mpg
mtcars
