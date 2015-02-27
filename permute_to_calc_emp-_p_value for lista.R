# permute_to_calc_emp-_p_value
library(gdata)
print("")
print("enter first a list of all possible genes, then the differentially expressed gene list, then the biological gene set list")
print("")

setwd("/home/ben/workspace/stemAssay/code")

print("read in one column tab delim text fie with  COLUMN HEADER  of all genes that could have been in the gene list - eg. the file named possible_targets.txt")
FileChoice<-file.choose() # choose input file;
allGenes <-read.delim(FileChoice, skip=0, sep="\t", as.is=TRUE)
allGenes <-as.character(allGenes[,1])
allGenes

print("read in one column tab delim text fie with  COLUMN HEADER of differentially expressed genes e.g. the file named actual5.txt")
FileChoice<-file.choose() # choose input file;
deGenes <- read.delim(FileChoice, skip=0, sep=",", as.is=TRUE)
deGenes <-as.character( deGenes[,1]) 
deGenes

print("read in one column tab delim text fie with  COLUMN HEADER  of probe sets in the biological gene set list e.g. the file named canonical_targets.txt")
FileChoice<-file.choose() # choose input file;
goGeneSet <- read.delim(FileChoice, skip=0, sep="\t", as.is=TRUE)
goGeneSet <-as.character(goGeneSet[,1])
goGeneSet

itterations <-10000 # use at least 10,000 itterations or this will not work well
numberOfDEGenes <-length(deGenes)

o<-vector("numeric",length=itterations)
o

?resample
for (cyr in 1:itterations) {
	print(paste("perm itter count = ",cyr))
	L1<-resample(allGenes, numberOfDEGenes ,replace=FALSE)	
  L1
	o[cyr]<-length(L1[L1 %in% goGeneSet])
  
}

Ly <- deGenes
Ly
actu <- length(Ly[Ly %in% goGeneSet])
actu
#actu<-length(Ly) #for debugging
		
emp_p <-length(o[o>=actu])/itterations
emp_p

emp_p_str<-paste("Empirical p ≤",emp_p)
emp_p_str


if(emp_p == 0){emp_p_str<-paste("empirical p ≤",1/itter)}

hh<-hist(o,breaks=itterations/800)
plot(hh,col="red", xlim=c(0,max(o,actu)+.1*max(o, actu)), ylim=c(0,max(hh$counts)+.1*max(hh$counts)),xlab="number of RNAs selected by chance in each of 10,000 random gene lists",ylab="frequency",main="")
arrows(actu,max(hh$counts), actu,.1*max(hh$counts),col="green",lwd=2)
abline(v=quantile(o , probs = c(0.05,0.95)),col="blue",lty=2)
