library(stats)
library(limma)

# find where genes in a set of gene sets lie your experimental results and plot barcode plots based on experimental t-values

print("How many degrees of freedom?")
df<-as.numeric(readLines(con = stdin(), n = 1, ok = TRUE))
# define cutoffs for significance in barcodeplots
co<-qt(c(0.05,0.95),df) # from t -distribution


# Identify gene sets relevant to your question
# This could be through the gene set search function of GeneSetDB at http://genesetdb.auckland.ac.nz
# Then load the gene sets - format to be as in Gene SetDB download format, tab delimited text file, has header row, column 1 is Class of gene set (e.g. Disease/Phenotype, GO, etc.), column 2 is Set Name, column 3 is SOurce DB (e.g. GP_BP, MPO), column 4 is Gene # (how many genes in set), and column 5 is Gene Names (of genes in set, as human OGS). 
print("Load file of Gene Sets you want to search")
FileChoice <- file.choose() # choose input file;
GeneSets <- read.delim(FileChoice, skip = 0, sep = "\t", as.is = TRUE)
GeneSets
nSets<-nrow(GeneSets)
nSets

# Load in the experimental results you want to querry - this is a t statistic - two columns - the first column is gene IDs converted to human OGS, the second column is the statistic
print("Load file of experimental results you want to use")
FileChoice <- file.choose() # choose input file;
ExpRes <- read.delim(FileChoice, skip = 0, sep = "\t", as.is = TRUE)
head(ExpRes)
ExpRes$SYMBOL <- tolower(ExpRes$SYMBOL)
colnames(ExpRes)[2]<-"T"
colnames(ExpRes)[1]<-"SYMBOL"

?lapply
sqrt()
ExpRes[]
ExpRes[4079,]
ExpRes[,2]

replacement <- sapply(ExpRes[,2], convertWaldtoT)
replacement
ExpRes <-ExpRes1
head(ExpRes1)
head(ExpRes)
ExpRes1$T <- replacement
ExpRes1 <- data.frame(ExpRes[2], lapply(ExpRes[,2], bob))


?data.frame
convertWaldtoT <- function(x) {
  if (x < 0) {
    -sqrt(abs(x))
  } else {
    sqrt(x)
  }
}

}
?sqrt


pdf(file=paste("Gene Set barcodeplots",date(),".pdf",sep=""),width=12,height=6)

write.table(ExpRes[,1],file="permutations/goLists/allPossible.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)


for (counter in 1:nSets) {
  fileName <- gsub("[^0-9]","",GeneSets[counter,2])
  print(fileName)
  #name <- grep(pattern="GO:(\d+)",)
  print(paste("Gene Set is",GeneSets[counter,3],": ",GeneSets[counter,2]))
  #lowSet<-	 tolower(unlist(strsplit(GeneSets[1, 5], ",")))
	lowSet<-   tolower(unlist(strsplit(GeneSets[counter, 5], ",")))
  m<-match(lowSet,ExpRes$SYMBOL)
  m <- m[!is.na(m)]

  subset<-ExpRes[m,]
  sig_down_in_set<-subset[subset$T<=co[1],1]
  sig_up_in_set<-subset[subset$T>=co[2],1]
  
  if (round(geneSetTest(match(lowSet,ExpRes$SYMBOL), ExpRes[,2], alternative="mixed", type="auto", ranks.only=FALSE, nsim=9999),4) < 0.05){
    write.table(subset[,1],file=paste("permutations/goLists/",fileName,".txt",sep=""),quote=FALSE, row.names=FALSE)
    sub= paste("Genes in set sig up =", paste(sig_up_in_set, collapse=", "),"\n","Genes in set sig down =", paste(sig_down_in_set, collapse=", "))
    #main=paste(GeneSets[1,3],": ",GeneSets[1,2],"\n \n Is the set of genes statistically highly ranked relative to other genes?    p =", round(geneSetTest(match(lowSet,ExpRes$SYMBOL), ExpRes[,2], alternative="mixed", type="auto", ranks.only=FALSE, nsim=9999),4))
    main=paste(GeneSets[counter,3],": ",GeneSets[counter,2],"\n \n Is the set of genes statistically highly ranked relative to other genes?    p =", round(geneSetTest(match(lowSet,ExpRes$SYMBOL), ExpRes[,2], alternative="mixed", type="auto", ranks.only=FALSE, nsim=9999),4))
    barcodeplot(ExpRes[,2], match(lowSet,ExpRes$SYMBOL), index2=NULL, quantiles=co, main=main, cex.main=0.8,,labels=c("highest  t  statistic","lowest t statistic"),sub=sub, cex.sub=0.5)
  }
}

dev.off()

