biocLite("gprege")
library(gprege)

datafile<-file.choose()
rawCounts = read.table( datafile,sep=",", header=TRUE,row.names=1 , stringsAsFactors = FALSE )
rawCounts<-cleanDataNames(rawCounts)

head(rawCounts)
design <- getDataRange(rawCounts,24,120)

timecolumns<-NULL
for(j in design[design$condition=="psa",2]){
  timecolumns<-append(timecolumns,design[design$time==j & design$condition=="psa",3:5])
}
for(j in design[design$condition=="control",2]){
  timecolumns<-append(timecolumns,design[design$time==j & design$condition=="control",3:5])
}
names(timecolumns) <- NULL
counts<-rawCounts[,unlist(timecolumns)]
head(counts)
M3 <- as.matrix(counts[1:39039,])



con <- url("http://gprege.googlecode.com/svn/trunk/DellaGattaData.RData")
while(!exists("DGdata")) try(load(con),TRUE); close.connection(con)
# Timepoints / GP inputs.
tTrue = matrix(seq(0,240,by=20), ncol=1)
?matrix
gpregeOptions <- list()
# Explore individual profiles in interactive mode.
gpregeOptions$explore <- FALSE
# Exhaustive plot resolution of the LML function.
gpregeOptions$exhaustPlotRes <- 30
# Exhaustive plot contour levels.
gpregeOptions$exhaustPlotLevels <- 10
# Exhaustive plot maximum lengthscale.
gpregeOptions$exhaustPlotMaxWidth <- 100
# Noisy ground truth labels: which genes are in the top 786 ranks of the TSNI ranking.
gpregeOptions$labels <- DGatta_labels_byTSNItop100
# SCG optimisation: maximum number of iterations.
gpregeOptions$iters <- 100
# SCG optimisation: no messages.
gpregeOptions$display <- FALSE
# Matrix of different hyperparameter configurations as rows:
# [inverse-lengthscale percent-signal-variance percent-noise-variance].
gpregeOptions$inithypers <-  + matrix( c( 1/1000, 1e-3, 0.999, 1/8, 0.999, 1e-3, 1/80, 2/3, 1/3 ), ncol=3, byrow=TRUE)
gpregeOptions
exprs_tp63_RMA
gpregeOutput<-gprege(data=exprs_tp63_RMA,inputs=tTrue,gpregeOptions=gpregeOptions)
