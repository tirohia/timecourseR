
library("ShortRead")

#softTrim
#adapted from Jeremy Leipzig http://jermdemo.blogspot.co.nz/2010/03/soft-trimming-in-r-using-shortread-and.html
#and http://manuals.bioinformatics.ucr.edu/home/ht-seq#TOC-Quality-Reports-of-FASTQ-Files-
#trim first position lower than minQuality and all subsequent positions
#omit sequences that after trimming are shorter than minLength or longer than maxLength
#left trim to firstBase, (1 implies no left trim)
#input: ShortReadQ reads
#       integer minQuality
#       integer firstBase
#       integer minLength
#       integer maxLength
#output: ShortReadQ trimmed reads
softTrim<-function(reads,minQuality,firstBase=1,minLength=5,maxLength=900){
  #qualMat<-as(FastqQuality(quality(quality(reads))),'matrix')
  qualMat<-as(SFastqQuality(quality(quality(reads))),'matrix')
  qualList<-split(qualMat,row(qualMat))
  ends<-as.integer(lapply(qualList,function(x){which(x < minQuality)[1]-1}))
  #length=end-start+1, so set start to no more than length+1 to avoid negative-length
  starts<-as.integer(lapply(ends,function(x){min(x+1,firstBase)}))
  #use whatever QualityScore subclass is sent
  newQ<-ShortReadQ(sread=subseq(sread(reads),start=starts,end=ends),
                   quality=new(Class=class(quality(reads)),quality=subseq(quality(quality(reads)),start=starts,end=ends)),
                   id=id(reads))
  #apply minLength using srFilter
  minlengthFilter <- srFilter(function(x) {width(x)>=minLength},name="minimum length cutoff")
  trimmedReads = newQ[minlengthFilter(newQ)]
  maxlengthFilter <- srFilter(function(x) {width(x)<=maxLength},name="maximum length cutoff")
  trimmedReads = trimmedReads[maxlengthFilter(trimmedReads)]
  return(trimmedReads)
}

#readnumtrim
#call softTrim and plot number of reads passing filter for different quality/length parameters
#randomly sample "nsamples" reads in fastq file and perform "nrep" replicates of filtering+counting
#create a plot if "do_plot" and write it to "pdfout"_max and "pdfout"_min files
#use quantiles of length and quality distribution for axis ticks if "quant"
#if "maxth" calculate for both maximum length threshold and minimum length threshold, otherwise only do minimum length 
#input: fastq file format
#       integer nsamples
#       integer nrep
#       boolean do_plot
#       character pdfout
#       boolean quant
#output: matrix of number of reads for each quantile of the distribution of quality scores and read lengths
datafile<-file.choose()
readnumtrim(datafile)
readnumtrim<-function(fastqfile,nsamples=100,nrep=100,do_plot=TRUE,pdfout="",quant=FALSE,maxth=FALSE){
  tnr = as.numeric(system(paste("cat",fastqfile,"|wc -l"),intern=TRUE))/4
  cat("total number of reads:",tnr,"\n")
  
  # use max 1e6 samples to estimate distributions of read length and quality scores
  if (tnr<1e6){
    reads <- readFastq(fastqfile, qualityType="Auto")
  }else{
    f <- FastqSampler(fastqfile,1e6)
    reads <- yield(f)
    close(f)
  }
  #qual = FastqQuality(quality(quality(reads))) # get quality scores
  qual = SFastqQuality(quality(quality(reads))) # get quality scores
  readM = as(qual,"matrix")
  max_qual = max(readM,na.rm=TRUE)
  max_length = max(width(reads))
  
  if (quant){
    quantile_seq = seq(0,1,length.out=10)
    length.ticks = round(quantile(width(reads),quantile_seq)) # get read lengths
    qual.ticks = round(quantile(as.numeric(readM),quantile_seq,na.rm=TRUE))
  } else {
    length.ticks = round(seq(0,max_length,length.out=10))
    qual.ticks = round(seq(0,max_qual,length.out=10))
  }
  rm(reads)
  rm(qual)
  rm(readM)
  
  # get subsamples to estimate number of reads for different pairs of quality and length threshold
  f <- FastqSampler(fastqfile,nsamples)
  
  numreadM = array(0, dim=c(10,10,nrep))
  numreadm = array(0, dim=c(10,10,nrep))
  for (n in 1:nrep){
    reads <- yield(f)
    
    if (maxth){
      #### QUALITY VS MAXIMUM LENGTH
      for (lp in 1:length(length.ticks)){
        for (qp in 1:length(qual.ticks)){
          tr = softTrim(reads=reads,
                        minQuality=as.numeric(qual.ticks[qp]),
                        firstBase=1,
                        minLength=1,
                        maxLength=as.numeric(length.ticks[lp]))
          numreadM[lp,qp,n] = length(tr)
        }
      }
    }
    
    #### QUALITY VS MINIMUM LENGTH
    #following does not work because need to vectorise softTrim() see http://stackoverflow.com/questions/5554305/simple-question-regarding-the-use-of-outer-and-user-defined-functions
    # call_softTrim <- function(x,y) {
    #   tr=softTrim(reads=reads,minQuality=as.numeric(qual.ticks[x]),firstBase=1,minLength=as.numeric(length.ticks[y]),maxLength=max(width(reads)))
    #   return(length(tr))
    # }
    # numread = outer(1:length(length.ticks), 1:length(qual.ticks),call_softTrim)
    for (lp in 1:length(length.ticks)){
      for (qp in 1:length(qual.ticks)){
        tr = softTrim(reads=reads,
                    minQuality=as.numeric(qual.ticks[qp]),
                    firstBase=1,
                    minLength=as.numeric(length.ticks[lp]),
                    maxLength=max_length)
        numreadm[lp,qp,n] = length(tr)
      }
    }
  }
  close(f)
  
  
  # average over replicates
  anumreadM = apply(numreadM,c(1,2),mean)
  anumreadm = apply(numreadm,c(1,2),mean)
  
  if (do_plot){# plot with colours
    cpalette = colorRampPalette(c("white","blue"))
    if (maxth){
      # maximum length threshold
      cpalette = colorRampPalette(c("green","red"))
      if (nchar(pdfout)>0){
        pdf(paste(pdfout,"_max.pdf",sep=""))
      }else{
        x11()
      }
      filled.contour2(seq(0,1,length.out=10),seq(0,1,length.out=10),t(anumreadM),
                      axes=FALSE,xlab="minimum quality",ylab="maximum length",color.palette=cpalette)
      title(paste("maximum length threshold vs base quality\n",fastqfile,"\n","total reads",tnr),cex.main=0.7)
      contour(seq(0,1,length.out=10),seq(0,1,length.out=10),t(anumreadM),axes=FALSE,add=T,levels=seq(0,nsamples,length.out=10),
              labels=paste(as.character(round(seq(0,tnr,length.out=10)*100/tnr)),"% - ",as.character(round(seq(0,tnr,length.out=10))),sep=""))
      #axis(1,at=seq(0,1,length.out=10),label=qual.ticks)
      axis(1,at=seq(0,1,length.out=10),label=paste(qual.ticks,round(exp(qual.ticks/(-10)),digits=3),sep="\n"),padj=.5)
      axis(2,at=seq(0,1,length.out=10),label=length.ticks)
      if (nchar(pdfout)>0){
        dev.off()
      }
    }
    # minimum length threshold
    if (nchar(pdfout)>0){
      pdf(paste(pdfout,"_min.pdf",sep=""))
    }else{
      x11()
    }
    filled.contour2(seq(0,1,length.out=10),seq(0,1,length.out=10),t(anumreadm),
                    axes=FALSE,xlab="minimum quality",ylab="minimum length",color.palette=cpalette)
    title(paste("minimum length threshold vs base quality\n",fastqfile,"\n","total reads",tnr),cex.main=0.7)
    contour(seq(0,1,length.out=10),seq(0,1,length.out=10),t(anumreadm),axes=FALSE,add=T,levels=seq(0,nsamples,length.out=10),
            labels=paste(as.character(round(seq(0,tnr,length.out=10)*100/tnr)),"% - ",as.character(round(seq(0,tnr,length.out=10))),sep=""))  
    #axis(1,at=seq(0,1,length.out=10),label=qual.ticks)
    axis(1,at=seq(0,1,length.out=10),label=paste(qual.ticks,round(exp(qual.ticks/(-10)),digits=3),sep="\n"),padj=.5)
    axis(2,at=seq(0,1,length.out=10),label=length.ticks)
    if (nchar(pdfout)>0){
      dev.off()
    }
  }
  
  #writeFastq(trimmedReads,file="trimmed.fastq")
  return(list(anumreadM,anumreadm))
}

# allow color controur plot with levels overplotted
filled.contour2<-function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...) 
  {
    # modification by Ian Taylor of the filled.contour function
    # to remove the key and facilitate overplotting with contour()
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
      stop("increasing 'x' and 'y' values expected")
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2]) * par("csi") * 2.54
    par(las = las)
    mar <- mar.orig
    plot.new()
    par(mar=mar)
    print(paste(xlim,ylim))
    plot.window(xlim = xlim, ylim = ylim, log = "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
      stop("no proper 'z' matrix specified")
    if (!is.double(z)) 
      storage.mode(z) <- "double"
    if (getRversion()<3){
      .Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels), col = col))
    }else{
      .filled.contour(as.double(x), as.double(y), z, as.double(levels), col = col) # fix for R3
    }
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(x, side = 1)
        Axis(y, side = 2)
      }
    }
    else plot.axes
    if (frame.plot) 
      box()
    if (missing(plot.title)) 
      title(...)
    else plot.title
    invisible()
}

