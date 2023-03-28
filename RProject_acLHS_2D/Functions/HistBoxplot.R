# HistBoxplot

# Description
# Boxplot over a histogram.

# Arguments
# "x" must be a numeric vector of data
# "nbins" is the number of bins of the histogram
# "mean"= mean of the x. "mean" could also be a vector of mean values when a mixture of distributions fit is to be plotted over the histogram. If it is not desired to display the mean do the same: median = NULL
# "median" = median of x. If it is not desired to display the median do the same: median = NULL
# "MedianLinePar" and "MeanLinePar" are 4-length vectos specifying line plotting parameterse: Name, lty, col,lwd
# "..." further arguments and graphical parameters passed to "hist2". breaks, col,mean, median, MedianLinePar, and MeanLinePar cannot be passed
# "col" is the color of the boxplot and of the histogram bars
# "plot" is the plot desired or only the breaks?

# Value
# the same as the base function hist.default

# Functions dependencies
# hist2 

HistBoxplot <- 
  function(x, nbins = 9, mean = NULL, median = NULL,
           MeanLinePar=list(name="Mean",col="Red",lty=1,lwd=2),
           MedianLinePar=list(name="Median",col="Blue",lty=3,lwd=2),
           col= "lightgray",plot=TRUE,...)
  {
    x <- x[is.finite(x)]
    Xmax <- max(x)
    Xmin <- min(x)
    Range <- Xmax-Xmin
    if (length(nbins) == 1) {
      BreaksL <- seq(Xmin-0.005*Range, Xmax+0.005*Range, length.out = nbins+1)
    } else BreaksL <- nbins
    if (plot) {
      orden<-matrix(c(2,1),ncol=1,nrow=2,byrow=T) # c(2,1) means that the second figure to be ploted will display in the 1st        frame of the plot window
      div<-layout(orden, widths= 7, heights=c(0.7, 3)) # "layout" Specify Complex Plot Arrangements
      par(mar=c(6,6,0,3)) 
      r <- hist2(x, breaks = BreaksL, col= col ,mean = mean, median = median, MedianLinePar=MedianLinePar,MeanLinePar=MeanLinePar,...)
      par(mar=c(0,6,3,3))
      plot(0,0,type="n", xlim=c(BreaksL[1],max(BreaksL)), ylim=c(0,1.54), xaxt='n',yaxt='n', xlab="",ylab="") 
      boxplot(x, range=1.5, horizontal= TRUE, col= col, pch= 21, axes= FALSE, add=TRUE ,outline = TRUE, varwidth=TRUE) 
      abline(v = mean, lty= MeanLinePar$lty,col = MeanLinePar$col, lwd=MeanLinePar$lwd )
      abline(v = median, lty = MedianLinePar$lty,col = MedianLinePar$col, lwd=MedianLinePar$lwd)
      invisible(r)
    } else return(BreaksL)
  }


#### EXAMPLE
# x<-rnorm(230)
# HistBoxplot(x=x,mean = mean(x), median = median(x), main ="", PercentFreq = FALSE )