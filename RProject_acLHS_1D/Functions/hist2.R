# hist2

# Description
# plot a histogram

# Arguments
# See the base function hist.default for the parameters in {{{----}}}
# mean = mean of the x. "mean" could also be a vector of mean values when a mixture of distributions fit is to be plotted over the histogram. If it is not desired to display the mean do the same: median = NULL
# median = median of x. If it is not desired to display the median do the same: median = NULL
# MedianLinePar, MeanLinePar and FitLinePar are 4-length vectos specifying line plotting parameterse: Name, lty, col,lwd
# AbsFreq = logical; if TRUE, Absolute frequencies are ploted over the histogram bins. Defaults to TRUE.
# PercentFreq = logical; if TRUE, percent frequencies are ploted inside the histogram bins. Defaults to TRUE.
# "legendPosition" a list as list(x="bottomright",y = NULL) or NULL (default) to automatic positioning
# fit = a 2-columns vector of values of a x-axis in the first column and fitting values in the second. This will be plotted as a black line over the histogram. It can also be a string character as "lnorm" or "exp" to mean that it is a known distribution function
# At = upscaling value to match the fitting curve when "freq = TRUE"
# "duplicate" is the histogram to be duplicated?
# "period" if data is periodic, what is the period. 2pi for circular data, for intance

# Value
# the same as the base function hist.default


hist2 <- function(
  #--------------{{{{{{-------------------------
  # see the documentation of the base function hist.default for the next parameters:
  x, breaks =  "Sturges", right = FALSE,
  freq = NULL,
  col = "lightgray", main = "Histogram", 
  ylim = NULL, xlab = "X", ylab = "",  
  #--------------}}}}}}-------------------------
  mean = NULL, median = NULL, fitMean= NULL,
  MeanLinePar = list(name="Media",col="Red",lty=1,lwd=2),
  MedianLinePar = list(name="Mediana",col="Blue",lty=3,lwd=2),
  FitLinePar = list(name="Ajuste",col="Black",lty=1,lwd=2),
  fitMeanLinePar = list(name="Media del Ajuste",col="Black",lty=3,lwd=2),
  AbsFreq = NULL, PercentFreq = NULL, 
  legendPosition = NULL, fit = NULL, At = NULL, duplicate = FALSE, period = NULL,   
  ...)
{
  n <- length(x <- x[is.finite(x)]) # "is.finite" returns a vector of the same length as x whose jth element of which is TRUE if x[j] is finite (i.e., it is not one of the values NA, NaN, Inf or -Inf) and FALSE otherwise.

  #res <- hist.default(x = x,right = right, breaks = breaks, plot = FALSE,...)
  res <- hist.default(x=x,right = right,breaks = breaks,plot = FALSE)
  Breaks <- res$breaks
  nB <- length(Breaks)
  counts <- res$counts
  dens <- counts/sum(counts)*100 #res$density
  mids <- res$mids

  if (duplicate == TRUE){
    if (!is.null(period) && Breaks[nB]!=period) { # this "if" is for ranges that do not span almost al the "period"
      fit<-rbind(fit,cbind(fit[-1,1]+period,fit[-1,2]))
    }
    if (!is.null(period) && Breaks[nB]==period) {
      counts <-c(counts,counts)
      Breaks <-c(Breaks,Breaks[-1]+period)
      x <- c(x, x+period)
      if (!is.null(fit)) fit<-rbind(fit,cbind(fit[-1,1]+period,fit[-1,2]))
    }
    else if (!is.null(period)) print("Maximum value of 'breaks' does not equal the period of data")
  }

  if (AbsFreq == TRUE) {
    ylim<-c(0,max(counts)*1.15) # From Javier?s function "BasicS". maxY+maxY*0.15 = maxY *(1+0.15)=  maxY*1.15
    if (is.null(freq) & !is.null(fit)) {
      Ymax<-max(c(counts,At*fit[,2]))
      ylim<-c(0,1.15*Ymax)
    }
    hist.default(x = x, breaks = Breaks, right = right, col = col, xlab = xlab, ylab = ylab, ylim = ylim, main = main, ...)
  }
  if (PercentFreq == TRUE) {
    ylim<-c(0,max(dens)*1.15) # From Javier?s function "BasicS". maxY+maxY*0.15 = maxY *(1+0.15)=  maxY*1.15
    if (is.null(freq) & !is.null(fit)) {
      Ymax<-max(c(dens,At*fit[,2]))
      ylim<-c(0,1.15*Ymax)
    }
    h = hist.default(x = x,plot = FALSE, breaks = Breaks) # , right = right, col = col, xlab = xlab, ylab = ylab, ylim = ylim, main = main
    h$density = h$counts/sum(h$counts)*100
    plot(h,freq=FALSE, col = col, xlab = xlab, ylab = ylab, ylim = ylim, main = main,...) # right = right,
  }    
  
  if (is.null(legendPosition)) {
    Skewness<-(sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
    if (Skewness > 0) {
      legendPosition <- list(x = "topright", y = NULL)
    } else {
      legendPosition <- list(x = "topleft", y = NULL)
    }
  } else if (is.character(legendPosition)) legendPosition <- list(x = legendPosition, y = NULL)
  abline(v = mean, lty= as.integer(MeanLinePar$lty),col = MeanLinePar$col, lwd=as.integer(MeanLinePar$lwd) )
  abline(v = median, lty = as.integer(MedianLinePar$lty),col = MedianLinePar$col, lwd=as.integer(MedianLinePar$lwd))
  if (is.null(fit)) {
    if (is.null(median) & !is.null(mean))
      legend(x = legendPosition$x, y = legendPosition$y,MeanLinePar$name, lty=as.integer(MeanLinePar$lty),col= MeanLinePar$col, bty="n",lwd=as.integer(MeanLinePar$lwd))
    else  if (!is.null(median) & !is.null(mean))
      legend(x = legendPosition$x, y = legendPosition$y,c(MeanLinePar$name,MedianLinePar$name), lty=as.integer(c(MeanLinePar$lty,MedianLinePar$lty)),col= c(MeanLinePar$col,MedianLinePar$col), bty="n",lwd=c(MeanLinePar$lwd,MedianLinePar$lwd))
    if (!is.null(median) & is.null(mean))
      legend(x = legendPosition$x, y = legendPosition$y,MedianLinePar$name, lty=as.integer(MedianLinePar$lty),col= MedianLinePar$col, bty="n",lwd=MedianLinePar$lwd)
  }
  if (is.null(freq) && !is.null(fit)) {
    segments(fit[-nrow(fit),1],At*fit[-nrow(fit),2],fit[-1,1],At*fit[-1,2],lwd=3)
    legend(x = legendPosition$x, y = legendPosition$y,c(MeanLinePar$name,FitLinePar$name), lty=as.integer(c(MeanLinePar$lty,FitLinePar$lty)),col= c(MeanLinePar$col,FitLinePar$col), bty="n",lwd=c(MeanLinePar$lwd,FitLinePar$lwd))
  }
  if (!is.null(fitMean)) {
    abline(v = fitMean, lty= as.integer(fitMeanLinePar$lty),col = fitMeanLinePar$col, lwd=fitMeanLinePar$lwd)
    legend(x=legendPosition$x, y = legendPosition$y,fitMeanLinePar$name, lty=as.integer(fitMeanLinePar$lty),col= fitMeanLinePar$col, bty="n",lwd=fitMeanLinePar$lwd)
  }
  if (AbsFreq) 
    text(mids,counts,labels = counts,pos =3)
  if (PercentFreq) 
    text(mids,dens,labels =round(dens*100/sum(dens)),pos =3, col= "black")
  if(!is.null(period)) axis(3, at=c(0,180,360))
  box()
  abline(h = 0)
  if (!is.null(period) && Breaks[nB]==period) {
    axis(2, pos=period)
    abline(v=period)
  }

  invisible(res)
}

############# EXAMPLE
# x<-rnorm(350,4,0.7) # 350 data with mean 4 and standard deviation 0.7
# hist2(x,       median = 3.5, )
# hist2(x,mean= 4,             )
# x <- runif(100)
# hist2(x,mean= 0.5, duplicate = TRUE, period = 1, breaks = seq(0,1,0.1))
# hist2(x,mean= 0.5, duplicate = TRUE, period = 1, xlab = "asfdas")
# 
# library(circular)
### DATA IN RADIANS
# set.seed(123)
# x <- rvonmises(n=1000, mu=circular( 0,units="degrees", modulo="pi", template="geographics"), kappa=10, control.circular=list(units="degrees", modulo="pi", template="geographics"))
# x <- rmixedvonmises(n=1000, mu1=circular(pi/4,modulo="2pi"), mu2=circular(3*pi/2,modulo="2pi"), kappa1=4, kappa2=4, prop=0.3, control.circular = list(modulo="2pi"))
# x <- as.numeric(x)
# hist2(x, duplicate = TRUE, period = 2*pi, breaks = seq(0, 2*pi, length.out = 24),
#       mean = 1, median = 1.5, AbsFreq = FALSE, PercentFreq = FALSE)
# 
# x <- x/2
# hist2(x, duplicate = TRUE, period = pi, breaks = seq(0, pi, length.out = 12),
#       mean = 1, median = 1.5)
# 
# x <- x/2
# hist2(x, duplicate = TRUE, period = pi/2, breaks = seq(0, pi/2, length.out = 12),
#       mean = 1, median = 1.5)