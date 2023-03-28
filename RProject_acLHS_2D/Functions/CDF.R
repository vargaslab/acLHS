# DESCRIPTION
# CDF plots the empirical and theoretical cumulative distribution function of a vector of data

# ARGUMENTS (INPUT)
# "x" a vector of values for which the histogram is desired.
# "distr" a character string of the distribution name. See the examples below and ?distributions.
# "para" A named list giving the parameters of the named distribution.
# "xlim" the x limits (x1, x2) of the plot. x1 > x2.
# "main" a main title for the plot, see also title. Default to "CDF".
# "xlab" a label for the x axis, defaults to a description of x.
# "ylab" a label for the y axis, defaults to "CDF".
# "col" a character string of the points color. See ?par
# "lcol" a character string of the line color of the model curve.  See "col" in ?par 
# "lty" a character string of the line type  of the model curve. See ?par
# "lwd" a character string of the line width of the model curve. See ?par
# "..." further arguments passed to histogram.default

# DETAILS
#  This function is based on the function plotdist from the package fitdistrplus

# VALUE
# a list with two elements:
# 1) a data.frame with two columns: one of observations and one with the CDF values
# 2) a vector with two elements of the "x" limits of the plot

CDF <- function (x, distr, para, xlim = NULL,
                 main = "CDF", xlab = NULL, ylab = NULL, col = "lightgray", pch = 20, cex = 1.8,
                 lcol = par("col"), lty = par("lty"), lwd = par("lwd"), ...)
{
  
  if (is.null(xlim))
    xlim <- c(min(x),max(x))
  
  xlab <- ifelse(!is.null(xlab),
                 xlab,
                 deparse(substitute(x)))
  
  ylab <- ifelse(!is.null(ylab),
                 ylab,
                 "CDF")
  
  s <- sort(x) # This sorting is needed for the CDF plot below. 
  obsp <- ecdf(s)(s) # "obsp" are the CDF values of "data"
  if (!is.character(distr)) {
    distname <- substring(as.character(match.call()$distr), 2) # gets the distribution name from the argument "distr" (by deleting the first letter and converting to character) in case is given as "distr=dnorm", for example
  } else distname <- distr
  if (!is.list(para)) 
    stop("'para' must be a named list")
  pdistname <- paste("p", distname, sep = "") # builds up the "p" probability distribution function name as a character string
  if (!exists(pdistname, mode = "function")) 
    stop(paste("The ", pdistname, " function must be defined"))
  
  plot.default(s, obsp, xlim = xlim,
               main = main, xlab = xlab, ylab = ylab, col = col, pch = pch, cex = cex,...) # plots the CDF as points with coordinates (s,obsp). Both "s" and "obsp" are given above
  sfin <- seq(xlim[1], xlim[2], by = (xlim[2]-xlim[1])/100)
  theopfin <- do.call(pdistname, c(list(q = sfin), as.list(para)))
  lines.default(sfin, theopfin, col = lcol, lty = lty, lwd = lwd) # plot of the theoretical CDF function over the empirical CDF
  
  r <- list(CDF = data.frame(s,obsp), xlim = xlim)
  invisible(r)
  
}

##### EXAMPLES
## EXAMPLE 1: NORMAL DISTRIBUTION
# set.seed(123)
# MEAN <- 2
# SD <- 0.5
# X <- rnorm(n = 300, mean = MEAN, sd = SD)
# PARA <- list(mean = MEAN, sd = SD)
# Results <- CDF(x = X, distr = "norm", para = PARA)
# CDF(x = X, distr = "norm", para = PARA, col = "orange")
# CDF(x = X, distr = "norm", para = PARA, col = "orange", main = "Test Main")
# CDF(x = X, distr = "norm", para = PARA, col = "orange", main = "Test Main", pch = "i")
# CDF(x = X, distr = "norm", para = PARA, col = "orange", main = "Test Main", pch = "i", lcol = "blue", lwd = 3)
# 
## EXAMPLE 2: EXPONENTIAL DISTRIBUTION
# set.seed(123)
# RATE <- 0.5
# X <- rexp(300, RATE)
# PARA <- list(rate = RATE)
# Results <- CDF(x = X, distr = "exp", para = PARA)