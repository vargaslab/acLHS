# PPplot

# Description
# PPplot plots the empirical and theoretical cumulative distribution function of a vector of data

# Arguments
# "x" a vector of values for which the histogram is desired.
# "distr" a character string of the distribution name. See the examples below and ?distributions.
# "para" A named list giving the parameters of the named distribution.
# "main" a main title for the plot, see also title. Default to "P-P plot".
# "xlab" a label for the x axis, defaults to a description of x.
# "ylab" a label for the y axis, defaults to "Probabilidades emp?ricas".
# "col" a character string of the points color. See ?par
# "pch" plotting "character", i.e., symbol to use
# "cex" the size of the points
# "lcol" a character string of the line color of the model curve.  See "col" in ?par 
# "lty" a character string of the line type  of the model curve. See ?par
# "lwd" a character string of the line width of the model curve. See ?par
# "..." further arguments passed to histogram.default

# Details
#  This function is based on the function plotdist from the package fitdistrplus

# Value
# a Percentile-Percentile plot

PPplot <- function (x, distr, para,
                    main = "P-P plot", xlab = NULL, ylab = NULL, col = "lightgray", pch = 20, cex = 1.8,
                    lcol = par("col"), lty = par("lty"), lwd = 2, ...)
{
  
  xlab <- ifelse(!is.null(xlab),
                 xlab,
                 "Probabilidades teoricas")
  
  ylab <- ifelse(!is.null(ylab),
                 ylab,
                 "Probabilidades empiricas")
  
  s <- sort(x) # sort data increasingly
  obsp <- ecdf(s)(s) # "obsp" are the CDF values of "data"
  if (!is.character(distr)) {
    distname <- substring(as.character(match.call()$distr), 2) # gets the distribution name from the argument "distr" (by deleting the first letter and converting to character) in case is given as "distr=dnorm", for example
  } else distname <- distr
  if (!is.list(para)) 
    stop("'para' must be a named list")
  pdistname <- paste("p", distname, sep = "") # builds up the "p" probability distribution function name as a character string
  if (!exists(pdistname, mode = "function")) 
    stop(paste("The ", pdistname, " function must be defined"))
  
  theop <- do.call(pdistname, c(list(q = s), as.list(para))) # "theop" are the theoretical probabilities
  
  plot.default(theop, obsp,
               main = main, xlab = xlab, ylab = ylab, col = col, pch = pch, cex = cex, ...)
  
  abline(0, 1,col = lcol, lty = lty, lwd = lwd) # adds a line (with slope=1) of the form y=a+bx, where a=0, b=1
  
  r <- list(pp = data.frame(theop,obsp))
  invisible(r)
  
}

##### EXAMPLES
# # EXAMPLE 1: NORMAL DISTRIBUTION
# set.seed(123)
# MEAN <- 2
# SD <- 0.5
# X <- rnorm(n = 300, mean = MEAN, sd = SD)
# PARA <- list(mean = MEAN, sd = SD)
# Results <- PPplot(x = X, distr = "norm", para = PARA)
# PPplot(x = X, distr = "norm", para = PARA, col = "orange")
# PPplot(x = X, distr = "norm", para = PARA, col = "orange", main = "Test Main")
# PPplot(x = X, distr = "norm", para = PARA, col = "orange", main = "Test Main", pch = "i")
# PPplot(x = X, distr = "norm", para = PARA, col = "orange", main = "Test Main", pch = "i", lcol = "blue", lwd = 3)
# # 
# ## EXAMPLE 2: EXPONENTIAL DISTRIBUTION
# set.seed(123)
# RATE <- 0.5
# X <- rexp(300, RATE)
# PARA <- list(rate = RATE)
# Results <- PPplot(x = X, distr = "exp", para = PARA)