# QQplot

# Description
# QQplot plots the empirical and theoretical quantiles of a vector of data.

# Arguments
# "x" a vector of values for which the histogram is desired.
# "distr" a character string of the distribution name. See the examples below and ?distributions.
# "para" A named list giving the parameters of the named distribution.
# "xlim" the x limits (x1, x2) of the plot. x1 > x2.
# "main" a main title for the plot, see also title. Default to "QQ-plot".
# "xlab" a label for the x axis, defaults to a description of x.
# "ylab" a label for the y axis, defaults to "Cuantiles".
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
# a Quantile-Quantile plot

QQplot <- function (x, distr, para, xlim = NULL,
                    main = "Q-Q plot", xlab = NULL, ylab = NULL, col = "lightgray", pch = 20, cex = 1.8,
                    lcol = par("col"), lty = par("lty"), lwd = 2, ...)
{
  
  if (is.null(xlim))
    xlim <- c(min(x),max(x))
  
  xlab <- ifelse(!is.null(xlab),
                 xlab,
                 "Cuantiles teoricos")
  
  ylab <- ifelse(!is.null(ylab),
                 ylab,
                 "Cuantiles muestrales")
  
  s <- sort(x) # sort data increasingly 
  if (!is.character(distr)) {
    distname <- substring(as.character(match.call()$distr), 2) # gets the distribution name from the argument "distr" (by deleting the first letter and converting to character) in case is given as "distr=dnorm", for example
  } else distname <- distr
  if (!is.list(para)) 
    stop("'para' must be a named list")
  qdistname <- paste("q", distname, sep = "") # builds up the "p" probability distribution function name as a character string
  if (!exists(qdistname, mode = "function")) 
    stop(paste("The ", qdistname, " function must be defined"))
  
  n <- length(x)
  theoq <- do.call(qdistname, c(list(p = ppoints(n)), as.list(para)))  # "theoq" are the quantiles
  plot.default(theoq, s, xlim = xlim,
               main = main, xlab = xlab, ylab = ylab, col = col, pch = pch, cex = cex, ...) # plots the Q-Q plot
  abline(0, 1,col = lcol, lty = lty, lwd = lwd) # adds a line (with slope=1) of the form y=a+bx, where a=0, b=1 to the Q-Q plot
  
  r <- list(qq = data.frame(theoq,s), xlim = xlim)
  invisible(r)
  
}

##### EXAMPLES
## EXAMPLE 1: NORMAL DISTRIBUTION
# set.seed(123)
# MEAN <- 2
# SD <- 0.5
# X <- rnorm(n = 300, mean = MEAN, sd = SD)
# PARA <- list(mean = MEAN, sd = SD)
# Results <- QQplot(x = X, distr = "norm", para = PARA)
# QQplot(x = X, distr = "norm", para = PARA, col = "orange")
# QQplot(x = X, distr = "norm", para = PARA, col = "orange", main = "Test Main")
# QQplot(x = X, distr = "norm", para = PARA, col = "orange", main = "Test Main", pch = "i")
# QQplot(x = X, distr = "norm", para = PARA, col = "orange", main = "Test Main", pch = "i", lcol = "blue", lwd = 3)
# # 
# ## EXAMPLE 2: EXPONENTIAL DISTRIBUTION
# set.seed(123)
# RATE <- 0.5
# X <- rexp(300, RATE)
# PARA <- list(rate = RATE)
# Results <- QQplot(x = X, distr = "exp", para = PARA)