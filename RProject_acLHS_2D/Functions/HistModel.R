# HistModel

# Description
# HistModel plots theoretical probability density models over a histogram.

# Arguments
# "x" a vector of values for which the histogram is desired.
# "distr" a character string of the distribution name. See the examples below and ?distributions.
# "para" A named list giving the parameters of the named distribution.
# "freq" Logical. see the histogram.default function documentation.
# "colCurve" a character string of the line color of the model curve.
# "ltyCurve" a character string of the line type  of the model curve.
# "lwdCurve" a character string of the line width of the model curve.
# "..." further arguments passed to histogram.default

# Details
#  This function is based on the function plotdist from the package fitdistrplus and in the function hist.default

# Value
# The same as in the function hist.default


HistModel <- function(x, distr, para, freq = TRUE,
                      colCurve =  "black", ltyCurve = "solid", lwdCurve = 2,
                      ...)
{
  library(poweRlaw)
  
  # Setting the histogram parameters
  Breaks <- list(...) $breaks
  breaks <- ifelse(!is.null(Breaks),Breaks, "Sturges") 
  Hist <- hist.default(x = x, plot = FALSE, breaks = breaks)
  breaks   <- Hist$breaks
  equidist <- Hist$equidist
  counts   <- Hist$counts
  density  <- Hist$density
  
  if(!is.null(distr)){ # if FALSE, then go to the Power-law case
    xhist <- seq(min(breaks), max(breaks), length = 1000) # "xhist" are the values for which the probability density will be computed
    if (!is.list(para)) 
      stop("'para' must be a named list")    
    if (distr != "pl") {
      if (!is.character(distr)){ 
        distname <- substring(as.character(match.call()$distr), 2) # gets the distribution name from the argument "distr" (by deleting the first letter and converting to character) in case is given as "distr=dnorm", for example
      } else distname <- distr
      ddistname <- paste("d", distname, sep = "") # builds up the "d" density function name as a character string
      if (!exists(ddistname, mode = "function")) # checks if the density function is already defined in R
        stop(paste("The ", ddistname, " function must be defined"))
      densfun <- get(ddistname, mode = "function") # converts "densfun"  to a function
      nm <- names(para) # gets the names (as characters) given by the user
      f <- formals(densfun) # Access to and Manipulation of the Formal Arguments. "f" is a named list with the names of the arguments of the function (densfun in this case). run formals(dnorm) for instance
      args <- names(f) # gets the names of the list "f" above. "args" is a character vector
      m <- match(nm, args) # m stores the elements which are the same in both "nm" and "args"
      if (any(is.na(m))) 
        stop(paste("'para' specifies names which are not arguments to ", ddistname))
      yhist <- do.call(ddistname, c(list(x = xhist), as.list(para))) # "yhist" are the theoretical probabilities densities
    } else { # Powew-law case:
      library(poweRlaw)
      m <- conpl$new() # Set the data as continous power-law object
      m$setXmin(para$xmin)
      m$setPars(para$alpha)
      yhist <- dist_pdf(m = m, q = xhist)
    }
  }
  
  if(equidist){
    BarsHeight <- ifelse(!freq, density, counts)
    if (freq)
      yhist <- yhist* counts[1]/density[1]# upscale to fit the histogram when freq = TRUE
  }
  ymax <- ifelse(is.finite(max(yhist)),
                 max(max(BarsHeight), max(yhist)),
                 max(BarsHeight))
  hist.default(x = x, freq = freq, ylim = c(0, ymax),...)
box()
  lines.default(xhist, yhist, col = colCurve, lty = ltyCurve, lwd = lwdCurve) # plot of the density function over the histogram
  
  invisible(Hist)
}

# ##### EXAMPLES
## EXAMPLE 1: NORMAL DISTRIBUTION
# MEAN <- 2
# SD <- 0.5
# X <- rnorm(n = 300, mean = MEAN, sd = SD)
# PARA <- list(mean = MEAN, sd = SD)
# HistModel(x = X, distr = "norm", para = PARA)
# HistModel(x = X, distr = "norm", para = PARA, breaks = 30)
# HistModel(x = X, distr = "norm", para = PARA, breaks = 30, freq = FALSE)
# 
# ## EXAMPLE 2: POWER-LAW DISTRIBUTION
# library(poweRlaw)
# data(native_american)
# X <- native_american$Cas
# XMIN <- 27
# ALPHA <- 2.243601
# PARA <- list( xmin = XMIN, alpha = ALPHA)
# HistModel(x = X, distr = "pl", para = PARA,
#           col = "green", main = "Histogram of Power-law data",
#           colCurve ="blue")
# 
