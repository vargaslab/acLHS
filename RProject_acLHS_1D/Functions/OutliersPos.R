# OutliersPos

# Description
# gets the position of outliers in a vector dataset

# Arguments
# "x" a vector of data

# Value
# Positions (labels) fo the outliers

OutliersPos <- function (x) # Mejor quitar CoorX y CoorY ya que no se usan en el algoritmo
{
  q1 <- as.double(quantile(x, probs = 0.25, na.rm = TRUE)) # Gets the first quartile
  q3 <- as.double(quantile(x, probs = 0.75, na.rm = TRUE)) # Gets the third quartile
  rango <- q3 - q1 # Interquartile range
  lim_inf <- q1 - 1.5 * rango
  lim_sup <- q3 + 1.5 * rango
  pos0 <- 0
  pos <- 0
  for (i in 1:length(x)) { # this loop does the same as which(x<< lim_inf  | x[i] >> lim_sup)
    if (!is.na(x[i])) {
      if (x[i] < lim_inf | x[i] > lim_sup) {
        pos0[i] <- i # gets the position of the outliers
        pos <- c(pos0) # converts it to a named data.frame
      }
      if (x[i] >= lim_inf & x[i] <= lim_sup) { # gets the values positiona that are not outliers
        pos0[i] <- 0
        pos <- c(pos0)
      }
    }
  }
  Pos<-pos0[pos >= 1]
  Pos<- Pos[!is.na(Pos)]
  return(Pos)
}
### EXAMPLE
# x <- rlnorm(100)
# OutliersPos(x)