Outliers<-function (X, Y, Z) 
{
  cuartil1 <- as.double(quantile(Z, probs = 0.25))
  cuartil3 <- as.double(quantile(Z, probs = 0.75))
  rango <- cuartil3 - cuartil1
  lim_inf <- cuartil1 - 1.5 * rango
  lim_sup <- cuartil3 + 1.5 * rango
  cual <- 0
  cuales <- 0
  for (i in 1:length(Z)) {
    if (Z[i] < lim_inf | Z[i] > lim_sup) {
      cual[i] <- i
      cuales <- c(cual)
    }
    if (Z[i] >= lim_inf & Z[i] <= lim_sup) {
      cual[i] <- 0
      cuales <- c(cual)
    }
  }
  Datos <- cbind(X, Y, Z)
  if (sum(cuales) != 0) {
    Datos <- Datos[-cuales, ]
  }
  
#  names(Datos)<-c("X","Y", "Z")
  return(Datos)
}