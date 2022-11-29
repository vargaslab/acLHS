OutliersCount<-function (CoorX, CoorY, propiedad) 
{
  cuartil1 <- as.real(quantile(propiedad, probs = 0.25))
  cuartil3 <- as.real(quantile(propiedad, probs = 0.75))
  rango <- cuartil3 - cuartil1
  lim_inf <- cuartil1 - 1.5 * rango
  lim_sup <- cuartil3 + 1.5 * rango
  cual <- 0
  cuales <- 0
  for (i in 1:length(propiedad)) {
    if (propiedad[i] < lim_inf | propiedad[i] > lim_sup) {
      cual[i] <- i
      cuales <- c(cual)
    }
    if (propiedad[i] >= lim_inf & propiedad[i] <= lim_sup) {
      cual[i] <- 0
      cuales <- c(cual)
    }
  }
  count <- length(cuales[cuales >= 1])
  return(count)
}