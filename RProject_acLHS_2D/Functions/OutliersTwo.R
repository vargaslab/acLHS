OutliersTwo<-function (CoorX, CoorY, P1, P2, XLAB, YLAB, Quitar) 
{
  require(reshape)
  DatosN <- as.data.frame(cbind(CoorX, CoorY, P1, P2))
  DatosOrig <- DatosN
  loess.P1 <- lm(P1 ~ P2)
  Estim <- predict(loess.P1)
  Error <- c(P1 - Estim)
  Error <- round(Error, digits = 3)
  cuartil1 <- as.real(quantile(Error, probs = 0.2))
  cuartil3 <- as.real(quantile(Error, probs = 0.8))
  rango <- cuartil3 - cuartil1
  lim_inf <- cuartil1 - 3 * rango
  lim_sup <- cuartil3 + 3 * rango
  cual <- 0
  cuales <- 0
  for (i in 1:length(Error)) {
    if (Error[i] < lim_inf | Error[i] > lim_sup) {
      cual[i] <- i
      cuales <- c(cual)
    }
    if (Error[i] >= lim_inf & Error[i] <= lim_sup) {
      cual[i] <- 0
      cuales <- c(cual)
    }
  }
  if (sum(cuales) > 0) {
    control <- cuales[cuales >= 1]
    Errores <- Error[cuales]
    ErrorAbs <- abs(Errores)
    ErrorU <- as.data.frame(cbind(control, ErrorAbs))
    ErrorOrd <- sort_df(ErrorU, vars = "ErrorAbs")
    orden <- matrix(1, ncol = 1, nrow = 1, byrow = T)
    div <- layout(orden, heights = 10, widths = 10, TRUE)
    layout.show(div)
    par(mar = c(5, 5, 5, 5))
    plot(DatosN[, 3], DatosN[, 4], xlab = XLAB, ylab = YLAB, 
         pch = 19, col = "black", cex.lab = 0.8, cex.axis = 0.8)
    NSup <- length(ErrorOrd[, 1])
    NInf <- length(ErrorOrd[, 1]) - Quitar + 1
    CualQuitar <- ErrorOrd[c(NSup:NInf), 1]
    DatosN <- as.matrix(DatosN[-CualQuitar, ])
    DatosN1 <- DatosOrig[CualQuitar, ]
    matplot(DatosN1[, 3], DatosN1[, 4], add = T, pch = 13, 
            col = "red")
  }
  DatosN <- as.matrix(DatosN)
  return(DatosN)
}
