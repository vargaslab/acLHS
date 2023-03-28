CrossVariograma<-function (CoorX, CoorY, P1, P2, NInt, lags, Direccion, Tol, NomP1, 
          NomP2, NomP1P2) 
{

library(RandomFields)

  Datos <- as.data.frame(cbind(CoorX, CoorY, P1, P2))
  colnames(Datos) <- c("X", "Y", "Prop1", "Prop2")
  coordinates(Datos) <- ~X + Y
  ambas <- gstat(NULL, id = "Prop1", formula = Prop1 ~ 1, data = Datos)
  ambas <- gstat(ambas, id = "Prop2", formula = Prop2 ~ 1, 
                 data = Datos)
  v.cross <- variogram(ambas, cutoff = NInt * lags, width = lags, 
                       alpha = Direccion, tol.hor = Tol)
  valsemi <- v.cross$gamma
  minsemi <- min(valsemi)
  if (minsemi == 0) {
    cualzero <- 0
    cualeszero <- 0
    nsemi <- length(valsemi)
    for (i in 1:nsemi) {
      if (valsemi[i] == 0) {
        cualzero[i] <- i
        cualeszero <- c(cualzero)
      }
    }
    cualeszero <- na.omit(cualeszero)
    v.cross <- v.cross[-c(cualeszero), ]
  }
  nvar <- length(v.cross[, 1])
  k0 <- 0
  k1 <- 0
  k2 <- 0
  for (i in 1:nvar) {
    if (v.cross[i, 6] == "Prop1.Prop2") {
      k0 <- k0 + 1
    }
    if (v.cross[i, 6] == "Prop2") {
      k1 <- k1 + 1
    }
    if (v.cross[i, 6] == "Prop1") {
      k2 <- k2 + 1
    }
  }
  MaxX1 <- max((v.cross$dist[(k0 + k1 + 1):(k0 + k1 + k2)])) + 
    0.15 * (max((v.cross$dist[(k0 + k1 + 1):(k0 + k1 + k2)])))
  MaxY1 <- max((v.cross$gamma[(k0 + k1 + 1):(k0 + k1 + k2)]), 
               var(P1)) + 0.15 * (max((v.cross$gamma[(k0 + k1 + 1):(k0 + 
                                                                      k1 + k2)]), var(P1)))
  MaxX2 <- max((v.cross$dist[1:k0])) + 0.15 * (max((v.cross$dist[1:k0])))
  MaxY2 <- max((v.cross$gamma[1:k0])) + 0.15 * (max((v.cross$gamma[1:k0])))
  MaxX3 <- max((v.cross$dist[(k0 + 1):(k0 + k1)])) + 0.15 * 
    (max((v.cross$dist[(k0 + 1):(k0 + k1)])))
  MaxY3 <- max((v.cross$gamma[(k0 + 1):(k0 + k1)]), var(P2)) + 
    0.15 * (max((v.cross$gamma[(k0 + 1):(k0 + k1)]), var(P2)))
  orden <- matrix(c(1, 3, 2), ncol = 1, nrow = 3, byrow = T)
  div <- layout(orden, widths = 8, heights = c(3, 3, 3), TRUE)
  layout.show(div)
  par(mar = c(5, 4, 2, 1))
  plot(v.cross$dist[(k0 + k1 + 1):(k0 + k1 + k2)], v.cross$gamma[(k0 + 
                                                                    k1 + 1):(k0 + k1 + k2)], xlab = "Distance", ylab = "Semivariance", 
       xlim = c(0, MaxX1), ylim = c(0, MaxY1), pch = 21, col = "red", 
       bg = "blue", cex.lab = 1.3, cex.axis = 1.2)
  grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
  par(new=TRUE)
  plot(v.cross$dist[(k0 + k1 + 1):(k0 + k1 + k2)], v.cross$gamma[(k0 + 
                                                                    k1 + 1):(k0 + k1 + k2)], xlab = "Distance", ylab = "Semivariance", 
       xlim = c(0, MaxX1), ylim = c(0, MaxY1), pch = 21, col = "red", 
       bg = "blue", cex.lab = 1.3, cex.axis = 1.2)
  title(NomP1, cex.main = 1.5)
  par(mar = c(5, 4, 2, 1))
  plot(v.cross$dist[1:k0], v.cross$gamma[1:k0], xlab = "Distance", 
       ylab = "Semivariance", xlim = c(0, MaxX2), ylim = c(0, 
                                                           MaxY2), pch = 21, col = "red", bg = "blue", cex.lab = 1.3, 
       cex.axis = 1.2)
  grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
  par(new=TRUE)
  plot(v.cross$dist[1:k0], v.cross$gamma[1:k0], xlab = "Distance", 
       ylab = "Semivariance", xlim = c(0, MaxX2), ylim = c(0, 
                                                           MaxY2), pch = 21, col = "red", bg = "blue", cex.lab = 1.3, 
       cex.axis = 1.2)
  title(NomP1P2, cex.main = 1.5)
  par(mar = c(5, 4, 2, 1))
  plot(v.cross$dist[(k0 + 1):(k0 + k1)], v.cross$gamma[(k0 + 
                                                          1):(k0 + k1)], xlab = "Distance", ylab = "Semivariance", 
       xlim = c(0, MaxX3), ylim = c(0, MaxY3), pch = 21, col = "red", 
       bg = "blue", cex.lab = 1.3, cex.axis = 1.2)
  grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
  par(new=TRUE)
  plot(v.cross$dist[(k0 + 1):(k0 + k1)], v.cross$gamma[(k0 + 
                                                          1):(k0 + k1)], xlab = "Distance", ylab = "Semivariance", 
       xlim = c(0, MaxX3), ylim = c(0, MaxY3), pch = 21, col = "red", 
       bg = "blue", cex.lab = 1.3, cex.axis = 1.2)
  title(NomP2, cex.main = 1.5)
  Salida <- as.matrix(v.cross)
  return(Salida)
}