ModelVariogram<-function (CoorX, CoorY, P1, P2, NInt, lags, Direccion, Tol, Modelo, 
          Sill1, Sill2, Sill3, Nugget1, Nugget2, Nugget3, Alcance, 
          NomP1, NomP2, NomP1P2) 
{
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
  Modelo1 <- Modelo
  if (Modelo == 1) {
    Modelo <- "Exp"
  }
  if (Modelo == 2) {
    Modelo <- "Sph"
  }
  if (Modelo == 3) {
    Modelo <- "Gau"
  }
  if (Modelo1 == 1) {
    Modelo1 <- "exp"
  }
  if (Modelo1 == 2) {
    Modelo1 <- "sph"
  }
  if (Modelo1 == 3) {
    Modelo1 <- "gau"
  }
  modelo <- gstat(ambas, id = "Prop1", model = vgm(Sill1, Modelo, 
                                                   Alcance, Nugget1), nmax = Inf, nmin = 0, maxdist = NInt * 
                    lags)
  modelo <- gstat(modelo, id = "Prop2", model = vgm(Sill2, 
                                                    Modelo, Alcance, Nugget2), nmax = Inf, nmin = 0, maxdist = NInt * 
                    lags)
  modelo <- gstat(modelo, id = c("Prop1", "Prop2"), model = vgm(Sill3, 
                                                                Modelo, Alcance, Nugget3), nmax = Inf, nmin = 0, maxdist = NInt * 
                    lags)
  MaxX1 <- max((v.cross$dist[(k0 + k1 + 1):(k0 + k1 + k2)]), 
               Alcance) + 0.22 * (max((v.cross$dist[(k0 + k1 + 1):(k0 + 
                                                                     k1 + k2)]), Alcance))
  MaxY1 <- max((v.cross$gamma[(k0 + k1 + 1):(k0 + k1 + k2)]), 
               var(P1), Sill1 + Nugget1) + 0.22 * (max((v.cross$gamma[(k0 + 
                                                                         k1 + 1):(k0 + k1 + k2)]), var(P1), Sill1 + Nugget1))
  MaxX2 <- max((v.cross$dist[1:k0]), Alcance) + 0.22 * (max((v.cross$dist[1:k0]), 
                                                            Alcance))
  MaxY2 <- max((v.cross$gamma[1:k0]), Sill3 + Nugget3) + 0.22 * 
    (max((v.cross$gamma[1:k0]), Sill3 + Nugget3))
  MaxX3 <- max((v.cross$dist[(k0 + 1):(k0 + k1)]), Alcance) + 
    0.22 * (max((v.cross$dist[(k0 + 1):(k0 + k1)]), Alcance))
  MaxY3 <- max((v.cross$gamma[(k0 + 1):(k0 + k1)]), var(P2), 
               Sill2 + Nugget2) + 0.22 * (max((v.cross$gamma[(k0 + 1):(k0 + 
                                                                         k1)]), var(P2), Sill2 + Nugget2))
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
  lines.variomodel(cov.model = Modelo1, cov.pars = c(Sill1, 
                                                     Alcance), nug = Nugget1, max.dist = NInt * lags, col = "black", 
                   lwd = 2)
  abline(h = Sill1 + Nugget1, lty = 3, col = "Blue", lwd = 2)
  abline(h = var(P1), lty = 1, col = "Red", lwd = 2)
  abline(v = Alcance, lty = 5, col = "green", lwd = 2)
  legend(0.7 * MaxX1, 0.5 * MaxY1, c("Sill+Nugget", "Range", 
                                     "Variance"), col = c("Blue", "green", "Red"), lty = c(3, 
                                                                                           5, 1), bty = "n")
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
  lines.variomodel(cov.model = Modelo1, cov.pars = c(Sill3, 
                                                     Alcance), nug = Nugget3, max.dist = NInt * lags, col = "black", 
                   lwd = 2)
  abline(h = Sill3 + Nugget3, lty = 3, col = "Blue", lwd = 2)
  abline(v = Alcance, lty = 5, col = "green", lwd = 2)
  legend(0.7 * MaxX2, 0.5 * MaxY2, c("Sill+Nugget", "Range"), 
         col = c("Blue", "green"), lty = c(3, 5), bty = "n")
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
  lines.variomodel(cov.model = Modelo1, cov.pars = c(Sill2, 
                                                     Alcance), nug = Nugget2, max.dist = NInt * lags, col = "black", 
                   lwd = 2)
  abline(h = Sill2 + Nugget2, lty = 3, col = "Blue", lwd = 2)
  abline(h = var(P2), lty = 1, col = "Red", lwd = 2)
  abline(v = Alcance, lty = 5, col = "green", lwd = 2)
  legend(0.7 * MaxX3, 0.5 * MaxY3, c("Sill+Nugget", "Range", 
                                     "Variance"), col = c("Blue", "green", "Red"), lty = c(3, 
                                                                                           5, 1), bty = "n")
  title(NomP2, cex.main = 1.5)
}