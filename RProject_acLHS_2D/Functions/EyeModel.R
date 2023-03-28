EyeModel<- function (CoorX, CoorY, Variable, Direccion, Tol, NIntervalos, 
          Lags, Npares, Modelo, Nugget, SillYNugget, Alcance, MainTitle) 
{
  Sill <- SillYNugget - Nugget
  DatosGeo <- as.geodata(cbind(CoorX, CoorY, Variable), header = TRUE, 
                         coords.col = 1:2, data.col = 3)
  Nuevo <- variog(DatosGeo, breaks = c(seq(Lags/2, Lags * (NIntervalos+1), 
                                           Lags)), trend = "cte", lambda = 1, estimator.type = "classical", 
                  nugget.tolerance = 0, direction = Direccion, tolerance = Tol, 
                  unit.angle = "degrees", pairs.min = Npares)
  Modelo1 <- Modelo
  if (Modelo == 1) {
    Modelo <- "exp"
  }
  if (Modelo == 2) {
    Modelo <- "sph"
  }
  if (Modelo == 3) {
    Modelo <- "gau"
  }
  if (Modelo1 == 1) {
    ModeloA <- "exponential"
  }
  if (Modelo1 == 2) {
    ModeloA <- "spherical"
  }
  if (Modelo1 == 3) {
    ModeloA <- "gaussian"
  }
  VMod <- (SillYNugget) - cov.spatial(Nuevo$u, cov.model = ModeloA, 
                                      cov.pars = c(Sill, Alcance))
  Error2 <- (Nuevo$v - VMod)^2
  SumaError <- sum(Error2)
  MaxX <- max(Nuevo$u, Alcance) + 0.22 * (max(Nuevo$u, Alcance))
  MaxY <- max(Nuevo$v, Nuevo$var.mark, Sill + Nugget) + 0.22 * 
    (max(Nuevo$v, Nuevo$var.mark, Sill + Nugget))
  orden <- matrix(c(1, 2), ncol = 1, nrow = 2, byrow = T)
  div <- layout(orden, widths = 9, heights = c(3.5, 0.8))
  layout.show(div)
  par(mar = c(5, 5, 2, 2))
  plot(Nuevo$u, Nuevo$v, pch = 21, col = "black", bg = "blue", 
       xlim = c(0, MaxX), ylim = c(0, MaxY), xlab = "Time [days]", 
       ylab = "Semivariance", cex.lab = 1.5, cex.axis = 1.2, 
       main = MainTitle)
  grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
  par(new=TRUE)
  plot(Nuevo$u, Nuevo$v, pch = 21, col = "black", bg = "blue", 
       xlim = c(0, MaxX), ylim = c(0, MaxY), xlab = "Time [days]", 
       ylab = "Semivariance", cex.lab = 1.5, cex.axis = 1.2, 
       main = MainTitle)
  lines.variomodel(cov.model = Modelo, cov.pars = c(Sill, Alcance), 
                   nug = Nugget, col = "Black", lwd = 2, max.dist = Lags * 
                     NIntervalos)
  abline(h = Nuevo$var.mark, lty = 5, lwd = 2, col = "Red")
  abline(h = Sill + Nugget, lty = 4, lwd = 2, col = "darkviolet")
  abline(v = Alcance, lty = 3, lwd = 2, col = "Green")
  legend("topleft", c(Modelo, "Variance", "Sill+Nugget", "Range"), # 0, MaxY,
         col = c("black", "red", "darkviolet", "green"), lty = c(1, 
                                                                 5, 4, 3), lwd = 2, bty = "n")
  ModEle <- cbind(Nugget, Sill + Nugget, Alcance)
  rownames(ModEle) <- c(Modelo)
  colnames(ModEle) <- c("Nugget", "Meseta+Nugget", "Alcance")
  par(mar = c(0.5, 0.5, 0, 0.5))
  plot(0, 0, type = "n", xlim = c(0, 60), ylim = c(0, 7), xaxt = "n", 
       yaxt = "n", xlab = "", ylab = "")
  text(4, 5, labels = "Model", cex = 1.1)
  text(17, 5, labels = "Nugget", cex = 1.1)
  text(32, 5, labels = "Sill+Nugget", cex = 1.1)
  text(47, 5, labels = "Range", cex = 1.1)
  text(58, 5, labels = "MSE", cex = 1.1)
  text(4, 3, labels = Modelo, cex = 1.1)
  text(17, 3, labels = sprintf("%.4f", Nugget), cex = 1.1)
  text(32, 3, labels = sprintf("%.4f", Sill + Nugget), cex = 1.1)
  text(47, 3, labels = sprintf("%.4f", Alcance), cex = 1.1)
  text(58, 3, labels = sprintf("%.4f", SumaError), cex = 1.1)
  return(ModEle)
}