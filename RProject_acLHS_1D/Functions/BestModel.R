BestModel<- function (CoorX, CoorY, Variable, Direccion, Tol, NIntervalos, 
          Lags, Npares, MainTitle) 
{
  DatosGeo <- as.geodata(cbind(CoorX, CoorY, Variable), header = TRUE, 
                         coords.col = 1:2, data.col = 3)
  Nuevo <- variog(DatosGeo, breaks = c(seq(Lags/2, Lags * (NIntervalos+1), 
                                           Lags)), trend = "cte", lambda = 1, estimator.type = "classical", 
                  nugget.tolerance = 0, direction = Direccion, tolerance = Tol, 
                  unit.angle = "degrees", pairs.min = Npares)
  max.pepita <- max(Nuevo$v)
  for (i in 1:length(Nuevo$v)) {
    if (Nuevo$v[i] == max.pepita) {
      rango <- Nuevo$u[i]
      rango <- trunc(rango)
    }
  }
  modelo_Exp00 <- variofit(Nuevo, ini = c(Nuevo$var.mark, rango), 
                           cov.model = "exponential", fix.nug = FALSE)
  modelo_Sph00 <- variofit(Nuevo, ini = c(Nuevo$var.mark, rango), 
                           cov.model = "spherical", fix.nug = FALSE)
  modelo_Gau00 <- variofit(Nuevo, ini = c(Nuevo$var.mark, rango), 
                           cov.model = "gaussian", fix.nug = FALSE)
  modelo_Exp0 <- variofit(Nuevo, ini = c(modelo_Exp00$cov.pars[1], 
                                         modelo_Exp00$cov.pars[2]), cov.model = "exponential", 
                          fix.nug = FALSE)
  modelo_Sph0 <- variofit(Nuevo, ini = c(modelo_Sph00$cov.pars[1], 
                                         modelo_Sph00$cov.pars[2]), cov.model = "spherical", fix.nug = FALSE)
  modelo_Gau0 <- variofit(Nuevo, ini = c(modelo_Gau00$cov.pars[1], 
                                         modelo_Gau00$cov.pars[2]), cov.model = "gaussian", fix.nug = FALSE)
  modelo_Exp <- variofit(Nuevo, ini = c(modelo_Exp0$cov.pars[1], 
                                        modelo_Exp0$cov.pars[2]), cov.model = "exponential", 
                         fix.nug = FALSE)
  modelo_Sph <- variofit(Nuevo, ini = c(modelo_Sph0$cov.pars[1], 
                                        modelo_Sph0$cov.pars[2]), cov.model = "spherical", fix.nug = FALSE)
  modelo_Gau <- variofit(Nuevo, ini = c(modelo_Gau0$cov.pars[1], 
                                        modelo_Gau0$cov.pars[2]), cov.model = "gaussian", fix.nug = FALSE)
  VMod_Exp <- (modelo_Exp$cov.pars[1] + modelo_Exp$nugget) - 
    cov.spatial(Nuevo$u, cov.model = modelo_Exp$cov.model, 
                cov.pars = c(modelo_Exp$cov.pars[1], modelo_Exp$cov.pars[2]))
  Error2_Exp <- (Nuevo$v - VMod_Exp)^2
  SumaError_Exp <- sum(Error2_Exp)
  VMod_Sph <- (modelo_Sph$cov.pars[1] + modelo_Sph$nugget) - 
    cov.spatial(Nuevo$u, cov.model = modelo_Sph$cov.model, 
                cov.pars = c(modelo_Sph$cov.pars[1], modelo_Sph$cov.pars[2]))
  Error2_Sph <- (Nuevo$v - VMod_Sph)^2
  SumaError_Sph <- sum(Error2_Sph)
  VMod_Gau <- (modelo_Sph$cov.pars[1] + modelo_Gau$nugget) - 
    cov.spatial(Nuevo$u, cov.model = modelo_Gau$cov.model, 
                cov.pars = c(modelo_Gau$cov.pars[1], modelo_Gau$cov.pars[2]))
  Error2_Gau <- (Nuevo$v - VMod_Gau)^2
  SumaError_Gau <- sum(Error2_Gau)
  SCerror <- min(c(SumaError_Exp, SumaError_Sph, SumaError_Gau))
  if (SCerror == SumaError_Exp) {
    modelo <- modelo_Exp
    modelog <- "Exp"
  }
  if (SCerror == SumaError_Sph) {
    modelo <- modelo_Sph
    modelog <- "Sph"
  }
  if (SCerror == SumaError_Gau) {
    modelo <- modelo_Gau
    modelog <- "Gau"
  }
  modeloele <- variofit(Nuevo, ini = c(modelo$cov.pars[1], 
                                       modelo$cov.pars[2]), cov.model = modelo$cov.model, nugget = modelo$nugget, 
                        fix.nug = TRUE)
  VMod_ele <- (modeloele$cov.pars[1] + modeloele$nugget) - 
    cov.spatial(Nuevo$u, cov.model = modeloele$cov.model, 
                cov.pars = c(modeloele$cov.pars[1], modeloele$cov.pars[2]))
  Error2_ele <- (Nuevo$v - VMod_ele)^2
  SumaError <- sum(Error2_ele)
  MaxX <- max(Nuevo$u, modeloele$cov.pars[2]) + 0.22 * (max(Nuevo$u, 
                                                            modeloele$cov.pars[2]))
  MaxY <- max(Nuevo$v, Nuevo$var.mark, modeloele$cov.pars[1] + 
                modeloele$nugget) + 0.22 * (max(Nuevo$v, Nuevo$var.mark, 
                                                modeloele$cov.pars[1] + modeloele$nugget))
  orden <- matrix(1, ncol = 1, nrow = 1, byrow = T)
  div <- layout(orden, TRUE)
  layout.show(div)
  par(mar = c(6, 6, 3, 3))
  plot(Nuevo$u, Nuevo$v, type = "p", pch = 21, col = "black", 
       bg = "blue", xlim = c(0, MaxX), ylim = c(0, MaxY), xlab = "Distancia [m]", 
       ylab = "Semivarianza", cex.lab = 1.5, cex.axis = 1.2, 
       main = MainTitle)
  grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
  par(new=TRUE)
  plot(Nuevo$u, Nuevo$v, type = "p", pch = 21, col = "black", 
       bg = "blue", xlim = c(0, MaxX), ylim = c(0, MaxY), xlab = "Distancia [m]", 
       ylab = "Semivarianza", cex.lab = 1.5, cex.axis = 1.2, 
       main = MainTitle)
  lines(modeloele, lty = 1, lwd = 2, col = "black")
  abline(h = Nuevo$var.mark, lty = 5, lwd = 2, col = "Red")
  abline(h = modeloele$cov.pars[1] + modeloele$nugget, lty = 4, 
         lwd = 2, col = "darkviolet")
  abline(v = modeloele$cov.pars[2], lty = 3, lwd = 2, col = "Green")
  legend(0, MaxY, c(modeloele$cov.model, "Varianza", "Sill+Nugget", 
                    "Range"), col = c("black", "red", "darkviolet", "green"), 
         lty = c(1, 5, 4, 3), lwd = 2, bty = "n")
  ModMejor <- cbind(modeloele$nugget, modeloele$cov.pars[1] + 
                      modeloele$nugget, modeloele$cov.pars[2], SumaError, max(Nuevo$v), 
                    min(Nuevo$v))
  rownames(ModMejor) <- c(modeloele$cov.model)
  colnames(ModMejor) <- c("Nugget", "Meseta+Nugget", "Alcance", 
                          "SCE", "MaxY", "MinY")
  return(ModMejor)
}