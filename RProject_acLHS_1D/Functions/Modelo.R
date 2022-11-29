Modelo<-function (datos, lagUsada, nUsada) 
{
  Nuevo <- Variograma(datos, lagUsada, nUsada)
  max.pepita <- max(Nuevo$v)
  for (i in 1:length(Nuevo$v)) {
    if (Nuevo$v[i] == max.pepita) {
      rango <- Nuevo$u[i]
      rango <- trunc(rango)
    }
  }
  modelo_Exp <- variofit(Nuevo, ini = c(Nuevo$var.mark, rango), 
                         cov.model = "exponential", fix.nug = FALSE)
  modelo_Sph <- variofit(Nuevo, ini = c(Nuevo$var.mark, rango), 
                         cov.model = "spherical", fix.nug = FALSE)
  modelo_Gau <- variofit(Nuevo, ini = c(Nuevo$var.mark, rango), 
                         cov.model = "gaussian", fix.nug = FALSE)
  SCerror <- min(c(modelo_Exp$value, modelo_Sph$value, modelo_Gau$value))
  SCerror
  if (SCerror == modelo_Exp$value) {
    modelo <- modelo_Exp
    modelog <- "Exp"
  }
  if (SCerror == modelo_Sph$value) {
    modelo <- modelo_Sph
    modelog <- "Sph"
  }
  if (SCerror == modelo_Gau$value) {
    modelo <- modelo_Gau
    modelog <- "Gau"
  }
  modeloele <- variofit(Nuevo, ini = c(modelo$cov.pars[1], 
                                       modelo$cov.pars[2]), cov.model = modelo$cov.model, nugget = modelo$nugget, 
                        fix.nug = TRUE)
  orden <- matrix(c(1, 2), ncol = 2, nrow = 1, byrow = T)
  div <- layout(orden, TRUE)
  layout.show(div)
  par(mar = c(5, 5, 5, 5))
  plot(Nuevo$u, Nuevo$v, type = "p", pch = 11, col = "red", 
       xlim = c(0, max(Nuevo$u)), ylim = c(0, max(Nuevo$v)), 
       xlab = "Distancias (Km)", ylab = "Semivarianzas", cex.lab = 1.5, 
       cex.axis = 1.2)
  lines(modelo_Exp, lty = 1, lwd = 2.5, col = "darkorange")
  lines(modelo_Sph, lty = 2, lwd = 2.5, col = "darkblue")
  lines(modelo_Gau, lty = 3, lwd = 2.5, col = "darkgreen")
  abline(h = Nuevo$var.mark, lty = 5, col = "Red")
  legend(20, 0.6, c("Exponencial", "Esférico"   , "Gaussiano"), 
         col = c("darkorange", "darkblue", "darkgreen"), lty = c(1, 
                                                                 2, 3), lwd = 2.5, bty = "n", cex = 1.5)
  legend(0, 1.13, "Varianza", bty = "n", cex = 1.5)
  par(mar = c(5, 5, 5, 5))
  plot(Nuevo$u, Nuevo$v, type = "p", pch = 11, col = "red", 
       xlim = c(0, max(Nuevo$u)), ylim = c(0, max(Nuevo$v)), 
       xlab = "Distancias (Km)", ylab = "Semivarianzas", cex.lab = 1.5, 
       cex.axis = 1.2)
  lines(modeloele)
  abline(h = modeloele$cov.pars[1] + modeloele$nugget, lty = 3, 
         col = "Blue")
  abline(h = Nuevo$var.mark, lty = 5, col = "Red")
  abline(v = modeloele$cov.pars[2], lty = 1, col = "Green")
  legend(20, 0.4, c("Meseta+Nugget", "Varianza", "Alcance"), 
         col = c("Blue", "Red", "Green"), lty = c(3, 5, 1), bty = "n", 
         cex = 1.5)
  print("Modelos Probados")
  print(modelo_Exp)
  print(modelo_Sph)
  print(modelo_Gau)
  print("Modelo Elegido")
  print(modeloele)
}