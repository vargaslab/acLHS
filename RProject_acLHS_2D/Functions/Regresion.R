Regresion<-function (CoorX, CoorY, Dependiente, Independiente) 
{
  regresion <- lm(Dependiente ~ Independiente)
  estimadores <- lsfit(Independiente, Dependiente, intercept = TRUE)
  resultados <- ls.diag(estimadores)
  b0 <- as.real(regresion$coefficients[1])
  b1 <- as.real(regresion$coefficients[2])
  Beta0 <- c(b0, rep(0, (length(regresion$residuals) - 1)))
  Beta1 <- c(b1, rep(0, (length(regresion$residuals) - 1)))
  Salida <- as.matrix(cbind(Beta0, Beta1, CoorX, CoorY, Dependiente, 
                            regresion$fitted.values, regresion$residuals))
  return(Salida)
}