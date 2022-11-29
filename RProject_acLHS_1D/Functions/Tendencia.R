Tendencia<-function (CoorX, CoorY, Prop1, grado) 
{
  DatosOrg <- as.data.frame(cbind(CoorX, CoorY, Prop1))
  colnames(DatosOrg) <- c("X", "Y", "V1")
  if (grado == 1) {
    Tendencia1 <- lm(V1 ~ X + Y, DatosOrg)
    Salida <- as.matrix(cbind(DatosOrg[, 1], DatosOrg[, 2], 
                              Tendencia1$residuals, Tendencia1$coefficients[1], 
                              Tendencia1$coefficients[2], Tendencia1$coefficients[3]))
    return(Salida)
  }
  if (grado == 2) {
    Tendencia1 <- lm(V1 ~ I(X^2) + I(Y^2) + I(X * Y) + X + 
                       Y, DatosOrg)
    Salida <- as.matrix(cbind(DatosOrg[, 1], DatosOrg[, 2], 
                              Tendencia1$residuals, Tendencia1$coefficients[1], 
                              Tendencia1$coefficients[5], Tendencia1$coefficients[6], 
                              Tendencia1$coefficients[4], Tendencia1$coefficients[2], 
                              Tendencia1$coefficients[3]))
    return(Salida)
  }
}