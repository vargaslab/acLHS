# Polynomial Trend Estimation
# Returns the residual (detrended) values and estimated coefficients
# degree = 1: R = Z - (c1 + c2*X + c3*Y)
# degree = 2: R = Z - (c1 + c2*X**2 + c3*Y**2 + c4*X*Y + c5*X + c6*Y)


Trend<-function (X, Y, Z, degree) 
{
  DatosOrg <- as.data.frame(cbind(X, Y, Z))
  colnames(DatosOrg) <- c("X", "Y", "V1")
  if (degree == 1) {
    Tendencia1 <- lm(V1 ~ X + Y, DatosOrg)
    Salida <- as.data.frame(cbind(DatosOrg[, 1], DatosOrg[, 2], 
                              Tendencia1$residuals, Tendencia1$coefficients[1], 
                              Tendencia1$coefficients[2], Tendencia1$coefficients[3]))
    names(Salida)<-c("X","Y","R","C1","C2","C3")
    return(Salida)
  }
  if (degree == 2) {
    Tendencia1 <- lm(V1 ~ I(X^2) + I(Y^2) + I(X * Y) + X + 
                       Y, DatosOrg)
    Salida <- as.data.frame(cbind(DatosOrg[, 1], DatosOrg[, 2], 
                              Tendencia1$residuals, Tendencia1$coefficients[1], 
                              Tendencia1$coefficients[5], Tendencia1$coefficients[6], 
                              Tendencia1$coefficients[4], Tendencia1$coefficients[2], 
                              Tendencia1$coefficients[3]))
    names(Salida)<-c("X","Y","R","C1","C2","C3","C4","C5","C6")
    return(Salida)
  }
}