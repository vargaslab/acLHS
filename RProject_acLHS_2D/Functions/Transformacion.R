Transformacion<-function (Datos, Variable, x) 
{
  if (x == 1) {
    Transformados <- log(Variable)
  }
  if (x == 2) {
    Transformados <- sqrt(Variable)
  }
  if (x == 3) {
    Transformados <- 1/(Variable)
  }
  Datos <- as.matrix(cbind(Datos, Transformados))
  return(Datos)
}