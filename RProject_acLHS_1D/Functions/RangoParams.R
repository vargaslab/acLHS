RangoParams<-function (CoorX, CoorY, Variable, Direccion, Tol, NIntervalos, 
          Lags, Npares) 
{
  DatosGeo <- as.geodata(cbind(CoorX, CoorY, Variable), header = TRUE, 
                         coords.col = 1:2, data.col = 3)
  Nuevo <- variog(DatosGeo, breaks = c(seq(0, Lags * NIntervalos, 
                                           Lags)), trend = "cte", lambda = 1, estimator.type = "classical", 
                  nugget.tolerance = 0, direction = Direccion, tolerance = Tol, 
                  unit.angle = "degrees", pairs.min = Npares)
  MinParams <- min(Nuevo$v)
  MaxParams <- max(Nuevo$v)
  parametros <- as.matrix(cbind(MinParams, MaxParams))
  return(parametros)
}