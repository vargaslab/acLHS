ValidacionCross<-function (CoorX, CoorY, P1, P2, NInt, lags, Modelo, Sill1, Sill2, 
          Sill3, Nugget1, Nugget2, Nugget3, Alcance) 
{
  Datos <- as.data.frame(cbind(CoorX, CoorY, P1, P2))
  colnames(Datos) <- c("X", "Y", "Prop1", "Prop2")
  coordinates(Datos) <- ~X + Y
  ambas <- gstat(NULL, id = "Prop1", formula = Prop1 ~ 1, data = Datos)
  ambas <- gstat(ambas, id = "Prop2", formula = Prop2 ~ 1, 
                 data = Datos)
  v.cross <- variogram(ambas, cutoff = NInt * lags, width = lags)
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
  if (Modelo == 1) {
    Modelo <- "Exp"
  }
  if (Modelo == 2) {
    Modelo <- "Sph"
  }
  if (Modelo == 3) {
    Modelo <- "Gau"
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
  validar <- gstat.cv(modelo, nfold = length(P1))
  validar
  validar <- as.data.frame(validar)
  nn <- length(validar[, 1])
  secuencia <- seq(1:nn)
  unir <- cbind(secuencia, validar[, 1])
  Quedan <- na.omit(unir)
  cuales <- Quedan[, 1]
  Salida <- cbind(CoorX, CoorY, P1, validar[, 1], validar[, 
                                                          4])
  Salida1 <- Salida[c(cuales), ]
  return(Salida1)
}