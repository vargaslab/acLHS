# Cross Validation


CrossValidation<-
function (CoorX, CoorY, Variable, Modelo, Nugget, SillYNugget, Alcance, MaxAnis, proporcion) 
{
library(RandomFields)

  Sill <- SillYNugget - Nugget
  DatosOrg <- as.data.frame(cbind(CoorX, CoorY, Variable))
  colnames(DatosOrg) <- c("X", "Y", "V1")
  coordinates(DatosOrg) <- ~X + Y
  if (Modelo == 1) {
    Modelo <- "Exp"
  }
  if (Modelo == 2) {
    Modelo <- "Sph"
  }
  if (Modelo == 3) {
    Modelo <- "Gau"
  }
  validar <- krige.cv(V1 ~ 1, DatosOrg, vgm(Sill, Modelo, Alcance, 
                                            Nugget, anis = c(MaxAnis, proporcion)), maxdist = Alcance)
  validar <- as.data.frame(validar)
  nn <- length(validar[, 1])
  secuencia <- seq(1:nn)
  unir <- cbind(secuencia, validar[, 1])
  Quedan <- na.omit(unir)
  cuales <- Quedan[, 1]
  Salida <- as.data.frame(cbind(CoorX, CoorY, Variable, validar[, 1], validar[, 
                                                                4]))
  names(Salida)<-c("X","Y","Z", "Z*","Z-Z*" )
  Salida1 <- Salida[c(cuales), ]
  
  return(Salida1)
}