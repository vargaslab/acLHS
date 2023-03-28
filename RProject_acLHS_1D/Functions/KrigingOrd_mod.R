

KrigingOrd_mod<-function (CoorX, CoorY, Prop1, Modelo, Nugget, SillYNugget, Alcance, 
    minPar, maxPar, Xmin, Xmax, Ymin, Ymax, TX, TY, InvT, NameX, 
    NameY, Titulo1, Titulo2) 
{
    DatosOrg <- as.data.frame(cbind(CoorX, CoorY, Prop1))
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
    localizacion <- as.data.frame(expand.grid(seq(Xmin, Xmax, 
        by = TX), seq(Ymin, Ymax, by = TY)))
    gridded(localizacion) <- ~Var1 + Var2
    Sill <- SillYNugget - Nugget
    modelocal <- vgm(Sill, Modelo, Alcance, Nugget)
    kriging <- krige(V1 ~ 1, DatosOrg, localizacion, model = modelocal, 
        nmax = maxPar, nmin = minPar, maxdist = Alcance, indicators = FALSE)

    kriging <- as.data.frame(kriging)
    nn <- length(kriging[, 1])
    secuencia <- seq(1:nn)
    unir <- cbind(secuencia, kriging[, 3])
    Quedan <- na.omit(unir)
    cuales <- Quedan[, 1]
    kriging <- kriging[c(cuales), ]

    if (InvT == 0) {
        Salida <- as.data.frame(cbind(kriging$Var1, kriging$Var2, kriging$var1.pred, 
            sqrt(kriging$var1.var)))
 	  Salida1<-Salida
       gridded(Salida) <- ~V1+V2
  
	  b1<-spplot(Salida, c("V3"),main=Titulo1, xlab = NameX , ylab = NameY , col.regions=topo.colors,
	             at=seq(min(Salida$V3), max(Salida$V3), 2),
	             do.log = TRUE,
	             key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
	             scales=list(draw = TRUE),
	             sp.layout = list("sp.points", DatosOrg, pch = 20, col = "black"))  
    b2<-spplot(Salida, c("V4"),main=Titulo2, xlab = NameX , ylab = NameY , col.regions=topo.colors,
                   do.log = TRUE,
                   key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
                   scales=list(draw = TRUE),
                   sp.layout = list("sp.points", DatosOrg, pch = 20, col = "black"))
	  print(b1, split = c(1,1,2,1), more = T)
	  print(b2, split = c(2,1,2,1), more = F)

        Salida1 <- as.matrix(Salida1)
        return(Salida1)

    }
    if (InvT == 1) {
        Salida <- as.data.frame(cbind(kriging$Var1, kriging$Var2, exp(kriging$var1.pred), 
            sqrt(kriging$var1.var)))
	 Salida1<-Salida
       gridded(Salida) <- ~V1+V2
  
	  b1<-spplot(Salida, c("V3"),main=Titulo1, xlab = NameX , ylab = NameY , col.regions=topo.colors,
	             at=seq(min(Salida$V3), max(Salida$V3), 2),
	             do.log = TRUE,
	             key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
	             scales=list(draw = TRUE),
	             sp.layout = list("sp.points", DatosOrg, pch = 20, col = "black"))  
    b2<-spplot(Salida, c("V4"),main=Titulo2, xlab = NameX , ylab = NameY , col.regions=topo.colors,
               do.log = TRUE,
               key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
               scales=list(draw = TRUE),
               sp.layout = list("sp.points", DatosOrg, pch = 20, col = "black"))
	  print(b1, split = c(1,1,2,1), more = T)
	  print(b2, split = c(2,1,2,1), more = F)

        Salida1 <- as.matrix(Salida1)
        return(Salida1)

    }
    if (InvT == 2) {
        Salida <- as.data.frame(cbind(kriging$Var1, kriging$Var2, (kriging$var1.pred)^2, 
            sqrt(kriging$var1.var)))
 	  Salida1<-Salida
       gridded(Salida) <- ~V1+V2
  
	  b1<-spplot(Salida, c("V3"),main=Titulo1, xlab = NameX , ylab = NameY ,col.regions=topo.colors,
	             at=seq(min(Salida$V3), max(Salida$V3), 2),
	             do.log = TRUE,
	             key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
	             scales=list(draw = TRUE),
	             sp.layout = list("sp.points", DatosOrg, pch = 20, col = "black"))  
    b2<-spplot(Salida, c("V4"),main=Titulo2, xlab = NameX , ylab = NameY , col.regions=topo.colors,
               do.log = TRUE,
               key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
               scales=list(draw = TRUE),
               sp.layout = list("sp.points", DatosOrg, pch = 20, col = "black"))
	  print(b1, split = c(1,1,2,1), more = T)
	  print(b2, split = c(2,1,2,1), more = F)

        Salida1 <- as.matrix(Salida1)
        return(Salida1)

    }
    if (InvT == 3) {
        Salida <- as.data.frame(cbind(kriging$Var1, kriging$Var2, 1/(kriging$var1.pred), 
            sqrt(kriging$var1.var)))
 	  Salida1<-Salida
       gridded(Salida) <- ~V1+V2
  
	  b1<-spplot(Salida, c("V3"),main=Titulo1, xlab = NameX , ylab = NameY , col.regions=topo.colors,
	             at=seq(min(Salida$V3), max(Salida$V3), 2),
	             do.log = TRUE,
	             key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
	             scales=list(draw = TRUE),
	             sp.layout = list("sp.points", DatosOrg, pch = 20, col = "black"))  
    b2<-spplot(Salida, c("V4"),main=Titulo2, xlab = NameX , ylab = NameY , col.regions=topo.colors,
               do.log = TRUE,
               key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
               scales=list(draw = TRUE),
               sp.layout = list("sp.points", DatosOrg, pch = 20, col = "black"))
	  print(b1, split = c(1,1,2,1), more = T)
	  print(b2, split = c(2,1,2,1), more = F)

        Salida1 <- as.matrix(Salida1)
        return(Salida1)

    }
}

























