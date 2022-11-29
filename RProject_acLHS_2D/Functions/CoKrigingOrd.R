

CoKrigingOrd<-function (CoorX, CoorY, P1, P2, Modelo, Sill1, Sill2, Sill3, 
    Nugget1, Nugget2, Nugget3, Alcance, minPar, maxPar, Xmin, 
    Xmax, Ymin, Ymax, TX, TY, InvT, NameX, NameY, Titulo1, Titulo2) 
{
    Datos <- as.data.frame(cbind(CoorX, CoorY, P1, P2))
    colnames(Datos) <- c("X", "Y", "Prop1", "Prop2")
    coordinates(Datos) <- ~X + Y
    ambas <- gstat(NULL, id = "Prop1", formula = Prop1 ~ 1, data = Datos)
    ambas <- gstat(ambas, id = "Prop2", formula = Prop2 ~ 1, 
        data = Datos)
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
        Alcance, Nugget1), nmax = maxPar, nmin = minPar, maxdist = Alcance)
    modelo <- gstat(modelo, id = "Prop2", model = vgm(Sill2, 
        Modelo, Alcance, Nugget2), nmax = maxPar, nmin = minPar, 
        maxdist = Alcance)
    modelo <- gstat(modelo, id = c("Prop1", "Prop2"), model = vgm(Sill3, 
        Modelo, Alcance, Nugget3), nmax = maxPar, nmin = minPar, 
        maxdist = Alcance)
    localizacion <- as.data.frame(expand.grid(seq(Xmin, Xmax, 
        by = TX), seq(Ymin, Ymax, by = TY)))
    coordinates(localizacion) <- ~Var1 + Var2
    kriging <- predict(modelo, localizacion)
    kriging <- as.data.frame(kriging)
    nn <- length(kriging[, 1])
    secuencia <- seq(1:nn)
    unir <- cbind(secuencia, kriging[, 3])
    Quedan <- na.omit(unir)
    cuales <- Quedan[, 1]
    kriging <- kriging[c(cuales), ]
    if (InvT == 0) {
     
        Salida <- as.data.frame(cbind(kriging$Var1, kriging$Var2, kriging$Prop1.pred, 
            sqrt(kriging$Prop1.var)))
        Salida1<-Salida
       gridded(Salida) <- ~V1+V2
  
	  b1<-spplot(Salida, c("V3"),main=Titulo1, xlab = NameX , ylab = NameY , col.regions=topo.colors,
	             do.log = TRUE,
	             key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
	             scales=list(draw = TRUE),
	             sp.layout = list("sp.points", Datos, pch = 20, col = "black"))  
    b2<-spplot(Salida, c("V4"),main=Titulo2, xlab = NameX , ylab = NameY , col.regions=topo.colors,
                   do.log = TRUE,
                   key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
                   scales=list(draw = TRUE),
                   sp.layout = list("sp.points", Datos, pch = 20, col = "black"))
	  print(b1, split = c(1,1,2,1), more = T)
	  print(b2, split = c(2,1,2,1), more = F)

        Salida1 <- as.matrix(Salida1)
        return(Salida1)


    }
    if (InvT == 1) {

        Salida <- as.data.frame(cbind(kriging$Var1, kriging$Var2, exp(kriging$Prop1.pred), 
            sqrt(kriging$Prop1.var)))
	Salida1<-Salida
       gridded(Salida) <- ~V1+V2
  
	  b1<-spplot(Salida, c("V3"),main=Titulo1, xlab = NameX , ylab = NameY , col.regions=topo.colors,
	             do.log = TRUE,
	             key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
	             scales=list(draw = TRUE),
	             sp.layout = list("sp.points", Datos, pch = 20, col = "black"))  
    b2<-spplot(Salida, c("V4"),main=Titulo2, xlab = NameX , ylab = NameY , col.regions=topo.colors,
               do.log = TRUE,
               key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
               scales=list(draw = TRUE),
               sp.layout = list("sp.points", Datos, pch = 20, col = "black"))
	  print(b1, split = c(1,1,2,1), more = T)
	  print(b2, split = c(2,1,2,1), more = F)

        Salida1 <- as.matrix(Salida1)
        return(Salida1)
    }
    if (InvT == 2) {
       
        Salida <- as.data.frame(cbind(kriging$Var1, kriging$Var2, (kriging$Prop1.pred)^2, 
            sqrt(kriging$Prop1.var)))
       Salida1<-Salida
       gridded(Salida) <- ~V1+V2
  
	  b1<-spplot(Salida, c("V3"),main=Titulo1, xlab = NameX , ylab = NameY , col.regions=topo.colors,
	             do.log = TRUE,
	             key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
	             scales=list(draw = TRUE),
	             sp.layout = list("sp.points", Datos, pch = 20, col = "black"))  
    b2<-spplot(Salida, c("V4"),main=Titulo2, xlab = NameX , ylab = NameY , col.regions=topo.colors,
               do.log = TRUE,
               key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
               scales=list(draw = TRUE),
               sp.layout = list("sp.points", Datos, pch = 20, col = "black"))
	  print(b1, split = c(1,1,2,1), more = T)
	  print(b2, split = c(2,1,2,1), more = F)

        Salida1 <- as.matrix(Salida1)
        return(Salida1)

    }
    if (InvT == 3) {
       
        Salida <- as.data.frame(cbind(kriging$Var1, kriging$Var2, 1/(kriging$Prop1.pred), 
            sqrt(kriging$Prop1.var)))
        Salida1<-Salida
       gridded(Salida) <- ~V1+V2
  
	  b1<-spplot(Salida, c("V3"),main=Titulo1, xlab = NameX , ylab = NameY , col.regions=topo.colors,
	             do.log = TRUE,
	             key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
	             scales=list(draw = TRUE),
	             sp.layout = list("sp.points", Datos, pch = 20, col = "black"))  
    b2<-spplot(Salida, c("V4"),main=Titulo2, xlab = NameX , ylab = NameY , col.regions=topo.colors,
               do.log = TRUE,
               key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
               scales=list(draw = TRUE),
               sp.layout = list("sp.points", Datos, pch = 20, col = "black"))
	  print(b1, split = c(1,1,2,1), more = T)
	  print(b2, split = c(2,1,2,1), more = F)

        Salida1 <- as.matrix(Salida1)
        return(Salida1)


    }
}








