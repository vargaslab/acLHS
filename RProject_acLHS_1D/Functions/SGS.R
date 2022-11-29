SGS <-function (CoorX, CoorY, Prop1, Modelo, Nugget, SillYNugget, Alcance, 
                          minPar, maxPar, Xmin, Xmax, Ymin, Ymax, TX, TY, InvT, NameX, 
                          NameY, Titulo1, Titulo2, n_sim) 
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
                   nmax = maxPar, nmin = minPar, maxdist = Alcance, nsim = n_sim)
  
  kriging <- as.data.frame(kriging)
  # sim1 <- kriging$sim1
  # for (i in length(sim1)) {
  #     if (sim1[i] <= 0) sim1[i] = 0
  # }
  nn <- length(kriging[, 1])
  secuencia <- seq(1:nn)
  unir <- cbind(secuencia, kriging[, 3])
  Quedan <- na.omit(unir)
  cuales <- Quedan[, 1]
  kriging <- kriging[c(cuales), ]
  
  if (InvT == 0) {
    Salida <- as.data.frame(cbind(kriging))
    Salida1<-Salida
    gridded(Salida) <- ~Var1+Var2
    b1<-spplot(Salida, c("sim1"),main=Titulo1, xlab = NameX , ylab = NameY , col.regions=topo.colors, 
               do.log = TRUE,
               key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
               scales=list(draw = TRUE),
               sp.layout = list("sp.points", DatosOrg, pch = 20, col = "black"))
    # b2<-spplot(Salida, c("V4"),main=Titulo2, xlab = NameX , ylab = NameY , col.regions=topo.colors,
    #            do.log = TRUE,
    #            key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
    #            scales=list(draw = TRUE),
    #            sp.layout = list("sp.points", DatosOrg, pch = 20, col = "black"))
    
    # df1 = data.frame(cbind(kriging$Var1, kriging$Var2, kriging$var1.pred))
    # df2 = data.frame(cbind(CoorX, CoorY, Prop1))
    # myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
    # ggplot(df1, aes(x = X1, y = X2)) + geom_tile(aes(fill = X3)) + scale_colour_gradientn(colours = terrain.colors(10)) +
    #   geom_point(data = df2, aes(color = Prop1), shape = 19, size = 2) +
    #   scale_color_gradient(low = "grey", high = "blue")
    
    print(b1)
    #print(b2, split = c(2,1,2,1), more = F)
    
    Salida1 <- as.matrix(Salida1)
    return(Salida1)
    
  }
  if (InvT == 1) {
    Salida <- as.data.frame(cbind(kriging$Var1, kriging$Var2, exp(kriging$sim1)))
    Salida1<-Salida
    gridded(Salida) <- ~V1+V2
    
    b1<-spplot(Salida, c("V3"),main=Titulo1, xlab = NameX , ylab = NameY , col.regions=topo.colors,
               do.log = TRUE,
               key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
               scales=list(draw = TRUE),
               sp.layout = list("sp.points", DatosOrg, pch = 20, col = "black"))  
    # b2<-spplot(Salida, c("V4"),main=Titulo2, xlab = NameX , ylab = NameY , col.regions=topo.colors,
    #            do.log = TRUE,
    #            key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
    #            scales=list(draw = TRUE),
    #            sp.layout = list("sp.points", DatosOrg, pch = 20, col = "black"))
    print(b1) #, split = c(1,1,2,1), more = T)
    #print(b2, split = c(2,1,2,1), more = F)
    
    Salida1 <- as.matrix(Salida1)
    return(Salida1)
    
  }
  if (InvT == 2) {
    Salida <- as.data.frame(cbind(kriging$Var1, kriging$Var2, (kriging$sim1)^2))
    Salida1<-Salida
    gridded(Salida) <- ~V1+V2
    
    b1<-spplot(Salida, c("V3"),main=Titulo1, xlab = NameX , ylab = NameY ,col.regions=topo.colors,
               do.log = TRUE,
               key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
               scales=list(draw = TRUE),
               sp.layout = list("sp.points", DatosOrg, pch = 20, col = "black"))  
    # b2<-spplot(Salida, c("V4"),main=Titulo2, xlab = NameX , ylab = NameY , col.regions=topo.colors,
    #            do.log = TRUE,
    #            key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
    #            scales=list(draw = TRUE),
    #            sp.layout = list("sp.points", DatosOrg, pch = 20, col = "black"))
    print(b1) #, split = c(1,1,2,1), more = T)
    #print(b2, split = c(2,1,2,1), more = F)
    
    Salida1 <- as.matrix(Salida1)
    return(Salida1)
    
  }
  if (InvT == 3) {
    Salida <- as.data.frame(cbind(kriging$Var1, kriging$Var2, 1/(kriging$sim1)))
    Salida1<-Salida
    gridded(Salida) <- ~V1+V2
    
    b1<-spplot(Salida, c("V3"),main=Titulo1, xlab = NameX , ylab = NameY , col.regions=topo.colors,
               do.log = TRUE,
               key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
               scales=list(draw = TRUE),
               sp.layout = list("sp.points", DatosOrg, pch = 20, col = "black"))  
    # b2<-spplot(Salida, c("V4"),main=Titulo2, xlab = NameX , ylab = NameY , col.regions=topo.colors,
    #            do.log = TRUE,
    #            key.space=list(x = 0.1, y = 0.95, corner = c(0, 1)),
    #            scales=list(draw = TRUE),
    #            sp.layout = list("sp.points", DatosOrg, pch = 20, col = "black"))
    print(b1)#, split = c(1,1,2,1), more = T)
    #print(b2, split = c(2,1,2,1), more = F)
    
    Salida1 <- as.matrix(Salida1)
    return(Salida1)
  }
}