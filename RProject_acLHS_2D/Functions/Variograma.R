## Javier?s Function taken from Rgeoestad

Variograma<-
function (CoorX, CoorY, Variable, Direccion=0, Tol=90, NIntervalos=10, 
    Lags, Npares=1, MainTitle="Variograma", xlab="Distancia", ylab = "Semivarianza") 
{
    library(geoR)

    DatosGeo <- as.geodata(cbind(CoorX, CoorY, Variable), header = TRUE, 
        coords.col = 1:2, data.col = 3)
    variograma <- variog(DatosGeo, breaks = c(seq(Lags/2, Lags * (NIntervalos+1), 
        Lags)), trend = "cte", lambda = 1, estimator.type = "classical", 
        nugget.tolerance = 0, direction = Direccion, tolerance = Tol, 
        unit.angle = "degrees", pairs.min = Npares)
    variograma
    MaxX <-1.2*max(variograma$u) # MaxX <- max(variograma$u) + 0.2 * max(variograma$u)
    MaxY <-1.2*max(variograma$v, variograma$var.mark) # MaxY <- max(variograma$v, variograma$var.mark) + 0.2 * (max(variograma$v, variograma$var.mark))
    orden <- matrix(1, ncol = 1, nrow = 1, byrow = T)
    div <- layout(orden, TRUE)
    layout.show(div)
    par(mar = c(5, 5, 5, 5))
    plot.default(variograma$u, variograma$v, type = "p", pch = 21, col = "black", 
        bg = "green", xlim = c(0, MaxX), ylim = c(0, MaxY), xlab = xlab, 
        ylab = ylab, cex.lab = 1.2, cex.axis = 1.2, cex.main=1.5,
        main = MainTitle)
    grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
    par(new=TRUE)
    plot.default(variograma$u, variograma$v, type = "p", pch = 21, col = "black", 
                 bg = "green", xlim = c(0, MaxX), ylim = c(0, MaxY), xlab = xlab, 
                 ylab = ylab, cex.lab = 1.2, cex.axis = 1.2, cex.main=1.5,
                 main = MainTitle)
    abline(h = variograma$var.mark, lty = 5, lwd = 2, col = "red")
    text(variograma$u, variograma$v, as.character(variograma$n), cex=1, pos=3, col="blue") # 
    legend("right", "Variance", col = "red", # x = "topleft", y = NULL
        lty = 5, lwd = 2, bty = "n")
    #text(1,10,"- Varianza Total                                ", col="red", cex=1)
    Salida <- as.data.frame(cbind(variograma$n, variograma$u, variograma$v))
    names(Salida)<- c("Npares","Lags","Semivarianzas")
    return(Salida)
}

### EXAMPLES
#  variograma <- variog(DatosGeo, breaks = c(seq(0, Lags * NIntervalos, 
#          Lags)), trend = "cte", lambda = 1, estimator.type = "classical", 
#          nugget.tolerance = 0, direction = Direccion, tolerance = Tol, 
#          unit.angle = "degrees", pairs.min = Npares)
   
## converting the data-set "topo" from the package MASS (VR?s bundle)
## to the geodata format:
#  topo
#  plot(topo)
#  DistMin<-min(dist(topo[,1:2])) # Minimum distance in data
#  Variogram<-Variograma(CoorX=topo[,1], CoorY=topo[,2], Variable=topo[,3], Direccion=0, Tol=90, NIntervalos=10, #  #  Lags=DistMin, Npares=1, MainTitle="Semivariogram Topo data")

# To check
#  CoorX=topo[,1]; CoorY=topo[,2]; Variable=topo[,3]; Direccion=0; Tol=90; NIntervalos=10; Lags=DistMin; Npares=1; #  #  MainTitle="Semivariogram Topo data"
#  topogeo <- as.geodata(topo)
# plot.geodata(topogeo)
#  variograma <- variog(topogeo,, breaks = c(seq(0, Lags * NIntervalos, 
#          Lags)), trend = "cte", lambda = 1, estimator.type = "classical", 
#          nugget.tolerance = 0, direction = Direccion, tolerance = Tol, 
#          unit.angle = "degrees", pairs.min = Npares)
#  plot.default(variograma$u,variograma$v)
#  all.equal(variograma$v, Variogram[,3])