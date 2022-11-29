## Javier?s Function taken from Rgeoestad
function (CoorX, CoorY, Variable, CoorX_sub, CoorY_sub, Variable_sub, Direccion=0, Tol=90, NIntervalos=10, 
    Lags, Npares=1, MainTitle="Variograma", xlab="Distancia", ylab = "Semivarianza") 
{
    # CoorX <- t
    # CoorY <- ti
    # Variable <- CO2
    # CoorX_sub <-  cbind(t_sub_CO2,t_temporal_sub)
    # CoorY_sub <-  cbind(ti_sub_CO2,ti_temporal_sub)
    # Variable_sub <- cbind(CO2_temp_opt_subsample,CO2_temporal_subsample)
    # Direccion=0
    # Tol=90
    # NIntervalos= N_lags 
    # Lags = lag_value
    # Npares=1
    # MainTitle = "Empirical variogram"
    # xlab = "Time (days)"
    # ylab = "Semivariogram"
    # 
    library(geoR)

    DatosGeo <- as.geodata(cbind(CoorX, CoorY, Variable), header = TRUE, 
        coords.col = 1:2, data.col = 3)
    variograma <- variog(DatosGeo, breaks = c(seq(0, Lags * NIntervalos, 
        Lags)), trend = "cte", lambda = 1, estimator.type = "classical", 
        nugget.tolerance = 0, direction = Direccion, tolerance = Tol, 
        unit.angle = "degrees", pairs.min = Npares)
    DatosGeo_sub <- as.geodata(cbind(CoorX_sub[,1], CoorY_sub[,1], Variable_sub[,1]), header = TRUE, 
                               coords.col = 1:2, data.col = 3)
    variograma_sub <- variog(DatosGeo_sub, breaks = c(seq(0, Lags * NIntervalos, 
                                                          Lags)), trend = "cte", lambda = 1, estimator.type = "classical", 
                             nugget.tolerance = 0, direction = Direccion, tolerance = Tol, 
                             unit.angle = "degrees", pairs.min = Npares)
    DatosGeo_sub2 <- as.geodata(cbind(CoorX_sub[,2], CoorY_sub[,2], Variable_sub[,2]), header = TRUE, 
                                coords.col = 1:2, data.col = 3)
    variograma_sub2 <- variog(DatosGeo_sub2, breaks = c(seq(0, Lags * NIntervalos, 
                                                            Lags)), trend = "cte", lambda = 1, estimator.type = "classical", 
                              nugget.tolerance = 0, direction = Direccion, tolerance = Tol, 
                              unit.angle = "degrees", pairs.min = Npares)
    MaxX <-1.5*max(variograma$u, variograma_sub$u, variograma_sub2$u) # MaxX <- max(variograma$u) + 0.2 * max(variograma$u)
    MaxY <-1.5*max(variograma$v, variograma$var.mark, variograma_sub$v, variograma_sub$var.mark,
                   variograma_sub2$v, variograma_sub2$var.mark) # MaxY <- max(variograma$v, variograma$var.mark) + 0.2 * (max(variograma$v, variograma$var.mark))
    orden <- matrix(1, ncol = 1, nrow = 1, byrow = T)
    div <- layout(orden, TRUE)
    layout.show(div)
    par(mar = c(5, 5, 5, 5))
    plot.default(variograma$u, variograma$v, type = "p", pch = 19, col = "black", 
        bg = "black", xlim = c(0, MaxX), ylim = c(0, MaxY), xlab = xlab, 
        ylab = ylab, cex.lab = 1.2, cex.axis = 1.2, cex.main=1.5,
        main = MainTitle)
    grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
    par(new=TRUE)
    plot.default(variograma$u, variograma$v, type = "p", pch = 19, col = "black", 
                 bg = "black", xlim = c(0, MaxX), ylim = c(0, MaxY), xlab = xlab, 
                 ylab = ylab, cex.lab = 1.2, cex.axis = 1.2, cex.main=1.5,
                 main = MainTitle)
    abline(h = variograma$var.mark, lty = 5, lwd = 2, col = "black")
    par(new=TRUE)
    plot.default(variograma_sub$u, variograma_sub$v, type = "p", pch = 19, col = "red", 
                     bg = "red", xlim = c(0, MaxX), ylim = c(0, MaxY), xlab = xlab, 
                     ylab = ylab, cex.lab = 1.2, cex.axis = 1.2, cex.main=1.5,
                     main = MainTitle)
    abline(h = variograma_sub$var.mark, lty = 5, lwd = 2, col = "red")
    par(new=TRUE)
    plot.default(variograma_sub2$u, variograma_sub2$v, type = "p", pch = 19, col = "blue", 
                 bg = "blue", xlim = c(0, MaxX), ylim = c(0, MaxY), xlab = xlab, 
                 ylab = ylab, cex.lab = 1.2, cex.axis = 1.2, cex.main=1.5,
                 main = MainTitle)
    abline(h = variograma_sub2$var.mark, lty = 5, lwd = 2, col = "blue")
    #text(variograma$u, variograma$v, as.character(variograma$n), cex=1, pos=3, col="blue") # 
    legend(x = "topleft", y = NULL, c("data ","sps", "sts", "data variance","sps variance",
                                      "sts variance"), col = c("black","red","blue","black","red","blue"), 
        lty = c(NA,NA,NA,5,5,5), lwd = c(NA,NA,NA,2,2,2), pch = c(19,19,19,NA,NA,NA),  bty = "n")
    #text(1,10,"- Varianza Total                                ", col="red", cex=1)
    Salida <- as.data.frame(cbind(variograma$n, variograma$u, variograma$v))
    names(Salida)<-c("Npares","Lags","Semivarianzas")
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