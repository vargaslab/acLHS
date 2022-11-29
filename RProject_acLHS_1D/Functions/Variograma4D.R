## Modified from Javier?s Function taken from Rgeoestad

Variograma4D<-
function (CoorX, CoorY, Variable, D1 = 0, D2 = 45, D3 = 90, D4 = 135, Tol = 22.5, NIntervalos = 10, 
    Lags, Npares = 1, MainTitle, xlab="Distancia", ylab = "Semivarianza") 
{

    library(geoR)

    DatosGeo <- as.geodata(cbind(CoorX, CoorY, Variable), header = TRUE, 
        coords.col = 1:2, data.col = 3)
    variogramaD1 <- variog(DatosGeo, breaks = c(seq(0, Lags * 
        NIntervalos, Lags)), trend = "cte", lambda = 1, estimator.type = "classical", 
        nugget.tolerance = 0, direction = D1, tolerance = Tol, 
        unit.angle = "degrees", pairs.min = Npares)
    variogramaD2 <- variog(DatosGeo, breaks = c(seq(0, Lags * 
        NIntervalos, Lags)), trend = "cte", lambda = 1, estimator.type = "classical", 
        nugget.tolerance = 0, direction = D2, tolerance = Tol, 
        unit.angle = "degrees", pairs.min = Npares)
    variogramaD3 <- variog(DatosGeo, breaks = c(seq(0, Lags * 
        NIntervalos, Lags)), trend = "cte", lambda = 1, estimator.type = "classical", 
        nugget.tolerance = 0, direction = D3, tolerance = Tol, 
        unit.angle = "degrees", pairs.min = Npares)
    variogramaD4 <- variog(DatosGeo, breaks = c(seq(0, Lags * 
        NIntervalos, Lags)), trend = "cte", lambda = 1, estimator.type = "classical", 
        nugget.tolerance = 0, direction = D4, tolerance = Tol, 
        unit.angle = "degrees", pairs.min = Npares)
    MaxY <- max(variogramaD1$v, variogramaD2$v, variogramaD3$v, 
        variogramaD4$v, variogramaD1$var.mark) + 0.2 * (max(variogramaD1$v, 
        variogramaD2$v, variogramaD3$v, variogramaD4$v, variogramaD1$var.mark))
    MaxX <- max(variogramaD1$u, variogramaD2$u, variogramaD3$u, 
        variogramaD4$u) + 0.2 * (max(variogramaD1$u, variogramaD2$u, 
        variogramaD3$u, variogramaD4$u))
    orden <- matrix(1, ncol = 1, nrow = 1, byrow = T)
    div <- layout(orden, TRUE)
    layout.show(div)
    par(mar = c(5, 5, 5, 5))
    plot(variogramaD1$u, variogramaD1$v, type = "b", pch = 1, 
        col = "black", xlim = c(0, MaxX), ylim = c(0, MaxY), 
        xlab = xlab , ylab = ylab , cex.lab = 1.8, 
        cex.axis = 1.2, lwd = 2, lty = 1, main = MainTitle)
    grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
    par(new=TRUE)
    plot(variogramaD1$u, variogramaD1$v, type = "b", pch = 1, 
         col = "black", xlim = c(0, MaxX), ylim = c(0, MaxY), 
         xlab = xlab , ylab = ylab , cex.lab = 1.8, 
         cex.axis = 1.2, lwd = 2, lty = 1, main = MainTitle)
    matplot(variogramaD2$u, variogramaD2$v, type = "b", pch = 5, 
        col = "blue", xlim = c(0, MaxX), ylim = c(0, MaxY), cex.lab = 1.8, 
        cex.axis = 1.2, add = T, lwd = 2, lty = 1)
    matplot(variogramaD3$u, variogramaD3$v, type = "b", pch = 2, 
        col = "green", xlim = c(0, MaxX), ylim = c(0, MaxY), 
        cex.lab = 1.8, cex.axis = 1.2, add = T, lwd = 2, lty = 1)
    matplot(variogramaD4$u, variogramaD4$v, type = "b", pch = 8, 
        col = "darkviolet", xlim = c(0, MaxX), ylim = c(0, MaxY), 
        cex.lab = 1.8, cex.axis = 1.2, add = T, lwd = 2, lty = 1)
    abline(h = variogramaD1$var.mark, pch = 20, lty = 5, lwd = 2, 
        col = "red")
    legend(x = "topleft", y = NULL, c(D1, D2, D3, D4, "Varianza"), col = c("black", 
        "blue", "green", "darkviolet", "red"), lty = c(1, 1, 1, 1, 5), lwd = 2, bty = "n")
    
    #text(3,110,"- Variograma Direccion   0 grados", col="black", cex=1)
    #text(3,100,"- Variograma Direccion  45 grados", col="blue", cex=1)
    #text(3,90,"- Variograma Direccion  90 grados", col="green", cex=1)
    #text(3,80,"- Variograma Direccion 135 grados", col="darkviolet", cex=1)
    #text(3,70,"- Varianza Total                                ", col="red", cex=1)
    
    Salida1 <- cbind(variogramaD1$n, variogramaD1$u, variogramaD1$v)
    Salida2 <- cbind(variogramaD2$n, variogramaD2$u, variogramaD2$v)
    Salida3 <- cbind(variogramaD3$n, variogramaD3$u, variogramaD3$v)
    Salida4 <- cbind(variogramaD4$n, variogramaD4$u, variogramaD4$v)
###------------{{{{{{{{{{  START MODIFICATION
    Salida <- list(Zero=Salida1, FortyFive=Salida2, Ninety= Salida3, OneThertyFive=Salida4)
###------------}}}}}}}}}}   END MODIFICATION
    return(Salida)
}
