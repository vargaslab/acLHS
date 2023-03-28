GDEspacial<-function (datos) 
{
  require(sp)
  require(geoR)
  DatosGeo <- as.geodata(datos, header = TRUE, coords.col = 1:2, 
                         data.col = 3)
  orden <- matrix(1, ncol = 1, nrow = 1, byrow = T)
  div <- layout(orden, heights = 10, widths = 10, TRUE)
  layout.show(div)
  par(mar = c(5, 5, 5, 5))
  plot(DatosGeo$coord[, 2] ~ DatosGeo$coord[, 1], pch = 22, 
       col = "black", bg = "aquamarine", ylab = "Coordenada Y", 
       xlab = "Coordanada X", cex.lab = 1.4, cex.axis = 1.2)
  title("Distribución Espacial de la Información"      )
}