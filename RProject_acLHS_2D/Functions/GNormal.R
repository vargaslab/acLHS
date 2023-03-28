GNormal<-function (Data, Property, CooX, CooY, Title) 
{
  DatosGeo <- as.geodata(Data, coords.col = CooX:CooY, data.col = Property)
  orden <- matrix(1, ncol = 1, nrow = 1, byrow = T)
  div <- layout(orden, heights = 10, widths = 10, TRUE)
  layout.show(div)
  par(mar = c(5, 5, 5, 5))
  qqnorm(DatosGeo$data, xlab = "Theoric Quantil", ylab = "Sample Quantil", 
         cex.lab = 1.3, cex.axis = 1, main = Title)
  qqline(DatosGeo$data)
}