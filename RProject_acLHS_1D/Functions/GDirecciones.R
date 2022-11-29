GDirecciones<-function (CoorX, CoorY, Variable, x_lab="X Coordinate", y_lab="Y Coordinate", var_lab="Variable Data",
                        main = "Median Regression Analysis in X and Y directions") 
{
  # orden <- matrix(c(1, 2), ncol = 2, nrow = 1, byrow = T)
  # div <- layout(orden, heights = 7.5, widths = c(4.5, 4.5), 
  #               TRUE)
  # layout.show(div)
  # par(mar = c(4.5, 4, 6, 1))
  # plot(Variable ~ CoorX, pch = 21, col = "black", bg = "green", 
  #      ylab = var_lab, xlab = x_lab, cex.lab = 1.4, 
  #      cex.axis = 1.2)
  # abline(h = mean(Variable), lty = 5, col = "red")
  # lines(loess.smooth(CoorX, Variable))
  # legend(min(CoorX), max(Variable), "Media", cex = 1.5, col = "red", 
  #        lty = 5, bty = "n")
  # par(mar = c(4.5, 4, 6, 1))
  # plot(Variable ~ CoorY, pch = 21, col = "black", bg = "green", 
  #      ylab = var_lab, xlab = y_lab, cex.lab = 1.4, 
  #      cex.axis = 1.2)
  # abline(h = mean(Variable), lty = 5, col = "Red")
  # lines(loess.smooth(CoorY, Variable))
  # legend(min(CoorY), max(Variable), "Media", cex = 1.5, col = "Red", 
  #        lty = 5, bty = "n")
  orden <- matrix(c(1, 2), ncol = 2, nrow = 1, byrow = T)
  div <- layout(orden, heights = 7.5, widths = c(4.5, 4.5), 
                TRUE)
  layout.show(div)
  par(mar = c(4.5, 4, 6, 1))
  plot(Variable ~ CoorX, pch = 21, col = "black", bg = "green", 
       ylab = var_lab, xlab = x_lab, cex.lab = 1.4, 
       cex.axis = 1.2)
  grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
  par(new=TRUE)
  plot(Variable ~ CoorX, pch = 21, col = "black", bg = "green", 
       ylab = var_lab, xlab = x_lab, cex.lab = 1.4, 
       cex.axis = 1.2)
  abline(h = mean(Variable), lty = 5, col = "red")
  lines(loess.smooth(CoorX, Variable))
  legend(min(CoorX), max(Variable), legend = c("Media","Mediana"), cex = 1.5, col = c("red", "black"), 
         lty = c(5,1), bty = "n")
  par(mar = c(4.5, 4, 6, 1))
  plot(Variable ~ CoorY, pch = 21, col = "black", bg = "green", 
       ylab = var_lab, xlab = y_lab, cex.lab = 1.4, 
       cex.axis = 1.2)
  grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
  par(new=TRUE)
  plot(Variable ~ CoorY, pch = 21, col = "black", bg = "green", 
       ylab = var_lab, xlab = y_lab, cex.lab = 1.4, 
       cex.axis = 1.2)
  abline(h = mean(Variable), lty = 5, col = "Red")
  lines(loess.smooth(CoorY, Variable))
  legend(min(CoorY), max(Variable), legend = c("Media","Mediana"), cex = 1.5, col = c("red", "black"), 
         lty = c(5,1), bty = "n")
  mtext(main, side = 3, line = -7, outer = TRUE, cex=2,font=2)
}