

BasicStats<-function (Data, Property, BREAKS, MainTitle, XLAB, Xmin, Xmax) 
{
    
Property<-na.omit(Property)
Xmax<-Xmax
Xmin<-Xmin

difhis<-Xmax-Xmin+0.002*(Xmax-Xmin)
numclases<- BREAKS
tamint<-difhis/numclases
valorhis<-0
vectorhis<-0
for (i in 1:numclases)
  {
   valorhis[i]<-Xmin+i*tamint-0.001*(Xmax-Xmin)
   vectorhis<-c(Xmin-0.001*(Xmax-Xmin), valorhis) 
  }
vectorhis
    cual <- 0
    cuales <- 0
    for (i in 1:length(Property)) {
        if (Property[i] < min(vectorhis) | Property[i] > max(vectorhis)) {
            cual[i] <- i
            cuales <- c(cual)
        }
        if (Property[i] >= min(vectorhis) & Property[i] <= max(vectorhis)) {
            cual[i] <- 0
            cuales <- c(cual)
        }
    }
    if (sum(cuales) != 0) {
        Property <- Property[-cuales]
    }
    Elementos <- 0
    cuantos <- 0
    for (i in 1:length(vectorhis) - 1) {
        uno <- Property[vectorhis[i] <= Property]
        uno
        Elementos[i] <- length(uno[uno < vectorhis[i + 1]])
        cuantos <- c(Elementos)
    }
    maxY <- max(cuantos)
    orden <- matrix(c(2, 3, 1, 4), ncol = 1, nrow = 4, byrow = T)
    div <- layout(orden, widths = 9, heights = c(3.5, 0.5, 2.5, 
        2.8))
    layout.show(div)
    par(mar = c(2.5, 6, 0, 3))
    plot(Property, seq(1, length(Property), 1), pch = 22, xlim = c(min(vectorhis), 
        max(vectorhis)), col = "black", bg = "#AAFF00FF", ylab = "Number of Observations", 
        xlab = XLAB, cex.lab = 1.3, cex.axis = 1.2)
    abline(v = mean(Property), lty = 4, col = "Red")
    abline(v = median(Property), lty = 5, col = "Blue")
    legend((max(Property) * 0.6), length(Property), c("Mean", 
        "Median"), cex = 1.5, col = c("Red", "Blue"), lty = c(4, 
        5), bty = "n")
    par(mar = c(0, 6, 3, 3))
    histograma <- hist(Property, breaks = vectorhis, xaxt = "n", 
        ylim = c(0, maxY + 4), col = "#FFAA00FF", ylab = "Frequency", 
        cex.lab = 1.3, cex.axis = 1.2, main = MainTitle, labels = TRUE, 
        cex.main = 2, cex = 3)
    abline(v = mean(Property), lty = 3, col = "Red")
    abline(v = median(Property), lty = 5, col = "Blue")
    par(mar = c(0, 6, 0, 3))
    plot(0, 0, type = "n", xlim = c(min(vectorhis), max(vectorhis)), 
        ylim = c(0, 1.54), xaxt = "n", yaxt = "n", xlab = "", 
        ylab = "")
    boxplot(Property, range = 1.5, ylim = c(min(vectorhis), max(vectorhis)), 
        horizontal = TRUE, col = "#AAFF00FF", pch = 22, axes = FALSE, 
        add = TRUE, outline = TRUE, varwidth = TRUE)
    abline(v = mean(Property), lty = 3, col = "Red")
    abline(v = median(Property), lty = 5, col = "Blue")
    summaryP1 <- t(t(summary(Property)))
    par(mar = c(0.5, 6, 0, 3))
    plot(0, 0, type = "n", xlim = c(0, 60), ylim = c(0, 7), xaxt = "n", 
        yaxt = "n", xlab = "", ylab = "")
    text(20, 5.5, labels = "            Min", cex = 1.5)
    text(20, 4.5, labels = "1° Quartile", cex = 1.5)
    text(20, 3.5, labels = "     Median", cex = 1.5)
    text(20, 2.5, labels = "        Mean", cex = 1.5)
    text(20, 1.5, labels = "3° Quartile", cex = 1.5)
    text(20, 0.5, labels = "          Max", cex = 1.5)
    text(32, 6.9, labels = XLAB, cex = 1.5)
    text(34, 5.5, labels = sprintf("%.5f", summaryP1[1, 1]), 
        cex = 1.5)
    text(34, 4.5, labels = sprintf("%.5f", summaryP1[2, 1]), 
        cex = 1.5)
    text(34, 3.5, labels = sprintf("%.5f", summaryP1[3, 1]), 
        cex = 1.5)
    text(34, 2.5, labels = sprintf("%.5f", summaryP1[4, 1]), 
        cex = 1.5)
    text(34, 1.5, labels = sprintf("%.5f", summaryP1[5, 1]), 
        cex = 1.5)
    text(34, 0.5, labels = sprintf("%.5f", summaryP1[6, 1]), 
        cex = 1.5)
    box("outer", lty = "solid", col = "black")
}






