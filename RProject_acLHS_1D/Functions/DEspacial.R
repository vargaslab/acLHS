DEspacial <-
function (CoorX, CoorY, P1, XLAB="X", YLAB="Y", NameP1= deparse(substitute(P1)), MainTitle = "Distribución espacial",
          cex.lab = 1.8, cex.axis = 1.5, cex.main= 1.8,
          AbsFreq = TRUE,Id = FALSE, IdCol = "black",
          # col = gray.colors(64, start = 0.3, end = 0.9, gamma = 2.2),
          Grid=NULL, breaks = NULL, TextPar=list(col="black", cex=1),
          win = NULL)
{
  library(reshape) # for the function "sort_df"
  library(fields) # for the function "as.image"
  
  Datos <- as.data.frame(cbind(CoorX, CoorY, P1))
  DatosOrd <- sort_df(Datos, vars = "P1") # "sort_df" sorts "Datos" in ascending order for "P1"
  
  if (is.null(breaks)) {
    nbins = 9
    Xmax <- max(P1)
    Xmin <- min(P1)
    Range <- Xmax-Xmin
    valhist <- seq(Xmin-0.005*Range, Xmax+0.005*Range, length.out = nbins+1)
  } else {
    valhist <-hist(P1,breaks = breaks,plot=F)$breaks
  }
  nbins<-length(valhist)-1
  ########################### Getting the elements of every bin in the histogram
  # Next, "Datoi" stores the position in "DatosOrd" of the data that belongs to the bin number i
  DATO<-NULL
  DATOS<-NULL
  AbsFreqVarName <- NULL
  for (i in 1:nbins) {
    DATO<-c(DATO,paste("Dato",i, sep=""))
    assign(DATO[i],NULL)
    DATOS<-c(DATOS,paste("Datos",i, sep=""))
    assign(DATOS[i],NULL)
    AbsFreqVarName<-c(AbsFreqVarName ,paste("AbsFreq",i, sep=""))
    assign(AbsFreqVarName[i],0)
  }
  
  N <- length(DatosOrd[,3])
  
  for (i in 1:N) {
    Jmax<-TRUE
    if (valhist[1] <= DatosOrd[i, 3] & DatosOrd[i, 3] <= 
          valhist[1+1]) { # checks if the datum "i-th" is in this interval
      Dato1 <- c(Dato1,i)
      # Datos1 <- c(Dato1)
      DatosOrd[i,4]<-1
      AbsFreq1 <- AbsFreq1+1 # This is the same as i1<-i1+1
      Jmax<-FALSE # "Jmax" indicates if the search must continue (TRUE)for this value of "i"
    }
    if (Jmax) {
      continue<-TRUE # To start the search for values of "j" greater than one
      j<-2
      while(continue){
        if (valhist[j] < DatosOrd[i, 3] & DatosOrd[i, 3] <= 
              valhist[j+1]) {
          assign(DATO[j],c(get(DATO[j]),i))
          DatosOrd[i,4]<-j
          assign(AbsFreqVarName[j],get(AbsFreqVarName[j]) + 1)
          continue<-FALSE
        }
        j<-j+1
      }
    }
  }
  ##############
  #Colors<-gray.colors(nbins, start = 0.3, end = 0.9, gamma = 2.2)
  Colors<- rainbow(nbins,start=0, end=4/6)
  
  #Colors<-c("#0000FFFF","green","violet","darkorange","aquamarine",
  #"darkgreen","gray40","#AA00FFFF","#FF0000FF") 
  # Run "colors()" for other colours. These are the colours in the Javiier's DEspacial function
  
  ###########
  orden <- matrix(c(1, 1, 3, 1, 1, 3, 1, 1, 2), ncol = 3, nrow = 3, 
                  byrow = T)
  div <- layout(orden, widths = c(3, 3, 3), heights = c(0.3, 
                                                        3, 3), TRUE)
  layout.show(div)
  par(mar = c(5, 5, 5, 0))
  
  
  plot.default(DatosOrd[, 1], DatosOrd[, 2], pch = 19, col = "transparent", 
               cex.lab = 1.8, cex.axis = 1.5, cex = 2, xlab = XLAB, 
               ylab = YLAB, main = MainTitle, cex.main = 1.8) # , asp=1
  grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
  par(new=TRUE)
  # Next, ploting the spatial distribution of every bin that has absolute frequency greater than one
  for (i in 1:nbins) {
    if ( get(AbsFreqVarName[i]) > 0) {
      points.default(DatosOrd[get(DATO[i]), 1], DatosOrd[get(DATO[i]), 2], 
                     pch = 20, col = Colors[i], cex = 3)
    }
  }
  ########################## Next, plot as image the spatial distribution. Use this instead of the spatial distribution plotting above
  #    if (is.null(Grid)) {
  #        m<-unique(CoorX)
  #        n<-unique(CoorY)
  #        Grid<-list(x = sort(m),y = sort(n))
  #    }
  #    colnames(DatosOrd)<-c("x","y","z","binId")
  #    imageData<-as.image(Z=DatosOrd[,4],x = DatosOrd[,1:2],grid=Grid) # "DatosOrd[,4]" stores the bins which a given value of Data it belongs
  #   # head(cbind(1:152,DatosOrd,imageData$xd,imageData$ind), n=20)
  #   # imageData$z
  #    image.default(x = imageData, asp = 1, col=Colors, xlab = XLAB, ylab = XLAB, breaks = seq(0.5,nbins+0.5),
  #        main=MainTitle,cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
  #    if (!is.null(win)) plot.owin(win, add=T, main="")
  #    if (Id)
  #       text(imageData$xd,labels = imageData$z[imageData$ind], col=IdCol)
  #
  ############# plotting the histogram
  par(mar = c(5, 5, 1, 2))
  hist2(DatosOrd[, 3], freq = TRUE, breaks = valhist, xlab = NameP1, ylab= "Frecuencia (conteo)",
        main = "",
        cex.lab = 0.8*cex.lab, cex.axis = 0.8*cex.axis,
        AbsFreq = AbsFreq, PercentFreq = F) # "AbsFreq" Is the absolute frequency to be plotted?
  for (i in 1:nbins) {
    if (Id)
      text((valhist[i]+valhist[i+1])/2,2,labels = i, pos=1, col=IdCol)
    rect(xleft = valhist[i], ybottom =0 , xright = valhist[i+1], ytop = get(AbsFreqVarName[i]) ,col=Colors[i])
  }
  ########## PLOTTING THE TEXT
  Id = TRUE
  par(mar = c(2, 0, 5, 2))
  plot.default(0, 0, type = "n", xlim = c(0, 40), ylim = c(0, 42), 
               xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  if (Id){
    text(0,40,labels = "No", cex = 1.3*TextPar$cex, font=2, pos=4)
  }
  text(10, 40, labels = "Min", cex = 1.3*TextPar$cex, font=2)
  text(20, 40, labels = "Bar", cex = 1.3*TextPar$cex, font=2)
  text(30, 40, labels = "Max", cex = 1.3*TextPar$cex, font=2)
  Pos<-seq(36,4,-4)
  for (i in 1:nbins) {
    if (Id){
      text(0,Pos[i],labels = sprintf("%i:",i), pos=4, font=2, col=IdCol)
    }
    text(10, Pos[i], labels = round(valhist[i], 3), cex = TextPar$cex, col = TextPar$col)
    points.default(20, Pos[i],pch = 19, col = Colors[i] , bg= Colors[i], cex = 2)
    text(30, Pos[i], labels = round(valhist[i+1], 3), cex = TextPar$cex, col = TextPar$col)
  }
  #    r<-list(imageData=imageData,breaks=valhist,Data=DatosOrd)
  #    invisible(r)
}