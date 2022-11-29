# "datos" is
# the input data must be an object of class data.frame
# "Coordinates" is a data.frame or a vector object
# "Round" is the number of decimals to be rounded 
# "Signif" rounds the values in its first argument to the specified number of significant digits. See ?signif
Estadisticas<-function(datos, Coordinates = NULL, Round = 5, Signif = NULL){

    # require(moments)
    # library(e1071)

    datos<-na.omit(datos)
    dimension<-length(datos)
    muestra<-dimension[1]
    muestras<-c(muestra)
    minimo<-min(datos)
    minimos<-c(minimo)
    cuantiles1<-as.double(quantile(datos,0.25))
    mediana<-median(datos)
    medianas<-c(mediana)
    media<-mean(datos)
    medias<-c(media)
    cuantiles3<-as.double(quantile(datos,0.75))
    maximo<-max(datos)
    maximos<-c(maximo)
    rango<-max(datos)-min(datos)
    rangos<-c(rango)
    rangoInt<-IQR(datos)
    rangosInt<-c(rangoInt)
    varianza<-var(datos)
    varianzas<-c(varianza)
    desv<-sd(datos)
    desvs<-c(desv)
    CV<-desv/media
    CVs<-c(CV)
    simetria<-skewness(datos)
    simetrias<-c(simetria)
    curtosis<-kurtosis(datos)
    curtosiss<-c(curtosis)
# var(as.data.frame(cbind(rnorm(100,1,0.7),rnorm(100,1,4),rnorm(100,5,2))))

    if(!is.null(Coordinates)) {
        DistMin<-min(as.vector(dist(Coordinates))) # Minimum distance in data
        DistMax<-max(as.vector(dist(Coordinates))) # Maximum distance in data
        Valores<-as.data.frame(rbind(muestras,DistMax,DistMin,minimos,cuantiles1,medianas,medias,cuantiles3,maximos,rangos,rangosInt,varianzas,desvs,CVs,simetrias,curtosiss))
        Estadigrafos<-c("n","Max. Distance","Min. Distance", "Minimum", "1st. Quartile","Median", "Mean","3rd. Quartile", "Maximum", "Rank","Interquartile Rank", "Variance", "Standard Deviation", "Variation Coeff." ,"Skewness", "Kurtosis")
    } else {
        Valores<-as.data.frame(rbind(muestras,                minimos,cuantiles1,medianas,medias,cuantiles3,maximos,rangos,rangosInt,varianzas,desvs,CVs,simetrias,curtosiss))
        Estadigrafos<-c("n",                                       "Minimum", "1st. Quartile","Median", "Mean","3rd. Quartile", "Maximum", "Rank","Interquartile Rank", "Variance", "Standard Deviation", "Variation Coeff." ,"Skewness", "Kurtosis")
    }

    if(!is.null(Signif)) Valores<-signif(Valores,Signif)
    if(!is.null(Round)) Valores<-round(Valores,Round)
    Estad<-cbind(Estadigrafos,Valores)
    colnames(Estad)<-c("Statistics","Values")

    return(Estad)
}


##### IMPROVEMENTS
# Missing the argument  "ColNames = NULL" as in the function EstadisticasCirculares1
# Display N (the number of data) as an integer, not with zero (0) decimal values
# Include the maximum and minimum distances when position coordinates are given
# # Accept the number of decimals values of the desired table