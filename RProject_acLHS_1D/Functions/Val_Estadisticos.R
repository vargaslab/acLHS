#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
#     NUEVAS METODOLOG?AS Y HERRAMIENTAS DE CARACTERIZACI?N EST?TICA Y DIN?MICA
#        CONSIDERANDO LAS PROPIEDADES FRACTALES DE LOS YACIMIENTOS PETROLEROS 
#                      Proyecto Fondos SENER-CONACYT  No. 143935 (Y-00114)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
#         Esta funci?n obtiene una tabla con todos los valores estadisticos de la variable(s) en estudio
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                        Autor: Javier M?ndez Venegas
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                          Fecha: Marzo de 2014
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                   C?digo fuente de la aplicaci?n:



Val_Estadisticos<-function(datos)
{
    
    require(moments)
    
    tam<-dim(datos)
    
    
    media<-NULL
    Media<-NULL
    mediana<-NULL
    Mediana<-NULL
    cuartil_1er<-NULL
    Cuartil_1er<-NULL
    cuartil_3er<-NULL
    Cuartil_3er<-NULL
    minimo<-NULL
    Minimo<-NULL
    maximo<-NULL
    Maximo<-NULL
    varianza<-NULL
    Varianza<-NULL
    desv<-NULL
    Desv_Estandar<-NULL
    simetria<-NULL
    Simetria<-NULL
    curtosis<-NULL
    Curtosis<-NULL
    rango<-NULL
    Rango<-NULL
    rangoInt<-NULL
    Rango_Intercuartil<-NULL
    muestra<-NULL
    dimension<-NULL
    No_muestras<-NULL
    
    for (i in 1:tam[2])
    {
        media[i]<-round(mean(na.omit(datos[,i])),5)
        Media<-c(media)
        mediana[i]<-round(median(na.omit(datos[,i])),5)
        Mediana<-c(mediana)
        minimo[i]<-round(min(na.omit(datos[,i])),5)
        Minimo<-c(minimo)
        maximo[i]<-round(max(na.omit(datos[,i])),5)
        Maximo<-c(maximo)
        cuartil_1er[i]<-round((quantile(na.omit(datos[,i]),0.25)),5)
        Cuartil_1er<-c(cuartil_1er)
        cuartil_3er[i]<-round((quantile(na.omit(datos[,i]),0.75)),5)
        Cuartil_3er<-c(cuartil_3er)
        varianza[i]<-round(var(na.omit(datos[,i])),5)
        desv[i]<-round(sd(na.omit(datos[,i])),5)
        Varianza<-c(varianza)
        Desv_Estandar<-c(desv)
        simetria[i]<-round(skewness(na.omit(datos[,i])),5)
        Simetria<-c(simetria)
        curtosis[i]<-round(kurtosis(na.omit(datos[,i])),5)
        Curtosis<-c(curtosis)
        rango[i]<-round((max(na.omit(datos[,i]))-min(na.omit(datos[,i]))),5)
        Rango<-c(rango)
        rangoInt[i]<-round(IQR(na.omit(datos[,i])),5)
        Rango_Intercuartil<-c(rangoInt)
        dimension<-length(na.omit(datos[,i]))
        muestra[i]<-dimension[1]
        No_muestras<-c(muestra)
    }
    
    
    Estad<-as.data.frame(rbind(No_muestras,Minimo,Cuartil_1er,Mediana,Media,Cuartil_3er,
                               Maximo,Rango,Rango_Intercuartil,Varianza,Desv_Estandar,Simetria, Curtosis))
    
    names(Estad)<-c(names(datos))
    
    return(Estad)
}





