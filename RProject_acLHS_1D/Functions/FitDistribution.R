


# FitDistribution

# Description
# fits and plots statistical distributions

# Arguments
# "data" A numeric vector
# "DISTR" A character string "name",Exponential    (exp), Log Normal (lnorm),
# Normal (norm), Uniform (unif), Weibull (weibull).
# "BREAKS" Number of classes or "scott",  "FD","Sturges"
# "col" the color of the histogram bins
# "DistName" the name (in spanish) of the distribution


# Value
# a list with three elements:
# 1) "x" a data.frame with some hyphothesis test
# 2) "y" the fit parameters
# 3) "DistName" the name (in spanish) of the distribution


FitDistribution<- function(data, DISTR, BREAKS, col, DistName)
{
    require(fitdistrplus)
    require(actuar)
    require(ADGofTest)
    
    BREAK<-BREAKS
    
    if (class(BREAKS)=="numeric")
    {
        
        Property<-na.omit(data)
        Xmax<-max(Property)
        Xmin<-min(Property)
        difhis<-Xmax-Xmin+0.002*(Xmax-Xmin)
        numclases<-BREAKS
        tamint<-difhis/numclases
        valorhis<-0
        vectorhis<-0
        for (i in 1:numclases)
        {
            valorhis[i]<-Xmin+i*tamint-0.001*(Xmax-Xmin)
            vectorhis<-c(Xmin-0.001*(Xmax-Xmin), valorhis) 
        }
        BREAK<-vectorhis
    }   
    
    
    
    DISTRIBUTION<-NULL
    if (DISTR=="exp")
    {DISTRIBUTION<-pexp}
    if (DISTR=="gamma")
    {DISTRIBUTION<-pgamma}
    if (DISTR=="lnorm")
    {DISTRIBUTION<-plnorm}
    if (DISTR=="norm")
    {DISTRIBUTION<-pnorm}
    if (DISTR=="unif")
    {DISTRIBUTION<-punif}
    if (DISTR=="weibull")
    {DISTRIBUTION<-pweibull}
    DISTRIBUTION
    
    
    Fit<-fitdist(data, DISTR, method = "mle")
    plot(Fit, breaks=BREAK, col=col)
    
    if(DISTR=="exp")
        
    {    Kolm<-ks.test(data,DISTRIBUTION,Fit$estimate[1])
         Anderson<-ad.test(data, DISTRIBUTION,Fit$estimate[1])
    }
    if(DISTR!="exp")
    {
        Kolm<-ks.test(data,DISTRIBUTION,Fit$estimate[1],Fit$estimate[2])
        Anderson<-ad.test(data, DISTRIBUTION,Fit$estimate[1],Fit$estimate[2])
    }
    
    p_value<-cbind(Kolm[2],Anderson[2])
    Estadistico<-cbind(Kolm[1],Anderson[1])
    
    Decision<-NULL
    Decisiones<-NULL
    
    for (i in 1:length(p_value))
    {
        if (p_value[i]<0.05)
        {Decision<-"Rechazo H0"
         Decisiones<-c(Decisiones,Decision)}
        if (p_value[i]>=0.05)
        {Decision<-"No rechazo H0"
         Decisiones<-c(Decisiones,Decision)}
    }
    
    Pre<-NULL
    TestNames <- c("Kolmogorov-Smirnov","Anderson-Darling")
    Pre<-as.data.frame(cbind(TestNames,c(0.05,0.05),format(p_value,digits=4),format(Estadistico, digits=4),Decisiones))
    colnames(Pre)<-c("Nombre","Nivel de significancia","P-valor","Estadístico","Decisión")
    
    Pre2<-NULL
    Pre2<-as.data.frame(rbind(as.numeric(Fit$estimate[1]),as.numeric(Fit$estimate[2]),as.numeric(Fit$loglik),as.numeric(Fit$aic)))
    rownames(Pre2)<-c("Media","Desviación estándar","Máxima Verosimilitud","AICC")
    colnames(Pre2)<-c(DistName)
    
    r <-list(x=Pre,y=Pre2)
    return(r)   
    
}





# #Ejemplo
# Familia<- readtable2(,header= TRUE,sep="") 
# dim(Familia)
# names(Familia)
# 
# FitDistribution(Familia$length, "norm",  15,  "blue",  "Normal")
# 
# 
# FitDistribution(Familia$length, "lnorm",  "Sturges",  "blue",  "Normal")
# 
# 
# FitDistribution(Familia$length, "exp",  "FD",  "blue",  "Exponencial")
# 
# 
# FitDistribution(Familia$length, "weibull",  "scott",  "blue",  "Weibull")






