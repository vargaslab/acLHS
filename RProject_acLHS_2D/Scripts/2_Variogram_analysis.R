#--------------------------------------------------------------#
#                   Variography Analysis                       #
#--------------------------------------------------------------#
root_dir<-getwd()

# Creates a folder to store results for a case
# dir.create(paste(getwd(),"/Results/spatial_2D", sep=""))
# result_dir<-paste(root_dir,"/Results/spatial_2D",sep="")
# Creates a folder to store results for Variogram Analysis (VA)
dir.create(paste(getwd(),"/Results/spatial_2D/VA", sep=""))
VA_dir<-paste(root_dir,"/Results/spatial_2D/VA",sep="")

#--------------------------------------------------------------#
#                     Spatial Distribution                     #
#--------------------------------------------------------------#
png(paste(VA_dir,"/Rs_Spatial_Distr.png",sep=""), 
    bg = "white", width = 1500, height = 1000, res = 100)
DEspacial(X, Y, Rs, 
          Xlabel , Ylabel, Rslabel, 'Rs Spatial Distribution')
dev.off()


# Estimation of the experimental variogram
X_rng<-X_Stat[8,2]
Y_rng<-Y_Stat[8,2]
N_lags<- 10
lag_value <- sqrt(X_rng*X_rng+Y_rng*Y_rng)/(2*N_lags)
DistMin<-min(dist(Data_Frame[,1:2])) # Minimum distance in data
DistMax<-max(dist(Data_Frame[,1:2])) # Maximum distance in data
lag_value<- DistMin#max((DistMax/2)/N_lags, DistMin)
#lag_value <- 100.0
png(paste(VA_dir,"/Rs_VarioEstimation_10lags.png",sep=""), bg = "white", width = 1500, height = 1000, res = 150)
Rs_VarioEstimation<-Variograma(X, Y, 
                               Rs, 0, 90, N_lags, lag_value, 1, "Variograma adireccionalde Rs", "Distancia [m]")
dev.off()



LR = lm(Rs_VarioEstimation$Semivarianzas~Rs_VarioEstimation$Lags)
png(paste(VA_dir,"/Rs_VarioEyeEstimation.png",sep=""),  bg = "white", width = 1500, height = 1500, res = 200)
plot(Rs_VarioEstimation[,c(2,3)], xlim = c(0,1.1*max(Rs_VarioEstimation$Lags)),
     pch = 19, ylim= c(0,1.1*max(Rs_VarioEstimation$Semivarianzas)))
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Rs_VarioEstimation[,c(2,3)], xlim = c(0,1.1*max(Rs_VarioEstimation$Lags)), 
     pch = 19, ylim= c(0,1.1*max(Rs_VarioEstimation$Semivarianzas)))
abline(LR)
dev.off()
