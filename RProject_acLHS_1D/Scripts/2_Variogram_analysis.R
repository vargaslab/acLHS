#--------------------------------------------------------------#
#                   Variogram Analysis                         #
#--------------------------------------------------------------#

root_dir<-getwd()

# Creates a folder to store results for a case
# dir.create(paste(getwd(),"/Results/Time_Series", sep=""))
# result_dir<-paste(root_dir,"/Results/Time_Series",sep="")
# Create a folder to store results for  Variogram Analysis (VA)
dir.create(paste(getwd(),"/Results/Time_Series/VA", sep=""))
VA_dir<-paste(root_dir,"/Results/Time_Series/VA",sep="")

#--------------------------------------------------------------#
#                    Temporal Distribution                     #
#--------------------------------------------------------------#
log_list=list(Temp,CO2)
n_logs= length(log_list)
log_xlabel_list=list(Templabel,CO2label)
plot_mean=TRUE
plot_median=TRUE
fontproportion=1.0
png(paste(VA_dir,"/Temp_CO2_Temporal_Distr.png",sep=""), 
    bg = "white", width = 1500, height = 1000, res = 150 )
par(mfrow = c(n_logs,1),mar = c(4, 6, 2, 2), cex.lab= 1.0, cex.axis = 1.0)
for (i in 1:n_logs) {
  p=unlist(log_list[i])
  log_label=unlist(log_xlabel_list[i])
  s=summary(p)
  plot(Time, p,  type = "b", pch = 19,
       ylab = log_label , xlab = Timelabel, bty="o")
  grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
  par(new=TRUE)
  plot(Time, p, type = "b", pch = 19,
       ylab =log_label , xlab = Timelabel, bty="o")
  if (plot_mean) {
    abline(h = s[4], col = "Red")
  }  
  if (plot_median) {
    abline(h = s[3], col = "Blue", lty = 2)
  }
  legend("topleft", legend=c("Mean", "Median"),
         col=c("red", "blue"), lty=1:2, box.lty=0)
  box()
}
dev.off()

#--------------------------------------------------------------#
#               Trend (Stationarity) Analysis                  #
#--------------------------------------------------------------#

# Estimation of the experimental variogram
N_lags<- 8       
lag_value <- 1   # min(dist(Time)) # or delta t
ti = numeric(length(Time))
png(paste(VA_dir,"/CO2_VarioEstimation.png",sep=""), bg = "white", width = 1500, height = 1000, res = 200)
CO2_VarioEstimation<-Variograma(Time, ti, 
                                         CO2, 0, 90, N_lags, lag_value, 1, "", Timelabel) # CO2 temporal variogram
dev.off()

#--------------------------------------------------------------#
#       Univariate Variography (Structural)  Modeling          #
#--------------------------------------------------------------#

# Manual Variogram Model Fitting
# Variogram models (1- exponential, 2- spherical, 3- gaussian)
N_lags<- 8 #length(t)/2.5
lag_value <- 1 # min(dist(t)) # or delta t
CO2_vario_model<- 2 #3
CO2_nugget<- 0.57 #0.1 #0.57 #0.4 #1
CO2_sill_and_nugget<- 1.95 #2.4 #1.95 #1.9 #29
CO2_rank <- 6.21 #6 #6.21 #5
png(paste(VA_dir,"/CO2_VarioEyeEstimation.png",sep=""),   bg = "white", width = 1500, height = 1000, res = 150)
CO2_EyeModelVarioFit<-EyeModel(Time, ti, 
                                           CO2, 0, 90, N_lags, lag_value, 1, 
                                    CO2_vario_model, CO2_nugget, CO2_sill_and_nugget, CO2_rank,
                                    "") # A fitted temporal variogram model
dev.off()

#Variogram models (1- exponential, 2- spherical, 3- gaussian)
N_lags<- 8
lag_value <- 1 
CO2_aclhs_sub48_vario_model<- 2 
CO2_aclhs_sub48_nugget<- 0.57 
CO2_aclhs_sub48_sill_and_nugget<- 1.95 
CO2_aclhs_sub48_range <- 6.21
png(paste(VA_dir,"/CO2_Variogram.png",sep=""), bg = "white",  width = 1500, height = 1000, res = 200)
par(mfrow = c(1,1),mar = c(4, 6, 2, 2), cex.lab= 1.2, cex.axis = 1.2)
plot(CO2_VarioEstimation[,c(2,3)], pch = 1, xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*1.7)),
     xlab = Timelabel ,ylab = Variolabel)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(CO2_VarioEstimation[,c(2,3)], pch = 1, xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*1.7)),
     xlab = Timelabel ,ylab = Variolabel )
# par(new=TRUE)
# plot(CO2_fixed_sub48_VarioEstimation[,c(2,3)], pch = 15, cex = 1.2, col = "green", xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*1.7)),
#      xlab = Timelabel ,ylab = Variolabel)
# par(new=TRUE)
# plot(CO2_clhs_sub48_VarioEstimation[,c(2,3)], pch = 17, cex = 1.2, col = "red",  xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*1.7)),
#      xlab = Timelabel , ylab = Variolabel)
# par(new=TRUE)
# plot(CO2_aclhs_sub48_VarioEstimation[,c(2,3)], pch = 18, cex = 1.5, col = "blue", xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*1.7)),
#      xlab = Timelabel , ylab = Variolabel)
par(new=TRUE)
lines.variomodel(cov.model = "sph", cov.pars = c(CO2_aclhs_sub48_sill_and_nugget-CO2_aclhs_sub48_nugget, CO2_aclhs_sub48_range)
                 , nug = CO2_aclhs_sub48_nugget, col = "Black", lwd = 3, max.dist = 1*8)
legend("topleft",legend=c("Empirical variogram", "Variogram model") # ,
       ,col=c("black", "black"),lty = c(1,1),lwd = c(NA, 3), pch = c(1, NA), box.lty=0) # , title ="Legend"
# legend("topleft",legend=c("Data", "Model",
#                            "FTS sample","cLHS sample","acLHS sample") # ,
#        ,col=c("black", "black", "green", "red", "blue"),lty = c(1,1,1,1,1),lwd = c(NA, 3, NA, NA, NA), pch = c(19, NA, 19, 19,19), box.lty=0) # , title ="Legend"
box()
dev.off()
