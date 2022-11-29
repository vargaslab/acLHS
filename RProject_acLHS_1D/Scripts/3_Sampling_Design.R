# Create a folder to store results for Sampling Methods (SM)
dir.create(paste(getwd(),"/Results/Time_Series/SM", sep=""))
SM_dir<-paste(root_dir,"/Results/Time_Series/SM",sep="")
#----------------------------------Fixed Temporal Sampling-------------------------------------------#
fixed_sub48 <- round(seq(1,length(Time),length.out = 48))
Time_fixed_sub48 <- Time[fixed_sub48]
Temp_fixed_sub48 <- Temp[fixed_sub48]
CO2_fixed_sub48 <- CO2[fixed_sub48]
# Stats
CO2_fixed_sub48_Stat <- Estadisticas(CO2_fixed_sub48)
Temp_fixed_sub48_Stat <- Estadisticas(Temp_fixed_sub48)
# fixed_Temporal_Sample <- cbind(fixed_sub48,Time_fixed_sub48,Temp_fixed_sub48,CO2_fixed_sub48)
# write.csv(fixed_Temporal_Sample , file = paste(SM_dir,"/fixed_Temporal_Sample.csv",sep=""))
#-----------------------------------------cLHS-------------------------------------------------------#
n_sub = 48
df = data.frame(Temp,CO2)
# Returning the indices of the sampled points
# res <- clhs(df, CoordXY = cbind(Time,ti), size = n_sub, use.cpp = FALSE, progress = FALSE, simple = TRUE)
# str(res)
# Returning a clhs_result object for plotting using C++
system.time(res_ <- clhs(df, size = n_sub, use.cpp = TRUE, iter = 500000, progress = TRUE, simple = FALSE))
res_
str(res_)
plot(res_)
#clhs_Temporal_Sample <- read.csv(file=file.choose(),header=T,na.strings="-999.25")
#clhs_sub48 <- clhs_Temporal_Sample$clhs_sub48
#clhs_sub48 <- res_$index_samples
Time_clhs_sub48 <- Time[clhs_sub48]
Temp_clhs_sub48 <- Temp[clhs_sub48]
CO2_clhs_sub48 <- CO2[clhs_sub48]
# Stats
CO2_clhs_sub48_Stat <- Estadisticas(CO2_clhs_sub48)
Temp_clhs_sub48_Stat <- Estadisticas(Temp_clhs_sub48)
# clhs_Temporal_Sample <- cbind(clhs_sub48,Time_clhs_sub48,Temp_clhs_sub48,CO2_clhs_sub48)
# write.csv(clhs_Temporal_Sample , file = paste(SM_dir,"/clhs_Temporal_Sample.csv",sep=""))

#-------------------------------------------------acLHS---------------------------------------------#
n_sub = 48
Fn_limit <- seq(0,1,length.out = n_sub+1)
Temp_limit <- quantile(Temp,Fn_limit)
CO2_limit <- quantile(CO2,Fn_limit)
Time_limit <- quantile(Time,Fn_limit)
# OF1
FuncO1  <- function(var_sub) {
  pos_sub <- NULL
  for (i in 1:n_sub) {
    pos_sub <- rbind(pos_sub,which.min(abs(CO2-var_sub[i])))
  }
  Temp_sub <- Temp[pos_sub]
  hist_Temp_sub <- hist(Temp_sub, breaks = Temp_limit, plot = FALSE)$counts
  frecuencia <- rep(1,n_sub)
  return (sum(abs(hist_Temp_sub-frecuencia))) # initial = 8 
}

corP_Temp_CO2 <- round(cor(Temp,CO2, method = "pearson"),3)
corS_Temp_CO2 <- round(cor(Temp,CO2, method = "spearman"),3)
corK_Temp_CO2 <- round(cor(Temp,CO2, method = "kendall"),3)
# OF2
FuncO2  <- function(var_sub) {
  pos_sub <- NULL
  for (i in 1:n_sub) {
    pos_sub <- rbind(pos_sub,which.min(abs(CO2-var_sub[i])))
  }
  corP_Temp_Vsub <- round(cor(Temp[pos_sub],CO2[pos_sub], method = "pearson"),3)
  corS_Temp_Vsub <- round(cor(Temp[pos_sub],CO2[pos_sub], method = "spearman"),3)
  corK_Temp_Vsub <- round(cor(Temp[pos_sub],CO2[pos_sub], method = "kendall"),3)
  return (abs(corP_Temp_Vsub-corP_Temp_CO2) + abs(corS_Temp_Vsub-corS_Temp_CO2) + abs(corK_Temp_Vsub-corK_Temp_CO2)) 
}

#OF3
FuncO3 <- function(var_sub) {
  N_lags<- 8 
  lag_value <- 1 # min(dist(Time))
  pos_sub <- NULL
  for (i in 1:n_sub) {
    pos_sub <- rbind(pos_sub,which.min(abs(CO2-var_sub[i])))
  }
  t_sub <- Time[pos_sub]
  ti_sub <- numeric(length(t_sub))
  Vario_var_sub <- variog(as.geodata(cbind(t_sub, ti_sub, var_sub), 
                                     coords.col = 1:2, data.col = 3),  breaks = c(seq(lag_value/2, lag_value * (N_lags+1), lag_value)),
                          trend = "cte", lambda = 1, estimator.type = "classical", nugget.tolerance = 0, direction = 0, 
                          tolerance = 90, unit.angle = "degrees", pairs.min = 1)
  return (sum(abs(Vario_var_sub$v-CO2_VarioEstimation[,3]))) # # initial = 11.20084
}

# OF = w1 * OF1 + w2 * OF2 + w3 * OF3
w1 = 0.3; w2 = 100; w3 = 1 #w1 = 0.2; w2 = 10; w3 = 1
FuncConj  <- function(var_sub) {
  total = w1*FuncO1(var_sub) + w2*FuncO2(var_sub) + w3*FuncO3(var_sub) 
  #           4                      0.009             0.5336005
  #           8                      0.001             0.4868278
  #           8                      0.002             1.451617
  return (total) 
}

# optimization
Fn_opt_sub <- seq(1/n_sub,1,length.out = n_sub)
CO2_optimal_sub <- unname(quantile(CO2,Fn_opt_sub)) 
CO2_order <- CO2[order(CO2)]
pos_otp_order <- NULL
for (i in 1:n_sub) {
  pos_otp_order <- rbind(pos_otp_order,which.min(abs(CO2_order-CO2_optimal_sub[i])))
}
tol_sub <-  round(length(CO2)/n_sub)
Pop_Inicial <- matrix(data=NA,nrow=tol_sub,ncol=n_sub)
for (i in 1:(n_sub)) {
  for (j in 1:tol_sub) {
    pos <- pos_otp_order[i]
    Pop_Inicial[j,i] <- CO2_order[pos-j]
  }
}
Fn_sub_limit <- seq(0,1,1/n_sub)
CO2_sub_limit <- unname(quantile(CO2, probs = Fn_sub_limit))
subvector_min = CO2_sub_limit[-(n_sub+1)]
subvector_max = CO2_sub_limit[-(1)]
system.time(aclhs_sub_outDEoptim_CO2 <- DEoptim(fn=FuncConj, lower= subvector_min,
                                                  upper= subvector_max,
                                                  control = DEoptim.control(VTR = 0.000001,strategy = 3, 
                                                   itermax =10000, reltol = 1e-8, CR = 0.5, F = 0.8,
                                                  NP= nrow(Pop_Inicial),initialpop = Pop_Inicial))) #, 

# user  system elapsed
# 148.531   9.554 158.825   itermax =5000,  n_sub = 48, 3functions=1+1+1,  aclhs_sub_outDEoptim_CO2$optim$bestval = 0.9856357
# 363.728  45.779 412.902   itermax =10000,  n_sub = 48, 3functions=0.2+100+1,  aclhs_sub_outDEoptim_CO2$optim$bestval = 2.186828
# 492.715  64.620 548.035   itermax =15000,  n_sub = 48, 3functions=0.2+100+1,  aclhs_sub_outDEoptim_CO2$optim$bestval =  2.709265
# 679.072  89.088 753.199   itermax =20000,  n_sub = 48, 3functions=0.2+100+1,  aclhs_sub_outDEoptim_CO2$optim$bestval =  1.913174
# 396.71   49.35  455.25    itermax =10000,  n_sub = 48, 3functions=0.3+100+1,  aclhs_sub_outDEoptim_CO2$optim$bestval =  2.954214
aclhs_sub_outDEoptim_CO2$optim$bestval 
Ite <- seq(1,length(aclhs_sub_outDEoptim_CO2$member$bestvalit),1)
plot(Ite,aclhs_sub_outDEoptim_CO2$member$bestvalit, "l", xlab ="Iterations", ylab ="Objective Function")
var_sub <- as.numeric(aclhs_sub_outDEoptim_CO2$optim$bestmem) 

Estadisticas(var_sub)
pos_sub_CO2 <- NULL
for (i in 1:n_sub) {
  pos_sub_CO2 <- rbind(pos_sub_CO2,which.min(abs(CO2-var_sub[i])))
}

# aclhs_Temporal_Sample <- read.csv(file=file.choose(),header=T,na.strings="-999.25")
# aclhs_sub48 <- aclhs_Temporal_Sample[,2]
aclhs_sub48 <- pos_sub_CO2
Time_aclhs_sub48 <- Time[aclhs_sub48]
Temp_aclhs_sub48 <- Temp[aclhs_sub48]
CO2_aclhs_sub48 <- CO2[aclhs_sub48]
# Stats
CO2_aclhs_sub48_Stat <- Estadisticas(CO2_aclhs_sub48)
Temp_aclhs_sub48_Stat <- Estadisticas(Temp_aclhs_sub48)
# aclhs_Temporal_Sample <- cbind(aclhs_sub48,Time_aclhs_sub48,Temp_aclhs_sub48,CO2_aclhs_sub48)
# write.csv(aclhs_Temporal_Sample , file = paste(SM_dir,"/aclhs_Temporal_Sample.csv",sep=""))
#--------------------Figure 2------------------------------#
png(paste(SM_dir,"/Process_optimization.png",sep=""), bg = "white",  width = 1500, height = 1500, res = 250)
plot(Ite,aclhs_sub_outDEoptim_CO2$member$bestvalit,
     "l", xlab = expression(bold("Iterations")), ylab = expression(bold("Objective Function")))
dev.off()

#--------------------------------Figure 3--------------------------------------#
log_list=list(Temp,CO2)
n_logs= length(log_list)
log_xlabel_list=list(Templabel,CO2label)
plot_mean=TRUE
plot_median=TRUE
fontproportion=1.0

#n_logs= 4 #length(log_list)
plot_mean=TRUE
plot_median=TRUE
fontproportion=1.0
png(paste(SM_dir,"/CO2_Temporal_Distr_Sampling_days.png",sep=""), 
    bg = "white", width = 2000, height = 1000, res = 150 )
par(mfrow = c(2,2),mar = c(4, 6, 2, 2), cex.lab= 1.0, cex.axis = 1.0)
plot(Time, CO2,  type = "p", lwd = 2, pch = 1, 
     ylab = CO2label , xlab = Timelabel, bty="o")
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Time, CO2, type = "p", lwd = 2, pch = 1, 
       ylab =CO2label , xlab = Timelabel, bty="o")
#
plot(Time, CO2,  type = "p", lwd = 2, pch = 1,
     ylab = CO2label , xlab = Timelabel, bty="o", col = "lightgray")
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Time_fixed_sub48,CO2_fixed_sub48, pch = 15, cex = 1.2, bty="o", col = "green", xlim = c(min(Time), max(Time)), ylim = c(min(CO2), max(CO2)), xlab =Timelabel, ylab = CO2label)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Time_fixed_sub48,CO2_fixed_sub48, pch = 15, cex = 1.2, bty="o", col = "green", xlim = c(min(Time), max(Time)), ylim = c(min(CO2), max(CO2)), xlab =Timelabel, ylab = CO2label)
#
plot(Time, CO2,  type = "p", lwd = 2, pch = 1,
     ylab = CO2label , xlab = Timelabel, bty="o", col = "lightgray")
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Time_clhs_sub48,CO2_clhs_sub48, pch = 17, cex = 1.2, bty="o", col = "red", xlim = c(min(Time), max(Time)), ylim=c(min(CO2), max(CO2)), xlab =Timelabel, ylab = CO2label)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Time_clhs_sub48,CO2_clhs_sub48, pch = 17, cex = 1.2, bty="o", col = "red", xlim = c(min(Time), max(Time)), ylim=c(min(CO2), max(CO2)), xlab =Timelabel, ylab = CO2label)
#
plot(Time, CO2,  type = "p", lwd = 2, pch = 1,
     ylab = CO2label , xlab = Timelabel, bty="o", col = "lightgray")
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Time_aclhs_sub48,CO2_aclhs_sub48, pch = 18,  cex = 1.5,  bty="o", col = "blue", xlim = c(min(Time), max(Time)), ylim=c(min(CO2), max(CO2)), xlab =Timelabel, ylab = CO2label)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Time_aclhs_sub48,CO2_aclhs_sub48, pch = 18,  cex = 1.5, bty="o", col = "blue", xlim = c(min(Time), max(Time)), ylim=c(min(CO2), max(CO2)), xlab =Timelabel, ylab = CO2label)

# legend("topright", legend=c("Data","FTS sample", "cLHS sample", "acLHS sample"),
#          col=c("black","green","red", "blue"), pch = c(19,15,17,18), box.lty=0)
box()
dev.off()

#---------------------------------------------Figure4_Univariate pdfs--------------------------------#
#Fn_Temp <- estandarizar(cbind(Temp,CO2))[,1]
#Fn_CO2 <- estandarizar(cbind(Temp,CO2))[,2]
U_limit <- seq(0,1,by = 0.2)
png(paste(SM_dir,"/Temp_univariate_pdfs.png",sep=""),   bg = "white", width = 1200, height = 1000, res = 120)
par(mfrow = c(2,2),mar = c(4, 6, 2, 2), cex.lab= 1.2, cex.axis = 1.2)
#-------------------Figure a----------------#
plot(ecdf(Temp), pch = 1, yaxt = "n",  main = " ", xlab= Templabel, 
     ylab = expression(bold(paste(Fn(Temperature)))),xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
axis(2, at=U_limit)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ecdf(Temp), pch = 1, yaxt = "n",  main = " ", xlab= Templabel, 
     ylab = expression(bold(paste(Fn(Temperature)))),xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
axis(2, at=U_limit)
box()

#-------------------Figure b----------------#
plot(ecdf(Temp), pch = 1, yaxt = "n",  main = " ", xlab= Templabel,  col = "lightgray",
     ylab = expression(bold(paste(Fn(Temperature)))),xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
axis(2, at=U_limit)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ecdf(Temp), pch = 1, yaxt = "n",  main = " ", xlab= Templabel, col = "lightgray",
     ylab = expression(bold(paste(Fn(Temperature)))),xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
axis(2, at=U_limit)
par(new=TRUE)
plot(ecdf(Temp_fixed_sub48),pch = 15, col = "green", cex = 1.2,
     yaxt = "n",  main = " ", xlab= Templabel, ylab = expression(bold(paste(Fn(Temperature)))),
     xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
box()

#-------------------Figure c----------------#
plot(ecdf(Temp), pch = 1, yaxt = "n",  main = " ", xlab= Templabel, col = "lightgray",
     ylab = expression(bold(paste(Fn(Temperature)))),xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
axis(2, at=U_limit)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ecdf(Temp), pch = 1, yaxt = "n",  main = " ", xlab= Templabel, col = "lightgray",
     ylab = expression(bold(paste(Fn(Temperature)))),xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
axis(2, at=U_limit)
par(new=TRUE)
plot(ecdf(Temp_clhs_sub48),pch = 17, col = "red", cex = 1.2,
     yaxt = "n",  main = " ", xlab= Templabel, ylab = expression(bold(paste(Fn(Temperature)))),
     xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
box()

#-------------------Figure d----------------#
plot(ecdf(Temp), pch = 1, yaxt = "n",  main = " ", xlab= Templabel, col = "lightgray",
     ylab = expression(bold(paste(Fn(Temperature)))),xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
axis(2, at=U_limit)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ecdf(Temp), pch = 1, yaxt = "n",  main = " ", xlab= Templabel, col = "lightgray",
     ylab = expression(bold(paste(Fn(Temperature)))),xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
axis(2, at=U_limit)
par(new=TRUE)
plot(ecdf(Temp_aclhs_sub48),pch = 18, col = "blue", cex = 1.5,
     yaxt = "n",  main = " ", xlab= Templabel, ylab = expression(bold(paste(Fn(Temperature)))),
     xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
box()
dev.off()

#---------------------------------------------Figure5_Univariate pdfs--------------------------------#
#Fn_Temp <- estandarizar(cbind(Temp,CO2))[,1]
#Fn_CO2 <- estandarizar(cbind(Temp,CO2))[,2]
U_limit <- seq(0,1,by = 0.2)
png(paste(SM_dir,"/CO2_univariate_pdfs.png",sep=""),   bg = "white", width = 1200, height = 1000, res = 120)
par(mfrow = c(2,2),mar = c(4, 6, 2, 2), cex.lab= 1.2, cex.axis = 1.2)
#-------------------Figure a----------------#
plot(ecdf(CO2), pch = 1, yaxt = "n",  main = " ", xlab= CO2label, 
     ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),xlim =c(min(CO2),max(CO2)), ylim = c(0,1))
axis(2, at=U_limit)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ecdf(CO2), pch = 1, yaxt = "n",  main = " ", xlab= CO2label, 
     ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),xlim =c(min(CO2),max(CO2)), ylim = c(0,1))
axis(2, at=U_limit)
box()

#-------------------Figure b----------------#
plot(ecdf(CO2), pch = 1, yaxt = "n",  main = " ", xlab= CO2label, col = "lightgray",
     ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),xlim =c(min(CO2),max(CO2)), ylim = c(0,1))
axis(2, at=U_limit)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ecdf(CO2), pch = 1, yaxt = "n",  main = " ", xlab= CO2label, col = "lightgray",
     ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),xlim =c(min(CO2),max(CO2)), ylim = c(0,1))
axis(2, at=U_limit)
par(new=TRUE)
plot(ecdf(CO2_fixed_sub48),pch = 15, col = "green", cex = 1.2,
     yaxt = "n",  main = " ", xlab= CO2label, ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),
     xlim =c(min(CO2),max(CO2)), ylim = c(0,1))
box()

#-------------------Figure c----------------#
plot(ecdf(CO2), pch = 1, yaxt = "n",  main = " ", xlab= CO2label, col = "lightgray",
     ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),xlim =c(min(CO2),max(CO2)), ylim = c(0,1))
axis(2, at=U_limit)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ecdf(CO2), pch = 1, yaxt = "n",  main = " ", xlab= CO2label, col = "lightgray",
     ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),xlim =c(min(CO2),max(CO2)), ylim = c(0,1))
axis(2, at=U_limit)
par(new=TRUE)
plot(ecdf(CO2_clhs_sub48),pch = 17, col = "red", cex = 1.2,
     yaxt = "n",  main = " ", xlab= CO2label, ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),
     xlim =c(min(CO2),max(CO2)), ylim = c(0,1))
box()

#-------------------Figure d----------------#
plot(ecdf(CO2), pch = 1, yaxt = "n",  main = " ", xlab= CO2label, col = "lightgray",
     ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),xlim =c(min(CO2),max(CO2)), ylim = c(0,1))
axis(2, at=U_limit)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ecdf(CO2), pch = 1, yaxt = "n",  main = " ", xlab= CO2label, col = "lightgray",
     ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),xlim =c(min(CO2),max(CO2)), ylim = c(0,1))
axis(2, at=U_limit)
par(new=TRUE)
plot(ecdf(CO2_aclhs_sub48),pch = 18, col = "blue", cex = 1.5,
     yaxt = "n",  main = " ", xlab= CO2label, ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),
     xlim =c(min(CO2),max(CO2)), ylim = c(0,1))
box()
dev.off()


#---------------------------------------------Figure6_scatterplots--------------------------------#
Fn_Temp <- estandarizar(cbind(Temp,CO2))[,1]
Fn_CO2 <- estandarizar(cbind(Temp,CO2))[,2]
U_limit <- seq(0,1,by = 0.2)
png(paste(SM_dir,"/CO2_Temp_scatterplots.png",sep=""),   bg = "white", width = 1200, height = 1000, res = 120)
par(mfrow = c(2,2),mar = c(4, 6, 2, 2), cex.lab= 1.2, cex.axis = 1.2)
#-------------------Figure a----------------#
plot(Temp, CO2, pch = 1, xlab= Templabel, ylab = CO2label,  
     xlim = c(min(Temp), max(Temp)), 
     ylim = c(min(CO2), max(CO2)))
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Temp, CO2, pch = 1, xlab= Templabel, ylab = CO2label, 
     xlim = c(min(Temp), max(Temp)), 
     ylim = c(min(CO2), max(CO2)))
box()

#----------------- Figure b----------------#
plot(Temp, CO2, pch = 1, xlab= Templabel, ylab = CO2label,  
     xlim = c(min(Temp), max(Temp)), col = "lightgray",
     ylim = c(min(CO2), max(CO2)))
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Temp, CO2, pch = 1, xlab= Templabel, ylab = CO2label, 
     xlim = c(min(Temp), max(Temp)),  col = "lightgray",
     ylim = c(min(CO2), max(CO2)))
par(new=TRUE)
plot(Temp_fixed_sub48, CO2_fixed_sub48,  xlab= Templabel, ylab = CO2label,
     pch = 15, col = "green", cex = 1.2,
     xlim = c(min(Temp), max(Temp)),  ylim = c(min(CO2), max(CO2)))
box()

#----------------- Figure c----------------#
plot(Temp, CO2, pch = 1, xlab= Templabel, ylab = CO2label,  
     xlim = c(min(Temp), max(Temp)), col = "lightgray",
     ylim = c(min(CO2), max(CO2)))
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Temp, CO2, pch = 1, xlab= Templabel, ylab = CO2label, 
     xlim = c(min(Temp), max(Temp)), col = "lightgray",
     ylim = c(min(CO2), max(CO2)))
par(new=TRUE)
plot(Temp_clhs_sub48, CO2_clhs_sub48,  xlab= Templabel, ylab = CO2label,
     pch = 17, col = "red", cex = 1.2,
     xlim = c(min(Temp), max(Temp)),  ylim = c(min(CO2), max(CO2)))
box()

#----------------- Figure d----------------#
plot(Temp, CO2, pch = 1, xlab= Templabel, ylab = CO2label,  
     xlim = c(min(Temp), max(Temp)), col = "lightgray",
     ylim = c(min(CO2), max(CO2)))
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Temp, CO2, pch = 1, xlab= Templabel, ylab = CO2label, 
     xlim = c(min(Temp), max(Temp)), col = "lightgray",
     ylim = c(min(CO2), max(CO2)))
par(new=TRUE)
plot(Temp_aclhs_sub48, CO2_aclhs_sub48,  xlab= Templabel, ylab = CO2label,
     pch = 18, col = "blue", cex = 1.5,
     xlim = c(min(Temp), max(Temp)),  ylim = c(min(CO2), max(CO2)))
box()
dev.off()


#------------------------------- Figure 7 Variograms--------------------------------------#
Variolabel <- expression(bold("Semivariance"))
N_lags<- 1 
lag_value <- 8 
CO2_fixed_sub48_VarioEstimation<-Variograma(Time_fixed_sub48, numeric(length(Time_fixed_sub48)), 
                                            CO2_fixed_sub48, 0, 90, N_lags, lag_value, 1, "", Timelabel) # CO2 temporal variogram
N_lags<- 8 
lag_value <- 1 
CO2_clhs_sub48_VarioEstimation<-Variograma(Time_clhs_sub48, numeric(length(Time_clhs_sub48)), 
                                            CO2_clhs_sub48, 0, 90, N_lags, lag_value, 1, "", Timelabel) # CO2 temporal variogram
CO2_aclhs_sub48_VarioEstimation<-Variograma(Time_aclhs_sub48, numeric(length(Time_aclhs_sub48)), 
                                            CO2_aclhs_sub48, 0, 90, N_lags, lag_value, 1, "", Timelabel) # CO2 temporal variogram

#Variogram models (1- exponential, 2- spherical, 3- gaussian)
N_lags<- 8
lag_value <- 1 
CO2_aclhs_sub48_vario_model<- 2 
CO2_aclhs_sub48_nugget<- 0.57 
CO2_aclhs_sub48_sill_and_nugget<- 1.95 
CO2_aclhs_sub48_range <- 6.21
png(paste(SM_dir,"/Variogram_comparison.png",sep=""), bg = "white",  width = 1500, height = 1000, res = 200)
par(mfrow = c(1,1),mar = c(4, 6, 2, 2), cex.lab= 1.2, cex.axis = 1.2)
plot(CO2_VarioEstimation[,c(2,3)], pch = 1, xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*1.7)),
     xlab = Timelabel ,ylab = Variolabel)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(CO2_VarioEstimation[,c(2,3)], pch = 1, xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*1.7)),
     xlab = Timelabel ,ylab = Variolabel )
par(new=TRUE)
plot(CO2_fixed_sub48_VarioEstimation[,c(2,3)], pch = 15, cex = 1.2, col = "green", xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*1.7)),
     xlab = Timelabel ,ylab = Variolabel)
par(new=TRUE)
plot(CO2_clhs_sub48_VarioEstimation[,c(2,3)], pch = 17, cex = 1.2, col = "red",  xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*1.7)),
     xlab = Timelabel , ylab = Variolabel)
par(new=TRUE)
plot(CO2_aclhs_sub48_VarioEstimation[,c(2,3)], pch = 18, cex = 1.5, col = "blue", xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*1.7)),
     xlab = Timelabel , ylab = Variolabel)
par(new=TRUE)
lines.variomodel(cov.model = "sph", cov.pars = c(CO2_aclhs_sub48_sill_and_nugget-CO2_aclhs_sub48_nugget, CO2_aclhs_sub48_range)
                 , nug = CO2_aclhs_sub48_nugget, col = "Black", lwd = 3, max.dist = 1*8)
box()
dev.off()


N_lags<- 8
lag_value <- 1 
CO2_aclhs_sub48_vario_model<- 2 
CO2_aclhs_sub48_nugget<- 0.57 
CO2_aclhs_sub48_sill_and_nugget<- 1.95 
CO2_aclhs_sub48_range <- 6.21
png(paste(SM_dir,"/Variogram_comparison_.png",sep=""), bg = "white",  width = 1500, height = 1000, res = 200)
par(mfrow = c(1,1),mar = c(4, 6, 2, 2), cex.lab= 1.2, cex.axis = 1.2)
plot(CO2_VarioEstimation[,c(2,3)], pch = 1, xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*1.2)),
     xlab = Timelabel ,ylab = Variolabel)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(CO2_VarioEstimation[,c(2,3)], pch = 1, xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*1.2)),
     xlab = Timelabel ,ylab = Variolabel )
par(new=TRUE)
plot(CO2_aclhs_sub48_VarioEstimation[,c(2,3)], pch = 18, cex = 1.5, col = "blue", xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*1.2)),
     xlab = Timelabel , ylab = Variolabel)
par(new=TRUE)
lines.variomodel(cov.model = "sph", cov.pars = c(CO2_aclhs_sub48_sill_and_nugget-CO2_aclhs_sub48_nugget, CO2_aclhs_sub48_range)
                 , nug = CO2_aclhs_sub48_nugget, col = "Black", lwd = 3, max.dist = 1*8)
box()
dev.off()

#----------------------Tables------------------#
# Statistical Properties
cbind(CO2_Stat,CO2_fixed_sub48_Stat,CO2_clhs_sub48_Stat,CO2_aclhs_sub48_Stat)
cbind(Temp_Stat,Temp_fixed_sub48_Stat,Temp_clhs_sub48_Stat,Temp_aclhs_sub48_Stat)

# Kolmogorov Smirnov test
ks_test_fixed_CO2 <- ks.test(CO2_fixed_sub48, CO2, conf.level = 0.95)
ks_test_clhs_CO2 <- ks.test(CO2_clhs_sub48, CO2, conf.level = 0.95)
ks_test_aclhs_CO2 <- ks.test(CO2_aclhs_sub48, CO2, conf.level = 0.95)
ks_test_fixed_CO2
ks_test_clhs_CO2
ks_test_aclhs_CO2
ks_test_fixed_Temp <- ks.test(Temp_fixed_sub48, Temp, conf.level = 0.95)
ks_test_clhs_Temp <- ks.test(Temp_clhs_sub48, Temp, conf.level = 0.95)
ks_test_aclhs_Temp <- ks.test(Temp_aclhs_sub48, Temp, conf.level = 0.95)
ks_test_fixed_Temp
ks_test_clhs_Temp
ks_test_aclhs_Temp

# Data
(corP_Temp_CO2 <- round(cor(Temp,CO2, method = "pearson"),3))
(corS_Temp_CO2 <- round(cor(Temp,CO2, method = "spearman"),3))
(corK_Temp_CO2 <- round(cor(Temp,CO2, method = "kendall"),3))
# fixed
(corP_Temp_CO2_fixed_sub48 <- round(cor(Temp_fixed_sub48,CO2_fixed_sub48, method = "pearson"),3))
(corS_Temp_CO2_fixed_sub48 <- round(cor(Temp_fixed_sub48,CO2_fixed_sub48, method = "spearman"),3))
(corK_Temp_CO2_fixed_sub48 <- round(cor(Temp_fixed_sub48,CO2_fixed_sub48, method = "kendall"),3))
# clhs
(corP_Temp_CO2_clhs_sub48 <- round(cor(Temp_clhs_sub48,CO2_clhs_sub48, method = "pearson"),3))
(corS_Temp_CO2_clhs_sub48 <- round(cor(Temp_clhs_sub48,CO2_clhs_sub48, method = "spearman"),3))
(corK_Temp_CO2_clhs_sub48 <- round(cor(Temp_clhs_sub48,CO2_clhs_sub48, method = "kendall"),3))
#aclhs
(corP_Temp_CO2_aclhs_sub48 <- round(cor(Temp_aclhs_sub48,CO2_aclhs_sub48, method = "pearson"),3))
(corS_Temp_CO2_aclhs_sub48 <- round(cor(Temp_aclhs_sub48,CO2_aclhs_sub48, method = "spearman"),3))
(corK_Temp_CO2_aclhs_sub48 <- round(cor(Temp_aclhs_sub48,CO2_aclhs_sub48, method = "kendall"),3))

round(abs(corP_Temp_CO2_fixed_sub48-corP_Temp_CO2),3)
round(abs(corS_Temp_CO2_fixed_sub48-corS_Temp_CO2),3)
round(abs(corK_Temp_CO2_fixed_sub48-corK_Temp_CO2),3)
round(abs(corP_Temp_CO2_clhs_sub48-corP_Temp_CO2),3)
round(abs(corS_Temp_CO2_clhs_sub48-corS_Temp_CO2),3)
round(abs(corK_Temp_CO2_clhs_sub48-corK_Temp_CO2),3)
round(abs(corP_Temp_CO2_aclhs_sub48-corP_Temp_CO2),3)
round(abs(corS_Temp_CO2_aclhs_sub48-corS_Temp_CO2),3)
round(abs(corK_Temp_CO2_aclhs_sub48-corK_Temp_CO2),3)

#
sum(abs(CO2_clhs_sub48_VarioEstimation[,3]-CO2_VarioEstimation[,3]))
sum(abs(CO2_aclhs_sub48_VarioEstimation[,3]-CO2_VarioEstimation[,3]))
