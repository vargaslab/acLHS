# Creates a folder to store results for Sampling Methods (SM)
dir.create(paste(getwd(),"/Results/spatial_2D/SM", sep=""))
SM_dir<-paste(root_dir,"/Results/spatial_2D/SM",sep="")
#----------------------------------Fixed Temporal Sampling-------------------------------------------#
Data_Frame_gridded <-Data_Frame
gridded(Data_Frame_gridded) = ~X+Y
image(Data_Frame_gridded, attr=2 )
#fixed_sub50  <- spsample(Data_Frame_gridded,n=50,type= "regular") 
# type="random", "stratified", "nonaligned", "regular", "Fibonacci" , "hexagonal"
fixed_sub50
str(fixed_sub50)
points(fixed_sub50, pch=3, cex=.5) 
# str(fixed_sub50)
Xfixed_sub50 <- fixed_sub50@coords[,1]
Yfixed_sub50 <- fixed_sub50@coords[,2]
pos_sub_fixed <- NULL
for (i in 1:length(Xfixed_sub50)) {
  pos_sub_fixed0 <- which.min(abs(X-Xfixed_sub50[i])+abs(Y-Yfixed_sub50[i]))
  pos_sub_fixed <- rbind(pos_sub_fixed, pos_sub_fixed0)
}
X_fixed_sub50 <- X[pos_sub_fixed]
Y_fixed_sub50 <- Y[pos_sub_fixed]
Temp_fixed_sub50 <- Temp[pos_sub_fixed]
Rs_fixed_sub50 <- Rs[pos_sub_fixed]
plot(Y,X)
plot(X,Y, xlim = c(min(X),max(X)), ylim = c(min(Y),max(Y)), pch = 19)
par(new=TRUE)
plot(X_fixed_sub50,Y_fixed_sub50, col= "red",xlim = c(min(X),max(X)), ylim = c(min(Y),max(Y)), pch = 19)
# Stats
Rs_fixed_sub50_Stat <- Estadisticas(Rs_fixed_sub50)
Temp_fixed_sub50_Stat <- Estadisticas(Temp_fixed_sub50)
#plot(Data_Frame[1:nrow(Data_Frame),c(1,2)], xlim = c(min(X),max(X)), ylim = c(min(Y),max(Y)))
# plot(Data_Frame[,c(1,2)],pch = 19)
# points(spsample(SpatialPoints(Data_Frame[,c(1,2)]), n = 100, "stratified"),pch = 19, col = "red")

# fixed_Sample <- cbind(X_fixed_sub50,Y_fixed_sub50,Temp_fixed_sub50,Rs_fixed_sub50)
# write.csv(fixed_Sample , file = paste(SM_dir,"/fixed_Sample.csv",sep=""))

#-----------------------------------------cLHS-------------------------------------------------------#
n_sub = 50
df = data.frame(Temp,Rs)
# Returning the indices of the sampled points
# res <- clhs(df, CoordXY = cbind(X,ti), size = n_sub, use.cpp = FALSE, progress = FALSE, simple = TRUE)
# str(res)
# Returning a clhs_result object for plotting using C++
system.time(res_ <- clhs(df, size = n_sub, use.cpp = TRUE, iter = 100000, progress = TRUE, simple = FALSE))
# user  system elapsed 
#33.478   0.000  33.618 
res_
#str(res_)
plot(res_)
#clhs_sub50 <- res_$index_samples
#clhs_sub50_ <- read.csv(file=file.choose(),header=T,na.strings="-999.25")
#clhs_sub50  <- clhs_sub50_$clhs_sub50
X_clhs_sub50 <- X[clhs_sub50]
Y_clhs_sub50 <- Y[clhs_sub50]
Temp_clhs_sub50 <- Temp[clhs_sub50]
Rs_clhs_sub50 <- Rs[clhs_sub50]
# Stats
Rs_clhs_sub50_Stat <- Estadisticas(Rs_clhs_sub50)
Temp_clhs_sub50_Stat <- Estadisticas(Temp_clhs_sub50)
plot(X,Y, xlim = c(min(X),max(X)), ylim = c(min(Y),max(Y)), pch = 19)
par(new=TRUE)
plot(X_clhs_sub50,Y_clhs_sub50, col= "red",xlim = c(min(X),max(X)), ylim = c(min(Y),max(Y)), pch = 19)

# clhs_Sample <- cbind(X_clhs_sub50,Y_clhs_sub50,Temp_clhs_sub50,Rs_clhs_sub50)
# write.csv(clhs_Sample , file = paste(SM_dir,"/clhs_Sample.csv",sep=""))

#-------------------------------------------------acLHS---------------------------------------------#
n_sub = 50
Fn_limit <- seq(0,1,length.out = n_sub+1)
Temp_limit <- quantile(Temp,Fn_limit)
Rs_limit <- quantile(Rs,Fn_limit)
X_limit <- quantile(X,Fn_limit)
Y_limit <- quantile(Y,Fn_limit)
# OF1
FuncO1  <- function(var_sub) {
  pos_sub <- NULL
  for (i in 1:n_sub) {
    pos_sub <- rbind(pos_sub,which.min(abs(Rs-var_sub[i])))
  }
  Temp_sub <- Temp[pos_sub]
  hist_Temp_sub <- hist(Temp_sub, breaks = Temp_limit, plot = FALSE)$counts
  frecuencia <- rep(1,n_sub)
  # X_sub <- X[pos_sub]
  # hist_X_sub <- hist(X_sub, breaks = X_limit, plot = FALSE)$counts
  # Y_sub <- Y[pos_sub]
  # hist_Y_sub <- hist(Y_sub, breaks = Y_limit, plot = FALSE)$counts
  return (sum(abs(hist_Temp_sub-frecuencia))) # + sum(abs(hist_X_sub-frecuencia)) + sum(abs(hist_Y_sub-frecuencia)) ) # initial = 8 
}                #0                               #16                             #48

# min = 64 for X,Y,Temp


corP_Temp_Rs <-  cor(Temp,Rs, method = "pearson")
corS_Temp_Rs <- cor(Temp,Rs, method = "spearman")
corK_Temp_Rs <- cor(Temp,Rs, method = "kendall")
# OF2
FunRs  <- function(var_sub) {
  pos_sub <- NULL
  for (i in 1:n_sub) {
    pos_sub <- rbind(pos_sub,which.min(abs(Rs-var_sub[i])))
  }
  corP_Temp_Vsub <- cor(Temp[pos_sub],Rs[pos_sub], method = "pearson")
  corS_Temp_Vsub <- cor(Temp[pos_sub],Rs[pos_sub], method = "spearman")
  corK_Temp_Vsub <- cor(Temp[pos_sub],Rs[pos_sub], method = "kendall")
  return (abs(corP_Temp_Vsub-corP_Temp_Rs) + abs(corS_Temp_Vsub-corS_Temp_Rs) + abs(corK_Temp_Vsub-corK_Temp_Rs)) 
}

#OF3
FuncO3 <- function(var_sub) {
  N_lags<- 10
  lag_value <- DistMin
  pos_sub <- NULL
  for (i in 1:n_sub) {
    pos_sub <- rbind(pos_sub,which.min(abs(Rs-var_sub[i])))
  }
  X_sub <- X[pos_sub]
  Y_sub <- Y[pos_sub]
  Vario_var_sub <- variog(as.geodata(cbind(X_sub, Y_sub, var_sub), 
                                     coords.col = 1:2, data.col = 3),  breaks = c(seq(lag_value/2, lag_value * (N_lags+1), lag_value)),
                          trend = "cte", lambda = 1, estimator.type = "classical", nugget.tolerance = 0, direction = 0, 
                          tolerance = 90, unit.angle = "degrees", pairs.min = 1)
  return (sum(abs(Vario_var_sub$v-Rs_VarioEstimation[,3]))) # # initial = 11.20084
}

# OF = w1 * OF1 + w2 * OF2 + w3 * OF3
w1 = 10; w2 = 1000 ; w3 = 0.001 #0.001 #w1 = 1; w2 = 100; w3 = 0.00012
FuncConj  <- function(var_sub) {
  total =  w1*FuncO1(var_sub) + w2*FunRs(var_sub) + w3*FuncO3(var_sub) 
  return (total)
}

# optimization
Fn_opt_sub <- seq(1/n_sub,1,length.out = n_sub)
Rs_optimal_sub <- unname(quantile(Rs,Fn_opt_sub)) 
Rs_order <- Rs[order(Rs)]
pos_otp_order <- NULL
for (i in 1:n_sub) {
  pos_otp_order <- rbind(pos_otp_order,which.min(abs(Rs_order-Rs_optimal_sub[i])))
}
tol_sub <-  round(length(Rs)/n_sub)
Pop_Inicial <- matrix(data=NA,nrow=tol_sub,ncol=n_sub)
for (i in 1:(n_sub)) {
  for (j in 1:tol_sub) {
    pos <- pos_otp_order[i]
    Pop_Inicial[j,i] <- Rs_order[pos-j]
  }
}
Fn_sub_limit <- seq(0,1,1/n_sub)
Rs_sub_limit <- unname(quantile(Rs, probs = Fn_sub_limit))
subvector_min = Rs_sub_limit[-(n_sub+1)]
subvector_max = Rs_sub_limit[-(1)]
system.time(aclhs_sub_outDEoptim_Rs <- DEoptim(fn=FuncConj, lower= subvector_min,
                                                  upper= subvector_max,
                                                  control = DEoptim.control(VTR = 0.000001,strategy = 3, 
                                                                            itermax =20000, reltol = 1e-8,
                                                                            CR = 0.5, F = 0.8, NP= nrow(Pop_Inicial),
                                                                            initialpop = Pop_Inicial))) #, 

# user  system elapsed
# 51.284   2.818  53.845    itermax =100,  n_sub = 50, 3functions=1+100+0.0001,  aclhs_sub_outDEoptim_Rs$optim$bestval = 26.07795
# 477.530  43.656 518.384   itermax =1000,  n_sub = 50, 3functions=1+100+0.0001,  aclhs_sub_outDEoptim_Rs$optim$bestval = 9.410707
# 990.644   90.114 1092.747 itermax =2000,  n_sub = 50, 3functions=1+100+0.0001,  aclhs_sub_outDEoptim_Rs$optim$bestval = 7.753369
# 1343.538  122.207 1431.503 itermax =3000,  n_sub = 50, 3functions=1+100+0.0001,  aclhs_sub_outDEoptim_Rs$optim$bestval = 5.076521
# 2245.971  202.853 2402.099  itermax =5000,  n_sub = 50, 3functions=1+100+0.0001,  aclhs_sub_outDEoptim_Rs$optim$bestval = 2.149056
# 1415.563  115.450 1505.200  itermax =3000,  n_sub = 50, 3functions=1+100+0.0001,  aclhs_sub_outDEoptim_Rs$optim$bestval = 2.164424
# 1401.718  120.080 1498.022  itermax =3000,  n_sub = 50, 3functions=1+100+0.0001,  aclhs_sub_outDEoptim_Rs$optim$bestval = 5.086694
# 2358.438  201.183 2509.314  itermax =5000,  n_sub = 50, 3functions=1+100+0.0001,  aclhs_sub_outDEoptim_Rs$optim$bestval = 4.524411
# 2569.010  215.139 2781.487  itermax =5000,  n_sub = 50, 3functions=1+100+0.001,  aclhs_sub_outDEoptim_Rs$optim$bestval = 24.14481
#--------------------------
# 894.832  98.565 977.312    itermax =10000,  n_sub = 50, 3functions=1+100+0.0001,  aclhs_sub_outDEoptim_Rs$optim$bestval = 2.682616
# 950.334  102.787 1037.541 itermax =10000,  n_sub = 50, 3functions=1+100+0.001,  aclhs_sub_outDEoptim_Rs$optim$bestval = 25.54513
# 953.261  101.977 1039.385  itermax =10000,  n_sub = 50, 3functions=1+100+0.0005,  aclhs_sub_outDEoptim_Rs$optim$bestval = 12.24948
# 583.380  60.488 675.809  itermax =5000,  n_sub = 50, 3functions=1+100+0.0003,  aclhs_sub_outDEoptim_Rs$optim$bestval = 12.29842
# 1139.254  117.554 1260.909 itermax =10000,  n_sub = 50, 3functions=1+100+0.0003,  aclhs_sub_outDEoptim_Rs$optim$bestval = 11.04133
# 909.259   99.741 1003.150  itermax =10000,  n_sub = 50, 3functions=1+100+0.0001,  aclhs_sub_outDEoptim_Rs$optim$bestval =  3.169948
# 1024.458  139.425 1109.602 itermax =20000,  n_sub = 50, 3functions=1+100+0.0001,  aclhs_sub_outDEoptim_Rs$optim$bestval =  3.686642
# 991.216  137.179 1072.770  itermax =20000,  n_sub = 50, 3functions=1+100+0.0002,  aclhs_sub_outDEoptim_Rs$optim$bestval =  7.831948
# 1234.100  160.064 1327.238 itermax =20000,  n_sub = 50, 3functions=1+100+0.0002,  aclhs_sub_outDEoptim_Rs$optim$bestval =  75.81003
#1259.091  162.862 1354.630  itermax =20000,  n_sub = 50, 3functions=1+100+0.0001,  aclhs_sub_outDEoptim_Rs$optim$bestval =  73.59492
# 634.308  81.747 683.351    itermax =10000,  n_sub = 50, 3functions=1+100+0.0001,  aclhs_sub_outDEoptim_Rs$optim$bestval =  13.0525
# 1266.020  163.959 1364.970 itermax =10000,  n_sub = 50, 3functions=1+100+0.0001,  aclhs_sub_outDEoptim_Rs$optim$bestval =  13.54867
# 499.485  68.481 540.304    itermax =10000,  n_sub = 50, 3functions=1+100+0.0001,  aclhs_sub_outDEoptim_Rs$optim$bestval =  3.89086
# 1022.832  138.820 1106.150 itermax =20000,  n_sub = 50, 3functions=1+100+0.00005,  aclhs_sub_outDEoptim_Rs$optim$bestval =  3.228897
#  996.394  135.832 1077.084 itermax =20000,  n_sub = 50, 3functions=10+1000+0.00001,  aclhs_sub_outDEoptim_Rs$optim$bestval =  0.8394497
# 1016.358  137.473 1099.743 itermax =20000,  n_sub = 50, 3functions=10+1000+0.001,  aclhs_sub_outDEoptim_Rs$optim$bestval =  44.32176
aclhs_sub_outDEoptim_Rs$optim$bestval 
Ite <- seq(1,length(aclhs_sub_outDEoptim_Rs$member$bestvalit),1)
plot(Ite,aclhs_sub_outDEoptim_Rs$member$bestvalit, "l", xlab ="Iterations", ylab ="Objective Function")
var_sub <- as.numeric(aclhs_sub_outDEoptim_Rs$optim$bestmem) 

#Estadisticas(var_sub)
pos_sub_Rs <- NULL
for (i in 1:n_sub) {
  pos_sub_Rs <- rbind(pos_sub_Rs,which.min(abs(Rs-var_sub[i])))
}

aclhs_sub50 <- pos_sub_Rs
X_aclhs_sub50 <- X[aclhs_sub50]
Y_aclhs_sub50 <- Y[aclhs_sub50]
Temp_aclhs_sub50 <- Temp[aclhs_sub50]
Rs_aclhs_sub50 <- Rs[aclhs_sub50]
# Stats
#Rs_aclhs_sub50_Stat <- Estadisticas(Rs_aclhs_sub50)
#Temp_aclhs_sub50_Stat <- Estadisticas(Temp_aclhs_sub50)

plot(X,Y, xlim = c(min(X),max(X)), ylim = c(min(Y),max(Y)), pch = 19)
par(new=TRUE)
plot(X_aclhs_sub50,Y_aclhs_sub50, xlab = "X", ylab = "Y",col= "red",xlim = c(min(X),max(X)), ylim = c(min(Y),max(Y)), pch = 19)

# aclhs_Sample <- cbind(X_aclhs_sub50,Y_aclhs_sub50,Temp_aclhs_sub50,Rs_aclhs_sub50)
# write.csv(aclhs_Sample , file = paste(SM_dir,"/aclhs_Sample.csv",sep=""))

#--------------------Figure 2------------------------------#
png(paste(SM_dir,"/2_Process_optimization_Temp_Rs.png",sep=""), bg = "white",  width = 1500, height = 1500, res = 250)
plot(Ite,aclhs_sub_outDEoptim_Rs$member$bestvalit,
     "l", xlab = expression(bold("Iterations")), ylab = expression(bold("Objective Function")))
dev.off()

png(paste(SM_dir,"/2_Temp_Rs_Spatial_Sampling_Temp_Rs_aclhs.png",sep=""),
    bg = "white", width = 1500, height = 1000, res = 250)
par(mfrow = c(1,1),mar = c(4, 6, 2, 2), cex.lab= 1.0, cex.axis = 1.0)
plot(X,Y, xlab = Xlabel, ylab = Ylabel, xlim = c(min(X),max(X)), ylim = c(min(Y),max(Y)), pch = 0)
par(new=TRUE)
plot(X_aclhs_sub50,Y_aclhs_sub50, pch = 18, col= "blue", xlab = Xlabel, ylab = Ylabel, xlim = c(min(X),max(X)), ylim = c(min(Y),max(Y)))
box()
dev.off()

#--------------------------------Figure 11--------------------------------------#
png(paste(SM_dir,"/2_Temp_Rs_Spatial_Sampling_Temp_Rs.png",sep=""),
    bg = "white", width = 1500, height = 1000, res = 250)
par(mfrow = c(1,1),mar = c(4, 6, 2, 2), cex.lab= 1.0, cex.axis = 1.0)
plot(X,Y, xlab = Xlabel, ylab = Ylabel, xlim = c(min(X),max(X)), ylim = c(min(Y),max(Y)), pch = 0)
par(new=TRUE)
plot(X_fixed_sub50,Y_fixed_sub50, pch = 15, col= "green", xlab = Xlabel, ylab = Ylabel, xlim = c(min(X),max(X)), ylim = c(min(Y),max(Y)))
par(new=TRUE)
plot(X_clhs_sub50,Y_clhs_sub50, pch = 17, col= "red", xlab = Xlabel, ylab = Ylabel, xlim = c(min(X),max(X)), ylim = c(min(Y),max(Y)))
par(new=TRUE)
plot(X_aclhs_sub50,Y_aclhs_sub50, pch = 18, col= "blue", xlab = Xlabel, ylab = Ylabel, xlim = c(min(X),max(X)), ylim = c(min(Y),max(Y)))
box()
dev.off()


png(paste(SM_dir,"/2_Rs_Spatial_Dist_Temp_Rs.png",sep=""),
    bg = "white", width = 1500, height = 1000, res = 250)
par(mfrow = c(1,1),mar = c(4, 6, 2, 2), cex.lab= 1.0, cex.axis = 1.0)
plot(X,Y, xlab = Xlabel, ylab = Ylabel, xlim = c(min(X),max(X)), ylim = c(min(Y),max(Y)), pch = 15)
dev.off()


#---------------------------------------------Figure4_Univariate pdfs--------------------------------#
#Fn_Temp <- estandarizar(cbind(Temp,CO2))[,1]
#Fn_CO2 <- estandarizar(cbind(Temp,CO2))[,2]
U_limit <- seq(0,1,by = 0.2)
png(paste(SM_dir,"/2_Temp_univariate_pdfs.png",sep=""),   bg = "white", width = 1200, height = 1000, res = 120)
par(mfrow = c(2,2),mar = c(4, 6, 2, 2), cex.lab= 1.2, cex.axis = 1.2)
#-------------------Figure a----------------#
plot(ecdf(Temp), pch = 0, yaxt = "n",  main = " ", xlab= Templabel, 
     ylab = expression(bold(paste(Fn(Temperature)))),xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
axis(2, at=U_limit)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ecdf(Temp), pch = 0, yaxt = "n",  main = " ", xlab= Templabel, 
     ylab = expression(bold(paste(Fn(Temperature)))),xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
axis(2, at=U_limit)
box()

#-------------------Figure b----------------#
plot(ecdf(Temp), pch = 0, yaxt = "n",  main = " ", xlab= Templabel,  col = "lightgray",
     ylab = expression(bold(paste(Fn(Temperature)))),xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
axis(2, at=U_limit)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ecdf(Temp), pch = 0, yaxt = "n",  main = " ", xlab= Templabel, col = "lightgray",
     ylab = expression(bold(paste(Fn(Temperature)))),xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
axis(2, at=U_limit)
par(new=TRUE)
plot(ecdf(Temp_fixed_sub50),pch = 15, col = "green", cex = 1.2,
     yaxt = "n",  main = " ", xlab= Templabel, ylab = expression(bold(paste(Fn(Temperature)))),
     xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
box()

#-------------------Figure c----------------#
plot(ecdf(Temp), pch = 0, yaxt = "n",  main = " ", xlab= Templabel, col = "lightgray",
     ylab = expression(bold(paste(Fn(Temperature)))),xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
axis(2, at=U_limit)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ecdf(Temp), pch = 0, yaxt = "n",  main = " ", xlab= Templabel, col = "lightgray",
     ylab = expression(bold(paste(Fn(Temperature)))),xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
axis(2, at=U_limit)
par(new=TRUE)
plot(ecdf(Temp_clhs_sub50),pch = 17, col = "red", cex = 1.2,
     yaxt = "n",  main = " ", xlab= Templabel, ylab = expression(bold(paste(Fn(Temperature)))),
     xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
box()

#-------------------Figure d----------------#
plot(ecdf(Temp), pch = 0, yaxt = "n",  main = " ", xlab= Templabel, col = "lightgray",
     ylab = expression(bold(paste(Fn(Temperature)))),xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
axis(2, at=U_limit)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ecdf(Temp), pch = 0, yaxt = "n",  main = " ", xlab= Templabel, col = "lightgray",
     ylab = expression(bold(paste(Fn(Temperature)))),xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
axis(2, at=U_limit)
par(new=TRUE)
plot(ecdf(Temp_aclhs_sub50),pch = 18, col = "blue", cex = 1.5,
     yaxt = "n",  main = " ", xlab= Templabel, ylab = expression(bold(paste(Fn(Temperature)))),
     xlim =c(min(Temp),max(Temp)), ylim = c(0,1))
box()
dev.off()


#---------------------------------------------Figure5_Univariate pdfs--------------------------------#
#Fn_Temp <- estandarizar(cbind(Temp,Rs))[,1]
#Fn_Rs <- estandarizar(cbind(Temp,Rs))[,2]
U_limit <- seq(0,1,by = 0.2)
png(paste(SM_dir,"/2_Rs_univariate_pdfs.png",sep=""),   bg = "white", width = 1200, height = 1000, res = 120)
par(mfrow = c(2,2),mar = c(4, 6, 2, 2), cex.lab= 1.2, cex.axis = 1.2)
#-------------------Figure a----------------#
plot(ecdf(Rs), pch = 0, yaxt = "n",  main = " ", xlab= Rslabel, 
     ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),xlim =c(min(Rs),max(Rs)), ylim = c(0,1))
axis(2, at=U_limit)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ecdf(Rs), pch = 0, yaxt = "n",  main = " ", xlab= Rslabel, 
     ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),xlim =c(min(Rs),max(Rs)), ylim = c(0,1))
axis(2, at=U_limit)
box()

#-------------------Figure b----------------#
plot(ecdf(Rs), pch = 0, yaxt = "n",  main = " ", xlab= Rslabel, col = "lightgray",
     ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),xlim =c(min(Rs),max(Rs)), ylim = c(0,1))
axis(2, at=U_limit)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ecdf(Rs), pch = 0, yaxt = "n",  main = " ", xlab= Rslabel, col = "lightgray",
     ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),xlim =c(min(Rs),max(Rs)), ylim = c(0,1))
axis(2, at=U_limit)
par(new=TRUE)
plot(ecdf(Rs_fixed_sub50),pch = 15, col = "green", cex = 1.2,
     yaxt = "n",  main = " ", xlab= Rslabel, ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),
     xlim =c(min(Rs),max(Rs)), ylim = c(0,1))
box()

#-------------------Figure c----------------#
plot(ecdf(Rs), pch = 0, yaxt = "n",  main = " ", xlab= Rslabel, col = "lightgray",
     ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),xlim =c(min(Rs),max(Rs)), ylim = c(0,1))
axis(2, at=U_limit)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ecdf(Rs), pch = 0, yaxt = "n",  main = " ", xlab= Rslabel, col = "lightgray",
     ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),xlim =c(min(Rs),max(Rs)), ylim = c(0,1))
axis(2, at=U_limit)
par(new=TRUE)
plot(ecdf(Rs_clhs_sub50),pch = 17, col = "red", cex = 1.2,
     yaxt = "n",  main = " ", xlab= Rslabel, ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),
     xlim =c(min(Rs),max(Rs)), ylim = c(0,1))
box()

#-------------------Figure d----------------#
plot(ecdf(Rs), pch = 0, yaxt = "n",  main = " ", xlab= Rslabel, col = "lightgray",
     ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),xlim =c(min(Rs),max(Rs)), ylim = c(0,1))
axis(2, at=U_limit)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(ecdf(Rs), pch = 0, yaxt = "n",  main = " ", xlab= Rslabel, col = "lightgray",
     ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),xlim =c(min(Rs),max(Rs)), ylim = c(0,1))
axis(2, at=U_limit)
par(new=TRUE)
plot(ecdf(Rs_aclhs_sub50),pch = 18, col = "blue", cex = 1.5,
     yaxt = "n",  main = " ", xlab= Rslabel, ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),
     xlim =c(min(Rs),max(Rs)), ylim = c(0,1))
box()
dev.off()


#---------------------------------------------Figure6_scatterplots--------------------------------#
Fn_Temp <- estandarizar(cbind(Temp,Rs))[,1]
Fn_Rs <- estandarizar(cbind(Temp,Rs))[,2]
U_limit <- seq(0,1,by = 0.2)
png(paste(SM_dir,"/2_Rs_Temp_scatterplots.png",sep=""),   bg = "white", width = 1200, height = 1000, res = 120)
par(mfrow = c(2,2),mar = c(4, 6, 2, 2), cex.lab= 1.2, cex.axis = 1.2)
#-------------------Figure a----------------#
plot(Temp, Rs, pch = 0, xlab= Templabel, ylab = Rslabel,  
     xlim = c(min(Temp), max(Temp)), 
     ylim = c(min(Rs), max(Rs)))
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Temp, Rs, pch = 0, xlab= Templabel, ylab = Rslabel, 
     xlim = c(min(Temp), max(Temp)), 
     ylim = c(min(Rs), max(Rs)))
box()

#----------------- Figure b----------------#
plot(Temp, Rs, pch = 0, xlab= Templabel, ylab = Rslabel,  
     xlim = c(min(Temp), max(Temp)), col = "lightgray",
     ylim = c(min(Rs), max(Rs)))
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Temp, Rs, pch = 0, xlab= Templabel, ylab = Rslabel, 
     xlim = c(min(Temp), max(Temp)),  col = "lightgray",
     ylim = c(min(Rs), max(Rs)))
# # abline(v=Temp_limit, col= 'green', lty =2, lwd=1)
# # abline(h=Rs_limit, col= 'green', lty =2, lwd=1)
par(new=TRUE)
plot(Temp_fixed_sub50, Rs_fixed_sub50,  xlab= Templabel, ylab = Rslabel,
     pch = 15, col = "green", cex = 1.2,
     xlim = c(min(Temp), max(Temp)),  ylim = c(min(Rs), max(Rs)))
box()

#----------------- Figure c----------------#
plot(Temp, Rs, pch = 0, xlab= Templabel, ylab = Rslabel,  
     xlim = c(min(Temp), max(Temp)), col = "lightgray",
     ylim = c(min(Rs), max(Rs)))
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Temp, Rs, pch = 0, xlab= Templabel, ylab = Rslabel, 
     xlim = c(min(Temp), max(Temp)), col = "lightgray",
     ylim = c(min(Rs), max(Rs)))
par(new=TRUE)
plot(Temp_clhs_sub50, Rs_clhs_sub50,  xlab= Templabel, ylab = Rslabel,
     pch = 17, col = "red", cex = 1.2,
     xlim = c(min(Temp), max(Temp)),  ylim = c(min(Rs), max(Rs)))
box()

#----------------- Figure d----------------#
plot(Temp, Rs, pch = 0, xlab= Templabel, ylab = Rslabel,  
     xlim = c(min(Temp), max(Temp)), col = "lightgray",
     ylim = c(min(Rs), max(Rs)))
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Temp, Rs, pch = 0, xlab= Templabel, ylab = Rslabel, 
     xlim = c(min(Temp), max(Temp)), col = "lightgray",
     ylim = c(min(Rs), max(Rs)))
par(new=TRUE)
plot(Temp_aclhs_sub50, Rs_aclhs_sub50,  xlab= Templabel, ylab = Rslabel,
     pch = 18, col = "blue", cex = 1.5,
     xlim = c(min(Temp), max(Temp)),  ylim = c(min(Rs), max(Rs)))
box()
dev.off()


#------------------------------- Figure 13--------------------------------------#

Variolabel <- expression(bold("Semivariance"))
N_lags<- 2 
lag_value <- 400
Rs_fixed_sub50_VarioEstimation<-Variograma(X_fixed_sub50, Y_fixed_sub50, 
                                           Rs_fixed_sub50, 0, 90, N_lags, lag_value, 1, "", Xlabel) # Rs Spatial variogram
N_lags<- 10
lag_value <- 100 
Rs_clhs_sub50_VarioEstimation<-Variograma(X_clhs_sub50, Y_clhs_sub50,  
                                          Rs_clhs_sub50, 0, 90, N_lags, lag_value, 1, "", Xlabel) # Rs Spatial variogram
Rs_aclhs_sub50_VarioEstimation<-Variograma(X_aclhs_sub50, Y_aclhs_sub50,
                                           Rs_aclhs_sub50, 0, 90, N_lags, lag_value, 1, "", Xlabel) # Rs Spatial variogram

#Variogram models (1- exponential, 2- spherical, 3- gaussian)
N_lags<- 10
lag_value <- 100
# Rs_aclhs_sub50_vario_model<- 2 
# Rs_aclhs_sub50_nugget<- 10000
# Rs_aclhs_sub50_sill_and_nugget<- 125000
# Rs_aclhs_sub50_range <- 2800
Dislabel <- "Distance [km]"
png(paste(SM_dir,"/2_Variogram_comparison_Temp_Rs.png",sep=""), bg = "white",  width = 1500, height = 1000, res = 200)
par(mfrow = c(1,1),mar = c(4, 6, 2, 2), cex.lab= 1.2, cex.axis = 1.2)
xzoom <- max(Rs_VarioEstimation[,2])*1.1
yzoom <- max(Rs_VarioEstimation[,3])*1.3
plot(Rs_VarioEstimation[,c(2,3)], pch = 0, xlim = c(0,xzoom),  ylim=c(0,yzoom),
     xlab = Dislabel,ylab = Variolabel)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
abline(lm(Rs_VarioEstimation$Semivarianzas~Rs_VarioEstimation$Lags), lwd = 2, lty = 1)
par(new=TRUE)
plot(Rs_VarioEstimation[,c(2,3)], pch = 0, xlim = c(0,xzoom), ylim=c(0,yzoom),
     xlab = Dislabel ,ylab = Variolabel )
par(new=TRUE)
plot(Rs_fixed_sub50_VarioEstimation[,c(2,3)], pch = 15, cex = 1.2, col = "green",  xlim = c(0,xzoom), ylim=c(0,yzoom),
     xlab = Dislabel ,ylab = Variolabel)
par(new=TRUE)
plot(Rs_clhs_sub50_VarioEstimation[,c(2,3)], pch = 17, cex = 1.2, col = "red",   xlim = c(0,xzoom), ylim=c(0,yzoom),
     xlab = Dislabel , ylab = Variolabel)
par(new=TRUE)
plot(Rs_aclhs_sub50_VarioEstimation[,c(2,3)], pch = 18, cex = 1.5, col = "blue",  xlim = c(0,xzoom), ylim=c(0,yzoom),
     xlab = Dislabel , ylab = Variolabel)
box()
dev.off()
sum(abs(Rs_aclhs_sub50_VarioEstimation[,3]-Rs_VarioEstimation[,3]))

png(paste(SM_dir,"/2_Variogram_comparison_Temp_Rs_aclhs.png",sep=""), bg = "white",  width = 1500, height = 1000, res = 200)
par(mfrow = c(1,1),mar = c(4, 6, 2, 2), cex.lab= 1.2, cex.axis = 1.2)
xzoom <- max(Rs_VarioEstimation[,2])*1.0
yzoom <- max(Rs_VarioEstimation[,3])*1.0
plot(Rs_VarioEstimation[,c(2,3)], pch = 0, xlim = c(100,xzoom),  ylim=c(100,yzoom),
     xlab = Dislabel,ylab = Variolabel)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
abline(lm(Rs_VarioEstimation$Semivarianzas~Rs_VarioEstimation$Lags), lwd = 2, lty = 1)
#abline(lm(Rs_aclhs_sub50_VarioEstimation$Semivarianzas~Rs_aclhs_sub50_VarioEstimation$Lags), lwd = 2, lty = 1)
par(new=TRUE)
plot(Rs_VarioEstimation[,c(2,3)], pch = 0, xlim = c(100,xzoom), ylim=c(100,yzoom),
     xlab = Dislabel ,ylab = Variolabel )
par(new=TRUE)
plot(Rs_aclhs_sub50_VarioEstimation[,c(2,3)], pch = 18, cex = 1.5, col = "blue",  xlim = c(100,xzoom), ylim=c(100,yzoom),
     xlab = Dislabel , ylab = Variolabel)
box()
dev.off()

#----------------------Tables------------------#
# Statistical Properties
cbind(Rs_Stat,Rs_fixed_sub50_Stat,Rs_clhs_sub50_Stat,Rs_aclhs_sub50_Stat)
cbind(Temp_Stat,Temp_fixed_sub50_Stat,Temp_clhs_sub50_Stat,Temp_aclhs_sub50_Stat)

# Kolmogorov Smirnov test
ks_test_fixed_Rs <- ks.test(Rs_fixed_sub50, Rs, conf.level = 0.95)
ks_test_clhs_Rs <- ks.test(Rs_clhs_sub50, Rs, conf.level = 0.95)
ks_test_aclhs_Rs <- ks.test(Rs_aclhs_sub50, Rs, conf.level = 0.95)
ks_test_fixed_Rs
ks_test_clhs_Rs
ks_test_aclhs_Rs
ks_test_fixed_Temp <- ks.test(Temp_fixed_sub50, Temp, conf.level = 0.95)
ks_test_clhs_Temp <- ks.test(Temp_clhs_sub50, Temp, conf.level = 0.95)
ks_test_aclhs_Temp <- ks.test(Temp_aclhs_sub50, Temp, conf.level = 0.95)
ks_test_fixed_Temp
ks_test_clhs_Temp
ks_test_aclhs_Temp

# Data
(corP_Temp_Rs <- round(cor(Temp,Rs, method = "pearson"),3))
(corS_Temp_Rs <- round(cor(Temp,Rs, method = "spearman"),3))
(corK_Temp_Rs <- round(cor(Temp,Rs, method = "kendall"),3))
# fixed
(corP_Temp_Rs_fixed_sub50 <- round(cor(Temp_fixed_sub50,Rs_fixed_sub50, method = "pearson"),3))
(corS_Temp_Rs_fixed_sub50 <- round(cor(Temp_fixed_sub50,Rs_fixed_sub50, method = "spearman"),3))
(corK_Temp_Rs_fixed_sub50 <- round(cor(Temp_fixed_sub50,Rs_fixed_sub50, method = "kendall"),3))
# clhs
(corP_Temp_Rs_clhs_sub50 <- round(cor(Temp_clhs_sub50,Rs_clhs_sub50, method = "pearson"),3))
(corS_Temp_Rs_clhs_sub50 <- round(cor(Temp_clhs_sub50,Rs_clhs_sub50, method = "spearman"),3))
(corK_Temp_Rs_clhs_sub50 <- round(cor(Temp_clhs_sub50,Rs_clhs_sub50, method = "kendall"),3))
#aclhs
(corP_Temp_Rs_aclhs_sub50 <- round(cor(Temp_aclhs_sub50,Rs_aclhs_sub50, method = "pearson"),3))
(corS_Temp_Rs_aclhs_sub50 <- round(cor(Temp_aclhs_sub50,Rs_aclhs_sub50, method = "spearman"),3))
(corK_Temp_Rs_aclhs_sub50 <- round(cor(Temp_aclhs_sub50,Rs_aclhs_sub50, method = "kendall"),3))

abs(corP_Temp_Rs_fixed_sub50-corP_Temp_Rs)
abs(corS_Temp_Rs_fixed_sub50-corS_Temp_Rs)
abs(corK_Temp_Rs_fixed_sub50-corK_Temp_Rs)
abs(corP_Temp_Rs_clhs_sub50-corP_Temp_Rs)
abs(corS_Temp_Rs_clhs_sub50-corS_Temp_Rs)
abs(corK_Temp_Rs_clhs_sub50-corK_Temp_Rs)
abs(corP_Temp_Rs_aclhs_sub50-corP_Temp_Rs)
abs(corS_Temp_Rs_aclhs_sub50-corS_Temp_Rs)
abs(corK_Temp_Rs_aclhs_sub50-corK_Temp_Rs)
#
sum(abs(Rs_clhs_sub50_VarioEstimation[,3]-Rs_VarioEstimation[,3]))
sum(abs(Rs_aclhs_sub50_VarioEstimation[,3]-Rs_VarioEstimation[,3]))
