# Creates a folder to store results for Simulations (Sims)
dir.create(paste(getwd(),"/Results/spatial_2D/Sims", sep=""))
Sims_dir<-paste(root_dir,"/Results/spatial_2D/Sims",sep="")

#------------------------------Simulations of Rs conditioned to Temp--------------------------#
#-------------------------------------------FS-----------------------------------#
Fixed_Sample <- data.frame(Temp_fixed_sub50, Rs_fixed_sub50) 
muestra <- Fixed_Sample
matriz.copem <- genmat.copem(muestra)
m <- 10 # m simulations per Temp point
Temp_aux = Temp
for (i in 1:length(Temp_aux)) {
  if (Temp_aux[i] > max(Temp_fixed_sub50) ){
    Temp_aux[i] = max(Temp_fixed_sub50)
  }
  if (Temp_aux[i] < min(Temp_fixed_sub50) ){
    Temp_aux[i] = min(Temp_fixed_sub50)
  }
}
system.time(sim.m_fixed <- mapply(simula.Bernshtein.condicional, Temp_aux, rep(m, length(Temp_aux))))
# sim.m_fixed is a matrix with m rows and length(Temp) columns
# be patient... 8 seconds approx 


#-----------Optimization------------------#
n_sim <- length(Rs)
Fn_limit <- seq(0,1,length.out = n_sub+1)
Rs_fixed_sub50_limit <- quantile(Rs_fixed_sub50,Fn_limit)

#--------------OF1------------#
Func01sim  <- function(var_sim) {
  hist_Rs_sim <- hist(var_sim, breaks = Rs_fixed_sub50_limit, plot = FALSE)$counts
  freq <- rep(n_sim/n_sub,n_sub)
  return (sum(abs(hist_Rs_sim-freq)))
}


corP_Temp_Rs_fixed_sub50 <- round(cor(Temp_fixed_sub50,Rs_fixed_sub50, method = "pearson"),3)
corS_Temp_Rs_fixed_sub50 <- round(cor(Temp_fixed_sub50,Rs_fixed_sub50, method = "spearman"),3)
corK_Temp_Rs_fixed_sub50 <- round(cor(Temp_fixed_sub50,Rs_fixed_sub50, method = "kendall"),3)
#------------OF2----------------#
Func02sim  <- function(var_sim) {
  corP_Temp_var_sim <- round(cor(Temp,var_sim, method = "pearson"),3)
  corS_Temp_var_sim <- round(cor(Temp,var_sim, method = "spearman"),3)
  corK_Temp_var_sim <- round(cor(Temp,var_sim, method = "kendall"),3)
  return (abs(corP_Temp_var_sim-corP_Temp_Rs_fixed_sub50) + abs(corS_Temp_var_sim-corS_Temp_Rs_fixed_sub50) + abs(corK_Temp_var_sim-corK_Temp_Rs_fixed_sub50)) 
}


#------------ OF = w1 * OF1 + w2 * OF2 + w3 * OF3--------#
ws1 = 1; ws2 = 1000;
FuncConjsim  <- function(var_sim) {
  total =  ws1*Func01sim(var_sim)  + ws2*Func02sim(var_sim) 
  #  Prior:              267.84    +  0.02707472          
  #  min invidual:       5.64      +  1.402266e-05        
  #  itermax =5000:     1*16.92    + 500 * 0.1070009      
  #  itermax =5000:     3*11.28    + 1000 * 0.02927405    
  return (total) 
}

#DE
Pop_Inicial_sim_fixed <- as.matrix(sim.m_fixed)
vector_min =  apply(Pop_Inicial_sim_fixed, 2, min)  
vector_max =  apply(Pop_Inicial_sim_fixed, 2, max)  
vector_min[pos_sub_fixed] = Rs_fixed_sub50
vector_max[pos_sub_fixed] = Rs_fixed_sub50
system.time(
  fixed_sim_outDEoptim_Rs <- DEoptim(fn=FuncConjsim, lower= vector_min,
                                     upper= vector_max,
                                     control = DEoptim.control(VTR = 0.0000001,strategy = 3, 
                                                               itermax =2000, reltol = 1e-8, 
                                                               CR = 0.5, F = 0.8, NP= nrow(Pop_Inicial_sim_fixed), initialpop = Pop_Inicial_sim_fixed)) #, , 
)
# user       system      elapsed 
# 233.858  13.889 245.408, itermax =2000, 3 functions: 2+1000, fixed_sim_outDEoptim_Rs$optim$bestval = 216
# 256.440   4.261 260.809, , itermax =2000, 3 functions: 1+1, fixed_sim_outDEoptim_Rs$optim$bestval = 111.383
fixed_sim_outDEoptim_Rs$optim$bestval
Ite <- seq(1,length(fixed_sim_outDEoptim_Rs$member$bestvalit),1)
plot(Ite,fixed_sim_outDEoptim_Rs$member$bestvalit, "l", xlab ="Iterations", ylab ="Objective Function")
var_sim <- as.numeric(fixed_sim_outDEoptim_Rs$optim$bestmem)
#var_sim <- Pop_Inicial_sim_fixed[1,]



#-------------------------------------------cLHS-----------------------------------#
clhs_Sample <- data.frame(Temp_clhs_sub50, Rs_clhs_sub50) 
muestra <- clhs_Sample
matriz.copem <- genmat.copem(muestra)
m <- 10 # m simulations per Temp point
Temp_aux = Temp
for (i in 1:length(Temp_aux)) {
  if (Temp_aux[i] > max(Temp_clhs_sub50) ){
    Temp_aux[i] = max(Temp_clhs_sub50)
  }
  if (Temp_aux[i] < min(Temp_clhs_sub50) ){
    Temp_aux[i] = min(Temp_clhs_sub50)
  }
}
system.time(sim.m_clhs <- mapply(simula.Bernshtein.condicional, Temp_aux, rep(m, length(Temp_aux))))
# sim.m_clhs is a matrix with m rows and length(Temp) columns
# be patient... 9 seconds approx 


#-----------Optimization------------------#
n_sim <- length(Rs)
Fn_limit <- seq(0,1,length.out = n_sub+1)
Rs_clhs_sub50_limit <- quantile(Rs_clhs_sub50,Fn_limit)

#--------------OF1------------#
Func01sim  <- function(var_sim) {
  hist_Rs_sim <- hist(var_sim, breaks = Rs_clhs_sub50_limit, plot = FALSE)$counts
  freq <- rep(n_sim/n_sub,n_sub)
  return (sum(abs(hist_Rs_sim-freq)))
}


corP_Temp_Rs_clhs_sub50 <- round(cor(Temp_clhs_sub50,Rs_clhs_sub50, method = "pearson"),3)
corS_Temp_Rs_clhs_sub50 <- round(cor(Temp_clhs_sub50,Rs_clhs_sub50, method = "spearman"),3)
corK_Temp_Rs_clhs_sub50 <- round(cor(Temp_clhs_sub50,Rs_clhs_sub50, method = "kendall"),3)
#------------OF2----------------#
Func02sim  <- function(var_sim) {
  corP_Temp_var_sim <- round(cor(Temp,var_sim, method = "pearson"),3)
  corS_Temp_var_sim <- round(cor(Temp,var_sim, method = "spearman"),3)
  corK_Temp_var_sim <- round(cor(Temp,var_sim, method = "kendall"),3)
  return (abs(corP_Temp_var_sim-corP_Temp_Rs_clhs_sub50) + abs(corS_Temp_var_sim-corS_Temp_Rs_clhs_sub50) + abs(corK_Temp_var_sim-corK_Temp_Rs_clhs_sub50)) 
}


#------------ OF = w1 * OF1 + w2 * OF2 + w3 * OF3--------#
ws1 = 1; ws2 = 1; 
FuncConjsim  <- function(var_sim) {
  total =  ws1*Func01sim(var_sim)  + ws2*Func02sim(var_sim) 
  #  Prior:              267.84    +  0.02707472          
  #  min invidual:       5.64      +  1.402266e-05        
  #  itermax =5000:     1*16.92    + 500 * 0.1070009      
  #  itermax =5000:     3*11.28    + 1000 * 0.02927405   
  return (total) 
}

#DE
Pop_Inicial_sim_clhs <- as.matrix(sim.m_clhs)
vector_min =  apply(Pop_Inicial_sim_clhs, 2, min)  
vector_max =  apply(Pop_Inicial_sim_clhs, 2, max)  
vector_min[clhs_sub50] = Rs_clhs_sub50
vector_max[clhs_sub50] = Rs_clhs_sub50
system.time(
  clhs_sim_outDEoptim_Rs <- DEoptim(fn=FuncConjsim, lower= vector_min,
                                    upper= vector_max,
                                    control = DEoptim.control(VTR = 0.0000001,strategy = 3, 
                                                              itermax =2000, reltol = 1e-8, 
                                                              CR = 0.5, F = 0.8, NP= nrow(Pop_Inicial_sim_clhs), initialpop = Pop_Inicial_sim_clhs)) #, , 
)
# user       system      elapsed 
# 233.858  13.889 245.408, itermax =2000, 3 functions: 2+1000, clhs_sim_outDEoptim_Rs$optim$bestval = 216
# 243.247  10.515 252.310, itermax =2000, 3 functions: 1+1, clhs_sim_outDEoptim_Rs$optim$bestval =  31.9712
clhs_sim_outDEoptim_Rs$optim$bestval
Ite <- seq(1,length(clhs_sim_outDEoptim_Rs$member$bestvalit),1)
plot(Ite,clhs_sim_outDEoptim_Rs$member$bestvalit, "l", xlab ="Iterations", ylab ="Objective Function")
var_sim <- as.numeric(clhs_sim_outDEoptim_Rs$optim$bestmem)
#var_sim <- Pop_Inicial_sim_clhs[1,]

#-------------------------------------------acLHS-----------------------------------#
aclhs_Sample <- data.frame(Temp_aclhs_sub50, Rs_aclhs_sub50) 
muestra <- aclhs_Sample
matriz.copem <- genmat.copem(muestra)
m <- 10 # simulations per Temp point
Temp_aux = Temp
for (i in 1:length(Temp_aux)) {
  if (Temp_aux[i] > max(Temp_aclhs_sub50) ){
    Temp_aux[i] = max(Temp_aclhs_sub50)
  }
  if (Temp_aux[i] < min(Temp_aclhs_sub50) ){
    Temp_aux[i] = min(Temp_aclhs_sub50)
  }
}
system.time(sim.m_aclhs <- mapply(simula.Bernshtein.condicional, Temp_aux, rep(m, length(Temp_aux))))
# sim.m_aclhs is a matrix with m rows and length(Temp) columns
# be patient... 29 segundos approx

#-----------Optimization------------------#
n_sim <- length(Rs)
Fn_limit <- seq(0,1,length.out = n_sub+1)
Rs_aclhs_sub50_limit <- quantile(Rs_aclhs_sub50,Fn_limit)

#--------------OF1------------#
Func01sim  <- function(var_sim) {
  hist_Rs_sim <- hist(var_sim, breaks = Rs_aclhs_sub50_limit, plot = FALSE)$counts
  freq <- rep(n_sim/n_sub,n_sub)
  return (sum(abs(hist_Rs_sim-freq)))
}


corP_Temp_Rs_aclhs_sub50 <- round(cor(Temp_aclhs_sub50,Rs_aclhs_sub50, method = "pearson"),3)
corS_Temp_Rs_aclhs_sub50 <- round(cor(Temp_aclhs_sub50,Rs_aclhs_sub50, method = "spearman"),3)
corK_Temp_Rs_aclhs_sub50 <- round(cor(Temp_aclhs_sub50,Rs_aclhs_sub50, method = "kendall"),3)
#------------OF2----------------#
Func02sim  <- function(var_sim) {
  corP_Temp_var_sim <- round(cor(Temp,var_sim, method = "pearson"),3)
  corS_Temp_var_sim <- round(cor(Temp,var_sim, method = "spearman"),3)
  corK_Temp_var_sim <- round(cor(Temp,var_sim, method = "kendall"),3)
  return (abs(corP_Temp_var_sim-corP_Temp_Rs_aclhs_sub50) + abs(corS_Temp_var_sim-corS_Temp_Rs_aclhs_sub50) + abs(corK_Temp_var_sim-corK_Temp_Rs_aclhs_sub50)) 
}

#------------OF3------------#
Func03sim <- function(var_sim) {
  N_lags<- 10
  lag_value <- 100
  X_sim <- X
  Y_sim <- Y
  Vario_var_sim <- variog(as.geodata(cbind(X_sim, Y_sim, var_sim), 
                                     coords.col = 1:2, data.col = 3), breaks = c(seq(lag_value/2, lag_value * (N_lags+1), lag_value)),
                          trend = "cte", lambda = 1, estimator.type = "classical", nugget.tolerance = 0, direction = 0, 
                          tolerance = 90, unit.angle = "degrees", pairs.min = 1)
  return (sum(abs(Vario_var_sim$v-Rs_VarioEstimation[,3]))) 
}
#------------ OF = w1 * OF1 + w2 * OF2 + w3 * OF3--------#
ws1 = 1; ws2 = 100; ws3 = 0.001#ws1 = 1; ws2 = 100; ws3 = 0.001 # ws1 = 1; ws2 = 100; ws3 = 0.001 #ws1 = 10; ws2 = 100; ws3 = 0.001 # ws1 = 0.1; ws2 = 10; ws3 = 5
FuncConjsim  <- function(var_sim) {
  total =  ws1*Func01sim(var_sim)  + ws2*Func02sim(var_sim) + ws3*Func03sim(var_sim) 
  #  Prior:              267.84    +  0.02707472          +  138102.8
  #  min invidual:       5.64      +  1.402266e-05        +  17102.48/3000
  #  itermax =5000:     1*16.92    + 500 * 0.1070009      +  0.01*4345.137
  #  itermax =5000:     3*11.28    + 1000 * 0.02927405    +  0.01*10590.71
  return (total) 
}

#DE
Pop_Inicial_sim <- as.matrix(sim.m_aclhs)
for (i in 1:length(aclhs_sub50)) {
  for (j in 1:nrow(Pop_Inicial_sim)) {
    Pop_Inicial_sim[j,aclhs_sub50[i]] = Rs_aclhs_sub50[i]
  }
}
vector_min =  apply(sim.m_aclhs, 2, min) #apply(sim.m_aclhs, 2, quantile, probs=c(0.001), na.rm=TRUE) # apply(Pop_Inicial_sim, 2, min) # apply(sim.m_aclhs, 2, min) 
vector_max =  apply(sim.m_aclhs, 2, max) #apply(sim.m_aclhs, 2, quantile, probs=c(0.999), na.rm=TRUE) #apply(Pop_Inicial_sim, 2, max) # apply(sim.m_aclhs, 2, max)
vector_min[aclhs_sub50] = Rs_aclhs_sub50
vector_max[aclhs_sub50] = Rs_aclhs_sub50
system.time(
aclhs_sim_outDEoptim_Rs <- DEoptim(fn=FuncConjsim, lower= vector_min,
                                      upper= vector_max,
                                      control = DEoptim.control(VTR = 0.0000001,strategy = 3, 
                                                                itermax =40000, reltol = 1e-8, 
                                                                CR = 0.5, F = 0.8,NP= nrow(Pop_Inicial_sim), initialpop = Pop_Inicial_sim)) # 
)
# user       system      elapsed 
# 1684.853   68.058 1752.046, itermax =7500, 3 functions: 1+500+0.01, aclhs_sim_outDEoptim_Rs$optim$bestval = 113.8718
# 1613.194   60.222 1662.723, itermax =7500, 3 functions: 1+1000+0.01, aclhs_sim_outDEoptim_Rs$optim$bestval = 84.11742
# 1709.974   79.104 1794.196, itermax =7500, 3 functions: 3+1000+0.01, aclhs_sim_outDEoptim_Rs$optim$bestval = 169.0212
# 1036.370   51.160 1074.323, itermax =5000, 3 functions: 3+1000+0.01, aclhs_sim_outDEoptim_Rs$optim$bestval = 354.3757
# 1540.810   61.189 1582.417, itermax =7500, 3 functions: 3+1000+0.005, aclhs_sim_outDEoptim_Rs$optim$bestval = 315.2039
# 2061.410   75.987 2111.391, itermax =10000, 3 functions: 3+1000+0.003, aclhs_sim_outDEoptim_Rs$optim$bestval = 366.7643
# 4157.773  147.721 4260.089, itermax =10000, 3 functions: 1+100+0.001, aclhs_sim_outDEoptim_Rs$optim$bestval = 28.15154
# 4113.713  367.773 4435.324, itermax =10000, 3 functions: 1+100+0.001, aclhs_sim_outDEoptim_Rs$optim$bestval = 35.47882 
# 2056.955   71.148 2105.685  itermax =10000, 3 functions: 3+100+0.001, aclhs_sim_outDEoptim_Rs$optim$bestval = 120.0764
# 4230.171  186.682 4372.547  itermax =20000, 3 functions: 3+100+0.001, aclhs_sim_outDEoptim_Rs$optim$bestval = 56.96843
# 4223.926  492.265 4669.924  itermax =20000, 3 functions: 3+100+0.0001, aclhs_sim_outDEoptim_Rs$optim$bestval = 35.35535
# 4297.320  166.357 4419.602  itermax =20000, 3 functions: 3+100+0.0001, aclhs_sim_outDEoptim_Rs$optim$bestval = 33.74357
#4122.170  175.512 4253.573   itermax =20000, 3 functions: 3+100+0.0001, aclhs_sim_outDEoptim_Rs$optim$bestval = 35.31221
# 4191.891  340.860 4488.657  itermax =20000, 3 functions: 3+100+0.001, aclhs_sim_outDEoptim_Rs$optim$bestval = 33.30092
# 1028.912   40.440 1058.339  itermax =5000, 3 functions: 3+100+0.0001, aclhs_sim_outDEoptim_Rs$optim$bestval = 39.1114
# 2060.739   82.005 2120.850  itermax =10000, 3 functions: 3+100+0.001, aclhs_sim_outDEoptim_Rs$optim$bestval = 121.9743
# 4098.191  140.493 4195.205  itermax =10000, 3 functions: 3+100+0.001, aclhs_sim_outDEoptim_Rs$optim$bestval = 65.34106
# 4111.823  150.098 4215.569  itermax =20000, 3 functions: 3+100+0.001, aclhs_sim_outDEoptim_Rs$optim$bestval = 126.0213
# 2107.729  229.694 2314.121  itermax =10000, 3 functions: 3+300+0.001, aclhs_sim_outDEoptim_Rs$optim$bestval = 204.1339
# 2117.044  103.480 2198.544  itermax =20000, 3 functions: 1+100+0.001, aclhs_sim_outDEoptim_Rs$optim$bestval = 115.3875
# 8420.257 1026.881 9351.671  itermax =40000, 3 functions: 1+100+0.001, aclhs_sim_outDEoptim_Rs$optim$bestval = 26.12571
aclhs_sim_outDEoptim_Rs$optim$bestval
Ite <- seq(1,length(aclhs_sim_outDEoptim_Rs$member$bestvalit),1)
plot(Ite,aclhs_sim_outDEoptim_Rs$member$bestvalit, "l", xlab ="Iterations", ylab ="Objective Function")
var_sim <- as.numeric(aclhs_sim_outDEoptim_Rs$optim$bestmem)
#var_sim <- Pop_Inicial_sim[1,]

#------------------- tables-------------------#
# Statistical Properties
Rs_fixed_pred <- fixed_sim_outDEoptim_Rs$optim$bestmem #sim.m_fixed #[1,]
Rs_fixed_pred[pos_sub_fixed] <- Rs_fixed_sub50
Rs_clhs_pred <- clhs_sim_outDEoptim_Rs$optim$bestmem #sim.m_clhs #[1,]
Rs_clhs_pred[clhs_sub50] <- Rs_clhs_sub50
Rs_aclhs_pred <- aclhs_sim_outDEoptim_Rs$optim$bestmem
Rs_aclhs_pred[aclhs_sub50] <- Rs_aclhs_sub50
Rs_fixed_pred_Stat <- Estadisticas(Rs_fixed_pred)
Rs_clhs_pred_Stat <- Estadisticas(Rs_clhs_pred)
Rs_aclhs_pred_Stat <- Estadisticas(Rs_aclhs_pred)
cbind(Rs_Stat,Rs_fixed_pred_Stat,Rs_clhs_pred_Stat,Rs_aclhs_pred_Stat)


ks_test_fixed_Rs_sim <- ks.test(Rs_fixed_pred, Rs, conf.level = 0.95)
ks_test_clhs_Rs_sim <- ks.test(Rs_clhs_pred, Rs, conf.level = 0.95)
#ks_test_aclhs_Rs_sim1 <- ks.test(sim.m_aclhs[1,], Rs, conf.level = 0.95)
ks_test_aclhs_Rs_sim2 <- ks.test(as.numeric(Rs_aclhs_pred), Rs, conf.level = 0.95)
ks_test_fixed_Rs_sim
ks_test_clhs_Rs_sim
#ks_test_aclhs_Rs_sim1
ks_test_aclhs_Rs_sim2

# Data
(corP_Temp_Rs <- round(cor(Temp,Rs, method = "pearson"),3))
(corS_Temp_Rs <- round(cor(Temp,Rs, method = "spearman"),3))
(corK_Temp_Rs <- round(cor(Temp,Rs, method = "kendall"),3))
# fixed
(corP_Temp_Rs_fixed_pred <- round(cor(Temp,Rs_fixed_pred, method = "pearson"),3))
(corS_Temp_Rs_fixed_pred <- round(cor(Temp,Rs_fixed_pred, method = "spearman"),3))
(corK_Temp_Rs_fixed_pred <- round(cor(Temp,Rs_fixed_pred, method = "kendall"),3))
abs(corP_Temp_Rs_fixed_sub50-corP_Temp_Rs_fixed_pred)
abs(corS_Temp_Rs_fixed_sub50-corS_Temp_Rs_fixed_pred)
abs(corK_Temp_Rs_fixed_sub50-corK_Temp_Rs_fixed_pred)
# clhs
(corP_Temp_Rs_clhs_pred <- round(cor(Temp,Rs_clhs_pred, method = "pearson"),3))
(corS_Temp_Rs_clhs_pred <- round(cor(Temp,Rs_clhs_pred, method = "spearman"),3))
(corK_Temp_Rs_clhs_pred <- round(cor(Temp,Rs_clhs_pred, method = "kendall"),3))
abs(corP_Temp_Rs_clhs_sub50-corP_Temp_Rs_clhs_pred)
abs(corS_Temp_Rs_clhs_sub50-corS_Temp_Rs_clhs_pred)
abs(corK_Temp_Rs_clhs_sub50-corK_Temp_Rs_clhs_pred)
#aclhs
(corP_Temp_Rs_aclhs_pred <- round(cor(Temp,Rs_aclhs_pred, method = "pearson"),3))
(corS_Temp_Rs_aclhs_pred <- round(cor(Temp,Rs_aclhs_pred, method = "spearman"),3))
(corK_Temp_Rs_aclhs_pred <- round(cor(Temp,Rs_aclhs_pred, method = "kendall"),3))
abs(corP_Temp_Rs_aclhs_sub50-corP_Temp_Rs_aclhs_pred)
abs(corS_Temp_Rs_aclhs_sub50-corS_Temp_Rs_aclhs_pred)
abs(corK_Temp_Rs_aclhs_sub50-corK_Temp_Rs_aclhs_pred)
#
abs(corP_Temp_Rs_fixed_pred-corP_Temp_Rs)
abs(corS_Temp_Rs_fixed_pred-corS_Temp_Rs)
abs(corK_Temp_Rs_fixed_pred-corK_Temp_Rs)
#
abs(corP_Temp_Rs_clhs_pred-corP_Temp_Rs)
abs(corS_Temp_Rs_clhs_pred-corS_Temp_Rs)
abs(corK_Temp_Rs_clhs_pred-corK_Temp_Rs)
# 
abs(corP_Temp_Rs_aclhs_pred-corP_Temp_Rs)
abs(corS_Temp_Rs_aclhs_pred-corS_Temp_Rs)
abs(corK_Temp_Rs_aclhs_pred-corK_Temp_Rs)

(sum(abs(Rs-as.numeric(Rs_fixed_pred)))/sum(Rs))*100
(sum(abs(Rs-as.numeric(Rs_clhs_pred)))/sum(Rs))*100
(sum(abs(Rs-as.numeric(Rs_aclhs_pred)))/sum(Rs))*100


sum(abs(Rs-as.numeric(Rs_fixed_pred)))
sum(abs(Rs-as.numeric(Rs_clhs_pred)))
sum(abs(Rs-as.numeric(Rs_aclhs_pred)))


# sum(Rs) [gCm^-2year^-1]  * 100 km * 100 km * 10^6 m^2 / 10^15 grams  
#sum(Rs)*100*100*1000000/1000000000000000
sum(abs(Rs-as.numeric(Rs_fixed_pred)))*100*100*1000000/1000000000000000
sum(abs(Rs-as.numeric(Rs_clhs_pred)))*100*100*1000000/1000000000000000
sum(abs(Rs-as.numeric(Rs_aclhs_pred)))*100*100*1000000/1000000000000000

#--------------------------------Figure 14--------------------------------------#

png(paste(Sims_dir,"/14a_Rs_Spatial_Distr.png",sep=""), 
    bg = "white", width = 3000, height = 2000, res = 250)
DEspacial(X, Y, Rs, 
          Xlabel , Ylabel, Rslabel, breaks = histf_Rs$breaks,  cex.lab = 2, cex.axis = 2)
dev.off()

png(paste(Sims_dir,"/14b_Rs_fixed_pred_Spatial_Distr.png",sep=""), 
    bg = "white", width = 3000, height = 2000, res = 250)
DEspacial(X, Y, Rs_fixed_pred, 
          Xlabel , Ylabel, Rslabel,  breaks = histf_Rs$breaks, cex.lab = 2, cex.axis = 2)
dev.off()

png(paste(Sims_dir,"/14c_Rs_clhs_pred_Spatial_Distr.png",sep=""), 
    bg = "white", width = 3000, height = 2000, res = 250)
DEspacial(X, Y, Rs_clhs_pred, 
          Xlabel , Ylabel, Rslabel, breaks = histf_Rs$breaks, cex.lab = 2, cex.axis = 2)
dev.off()

png(paste(Sims_dir,"/14d_Rs_aclhs_pred_Spatial_Distr.png",sep=""), 
    bg = "white", width = 3000, height = 2000, res = 250)
DEspacial(X, Y, Rs_aclhs_pred, 
          Xlabel , Ylabel, Rslabel, breaks = histf_Rs$breaks, cex.lab = 2, cex.axis = 2)
dev.off()


#---------------------------------------------Figure5_Univariate pdfs--------------------------------#
#Fn_Temp <- estandarizar(cbind(Temp,Rs))[,1]
#Fn_Rs <- estandarizar(cbind(Temp,Rs))[,2]
U_limit <- seq(0,1,by = 0.2)
png(paste(Sims_dir,"/2_Rs_univariate_pdfs_sims.png",sep=""),   bg = "white", width = 1200, height = 1000, res = 120)
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
plot(ecdf(Rs_fixed_pred),pch = 15, col = "green", cex = 1.2,
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
plot(ecdf(Rs_clhs_pred),pch = 17, col = "red", cex = 1.2,
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
plot(ecdf(Rs_aclhs_pred),pch = 18, col = "blue", cex = 1.5,
     yaxt = "n",  main = " ", xlab= Rslabel, ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),
     xlim =c(min(Rs),max(Rs)), ylim = c(0,1))
box()
dev.off()



#---------------------------------------------Figure6_scatterplots--------------------------------#
Fn_Temp <- estandarizar(cbind(Temp,Rs))[,1]
Fn_Rs <- estandarizar(cbind(Temp,Rs))[,2]
U_limit <- seq(0,1,by = 0.2)
png(paste(Sims_dir,"/2_Rs_Temp_scatterplots_sims.png",sep=""),   bg = "white", width = 1200, height = 1000, res = 120)
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
plot(Temp, Rs_fixed_pred,  xlab= Templabel, ylab = Rslabel,
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
plot(Temp, Rs_clhs_pred,  xlab= Templabel, ylab = Rslabel,
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
plot(Temp, Rs_aclhs_pred,  xlab= Templabel, ylab = Rslabel,
     pch = 18, col = "blue", cex = 1.5,
     xlim = c(min(Temp), max(Temp)),  ylim = c(min(Rs), max(Rs)))
box()

dev.off()


#------------------------------- Figure 16--------------------------------------#
Variolabel <- expression(bold("Semivariance"))
# Estimation of the experimental variogram
N_lags<- 10 #length(t)/2.5
lag_value <- 100 # min(dist(t)) # or delta t

#png(paste(VA_dir,"/Rs_fixed_sim_VarioEstimation.png",sep=""), bg = "white", width = 1500, height = 1000, res = 200)
Rs_fixed_sim_VarioEstimation<-Variograma(X, Y, 
                                          Rs_fixed_pred, 0, 90, N_lags, lag_value, 1, "", Dislabel) # Rs Temporal variogram
#dev.off()

#png(paste(VA_dir,"/Rs_clhs_sim_VarioEstimation.png",sep=""), bg = "white", width = 1500, height = 1000, res = 200)
Rs_clhs_sim_VarioEstimation<-Variograma(X, Y, 
                                         Rs_clhs_pred, 0, 90, N_lags, lag_value, 1, "", Dislabel) # Rs Temporal variogram
#dev.off()

#png(paste(VA_dir,"/Rs_aclhs_sim_VarioEstimation.png",sep=""), bg = "white", width = 1500, height = 1000, res = 200)
Rs_aclhs_sim_VarioEstimation<-Variograma(X, Y, 
                                          as.numeric(Rs_aclhs_pred), 0, 90, N_lags, lag_value, 1, "", Dislabel) # Rs Temporal variogram
#dev.off()

sum(abs(Rs_fixed_sim_VarioEstimation$Semivarianzas-Rs_VarioEstimation$Semivarianzas))
sum(abs(Rs_clhs_sim_VarioEstimation$Semivarianzas-Rs_VarioEstimation$Semivarianzas))
sum(abs(Rs_aclhs_sim_VarioEstimation$Semivarianzas-Rs_VarioEstimation$Semivarianzas))

png(paste(Sims_dir,"/2_Rs_Variogram_comparison_sims.png",sep=""), bg = "white",  width = 1500, height = 1000, res = 200)
par(mfrow = c(1,1),mar = c(4, 6, 2, 2), cex.lab= 1.2, cex.axis = 1.2)
xzoom = max(Rs_VarioEstimation[,2])*1.1
yzoom = max(Rs_VarioEstimation[,3])*1.3
plot(Rs_VarioEstimation[,c(2,3)], pch = 0,  xlim = c(0,xzoom), ylim=c(0,yzoom),
     xlab = Dislabel ,ylab = Variolabel)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Rs_VarioEstimation[,c(2,3)], pch = 0,  xlim = c(0,xzoom), ylim=c(0,yzoom),
     xlab = Dislabel ,ylab = Variolabel )
par(new=TRUE)
plot(Rs_fixed_sim_VarioEstimation[,c(2,3)], pch = 15, cex = 1.2, col = "green",  xlim = c(0,xzoom), ylim=c(0,yzoom),
     xlab = Dislabel ,ylab = Variolabel)
par(new=TRUE)
plot(Rs_clhs_sim_VarioEstimation[,c(2,3)], pch = 17, cex = 1.5, col = "red",   xlim = c(0,xzoom), ylim=c(0,yzoom),
     xlab = Dislabel , ylab = Variolabel)
par(new=TRUE)
plot(Rs_aclhs_sim_VarioEstimation[,c(2,3)], pch = 18, cex = 1.5, col = "blue",  xlim = c(0,xzoom), ylim=c(0,yzoom),
     xlab = Dislabel , ylab = Variolabel)
#par(new=TRUE)
abline(lm(Rs_VarioEstimation$Semivarianzas~Rs_VarioEstimation$Lags),lwd = 3)
box()
dev.off()


