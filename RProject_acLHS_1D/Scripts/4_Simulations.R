# Create a folder to store results for Sampling Methods (SM)
dir.create(paste(getwd(),"/Results/Time_Series/Sims", sep=""))
Sims_dir<-paste(root_dir,"/Results/Time_Series/Sims",sep="")
#------------------------------Simulations of CO2 conditioned to Temp--------------------------#
#-------------------------------------------FTS-----------------------------------#
Fixed_Sample <- data.frame(Temp_fixed_sub48, CO2_fixed_sub48) 
muestra <- Fixed_Sample
matriz.copem <- genmat.copem(muestra)
m <- 10 # m simulations per Temp point
Temp_aux = Temp
for (i in 1:length(Temp_aux)) {
  if (Temp_aux[i] > max(Temp_fixed_sub48) ){
    Temp_aux[i] = max(Temp_fixed_sub48)
  }
  if (Temp_aux[i] < min(Temp_fixed_sub48) ){
    Temp_aux[i] = min(Temp_fixed_sub48)
  }
}
system.time(sim.m_fixed <- mapply(simula.Bernshtein.condicional, Temp_aux, rep(m, length(Temp_aux))))
# sim.m_fixed is a matrix with m rows and length(Temp) columns
# be patient... 24 seconds approx 


#-----------Optimization------------------#
n_sim <- 365
Fn_limit <- seq(0,1,length.out = n_sub+1)
CO2_fixed_sub48_limit <- quantile(CO2_fixed_sub48,Fn_limit)

#--------------OF1------------#
FuncO1sim  <- function(var_sim) {
  hist_CO2_sim <- hist(var_sim, breaks = CO2_fixed_sub48_limit, plot = FALSE)$counts
  freq <- rep(n_sim/n_sub,n_sub)
  return (sum(abs(hist_CO2_sim-freq)))
}


corP_Temp_CO2_fixed_sub48 <- cor(Temp_fixed_sub48,CO2_fixed_sub48, method = "pearson")
corS_Temp_CO2_fixed_sub48 <- cor(Temp_fixed_sub48,CO2_fixed_sub48, method = "spearman")
corK_Temp_CO2_fixed_sub48 <- cor(Temp_fixed_sub48,CO2_fixed_sub48, method = "kendall")
#------------OF2----------------#
FuncO2sim  <- function(var_sim) {
  corP_Temp_var_sim <- cor(Temp,var_sim, method = "pearson")
  corS_Temp_var_sim <- cor(Temp,var_sim, method = "spearman")
  corK_Temp_var_sim <- cor(Temp,var_sim, method = "kendall")
  return (abs(corP_Temp_var_sim-corP_Temp_CO2_fixed_sub48) + abs(corS_Temp_var_sim-corS_Temp_CO2_fixed_sub48) + abs(corK_Temp_var_sim-corK_Temp_CO2_fixed_sub48)) 
}

#------------ OF = w1 * OF1 + w2 * OF2 --------#
ws1 = 1; ws2 = 1
FuncConjsim  <- function(var_sim) {
  total = ws1*FuncO1sim(var_sim) + ws2*FuncO2sim(var_sim) 
  return (total) 
  #             26.125                   6.72777e-06
}

#DE
Pop_Inicial <- as.matrix(sim.m_fixed)
vector_min =  apply(Pop_Inicial, 2, min)  
vector_max =  apply(Pop_Inicial, 2, max)  
vector_min[fixed_sub48] = CO2_fixed_sub48
vector_max[fixed_sub48] = CO2_fixed_sub48
system.time(
  fixed_sim_outDEoptim_CO2 <- DEoptim(fn=FuncConjsim, lower= vector_min,
                                      upper= vector_max,
                                      control = DEoptim.control(VTR = 0.0000001,strategy = 3, 
                                                                itermax =10000, reltol = 1e-8, 
                                                                CR = 0.5, F = 0.8, NP= nrow(Pop_Inicial), initialpop = Pop_Inicial)) #, , 
)
# user       system      elapsed 
# 278.630  12.494 290.836, itermax =10000, 2 functions: 1+1, fixed_sim_outDEoptim_CO2$optim$bestval = 23.75015
# 309.893  51.494 355.060,  itermax =10000, 2 functions: 1+1, clhs_sim_outDEoptim_CO2$optim$bestval = 23.75041
fixed_sim_outDEoptim_CO2$optim$bestval
Ite <- seq(1,length(fixed_sim_outDEoptim_CO2$member$bestvalit),1)
plot(Ite,fixed_sim_outDEoptim_CO2$member$bestvalit, "l", xlab ="Iterations", ylab ="Objective Function")
var_sim <- as.numeric(fixed_sim_outDEoptim_CO2$optim$bestmem)

abs(corP_Temp_var_sim-corP_Temp_CO2_fixed_sub48) 
abs(corP_Temp_var_sim-corP_Temp_CO2)
abs(corS_Temp_var_sim-corS_Temp_CO2_fixed_sub48)
abs(corK_Temp_var_sim-corK_Temp_CO2_fixed_sub48)


#-------------------------------------------cLHS-----------------------------------#
clhs_Sample <- data.frame(Temp_clhs_sub48, CO2_clhs_sub48) 
muestra <- clhs_Sample
matriz.copem <- genmat.copem(muestra)
m <- 10 # m simulations per Temp point
Temp_aux = Temp
for (i in 1:length(Temp_aux)) {
  if (Temp_aux[i] > max(Temp_clhs_sub48) ){
    Temp_aux[i] = max(Temp_clhs_sub48)
  }
  if (Temp_aux[i] < min(Temp_clhs_sub48) ){
    Temp_aux[i] = min(Temp_clhs_sub48)
  }
}
system.time(sim.m_clhs <- mapply(simula.Bernshtein.condicional, Temp_aux, rep(m, length(Temp_aux))))
# sim.m_clhs is a matrix with m rows and length(Temp) columns
# be patient... 25 seconds approx 


#-----------Optimization------------------#
n_sim <- 365
Fn_limit <- seq(0,1,length.out = n_sub+1)
CO2_clhs_sub48_limit <- quantile(CO2_clhs_sub48,Fn_limit)

#--------------OF1------------#
FuncO1sim  <- function(var_sim) {
  hist_CO2_sim <- hist(var_sim, breaks = CO2_clhs_sub48_limit, plot = FALSE)$counts
  freq <- rep(n_sim/n_sub,n_sub)
  return (sum(abs(hist_CO2_sim-freq)))
}


corP_Temp_CO2_clhs_sub48 <- cor(Temp_clhs_sub48,CO2_clhs_sub48, method = "pearson")
corS_Temp_CO2_clhs_sub48 <- cor(Temp_clhs_sub48,CO2_clhs_sub48, method = "spearman")
corK_Temp_CO2_clhs_sub48 <- cor(Temp_clhs_sub48,CO2_clhs_sub48, method = "kendall")
#------------OF2----------------#
FuncO2sim  <- function(var_sim) {
  corP_Temp_var_sim <- cor(Temp,var_sim, method = "pearson")
  corS_Temp_var_sim <- cor(Temp,var_sim, method = "spearman")
  corK_Temp_var_sim <- cor(Temp,var_sim, method = "kendall")
  return (abs(corP_Temp_var_sim-corP_Temp_CO2_clhs_sub48) + abs(corS_Temp_var_sim-corS_Temp_CO2_clhs_sub48) + abs(corK_Temp_var_sim-corK_Temp_CO2_clhs_sub48)) 
}

#------------ OF = w1 * OF1 + w2 * OF2 + w3 * OF3--------#
ws1 = 1; ws2 = 1
FuncConjsim  <- function(var_sim) {
  total = ws1*FuncO1sim(var_sim) + ws2*FuncO2sim(var_sim) 
  return (total) 
  #             26.125                   6.72777e-06
}

#DE
Pop_Inicial <- as.matrix(sim.m_clhs)
vector_min =  apply(Pop_Inicial, 2, min)  
vector_max =  apply(Pop_Inicial, 2, max)  
vector_min[clhs_sub48] = CO2_clhs_sub48
vector_max[clhs_sub48] = CO2_clhs_sub48
system.time(
  clhs_sim_outDEoptim_CO2 <- DEoptim(fn=FuncConjsim, lower= vector_min,
                                     upper= vector_max,
                                     control = DEoptim.control(VTR = 0.0000001,strategy = 3, 
                                                               itermax =10000, reltol = 1e-8, 
                                                               CR = 0.5, F = 0.8, NP= nrow(Pop_Inicial), initialpop = Pop_Inicial)) #, , 
)
# user       system      elapsed 
# 1274.096   94.597 1353.162,  itermax =10000, 2 functions: 1+1, clhs_sim_outDEoptim_CO2$optim$bestval = 22.96115

clhs_sim_outDEoptim_CO2$optim$bestval
Ite <- seq(1,length(clhs_sim_outDEoptim_CO2$member$bestvalit),1)
plot(Ite,clhs_sim_outDEoptim_CO2$member$bestvalit, "l", xlab ="Iterations", ylab ="Objective Function")
var_sim <- as.numeric(clhs_sim_outDEoptim_CO2$optim$bestmem)

#------------------------------------------------------------------acLHS-------------------------------------------------#
aclhs_Sample <- data.frame(Temp_aclhs_sub48, CO2_aclhs_sub48) 
muestra <- aclhs_Sample
matriz.copem <- genmat.copem(muestra)
m <- 10 # simulations per Temp point
Temp_aux = Temp
for (i in 1:length(Temp_aux)) {
  if (Temp_aux[i] > max(Temp_aclhs_sub48) ){
    Temp_aux[i] = max(Temp_aclhs_sub48)
  }
  if (Temp_aux[i] < min(Temp_aclhs_sub48) ){
    Temp_aux[i] = min(Temp_aclhs_sub48)
  }
}
system.time(sim.m_aclhs <- mapply(simula.Bernshtein.condicional, Temp_aux, rep(m, length(Temp_aux))))
# sim.m_aclhs is a matrix with m rows and length(Temp) columns
# be patient... 27 segundos approx

#-----------Optimization------------------#
n_sim <- 365
Fn_limit <- seq(0,1,length.out = n_sub+1)
CO2_aclhs_sub48_limit <- quantile(CO2_aclhs_sub48,Fn_limit)

#--------------OF1------------#
FuncO1sim  <- function(var_sim) {
  hist_CO2_sim <- hist(var_sim, breaks = CO2_aclhs_sub48_limit, plot = FALSE)$counts
  freq <- rep(n_sim/n_sub,n_sub)
  return (sum(abs(hist_CO2_sim-freq)))
}


corP_Temp_CO2_aclhs_sub48 <- cor(Temp_aclhs_sub48,CO2_aclhs_sub48, method = "pearson")
corS_Temp_CO2_aclhs_sub48 <- cor(Temp_aclhs_sub48,CO2_aclhs_sub48, method = "spearman")
corK_Temp_CO2_aclhs_sub48 <- cor(Temp_aclhs_sub48,CO2_aclhs_sub48, method = "kendall")
#------------OF2----------------#
FuncO2sim  <- function(var_sim) {
  corP_Temp_var_sim <- cor(Temp,var_sim, method = "pearson")
  corS_Temp_var_sim <- cor(Temp,var_sim, method = "spearman")
  corK_Temp_var_sim <- cor(Temp,var_sim, method = "kendall")
  return (abs(corP_Temp_var_sim-corP_Temp_CO2_aclhs_sub48) + abs(corS_Temp_var_sim-corS_Temp_CO2_aclhs_sub48) + abs(corK_Temp_var_sim-corK_Temp_CO2_aclhs_sub48)) 
}

#Variogram models (1- exponential, 2- spherical, 3- gaussian)
N_lags<- 8
lag_value <- 1 
CO2_aclhs_sub48_vario_model<- 2 
CO2_aclhs_sub48_nugget<- 0.57 
CO2_aclhs_sub48_sill_and_nugget<- 1.95 
CO2_aclhs_sub48_range <- 6.21
Modelo1 = CO2_aclhs_sub48_vario_model
if (Modelo1 == 1) {
  ModeloA <- "exponential"
}
if (Modelo1 == 2) {
  ModeloA <- "spherical"
}
if (Modelo1 == 3) {
  ModeloA <- "gaussian"
}
Sill = CO2_aclhs_sub48_sill_and_nugget - CO2_aclhs_sub48_nugget
Range = CO2_aclhs_sub48_range
VMod = (CO2_aclhs_sub48_sill_and_nugget) - cov.spatial(seq(lag_value,N_lags*lag_value,lag_value), cov.model = ModeloA, 
                                                       cov.pars = c(Sill, Range))
semi_prior = VMod
#------------OF3------------#
FunCO3sim <- function(var_sim) {
  N_lags<- 8 
  lag_value <- 1 
  t_sim <- Time
  ti_sim <- rep(0,365)
  Vario_var_sim <- variog(as.geodata(cbind(t_sim, ti_sim, var_sim), 
                                     coords.col = 1:2, data.col = 3), breaks = c(seq(lag_value/2, lag_value * (N_lags+1), lag_value)),
                          trend = "cte", lambda = 1, estimator.type = "classical", nugget.tolerance = 0, direction = 0, 
                          tolerance = 90, unit.angle = "degrees", pairs.min = 1)
  return (sum(c(10,rep(1,N_lags-1))*abs(Vario_var_sim$v-semi_prior))) 
}
#------------ OF = w1 * OF1 + w2 * OF2 + w3 * OF3--------#
ws1 = 0.2; ws2 = 100; ws3 = 1 # ws1 = 0.1; ws2 = 10; ws3 = 5
FuncConjsim  <- function(var_sim) {
  total = ws1*FuncO1sim(var_sim) + ws2*FuncO2sim(var_sim) + ws3*FunCO3sim(var_sim) 
  return (total) 
}

#DE
Pop_Inicial <- as.matrix(sim.m_aclhs)
vector_min =  apply(Pop_Inicial, 2, min)  
vector_max =  apply(Pop_Inicial, 2, max)  
vector_min[aclhs_sub48] = CO2_aclhs_sub48
vector_max[aclhs_sub48] = CO2_aclhs_sub48
# plot(Time, vector_min, "l", col = "red")
# par(new=TRUE)
# plot(Time, vector_max, "l")

system.time(
aclhs_sim_outDEoptim_CO2 <- DEoptim(fn=FuncConjsim, lower= vector_min,
                                      upper= vector_max,
                                      control = DEoptim.control(VTR = 0.0000001,strategy = 3, 
                                                                itermax =10000, reltol = 1e-8, 
                                                                CR = 0.5, F = 0.8, NP= nrow(Pop_Inicial), initialpop = Pop_Inicial)) #, , 
)
# user       system      elapsed 
# 1274.096   94.597 1353.162, 3 functions: 0.1+10+1, aclhs_sim_outDEoptim_CO2$optim$bestval = 3.542966
# 1487.643  135.245 1681.886, 3 functions: 1+1+1, aclhs_sim_outDEoptim_CO2$optim$bestval = 30.96094
# 856.976  14.504 879.880, m <- 10, itermax =10000, 3 functions: 0.1+10+1, aclhs_sim_outDEoptim_CO2$optim$bestval = 2.445662
# 935.065   79.601 1045.929 m <- 10, itermax =10000, 3 functions: 0.1+10+1, aclhs_sim_outDEoptim_CO2$optim$bestval = 2.79532
# 521.218  39.790 566.670, m <- 30, itermax =2000, 3 functions: 0.1+10+1, aclhs_sim_outDEoptim_CO2$optim$bestval = 3.712153
# 967.742   95.338 1080.286 m <- 10, itermax =10000, 3 functions: 0.1+100+1, aclhs_sim_outDEoptim_CO2$optim$bestval = 5.444891
# 902.744   92.822 1011.953 m <- 10, itermax =10000, 3 functions: 0.1+100+1, aclhs_sim_outDEoptim_CO2$optim$bestval = 5.408111
# 784.942  82.992 857.929 m <- 10, itermax =10000, 3 functions: 0.1+100+1, aclhs_sim_outDEoptim_CO2$optim$bestval = 14.59955 (weight for lags)
# 901.140  69.554 990.032 m <- 10, itermax =10000, 3 functions: 0.1+100+1, aclhs_sim_outDEoptim_CO2$optim$bestval = 10.43702 (weight for lags)
# 831.408  64.206 899.837 m <- 10, itermax =10000, 3 functions: 0.1+100+1, aclhs_sim_outDEoptim_CO2$optim$bestval = 5.368356 (weight for lags)
aclhs_sim_outDEoptim_CO2$optim$bestval
Ite <- seq(1,length(aclhs_sim_outDEoptim_CO2$member$bestvalit),1)
plot(Ite,aclhs_sim_outDEoptim_CO2$member$bestvalit, "l", xlab ="Iterations", ylab ="Objective Function")
var_sim <- as.numeric(aclhs_sim_outDEoptim_CO2$optim$bestmem)


#------------------- tables-------------------#
# Statistical Properties
CO2_fixed_pred <- fixed_sim_outDEoptim_CO2$optim$bestmem #sim.m_fixed[1,]
CO2_fixed_pred[fixed_sub48] <- CO2_fixed_sub48
CO2_clhs_pred <- clhs_sim_outDEoptim_CO2$optim$bestmem #sim.m_clhs[1,]
CO2_clhs_pred[clhs_sub48] <- CO2_clhs_sub48
CO2_aclhs_pred <- aclhs_sim_outDEoptim_CO2$optim$bestmem
CO2_aclhs_pred[aclhs_sub48] <- CO2_aclhs_sub48
CO2_fixed_pred_Stat <- Estadisticas(CO2_fixed_pred)
CO2_clhs_pred_Stat <- Estadisticas(CO2_clhs_pred)
CO2_aclhs_pred_Stat <- Estadisticas(CO2_aclhs_pred)
cbind(CO2_Stat,CO2_fixed_pred_Stat,CO2_clhs_pred_Stat,CO2_aclhs_pred_Stat)



ks.test(CO2_fixed_pred, CO2_fixed_sub48, conf.level = 0.95)
ks.test(CO2_clhs_pred, CO2_clhs_sub48, conf.level = 0.95)
ks.test(sim.m_aclhs[1,], CO2_aclhs_sub48, conf.level = 0.95)
ks.test(as.numeric(CO2_aclhs_pred), CO2_aclhs_sub48, conf.level = 0.95)

ks_test_fixed_CO2_sim <- ks.test(CO2_fixed_pred, CO2, conf.level = 0.95)
ks_test_clhs_CO2_sim <- ks.test(CO2_clhs_pred, CO2, conf.level = 0.95)
#ks_test_aclhs_CO2_sim1 <- ks.test(sim.m_aclhs[1,], CO2, conf.level = 0.95)
ks_test_aclhs_CO2_sim2 <- ks.test(as.numeric(CO2_aclhs_pred), CO2, conf.level = 0.95)
ks_test_fixed_CO2_sim
ks_test_clhs_CO2_sim
#ks_test_aclhs_CO2_sim1
ks_test_aclhs_CO2_sim2

# Data
(corP_Temp_CO2 <- round(cor(Temp,CO2, method = "pearson"),3))
(corS_Temp_CO2 <- round(cor(Temp,CO2, method = "spearman"),3))
(corK_Temp_CO2 <- round(cor(Temp,CO2, method = "kendall"),3))
# fixed
(corP_Temp_CO2_fixed_pred <- round(cor(Temp,CO2_fixed_pred, method = "pearson"),3))
(corS_Temp_CO2_fixed_pred <- round(cor(Temp,CO2_fixed_pred, method = "spearman"),3))
(corK_Temp_CO2_fixed_pred <- round(cor(Temp,CO2_fixed_pred, method = "kendall"),3))
# clhs
(corP_Temp_CO2_clhs_pred <- round(cor(Temp,CO2_clhs_pred, method = "pearson"),3))
(corS_Temp_CO2_clhs_pred <- round(cor(Temp,CO2_clhs_pred, method = "spearman"),3))
(corK_Temp_CO2_clhs_pred <- round(cor(Temp,CO2_clhs_pred, method = "kendall"),3))
#aclhs
(corP_Temp_CO2_aclhs_pred <- round(cor(Temp,CO2_aclhs_pred, method = "pearson"),3))
(corS_Temp_CO2_aclhs_pred <- round(cor(Temp,CO2_aclhs_pred, method = "spearman"),3))
(corK_Temp_CO2_aclhs_pred <- round(cor(Temp,CO2_aclhs_pred, method = "kendall"),3))

abs(corP_Temp_CO2_fixed_pred-corP_Temp_CO2)
abs(corS_Temp_CO2_fixed_pred-corS_Temp_CO2)
abs(corK_Temp_CO2_fixed_pred-corK_Temp_CO2)
abs(corP_Temp_CO2_clhs_pred-corP_Temp_CO2)
abs(corS_Temp_CO2_clhs_pred-corS_Temp_CO2)
abs(corK_Temp_CO2_clhs_pred-corK_Temp_CO2)
abs(corP_Temp_CO2_aclhs_pred-corP_Temp_CO2)
abs(corS_Temp_CO2_aclhs_pred-corS_Temp_CO2)
abs(corK_Temp_CO2_aclhs_pred-corK_Temp_CO2)

(sum(abs(CO2-as.numeric(CO2_fixed_pred)))/sum(CO2))*100
(sum(abs(CO2-as.numeric(CO2_clhs_pred)))/sum(CO2))*100
(sum(abs(CO2-as.numeric(CO2_aclhs_pred)))/sum(CO2))*100


sum(abs(CO2-as.numeric(CO2_fixed_pred)))
sum(abs(CO2-as.numeric(CO2_clhs_pred)))
sum(abs(CO2-as.numeric(CO2_aclhs_pred)))


#--------------------------------Figure 8--------------------------------------#
#n_logs= 4 #length(log_list)
plot_mean=TRUE
plot_median=TRUE
fontproportion=1.0
png(paste(Sims_dir,"/CO2_Temp_dist_Comparison.png",sep=""), 
    bg = "white", width = 2000, height = 1000, res = 150 )
par(mfrow = c(2,2),mar = c(4, 6, 2, 2), cex.lab= 1.0, cex.axis = 1.0)

# (a)
plot(Time, CO2,  type = "p", lwd = 2, pch = 1, 
     ylab = CO2label , xlab = Timelabel, bty="o")
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Time, CO2, type = "p", lwd = 2, pch = 1, 
     ylab =CO2label , xlab = Timelabel, bty="o")

# (b)
plot(Time, CO2,  type = "p", lwd = 2, pch = 1, xlim = c(min(Time), max(Time)), ylim = c(min(CO2), max(CO2)),
     ylab = CO2label , xlab = Timelabel, bty="o", col = "lightgray")
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Time, CO2, pch = 15, cex = 1.2, bty="o", col = "lightgray", xlim = c(min(Time), max(Time)), ylim = c(min(CO2), max(CO2)), xlab =Timelabel, ylab = CO2label)
par(new=TRUE)
plot(Time,CO2_fixed_pred, pch = 15, cex = 1.2, bty="o", col = "green", xlim = c(min(Time), max(Time)), ylim = c(min(CO2), max(CO2)), xlab =Timelabel, ylab = CO2label)

# (c)
plot(Time, CO2,  type = "p", lwd = 2, pch = 1, xlim = c(min(Time), max(Time)), ylim = c(min(CO2), max(CO2)),
     ylab = CO2label , xlab = Timelabel, bty="o", col = "lightgray")
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Time, CO2, pch = 17, cex = 1.2, bty="o", col = "lightgray", xlim = c(min(Time), max(Time)), ylim=c(min(CO2), max(CO2)), xlab =Timelabel, ylab = CO2label)
par(new=TRUE)
plot(Time,CO2_clhs_pred, pch = 17, cex = 1.2, bty="o", col = "red", xlim = c(min(Time), max(Time)), ylim=c(min(CO2), max(CO2)), xlab =Timelabel, ylab = CO2label)

# (d)
plot(Time, CO2,  type = "p", lwd = 2, pch = 1,
     ylab = CO2label , xlab = Timelabel, bty="o", col = "lightgray")
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(Time, CO2, pch = 18,  cex = 1.5,  bty="o", col = "lightgray", xlim = c(min(Time), max(Time)), ylim=c(min(CO2), max(CO2)), xlab =Timelabel, ylab = CO2label)
par(new=TRUE)
plot(Time,CO2_aclhs_pred, pch = 18,  cex = 1.5, bty="o", col = "blue", xlim = c(min(Time), max(Time)), ylim=c(min(CO2), max(CO2)), xlab =Timelabel, ylab = CO2label)

box()
dev.off()


#---------------------------------------------Figure9_Univariate pdfs--------------------------------#
#Fn_Temp <- estandarizar(cbind(Temp,CO2))[,1]
#Fn_CO2 <- estandarizar(cbind(Temp,CO2))[,2]
U_limit <- seq(0,1,by = 0.2)
png(paste(Sims_dir,"/CO2_univariate_pdfs_sim.png",sep=""),   bg = "white", width = 1200, height = 1000, res = 120)
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
plot(ecdf(CO2_fixed_pred),pch = 15, col = "green", cex = 1.2,
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
plot(ecdf(CO2_clhs_pred),pch = 17, col = "red", cex = 1.2,
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
plot(ecdf(CO2_aclhs_pred),pch = 18, col = "blue", cex = 1.5,
     yaxt = "n",  main = " ", xlab= CO2label, ylab = expression(bold(paste(Fn(CO[2] ~ efflux)))),
     xlim =c(min(CO2),max(CO2)), ylim = c(0,1))
box()
dev.off()



#---------------------------------------------Figure10_scatterplot_sims--------------------------------#
Fn_Temp <- estandarizar(cbind(Temp,CO2))[,1]
Fn_CO2 <- estandarizar(cbind(Temp,CO2))[,2]
U_limit <- seq(0,1,by = 0.2)
png(paste(Sims_dir,"/CO2_Temp_scatterplots_sims.png",sep=""),   bg = "white", width = 1200, height = 1000, res = 120)
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
plot(Temp, CO2_fixed_pred,  xlab= Templabel, ylab = CO2label,
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
plot(Temp, CO2_clhs_pred,  xlab= Templabel, ylab = CO2label,
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
plot(Temp, CO2_aclhs_pred,  xlab= Templabel, ylab = CO2label,
     pch = 18, col = "blue", cex = 1.5,
     xlim = c(min(Temp), max(Temp)),  ylim = c(min(CO2), max(CO2)))
box()
dev.off()


#------------------------------- Figure 8--------------------------------------#
Variolabel <- expression(bold("Semivariance"))
# Estimation of the experimental variogram
N_lags<- 8 #length(t)/2.5
lag_value <- 1 # min(dist(t)) # or delta t
ti = numeric(length(Time))
#png(paste(VA_dir,"/CO2_fixed_sim_VarioEstimation.png",sep=""), bg = "white", width = 1500, height = 1000, res = 200)
CO2_fixed_sim_VarioEstimation<-Variograma(Time, ti, 
                                          CO2_fixed_pred, 0, 90, N_lags, lag_value, 1, "", Timelabel) # CO2 temporal variogram
#dev.off()

#png(paste(VA_dir,"/CO2_clhs_sim_VarioEstimation.png",sep=""), bg = "white", width = 1500, height = 1000, res = 200)
CO2_clhs_sim_VarioEstimation<-Variograma(Time, ti, 
                                         CO2_clhs_pred, 0, 90, N_lags, lag_value, 1, "", Timelabel) # CO2 temporal variogram
#dev.off()

#png(paste(VA_dir,"/CO2_aclhs_sim_VarioEstimation.png",sep=""), bg = "white", width = 1500, height = 1000, res = 200)
CO2_aclhs_sim_VarioEstimation<-Variograma(Time, ti, 
                                          as.numeric(CO2_aclhs_pred), 0, 90, N_lags, lag_value, 1, "", Timelabel) # CO2 temporal variogram
#dev.off()

sum(abs(CO2_fixed_sim_VarioEstimation$Semivarianzas-CO2_VarioEstimation$Semivarianzas))
sum(abs(CO2_clhs_sim_VarioEstimation$Semivarianzas-CO2_VarioEstimation$Semivarianzas))
sum(abs(CO2_aclhs_sim_VarioEstimation$Semivarianzas-CO2_VarioEstimation$Semivarianzas))

png(paste(Sims_dir,"/Variogram_comparison_sims.png",sep=""), bg = "white",  width = 1500, height = 1000, res = 200)
par(mfrow = c(1,1),mar = c(4, 6, 2, 2), cex.lab= 1.2, cex.axis = 1.2)
plot(CO2_VarioEstimation[,c(2,3)], pch = 1, xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*3)),
     xlab = Timelabel ,ylab = Variolabel)
grid(col = "lightgray", lty = "dashed", lwd = par("lwd"), equilogs = TRUE)
par(new=TRUE)
plot(CO2_VarioEstimation[,c(2,3)], pch = 1, xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*3)),
     xlab = Timelabel ,ylab = Variolabel )
par(new=TRUE)
plot(CO2_fixed_sim_VarioEstimation[,c(2,3)], pch = 15, cex = 1.2, col = "green", xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*3)),
     xlab = Timelabel ,ylab = Variolabel)
par(new=TRUE)
plot(CO2_clhs_sim_VarioEstimation[,c(2,3)], pch = 17, cex = 1.2, col = "red",  xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*3)),
     xlab = Timelabel , ylab = Variolabel)
par(new=TRUE)
plot(CO2_aclhs_sim_VarioEstimation[,c(2,3)], pch = 18, cex = 1.5, col = "blue", xlim = c(0,9), ylim=c(0,max(CO2_VarioEstimation[,3]*3)),
     xlab = Timelabel , ylab = Variolabel)
par(new=TRUE)
lines.variomodel(cov.model = "sph", cov.pars = c(CO2_aclhs_sub48_sill_and_nugget-CO2_aclhs_sub48_nugget, CO2_aclhs_sub48_range)
                 , nug = CO2_aclhs_sub48_nugget, col = "Black", lwd = 3, max.dist = 1*8)
box()
dev.off()

