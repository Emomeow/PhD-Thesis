# B spline
######################################
# import packages and generate data
library(splines)
library(splines2)
library(MASS)
library(dplyr)
library(tidyr)
library(statmod)
library(Matrix)
library(matrixcalc)
# library(expm)
library(survival)

setwd("C:/Users/yuezh/OneDrive - University of Nebraska Medical Center/Research/Joint_model/joint_likelihood")
source('Data generation.R')
source('Numeric Diff.R')
mode = 3
source('Initial parameter.R')
source('Causal Mediation.R')
source('NR.R')
source('Dynamic Prediction.R')

# args = commandArgs(trailingOnly=TRUE)
# i=as.numeric(args[1])
# 
# iniseed=1111+i


set.seed(1234)


col_names <- colnames(parameters)
num_para <- length(parameters)


long.data = data.gen.vec(parameters,n.id,J,par,Sigma.x)
#################################
# Give initial guess of B-spline cumulative hazard
long.data.id = long.data[!duplicated(long.data$id),]
mean(long.data.id$status)
end_times <- pmin(long.data.id$etime, long.data.id$ctime)

# Combine all end times with only the start times that are different
time <- c(end_times, long.data.id$preetime[long.data.id$preetime != end_times])

knots = quantile(time, probs = seq(1, num_knots)/(num_knots+1), na.rm = TRUE)

# set boundary knots
# boundary_knots = c(0, max(max(long.data$visittime, na.rm = TRUE), knots))
boundary_knots = c(0, J-0.8)

# give initial guess of alpha
surv_data <- long.data.id %>%
  mutate(U = if_else(status == 1, etime, Inf))
# surv_data <- long.data.id %>%
#   mutate(time_cox = (preetime+pmin(etime,ctime))/2)
surv_object <- Surv(time = surv_data$preetime, time2 = surv_data$U, type = "interval2")
# surv_object <- Surv(time = surv_data$time_cox, event = surv_data$status, type = 'right')
# boundary_knots = c(0, max(long.data$visittime))

#hazard_fit <- estimate_bspline_hazard(surv_object, knots, degree)
hazard_fit <- estimate_mspline_hazard(surv_object, knots, degree)

alpha = data.frame(t(hazard_fit[["coefficients"]][["estimate"]]))
colnames(alpha) = paste("alpha",1:(num_knots+degree+1),sep = "_")
# test 
# test.time = seq(0, max(long.data$visittime), 0.01)
# 
# test.spline = evaluate_Ispline(test.time, knots, t(alpha), boundary_knots, degree)
# plot(test.time, test.spline)
# lines(test.time, (test.time/par[1])^par[2])
# test.hazard = evaluate_Mspline(test.time, knots, t(alpha), boundary_knots, degree)
# plot(test.time, test.hazard)
# lines(test.time, par[2]/par[1]*(test.time/par[1])^(par[2]-1))

# ini.parameters = cbind(parameters, alpha)
ini.parameters = cbind(parameters + runif(num_para, min=-0.2,max=0.2)*parameters, alpha)
###################################
## Give initial guess of other parameters
# library(lme4)
# library(JMbayes2)
# library(survival)
# # Fit 2 longitudinal models for 2 biomarkers
# data_longitudinal <- long.data |>
#   group_by(id) |>
#   mutate(gammaterm = as.numeric(visittime >= etime)*(visittime - (etime+preetime)/2)) |>
#   ungroup()
# 
# model1_long = lme(fixed = y ~ visittime + gammaterm + x,
#                   random = ~ 1|id,
#                   data = data_longitudinal)
# 
# model2_long = lme(fixed = stroopwo ~ visittime + gammaterm + age1 + educ_yrs + gender + base_cap,
#                   random = ~ 1|id,
#                   data = data_longitudinal)
# 
# summary(model1_long)
# summary(model2_long)

hat.parameters = ini.parameters
#############################
# test
# mode = 3
# source('Initial parameter.R')
# degree <- 1
# knots <- NULL
# num_knots = 0
# alpha <- c(-Inf, log(boundary_knots[2]/par[1]))
# alpha <- data.frame(t(alpha))
# colnames(alpha) = paste("alpha",1:(num_knots+degree+1),sep = "_")
# ini.parameters = cbind(parameters, alpha)
# 
# beta.0.test = seq(-0.2,0.2,by=0.01)
# likelihood.test = rep(0,length(beta.0.test))
# 
# for (i in 1:length(beta.0.test)){
#   parameters.test.i = ini.parameters + beta.0.test[i]*t(c(1,rep(0, length(ini.parameters)-1)))
#   likelihood.test[i]=likelihood.spline(long.data, parameters.test.i, knots)$ll
# }
# plot(beta.0.test,likelihood.test)
# 
# mode = 2
# source('Initial parameter.R')
# long.data.id = long.data[!duplicated(long.data$id),]
# time = c(long.data.id$preetime,pmin(long.data.id$etime,long.data.id$ctime))
# 
# d = c(0, quantile(time, probs = seq(1, knots)/knots, na.rm = TRUE))
# ini.hazard = data.frame(t(rep(1/par[1], knots)))
# colnames(ini.hazard) = paste("hazard",1:knots,sep = "_")
# ini.parameters = parameters # + runif(num_para, min=-0.5,max=0.5)*parameters
# ini.parameters = cbind(ini.parameters, ini.hazard)
# 
# beta.0.test = seq(-0.2,0.2,by=0.01)
# likelihood.test = rep(0,length(beta.0.test))
# 
# for (i in 1:length(beta.0.test)){
#   parameters.test.i = ini.parameters + beta.0.test[i]*t(c(1,rep(0, length(ini.parameters)-1)))
#   likelihood.test[i]=likelihood.piecewise(long.data, parameters.test.i)$ll
# }
# plot(beta.0.test,likelihood.test)
# likelihood has some problem
ini.likelihood = likelihood.spline3(long.data, ini.parameters, knots)$ll
##################################

hat.parameters = NR_spline(long.data, hat.parameters, knots)
coef = hat.parameters %>% dplyr::select(!contains("alpha"))
colnames(coef) = col_names

# Hessian.MLE = diff.likelihood.vec(long.data, hat.parameters)$Hessian
# stddev = data.frame(t(sqrt(-diag(Hessian.MLE))))
# stddev = stddev[,1:num_para]
# 
# colnames(stddev) = col_names
# CP = as.matrix((parameters >= coef-1.96*stddev)&(parameters <= coef+1.96*stddev))
# 
# colnames(coef) = col_names
# colnames(stddev) = col_names
# colnames(CP) = col_names

write.csv(coef,paste("est_coef_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)

# write.csv(bias,paste("est_bias_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)

knots_i = data.frame(t(knots))
colnames(knots_i) = paste("knots_",1:num_knots,sep = '')
alpha.hat = hat.parameters %>% dplyr::select(contains('alpha'))
#boundary = data.frame(boundary_knots[2])
#colnames(boundary) = paste('boundary')
hazard = cbind(knots_i, alpha.hat)
write.csv(hazard,paste("hazard_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)

# bootstrap for stddev
B = 2  #test
bootstrap = 1



#hat.hazard = hat.parameters %>% dplyr::select(contains("hazard"))
#piecewise_hazard = cbind(as.vector(d[-length(d)]),t(hat.hazard))
#colnames(piecewise_hazard) = c("knots","hazard")
#coef = hat.parameters %>% dplyr::select(!contains("hazard"))
# bias = coef-parameters
#colnames(coef) = col_names
t = 6
x1 = 0
x2 = 1
index = 1
NE_S = With_cov(long.data, t, x1, x2, hat.parameters, index)
NDE_S = matrix(NE_S$NDE)
rownames(NDE_S) = paste("NDE")
NIE_S = matrix(NE_S$NIE)
rownames(NIE_S) = paste("NIE")
causal_S = data.frame(t(NDE_S),t(NIE_S))
# 
# a.re = long.data.id %>% dplyr::select(contains("aval"))
true_NE_S = With_cov_true(long.data, t, x1, x2, parameters, index, par)
true_NDE_S = matrix(true_NE_S$NDE)
rownames(true_NDE_S) = paste("NDE")
true_NIE_S = matrix(true_NE_S$NIE)
rownames(true_NIE_S) = paste("NIE")
true_NE_S = data.frame(t(true_NDE_S),t(true_NIE_S))


write.csv(causal_S, paste("causal_S_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)
write.csv(true_NE_S, paste("true_NE_S_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)

# M2
NE_M = NE_M2(long.data, t, x1, x2, hat.parameters, index)
NDE_M = NE_M$NDE
colnames(NDE_M) = paste("NDE",1:K,sep='')
NIE_M = NE_M$NIE
colnames(NIE_M) = paste("NIE",1:K,sep='')
causal_M = data.frame(NDE_M,NIE_M)

true_NE_M = NE_M2_true(long.data, t, x1, x2, parameters, index, par)
true_NDE_M = true_NE_M$NDE
colnames(true_NDE_M) = paste("NDE",1:K, sep='')
true_NIE_M = true_NE_M$NIE
colnames(true_NIE_M) = paste("NIE",1:K, sep='')
true_NE_M = data.frame(true_NDE_M,true_NIE_M)

write.csv(causal_M, paste("causal_M_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)
write.csv(true_NE_M, paste("true_NE_M_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)

# Bootstrap code
if (bootstrap == 1){
  coef_bootstrap = matrix(NA, nrow = B, ncol = ncol(coef))
  Bootstrap_NE_S = matrix(NA, nrow = B, ncol = 2)
  Bootstrap_NE_M = matrix(NA, nrow = B, ncol = 2*K)
  for (b in 1:B){
    #try({
    bootstrap_data = bootstrap_sample(long.data)
    bootstrap_data$id = rep(1:n.id, each = J)
    long.data.id = bootstrap_data[!duplicated(bootstrap_data$id),]
    end_times <- pmin(long.data.id$etime, long.data.id$ctime)
    
    # Combine all end times with only the start times that are different
    time <- c(end_times, long.data.id$preetime[long.data.id$preetime != end_times])
    # time = c(long.data.id$preetime,pmin(long.data.id$etime,long.data.id$ctime))
    
    knots = quantile(time, probs = seq(1, num_knots)/(num_knots+1), na.rm = TRUE)
    
    # set boundary knots
    boundary_knots = c(0, J-0.8)
    
    # give initial guess of alpha
    # surv_data <- long.data.id %>%
    # mutate(U = if_else(status == 1, etime, Inf))
    
    #surv_object <- Surv(time = surv_data$preetime, time2 = surv_data$U, type = "interval2")
    
    #hazard_fit <- estimate_bspline_hazard(surv_object, knots, degree)
    
    #alpha = data.frame(t(hazard_fit[["coefficients"]][["estimate"]]))
    #colnames(alpha) = paste("alpha",1:(num_knots+degree+1),sep = "_")
    
    #hat.parameters =  cbind(coef, alpha)
    
    #Calculate MLE for bootstrap data
    bootstrap.parameters <- NR_spline(bootstrap_data, hat.parameters, knots)
    Bootstrap_NE_S_i = With_cov(bootstrap_data, t, x1, x2, bootstrap.parameters, index)
    Bootstrap_NE_S[b,1] = Bootstrap_NE_S_i$NDE
    Bootstrap_NE_S[b,2] = Bootstrap_NE_S_i$NIE
    Bootstrap_NE_M_i = NE_M2(bootstrap_data, t, x1, x2, bootstrap.parameters, index)
    Bootstrap_NE_M[b,1:K] = Bootstrap_NE_M_i$NDE
    Bootstrap_NE_M[b,(K+1):(2*K)] = Bootstrap_NE_M_i$NIE
    bootstrap.parameters <- as.matrix(bootstrap.parameters)
    coef_bootstrap[b, ] <- bootstrap.parameters[, 1:length(parameters)]
    #})
  }
  # stddev = t(as.matrix(apply(coef_bootstrap, 2, sd,na.rm=TRUE)))
}

if (bootstrap == 0) {
  # stddev = t(as.matrix(Est_sd(long.data, hat.parameters)))
  Hessian.MLE = diff.likelihood.vec(long.data, hat.parameters)$Hessian
  stddev = data.frame(t(sqrt(-diag(Hessian.MLE))))
  stddev = stddev[,1:num_para]
}

stddev = data.frame(t(apply(coef_bootstrap, 2, sd, na.rm=TRUE)))

colnames(stddev) = col_names
CP = as.matrix((parameters >= coef-1.96*stddev)&(parameters <= coef+1.96*stddev))
colnames(CP) = col_names

write.csv(stddev,paste("est_std_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)
write.csv(CP,paste("CP_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)

sd_NDE_S = matrix(sd(Bootstrap_NE_S[, 1]))
sd_NIE_S = matrix(sd(Bootstrap_NE_S[, 2]))
rownames(sd_NDE_S) = paste("NDE")
rownames(sd_NIE_S) = paste("NIE")
sd_NE_S = data.frame(t(sd_NDE_S),t(sd_NIE_S))

CP_NDE_S = as.matrix((NDE_S >= true_NDE_S-1.96*sd_NDE_S)&(NDE_S <= true_NDE_S+1.96*sd_NDE_S))
CP_NIE_S = as.matrix((NIE_S >= true_NIE_S-1.96*sd_NIE_S)&(NIE_S <= true_NIE_S+1.96*sd_NIE_S))
rownames(CP_NDE_S) = paste("NDE")
rownames(CP_NIE_S) = paste("NIE")
CP_NE_S = data.frame(t(CP_NDE_S),t(CP_NIE_S))

sd_NDE_M = matrix(apply(Bootstrap_NE_M[, 1:K], 2, sd, na.rm=TRUE))
sd_NIE_M = matrix(apply(Bootstrap_NE_M[, (K+1):(2:K)], 2, sd, na.rm=TRUE))
rownames(sd_NDE_M) = paste("NDE",1:K, sep='')
rownames(sd_NIE_M) = paste("NIE",1:K, sep='')
sd_NE_M = data.frame(t(sd_NDE_M),t(sd_NIE_M))

CP_NDE_M = as.matrix((NDE_M >= true_NDE_M-1.96*t(sd_NDE_M))&(NDE_M <= true_NDE_M+1.96*t(sd_NDE_M)))
CP_NIE_M = as.matrix((NIE_M >= true_NIE_M-1.96*t(sd_NIE_M))&(NIE_M <= true_NIE_M+1.96*t(sd_NIE_M)))
colnames(CP_NDE_M) = paste("NDE",1:K, sep='')
colnames(CP_NIE_M) = paste("NIE",1:K, sep='')
CP_NE_M = data.frame(CP_NDE_M,CP_NIE_M)

write.csv(sd_NE_S, paste("sd_NE_S_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)
write.csv(CP_NE_S, paste("CP_NE_S_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)

write.csv(sd_NE_M, paste("sd_NE_M_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)
write.csv(CP_NE_M, paste("CP_NE_M_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)

#########################################
## Dynamic Prediction
t = 6
data.test = data.gen.vec(parameters, 5, J, par, Sigma.x)
CP_DP = vector("list", 5)
for (i in 1:5){
  long.data.i = data.test[data.test$id == i, ]
  CP.i = DP(long.data.i, t, hat.parameters)
  CP_DP[[i]] <- CP.i
}

CP_DP = do.call(rbind, CP_DP)
###############################
# doparallel
library(doParallel)
library(foreach)
# --- 1. Set up the Parallel Backend ---
# Use one less than your total cores to keep your computer responsive
n_cores <- detectCores() - 1 

# Create a cluster
my_cluster <- makeCluster(n_cores) 

# Register the cluster for 'foreach' to use
registerDoParallel(my_cluster)

# --- 2. Run the Parallel Loop ---
# We set a seed for reproducible results in parallel
RNGkind("L'Ecuyer-CMRG") # Use a parallel-safe random number generator

# The 'coef_bootstrap' matrix is now created by foreach
coef_bootstrap <- foreach(
  b = 1:B, 
  .combine = 'rbind', # Combine results by row-binding them
  
  # List all packages your loop code depends on
  .packages = c('dplyr', 'survival','splines','splines2','Matrix','statmod','expm','matrixcalc','MASS'), 
  
  # List all functions and variables from your main environment that
  # the workers need to access.
  .export = c('bootstrap_sample', 'n.id', 'J', 'num_knots', 'degree',
              'estimate_bspline_hazard', 'coef', 'alpha', 'NR_spline', 
              'parameters', 'long.data') 
) %dopar% {
  
  # --- This is the code from *inside* your loop ---
  # (No changes are needed to the logic)
  try({
  bootstrap_data = bootstrap_sample(long.data)
  bootstrap_data$id = rep(1:n.id, each = J)
  long.data.id = bootstrap_data[!duplicated(bootstrap_data$id),]
  time = c(long.data.id$preetime,pmin(long.data.id$etime,long.data.id$ctime))
  
  knots = quantile(time, probs = seq(1, num_knots)/(num_knots+1), na.rm = TRUE)
  
  boundary_knots = c(0, max(time, na.rm = TRUE))
  
  # surv_data <- long.data.id %>%
  #   mutate(U = if_else(status == 1, etime, Inf))
  # 
  # surv_object <- Surv(time = surv_data$preetime, time2 = surv_data$U, type = "interval2")
  # 
  # hazard_fit <- estimate_bspline_hazard(surv_object, knots, degree)
  # 
  # alpha = data.frame(t(hazard_fit[["coefficients"]][["estimate"]]))
  # colnames(alpha) = paste("alpha",1:(num_knots+degree+1),sep = "_")
  # 
  # hat.parameters =  cbind(coef, alpha)
  
  bootstrap.parameters <- NR_spline(bootstrap_data, hat.parameters, knots)
  boostrap.parameters <- as.matrix(bootstrap.parameters)
  
  # --- Return the result for this iteration ---
  # 'foreach' automatically combines this output
  bootstrap.parameters[, 1:length(parameters)]
  })
} # --- End of %dopar% loop ---

stopCluster(my_cluster)

stddev = data.frame(t(apply(coef_bootstrap, 2, sd)))
# Hessian.MLE = diff.likelihood.vec(long.data, hat.parameters)$Hessian
# stddev = data.frame(t(sqrt(-diag(Hessian.MLE))))
# stddev = stddev[,-((length(stddev)-knots+1):length(stddev))]
colnames(stddev) = col_names
CP = as.matrix((parameters >= coef-1.96*stddev)&(parameters <= coef+1.96*stddev))
colnames(CP) = col_names

write.csv(stddev,paste("est_std_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)
write.csv(CP,paste("CP_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)

##########################################

# # test likelihood for piecewise and Bspline
# degree <- 1
# knots <- NULL
# num_knots = 0
# # alpha = rep(1/par[1], num_knots+degree+1)
# 
# 
# #use cumsum(exp(alpha))
# boundary_knots = c(0, max(long.data$visittime))
# alpha <- c(-Inf, log(boundary_knots[2]/par[1]))
# hazard_value2 = spline_hazard(test.time, knots, alpha, boundary_knots, degree)
# #hazard_value_matrix = evaluate_spline_matrix(test.time, knots, alpha, boundary_knots, degree)
# cumulative_hazard_value2 = spline_cumulative_hazard(test.time, knots, alpha, boundary_knots, degree)
# 
# alpha <- data.frame(t(alpha))
# colnames(alpha) = paste("alpha",1:(num_knots+degree+1),sep = "_")
# ini.parameters = cbind(parameters, alpha)
# hat.parameters = ini.parameters
# ll1_2 = likelihood.spline(long.data, hat.parameters, knots)$ll
# ll1_3 = likelihood.spline2(long.data, hat.parameters, knots)$ll
# mode = 2
# source('Initial parameter.R')
# d = c(0, quantile(time, probs = seq(1, knots)/knots, na.rm = TRUE))
# # ini.hazard = data.frame(t(1/knots/diff(d)))
# ini.hazard = data.frame(t(rep(1/par[1], knots)))
# colnames(ini.hazard) = paste("hazard",1:knots,sep = "_")
# ini.parameters = cbind(parameters, ini.hazard)
# hat.parameters = ini.parameters
# ll2 = likelihood.piecewise(long.data, hat.parameters)$ll
# 
# mode = 1
# source('Initial parameter.R')
# hat.parameters = unlist(parameters)
# ll3 = likelihood.vec2(long.data, hat.parameters)

# Now the three codes in constant hazard case are the same.

# verify that fitted hazard close to true
