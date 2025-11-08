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
library(expm)
library(survival)

setwd("C:/Users/yuezh/OneDrive - University of Nebraska Medical Center/Research/Joint_model/joint_likelihood")
source('Data generation.R')
source('Numeric Diff.R')
mode = 3
source('Initial parameter.R')
# source('Causal Mediation.R')
source('NR.R')


args = commandArgs(trailingOnly=TRUE)
i=as.numeric(args[1])
# 
iniseed=1111+i


set.seed(iniseed)


col_names <- colnames(parameters)
num_para <- length(parameters)


long.data = data.gen.vec(parameters,n.id,J,par,Sigma.x)
#################################
# Give initial guess of B-spline cumulative hazard
long.data.id = long.data[!duplicated(long.data$id),]
end_times <- pmin(long.data.id$etime, long.data.id$ctime)

# Combine all end times with only the start times that are different
time <- c(end_times, long.data.id$preetime[long.data.id$preetime != end_times])

knots = quantile(time, probs = seq(1, num_knots)/(num_knots+1), na.rm = TRUE)

# set boundary knots
boundary_knots = c(0, max(max(long.data$visittime, na.rm = TRUE), knots))

# give initial guess of alpha
surv_data <- long.data.id %>%
  mutate(U = if_else(status == 1, etime, Inf))

surv_object <- Surv(time = surv_data$preetime, time2 = surv_data$U, type = "interval2")

# boundary_knots = c(0, max(long.data$visittime))
estimate_bspline_hazard <- function(surv_obj, knots, degree) {
  
  # --- 1. Define the Negative Log-Likelihood Function ---
  # This is the function that `optim` will try to MINIMIZE.
  # It now models h(t) = S(t) with the constraint alpha > 0.
  neg_log_likelihood <- function(alpha, time1, time2, status, knots, degree) {
    total_ll <- 0
    # a. Calculate the log likelihood for right-censored data
    idx_right <- which(status == 0)
    if (length(idx_right) > 0) {
      times_L <- time1[idx_right]
      H_L <- spline_cumulative_hazard(times_L, knots, alpha, boundary_knots, degree)
      total_ll <- total_ll + sum(-H_L)
    }
    # b. Calculate the log likelihood for interval-censored data
    idx_interval <- which(status == 3)
    if (length(idx_interval) > 0) {
      times_L <- time1[idx_interval]
      times_R <- time2[idx_interval]
      
      H_L <- spline_cumulative_hazard(times_L, knots, alpha, boundary_knots, degree)
      H_R <- spline_cumulative_hazard(times_R, knots, alpha, boundary_knots, degree)
      
      S_L <- exp(-H_L)
      S_R <- exp(-H_R)
      
      total_ll <- total_ll + sum(log(S_L - S_R))
    }
    
    return(-total_ll)
  }
  
  # --- 2. Set Up and Run the Optimization ---
  
  # The number of alpha coefficients needed is determined by the spline parameters
  num_alphas <- degree + length(knots) + 1
  
  # Provide an initial guess for the coefficients
  initial_alphas <- rnorm(num_alphas, 0, 1)
  
  # Run the optimizer to find the alphas that minimize the negative log-likelihood
  optim_results <- optim(
    par = initial_alphas,
    fn = neg_log_likelihood,
    method = "L-BFGS-B", # Use a method that allows box constraints
    hessian = TRUE, # Ask for the Hessian matrix to calculate standard errors
    time1 = surv_obj[, 1],
    time2 = surv_obj[, 2],
    status = surv_obj[, 3],
    knots = knots,
    degree = degree
  )
  
  # --- 3. Format and Return the Results ---
  
  # Extract the estimated coefficients
  estimated_alphas <- optim_results$par
  
  # Calculate standard errors from the inverse of the Hessian matrix
  if (optim_results$convergence == 0) {
    fisher_info <- optim_results$hessian
    # Check for positive definiteness before inverting
    if (all(eigen(fisher_info)$values > 1e-6)) {
      vcov_matrix <- solve(fisher_info)
      std_errors <- sqrt(diag(vcov_matrix))
    } else {
      warning("Hessian matrix is not positive definite. Standard errors are unreliable.")
      std_errors <- rep(NA, num_alphas)
    }
  } else {
    warning("Optimization did not converge. Standard errors are not reliable.")
    std_errors <- rep(NA, num_alphas)
  }
  
  
  coefficients_df <- data.frame(
    alpha = paste0("alpha_", 1:num_alphas),
    estimate = estimated_alphas,
    std_error = std_errors
  )
  
  return(list(
    coefficients = coefficients_df,
    optim_results = optim_results,
    knots = knots,
    degree = degree
  ))
}

hazard_fit <- estimate_bspline_hazard(surv_object, knots, degree)

alpha = data.frame(t(hazard_fit[["coefficients"]][["estimate"]]))
colnames(alpha) = paste("alpha",1:(num_knots+degree+1),sep = "_")
# test 
# test.time = seq(0, max(long.data$visittime), 0.01)
# 
# test.spline = spline_cumulative_hazard(test.time, knots, t(alpha), boundary_knots, degree)
# plot(test.time, test.spline)
# lines(test.time, (test.time/par[1])^par[2])
# test.hazard = spline_hazard(test.time, knots, t(alpha), boundary_knots, degree)
# plot(test.time, test.hazard)
# lines(test.time, par[2]/par[1]*(test.time/par[1])^(par[2]-1))

# ini.parameters = cbind(parameters, alpha)
ini.parameters = cbind(parameters + runif(num_para, min=-0.5,max=0.5)*parameters, alpha)
###################################

hat.parameters = ini.parameters

# ini.likelihood = likelihood.spline2(long.data, ini.parameters, knots)$ll

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

knots = data.frame(t(knots))
colnames(knots) = paste("knots_",1:num_knots,sep = '')
alpha.hat = hat.parameters %>% dplyr::select(contains('alpha'))
boundary = data.frame(boundary_knots[2])
colnames(boundary) = paste('boundary')
hazard = cbind(knots, alpha.hat, boundary)
write.csv(hazard,paste("hazard_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)

# bootstrap for stddev
B = 50
bootstrap = 1

# Bootstrap code
if (bootstrap == 1){
  coef_bootstrap = matrix(NA, nrow = B, ncol = ncol(coef))
  
  for (b in 1:B){
    bootstrap_data = bootstrap_sample(long.data)
    bootstrap_data$id = rep(1:n.id, each = J)
    long.data.id = bootstrap_data[!duplicated(bootstrap_data$id),]
    time = c(long.data.id$preetime,pmin(long.data.id$etime,long.data.id$ctime))
    
    #d = c(0, quantile(time, probs = seq(1, knots)/knots, na.rm = TRUE))
    #ini.hazard = data.frame(t(1/knots/diff(d)))
    #colnames(ini.hazard) = paste("hazard",1:knots,sep = "_")
    knots = quantile(time, probs = seq(1, num_knots)/(num_knots+1), na.rm = TRUE)
    
    # set boundary knots
    boundary_knots = c(0, max(time, na.rm = TRUE))
    
    # give initial guess of alpha
    surv_data <- long.data.id %>%
      mutate(time_cox = (preetime + pmin(etime,ctime))/2)
    
    surv_object <- Surv(time = surv_data$time_cox, event = surv_data$status)
    
    hazard_fit <- estimate_bspline_hazard(surv_object, knots, degree)
    
    alpha = data.frame(t(hazard_fit[["coefficients"]][["estimate"]]))
    colnames(alpha) = paste("alpha",1:(num_knots+degree+1),sep = "_")
    
    hat.parameters =  cbind(coef, alpha)
    
    #Calculate MLE for bootstrap data
    hat.parameters <- NR_spline(bootstrap_data, hat.parameters, knots)
    hat.parameters <- as.matrix(hat.parameters)
    coef_bootstrap[b, ] <- hat.parameters[, 1:length(parameters)]
  }
  stddev = t(as.matrix(apply(coef_bootstrap, 2, sd,na.rm=TRUE)))
}

if (bootstrap == 0) {
  # stddev = t(as.matrix(Est_sd(long.data, hat.parameters)))
  Hessian.MLE = diff.likelihood.vec(long.data, hat.parameters)$Hessian
  stddev = data.frame(t(sqrt(-diag(Hessian.MLE))))
  stddev = stddev[,1:num_para]
}

# Hessian.MLE = diff.likelihood.vec(long.data, hat.parameters)$Hessian
# stddev = data.frame(t(sqrt(-diag(Hessian.MLE))))
# stddev = stddev[,-((length(stddev)-knots+1):length(stddev))]
colnames(stddev) = col_names
CP = as.matrix((parameters >= coef-1.96*stddev)&(parameters <= coef+1.96*stddev))
colnames(CP) = col_names

write.csv(stddev,paste("est_std_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)
write.csv(CP,paste("CP_i=",as.character(i),"K=",K,"p=",p,"n=",n.id,"seed=",iniseed,".csv",sep=""),row.names=FALSE)
