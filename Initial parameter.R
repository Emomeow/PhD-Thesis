## Initialize parameter
## sample size
n.id=400
## assume before t0, the baseline hazard=0
t0=0
# visit time centered at=0,1,2,3,...,J-1 (if not censored)
## visit label: t1=0, t2, ...,t_J=J-1
J=11

par=c(5,1.5)
# True values of parameters
# K domains for biomarkers
K=2
## p dims for covarite X
p=1
#### Define variance process model here
Sigma.x=0.8*diag(p)
Sigma.a=Sigma_CS(c(1:K), c(0.5,0.5,0.2))
Sigma.a.vec=matrix(Sigma.a[upper.tri(Sigma.a,diag = TRUE)])
pos_labels_a <- sapply(1:K, function(i) {
  sapply(i:K, function(j) paste0("Sigma.a_", i, "_", j))
})
pos_labels_a <- unlist(pos_labels_a)
rownames(Sigma.a.vec) <- pos_labels_a
Sigma.e=Sigma_CS(c(1:K), c(0.5,1,0.2))
Sigma.e.vec=matrix(Sigma.e[upper.tri(Sigma.e,diag = TRUE)])
pos_labels_e <- sapply(1:K, function(i) {
  sapply(i:K, function(j) paste0("Sigma.e_", i, "_", j))
})
pos_labels_e <- unlist(pos_labels_e)
rownames(Sigma.e.vec) <- pos_labels_e

## fix effect parameter in lm
## beta.0 intercept, beta.1 for slope (t), beta.2 for X
beta.0=matrix(rep(5,K),nrow=K)
rownames(beta.0)=paste("beta.0","_",1:K,sep="")
beta.1=matrix(rep(-0.2,K),nrow=K)
rownames(beta.1)=paste("beta.1","_",1:K,sep="")
beta.2=matrix(-0.3,nrow = K,ncol= p) 
index_grid <- expand.grid(row = 1:K, col = 1:p)
beta.2.vec=matrix(as.vector(beta.2))
rownames(beta.2.vec)=paste("beta.2",index_grid$row, index_grid$col,sep='_')

## parameter in hazard exponential
theta.x=matrix(rep(0.3,p),nrow=p)
rownames(theta.x)=paste("theta.x",1:p,sep="_")
## dim p*1 
theta.a=matrix(rep(0.1,K),nrow=K)
rownames(theta.a)=paste("theta.a",1:K,sep="_")
## theta.a=rep(0.2,K)
theta.m=matrix(rep(-0.2,K),nrow=K)
rownames(theta.m)=paste("theta.m",1:K,sep="_")

## parameter in the model for updating lm
## gamma.1 for slope (t), gamma.2 for X
gamma=matrix(-0.4,nrow=K)
rownames(gamma)=paste("gamma",1:K,sep="_")
# G-H
no.gauss = 20
out.gauss=gauss.quad(no.gauss, kind="hermite")
nodes=as.matrix(do.call(expand.grid,replicate(K,out.gauss$nodes,simplify = FALSE)))
if (K>1){
  weight=as.vector(do.call(outer,replicate(K,out.gauss$weights,simplify = FALSE)))
} else {
  weight=as.vector(out.gauss$weights)
}
# nodes_list  <- replicate(K, out.gauss$nodes, simplify = FALSE)
# weights_list <- replicate(K, out.gauss$weights, simplify = FALSE)
# 
# nodes  <- as.matrix(do.call(expand.grid, nodes_list))
# weight <- apply(do.call(expand.grid, weights_list), 1, prod)

## mode = 1 # constant hazard
# mode = 2 # piecewise constant
## mode = 3 # B-spline
if (mode == 1){

## baseline cumulative hazard parameter
lambda0 = matrix(1/par[1])
rownames(lambda0) = paste("lambda0")

## Combine all parameters into a dataframe
parameters = data.frame(t(beta.0),t(beta.1),t(beta.2.vec),t(gamma),t(lambda0),
                        t(theta.x),t(theta.a),t(theta.m),t(Sigma.a.vec),t(Sigma.e.vec))

}

##########################
if (mode == 2){
  knots = 3
  ## Combine all parameters into a dataframe
  parameters = data.frame(t(beta.0),t(beta.1),t(beta.2.vec),t(gamma),
                          t(theta.x),t(theta.a),t(theta.m),t(Sigma.a.vec),t(Sigma.e.vec))

}

##############################
if (mode == 3){
  num_knots = 6
  degree = 3
  parameters = data.frame(t(beta.0),t(beta.1),t(beta.2.vec),t(gamma),
                          t(theta.x),t(theta.a),t(theta.m),t(Sigma.a.vec),t(Sigma.e.vec))
  leg <- gauss.quad(no.gauss, kind = "legendre")
  z <- leg$nodes    # length = 20
  w <- leg$weights  # length = 20
}

estimate_bspline_hazard <- function(surv_obj, knots, degree) {
  
  # --- 1. Define the Negative Log-Likelihood Function ---
  # This is the function that `optim` will try to MINIMIZE.
  # It now models h(t) = S(t) with the constraint alpha > 0.
  neg_log_likelihood <- function(alpha, time1, time2, status, knots, degree) {
    total_ll <- 0
    # a. Calculate the log-hazard at the event times
    idx_right <- which(status == 0)
    if (length(idx_right) > 0) {
      times_L <- time1[idx_right]
      H_L <- spline_cumulative_hazard(times_L, knots, alpha, boundary_knots, degree)
      total_ll <- total_ll + sum(-H_L)
    }
    
    # idx_interval <- which(status == 1)
    # if (length(idx_interval) > 0) {
    #   # times_L <- time1[idx_interval]
    #   # times_R <- time2[idx_interval]
    #   times_M <- time[idx_interval]
    #   hazard <- spline_hazard(times_M, knots, alpha, boundary_knots, degree)
    #   H_M = spline_cumulative_hazard(times_M, knots, alpha, boundary_knots, degree)
    #   total_ll <- total_ll + sum(log(pmax(hazard,1e-9))) + sum(-H_M)
    # }
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
    # b. Calculate the cumulative hazard for ALL subjects at their follow-up time
    
    # We use sapply to integrate the hazard function from 0 to each subject's time
    # cumulative_hazard_values <- as.vector(ibs(times, knots = knots, degree = degree, intercept = TRUE) %*% alpha)
    #cumulative_hazard_values <- spline_cumulative_hazard(times, knots, alpha, boundary_knots, degree)
    # c. Calculate the total log-likelihood and negate it for minimization
    #log_likelihood <- sum(log_hazard_values) - sum(cumulative_hazard_values)
    
    # Return a large value if the log-likelihood is not finite, to guide the optimizer
    # if (!is.finite(log_likelihood)) {
    #   return(1e9)
    # }
    
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
    # lower = rep(-1, num_alphas), # Set a small positive lower bound for alpha
    # upper = rep(1, num_alphas),
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

estimate_mspline_hazard <- function(surv_obj, knots, degree) {
  
  # --- 1. Define the Negative Log-Likelihood Function ---
  # This is the function that `optim` will try to MINIMIZE.
  # It now models h(t) = S(t) with the constraint alpha > 0.
  neg_log_likelihood <- function(alpha, time1, time2, status, knots, degree) {
    total_ll <- 0
    # a. Calculate the log-hazard at the event times
    idx_right <- which(status == 0)
    if (length(idx_right) > 0) {
      times_L <- time1[idx_right]
      H_L <- evaluate_Ispline(times_L, knots, alpha, boundary_knots, degree)
      total_ll <- total_ll + sum(-H_L)
    }
    
    idx_interval <- which(status == 3)
    if (length(idx_interval) > 0) {
      times_L <- time1[idx_interval]
      times_R <- time2[idx_interval]
      
      H_L <- evaluate_Ispline(times_L, knots, alpha, boundary_knots, degree)
      H_R <- evaluate_Ispline(times_R, knots, alpha, boundary_knots, degree)
      
      S_L <- exp(-H_L)
      S_R <- exp(-H_R)
      
      total_ll <- total_ll + sum(log(S_L - S_R))
    }
    # b. Calculate the cumulative hazard for ALL subjects at their follow-up time
    
    # We use sapply to integrate the hazard function from 0 to each subject's time
    # cumulative_hazard_values <- as.vector(ibs(times, knots = knots, degree = degree, intercept = TRUE) %*% alpha)
    #cumulative_hazard_values <- spline_cumulative_hazard(times, knots, alpha, boundary_knots, degree)
    # c. Calculate the total log-likelihood and negate it for minimization
    #log_likelihood <- sum(log_hazard_values) - sum(cumulative_hazard_values)
    
    # Return a large value if the log-likelihood is not finite, to guide the optimizer
    # if (!is.finite(log_likelihood)) {
    #   return(1e9)
    # }
    
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
    # lower = rep(0, num_alphas), # Set a small positive lower bound for alpha
    # upper = rep(1, num_alphas),
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