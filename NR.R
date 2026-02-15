NR = function(long.data, hat.parameters, tol = 1e-3, maxloop = 200){
  #long.data.id = long.data[!duplicated(long.data$id),]
  #time = c(long.data.id$preetime,pmin(long.data.id$etime,long.data.id$ctime))
  
  #d = c(0, quantile(time, probs = seq(1, knots)/knots, na.rm = TRUE))
  
  error = 0
  diff.theta = 1000
  # tol = 1e-3
  loop = 1
  # maxloop = 200
  while (diff.theta>tol & loop<=maxloop) {
    old.parameters = hat.parameters
    error=0
    #while (1) {
    diff = diff.likelihood.vec(long.data, hat.parameters)
    score = diff$Score
    lr = 1
    # Hessian = diff$Hessian
    Hessian = diff$Hessian
    if (any(is.na(Hessian))||any(is.infinite(Hessian))){
      error = 1
    }
    step = as.vector(Hessian%*%score)
    likelihood.hat = likelihood.piecewise(long.data, hat.parameters)$ll
    while (1){
      est.parameters = hat.parameters - lr*t(step)
      hazard.est = as.numeric(est.parameters %>% dplyr::select(contains("hazard")))
      sigma.e.est = est.parameters %>% dplyr::select(contains(paste("Sigma.e",1:K,1:K,sep="_")))
      sigma.a.est = est.parameters %>% dplyr::select(contains(paste("Sigma.a",1:K,1:K,sep="_")))
      if (any(hazard.est<=0)|any(sigma.e.est<=0)|any(sigma.a.est<=0)){
        lr = lr/2
        next
      }
      likelihood.est = likelihood.piecewise(long.data, est.parameters)$ll
      # if (any(is.na(likelihood.est)) || any(is.infinite(likelihood.est))) {
      #   lr = lr/2
      #   next
      # }
      if (lr*max(abs(step))<1e-6){
        temp.parameters = hat.parameters
        hat.parameters = est.parameters
        break
      }
      if (likelihood.est<likelihood.hat){
        lr = lr/2
        next
      }
      else {
        temp.parameters = hat.parameters
        hat.parameters = est.parameters
        break
      }
    }
    # if (max(abs(temp.parameters-hat.parameters))<tol) {
    #   break
    # }
    #}
    if (error == 1){
      hat.parameters = matrix(NA, nrow=1, ncol=length(temp.parameters))
      coef = matrix(NA, nrow=1, ncol=length(temp.parameters))
      stddev = matrix(NA, nrow=1, ncol=length(parameters))
      CP = matrix(NA, nrow=1, ncol=length(parameters))
      break
    }
    diff.theta = max(abs(hat.parameters-old.parameters))
    loop = loop + 1
  }
  return(hat.parameters)
}


bootstrap_sample <- function(data) {
  boot_ids <- data |>
    distinct(id) |>
    sample_frac(replace = TRUE)
  
  boot_df <- boot_ids |>
    mutate(new_id = row_number()) |>
    left_join(data, by = "id", relationship = "many-to-many") |>
    arrange(new_id, visit) |>
    dplyr::select(-id) |>                # drop original id
    dplyr::rename(id = new_id) 
  return(boot_df)
}

NR_spline = function(long.data, hat.parameters, knots, degree=3, tol = 1e-3, maxloop = 200){
  #long.data.id = long.data[!duplicated(long.data$id),]
  #time = c(long.data.id$preetime,pmin(long.data.id$etime,long.data.id$ctime))
  
  #d = c(0, quantile(time, probs = seq(1, knots)/knots, na.rm = TRUE))
  
  error = 0
  diff.theta = 1000
  # tol = 1e-3
  loop = 1
  # maxloop = 200
  while (diff.theta>tol & loop<=maxloop) {
    old.parameters = hat.parameters
    error=0
    #while (1) {
    diff = diff.likelihood.vec(long.data, hat.parameters)
    score = diff$Score
    #if (loop<=10){lr = 1}else if (loop > 10){lr = 5}
    lr = 5
    # Hessian = diff$Hessian
    Hessian = diff$Hessian
    if (any(is.na(Hessian))||any(is.infinite(Hessian))){
      error = 1
    }
    step = as.vector(Hessian%*%score)
    #likelihood.hat = likelihood.spline2(long.data, hat.parameters, knots)$ll
    likelihood.hat = likelihood.spline3(long.data, hat.parameters, knots)$ll
    while (1){
      est.parameters = hat.parameters - lr*t(step)
      # alpha.est = t(est.parameters %>% dplyr::select(contains("alpha")))
      # Sigma.e.est = matrix(0, nrow=K,ncol=K)
      # sigma.e.est.vec = as.matrix(est.parameters %>% dplyr::select(contains("Sigma.e")))
      # Sigma.e.est[upper.tri(Sigma.e,diag = TRUE)] = sigma.e.est.vec
      # if (K>1){
      #   Sigma.e.est=Sigma.e.est+t(Sigma.e.est)-diag(diag(Sigma.e.est))
      # }
      # eigen_e = eigen(Sigma.e.est, symmetric = TRUE, only.values = TRUE)$values
      # 
      # Sigma.a.est = matrix(0, nrow=K,ncol=K)
      # sigma.a.est.vec = as.matrix(est.parameters %>% dplyr::select(contains("Sigma.a")))
      # Sigma.a.est[upper.tri(Sigma.a,diag = TRUE)] = sigma.a.est.vec
      # if (K>1){
      #   Sigma.a.est=Sigma.a.est+t(Sigma.a.est)-diag(diag(Sigma.a.est))
      # }
      # eigen_a = eigen(Sigma.a.est, symmetric = TRUE, only.values = TRUE)$values
      # if (any(eigen_e<=0)|any(eigen_a<=0)){
      #   lr = lr/2
      #   next
      # }

      # likelihood.est = likelihood.spline2(long.data, est.parameters, knots)$ll
      likelihood.est = likelihood.spline3(long.data, est.parameters, knots)$ll
      if (is.na(likelihood.est) || is.infinite(likelihood.est)) {
        lr = lr/2
        next
      }
      if (lr*max(abs(step))<1e-6){
        temp.parameters = hat.parameters
        hat.parameters = est.parameters
        break
      }
      if (likelihood.est<likelihood.hat){
        lr = lr/2
        next
      }
      else {
        temp.parameters = hat.parameters
        hat.parameters = est.parameters
        break
      }
    }
    # if (max(abs(temp.parameters-hat.parameters))<tol) {
    #   break
    # }
    #}
    if (error == 1){
      hat.parameters = matrix(NA, nrow=1, ncol=length(temp.parameters))
      coef = matrix(NA, nrow=1, ncol=length(temp.parameters))
      stddev = matrix(NA, nrow=1, ncol=length(parameters))
      CP = matrix(NA, nrow=1, ncol=length(parameters))
      break
    }
    # diff.theta = max(abs(hat.parameters-old.parameters))
    # if (any(old.parameters) == 0){
    #   diff.theta = max(abs(hat.parameters-old.parameters)/(abs(old.parameters)+1e-6))
    # }else{
    #   diff.theta = max(abs(hat.parameters-old.parameters)/abs(old.parameters))
    # }
    diff.theta = max(abs(hat.parameters-old.parameters)/pmax(abs(old.parameters),1e-3))
    max_grad <- max(abs(score))
    
    # 2. Max Absolute Step (The raw step vector magnitude)
    max_step <- max(abs(step))
    
    # 3. Print cleanly using sprintf for scientific notation
    # cat(sprintf("Iter:%3d | MaxGrad: %.2f | MaxStep: %.2f | LR: %.2f | Diff: %.2f | ll: %.2f\n", 
    #             loop, max_grad, max_step, lr, diff.theta, likelihood.est))
    loop = loop + 1
  }
  return(hat.parameters)
}

