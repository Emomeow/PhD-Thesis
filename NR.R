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
    lr = 1
    # Hessian = diff$Hessian
    Hessian = diff$Hessian
    if (any(is.na(Hessian))||any(is.infinite(Hessian))){
      error = 1
    }
    step = as.vector(Hessian%*%score)
    likelihood.hat = likelihood.spline2(long.data, hat.parameters, knots)$ll
    while (1){
      est.parameters = hat.parameters - lr*t(step)
      # alpha.est = t(est.parameters %>% dplyr::select(contains("alpha")))
      sigma.e.est = est.parameters %>% dplyr::select(contains(paste("Sigma.e",1:K,1:K,sep="_")))
      sigma.a.est = est.parameters %>% dplyr::select(contains(paste("Sigma.a",1:K,1:K,sep="_")))
      if (any(sigma.e.est<=0)|any(sigma.a.est<=0)){
        lr = lr/2
        next
      }
      # test.time = seq(0, max(long.data$visittime), 0.01)
      # test.spline = evaluate_spline(test.time, knots, alpha.est)
      # if (any(test.spline < 0)){
      #   lr = lr/2
      #   next
      # }
      likelihood.est = likelihood.spline2(long.data, est.parameters, knots)$ll
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
    # diff.theta = max(abs(hat.parameters-old.parameters))
    diff.theta = max(abs(hat.parameters-old.parameters)/abs(old.parameters))
    loop = loop + 1
  }
  return(hat.parameters)
}

