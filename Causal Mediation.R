#############################
## functions for causal mediation analysis


# function for only one predictor (one x and no v)
No_cov = function(long.data, t, x1, x2, parameters){
  # import coefficients
  n.id = length(unique(long.data$id))
  beta.0.hat =  t(as.matrix(parameters %>% dplyr::select(contains("beta.0"))))
  K = length(beta.0.hat)
  beta.1.hat = t(as.matrix(parameters %>% dplyr::select(contains("beta.1"))))
  beta.2.vec = data.frame(parameters %>% dplyr::select(contains("beta.2")))
  col_names <- colnames(beta.2.vec)
  if (length(col_names)==1){
    beta.2.hat=as.matrix(beta.2.vec)
    p = ncol(beta.2.hat)
  } else{
    indices <- do.call(rbind, strsplit(col_names, "_"))[, 2:3]
    indices <- apply(indices, 2, as.integer)
    beta.2.hat = matrix(0, nrow=max(indices[,1]),ncol=max(indices[,2]))
    p=ncol(beta.2.hat)
    for (i in seq_along(col_names)){
      row_idx <- indices[i, 1]
      col_idx <- indices[i, 2]
      beta.2.hat[row_idx, col_idx] <- beta.2.vec[[col_names[i]]]
    }
  }
  gamma.hat = t(as.matrix(parameters %>% dplyr::select(contains("gamma"))))
  
  lambda.hat = as.matrix(parameters %>% dplyr::select(contains("hazard")))
  theta.x.hat = as.matrix(parameters %>% dplyr::select(contains("theta.x")))
  theta.a.hat = as.matrix(parameters %>% dplyr::select(contains("theta.a")))
  theta.m.hat = as.matrix(parameters %>% dplyr::select(contains("theta.m")))
  
  Sigma.a.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_a = as.matrix(parameters %>% dplyr::select(contains("Sigma.a")))
  Sigma.a.hat[upper.tri(Sigma.a,diag = TRUE)] = upper_vec_a
  if (K>1){
    Sigma.a.hat=Sigma.a.hat+t(Sigma.a.hat)-diag(diag(Sigma.a.hat))
  }
  Sigma.e.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_e = as.matrix(parameters %>% dplyr::select(contains("Sigma.e")))
  Sigma.e.hat[upper.tri(Sigma.e,diag = TRUE)] = upper_vec_e
  if (K>1){
    Sigma.e.hat=Sigma.e.hat+t(Sigma.e.hat)-diag(diag(Sigma.e.hat))
  }
  # inv.Sigma.e.hat = solve(Sigma.e.hat)
  
  a.re = nodes%*%sqrtm(2*Sigma.a.hat)
  
  # a.re.full <- a.re[rep(seq_len(nrow(a.re)),each = J+1), ]
  
  # generate J time points
  J=100
  time1 <- c(0, rexp(J, rate = 1/t))
  time1 <- sort(time1)
  
  mu1 = rep(1,(J+1)) %*% t(beta.0.hat) + time1 %*% t(beta.1.hat) + rep(x1, (J+1)) %*% t(beta.2.hat)
  Z <- mvrnorm(J+1, rep(0,K), diag(K))
  L <- chol(Sigma.e.hat)
  M1_norand <- mu1 + Z %*% t(L)
  ## creating M with x1 including knots
  M1_with_time <- as.data.frame(cbind(time1, M1_norand))
  interval = sort(c(t, d[-length(d)]))
  fixed_grid <- expand.grid(time1 = interval)
  M1_expanded <- full_join(M1_with_time, fixed_grid, by = c("time1"))
  time1.full <- M1_expanded$time1
  M1_expanded <- M1_expanded %>%
    arrange(time1) %>%
    fill(everything(), .direction = "down")
  
  interval2 = d[-length(d)]
  hazard_intervals = data.frame(
    start = as.vector(interval2),
    end = c(interval2[-1], Inf),
    hazard = as.vector(lambda.hat)
  )
  M1_expanded_hazard <- M1_expanded %>%
    rowwise() %>%
    mutate(
      hazard = hazard_intervals$hazard[
        time1 >= hazard_intervals$start & time1 < hazard_intervals$end
      ]
    )
  hazard1 <- M1_expanded_hazard$hazard
  
  M1_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M1_expanded[,-1]))
  
  a.re.full <- a.re[rep(seq_len(nrow(a.re)),each = nrow(M1_expanded)), ]
  
  M1.full <- M1_norand.full + a.re.full
  # create M with x2
  time2 <- c(0, rexp(J, rate = 1/t))
  time2 <- sort(time2)
  mu2 = rep(1,(J+1)) %*% t(beta.0.hat) + time2 %*% t(beta.1.hat) + rep(x2, (J+1)) %*% t(beta.2.hat)
  Z <- mvrnorm(J+1, rep(0,K), diag(K))
  M2_norand <- mu2 + Z %*% t(L)
  
  M2_with_time <- as.data.frame(cbind(time2, M2_norand))
  fixed_grid <- expand.grid(time2 = interval)
  M2_expanded <- full_join(M2_with_time, fixed_grid, by = c("time2"))
  time2.full <- M2_expanded$time2
  M2_expanded <- M2_expanded %>%
    arrange(time2) %>%
    fill(everything(), .direction = "down")
  
  M2_expanded_hazard <- M2_expanded %>%
    rowwise() %>%
    mutate(
      hazard = hazard_intervals$hazard[
        time2 >= hazard_intervals$start & time2 < hazard_intervals$end
      ]
    )
  hazard2 <- M2_expanded_hazard$hazard
  
  M2_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M2_expanded[,-1]))
  M2.full <- M2_norand.full + a.re.full
  # reshape theta_m*M
  theta.mm1 <- M1.full %*% t(theta.m.hat)
  theta.mm1_expanded <- matrix(theta.mm1, nrow = nrow(M1_expanded))
  theta.mm2 <- M2.full %*% t(theta.m.hat)
  theta.mm2_expanded <- matrix(theta.mm2, nrow = nrow(M2_expanded))
  

  # calculating S(x1), S(x2)
  theta.xx1 <- x1 %*% t(theta.x.hat)
  theta.xx2 <- x2 %*% t(theta.x.hat)
  theta.xx1_expanded <- kronecker(matrix(1,nrow(M1_expanded),nrow(a.re)), theta.xx1)
  theta.xx2_expanded <- kronecker(matrix(1,nrow(M1_expanded),nrow(a.re)), theta.xx2)
  
  theta.aa <- a.re.full %*% t(theta.a.hat)
  theta.aa_expanded <- matrix(theta.aa, nrow = nrow(M1_expanded))
  
  It1 = as.numeric(time1.full < t)
  It2 = as.numeric(time2.full < t)
  
  Sx1m1_random = exp(-t(hazard1*diff(c(time1.full, 0))*exp(theta.xx1_expanded+theta.aa_expanded+theta.mm1_expanded))%*%It1)
  Sx1m1 = sum(Sx1m1_random*weight)/sqrt(pi)^K
  
  Sx2m1_random = exp(-t(hazard1*diff(c(time1.full, 0))*exp(theta.xx2_expanded+theta.aa_expanded+theta.mm1_expanded))%*%It1)
  Sx2m1 = sum(Sx2m1_random*weight)/sqrt(pi)^K
  
  Sx1m2_random = exp(-t(hazard2*diff(c(time2.full, 0))*exp(theta.xx1_expanded+theta.aa_expanded+theta.mm2_expanded))%*%It2)
  Sx1m2 = sum(Sx1m2_random*weight)/sqrt(pi)^K
  
  # Sx2m2_random = exp(-t(hazard2*diff(c(time2.full, 0))*exp(theta.xx2_expanded+theta.aa_expanded+theta.mm2_expanded))%*%It2)
  # Sx2m2 = sum(Sx2m2_random*weight)/sqrt(pi)^K
  
  NDE = Sx2m1 - Sx1m1
  NIE = Sx1m2 - Sx1m1
  
  return(list("NDE=",NDE,"NIE=",NIE))
}


# function for one predictor and other covariates (one x and some v)
# 11/23/2025
# I-spline version
With_cov = function(long.data, t, x1, x2, parameters, index){
  n.id = length(unique(long.data$id))
  beta.0.hat =  t(as.matrix(parameters %>% dplyr::select(contains("beta.0"))))
  K = length(beta.0.hat)
  beta.1.hat = t(as.matrix(parameters %>% dplyr::select(contains("beta.1"))))
  beta.2.vec = data.frame(parameters %>% dplyr::select(contains("beta.2")))
  col_names <- colnames(beta.2.vec)
  if (length(col_names)==1){
    beta.2.hat=as.matrix(beta.2.vec)
    p = ncol(beta.2.hat)
  } else{
    indices <- do.call(rbind, strsplit(col_names, "_"))[, 2:3]
    indices <- apply(indices, 2, as.integer)
    beta.2.hat = matrix(0, nrow=max(indices[,1]),ncol=max(indices[,2]))
    p=ncol(beta.2.hat)
    for (i in seq_along(col_names)){
      row_idx <- indices[i, 1]
      col_idx <- indices[i, 2]
      beta.2.hat[row_idx, col_idx] <- beta.2.vec[[col_names[i]]]
    }
  }
  gamma.hat = t(as.matrix(parameters %>% dplyr::select(contains("gamma"))))
  
  # lambda.hat = as.matrix(parameters %>% dplyr::select(contains("hazard")))
  alpha.hat = t(as.matrix(parameters %>% dplyr::select(contains("alpha"))))
  theta.x.hat = as.matrix(parameters %>% dplyr::select(contains("theta.x")))
  theta.a.hat = as.matrix(parameters %>% dplyr::select(contains("theta.a")))
  theta.m.hat = as.matrix(parameters %>% dplyr::select(contains("theta.m")))
  
  Sigma.a.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_a = as.matrix(parameters %>% dplyr::select(contains("Sigma.a")))
  Sigma.a.hat[upper.tri(Sigma.a.hat,diag = TRUE)] = upper_vec_a
  if (K>1){
    Sigma.a.hat=Sigma.a.hat+t(Sigma.a.hat)-diag(diag(Sigma.a.hat))
  }
  Sigma.e.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_e = as.matrix(parameters %>% dplyr::select(contains("Sigma.e")))
  Sigma.e.hat[upper.tri(Sigma.e.hat,diag = TRUE)] = upper_vec_e
  if (K>1){
    Sigma.e.hat=Sigma.e.hat+t(Sigma.e.hat)-diag(diag(Sigma.e.hat))
  }
  # inv.Sigma.e.hat = solve(Sigma.e.hat)
  U = chol(Sigma.a.hat)
  # a.re = nodes%*%sqrtm(2*Sigma.a.hat)
  a.re = sqrt(2)*nodes%*%t(U)
  # create X, v with x1 x2
  long.data.id <- long.data[!duplicated(long.data$id), ]
  x = long.data.id %>% dplyr::select(contains("xval"))
  X1 = x
  X1[, paste0("xval_", index)] = rep(x1, nrow(x))
  X2 = x
  X2[, paste0("xval_", index)] = rep(x2, nrow(x))
  
  
  J=100
  id = as.factor(kronecker(long.data.id$id, matrix(1,J+1,1)))
  id.matrix <- model.matrix(~ id - 1)
  sm_expanded = Matrix(t(id.matrix),sparse = T)
  time1 <- seq(0, t, length.out = J+1)  # fixed time points
  # generate M1, M2
  mu1 = rep(1,(J+1)*n.id) %*% t(beta.0.hat) + rep(time1, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X1), matrix(1,J+1,1)) %*% t(beta.2.hat)
  time2 <- seq(0, t, length.out = J+1)
  mu2 = rep(1,(J+1)*n.id) %*% t(beta.0.hat) + rep(time1, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X2), matrix(1,J+1,1)) %*% t(beta.2.hat)
  a.re.full <- a.re[rep(seq_len(nrow(a.re)),each = nrow(mu1)), ]
  L <- chol(Sigma.e.hat)
  
  theta.xx1 <- as.matrix(X1) %*% t(theta.x.hat)
  theta.xx2 <- as.matrix(X2) %*% t(theta.x.hat)
  theta.xx1_expanded <- kronecker(theta.xx1, matrix(1,nrow(mu1)/n.id,nrow(a.re)))
  theta.xx2_expanded <- kronecker(theta.xx2, matrix(1,nrow(mu2)/n.id,nrow(a.re)))
  theta.aa <- a.re.full %*% t(theta.a.hat)
  theta.aa_expanded <- matrix(theta.aa, nrow = nrow(mu1))
  factor.ax1 = exp(theta.xx1_expanded+theta.aa_expanded)
  factor.ax2 = exp(theta.xx2_expanded+theta.aa_expanded)
  # hazard
  hazard = spline_cumulative_hazard(time1, knots, alpha.hat, boundary_knots = boundary_knots, degree=degree)
  hazard.full = rep(hazard, n.id)
  It1 = as.numeric(time1 < t)
  It2 = as.numeric(time2 < t)
  It1.full = rep(It1, n.id)
  It2.full = rep(It2, n.id)
  ## generate M_1(X1), M_1(X2) by MC
  MC_sample = 100
  Sx1m1_random = 0
  Sx2m1_random = 0
  Sx2m2_random = 0
  for (i in 1:MC_sample){
  Z1 <- mvrnorm((J+1)*n.id, rep(0,K), diag(K))
  # L <- chol(Sigma.e.hat)
  M1_norand <- mu1 + Z1 %*% t(L)
  # M1_with_time <- data.frame(id = id, time1 = rep(time1,n.id), M1_norand)
  
  # interval = sort(c(t, d[(d<t)]))
  # fixed_grid <- expand.grid(id = unique(long.data$id), time1 = interval)
  # M1_expanded <- full_join(M1_with_time, fixed_grid, by = c("id", "time1"))
  # M1_expanded <- M1_expanded %>%
  #   arrange(id, time1) %>%
  #   fill(everything(), .direction = "down")
  # time1.full <- M1_expanded$time1
  # # create hazard
  # interval2 = d[-length(d)]
  # hazard_intervals = data.frame(
  #   start = as.vector(interval2),
  #   end = c(interval2[-1], Inf),
  #   hazard = as.vector(lambda.hat)
  # )
  # M1_expanded_hazard <- M1_expanded %>%
  #   rowwise() %>%
  #   mutate(
  #     hazard = hazard_intervals$hazard[
  #       time1 >= hazard_intervals$start & time1 < hazard_intervals$end
  #     ]
  #   )
  # hazard1 <- M1_expanded_hazard$hazard
  
  M1_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M1_norand))
  # a.re.full <- a.re[rep(seq_len(nrow(a.re)),each = nrow(M1_norand)), ]
  M1.full <- M1_norand.full + a.re.full
  
  
  # time2 <- seq(0, t, length.out = J+1)
  # mu2 = rep(1,(J+1)*n.id) %*% t(beta.0.hat) + rep(time1, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X2), matrix(1,J+1,1)) %*% t(beta.2.hat)
  Z2 <- mvrnorm((J+1)*n.id, rep(0,K), diag(K))
  M2_norand <- mu2 + Z2 %*% t(L)
  
  # M2_with_time <- data.frame(id = id, time2 = rep(time1, n.id), M2_norand)
  # fixed_grid <- expand.grid(id = unique(long.data$id), time2 = interval)
  # M2_expanded <- full_join(M2_with_time, fixed_grid, by = c("id", "time2"))
  # M2_expanded <- M2_expanded %>%
  #   arrange(id, time2) %>%
  #   fill(everything(), .direction = "down")
  # time2.full <- M2_expanded$time2
  # 
  # M2_expanded_hazard <- M2_expanded %>%
  #   rowwise() %>%
  #   mutate(
  #     hazard = hazard_intervals$hazard[
  #       time2 >= hazard_intervals$start & time2 < hazard_intervals$end
  #     ]
  #   )
  # hazard2 <- M2_expanded_hazard$hazard
  
  M2_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M2_norand))
  M2.full <- M2_norand.full + a.re.full
  # calculate survival factor theta_xx, theta_aa, theta_mm
  theta.mm1 <- M1.full %*% t(theta.m.hat)
  theta.mm1_expanded <- matrix(theta.mm1, nrow = nrow(M1_norand))
  theta.mm2 <- M2.full %*% t(theta.m.hat)
  theta.mm2_expanded <- matrix(theta.mm2, nrow = nrow(M2_norand))
  # theta.xx1 <- as.matrix(X1) %*% t(theta.x.hat)
  # theta.xx2 <- as.matrix(X2) %*% t(theta.x.hat)
  # theta.xx1_expanded <- kronecker(matrix(1,nrow(M1_norand)/n.id,nrow(a.re)), theta.xx1)
  # theta.xx2_expanded <- kronecker(matrix(1,nrow(M1_norand)/n.id,nrow(a.re)), theta.xx2)
  # theta.aa <- a.re.full %*% t(theta.a.hat)
  # theta.aa_expanded <- matrix(theta.aa, nrow = nrow(M1_norand))
  
  # # hazard
  # hazard1 = evaluate_Ispline(time1, knots, alpha.hat, boundary_knots = boundary_knots, degree=degree)
  # hazard1.full = rep(hazard1, n.id)
  # # sum by id
  # M1_expanded$id <- as.factor(M1_expanded$id)
  # id.matrix <- model.matrix(~ id - 1, data = M1_expanded)
  # sm_expanded = Matrix(t(id.matrix),sparse = T)
  
  # It1 = as.numeric(time1.full < t)
  # It2 = as.numeric(time2.full < t)
  Sx1m1_random_MC = exp(-sm_expanded%*%(It1.full*diff(c(hazard.full, 0))*exp(theta.mm1_expanded)*factor.ax1))
  Sx2m1_random_MC = exp(-sm_expanded%*%(It1.full*diff(c(hazard.full, 0))*exp(theta.mm1_expanded)*factor.ax2))
  Sx2m2_random_MC = exp(-sm_expanded%*%(It2.full*diff(c(hazard.full, 0))*exp(theta.mm2_expanded)*factor.ax2))
  Sx1m1_random = Sx1m1_random + Sx1m1_random_MC/MC_sample
  Sx2m1_random = Sx2m1_random + Sx2m1_random_MC/MC_sample
  Sx2m2_random = Sx2m2_random + Sx2m2_random_MC/MC_sample
  gc()
  }
  # Sx1m1_random = exp(-sm_expanded%*%(It1*diff(c(hazard1, 0))*exp(theta.xx1_expanded+theta.aa_expanded+theta.mm1_expanded)))
  Sx1m1 = mean(Sx1m1_random%*%weight)/sqrt(pi)^K
  
  # Sx2m1_random = exp(-sm_expanded%*%(It1*hazard1*diff(c(time1.full, 0))*exp(theta.xx2_expanded+theta.aa_expanded+theta.mm1_expanded)))
  Sx2m1 = mean(Sx2m1_random%*%weight)/sqrt(pi)^K
  
  # Sx2m2_random = exp(-sm_expanded%*%(It2*hazard2*diff(c(time2.full, 0))*exp(theta.xx2_expanded+theta.aa_expanded+theta.mm2_expanded)))
  Sx2m2 = mean(Sx2m2_random%*%weight)/sqrt(pi)^K
  
  NDE = Sx2m1 - Sx1m1
  NIE = Sx2m2 - Sx2m1
  
  return(list(NDE = NDE,NIE = NIE))
}

With_cov_true = function(long.data, t, x1, x2, parameters, index, par){
  n.id = length(unique(long.data$id))
  beta.0.hat =  t(as.matrix(parameters %>% dplyr::select(contains("beta.0"))))
  K = length(beta.0.hat)
  beta.1.hat = t(as.matrix(parameters %>% dplyr::select(contains("beta.1"))))
  beta.2.vec = data.frame(parameters %>% dplyr::select(contains("beta.2")))
  col_names <- colnames(beta.2.vec)
  if (length(col_names)==1){
    beta.2.hat=as.matrix(beta.2.vec)
    p = ncol(beta.2.hat)
  } else{
    indices <- do.call(rbind, strsplit(col_names, "_"))[, 2:3]
    indices <- apply(indices, 2, as.integer)
    beta.2.hat = matrix(0, nrow=max(indices[,1]),ncol=max(indices[,2]))
    p=ncol(beta.2.hat)
    for (i in seq_along(col_names)){
      row_idx <- indices[i, 1]
      col_idx <- indices[i, 2]
      beta.2.hat[row_idx, col_idx] <- beta.2.vec[[col_names[i]]]
    }
  }
  gamma.hat = t(as.matrix(parameters %>% dplyr::select(contains("gamma"))))
  
  #lambda.hat = as.matrix(parameters %>% dplyr::select(contains("hazard")))
  theta.x.hat = as.matrix(parameters %>% dplyr::select(contains("theta.x")))
  theta.a.hat = as.matrix(parameters %>% dplyr::select(contains("theta.a")))
  theta.m.hat = as.matrix(parameters %>% dplyr::select(contains("theta.m")))
  
  # Sigma.a.hat = matrix(0, nrow=K,ncol=K)
  # upper_vec_a = as.matrix(parameters %>% dplyr::select(contains("Sigma.a")))
  # Sigma.a.hat[upper.tri(Sigma.a,diag = TRUE)] = upper_vec_a
  # if (K>1){
  #   Sigma.a.hat=Sigma.a.hat+t(Sigma.a.hat)-diag(diag(Sigma.a.hat))
  # }
  Sigma.e.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_e = as.matrix(parameters %>% dplyr::select(contains("Sigma.e")))
  Sigma.e.hat[upper.tri(Sigma.e,diag = TRUE)] = upper_vec_e
  if (K>1){
    Sigma.e.hat=Sigma.e.hat+t(Sigma.e.hat)-diag(diag(Sigma.e.hat))
  }
  # inv.Sigma.e.hat = solve(Sigma.e.hat)
  
  # a.re = nodes%*%sqrtm(2*Sigma.a.hat)
  # create X, v with x1 x2
  long.data.id <- long.data[!duplicated(long.data$id), ]
  a.re = as.matrix(long.data.id %>% dplyr::select(contains("aval")))
  x = long.data.id %>% dplyr::select(contains("xval"))
  X1 = x
  X1[, paste0("xval_", index)] = rep(x1, nrow(x))
  X2 = x
  X2[, paste0("xval_", index)] = rep(x2, nrow(x))
  
  J=100
  id = as.factor(kronecker(long.data.id$id, matrix(1,J+1,1)))
  id.matrix <- model.matrix(~ id - 1)
  sm_expanded = Matrix(t(id.matrix),sparse = T)
  time1 <- seq(0, t, length.out = J+1)  # fixed time points
  # generate M1, M2
  mu1 = rep(1,(J+1)*n.id) %*% t(beta.0.hat) + rep(time1, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X1), matrix(1,J+1,1)) %*% t(beta.2.hat)
  time2 <- seq(0, t, length.out = J+1)
  mu2 = rep(1,(J+1)*n.id) %*% t(beta.0.hat) + rep(time1, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X2), matrix(1,J+1,1)) %*% t(beta.2.hat)
  a.re.full <- a.re[rep(seq_len(nrow(a.re)),each = J+1), ]
  L <- chol(Sigma.e.hat)
  
  theta.xx1 <- as.matrix(X1) %*% t(theta.x.hat)
  theta.xx2 <- as.matrix(X2) %*% t(theta.x.hat)
  theta.xx1_expanded <- kronecker(theta.xx1, matrix(1,J+1,1))
  theta.xx2_expanded <- kronecker(theta.xx2, matrix(1,J+1,1))
  theta.aa <- a.re.full %*% t(theta.a.hat)
  # theta.aa_expanded <- matrix(theta.aa, nrow = nrow(mu1))
  factor.ax1 = exp(theta.xx1_expanded+theta.aa)
  factor.ax2 = exp(theta.xx2_expanded+theta.aa)
  # hazard
  hazard = (time1/par[1])^par[2]
  hazard.full = rep(hazard, n.id)
  It1 = as.numeric(time1 < t)
  It2 = as.numeric(time2 < t)
  It1.full = rep(It1, n.id)
  It2.full = rep(It2, n.id)
  ## generate M_1(X1), M_1(X2) by MC
  MC_sample = 100
  Sx1m1_random = 0
  Sx2m1_random = 0
  Sx2m2_random = 0
  for (i in 1:MC_sample){
    Z1 <- mvrnorm((J+1)*n.id, rep(0,K), diag(K))
    # L <- chol(Sigma.e.hat)
    M1_norand <- mu1 + Z1 %*% t(L)
    id = as.factor(kronecker(long.data.id$id, matrix(1,J+1,1)))
    id.matrix <- model.matrix(~ id - 1)
    sm_expanded = Matrix(t(id.matrix),sparse = T)
    # M1_with_time <- data.frame(id = id, time1 = rep(time1,n.id), M1_norand)
    
    # interval = sort(c(t, d[(d<t)]))
    # fixed_grid <- expand.grid(id = unique(long.data$id), time1 = interval)
    # M1_expanded <- full_join(M1_with_time, fixed_grid, by = c("id", "time1"))
    # M1_expanded <- M1_expanded %>%
    #   arrange(id, time1) %>%
    #   fill(everything(), .direction = "down")
    # time1.full <- M1_expanded$time1
    # # create hazard
    # interval2 = d[-length(d)]
    # hazard_intervals = data.frame(
    #   start = as.vector(interval2),
    #   end = c(interval2[-1], Inf),
    #   hazard = as.vector(lambda.hat)
    # )
    # M1_expanded_hazard <- M1_expanded %>%
    #   rowwise() %>%
    #   mutate(
    #     hazard = hazard_intervals$hazard[
    #       time1 >= hazard_intervals$start & time1 < hazard_intervals$end
    #     ]
    #   )
    # hazard1 <- M1_expanded_hazard$hazard
    
    # M1_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M1_norand))
    # a.re.full <- a.re[rep(seq_len(nrow(a.re)),each = nrow(M1_norand)), ]
    M1.full <- M1_norand + a.re.full
    
    
    # time2 <- seq(0, t, length.out = J+1)
    # mu2 = rep(1,(J+1)*n.id) %*% t(beta.0.hat) + rep(time1, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X2), matrix(1,J+1,1)) %*% t(beta.2.hat)
    Z2 <- mvrnorm((J+1)*n.id, rep(0,K), diag(K))
    M2_norand <- mu2 + Z2 %*% t(L)
    
    # M2_with_time <- data.frame(id = id, time2 = rep(time1, n.id), M2_norand)
    # fixed_grid <- expand.grid(id = unique(long.data$id), time2 = interval)
    # M2_expanded <- full_join(M2_with_time, fixed_grid, by = c("id", "time2"))
    # M2_expanded <- M2_expanded %>%
    #   arrange(id, time2) %>%
    #   fill(everything(), .direction = "down")
    # time2.full <- M2_expanded$time2
    # 
    # M2_expanded_hazard <- M2_expanded %>%
    #   rowwise() %>%
    #   mutate(
    #     hazard = hazard_intervals$hazard[
    #       time2 >= hazard_intervals$start & time2 < hazard_intervals$end
    #     ]
    #   )
    # hazard2 <- M2_expanded_hazard$hazard
    
    # M2_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M2_norand))
    M2.full <- M2_norand + a.re.full
    # calculate survival factor theta_xx, theta_aa, theta_mm
    theta.mm1 <- M1.full %*% t(theta.m.hat)
    #theta.mm1_expanded <- matrix(theta.mm1, nrow = nrow(M1_norand))
    theta.mm2 <- M2.full %*% t(theta.m.hat)
    #theta.mm2_expanded <- matrix(theta.mm2, nrow = nrow(M2_norand))
    # theta.xx1 <- as.matrix(X1) %*% t(theta.x.hat)
    # theta.xx2 <- as.matrix(X2) %*% t(theta.x.hat)
    # theta.xx1_expanded <- kronecker(matrix(1,nrow(M1_norand)/n.id,nrow(a.re)), theta.xx1)
    # theta.xx2_expanded <- kronecker(matrix(1,nrow(M1_norand)/n.id,nrow(a.re)), theta.xx2)
    # theta.aa <- a.re.full %*% t(theta.a.hat)
    # theta.aa_expanded <- matrix(theta.aa, nrow = nrow(M1_norand))
    
    # # hazard
    # hazard1 = evaluate_Ispline(time1, knots, alpha.hat, boundary_knots = boundary_knots, degree=degree)
    # hazard1.full = rep(hazard1, n.id)
    # # sum by id
    # M1_expanded$id <- as.factor(M1_expanded$id)
    # id.matrix <- model.matrix(~ id - 1, data = M1_expanded)
    # sm_expanded = Matrix(t(id.matrix),sparse = T)
    
    # It1 = as.numeric(time1.full < t)
    # It2 = as.numeric(time2.full < t)
    Sx1m1_random_MC = exp(-sm_expanded%*%(It1.full*diff(c(hazard.full, 0))*exp(theta.mm1)*factor.ax1))
    Sx2m1_random_MC = exp(-sm_expanded%*%(It1.full*diff(c(hazard.full, 0))*exp(theta.mm1)*factor.ax2))
    Sx2m2_random_MC = exp(-sm_expanded%*%(It2.full*diff(c(hazard.full, 0))*exp(theta.mm2)*factor.ax2))
    Sx1m1_random = Sx1m1_random + Sx1m1_random_MC/MC_sample
    Sx2m1_random = Sx2m1_random + Sx2m1_random_MC/MC_sample
    Sx2m2_random = Sx2m2_random + Sx2m2_random_MC/MC_sample
    gc()
  }
  # Sx1m1_random = exp(-sm_expanded%*%(It1*hazard1*diff(c(time1.full, 0))*exp(theta.xx1_expanded+theta.aa+theta.mm1)))
  Sx1m1 = mean(Sx1m1_random)
  
  # Sx2m1_random = exp(-sm_expanded%*%(It1*hazard1*diff(c(time1.full, 0))*exp(theta.xx2_expanded+theta.aa+theta.mm1)))
  Sx2m1 = mean(Sx2m1_random)
  
  # Sx2m2_random = exp(-sm_expanded%*%(It2*hazard2*diff(c(time2.full, 0))*exp(theta.xx2_expanded+theta.aa+theta.mm2)))
  Sx2m2 = mean(Sx2m2_random)
  
  NDE = Sx2m1 - Sx1m1
  NIE = Sx2m2 - Sx2m1
  
  return(list(NDE = NDE,NIE = NIE))
}


################################
# NDE NIE on M2
NE_M2 = function(long.data, t, x1, x2, parameters, index){
  n.id = length(unique(long.data$id))
  beta.0.hat =  t(as.matrix(parameters %>% dplyr::select(contains("beta.0"))))
  K = length(beta.0.hat)
  beta.1.hat = t(as.matrix(parameters %>% dplyr::select(contains("beta.1"))))
  beta.2.vec = data.frame(parameters %>% dplyr::select(contains("beta.2")))
  col_names <- colnames(beta.2.vec)
  if (length(col_names)==1){
    beta.2.hat=as.matrix(beta.2.vec)
    p = ncol(beta.2.hat)
  } else{
    indices <- do.call(rbind, strsplit(col_names, "_"))[, 2:3]
    indices <- apply(indices, 2, as.integer)
    beta.2.hat = matrix(0, nrow=max(indices[,1]),ncol=max(indices[,2]))
    p=ncol(beta.2.hat)
    for (i in seq_along(col_names)){
      row_idx <- indices[i, 1]
      col_idx <- indices[i, 2]
      beta.2.hat[row_idx, col_idx] <- beta.2.vec[[col_names[i]]]
    }
  }
  gamma.hat = t(as.matrix(parameters %>% dplyr::select(contains("gamma"))))
  
  #lambda.hat = as.matrix(parameters %>% dplyr::select(contains("hazard")))
  alpha.hat = t(as.matrix(parameters %>% dplyr::select(contains("alpha"))))
  theta.x.hat = as.matrix(parameters %>% dplyr::select(contains("theta.x")))
  theta.a.hat = as.matrix(parameters %>% dplyr::select(contains("theta.a")))
  theta.m.hat = as.matrix(parameters %>% dplyr::select(contains("theta.m")))
  
  Sigma.a.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_a = as.matrix(parameters %>% dplyr::select(contains("Sigma.a")))
  Sigma.a.hat[upper.tri(Sigma.a.hat,diag = TRUE)] = upper_vec_a
  if (K>1){
    Sigma.a.hat=Sigma.a.hat+t(Sigma.a.hat)-diag(diag(Sigma.a.hat))
  }
  Sigma.e.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_e = as.matrix(parameters %>% dplyr::select(contains("Sigma.e")))
  Sigma.e.hat[upper.tri(Sigma.e.hat,diag = TRUE)] = upper_vec_e
  if (K>1){
    Sigma.e.hat=Sigma.e.hat+t(Sigma.e.hat)-diag(diag(Sigma.e.hat))
  }
  # inv.Sigma.e.hat = solve(Sigma.e.hat)
  
  #a.re = nodes%*%sqrtm(2*Sigma.a.hat)
  U = chol(Sigma.a.hat)
  a.re = sqrt(2)*nodes%*%t(U)
  # create X, v with x1 x2
  long.data.id <- long.data[!duplicated(long.data$id), ]
  x = long.data.id %>% dplyr::select(contains("xval"))
  X1 = x
  X1[, paste0("xval_", index)] = rep(x1, nrow(x))
  X2 = x
  X2[, paste0("xval_", index)] = rep(x2, nrow(x))
  
  max_t = max(pmin(long.data.id$etime,long.data.id$ctime))
  J=100
  time1 <- sort(c(t,seq(0, max_t, length.out = J+1)))  # fixed time points
  time1 <- unique(time1)
  time1.full = rep(time1, n.id*no.gauss^K)
  # generate M1, M2
  mu1 = rep(1,length(time1)*n.id) %*% t(beta.0.hat) + rep(time1, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X1), matrix(1,length(time1),1)) %*% t(beta.2.hat)
  mu1.full = kronecker(matrix(1,nrow(a.re),1),as.matrix(mu1))
  time2 <- sort(c(t,seq(0, max_t, length.out = J+1)))
  time2 <- unique(time2)
  time2.full = rep(time2, n.id*no.gauss^K)
  mu2 = rep(1,length(time2)*n.id) %*% t(beta.0.hat) + rep(time2, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X2), matrix(1,length(time2),1)) %*% t(beta.2.hat)
  mu2.full = kronecker(matrix(1,nrow(a.re),1),as.matrix(mu2))
  a.re.full <- a.re[rep(seq_len(nrow(a.re)),each = nrow(mu1)), ]
  noerror1.full = mu1.full + a.re.full
  noerror2.full = mu2.full + a.re.full
  
  L <- chol(Sigma.e.hat)
  
  theta.xx1 <- as.matrix(X1) %*% t(theta.x.hat)
  theta.xx2 <- as.matrix(X2) %*% t(theta.x.hat)
  theta.xx1_expanded <- kronecker(theta.xx1, matrix(1,nrow(mu1)/n.id,nrow(a.re)))
  theta.xx2_expanded <- kronecker(theta.xx2, matrix(1,nrow(mu2)/n.id,nrow(a.re)))
  theta.aa <- a.re.full %*% t(theta.a.hat)
  theta.aa_expanded <- matrix(theta.aa, nrow = nrow(mu1))
  factor.ax1 = exp(theta.xx1_expanded+theta.aa_expanded)
  factor.ax2 = exp(theta.xx2_expanded+theta.aa_expanded)
  # hazard
  hazard = spline_cumulative_hazard(time1, knots, alpha.hat, boundary_knots = boundary_knots, degree=degree)
  hazard.full = rep(hazard, n.id)
  diff_hazard = c(0, diff(hazard.full))
  diff_hazard[diff_hazard<0] = 0
  # It1 = as.numeric(time1 < t)
  # It2 = as.numeric(time2 < t)
  # It1.full = rep(It1, n.id)
  # It2.full = rep(It2, n.id)
  ## generate M_1(X1), M_1(X2) by MC
  MC_sample = 100
  M211_full = 0
  M221_full = 0
  M222_full = 0
  
  # 1. Create the small building block *once*
  # This is the lower-triangular matrix for a single individual
  block_matrix <- tril(matrix(1, length(time1), length(time1)))
  # 2. Create a list of these blocks, one for each individual
  list_of_blocks <- rep(list(block_matrix), n.id)
  # 3. Combine them into a sparse block-diagonal matrix directly
  # This avoids the huge dense intermediate object.
  cumsum_mat <- bdiag(list_of_blocks)
  
  for (i in 1:MC_sample){
    Z1 <- mvrnorm(length(time1)*n.id*no.gauss^K, rep(0,K), diag(K))
    # L <- chol(Sigma.e.hat)
    # M1_norand <- mu1 + Z1 %*% t(L)
    # id = as.factor(kronecker(long.data.id$id, matrix(1,J+1,1)))
    # id.matrix <- model.matrix(~ id - 1)
    # sm_expanded = Matrix(t(id.matrix),sparse = T)
    # M1_with_time <- data.frame(id = id, time1 = rep(time1,n.id), M1_norand)
    
    # interval = sort(c(t, d[(d<t)]))
    # fixed_grid <- expand.grid(id = unique(long.data$id), time1 = interval)
    # M1_expanded <- full_join(M1_with_time, fixed_grid, by = c("id", "time1"))
    # M1_expanded <- M1_expanded %>%
    #   arrange(id, time1) %>%
    #   fill(everything(), .direction = "down")
    # time1.full <- M1_expanded$time1
    # # create hazard
    # interval2 = d[-length(d)]
    # hazard_intervals = data.frame(
    #   start = as.vector(interval2),
    #   end = c(interval2[-1], Inf),
    #   hazard = as.vector(lambda.hat)
    # )
    # M1_expanded_hazard <- M1_expanded %>%
    #   rowwise() %>%
    #   mutate(
    #     hazard = hazard_intervals$hazard[
    #       time1 >= hazard_intervals$start & time1 < hazard_intervals$end
    #     ]
    #   )
    # hazard1 <- M1_expanded_hazard$hazard
    
    # M1_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M1_norand))
    # a.re.full <- a.re[rep(seq_len(nrow(a.re)),each = nrow(M1_norand)), ]
    M1.full <- noerror1.full + Z1 %*% t(L)
    
    
    # time2 <- seq(0, t, length.out = J+1)
    # mu2 = rep(1,(J+1)*n.id) %*% t(beta.0.hat) + rep(time1, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X2), matrix(1,J+1,1)) %*% t(beta.2.hat)
    Z2 <- mvrnorm(length(time1)*n.id*no.gauss^K, rep(0,K), diag(K))
    # M2_norand <- mu2 + Z2 %*% t(L)
    
    # M2_with_time <- data.frame(id = id, time2 = rep(time1, n.id), M2_norand)
    # fixed_grid <- expand.grid(id = unique(long.data$id), time2 = interval)
    # M2_expanded <- full_join(M2_with_time, fixed_grid, by = c("id", "time2"))
    # M2_expanded <- M2_expanded %>%
    #   arrange(id, time2) %>%
    #   fill(everything(), .direction = "down")
    # time2.full <- M2_expanded$time2
    # 
    # M2_expanded_hazard <- M2_expanded %>%
    #   rowwise() %>%
    #   mutate(
    #     hazard = hazard_intervals$hazard[
    #       time2 >= hazard_intervals$start & time2 < hazard_intervals$end
    #     ]
    #   )
    # hazard2 <- M2_expanded_hazard$hazard
    
    # M2_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M2_norand))
    M2.full <- noerror2.full + Z2 %*% t(L)
    # calculate survival factor theta_xx, theta_aa, theta_mm
    theta.mm1 <- M1.full %*% t(theta.m.hat)
    theta.mm1_expanded <- matrix(theta.mm1, nrow = nrow(mu1))
    theta.mm2 <- M2.full %*% t(theta.m.hat)
    theta.mm2_expanded <- matrix(theta.mm2, nrow = nrow(mu2))
  
  # generate survival time T
  
  # 1. Create the small building block *once*
  # This is the lower-triangular matrix for a single individual
  # block_matrix <- tril(matrix(1, length(time1), length(time1)))
  # 2. Create a list of these blocks, one for each individual
  # list_of_blocks <- rep(list(block_matrix), n.id)
  
  # 3. Combine them into a sparse block-diagonal matrix directly
  # This avoids the huge dense intermediate object.
  # cumsum_mat <- bdiag(list_of_blocks)
  # cumsum_mat = Matrix(kronecker(diag(n.id), lower.tri(matrix(1, nrow(M1_expanded)/n.id, nrow(M1_expanded)/n.id), diag = TRUE)+0), sparse = T)
  # S1_matrix = diff_hazard*exp(theta.mm1_expanded)*factor.ax1
  # S1_rmatrix = matrix(diff_hazard*exp(theta.mm1_expanded)*factor.ax1, nrow = length(time1), ncol = n.id*no.gauss^K)
  # S1_ijb = matrix(exp(-apply(S1_rmatrix,2,cumsum)), nrow = n.id*length(time1), ncol=no.gauss^K)
  S1_ijb = exp(-cumsum_mat%*%(diff_hazard*exp(theta.mm1_expanded)*factor.ax1))
  S2_ijb = exp(-cumsum_mat%*%(diff_hazard*exp(theta.mm2_expanded)*factor.ax2))
  u1 = matrix(rep(runif(n.id*no.gauss^K), each = length(time1)), nrow = n.id*length(time1), ncol = no.gauss^K)
  u2 = matrix(rep(runif(n.id*no.gauss^K), each = length(time2)), nrow = n.id*length(time2), ncol = no.gauss^K)
  # time_full = matrix(time1.full, ncol = no.gauss^K)
  reformat_S1 = matrix(abs(1-S1_ijb-u1), nrow = length(time1), ncol = n.id*no.gauss^K)
  reformat_S2 = matrix(abs(1-S2_ijb-u2), nrow = length(time2), ncol = n.id*no.gauss^K)
  reformat_t = matrix(time1.full, nrow = length(time1), ncol = n.id*no.gauss^K)
  T1_ijb = matrix(reformat_t[apply(reformat_S1, 2, which.min)], nrow = n.id, ncol = no.gauss^K)
  T2_ijb = matrix(reformat_t[apply(reformat_S2, 2, which.min)], nrow = n.id, ncol = no.gauss^K)
  
  #Calculate M2
  M1_t = M1.full[time1.full == t, ]
  M2_t = M2.full[time2.full == t, ]
  # reformat
  row_index <- seq_len(nrow(M1_t))
  groups <- rep(seq_len(no.gauss^K), each = n.id)
  
  # Split *indices*, then subset x → preserves row structure
  blocks1 <- lapply(
    split(row_index, groups),
    function(i) M1_t[i, , drop = FALSE]   # each block: 400 × 2
  )
  blocks2 <- lapply(
    split(row_index, groups),
    function(i) M2_t[i, , drop = FALSE]   # each block: 400 × 2
  )
  # Column-bind all 400 blocks → 400 × 800
  M1_t.full <- do.call(cbind, blocks1)
  M2_t.full <- do.call(cbind, blocks2)
  # M1_t.full = kronecker(as.matrix(M1_t), t(rep(1,no.gauss^K)))
  # M2_t.full = kronecker(as.matrix(M2_t), t(rep(1,no.gauss^K)))

  M211_full_MC = M1_t.full + kronecker((t > T1_ijb)*(t - T1_ijb), t(gamma.hat))
  M221_full_MC = M2_t.full + kronecker((t > T1_ijb)*(t - T1_ijb), t(gamma.hat))
  M222_full_MC = M2_t.full + kronecker((t > T2_ijb)*(t - T2_ijb), t(gamma.hat))
  
  M211_full = M211_full + M211_full_MC/MC_sample
  M221_full = M221_full + M221_full_MC/MC_sample
  M222_full = M222_full + M222_full_MC/MC_sample
  # gc()
  }
  # NDE and NIE
  NIE = apply((M222_full - M221_full)%*%kronecker(weight, diag(K)), 2, mean)/(sqrt(pi)^K)
  NDE = apply((M221_full - M211_full)%*%kronecker(weight, diag(K)), 2, mean)/(sqrt(pi)^K)

  
  return(list(NDE = t(NDE),NIE = t(NIE)))
}

NE_M2_true = function(long.data, t, x1, x2, parameters, index, par){
  n.id = length(unique(long.data$id))
  beta.0.hat =  t(as.matrix(parameters %>% dplyr::select(contains("beta.0"))))
  K = length(beta.0.hat)
  beta.1.hat = t(as.matrix(parameters %>% dplyr::select(contains("beta.1"))))
  beta.2.vec = data.frame(parameters %>% dplyr::select(contains("beta.2")))
  col_names <- colnames(beta.2.vec)
  if (length(col_names)==1){
    beta.2.hat=as.matrix(beta.2.vec)
    p = ncol(beta.2.hat)
  } else{
    indices <- do.call(rbind, strsplit(col_names, "_"))[, 2:3]
    indices <- apply(indices, 2, as.integer)
    beta.2.hat = matrix(0, nrow=max(indices[,1]),ncol=max(indices[,2]))
    p=ncol(beta.2.hat)
    for (i in seq_along(col_names)){
      row_idx <- indices[i, 1]
      col_idx <- indices[i, 2]
      beta.2.hat[row_idx, col_idx] <- beta.2.vec[[col_names[i]]]
    }
  }
  gamma.hat = t(as.matrix(parameters %>% dplyr::select(contains("gamma"))))
  
  #lambda.hat = as.matrix(parameters %>% dplyr::select(contains("hazard")))
  theta.x.hat = as.matrix(parameters %>% dplyr::select(contains("theta.x")))
  theta.a.hat = as.matrix(parameters %>% dplyr::select(contains("theta.a")))
  theta.m.hat = as.matrix(parameters %>% dplyr::select(contains("theta.m")))
  
  Sigma.a.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_a = as.matrix(parameters %>% dplyr::select(contains("Sigma.a")))
  Sigma.a.hat[upper.tri(Sigma.a,diag = TRUE)] = upper_vec_a
  if (K>1){
    Sigma.a.hat=Sigma.a.hat+t(Sigma.a.hat)-diag(diag(Sigma.a.hat))
  }
  Sigma.e.hat = matrix(0, nrow=K,ncol=K)
  upper_vec_e = as.matrix(parameters %>% dplyr::select(contains("Sigma.e")))
  Sigma.e.hat[upper.tri(Sigma.e,diag = TRUE)] = upper_vec_e
  if (K>1){
    Sigma.e.hat=Sigma.e.hat+t(Sigma.e.hat)-diag(diag(Sigma.e.hat))
  }
  # inv.Sigma.e.hat = solve(Sigma.e.hat)
  
  # a.re = nodes%*%sqrtm(2*Sigma.a.hat)
  # create X, v with x1 x2
  long.data.id <- long.data[!duplicated(long.data$id), ]
  a.re = as.matrix(long.data.id %>% dplyr::select(contains("aval")))
  x = long.data.id %>% dplyr::select(contains("xval"))
  X1 = x
  X1[, paste0("xval_", index)] = rep(x1, nrow(x))
  X2 = x
  X2[, paste0("xval_", index)] = rep(x2, nrow(x))
  
  max_t = max(long.data.id$etime)
  J=100
  time1 <- sort(c(t,seq(0, max_t, length.out = J+1)))  # fixed time points
  time1 <- unique(time1)
  time1.full = rep(time1, n.id)
  # generate M1, M2
  mu1 = rep(1,length(time1)*n.id) %*% t(beta.0.hat) + rep(time1, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X1), matrix(1,length(time1),1)) %*% t(beta.2.hat)
  time2 <- sort(c(t,seq(0, max_t, length.out = J+1)))
  time2 <- unique(time2)
  time2.full = rep(time2, n.id)
  mu2 = rep(1,length(time2)*n.id) %*% t(beta.0.hat) + rep(time2, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X2), matrix(1,length(time2),1)) %*% t(beta.2.hat)
  # a.re.full <- a.re[rep(seq_len(nrow(a.re)),each = nrow(mu1)), ]
  a.re.full <- a.re[rep(seq_len(nrow(a.re)),each = length(time1)), ]
  L <- chol(Sigma.e.hat)
  
  theta.xx1 <- as.matrix(X1) %*% t(theta.x.hat)
  theta.xx2 <- as.matrix(X2) %*% t(theta.x.hat)
  theta.xx1_expanded <- kronecker(theta.xx1, matrix(1,length(time1),1))
  theta.xx2_expanded <- kronecker(theta.xx2, matrix(1,length(time2),1))
  theta.aa <- a.re.full %*% t(theta.a.hat)
  # theta.aa_expanded <- matrix(theta.aa, nrow = nrow(mu1))
  factor.ax1 = exp(theta.xx1_expanded+theta.aa)
  factor.ax2 = exp(theta.xx2_expanded+theta.aa)
  # hazard
  hazard = (time1/par[1])^par[2]
  hazard.full = rep(hazard, n.id)
  
  diff_hazard = c(0, diff(hazard.full))
  diff_hazard[diff_hazard<0] = 0
  # It1 = as.numeric(time1 < t)
  # It2 = as.numeric(time2 < t)
  # It1.full = rep(It1, n.id)
  # It2.full = rep(It2, n.id)
  ## generate M_1(X1), M_1(X2) by MC
  MC_sample = 100
  M211_full = 0
  M221_full = 0
  M222_full = 0
  
  # 1. Create the small building block *once*
  # This is the lower-triangular matrix for a single individual
  block_matrix <- tril(matrix(1, length(time1), length(time1)))
  # 2. Create a list of these blocks, one for each individual
  list_of_blocks <- rep(list(block_matrix), n.id)
  # 3. Combine them into a sparse block-diagonal matrix directly
  # This avoids the huge dense intermediate object.
  cumsum_mat <- bdiag(list_of_blocks)
  
  for (i in 1:MC_sample){
    Z1 <- mvrnorm(length(time1)*n.id, rep(0,K), diag(K))
    # L <- chol(Sigma.e.hat)
    M1_norand <- mu1 + Z1 %*% t(L)
    # id = as.factor(kronecker(long.data.id$id, matrix(1,J+1,1)))
    # id.matrix <- model.matrix(~ id - 1)
    # sm_expanded = Matrix(t(id.matrix),sparse = T)
    # M1_with_time <- data.frame(id = id, time1 = rep(time1,n.id), M1_norand)
    
    # interval = sort(c(t, d[(d<t)]))
    # fixed_grid <- expand.grid(id = unique(long.data$id), time1 = interval)
    # M1_expanded <- full_join(M1_with_time, fixed_grid, by = c("id", "time1"))
    # M1_expanded <- M1_expanded %>%
    #   arrange(id, time1) %>%
    #   fill(everything(), .direction = "down")
    # time1.full <- M1_expanded$time1
    # # create hazard
    # interval2 = d[-length(d)]
    # hazard_intervals = data.frame(
    #   start = as.vector(interval2),
    #   end = c(interval2[-1], Inf),
    #   hazard = as.vector(lambda.hat)
    # )
    # M1_expanded_hazard <- M1_expanded %>%
    #   rowwise() %>%
    #   mutate(
    #     hazard = hazard_intervals$hazard[
    #       time1 >= hazard_intervals$start & time1 < hazard_intervals$end
    #     ]
    #   )
    # hazard1 <- M1_expanded_hazard$hazard
    
    # M1_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M1_norand))
    # a.re.full <- a.re[rep(seq_len(nrow(a.re)),each = nrow(M1_norand)), ]
    M1.full <- M1_norand + a.re.full
    
    
    # time2 <- seq(0, t, length.out = J+1)
    # mu2 = rep(1,(J+1)*n.id) %*% t(beta.0.hat) + rep(time1, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X2), matrix(1,J+1,1)) %*% t(beta.2.hat)
    Z2 <- mvrnorm(length(time1)*n.id, rep(0,K), diag(K))
    M2_norand <- mu2 + Z2 %*% t(L)
    
    # M2_with_time <- data.frame(id = id, time2 = rep(time1, n.id), M2_norand)
    # fixed_grid <- expand.grid(id = unique(long.data$id), time2 = interval)
    # M2_expanded <- full_join(M2_with_time, fixed_grid, by = c("id", "time2"))
    # M2_expanded <- M2_expanded %>%
    #   arrange(id, time2) %>%
    #   fill(everything(), .direction = "down")
    # time2.full <- M2_expanded$time2
    # 
    # M2_expanded_hazard <- M2_expanded %>%
    #   rowwise() %>%
    #   mutate(
    #     hazard = hazard_intervals$hazard[
    #       time2 >= hazard_intervals$start & time2 < hazard_intervals$end
    #     ]
    #   )
    # hazard2 <- M2_expanded_hazard$hazard
    
    # M2_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M2_norand))
    M2.full <- M2_norand + a.re.full
    # calculate survival factor theta_xx, theta_aa, theta_mm
    theta.mm1 <- M1.full %*% t(theta.m.hat)
    # theta.mm1_expanded <- matrix(theta.mm1, nrow = nrow(M1_norand))
    theta.mm2 <- M2.full %*% t(theta.m.hat)
    # theta.mm2_expanded <- matrix(theta.mm2, nrow = nrow(M2_norand))
    
    # generate survival time T
    
    # 1. Create the small building block *once*
    # This is the lower-triangular matrix for a single individual
    # block_matrix <- tril(matrix(1, length(time1), length(time1)))
    # 2. Create a list of these blocks, one for each individual
    # list_of_blocks <- rep(list(block_matrix), n.id)
    
    # 3. Combine them into a sparse block-diagonal matrix directly
    # This avoids the huge dense intermediate object.
    # cumsum_mat <- bdiag(list_of_blocks)
    # cumsum_mat = Matrix(kronecker(diag(n.id), lower.tri(matrix(1, nrow(M1_expanded)/n.id, nrow(M1_expanded)/n.id), diag = TRUE)+0), sparse = T)
    S1_ijb = exp(-cumsum_mat%*%(diff_hazard*exp(theta.mm1)*factor.ax1))
    S2_ijb = exp(-cumsum_mat%*%(diff_hazard*exp(theta.mm2)*factor.ax2))
    u1 = rep(runif(n.id), each = length(time1))
    u2 = rep(runif(n.id), each = length(time2))
    # time_full = matrix(time1.full, ncol = no.gauss^K)
    reformat_S1 = matrix(abs(1-S1_ijb-u1), nrow = length(time1), ncol = n.id)
    reformat_S2 = matrix(abs(1-S2_ijb-u2), nrow = length(time2), ncol = n.id)
    reformat_t = matrix(time1.full, nrow = length(time1), ncol = n.id)
    T1_ijb = matrix(reformat_t[apply(reformat_S1, 2, which.min)], nrow = n.id, ncol = 1)
    T2_ijb = matrix(reformat_t[apply(reformat_S2, 2, which.min)], nrow = n.id, ncol = 1)
    
    #Calculate M2
    M1_t = M1.full[time1.full == t, ]
    M2_t = M2.full[time2.full == t, ]
    # # reformat
    # row_index <- seq_len(nrow(M1_t))
    # groups <- rep(seq_len(no.gauss^K), each = n.id)
    # 
    # # Split *indices*, then subset x → preserves row structure
    # blocks1 <- lapply(
    #   split(row_index, groups),
    #   function(i) M1_t[i, , drop = FALSE]   # each block: 400 × 2
    # )
    # blocks2 <- lapply(
    #   split(row_index, groups),
    #   function(i) M2_t[i, , drop = FALSE]   # each block: 400 × 2
    # )
    # # Column-bind all 400 blocks → 400 × 800
    # M1_t.full <- do.call(cbind, blocks1)
    # M2_t.full <- do.call(cbind, blocks2)
    # M1_t.full = kronecker(as.matrix(M1_t), t(rep(1,no.gauss^K)))
    # M2_t.full = kronecker(as.matrix(M2_t), t(rep(1,no.gauss^K)))
    
    M211_full_MC = M1_t + kronecker((t > T1_ijb)*(t - T1_ijb), t(gamma.hat))
    M221_full_MC = M2_t + kronecker((t > T1_ijb)*(t - T1_ijb), t(gamma.hat))
    M222_full_MC = M2_t + kronecker((t > T2_ijb)*(t - T2_ijb), t(gamma.hat))
    
    M211_full = M211_full + M211_full_MC/MC_sample
    M221_full = M221_full + M221_full_MC/MC_sample
    M222_full = M222_full + M222_full_MC/MC_sample
    gc()
  }
  # NDE and NIE
  NIE = apply((M222_full - M221_full), 2, mean)
  NDE = apply((M221_full - M211_full), 2, mean)
  
  return(list(NDE = t(NDE),NIE = t(NIE)))
}
