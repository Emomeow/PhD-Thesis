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
  
  lambda.hat = as.matrix(parameters %>% dplyr::select(contains("hazard")))
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
  
  a.re = nodes%*%sqrtm(2*Sigma.a.hat)
  # create X, v with x1 x2
  long.data.id <- long.data[!duplicated(long.data$id), ]
  x = long.data.id %>% dplyr::select(contains("xval"))
  X1 = x
  X1[, paste0("xval_", index)] = rep(x1, nrow(x))
  X2 = x
  X2[, paste0("xval_", index)] = rep(x2, nrow(x))
  
  J=100
  time1 <- seq(0, t, length.out = J+1)  # fixed time points
  # generate M1, M2
  mu1 = rep(1,(J+1)*n.id) %*% t(beta.0.hat) + rep(time1, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X1), matrix(1,J+1,1)) %*% t(beta.2.hat)
  Z <- mvrnorm((J+1)*n.id, rep(0,K), diag(K))
  L <- chol(Sigma.e.hat)
  M1_norand <- mu1 + Z %*% t(L)
  id = kronecker(long.data.id$id, matrix(1,J+1,1))
  M1_with_time <- data.frame(id = id, time1 = rep(time1,n.id), M1_norand)
  
  interval = sort(c(t, d[(d<t)]))
  fixed_grid <- expand.grid(id = unique(long.data$id), time1 = interval)
  M1_expanded <- full_join(M1_with_time, fixed_grid, by = c("id", "time1"))
  M1_expanded <- M1_expanded %>%
    arrange(id, time1) %>%
    fill(everything(), .direction = "down")
  time1.full <- M1_expanded$time1
  # create hazard
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
  
  M1_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M1_expanded[,-c(1,2)]))
  a.re.full <- a.re[rep(seq_len(nrow(a.re)),each = nrow(M1_expanded)), ]
  M1.full <- M1_norand.full + a.re.full
  
  
  time2 <- seq(0, t, length.out = J+1)
  mu2 = rep(1,(J+1)*n.id) %*% t(beta.0.hat) + rep(time1, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X2), matrix(1,J+1,1)) %*% t(beta.2.hat)
  Z <- mvrnorm((J+1)*n.id, rep(0,K), diag(K))
  M2_norand <- mu2 + Z %*% t(L)
  
  M2_with_time <- data.frame(id = id, time2 = rep(time1, n.id), M2_norand)
  fixed_grid <- expand.grid(id = unique(long.data$id), time2 = interval)
  M2_expanded <- full_join(M2_with_time, fixed_grid, by = c("id", "time2"))
  M2_expanded <- M2_expanded %>%
    arrange(id, time2) %>%
    fill(everything(), .direction = "down")
  time2.full <- M2_expanded$time2
  
  M2_expanded_hazard <- M2_expanded %>%
    rowwise() %>%
    mutate(
      hazard = hazard_intervals$hazard[
        time2 >= hazard_intervals$start & time2 < hazard_intervals$end
      ]
    )
  hazard2 <- M2_expanded_hazard$hazard
  
  M2_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M2_expanded[,-c(1,2)]))
  M2.full <- M2_norand.full + a.re.full
  # calculate survival factor theta_xx, theta_aa, theta_mm
  theta.mm1 <- M1.full %*% t(theta.m.hat)
  theta.mm1_expanded <- matrix(theta.mm1, nrow = nrow(M1_expanded))
  theta.mm2 <- M2.full %*% t(theta.m.hat)
  theta.mm2_expanded <- matrix(theta.mm2, nrow = nrow(M2_expanded))
  theta.xx1 <- as.matrix(X1) %*% t(theta.x.hat)
  theta.xx2 <- as.matrix(X2) %*% t(theta.x.hat)
  theta.xx1_expanded <- kronecker(matrix(1,nrow(M1_expanded)/n.id,nrow(a.re)), theta.xx1)
  theta.xx2_expanded <- kronecker(matrix(1,nrow(M1_expanded)/n.id,nrow(a.re)), theta.xx2)
  theta.aa <- a.re.full %*% t(theta.a.hat)
  theta.aa_expanded <- matrix(theta.aa, nrow = nrow(M1_expanded))
  
  # sum by id
  M1_expanded$id <- as.factor(M1_expanded$id)
  id.matrix <- model.matrix(~ id - 1, data = M1_expanded)
  sm_expanded = Matrix(t(id.matrix),sparse = T)
  
  It1 = as.numeric(time1.full < t)
  It2 = as.numeric(time2.full < t)
  
  Sx1m1_random = exp(-sm_expanded%*%(It1*hazard1*diff(c(time1.full, 0))*exp(theta.xx1_expanded+theta.aa_expanded+theta.mm1_expanded)))
  Sx1m1 = mean(Sx1m1_random%*%weight)/sqrt(pi)^K
  
  Sx2m1_random = exp(-sm_expanded%*%(It1*hazard1*diff(c(time1.full, 0))*exp(theta.xx2_expanded+theta.aa_expanded+theta.mm1_expanded)))
  Sx2m1 = mean(Sx2m1_random%*%weight)/sqrt(pi)^K
  
  Sx2m2_random = exp(-sm_expanded%*%(It2*hazard2*diff(c(time2.full, 0))*exp(theta.xx2_expanded+theta.aa_expanded+theta.mm2_expanded)))
  Sx2m2 = mean(Sx2m2_random%*%weight)/sqrt(pi)^K
  
  NDE = Sx2m1 - Sx1m1
  NIE = Sx2m2 - Sx2m1
  
  return(list(NDE = NDE,NIE = NIE))
}

With_cov_true = function(long.data, t, x1, x2, parameters, index, par, a.re){
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
  x = long.data.id %>% dplyr::select(contains("xval"))
  X1 = x
  X1[, paste0("xval_", index)] = rep(x1, nrow(x))
  X2 = x
  X2[, paste0("xval_", index)] = rep(x2, nrow(x))
  
  J=100
  time1 <- seq(0, t, length.out = J+1)  # fixed time points
  # generate M1, M2
  mu1 = rep(1,(J+1)*n.id) %*% t(beta.0.hat) + rep(time1, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X1), matrix(1,J+1,1)) %*% t(beta.2.hat)
  Z <- mvrnorm((J+1)*n.id, rep(0,K), diag(K))
  L <- chol(Sigma.e.hat)
  M1_norand <- mu1 + Z %*% t(L)
  id = kronecker(long.data.id$id, matrix(1,J+1,1))
  M1_with_time <- data.frame(id = id, time1 = rep(time1,n.id), M1_norand)
  
  interval = sort(c(t, d[d<=t]))
  fixed_grid <- expand.grid(id = unique(long.data$id), time1 = interval)
  M1_expanded <- full_join(M1_with_time, fixed_grid, by = c("id", "time1"))
  M1_expanded <- M1_expanded %>%
    arrange(id, time1) %>%
    fill(everything(), .direction = "down")
  time1.full <- M1_expanded$time1
  # interval2 = d[-length(d)]
  # hazard_intervals = data.frame(
  #   start = as.vector(interval2),
  #   end = c(interval2[-1], Inf),
  #   hazard = as.vector(lambda.hat)
  # )
  M1_expanded_hazard <- M1_expanded %>%
    rowwise() %>%
    mutate(
      hazard = par[2]/par[1]*(time1/par[1])^(par[2]-1)
    )
  hazard1 <- M1_expanded_hazard$hazard
  
  # M1_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M1_expanded[,-c(1,2)]))
  M1_norand.full <- as.matrix(M1_expanded[,-c(1,2)])
  a.re.full <- as.matrix(a.re[rep(seq_len(nrow(a.re)),each = nrow(M1_expanded_hazard)/n.id), ])
  M1.full <- M1_norand.full + a.re.full
  
  
  time2 <- seq(0, t, length.out = J+1)
  mu2 = rep(1,(J+1)*n.id) %*% t(beta.0.hat) + rep(time2, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X2), matrix(1,J+1,1)) %*% t(beta.2.hat)
  Z <- mvrnorm((J+1)*n.id, rep(0,K), diag(K))
  M2_norand <- mu2 + Z %*% t(L)
  
  M2_with_time <- data.frame(id = id, time2 = rep(time2, n.id), M2_norand)
  fixed_grid <- expand.grid(id = unique(long.data$id), time2 = interval)
  M2_expanded <- full_join(M2_with_time, fixed_grid, by = c("id", "time2"))
  time2.full <- M2_expanded$time2
  
  M2_expanded <- M2_expanded %>%
    arrange(id, time2) %>%
    fill(everything(), .direction = "down")
  
  M2_expanded_hazard <- M2_expanded %>%
    rowwise() %>%
    mutate(
      hazard = par[2]/par[1]*(time2/par[1])^(par[2]-1)
    )
  hazard2 <- M2_expanded_hazard$hazard
  
  # M2_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M2_expanded[,-c(1,2)]))
  M2_norand.full <- as.matrix(M2_expanded[,-c(1,2)])
  M2.full <- M2_norand.full + a.re.full
  # calculate survival factor theta_xx, theta_aa, theta_mm
  theta.mm1 <- M1.full %*% t(theta.m.hat)
  # theta.mm1_expanded <- matrix(theta.mm1, nrow = nrow(M1_expanded))
  theta.mm2 <- M2.full %*% t(theta.m.hat)
  # theta.mm2_expanded <- matrix(theta.mm2, nrow = nrow(M2_expanded))
  theta.xx1 <- as.matrix(X1) %*% t(theta.x.hat)
  theta.xx2 <- as.matrix(X2) %*% t(theta.x.hat)
  # theta.xx1_expanded <- kronecker(matrix(1,nrow(M1_expanded_hazard)/n.id,1), theta.xx1)
  theta.xx1_expanded <- rep(theta.xx1, nrow(M1_expanded_hazard)/n.id)
  # theta.xx2_expanded <- kronecker(matrix(1,nrow(M1_expanded_hazard)/n.id,1), theta.xx2)
  theta.xx2_expanded <- rep(theta.xx2, nrow(M2_expanded_hazard)/n.id)
  theta.aa <- a.re.full %*% t(theta.a.hat)
  # theta.aa_expanded <- matrix(theta.aa, nrow = nrow(M1_expanded))
  
  # sum by id
  M1_expanded$id <- as.factor(M1_expanded$id)
  id.matrix <- model.matrix(~ id - 1, data = M1_expanded)
  sm_expanded = Matrix(t(id.matrix),sparse = T)
  
  It1 = as.numeric(time1.full < t)
  It2 = as.numeric(time2.full < t)
  
  # weight2 = rep(1/(no.gauss^K), no.gauss^K)
  Sx1m1_random = exp(-sm_expanded%*%(It1*hazard1*diff(c(time1.full, 0))*exp(theta.xx1_expanded+theta.aa+theta.mm1)))
  Sx1m1 = mean(Sx1m1_random)
  
  Sx2m1_random = exp(-sm_expanded%*%(It1*hazard1*diff(c(time1.full, 0))*exp(theta.xx2_expanded+theta.aa+theta.mm1)))
  Sx2m1 = mean(Sx2m1_random)
  
  Sx2m2_random = exp(-sm_expanded%*%(It2*hazard2*diff(c(time2.full, 0))*exp(theta.xx2_expanded+theta.aa+theta.mm2)))
  Sx2m2 = mean(Sx2m2_random)
  
  NDE = Sx2m1 - Sx1m1
  NIE = Sx2m2 - Sx2m1
  
  return(list(NDE = NDE,NIE = NIE))
}
# 4.5 seconds


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
  
  lambda.hat = as.matrix(parameters %>% dplyr::select(contains("hazard")))
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
  
  a.re = nodes%*%sqrtm(2*Sigma.a.hat)
  # create X, v with x1 x2
  long.data.id <- long.data[!duplicated(long.data$id), ]
  x = long.data.id %>% dplyr::select(contains("xval"))
  X1 = x
  X1[, paste0("xval_", index)] = rep(x1, nrow(x))
  X2 = x
  X2[, paste0("xval_", index)] = rep(x2, nrow(x))
  
  max_t = max(long.data.id$etime)
  J=100
  time1 <- seq(0, max_t, length.out = J+1)  # fixed time points
  # generate M1, M2
  mu1 = rep(1,(J+1)*n.id) %*% t(beta.0.hat) + rep(time1, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X1), matrix(1,J+1,1)) %*% t(beta.2.hat)
  Z <- mvrnorm((J+1)*n.id, rep(0,K), diag(K))
  L <- chol(Sigma.e.hat)
  M1_norand <- mu1 + Z %*% t(L)
  id = kronecker(long.data.id$id, matrix(1,J+1,1))
  M1_with_time <- data.frame(id = id, time1 = rep(time1,n.id), M1_norand)
  
  interval = sort(c(t, d))
  fixed_grid <- expand.grid(id = unique(long.data$id), time1 = interval)
  M1_expanded <- full_join(M1_with_time, fixed_grid, by = c("id", "time1"))
  M1_expanded <- M1_expanded %>%
    arrange(id, time1) %>%
    fill(everything(), .direction = "down")
  time1.full <- M1_expanded$time1
  # create hazard
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
  
  M1_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M1_expanded[,-c(1,2)]))
  a.re.full <- a.re[rep(seq_len(nrow(a.re)),each = nrow(M1_expanded)), ]
  M1.full <- M1_norand.full + a.re.full
  
  
  # time2 <- seq(0, t, length.out = J+1)
  mu2 = rep(1,(J+1)*n.id) %*% t(beta.0.hat) + rep(time1, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X2), matrix(1,J+1,1)) %*% t(beta.2.hat)
  Z <- mvrnorm((J+1)*n.id, rep(0,K), diag(K))
  M2_norand <- mu2 + Z %*% t(L)
  
  M2_with_time <- data.frame(id = id, time2 = rep(time1, n.id), M2_norand)
  fixed_grid <- expand.grid(id = unique(long.data$id), time2 = interval)
  M2_expanded <- full_join(M2_with_time, fixed_grid, by = c("id", "time2"))
  M2_expanded <- M2_expanded %>%
    arrange(id, time2) %>%
    fill(everything(), .direction = "down")
  time2.full <- M2_expanded$time2
  
  M2_expanded_hazard <- M2_expanded %>%
    rowwise() %>%
    mutate(
      hazard = hazard_intervals$hazard[
        time2 >= hazard_intervals$start & time2 < hazard_intervals$end
      ]
    )
  hazard2 <- M2_expanded_hazard$hazard
  
  M2_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M2_expanded[,-c(1,2)]))
  M2.full <- M2_norand.full + a.re.full
  # calculate survival factor theta_xx, theta_aa, theta_mm
  theta.mm1 <- M1.full %*% t(theta.m.hat)
  theta.mm1_expanded <- matrix(theta.mm1, nrow = nrow(M1_expanded))
  theta.mm2 <- M2.full %*% t(theta.m.hat)
  theta.mm2_expanded <- matrix(theta.mm2, nrow = nrow(M2_expanded))
  theta.xx1 <- as.matrix(X1) %*% t(theta.x.hat)
  theta.xx2 <- as.matrix(X2) %*% t(theta.x.hat)
  theta.xx1_expanded <- kronecker(matrix(1,nrow(M1_expanded)/n.id,nrow(a.re)), theta.xx1)
  theta.xx2_expanded <- kronecker(matrix(1,nrow(M1_expanded)/n.id,nrow(a.re)), theta.xx2)
  theta.aa <- a.re.full %*% t(theta.a.hat)
  theta.aa_expanded <- matrix(theta.aa, nrow = nrow(M1_expanded))
  
  # generate survival time T
  diff_time = c(0, diff(time1.full))
  diff_time[diff_time<0] = 0
  
  # 1. Create the small building block *once*
  # This is the lower-triangular matrix for a single individual
  block_matrix <- tril(matrix(1, nrow(M1_expanded)/n.id, nrow(M1_expanded)/n.id))
  # 2. Create a list of these blocks, one for each individual
  list_of_blocks <- rep(list(block_matrix), n.id)
  
  # 3. Combine them into a sparse block-diagonal matrix directly
  # This avoids the huge dense intermediate object.
  cumsum_mat <- bdiag(list_of_blocks)
  # cumsum_mat = Matrix(kronecker(diag(n.id), lower.tri(matrix(1, nrow(M1_expanded)/n.id, nrow(M1_expanded)/n.id), diag = TRUE)+0), sparse = T)
  S1_ijb = exp(-cumsum_mat%*%(hazard1*diff_time*exp(theta.xx1_expanded+theta.aa_expanded+theta.mm1_expanded)))
  S2_ijb = exp(-cumsum_mat%*%(hazard2*diff_time*exp(theta.xx2_expanded+theta.aa_expanded+theta.mm2_expanded)))
  u = rep(runif(n.id), each = nrow(M1_expanded)/n.id)
  time_full = kronecker(t(rep(1, no.gauss^K)), time1.full)
  reformat_S1 = matrix(abs(1-S1_ijb-u), nrow = nrow(M1_expanded)/n.id, ncol = n.id*no.gauss^K)
  reformat_S2 = matrix(abs(1-S2_ijb-u), nrow = nrow(M1_expanded)/n.id, ncol = n.id*no.gauss^K)
  reformat_t = matrix(time_full, nrow = nrow(M1_expanded)/n.id, ncol = n.id*no.gauss^K)
  T1_ijb = matrix(reformat_t[apply(reformat_S1, 2, which.min)], nrow = n.id, ncol = no.gauss^K)
  T2_ijb = matrix(reformat_t[apply(reformat_S2, 2, which.min)], nrow = n.id, ncol = no.gauss^K)
  
  #Calculate M2
  M1_t = M1_expanded[M1_expanded$time1 == t, -c(1,2)]
  M2_t = M2_expanded[M2_expanded$time2 == t, -c(1,2)]
  # reformat
  M1_t.full = kronecker(as.matrix(M1_t), t(rep(1,no.gauss^K)))
  M2_t.full = kronecker(as.matrix(M2_t), t(rep(1,no.gauss^K)))

  M211_full = M1_t.full + kronecker(t(gamma.hat), (t > T1_ijb)*(t - T1_ijb))
  M212_full = M1_t.full + kronecker(t(gamma.hat), (t > T2_ijb)*(t - T2_ijb))
  M222_full = M2_t.full + kronecker(t(gamma.hat), (t > T2_ijb)*(t - T2_ijb))
  # NDE and NIE
  NIE = apply((M212_full - M211_full)%*%kronecker(diag(K), weight), 2, mean)
  NDE = apply((M222_full - M212_full)%*%kronecker(diag(K), weight), 2, mean)

  
  return(list(NDE = t(NDE),NIE = t(NIE)))
}

NE_M2_true = function(long.data, t, x1, x2, parameters, index, par, a.re){
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
  x = long.data.id %>% dplyr::select(contains("xval"))
  X1 = x
  X1[, paste0("xval_", index)] = rep(x1, nrow(x))
  X2 = x
  X2[, paste0("xval_", index)] = rep(x2, nrow(x))
  
  max_t = max(long.data.id$etime)
  J=100
  time1 <- seq(0, max_t, length.out = J+1)  # fixed time points
  # generate M1, M2
  mu1 = rep(1,(J+1)*n.id) %*% t(beta.0.hat) + rep(time1, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X1), matrix(1,J+1,1)) %*% t(beta.2.hat)
  Z <- mvrnorm((J+1)*n.id, rep(0,K), diag(K))
  L <- chol(Sigma.e.hat)
  M1_norand <- mu1 + Z %*% t(L)
  id = kronecker(long.data.id$id, matrix(1,J+1,1))
  M1_with_time <- data.frame(id = id, time1 = rep(time1,n.id), M1_norand)
  
  interval = sort(c(t, d[d<=t]))
  fixed_grid <- expand.grid(id = unique(long.data$id), time1 = interval)
  M1_expanded <- full_join(M1_with_time, fixed_grid, by = c("id", "time1"))
  M1_expanded <- M1_expanded %>%
    arrange(id, time1) %>%
    fill(everything(), .direction = "down")
  time1.full <- M1_expanded$time1
  # interval2 = d[-length(d)]
  # hazard_intervals = data.frame(
  #   start = as.vector(interval2),
  #   end = c(interval2[-1], Inf),
  #   hazard = as.vector(lambda.hat)
  # )
  M1_expanded_hazard <- M1_expanded %>%
    rowwise() %>%
    mutate(
      hazard = par[2]/par[1]*(time1/par[1])^(par[2]-1)
    )
  hazard1 <- M1_expanded_hazard$hazard
  
  # M1_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M1_expanded[,-c(1,2)]))
  M1_norand.full <- as.matrix(M1_expanded[,-c(1,2)])
  a.re.full <- as.matrix(a.re[rep(seq_len(nrow(a.re)),each = nrow(M1_expanded_hazard)/n.id), ])
  M1.full <- M1_norand.full + a.re.full
  
  
  time2 <- seq(0, max_t, length.out = J+1)
  mu2 = rep(1,(J+1)*n.id) %*% t(beta.0.hat) + rep(time2, n.id) %*% t(beta.1.hat) + kronecker(as.matrix(X2), matrix(1,J+1,1)) %*% t(beta.2.hat)
  Z <- mvrnorm((J+1)*n.id, rep(0,K), diag(K))
  M2_norand <- mu2 + Z %*% t(L)
  
  M2_with_time <- data.frame(id = id, time2 = rep(time2, n.id), M2_norand)
  fixed_grid <- expand.grid(id = unique(long.data$id), time2 = interval)
  M2_expanded <- full_join(M2_with_time, fixed_grid, by = c("id", "time2"))
  time2.full <- M2_expanded$time2
  
  M2_expanded <- M2_expanded %>%
    arrange(id, time2) %>%
    fill(everything(), .direction = "down")
  
  M2_expanded_hazard <- M2_expanded %>%
    rowwise() %>%
    mutate(
      hazard = par[2]/par[1]*(time2/par[1])^(par[2]-1)
    )
  hazard2 <- M2_expanded_hazard$hazard
  
  # M2_norand.full <- kronecker(matrix(1,nrow(a.re),1),as.matrix(M2_expanded[,-c(1,2)]))
  M2_norand.full <- as.matrix(M2_expanded[,-c(1,2)])
  M2.full <- M2_norand.full + a.re.full
  # calculate survival factor theta_xx, theta_aa, theta_mm
  theta.mm1 <- M1.full %*% t(theta.m.hat)
  # theta.mm1_expanded <- matrix(theta.mm1, nrow = nrow(M1_expanded))
  theta.mm2 <- M2.full %*% t(theta.m.hat)
  # theta.mm2_expanded <- matrix(theta.mm2, nrow = nrow(M2_expanded))
  theta.xx1 <- as.matrix(X1) %*% t(theta.x.hat)
  theta.xx2 <- as.matrix(X2) %*% t(theta.x.hat)
  # theta.xx1_expanded <- kronecker(matrix(1,nrow(M1_expanded_hazard)/n.id,1), theta.xx1)
  theta.xx1_expanded <- rep(theta.xx1, nrow(M1_expanded_hazard)/n.id)
  # theta.xx2_expanded <- kronecker(matrix(1,nrow(M1_expanded_hazard)/n.id,1), theta.xx2)
  theta.xx2_expanded <- rep(theta.xx2, nrow(M2_expanded_hazard)/n.id)
  theta.aa <- a.re.full %*% t(theta.a.hat)
  # theta.aa_expanded <- matrix(theta.aa, nrow = nrow(M1_expanded))
  
  diff_time = c(0, diff(time1.full))
  diff_time[diff_time<0] = 0
  
  # 1. Create the small building block *once*
  # This is the lower-triangular matrix for a single individual
  block_matrix <- tril(matrix(1, nrow(M1_expanded)/n.id, nrow(M1_expanded)/n.id))
  # 2. Create a list of these blocks, one for each individual
  list_of_blocks <- rep(list(block_matrix), n.id)
  
  # 3. Combine them into a sparse block-diagonal matrix directly
  # This avoids the huge dense intermediate object.
  cumsum_mat <- bdiag(list_of_blocks)
  # cumsum_mat = Matrix(kronecker(diag(n.id), lower.tri(matrix(1, nrow(M1_expanded)/n.id, nrow(M1_expanded)/n.id), diag = TRUE)+0), sparse = T)
  S1_ijb = exp(-cumsum_mat%*%(hazard1*diff_time*exp(theta.xx1_expanded+theta.aa+theta.mm1)))
  S2_ijb = exp(-cumsum_mat%*%(hazard2*diff_time*exp(theta.xx2_expanded+theta.aa+theta.mm2)))
  
  u = rep(runif(n.id), each = nrow(M1_expanded)/n.id)
  reformat_S1 = matrix(abs(1-S1_ijb-u), nrow = nrow(M1_expanded)/n.id, ncol = n.id)
  reformat_S2 = matrix(abs(1-S2_ijb-u), nrow = nrow(M1_expanded)/n.id, ncol = n.id)
  reformat_t = matrix(time1.full, nrow = nrow(M1_expanded)/n.id, ncol = n.id)
  T1_ijb = matrix(reformat_t[apply(reformat_S1, 2, which.min)], nrow = n.id, ncol = 1)
  T2_ijb = matrix(reformat_t[apply(reformat_S2, 2, which.min)], nrow = n.id, ncol = 1)
  
  M1_t = M1_expanded[M1_expanded$time1 == t, -c(1,2)]
  M2_t = M2_expanded[M2_expanded$time2 == t, -c(1,2)]
  # reformat
  M1_t.full = kronecker(as.matrix(M1_t), t(rep(1,no.gauss^K)))
  M2_t.full = kronecker(as.matrix(M2_t), t(rep(1,no.gauss^K)))
  
  M211 = M1_t + kronecker(t(gamma.hat), (t > T1_ijb)*(t - T1_ijb))
  M212 = M1_t + kronecker(t(gamma.hat), (t > T2_ijb)*(t - T2_ijb))
  M222 = M2_t + kronecker(t(gamma.hat), (t > T2_ijb)*(t - T2_ijb))
  
  NDE = apply(M222-M212, 2, mean)
  NIE = apply(M212-M211, 2, mean)
  
  return(list(NDE = t(NDE),NIE = t(NIE)))
}

NE_S = function(long.data, t, x1, x2, parameters, index, trials = 100){
  NE_S_matrix = matrix(0, nrow = trials, ncol=2)
  for (i in 1:trials){
    NE_S_i = With_cov(long.data, t, x1, x2, parameters, index)
    NE_S_matrix[i,1] = NE_S_i$NDE
    NE_S_matrix[i,2] = NE_S_i$NIE
  }
  NE_S_avg = apply(NE_S_matrix, 2, mean)
  NDE = NE_S_avg[1]
  NIE = NE_S_avg[2]
  return(list(NDE = NDE,NIE = NIE))
}

NE_M = function(long.data, t, x1, x2, parameters, index, trials = 100){
  NE_M_matrix = matrix(0, nrow = trials, ncol=2*K)
  for (i in 1:trials){
    NE_M_i = NE_M2(long.data, t, x1, x2, parameters, index)
    NE_M_matrix[i,1:K] = NE_M_i$NDE
    NE_M_matrix[i,(K+1):(2*K)] = NE_M_i$NIE
  }
  NE_M_avg = apply(NE_M_matrix, 2, mean)
  NDE = NE_S_avg[1:K]
  NIE = NE_S_avg[(K+1):(2*K)]
  return(list(NDE = NDE,NIE = NIE))
}